import functools
import itertools
import multiprocessing
import operator
import inspect

import numpy as np
import pandas as pd
import pyranges as pr
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

import handygenome
import handygenome.tools as tools
import handygenome.logutils as logutils
import handygenome.refgenome.refgenome as refgenome
import handygenome.deco as deco
import handygenome.ucscdata as ucscdata

from handygenome.genomedf.genomedf_base import GenomeDataFrameBase
import handygenome.plot.gui as gui


DOT_ALPHA_CONST = 0.01 / 3e6
GAPREGION_PREFIX = '::GAP::'


######################
# alpha interpolator #
######################

def make_interp1d(data):
    data = sorted(data, key=(lambda x: x[0]))
    xs, ys = zip(*data)
    xs = np.array(xs)
    ys = np.array(ys)
    interp = scipy.interpolate.interp1d(
        xs, ys, bounds_error=False, fill_value=(ys[0], ys[-1]),
    )

    return interp

DEPTH_ALPHA_DATA = [
    (3101808, 0.01),
    (249250, 0.1),
    (59129, 0.3),
    (20000, 0.7),
    (10000, 1),
]

BAF_ALPHA_DATA = [
    (3010, 1),
    (11732, 0.7),
    (147650, 0.3),
    (438827, 0.08),
    (1966549, 0.02),
]

ALPHA_INTERP_DEPTH = make_interp1d(DEPTH_ALPHA_DATA)
ALPHA_INTERP_BAF = make_interp1d(BAF_ALPHA_DATA)


def calc_dot_alpha_depth(n_dots):
    return float(ALPHA_INTERP_DEPTH(n_dots))


def calc_dot_alpha_baf(n_dots):
    return float(ALPHA_INTERP_BAF(n_dots))


def calc_dot_alpha_old(
    n_dots, 
    max_n_dots=3_000_000,
    min_n_dots=5_000,
    max_alpha=1,
    min_alpha=0.01,
):
    result = max_alpha + (n_dots - min_n_dots) * (
        (min_alpha - max_alpha) 
        / (max_n_dots - min_n_dots)
    )
    result = np.clip(result, min_alpha, max_alpha)
    return result



################
# main classes #
################

class CoordinateConverter:
    ########
    # init #
    ########

    def __init__(
        self, 
        refver, 

        chroms=None, 
        start0s=None, 
        end0s=None, 
        weights=1,

        region_gdf=None, 

        region_gaps=0.1,

        xmin=0,
        xmax=1,
    ):
        """Args:
            region_gaps: fraction of sum of gap lengths over plotting region lengths
        """
        region_gaps = float(region_gaps)
        assert (
            (region_gaps >= 0)
            and (region_gaps < 1)
        )

        self.refver = refver
        self.chromdict = refgenome.get_chromdict(self.refver)
        self.xmin = xmin
        self.xmax = xmax

        # input region_gdf sanitycheck
        if region_gdf is not None:
            region_gdf = region_gdf.copy()
            if 'weight' not in region_gdf.annot_cols:
                region_gdf['weight'] = weights

        # set region_gaps
        #if region_gaps is None:
        #    region_gaps = 0.1
        self.region_gaps = region_gaps

        # prepare unmodified region_gdf
        if region_gdf is None:
            if chroms is None:
                region_gdf = GenomeDataFrameBase.all_regions(self.refver, assembled_only=True)
                region_gdf['weight'] = weights
            else:
                region_gdf = GenomeDataFrameBase.from_data(
                    refver, 
                    chroms=chroms, 
                    start0s=start0s, 
                    end0s=end0s, 
                    weight=weights,
                )

        self.trim_coords(region_gdf)

        # insert gap regions into region_gdf
        region_gdf = self.insert_gap_to_region_gdf(region_gdf, region_gaps)

        # set attributes for internal use
        self.set_params(region_gdf)

    @staticmethod
    def trim_coords(region_gdf):
        chromdict = refgenome.get_chromdict(region_gdf.refver)
        region_gdf['Start'] = np.maximum(0, region_gdf.start0s)
        region_gdf['End'] = np.minimum([chromdict[x] for x in region_gdf.chroms], region_gdf.end0s)

    @staticmethod
    def insert_gap_to_region_gdf(gdf, region_gaps, last=False):
        new_gdf = gdf.copy().choose_annots('weight')
        assert not np.char.startswith(new_gdf.chroms.astype(str), GAPREGION_PREFIX).any()

        # intercalate with gap rows
        do_make_gaps = (
            (last and (region_gaps != 0))
            or (
                (not last) 
                and (region_gaps != 0) 
                and (new_gdf.nrow > 1)
            )
        )
        if do_make_gaps:
        #if (region_gaps != 0) and (new_gdf.nrow > 1):
            weighted_length_mean = (new_gdf.lengths * new_gdf['weight']).sum() / new_gdf.nrow
            gap_length = weighted_length_mean * region_gaps

            #total_gap_length = weighted_length_sum * region_gaps
            #gap_weight = total_gap_length / (new_gdf.nrow - 1)
                # since raw length of every gap region is 1, weight is equal to weighted length

            nrow = new_gdf.nrow * 2 - 1
            num_gap = new_gdf.nrow - 1
            if last:
                nrow += 1
                num_gap += 1

            src_arr = np.empty((nrow,  new_gdf.ncol), dtype=object)
            src_arr[0::2, :] = new_gdf.df
            src_arr[1::2, 0] = [f'{GAPREGION_PREFIX}{x}' for x in range(num_gap)]  # gap region Chromosome
            src_arr[1::2, 1] = 0  # gap region Start
            src_arr[1::2, 2] = 1  # gap region End
            src_arr[1::2, 3] = gap_length  # gap region weight is equal to length because start and end is 0 and 1 for every gap region

            new_df = pd.DataFrame(src_arr, columns=new_gdf.columns)
            new_gdf.assign_frame(new_df, dtype={'weight': float})

        return new_gdf

    def set_params(self, region_gdf):
        # check internal overlap
        if region_gdf.check_self_overlap():
            raise Exception(f'Plot region dataframe must not have overlapping intervals.')

        # set totalregion_gdf
        totalregion_gdf = region_gdf
        totalregion_gdf['raw_region_length'] = totalregion_gdf.lengths

        plot_lengths = totalregion_gdf['raw_region_length'] * totalregion_gdf['weight']
        plot_lengths = (self.xmax - self.xmin) * (plot_lengths / plot_lengths.sum())
        totalregion_gdf['plot_region_length'] = plot_lengths

        cumsum = totalregion_gdf['plot_region_length'].cumsum()
        plot_end0s = cumsum + self.xmin
        plot_start0s = np.insert(plot_end0s[:-1], 0, self.xmin)
        totalregion_gdf['plot_region_end0'] = plot_end0s
        totalregion_gdf['plot_region_start0'] = plot_start0s

        # set dfs without gap regions
        gap_selector = np.char.startswith(totalregion_gdf.chroms.astype(str), GAPREGION_PREFIX)
        totalregion_gdf_wogap = totalregion_gdf.loc[~gap_selector, :]
        gapregion_gdf = totalregion_gdf.loc[gap_selector, :]

        # set chromosome-wise params
        chromwise_params = dict()
        for chrom, subgdf in totalregion_gdf.group_bychrom(sort=False).items():
            chromwise_params[chrom] = {
                'start0': subgdf.start0s,
                'end0': subgdf.end0s,
                'raw_region_length': subgdf['raw_region_length'],
                'plot_region_length': subgdf['plot_region_length'],
                'plot_region_start0': subgdf['plot_region_start0'],
            }

        # result
        self.totalregion_gdf = totalregion_gdf
        self.totalregion_gdf_wogap = totalregion_gdf_wogap
        self.gapregion_gdf = gapregion_gdf
        self.chromwise_params = chromwise_params

    ###################
    # utility methods #
    ###################

    def __eq__(self, other):
        return all(
            (
                refgenome.compare_refvers(self.refver, other.refver),
                self.xmin == other.xmin,
                self.xmax == other.xmax,
                self.totalregion_gdf == other.totalregion_gdf,
            )
        )

    def iter_totalregion_gdf(self, merge_same_chroms=True):
        if merge_same_chroms:
            chrom_idxs = self.totalregion_gdf_wogap.chromosome_indexes
            grouper = np.cumsum(np.concatenate([[0], np.diff(chrom_idxs)]))
        else:
            grouper = np.arange(self.totalregion_gdf_wogap.nrow)
        return iter(self.totalregion_gdf_wogap.df.groupby(grouper))

    @property
    def xlim(self):
        return (self.xmin, self.xmax)
        #start0 = self.totalregion_gdf['plot_region_start0'][0]
        #end0 = self.totalregion_gdf['plot_region_end0'][-1] - 1
        #return (start0, end0)

    @staticmethod
    @deco.get_deco_atleast1d(['pos0_list'])
    def genomic_to_plot(chromwise_params, chrom, pos0_list):
        """
        Example:
            Interval (0, 5)
            |-----|-----|-----|-----|-----|
            0     1     2     3     4     5

            coord of data to draw = 2
                        <--> : common offset
            <-----------> : data offset
            |-----|-----|--*--|-----|-----|
            0     1     2     3     4     5

            coord of data to draw = 0
            <--> : common offset
            |--*--|-----|-----|-----|-----|
            0     1     2     3     4     5
        """

        if chrom not in chromwise_params.keys():
            raise Exception(f'Input "chrom" argument is not included in the plotting region.')

        pos0_list = pos0_list[:, np.newaxis]
        params = chromwise_params[chrom]

        contains = np.logical_and(
            (pos0_list >= params['start0']), (pos0_list <= params['end0'])
        )
        idxs = np.stack(np.nonzero(contains), axis=1)
        sorted_idxs = idxs[np.argsort(idxs[:, 0])]
        selector = np.insert(np.diff(sorted_idxs[:, 0]), 0, 1).nonzero()
        dedup_idxs = sorted_idxs[selector]
        pos0s_indexes = dedup_idxs[:, 0]
        intv_indexes = dedup_idxs[:, 1]
        #pos0s_indexes, intv_indexes = np.where(contains)
            # np.ndarray composed of the indexes of the containing intervals
            # identical intervals can appear many times

        selected_plotregion_lengths = params['plot_region_length'][intv_indexes]
        selected_rawregion_lengths = params['raw_region_length'][intv_indexes]  # this is integer
        #common_offsets = 0.5 * (selected_plotregion_lengths / selected_rawregion_lengths)
        data_offsets = (
            selected_plotregion_lengths
            * (
                (pos0_list[pos0s_indexes, 0] - params['start0'][intv_indexes]) 
                / selected_rawregion_lengths
            )
        )
        #within_region_offsets = data_offsets + common_offsets
        within_region_offsets = data_offsets
        return params['plot_region_start0'][intv_indexes] + within_region_offsets

    #def genomic_to_plot(self, chrom, pos0_list):
    #    return genomic_to_plot(self.chromwise_params, chrom, pos0_list)

    #def genomic_to_plot_with_indexes(self, chrom, pos0_list, indexes):
    #    plot_coords = self.genomic_to_plot(self.chromwise_params, chrom, pos0_list)
    #    return (indexes, plot_coords)

    def plot_to_genomic(self, plotcoord_list, return_modified=False):
        """
        Input plot coordinate may lie within a gap region

        Example:
            Interval (0, 5)
            |-----|-----|-----|-----|-----|
            0     1     2     3     4     5

            Input plot coord is here => treated as "2"
                            |
            |-----|-----|---*-|-----|-----|
            0     1     2     3     4     5
        """
        plotcoord_list = np.atleast_1d(plotcoord_list)
        assert np.logical_and(
            plotcoord_list >= self.xlim[0],
            plotcoord_list <= self.xlim[1],
        ).all(), f'Input plot coordinates are out of plot limits'

        bins = np.append(
            self.totalregion_gdf['plot_region_start0'], 
            self.totalregion_gdf['plot_region_end0'][-1],
        )
        plotregion_indexes = np.digitize(plotcoord_list, bins, right=False)
        #assert not np.isin(plotregion_indexes, [0, len(bins)]).any()
        assert not (plotregion_indexes == 0).any()

        # shift indexes of plotcoord on the right margin
        plotregion_indexes[plotregion_indexes == len(bins)] -= 1

        # shift indexes of plotcoord on the left margin of gap regions
        for gapregion_start in self.gapregion_gdf['plot_region_start0']:
            plotregion_indexes[plotcoord_list == gapregion_start] -= 1

        # final shift
        plotregion_indexes -= 1

        # results
        subgdf = self.totalregion_gdf.iloc[plotregion_indexes, :]
        result_chroms = subgdf.chroms
        within_region_frac = (
            (plotcoord_list - subgdf['plot_region_start0'])
            / subgdf['plot_region_length']
        )
        offsets = np.floor(within_region_frac * subgdf.lengths).astype(int)
        result_pos0s = subgdf.start0s + offsets

        if return_modified:
            # move to nearest integer positions
            #wogap_regions = self.totalregion_gdf_wogap
            #region_genomic_lengths = wogap_regions.end0s - wogap_regions.start0s
            #region_plot_lengths = wogap_regions['plot_region_end0'] - wogap_regions['plot_region_start0']
            #onebase_lengths = region_plot_lengths / region_genomic_lengths
            onebase_lengths = (
                (subgdf['plot_region_end0'] - subgdf['plot_region_start0'])
                / subgdf.lengths
            )
            offset_ints = np.floor(
                (plotcoord_list - subgdf['plot_region_start0'])
                / onebase_lengths
            )
            modified_plotcoords = subgdf['plot_region_start0'] + (offset_ints * onebase_lengths)

            return result_chroms, result_pos0s, modified_plotcoords
        else:
            return result_chroms, result_pos0s

#    def plot_to_genomic_old(self, x):
#        plot_intvlist = self.data.index.get_level_values('plot_interval')
#        genome_intvlist = self.data.index.get_level_values('genome_interval')
#        contains = plot_intvlist.contains(x)
#
#        num_hit = contains.sum()
#        if num_hit == 0:
#            return None
#        elif num_hit > 1:
#            raise Exception(f'More than one intervals contains the input position.')
#
#        idx = np.where(contains)[0][0]
#
#        chrom = self.data.index.get_level_values('chromosome')[idx]
#
#        plot_intv = plot_intvlist[idx]
#        genome_intv = genome_intvlist[idx]
#        regional_offset_fraction = (x - plot_intv.left) / plot_intv.length
#        pos0 = int(np.rint(genome_intv.left + (genome_intv.length * regional_offset_fraction)))
#
#        return chrom, pos0

    ##################################################
    # pos1-style genomic_to_plot and plot_to_genomic #
    ##################################################

    @staticmethod
    @deco.get_deco_atleast1d(['pos0_list'])
    def genomic_to_plot_pos1style(chromwise_params, chrom, pos0_list):
        """
        Example:
            Interval (0, 5)
            |-----|-----|-----|-----|-----|
               0     1     2     3     4   

            coord of data to draw = 2
                        <--> : common offset
            <-----------> : data offset
            |-----|-----|--*--|-----|-----|
               0     1     2     3     4   

            coord of data to draw = 0
            <--> : common offset
            |--*--|-----|-----|-----|-----|
               0     1     2     3     4   
        """

        if chrom not in chromwise_params.keys():
            raise Exception(f'Input "chrom" argument is not included in the plotting region.')

        pos0_list = pos0_list[:, np.newaxis]
        params = chromwise_params[chrom]

        contains = np.logical_and(
            (pos0_list >= params['start0']), (pos0_list <= params['end0'])
        )
        idxs = np.stack(np.nonzero(contains), axis=1)
        sorted_idxs = idxs[np.argsort(idxs[:, 0])]
        selector = np.insert(np.diff(sorted_idxs[:, 0]), 0, 1).nonzero()
        dedup_idxs = sorted_idxs[selector]
        pos0s_indexes = dedup_idxs[:, 0]
        intv_indexes = dedup_idxs[:, 1]
        #pos0s_indexes, intv_indexes = np.where(contains)
            # np.ndarray composed of the indexes of the containing intervals
            # identical intervals can appear many times

        selected_plotregion_lengths = params['plot_region_length'][intv_indexes]
        selected_rawregion_lengths = params['raw_region_length'][intv_indexes]  # this is integer
        onebase_lengths = selected_plotregion_lengths / selected_rawregion_lengths
        common_offsets = 0.5 * onebase_lengths
        data_offsets = (
            selected_plotregion_lengths
            * (
                (pos0_list[pos0s_indexes, 0] - params['start0'][intv_indexes]) 
                / selected_rawregion_lengths
            )
        )
        within_region_offsets = data_offsets + common_offsets
        return params['plot_region_start0'][intv_indexes] + within_region_offsets

    def plot_to_genomic_pos1style(self, plotcoord_list, return_modified=False):
        """
        Input plot coordinate may lie within a gap region

        Example:
            Interval (0, 5)
            |-----|-----|-----|-----|-----|
               0     1     2     3     4   

            Input plot coord is here => treated as "2"
                            |
            |-----|-----|-----|-----|-----|
               0     1     2     3     4   
        """
        plotcoord_list = np.atleast_1d(plotcoord_list)
        assert np.logical_and(
            plotcoord_list >= self.xlim[0],
            plotcoord_list < self.xlim[1],
        ).all(), f'Input plot coordinates are out of plot limits'

        bins = np.append(
            self.totalregion_gdf['plot_region_start0'], 
            self.totalregion_gdf['plot_region_end0'][-1],
        )
        plotregion_indexes = np.digitize(plotcoord_list, bins, right=False)
        assert not np.isin(plotregion_indexes, [0, len(bins)]).any()
        plotregion_indexes = plotregion_indexes - 1

        # results
        subgdf = self.totalregion_gdf.iloc[plotregion_indexes, :]
        result_chroms = subgdf.chroms
        within_region_frac = (
            (plotcoord_list - subgdf['plot_region_start0'])
            / subgdf['plot_region_length']
        )
        offsets = np.floor(within_region_frac * subgdf.lengths).astype(int)
        result_pos0s = subgdf.start0s + offsets

        if return_modified:
            # move to nearest integer positions
            wogap_regions = self.totalregion_gdf_wogap
            region_genomic_lengths = wogap_regions.end0s - wogap_regions.start0s
            region_plot_lengths = wogap_regions['plot_region_end0'] - wogap_regions['plot_region_start0']
            onebase_lengths = region_plot_lengths / region_genomic_lengths

            offset_ints = np.floor(
                (plotcoord_list - wogap_regions['plot_region_start0'][plotregion_indexes])
                / onebase_lengths[plotregion_indexes]
            )
            modified_plotcoords = (
                wogap_regions['plot_region_start0'][plotregion_indexes] 
                + ((offset_ints + 0.5) * onebase_lengths[plotregion_indexes])
            )

            return result_chroms, result_pos0s, modified_plotcoords
        else:
            return result_chroms, result_pos0s

    # Axes modification
    def get_chrom_borders(self, merge_same_chroms=True):
        """Chromosome names starting with "-1" are omitted"""
        result = list()
        for key, subdf in self.iter_totalregion_gdf(merge_same_chroms=merge_same_chroms):
            chroms = set(subdf['Chromosome'])
            assert len(chroms) == 1
            chrom = chroms.pop()
            assert not chrom.startswith(GAPREGION_PREFIX)

            result.append(
                (
                    chrom, 
                    subdf['plot_region_start0'].iloc[0], 
                    subdf['plot_region_end0'].iloc[-1],
                )
            )
        return result

    #######################
    # plotdata generation #
    #######################

#    def make_plotdata_old(self, data, log_suffix=None):
#        if log_suffix is not None:
#            logutils.log(f'Beginning plotdata generation{log_suffix}')
#
#        #assert isinstance(data, GDF)
#        isec_gdf = data.intersect(self.totalregion_gdf_wogap)
#        if isec_gdf.is_empty:
#            return False
#
#        isec_gdf.sort()
#
#        result_start0s = list()
#        result_end0s = list()
#        ordered_chroms = tools.unique_keeporder(isec_gdf.chromosomes)
#
#        for chrom in ordered_chroms:
#            subgdf = isec_gdf.subset_chroms(chrom)
#            result_start0s.extend(
#                genomic_to_plot(self.chromwise_params, chrom, subgdf.starts)
#            )
#            result_end0s.extend(
#                genomic_to_plot(self.chromwise_params, chrom, subgdf.ends - 1) + 1
#            )
#
#        isec_gdf['plot_start0s'] = result_start0s
#        isec_gdf['plot_end0s'] = result_end0s
#
#        if log_suffix is not None:
#            logutils.log(f'Finished plotdata generation{log_suffix}')
#
#        return isec_gdf

    #@staticmethod
    @classmethod
    def make_plotdata_targetfunc(cls, partial_isec_gdf, chromwise_params):
        chrom = partial_isec_gdf.chroms[0]
        plot_start0s = cls.genomic_to_plot(chromwise_params, chrom, partial_isec_gdf.start0s)
        #plot_end0s = cls.genomic_to_plot(chromwise_params, chrom, partial_isec_gdf.end0s - 1) + 1
        plot_end0s = cls.genomic_to_plot(chromwise_params, chrom, partial_isec_gdf.end0s)
        return plot_start0s, plot_end0s

    def make_plotdata(self, data, log_suffix='', nproc=1, split_width=10000, verbose=True):
        if verbose:
            logutils.log(f'Beginning plotdata generation{log_suffix} (nproc={nproc})')

        #assert isinstance(data, GDF)
        if verbose:
            logutils.log(f'Beginning intersection')
        isec_gdf = data.intersect(self.totalregion_gdf_wogap, nproc=nproc)
        if verbose:
            logutils.log(f'Finished intersection')

        if isec_gdf.is_empty:
            return False

        if verbose:
            logutils.log(f'Beginning plot coordinate calculation')

        isec_gdf.sort()
        isec_gdf_splits = isec_gdf.equal_nrow_split_keepchrom(width=split_width)

        args = ((gdf, self.chromwise_params) for gdf in isec_gdf_splits)
        with multiprocessing.Pool(nproc) as pool:
            mp_result = pool.starmap(self.make_plotdata_targetfunc, args)

        plot_start0s, plot_end0s = zip(*mp_result)
        plot_start0s = np.concatenate(plot_start0s)
        plot_end0s = np.concatenate(plot_end0s)

        isec_gdf['plot_start0s'] = plot_start0s
        isec_gdf['plot_end0s'] = plot_end0s

        if verbose:
            logutils.log(f'Finished plot coordinate calculation')

        if verbose:
            logutils.log(f'Finished plotdata generation{log_suffix} (nproc={nproc})')

        return isec_gdf


class GenomePlotter:
    def __init__(self, refver, **kwargs):
        self.refver = refver
        self.cconv = CoordinateConverter(refver, **kwargs)

    def make_plotdata(self, data, log_suffix='', nproc=1, split_width=1000, verbose=True):
        return self.cconv.make_plotdata(
            data, 
            log_suffix=log_suffix, 
            nproc=nproc, 
            split_width=split_width,
            verbose=verbose,
        )

    def genomic_to_plot(self, *args, **kwargs):
        return self.cconv.genomic_to_plot(*args, **kwargs)

    def plot_to_genomic(self, *args, **kwargs):
        return self.cconv.plot_to_genomic(*args, **kwargs)

    @property
    def region_gdf(self):
        return self.cconv.totalregion_gdf_wogap

    def get_intersect(self, gdf):
        return gdf.intersect(self.region_gdf)

    #############################
    # common drawer & decorator #
    #############################

    def draw_common(
        self, 
        ax, 

        n_xlabel=None,

        split_spines=True,

        merge_same_chroms=True,

        chromlabel_kwargs=dict(), 
        draw_chromlabel=True,

        title=None,
        title_kwargs=dict(),
    ):
        """Should be done after data drawings are finished"""

        draw_common_artists = list()

        self.set_xlim(ax)

        # horizontal lines at the level of yticks
        ylims = ax.get_ylim()
        yticks = [
            x for x in ax.get_yticks()
            if (x > ylims[0]) and (x < ylims[1])
        ]
        linelist = self.draw_grids(
            ax, 
            ys=yticks, 
            line_params=dict(), 
            merge_same_chroms=merge_same_chroms,
        )
        draw_common_artists.extend(linelist)

        # xaxis label (genomic coordinates)
        if n_xlabel is not None:
            self.draw_genomecoord_labels(ax, n=n_xlabel)
        else:
            ax.set_xticks([])

        # centromere bgcolor
        ylims = ax.get_ylim()
        patchcol = self.draw_centromeres(ax, ymins=ylims[0], ymaxs=ylims[1])
        draw_common_artists.append(patchcol)
        #self.draw_centromeres_type2(ax)

        # spine modification
        if split_spines:
            subartists = self.fit_spines_to_regions(
                ax,
                ylims=ylims,
                draw_chromlabel=draw_chromlabel,
                merge_same_chroms=merge_same_chroms,
                chromlabel_kwargs=chromlabel_kwargs,
                title=title,
                title_kwargs=title_kwargs,
            )
            draw_common_artists.extend(subartists)

        return draw_common_artists

    def draw_decorator(func):
        sig = inspect.signature(func)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()

            # data, plotdata
            if set(['plotdata', 'data']).issubset(ba.arguments.keys()):
                assert set(['log_suffix', 'nproc', 'verbose']).issubset(ba.arguments.keys())
                if ba.arguments['plotdata'] is None:

                    #logutils.log(ba.arguments['verbose'], verbose_locstring=True)

                    gplotter = ba.arguments['self']
                    ba.arguments['plotdata'] = gplotter.make_plotdata(
                        ba.arguments['data'],
                        log_suffix=ba.arguments['log_suffix'], 
                        nproc=ba.arguments['nproc'], 
                        verbose=ba.arguments['verbose'],
                    )
                    del ba.arguments['data']

                if ba.arguments['plotdata'] is False:
                    nodraw = True
                else:
                    nodraw = False
            else:
                nodraw = True

            # ys
            if set(['ys', 'y_colname']).issubset(ba.arguments.keys()):
                if (ba.arguments['y_colname'] is not None) and (not nodraw):
                    ba.arguments['ys'] = ba.arguments['plotdata'][ba.arguments['y_colname']]
                    del ba.arguments['y_colname']

            # main drawing
            if nodraw:
                result = None
            else:
                result = func(**ba.arguments)

            # draw_common
            draw_common_artists = None
            if 'draw_common' in ba.arguments:
                if ba.arguments['draw_common']:
                    draw_common_artists = ba.arguments['self'].draw_common(
                        ba.arguments['ax'], 
                        **ba.arguments['draw_common_kwargs'], 
                    )

            return result, draw_common_artists

        return wrapper

    ###############################
    # components of common drawer #
    ###############################

    def set_xlim(self, ax, margin=0.01):
        #ax.set_xlim(*self.cconv.xlim)
        width = self.cconv.xlim[1] - self.cconv.xlim[0]
        margin = width * margin
        ax.set_xlim(self.cconv.xlim[0] - margin, self.cconv.xlim[1] + margin)

    def draw_genomecoord_labels(self, ax, n=10, pos1=False):
        """Should be done after data drawings are finished"""
        #xlim = self.cconv.xlim
        #plotcoords = np.linspace(xlim[0], xlim[1], num=n, endpoint=True)

        # generate labeling positions
        raw_plotcoords = tools.gapped_linspace(
            self.region_gdf['plot_region_start0'], 
            self.region_gdf['plot_region_end0'], 
            num=n,
            return_indexes=False,
            endpoint=True,
        )

        chroms, pos0s, plotcoords = self.cconv.plot_to_genomic(raw_plotcoords, return_modified=True)
        plotcoords, uniq_idxs = np.unique(plotcoords, return_index=True)
        chroms = chroms[uniq_idxs]
        pos0s = pos0s[uniq_idxs]

        #chroms = [
        #    x if x.startswith('chr') else ('chr' + x)
        #    for x in chroms
        #]
        if pos1:
            pos1s = pos0s + 1
            #pos1_strings = tools.shorten_int(pos1s)
            pos_strings = [f'{x:,}' for x in pos1s]
        else:
            pos_strings = [f'{x:,}' for x in pos0s]

        labels = [f'{x} : {y}' for x, y in zip(chroms, pos_strings)]
        ax.set_xticks(plotcoords, labels=labels, minor=False, rotation=90)

    def fit_spines_to_regions(
        self, ax, ylims,
        draw_chromlabel=True,
        prefix_with_chr=False,
        chromlabel_offset=0.01,
        chromlabel_kwargs=dict(), 
        line_kwargs=dict(),
        merge_same_chroms=True,
        title=None,
        title_kwargs=dict(),
    ):
        """What it does:
            1) draw spines (axes border lines) (top, bottom, left, right, region borders)
            2) draw region names above top spine (e.g. chr1, chr3:100-200)

        Should be done after data drawings are finished

        Args:
            chromlabel_offset: 
        """
        artists = list()

        # set plotting kwargs
        chromlabel_kwargs = (
            dict(ha='center', va='bottom', size=8)
            | chromlabel_kwargs
        )
        line_kwargs = (
            dict(color='black', linewidth=1)
            | line_kwargs
        )

        # get region border coordinates
        chrom_borders = self.cconv.get_chrom_borders(
            merge_same_chroms=merge_same_chroms,
        )

        # top and bottom spines
        ax.spines[['top', 'bottom']].set_visible(False)
        _, start0s, end0s = zip(*chrom_borders)
        linecol = ax.hlines(
            *np.broadcast_arrays(ylims[1], start0s, end0s),
            #np.repeat(ylims[1], len(start0s)), start0s, end0s, 
            color='black', linewidth=1,
        )
        artists.append(linecol)
        linecol = ax.hlines(
            *np.broadcast_arrays(ylims[0], start0s, end0s),
            #np.repeat(ylims[0], len(start0s)), start0s, end0s, 
            color='black', linewidth=1.5,
        )
        artists.append(linecol)

        # vertical spines - left and right margins
        #ax.spines[['left', 'right']].set_visible(False)
        #ax.vlines(ax.get_xlim(), ylims[0], ylims[1], **line_kwargs)

        # vertical spines - region borderlines
        ax.spines[['left', 'right']].set_visible(False)
        border_pos0s = set()
        #xlim = self.cconv.xlim
        for _, start0, end0 in chrom_borders:
            #if start0 != xlim[0]:
            border_pos0s.add(start0)
            #if end0 != xlim[1]:
            border_pos0s.add(end0)
        #ax.vlines(tuple(border_pos0s), ylims[0], ylims[1], **line_kwargs)
        linecol = ax.vlines(*np.broadcast_arrays(tuple(border_pos0s), ylims[0], ylims[1]), **line_kwargs)
        artists.append(linecol)

        # chromosome name texts
        if draw_chromlabel:
            chromlabel_y = ylims[1] + chromlabel_offset * (ylims[1] - ylims[0])
            for chrom, start0, end0 in chrom_borders:
                if chrom.startswith(GAPREGION_PREFIX):
                    continue
                if prefix_with_chr:
                    if not chrom.startswith('chr'):
                        chrom = 'chr' + chrom
                txt = ax.text(
                    (start0 + end0) / 2, 
                    chromlabel_y, 
                    chrom, 
                    **chromlabel_kwargs,
                )
                artists.append(txt)

        # axes title
        if title is not None:
            kwargs = (
                {'weight': 'bold', 'size': 20, 'y': 1.1} 
                | title_kwargs
            )
            ax.set_title(title, **kwargs)
        else:
            ax.set_title(None)

        return artists

    def draw_grids(
        self, 
        ax, 
        ys, 
        line_params=dict(), 
        merge_same_chroms=True,
    ):
        line_params = (
            dict(color='black', linewidth=0.2, alpha=0.5)
            | line_params
        )
        chroms, start0s, end0s = zip(
            *self.cconv.get_chrom_borders(
                merge_same_chroms=merge_same_chroms,
            )
        )

        linelist = list()
        for y in ys:
            linecol = ax.hlines(
                np.repeat(y, len(start0s)), start0s, end0s, 
                **line_params,
            )
            linelist.append(linecol)
        return linelist

    #############################
    # elemental drawing methods #
    #############################

    @staticmethod
    def get_point_plot_coord(plot_start0s, plot_end0s):
        return (plot_start0s + plot_end0s) / 2

    @draw_decorator
    @deco.get_deco_num_set_differently(('ys', 'y_colname'), 1)
    def draw_hlines(
        self, 
        ax, 
        *, 

        # plotdata generation
        data=None, 
        plotdata=None, 
        nproc=1,
        log_suffix='',
        verbose=True,

        ys=None,
        y_colname=None, 

        offset=None,
        plot_kwargs=dict(),

        draw_common=True, 
        draw_common_kwargs=dict(),
    ):
        plot_kwargs = (
            dict()
            | plot_kwargs
        )
        # draw data
        if offset is not None:
            ys = ys + offset

        if plotdata.nrow > 1:
            ys, xmins, xmaxs = self._merge_adjacent_data_new(
                genome_xmins=plotdata.start0s, 
                genome_xmaxs=plotdata.end0s, 
                plot_xmins=plotdata['plot_start0s'], 
                plot_xmaxs=plotdata['plot_end0s'], 
                ys=ys,
            )
        else:
            xmins = plotdata['plot_start0s'][0]
            xmaxs = plotdata['plot_end0s'][0]

        linecol = ax.hlines(ys, xmins, xmaxs, **plot_kwargs)
        return linecol

    @draw_decorator
    @deco.get_deco_num_set_differently(('ys', 'y_colname'), 1)
    def draw_dots(
        self, ax, 
        *, 

        # plotdata generation
        data=None, 
        plotdata=None,
        nproc=1,
        log_suffix='',
        verbose=True,

        ys=None,
        y_colname=None, 

        plot_kwargs=dict(),

        draw_common=True, 
        draw_common_kwargs=dict(),

        label=None,
    ):
        # kwargs handling
        plot_kwargs = (
            {
                'color': 'black',
                'marker': 'o',
                'linestyle': '',
            } | plot_kwargs
        )

        #xs = (plotdata['plot_start0s'] + (plotdata['plot_end0s'] - 1)) / 2
        xs = self.get_point_plot_coord(plotdata['plot_start0s'], plotdata['plot_end0s'])
        line2d, = ax.plot(xs, ys, label=label, **plot_kwargs)
        return line2d

    @draw_decorator
    @deco.get_deco_num_set_differently(('ys', 'y_colname'), 1)
    @deco.get_deco_num_set_differently(('color_colname', 'color_values'), 2, 'lt')
    def draw_dots_scatter(
        self, ax,
        *, 

        # plotdata generation
        data=None, 
        plotdata=None,
        nproc=1,
        log_suffix='',
        verbose=True,

        ys=None,
        y_colname=None, 

        color_colname=None,
        color_values=None,
        plot_kwargs=dict(),

        draw_common=True, 
        draw_common_kwargs=dict(),
    ):
        # kwargs handling
        plot_kwargs = (
            {
                'c': 'black',
                'marker': 'o',
                'markersize': 0.3, 
            } | plot_kwargs
        )

        # main
        #xs = (plotdata['plot_start0s'] + (plotdata['plot_end0s'] - 1)) / 2
        xs = self.get_point_plot_coord(plotdata['plot_start0s'], plotdata['plot_end0s'])
        if color_colname is not None:
            plot_kwargs['c'] = plotdata[color_colname]
        elif color_values is not None:
            plot_kwargs['c'] = color_values

        if 'color' in plot_kwargs:
            del plot_kwargs['color']
        if 'markersize' in plot_kwargs:
            plot_kwargs['s'] = plot_kwargs['markersize']
            del plot_kwargs['markersize']

        pathcol = ax.scatter(xs, ys, **plot_kwargs)
        return pathcol

    @draw_decorator
    @deco.get_deco_num_set_differently(('ys', 'y_colname'), 1)
    @deco.get_deco_num_set_differently(('bottom_values', 'bottom_colname'), 2, 'lt')
    def draw_bars(
        self, ax, 
        *, 

        # plotdata generation
        data=None, 
        plotdata=None,
        nproc=1,
        log_suffix='',
        verbose=True,

        ys=None,
        y_colname=None, 

        bottom_values=None,
        bottom_colname=None,

        plot_kwargs=dict(),

        draw_common=True, 
        draw_common_kwargs=dict(),
    ):
        plot_kwargs = (
            dict(alpha=0.5, color='tab:blue')
            | plot_kwargs
        )

        xs = (plotdata['plot_end0s'] + plotdata['plot_start0s']) / 2
        widths = plotdata['plot_end0s'] - plotdata['plot_start0s']

        if bottom_values is not None:
            bottoms = bottom_values
        elif bottom_colname is not None:
            bottoms = plotdata[bottom_colname]
        else:
            bottoms = 0 

        bars = ax.bar(xs, height=ys, width=widths, bottom=bottoms, **plot_kwargs)
        return bars

    @draw_decorator
    @deco.get_deco_num_set_differently(('texts', 'text_colname'), 1)
    def draw_texts(
        self, ax, 
        *,

        # plotdata generation
        data=None, 
        plotdata=None,
        nproc=1,
        log_suffix='',
        verbose=True,

        ys=None,
        y_colname=None, 

        texts=None,
        text_colname=None,

        text_kwargs=dict(),

        draw_common=True, 
        draw_common_kwargs=dict(),
    ):
        default_text_kwargs = dict(size=8)
        text_kwargs = default_text_kwargs | text_kwargs

        # prepare text positions
        #xs = (plotdata['plot_start0s'] + (plotdata['plot_end0s'] - 1)) / 2
        xs = self.get_point_plot_coord(plotdata['plot_start0s'], plotdata['plot_end0s'])
        if ys is None:
            ys = 1.05
            transform = ax.transAxes
        else:
            transform = ax.transData

        # prepare text values
        if texts is None:
            texts = plotdata[text_colname].to_numpy()
        else:
            texts = np.atleast_1d(texts)

        xs, ys, texts = np.broadcast_arrays(xs, ys, texts)

        # draw
        result = list()
        for t, x, y in zip(texts, xs, ys):
            result.append(
                ax.text(x, y, t, **text_kwargs)
            )
        return result

    @draw_decorator
    def draw_boxes(
        self, ax, 
        *,

        # plotdata generation
        data=None, 
        plotdata=None,
        nproc=1,
        log_suffix='',
        verbose=True,

        ymins=None,
        ymaxs=None,
        colors=None,

        #plot_kwargs=dict(),
        rect_kwargs=dict(),

        draw_common=True, 
        draw_common_kwargs=dict(),
    ):
        # x
        if plotdata.nrow > 1:
            _, xmins, xmaxs = self._merge_adjacent_data_new(
                genome_xmins=plotdata.start0s, 
                genome_xmaxs=plotdata.end0s, 
                plot_xmins=plotdata['plot_start0s'], 
                plot_xmaxs=plotdata['plot_end0s'], 
                ys=None,
            )
        else:
            xmins = np.atleast_1d(plotdata['plot_start0s'][0])
            xmaxs = np.atleast_1d(plotdata['plot_end0s'][0])

        widths = xmaxs - xmins

        # y
        ylims = ax.get_ylim()
        if ymins is None:
            ymins = np.repeat(ylims[0], len(widths))
        else:
            ymins = np.broadcast_to(ymins, widths.shape)

        if ymaxs is None:
            ymaxs = np.repeat(ylims[1], len(widths))
        else:
            ymaxs = np.broadcast_to(ymaxs, widths.shape)

        heights = ymaxs - ymins

        # color
        if colors is None:
            colors = 'yellow'
        colors = np.broadcast_to(colors, widths.shape)

        # setup rect_kwargs
        default_rect_kwargs = {'alpha': 0.1, 'zorder': 0}

        # final
        rects = list()
        for (xm, w, ym, h, col) in zip(xmins, widths, ymins, heights, colors):
            this_rect_kwargs = (
                default_rect_kwargs
                | {'facecolor': col, 'edgecolor': col} 
                | rect_kwargs
            )
            rects.append(
                Rectangle((xm, ym), width=w, height=h, **this_rect_kwargs)
            )

        patchcol = PatchCollection(rects, match_original=True)
        ax.add_collection(patchcol)
        return patchcol

#        # color
#        if np.isscalar(colors):
#            match_original = False
#            default_plot_kwargs['color'] = colors
#        else:
#            match_original = True
#
#        # final
#        boxes = [
#            Rectangle((xm, ym), width=w, height=h)
#            for (xm, w, ym, h) in zip(xmins, widths, ymins, heights)
#        ]
#        if match_original:
#            for col, box in zip(colors, boxes):
#                box.set(color=col)
#
#        ax.add_collection(
#            PatchCollection(
#                boxes, 
#                match_original=match_original, 
#                **default_plot_kwargs,
#            )
#        )


    #####################################
    # combinations of elemental drawers #
    #####################################

    @draw_decorator
    def draw_features(
        self, ax,
        *,

        # plotdata generation
        data=None,
        plotdata=None,
        nproc=1,
        log_suffix='',
        verbose=True,

        y_features=None,
        y_labels=None,
        #feature_as_dot=False,
        draw_label=True,

        text_kwargs=dict(),
        line_kwargs=dict(),

        draw_common=True, 
        draw_common_kwargs=dict(),
    ):
        # setup plot_kwargs
        text_kwargs = (
            dict()
            | text_kwargs
        )
        line_kwargs = (
            dict(linewidth=2, color='black')
            | line_kwargs
        )

        # main
        ylims = ax.get_ylim()
        if y_features is None:
            y_features = ylims[0] + 0.66 * (ylims[1] - ylims[0])
        if y_labels is None:
            y_labels = ylims[0] + 0.33 * (ylims[1] - ylims[0])

        # draw hline
#        if feature_as_dot:
#            self.draw_dots(
#                ax, 
#                y_val=y_features,
#                df_plotdata=df_plotdata,
#            )
#        else:
        self.draw_hlines(
            ax, 
            ys=y_features,
            plotdata=plotdata, 
            offset=0,
            plot_kwargs=line_kwargs,
            draw_common=False, 
        )

        # draw texts
        if draw_label:
            assert 'name' in plotdata.columns
            for name, subdf in plotdata.df.groupby('name'):
                subgdf = plotdata.spawn(subdf)
                subgdf.sort()
                #subdf = cnvmisc.sort_genome_df(subdf, refver=self.refver)
                if subgdf.nrow == 1:
                    xmins = subgdf['plot_start0s']
                    xmaxs = subgdf['plot_end0s']
                else:
                    _, xmins, xmaxs = self._merge_adjacent_data_new(
                        genome_xmins=subgdf.start0s, 
                        genome_xmaxs=subgdf.end0s, 
                        plot_xmins=subgdf['plot_start0s'], 
                        plot_xmaxs=subgdf['plot_end0s'], 
                        ys=None,
                    )
                text_xlist = 0.5 * (xmins + xmaxs - 1)
                for text_x in text_xlist:
                    ax.annotate(
                        name, xy=(text_x, y_features), xytext=(text_x, y_labels),
                        arrowprops=dict(arrowstyle='-'),
                        ha='center',
                        va='center',
                        **text_kwargs,
                    )

    def draw_ideogram(
        self, 
        ax,
    ):
        cytoband_gdf = ucscdata.get_cytoband(self.refver)
        colors = [ucscdata.CYTOBAND_COLORMAP[x] for x in cytoband_gdf['Stain']]
        self.draw_boxes(
            ax=ax, 
            data=cytoband_gdf, 
            rect_kwargs=dict(alpha=1),
            colors=colors,
            draw_common=False, 
            verbose=False,
        )

    def draw_centromeres(
        self, 
        ax, 
        ymins=None, 
        ymaxs=None,
        #draw_common=True, 
        #draw_common_kwargs=dict(),
    ):
        cytoband_gdf = ucscdata.get_cytoband(refver=self.refver)
        rects, _ = self.draw_boxes(
            ax, 
            data=cytoband_gdf.loc[cytoband_gdf['Stain'] == 'acen', :],
            #data=cytoband_gdf.loc[cytoband_gdf['Stain'].isin(['acen', 'gvar', 'stalk']), :],
            ymins=ymins,
            ymaxs=ymaxs,
            colors='red',
            rect_kwargs=dict(alpha=0.3, fill=None, hatch='//'),
            draw_common=False,
            verbose=False,
        )
        return rects

    @draw_decorator
    def draw_centromeres_type2(
        self, 
        ax,
        draw_common=True, 
        draw_common_kwargs=dict(),
    ):
        cytoband_gdf = ucscdata.get_cytoband(refver=self.refver)
        mapping = {
            'acen': 'red', 
            'gvar': 'green', 
            'stalk': 'blue',
        }

        def helper(bandname):
            self.draw_boxes(
                ax, 
                data=cytoband_gdf[cytoband_gdf['Stain'] == bandname, :], 
                colors=mapping[bandname],
                rect_kwargs=dict(alpha=0.3, linewidth=0),
                draw_common=False, 
                verbose=False,
            )

        helper('acen')
        helper('gvar')
        helper('stalk')

#    def prepare_plot_data_old(self, df, nproc=None):
#        # create isec between total region and input data
#        isec_gr, subgrs_bychrom = self._isec_trim_data_df(df)
#
#        # Add "End_minus1" columns; "End" columns cannot be used for plot coordinate calculation
#        for chrom, subgr in subgrs_bychrom.items():
#            subgr.End_minus1 = subgr.End - 1
#
#        xmins = self._get_ordered_plot_coords(subgrs_bychrom, 'Start', nproc=nproc)
#        xmaxs_minus1 = self._get_ordered_plot_coords(subgrs_bychrom, 'End_minus1', nproc=nproc)
#        xmaxs = xmaxs_minus1 + 1
#
#        return {
#            'isec_gr': isec_gr,
#            'subgrs_bychrom': subgrs_bychrom,
#            'xmins': xmins,
#            'xmaxs': xmaxs,
#        }

    ############################################

    @classmethod
    def _merge_adjacent_data_new(cls, genome_xmins, genome_xmaxs, plot_xmins, plot_xmaxs, ys=None):
        """Helper of draw_hlines"""
        if ys is None:
            ys = np.repeat(1, len(plot_xmins))

        ys = np.array(ys)
        plot_xmins = np.array(plot_xmins)
        plot_xmaxs = np.array(plot_xmaxs)
        genome_xmins = np.array(genome_xmins)
        genome_xmaxs = np.array(genome_xmaxs)

        flags = (ys[:-1] == ys[1:]) & (genome_xmaxs[:-1] == genome_xmins[1:])
        idx = 0
        indexes = list()
        if not flags[0]:
            indexes.append((0, 0))

        for key, subiter in itertools.groupby(flags):
            len_val = len(tuple(subiter))
            if key:
                start = idx
                end = idx + len_val
                indexes.append((start, end))
                idx += len_val
            else:
                indexes.extend((k, k) for k in range(idx + 1, idx + len_val))
                idx += len_val

        if not key:
            indexes.append((idx, idx))

        new_ys = ys[[x[0] for x in indexes]]
        new_plot_xmins = plot_xmins[[x[0] for x in indexes]]
        new_plot_xmaxs = plot_xmaxs[[x[1] for x in indexes]]

        return new_ys, new_plot_xmins, new_plot_xmaxs


#    @classmethod
#    def _merge_adjacent_data(cls, xmins, xmaxs, ys=None):
#        """Helper of draw_hlines"""
#        if ys is None:
#            ys = np.repeat(1, len(xmins))
#
#        ys = np.array(ys)
#        xmins = np.array(xmins)
#        xmaxs = np.array(xmaxs)
#
#        flags = (ys[:-1] == ys[1:]) & (xmaxs[:-1] == xmins[1:])
#        idx = 0
#        indexes = list()
#        if not flags[0]:
#            indexes.append((0, 0))
#
#        for key, subiter in itertools.groupby(flags):
#            len_val = len(tuple(subiter))
#            if key:
#                start = idx
#                end = idx + len_val
#                indexes.append((start, end))
#                idx += len_val
#            else:
#                indexes.extend((k, k) for k in range(idx + 1, idx + len_val))
#                idx += len_val
#
#        if not key:
#            indexes.append((idx, idx))
#
#        new_ys = ys[[x[0] for x in indexes]]
#        new_xmins = xmins[[x[0] for x in indexes]]
#        new_xmaxs = xmaxs[[x[1] for x in indexes]]
#
#        return new_ys, new_xmins, new_xmaxs

#    def _isec_trim_data_df(self, df, asis=False):
#        """helper of prepare_plot_data"""
#        assert '_index' not in df.columns
#
#        gr = cnvmisc.arg_into_gr(df)
#        if asis:
#            isec_gr = gr
#        else:
#            isec_gr = gr.intersect(self.cconv.totalregion_gr)
#
#        isec_gr._index = list(range(isec_gr.df.shape[0]))
#
#        subgrs_bychrom = dict()
#        for chrom in isec_gr.Chromosome.unique():
#            subgrs_bychrom[chrom] = isec_gr[chrom]
#
#        return isec_gr, subgrs_bychrom

#    def _get_ordered_plot_coords(self, subgrs_bychrom, pos0_colname, nproc=None):
#        """helper of prepare_plot_data"""
#        # Multiprocessing is slower than serial jobs!!
##        with multiprocessing.Pool(nproc) as pool:
##            result = pool.starmap(
##                self.genomic_to_plot_with_indexes, 
##                (
##                    (chrom, subdf[pos0_colname], subdf['_index']) 
##                    for chrom, subdf in subgrs_bychrom.items()
##                )
##            )
#
#        result = list()
#        for chrom, subgr in subgrs_bychrom.items():
#            result.append(
#                self.cconv.genomic_to_plot_with_indexes(
#                    chrom, getattr(subgr, pos0_colname), subgr._index
#                )
#            )
#
#        index_coord_pairs = itertools.chain.from_iterable(zip(*x) for x in result)
#        plot_coords = np.fromiter(
#            (x[1] for x in sorted(index_coord_pairs, key=operator.itemgetter(0))),
#            dtype=np.int_,
#        )
#
#        return plot_coords

#    def _get_ordered_plot_coords_woidx(self, subgrs_bychrom, pos0_colname):
#        """helper of prepare_plot_data"""
#        return np.fromiter(
#            itertools.chain.from_iterable(
#                self.cconv.genomic_to_plot(chrom, getattr(subgr, pos0_colname))
#                for chrom, subgr in subgrs_bychrom.items()
#            ),
#            dtype=int,
#        )


#class CNVPlotter:
#    def __init__(
#        self, 
#        refver=handygenome.OPTION['default_refver'], 
#        region_df=None, 
#        chroms=None, start0s=None, end0s=None, weights=None,
#        region_gaps=None,
#    ):
#        #super().__init__(refver=refver, region_df=region_df)
#        self.refver = refgenome.standardize_refver(refver)
#
#        region_df, region_gaps = handle_region_args(
#            refver, 
#            region_df=region_df, 
#            chroms=chroms, 
#            start0s=start0s, 
#            end0s=end0s, 
#            weights=weights,
#            region_gaps=region_gaps,
#        )
#
#        self.genomeplotter = GenomePlotter(
#            refver, region_df=region_df, region_gaps=region_gaps,
#        )
#        self.data = dict()
#        self.default_binsize = 100
#
#    def reset_genomeplotter(
#        self, 
#        region_df=None, 
#        chroms=None, start0s=None, end0s=None, weights=None,
#        region_gaps=None,
#    ):
#        region_df, region_gaps = handle_region_args(
#            refver=self.refver, 
#            region_df=region_df, 
#            chroms=chroms, 
#            start0s=start0s, 
#            end0s=end0s, 
#            weights=weights,
#            region_gaps=region_gaps,
#        )
#        self.genomeplotter = GenomePlotter(
#            self.refver, region_df=region_df, region_gaps=region_gaps,
#        )
#
#    def save_data(self, sampleid, outfile_path):
#        with open(outfile_path, 'wb') as outfile:
#            pickle.dump(self.data[sampleid], outfile)
#
#    def load_data(self, sampleid, infile_path):
#        with open(infile_path, 'rb') as infile:
#            self.data[sampleid] = pickle.load(infile)
#
#    ##############
#    # mainstream #
#    ##############
#
#    @deco.get_deco_num_set_differently(
#        ('normal_bam_path', 'normal_depth_path', 'normal_depth_df'), 2, 'lt',
#    )
#    @deco.get_deco_num_set_differently(
#        ('tumor_bam_path', 'tumor_depth_path', 'tumor_depth_df'), 2, 'lt',
#    )
#    @deco.get_deco_num_set_differently(
#        ('germline_vcf_path', 'germline_vafdf_path'), 2, 'lt',
#    )
#    def add_sample_file_new(
#        self, 
#        sampleid, 
#        is_female,
#
#        germline_vcf_path=None,
#        germline_vafdf_path=None,
#        vcf_sampleid_tumor=None,
#        vcf_sampleid_normal=None,
#
#        *,
#
#        norm_method='plain',
#
#        mode='wgs',
#        target_region=None,
#
#        normal_bam_path=None,
#        normal_depth_path=None, 
#        normal_depth_df=None, 
#
#        tumor_bam_path=None,
#        tumor_depth_path=None, 
#        tumor_depth_df=None, 
#
#        vcfload_nproc=1,
#        #mosdepth_postprocess_kwargs=dict(),
#
#        verbose=True,
#    ):
#        """Args:
#            *_depth_path: mosdepth output file
#            germline_vcf_path: germline variant vcf
#            vcf_sampleids: tuple of (normal sample id, tumor sample id)
#        """
#        def depth_loading_helper(
#            self, 
#            bam_path, 
#            depth_path, 
#            depth_df, 
#            sampleid, 
#            sampletype,
#        ):
#            assert sampletype in ('tumor', 'normal')
#
#            if bam_path is not None:
#                LOGGER_INFO.info(f'Loading {sampletype} depth - running mosdepth')
#                depth_df = libmosdepth.run_mosdepth(
#                    bam_path, 
#                    t=8, 
#                    use_median=False, 
#                    region_bed_path=None, 
#                    region_gr=(
#                        None
#                        if self.data[sampleid]['mode'] == 'wgs' else
#                        self.data[sampleid]['target_region']
#                    ), 
#                    window_size=self.default_binsize, 
#                    donot_subset_bam=True,
#                    as_gr=False, 
#                    load_perbase=False,
#                )
#            elif depth_path is not None:
#                LOGGER_INFO.info(f'Loading {sampletype} depth - reading mosdepth output file')
#                depth_df = libmosdepth.load_mosdepth_output(
#                    depth_path, depth_colname='mean_depth', as_gr=False,
#                )
#            elif depth_df is not None:
#                LOGGER_INFO.info(f'Loading {sampletype} depth - using the given depth dataframe')
#                assert isinstance(depth_df, pd.DataFrame)
#                assert set(depth_df.columns) == {'Chromosome', 'Start', 'End', 'mean_depth'}
#
#            if depth_df is not None:
#                self.data[sampleid][f'{sampletype}_depth'] = depth_df
#
#        def sanity_check(
#            mode, target_region, vcf_sampleid_tumor, vcf_sampleid_normal, germline_vcf_path,
#        ):
#            # target_region
#            assert mode in ('wgs', 'panel')
#            if (mode == 'panel') and (target_region is None):
#                raise Exception(f'"target_region" must be given when "mode" is "panel"')
#            elif (mode == 'wgs') and (target_region is not None):
#                raise Exception(f'"target_region" must not be given when "mode" is "wgs"')
#
#            # germline VCF file arguments
#            if germline_vcf_path is not None:
#                if vcf_sampleid_tumor is None:
#                    raise Exception(f'When "germline_vcf_path" is used, "vcf_sampleid_tumor" must be given.')
#                if (mode == 'wgs') and (vcf_sampleid_normal is None):
#                    raise Exception(f'"vcf_sampleid_normal" must be given when "mode" is "wgs"')
#
#        # main
#        sanity_check(mode, target_region, vcf_sampleid_tumor, vcf_sampleid_normal, germline_vcf_path)
#        self.data[sampleid] = dict()
#        self.data[sampleid]['is_female'] = is_female
#        self.data[sampleid]['mode'] = mode
#
#        # load normal
#        depth_loading_helper(
#            self, 
#            normal_bam_path, 
#            normal_depth_path, 
#            normal_depth_df,
#            sampleid, 
#            'normal', 
#        )
#
#        # load tumor
#        depth_loading_helper(
#            self, 
#            tumor_bam_path, 
#            tumor_depth_path, 
#            tumor_depth_df,
#            sampleid, 
#            'tumor', 
#        )
#
#        # set target_region
#        self.set_target_region(sampleid, mode, target_region)
#
#        # normal mean ploidy
#        self.set_normal_mean_ploidy(sampleid)
#
#        # load germline vcf
#        if germline_vcf_path is not None:
#            LOGGER_INFO.info('Loading tumor germline vcf')
#            self.load_germline_vcf(
#                sampleid=sampleid, 
#                vcf_path=germline_vcf_path, 
#                vcf_sampleid_tumor=vcf_sampleid_tumor,
#                vcf_sampleid_normal=vcf_sampleid_normal,
#                logging_lineno=50000,
#                nproc=vcfload_nproc,
#            )
#        elif germline_vafdf_path is not None:
#            LOGGER_INFO.info('Loading tumor germline vaf dataframe')
#            self.load_germline_vafdf(sampleid, germline_vafdf_path)
#
#        #self.postprocess_bafdf(sampleid)
#
#        # postprocess depths
#        self.postprocess_depth(sampleid, verbose=verbose, norm_method=norm_method)
#
#    def load_germline_vafdf(self, sampleid, vafdf_path):
#        self.data[sampleid]['original_baf'] = pd.read_csv(
#            vafdf_path,
#            sep='\t',
#            dtype={
#                'Chromosome': 'string',  
#                'Start': int,   
#                'End': int,     
#                'vaf_raw_tumor': float,   
#                'baf_raw_tumor': float,   
#                'vaf_raw_normal': float,  
#                'baf_raw_normal': float,
#            },
#        )
#
#    def postprocess_bafdf(self, sampleid):
#        modified_baf = self.data[sampleid]['original_baf'].copy()
#        #self.data[sampleid]['baf'] = modified_baf.loc[modified_baf['baf_raw_tumor'] > 0, :]
#
#    def add_bafpeak_to_segment(self, sampleid, bw=1):
#        # join segment and raw bafs
#        left = self.data[sampleid]['baf_segment']
#
#        right = self.data[sampleid]['original_baf']
#        right = right.loc[
#            right['baf_raw_tumor'] > 0, 
#            ['Chromosome', 'Start', 'End', 'baf_raw_tumor'],
#        ]
#
#        joined = pyranges_helper.join(
#            left, right, how='left', merge=None, sort=True, refver=self.refver,
#        )
#        # find peaks
#        groupkey = cnvmisc.genome_df_groupkey(joined, refver=self.refver)
#        peaks = list()
#        for k, v in joined.groupby(groupkey)[['baf_raw_tumor']]:
#            baf_values = v['baf_raw_tumor'].to_numpy()
#            peaks.append(
#                libbaf.infer_baf_density(baf_values, bw=bw, rmzero=False)
#            )
#        # assign values
#        assert len(peaks) == left.shape[0], (
#            f'The number of groupby groups and segment df row number are different'
#        )
#        left['baf_segment_peak'] = peaks
#        self.data[sampleid]['baf_segment'] = left
#
#    def postprocess_depth(self, sampleid, norm_method='plain', verbose=False):
#        logger = (LOGGER_DEBUG if verbose else LOGGER_INFO)
#
#        # get GC df
#        if norm_method == 'plain':
#            gc_df = None
#        else:
#            if self.data[sampleid]['mode'] == 'wgs':
#                logger.info(f'Getting gc fraction dataframe')
#                gc_df = libgcfraction.get_gc_df(
#                    self.refver, 
#                    self.default_binsize, 
#                    coords_as_index=True,
#                )
#            elif self.data[sampleid]['mode'] == 'panel':
#                gc_df = None
#
#        # main
#        for key in ('normal', 'tumor'):
#            if f'{key}_depth' not in self.data[sampleid].keys():
#                continue
#
#            logger.info(f'Beginning postprocess of {key} depth')
#            refver_arg = (
#                self.refver 
#                if self.data[sampleid]['mode'] == 'panel' else 
#                None
#            )
#            output_depth_df, gcbin_average_depths = cnvmisc.postprocess_depth_df(
#                self.data[sampleid][f'{key}_depth'],
#                refver=refver_arg,
#                gc_df=gc_df,
#                included_region=self.data[sampleid]['target_region'],
#                as_gr=False,
#                verbose=verbose,
#                norm_method=norm_method,
#                #add_norm_depth=True,
#            )
#
#            self.data[sampleid][f'{key}_depth'] = output_depth_df
#            self.data[sampleid][f'{key}_gcdata'] = gcbin_average_depths
#
#    def set_depthratio(self, sampleid):
#        depthratio_df = cnvmisc.make_depth_ratio(
#            self.data[sampleid]['tumor_depth'], 
#            self.data[sampleid]['normal_depth'],
#            #make_depthratio_mystyle=False,
#            #make_depthratio_plain=True,
#            as_gr=False,
#        )
#        depthratio_df.rename(
#            columns={
#                #'depth_ratio_sequenzastyle': 'depthratio_raw_seqzstyle',
#                #'depth_ratio_plain': 'depthratio_raw_plain',
#                'depthratio': 'depthratio_raw',
#            }, 
#            inplace=True,
#        )
#        self.data[sampleid]['depthratio'] = depthratio_df
#
#    def upscale_preprocessing(self, input_df):
#        result = input_df.copy()
#        annot_cols = cnvmisc.get_genome_df_annotcols(input_df)
#        result.loc[result['excluded'], annot_cols] = np.nan
#        result.drop('excluded', axis=1, inplace=True)
#        return result
#
#    def upscale_depthratio(self, sampleid, binsize=1000):
#        #input_df = self.upscale_preprocessing(self.data[sampleid]['depthratio'])
#        input_df = self.data[sampleid]['depthratio']
#        self.data[sampleid]['depthratio_upscaled'] = cnvmisc.upsize_depth_df_bin(
#            input_df, 
#            size=binsize, 
#            refver=self.refver,
#        )
#
#    def upscale_depth(self, sampleid, binsize=1000, do_normal=True, do_tumor=True):
#        if do_normal:
#            self.data[sampleid]['normal_depth_upscaled'] = cnvmisc.upsize_depth_df_bin(
#                self.data[sampleid]['normal_depth'], 
#                size=binsize, 
#                refver=self.refver,
#            )
#
#        if do_tumor:
#            self.data[sampleid]['tumor_depth_upscaled'] = cnvmisc.upsize_depth_df_bin(
#                self.data[sampleid]['tumor_depth'], 
#                size=binsize, 
#                refver=self.refver,
#            )
#
#    def make_segments(
#        self,
#        sampleid,
#        winsorize=False,
#        depthratio_gamma=None,
#        depthratio_kmin=None,
#        baf_gamma=100,
#        baf_kmin=None,
#        verbose=False,
#        segment_baf_cutoff=0.1,
#
#        bafcorrection_cutoff=None,
#        bafcorrector=None,
#
#        bw=1,
#    ):
#        depthratio_df = (
#            self.data[sampleid]['depthratio_upscaled']
#            if 'depthratio_upscaled' in self.data[sampleid] else
#            self.data[sampleid]['depthratio']
#        )
#        depth_segment, baf_segment = _make_segments_main(
#            depthratio_df=depthratio_df,
#            mode=self.data[sampleid]['mode'],
#            refver=self.refver,
#            winsorize=winsorize,
#            depthratio_gamma=depthratio_gamma,
#            depthratio_kmin=depthratio_kmin,
#            baf_gamma=baf_gamma,
#            baf_kmin=baf_kmin,
#            verbose=verbose,
#
#            baf_df=self.data[sampleid]['original_baf'],
#            target_region=self.data[sampleid]['target_region'],
#            baf_cutoff=segment_baf_cutoff,
#        )
#        self.data[sampleid]['depthratio_segment'] = depth_segment
#        self.data[sampleid]['baf_segment'] = baf_segment
#
#        self.make_segments_postprocess(
#            sampleid, 
#            bw=bw,
#            bafcorrection_cutoff=bafcorrection_cutoff,
#            bafcorrector=bafcorrector,
#        )
#
#    def make_segments_postprocess(
#        self, sampleid, bw=1,
#        bafcorrection_cutoff=None,
#        bafcorrector=None,
#    ):
#        # add baf peak
#        self.add_bafpeak_to_segment(sampleid, bw=bw)
#
#        # make merged segments
#        self.data[sampleid]['merged_segment'] = self.make_merged_segment(
#            sampleid,
#            self.data[sampleid]['depthratio_segment'],
#        )
#        
#        # add corrected baf
#        self.add_corrected_baf(
#            sampleid, 
#            round_cutoff=bafcorrection_cutoff,
#            bafcorrector=bafcorrector,
#        )
#
#    def add_corrected_baf(
#        self, 
#        sampleid, 
#        round_cutoff=None,
#        bafcorrector=None,
#    ):
#        if bafcorrector is None:
#            bafcorrector = libbaf.load_bafcorrect_func(x_cutoff=round_cutoff)
#
#        for df in (
#            self.data[sampleid]['baf_segment'], 
#            self.data[sampleid]['merged_segment'], 
#        ):
#            df['corrected_baf_segment_mean'] = bafcorrector(df['baf_segment_peak'].to_numpy())
#
#    def get_cp_from_twodata(
#        self, sampleid, depthratio1, CNt1, depthratio2, CNt2,
#    ):
#        return cnvmisc.get_cp_from_twodata(
#            depthratio1, 
#            CNt1, 
#            depthratio2, 
#            CNt2, 
#            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'], 
#            CNn=2, 
#            tumor_avg_depth_ratio=1, 
#            normal_avg_depth_ratio=1,
#        )
#
#    def calculate_tumor_ploidy(self, sampleid, segment_df):
#        if self.data[sampleid]['mode'] == 'wgs':
#            weights = segment_df['End'] - segment_df.start0s
#        else:
#            stripped_segment_df = segment_df.loc[:, ['Chromosome', 'Start', 'End']].copy()
#            all_indexes = list(range(stripped_segment_df.shape[0]))  # index of segments
#            stripped_segment_df['index'] = all_indexes
#
#            target_region = cnvmisc.arg_into_gr(self.data[sampleid]['target_region'])
#            index_annotated_targetregion_df = cnvmisc.annotate_region_with_segment(  
#                # each region is annotated with corresponding segment index
#                target_region[[]],
#                stripped_segment_df,
#                as_gr=False,
#            )
#
#            index_annotated_targetregion_df['length'] = (
#                index_annotated_targetregion_df['End']
#                - index_annotated_targetregion_df.start0s
#            )
#            weights_dict = index_annotated_targetregion_df.loc[
#                :, ['length', 'index']
#            ].groupby('index').sum().to_dict()['length']
#            weights = [
#                (weights_dict[x] if x in weights_dict else 0)
#                for x in all_indexes   
#            ]
#
#        return np.average(segment_df['CNt'], weights=weights)
#
#    @plotter_decorator
#    def plot_beforecp(
#        self, 
#        sampleid, 
#        figsize=None, 
#        hspace=None,
#        draw_invalid_regions=False, 
#        use_saved_plotdata=False,
#        use_merged_segment=True,
#
#        draw_depthratio_hist=True,
#        #rm_haploid_from_hist=True,
#
#        depthratio_hist_threshold=None,
#        depthratio_hist_bw=None,
#
#        draw_depth=False,
#        is_rawdepth=True,
#        depth_binsize=10000,
#        depth_ymax=None,
#
#        n_xlabel=None,
#        depthratio_ymax=None,
#
#        depthratio_dot_kwargs=dict(),
#        depthratio_line_segmean_kwargs=dict(),
#        #depthratio_line_predict_kwargs=dict(),
#
#        depthratio_hist_annotate_kwargs=dict(),
#        depthratio_hist_plot_kwargs=dict(),
#
#        depth_dot_kwargs=dict(),
#
#        baf_dot_kwargs=dict(),
#        baf_line_segmean_kwargs=dict(),
#        baf_line_corr_segmean_kwargs=dict(),
#        #baf_line_predict_kwargs=dict(),
#
#        draw_upscaled_depthratio=True,
#        draw_upscaled_depth=True,
#    ):
#        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
#        self.make_plotdata_basic(sampleid, draw_upscaled_depthratio)
#        if draw_depth:
#            self.make_plotdata_fordepth(sampleid, use_upscaled=draw_upscaled_depth, binsize=depth_binsize)
#        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')
#
#        fig, axd, bottom_axkey = self.make_axd(
#            figsize=figsize, 
#            hspace=hspace, 
#            draw_depthratio_hist=draw_depthratio_hist, 
#            draw_solution=False,
#            draw_tumor_baf=True,
#            draw_depthratio=True,
#            draw_normal_baf=draw_depth,
#            draw_tumor_depth=draw_depth,
#            draw_normal_depth=draw_depth,
#        )
#
#        fig.suptitle(
#            f'sample_id={sampleid}, is_female={self.data[sampleid]["is_female"]}',
#            fontsize=20,
#        )
#
#        self.draw_baf_ax(
#            sampleid, 
#            axd['baf'], 
#
#            use_merged_segment=True,
#
#            draw_predicted=False,
#            draw_corrected=True,
#            draw_segmean=True,
#            n_xlabel=(n_xlabel if 'baf' == bottom_axkey else None),
#
#            is_tumor=True,
#
#            mark_unfit_regions=False,
#
#            dot_kwargs=baf_dot_kwargs,
#            line_segmean_kwargs=baf_line_segmean_kwargs,
#            line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,
#        )
#        self.draw_depthratio_ax(
#            sampleid, 
#            axd['depthratio'], 
#            use_merged_segment=True,
#
#            draw_predicted=False,
#            draw_segmean=True,
#            draw_deviation=False,
#            draw_depthratio_peaks=False,
#
#            n_xlabel=(n_xlabel if 'depthratio' == bottom_axkey else None),
#            dot_kwargs=depthratio_dot_kwargs,
#            line_segmean_kwargs=depthratio_line_segmean_kwargs,
#            ymax=depthratio_ymax,
#
#            use_upscaled=draw_upscaled_depthratio,
#        )
#
#        if draw_depthratio_hist:
#            peak_depthratios = self.draw_depthratio_hist_ax(
#                sampleid,
#                axd['depthratio_hist'],
#                use_merged_segment=True, 
#                depth_ylim=axd['depthratio'].get_ylim(),
#                #rm_haploid=rm_haploid_from_hist,
#                peak_threshold=depthratio_hist_threshold,
#                bw=depthratio_hist_bw,
#                annotate_kwargs=depthratio_hist_annotate_kwargs,
#                plot_kwargs=depthratio_hist_plot_kwargs,
#            )
#
#            for y in peak_depthratios:
#                axd['depthratio'].axhline(y, color='orange', linewidth=1, alpha=0.6)
#
#        if draw_depth:
#            self.draw_depth_bundle(
#                sampleid, axd, n_xlabel, is_rawdepth, 
#                depth_dot_kwargs, depth_ymax, 
#                baf_dot_kwargs, baf_line_segmean_kwargs, baf_line_corr_segmean_kwargs,
#                bottom_axkey,
#                draw_upscaled_depth,
#            )
#
#        return fig, axd
#
#    @plotter_decorator
#    def plot_aftercp_freeccf(
#        self, 
#        sampleid, 
#        cellularity,
#        ploidy,
#
#        figsize=(30, 16), 
#        hspace=None,
#
#        n_xlabel=None,
#        depthratio_ymax=None,
#        CN_ymax=None,
#        subCN_ymax=None,
#
#        depthratio_std_factor=1,
#
#        draw_depth=False,
#        is_rawdepth=True,
#        depth_binsize=10000,
#        depth_ymax=None,
#
#        depthratio_dot_kwargs=dict(),
#        depthratio_line_segmean_kwargs=dict(),
#        depthratio_line_predict_kwargs=dict(),
#        depthratio_line_predict_clonal_kwargs=dict(),
#
#        baf_dot_kwargs=dict(),
#        baf_line_segmean_kwargs=dict(),
#        baf_line_corr_segmean_kwargs=dict(),
#        baf_line_predict_kwargs=dict(),
#        baf_line_predict_clonal_kwargs=dict(),
#
#        CN_line_CNt_kwargs=dict(),
#        subCN_line_CNt_kwargs=dict(),
#        ccf_bar_kwargs=dict(),
#
#        depth_ratio_diff=None,
#        baf_diff=0.05,
#
#        Bn=1,
#
#        min_N_CNt_candidates=5,
#        N_CNt_candidates_fraction=0.5,
#
#        limited_clonal=True,
#
#        draw_upscaled_depthratio=True,
#        draw_upscaled_depth=True,
#    ):
#        LOGGER_INFO.info(f'Beginning calculation of subclonal solution')
#        self.make_CN_solution_freeccf(
#            sampleid,
#            cellularity,
#            ploidy,
#            depth_ratio_diff=depth_ratio_diff,
#            baf_diff=baf_diff,
#            min_N_CNt_candidates=min_N_CNt_candidates,
#            N_CNt_candidates_fraction=N_CNt_candidates_fraction,
#            limited_clonal=limited_clonal,
#        )
#        self.add_freeccf_solution_to_segment(sampleid)
#
#        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
#        self.make_plotdata_basic(sampleid, draw_upscaled_depthratio)
#        if draw_depth:
#            self.make_plotdata_fordepth(sampleid, use_upscaled=draw_upscaled_depth, binsize=depth_binsize)
#        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')
#
#        fig, axd, bottom_axkey = self.make_axd(
#            figsize=figsize, 
#            hspace=hspace, 
#            draw_depthratio_hist=False, 
#            draw_solution=True,
#            draw_tumor_baf=True,
#            draw_depthratio=True,
#            draw_normal_baf=draw_depth,
#            draw_tumor_depth=draw_depth,
#            draw_normal_depth=draw_depth,
#        )
#
#        fig.suptitle(
#            ', '.join([
#                f'sample_id={sampleid}',
#                f'is_female={self.data[sampleid]["is_female"]}',
#                f'cellularity={round(cellularity, 3)}',
#                f'ploidy={round(ploidy, 3)}',
#            ]),
#            fontsize=20,
#        )
#
#        # depth
#        self.draw_depthratio_ax(
#            sampleid, 
#            axd['depthratio'], 
#            use_merged_segment=True,
#
#            draw_predicted=True,
#            draw_segmean=True,
#            draw_deviation=False,
#            draw_depthratio_peaks=False,
#
#            mark_unfit_regions=True,
#
#            cellularity=cellularity,
#            ploidy=ploidy,
#            draw_integerCN_lines=True,
#
#            std_factor=depthratio_std_factor,
#            n_xlabel=n_xlabel,
#            dot_kwargs=depthratio_dot_kwargs,
#            line_segmean_kwargs=depthratio_line_segmean_kwargs,
#            line_predict_kwargs=depthratio_line_predict_kwargs,
#            line_predict_clonal_kwargs=depthratio_line_predict_clonal_kwargs,
#            ymax=depthratio_ymax,
#
#            use_upscaled=draw_upscaled_depthratio,
#        )
#        if draw_depth:
#            self.draw_depth_bundle(
#                sampleid, axd, n_xlabel, is_rawdepth, 
#                depth_dot_kwargs, depth_ymax, 
#                baf_dot_kwargs, baf_line_segmean_kwargs, baf_line_corr_segmean_kwargs,
#                bottom_axkey,
#                draw_upscaled_depth,
#            )
#
#        # baf
#        self.draw_baf_ax(
#            sampleid, 
#            axd['baf'], 
#            use_merged_segment=True,
#            draw_predicted=True,
#            draw_corrected=True,
#            draw_segmean=True,
#            n_xlabel=n_xlabel,
#
#            is_tumor=True,
#            mark_unfit_regions=True,
#
#            dot_kwargs=baf_dot_kwargs,
#            line_segmean_kwargs=baf_line_segmean_kwargs,
#            line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,
#            line_predict_kwargs=baf_line_predict_kwargs,
#            line_predict_clonal_kwargs=baf_line_predict_clonal_kwargs,
#        )
#
#        # clonal CN
#        self.draw_CN_ax(
#            sampleid, 
#            axd['clonal_CN'],
#            n_xlabel=n_xlabel,
#            line_CNt_kwargs=CN_line_CNt_kwargs,
#            ymax=CN_ymax,
#            draw_CNt=True,
#            draw_A=False,
#            draw_B=True,
#        )
#        # subclonal CN
#        self.draw_subclonal_CN_ax(
#            sampleid, 
#            axd['subclonal_CN'],
#            n_xlabel=n_xlabel,
#            line_CNt_kwargs=subCN_line_CNt_kwargs,
#            ymax=subCN_ymax,
#            draw_CNt=True,
#            draw_A=False,
#            draw_B=True,
#        )
#        # ccf
#        self.draw_ccf_ax(
#            sampleid,
#            axd['ccf'],
#            n_xlabel=n_xlabel,
#            bar_kwargs=ccf_bar_kwargs,
#        )
#
#        return fig, axd
#
#    def show_ccfs(self, sampleid, bandwidth=0.1):
#        self.select_fixed_ccfs(sampleid, bandwidth=bandwidth)
#        ccf_plotdata = self.data[sampleid]['ccf_plotdata']
#        self.show_ccfs_main(
#            ccfs=ccf_plotdata['ccfs'],
#            lengths=ccf_plotdata['lengths'],
#            density=ccf_plotdata['density'],
#            peak_values=ccf_plotdata['peak_values'],
#        )
#
#    @plotter_decorator
#    def plot_aftercp_fixedccf(
#        self, 
#        sampleid, 
#        cellularity,
#        ploidy,
#
#        figsize=(30, 16), 
#        hspace=None,
#
#        n_xlabel=None,
#        depthratio_ymax=None,
#        CN_ymax=None,
#        subCN_ymax=None,
#
#        depthratio_std_factor=1,
#
#        draw_depth=False,
#        is_rawdepth=True,
#        depth_binsize=10000,
#        depth_ymax=None,
#
#        depthratio_dot_kwargs=dict(),
#        depthratio_line_segmean_kwargs=dict(),
#        depthratio_line_predict_kwargs=dict(),
#        depthratio_line_predict_clonal_kwargs=dict(),
#
#        baf_dot_kwargs=dict(),
#        baf_line_segmean_kwargs=dict(),
#        baf_line_corr_segmean_kwargs=dict(),
#        baf_line_predict_kwargs=dict(),
#        baf_line_predict_clonal_kwargs=dict(),
#
#        CN_line_CNt_kwargs=dict(),
#        subCN_line_CNt_kwargs=dict(),
#        ccf_bar_kwargs=dict(),
#
#        depth_ratio_diff=None,
#        baf_diff=0.05,
#
#        min_N_CNt_candidates=5,
#        N_CNt_candidates_fraction=0.5,
#
#        ccf_bw=0.1,
#
#        update_plotdata=False,
#        CNt_diff_factor=0.1,
#
#        mark_unfit_regions=False,
#
#        limited_clonal=True,
#
#        draw_upscaled_depthratio=True,
#        draw_upscaled_depth=True,
#
#        merge_same_chroms=True,
#    ):
#        LOGGER_INFO.info(f'Beginning calculation of subclonal solution')
#        self.make_CN_solution_after_ccfs(
#            sampleid,
#            cellularity,
#            ploidy,
#            min_N_CNt_candidates=min_N_CNt_candidates,
#            N_CNt_candidates_fraction=N_CNt_candidates_fraction,
#            CNt_diff_factor=CNt_diff_factor,
#            limited_clonal=limited_clonal,
#        )
#        self.add_fixedccf_solution_to_segment(sampleid)
#        LOGGER_INFO.info(f'Finished calculation of subclonal solution')
#
#        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
#        if update_plotdata:
#            self.add_solution_to_plotdata(sampleid)
#        else:
#            self.make_plotdata_basic(sampleid, draw_upscaled_depthratio)
#            if draw_depth:
#                self.make_plotdata_fordepth(sampleid, use_upscaled=draw_upscaled_depth, binsize=depth_binsize)
#        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')
#
#        fig, axd, bottom_axkey = self.make_axd(
#            figsize=figsize, 
#            hspace=hspace, 
#            draw_depthratio_hist=False, 
#            draw_solution=True,
#            draw_tumor_baf=True,
#            draw_depthratio=True,
#            draw_normal_baf=draw_depth,
#            draw_tumor_depth=draw_depth,
#            draw_normal_depth=draw_depth,
#        )
#
#        fig.suptitle(
#            ', '.join([
#                f'sample_id={sampleid}',
#                f'is_female={self.data[sampleid]["is_female"]}',
#                f'cellularity={round(cellularity, 3)}',
#                f'ploidy={round(ploidy, 3)}',
#            ]),
#            fontsize=20,
#        )
#
#        # depth
#        self.draw_depthratio_ax(
#            sampleid, 
#            axd['depthratio'], 
#            use_merged_segment=True,
#
#            draw_predicted=True,
#            draw_segmean=True,
#            draw_deviation=False,
#            draw_depthratio_peaks=False,
#
#            mark_unfit_regions=mark_unfit_regions,
#
#            cellularity=cellularity,
#            ploidy=ploidy,
#            draw_integerCN_lines=True,
#
#            std_factor=depthratio_std_factor,
#            n_xlabel=(n_xlabel if 'depthratio' == bottom_axkey else None),
#            dot_kwargs=depthratio_dot_kwargs,
#            line_segmean_kwargs=depthratio_line_segmean_kwargs,
#            line_predict_kwargs=depthratio_line_predict_kwargs,
#            line_predict_clonal_kwargs=depthratio_line_predict_clonal_kwargs,
#            ymax=depthratio_ymax,
#
#            use_upscaled=draw_upscaled_depthratio,
#
#            merge_same_chroms=merge_same_chroms,
#        )
#        if draw_depth:
#            self.draw_depth_bundle(
#                sampleid, axd, n_xlabel, is_rawdepth, 
#                depth_dot_kwargs, depth_ymax, 
#                baf_dot_kwargs, baf_line_segmean_kwargs, baf_line_corr_segmean_kwargs,
#                bottom_axkey,
#                draw_upscaled_depth,
#            )
#
#        # baf
#        self.draw_baf_ax(
#            sampleid, 
#            axd['baf'], 
#            use_merged_segment=True,
#            draw_predicted=True,
#            draw_corrected=True,
#            draw_segmean=True,
#            n_xlabel=(n_xlabel if 'baf' == bottom_axkey else None),
#
#            is_tumor=True,
#            mark_unfit_regions=mark_unfit_regions,
#
#            dot_kwargs=baf_dot_kwargs,
#            line_segmean_kwargs=baf_line_segmean_kwargs,
#            line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,
#            line_predict_kwargs=baf_line_predict_kwargs,
#            line_predict_clonal_kwargs=baf_line_predict_clonal_kwargs,
#
#            merge_same_chroms=merge_same_chroms,
#        )
#
#        # clonal CN
#        self.draw_CN_ax(
#            sampleid, 
#            axd['clonal_CN'],
#            n_xlabel=(n_xlabel if 'clonal_CN' == bottom_axkey else None),
#            line_CNt_kwargs=CN_line_CNt_kwargs,
#            ymax=CN_ymax,
#            draw_CNt=True,
#            draw_A=False,
#            draw_B=True,
#
#            merge_same_chroms=merge_same_chroms,
#        )
#        # subclonal CN
#        self.draw_subclonal_CN_ax(
#            sampleid, 
#            axd['subclonal_CN'],
#            n_xlabel=(n_xlabel if 'subclonal_CN' == bottom_axkey else None),
#            line_CNt_kwargs=subCN_line_CNt_kwargs,
#            ymax=subCN_ymax,
#            draw_CNt=True,
#            draw_A=False,
#            draw_B=True,
#
#            merge_same_chroms=merge_same_chroms,
#        )
#        # ccf
#        self.draw_ccf_ax(
#            sampleid,
#            axd['ccf'],
#            n_xlabel=(n_xlabel if 'ccf' == bottom_axkey else None),
#            bar_kwargs=ccf_bar_kwargs,
#
#            merge_same_chroms=merge_same_chroms,
#        )
#        fixed_ccfs = self.data[sampleid]['fixed_ccfs']
#        for y in fixed_ccfs:
#            axd['ccf'].axhline(y, color='red', linewidth=1)
#
#        return fig, axd
#
#    @plotter_decorator
#    def plot_custom(
#        self, 
#        sampleid, 
#
#        fig=None,
#        axd=None,
#        figsize=(30, 16), 
#        hspace=None,
#        height_ratios=None,
#        title=None,
#        title_size=20,
#        title_y=0.95,
#
#        draw_depthratio_hist=False, 
#        draw_tumor_baf=False,
#        draw_depthratio=False,
#        draw_normal_baf=False,
#        draw_normal_depth=False,
#        draw_tumor_depth=False,
#        draw_feature=False,
#
#        tumor_baf_ylabel=None,
#        depthratio_ylabel=None,
#        normal_baf_ylabel=None,
#        normal_depth_ylabel=None,
#        tumor_depth_ylabel=None,
#        feature_ylabel=None,
#        tumor_baf_ylabel_kwargs=dict(),
#        depthratio_ylabel_kwargs=dict(),
#        normal_baf_ylabel_kwargs=dict(),
#        normal_depth_ylabel_kwargs=dict(),
#        tumor_depth_ylabel_kwargs=dict(),
#        feature_ylabel_kwargs=dict(),
#
#        n_xlabel=None,
#        depthratio_ymax=None,
#        CN_ymax=None,
#        subCN_ymax=None,
#
#        is_rawdepth=True,
#        depth_binsize=1000,
#        depth_ymax=None,
#        use_upscaled_depthratio=True,
#        use_upscaled_tumor_depth=True,
#        use_upscaled_normal_depth=True,
#
#        depthratio_dot_kwargs=dict(),
#        depthratio_line_segmean_kwargs=dict(),
#        depthratio_line_predict_kwargs=dict(),
#        depthratio_line_predict_clonal_kwargs=dict(),
#
#        baf_dot_kwargs=dict(),
#        baf_line_segmean_kwargs=dict(),
#        baf_line_corr_segmean_kwargs=dict(),
#        baf_line_predict_kwargs=dict(),
#        baf_line_predict_clonal_kwargs=dict(),
#
#        feature_text_kwargs=dict(),
#        feature_line_kwargs=dict(),
#
#        CN_line_CNt_kwargs=dict(),
#        subCN_line_CNt_kwargs=dict(),
#        ccf_bar_kwargs=dict(),
#
#        normal_depth_dot_kwargs=dict(),
#        tumor_depth_dot_kwargs=dict(),
#
#        feature_df=None,
#        draw_feature_label=True,
#        #feature_as_dot=False,
#
#        #depth_ratio_diff=None,
#        #baf_diff=0.05,
#
#        #min_N_CNt_candidates=5,
#        #N_CNt_candidates_fraction=0.5,
#
#        #ccf_bw=0.1,
#
#        #update_plotdata=False,
#        #CNt_diff_factor=0.1,
#
#        #mark_unfit_regions=False,
#
#        split_spines=True,
#        merge_same_chroms=True,
#    ):
#        # make plotdata
#        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
#        if draw_tumor_baf:
#            self.make_tumor_baf_plotdata(sampleid)
#        if draw_depthratio:
#            self.make_depthratio_plotdata(sampleid, use_upscaled=use_upscaled_depthratio)
#        if draw_normal_baf:
#            self.make_normal_baf_plotdata(sampleid)
#        if draw_normal_depth:
#            self.make_depth_plotdata(
#                sampleid, 
#                is_tumor=False, 
#                use_upscaled=use_upscaled_normal_depth, 
#                binsize=depth_binsize,
#            )
#        if draw_tumor_depth:
#            self.make_depth_plotdata(
#                sampleid, 
#                is_tumor=True, 
#                use_upscaled=use_upscaled_tumor_depth, 
#                binsize=depth_binsize,
#            )
#        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')
#
#        # main
#        assert not (
#            (fig is None)
#            and (axd is not None)
#        )
#        fig_not_given = (fig is None)
#        axd_not_given = (axd is None)
#        if axd is None:
#            fig, axd, bottom_axkey = self.make_axd(
#                figsize=figsize, 
#                hspace=hspace, 
#                height_ratios=height_ratios,
#                draw_depthratio_hist=draw_depthratio_hist, 
#                draw_solution=False,
#                draw_tumor_baf=draw_tumor_baf,
#                draw_depthratio=draw_depthratio,
#                draw_normal_baf=draw_normal_baf,
#                draw_tumor_depth=draw_tumor_depth,
#                draw_normal_depth=draw_normal_depth,
#                draw_feature=draw_feature,
#                draw_xlabel=(n_xlabel is not None),
#                fig=fig,
#            )
#            
#        if fig_not_given:
#            if title is None:
#                title = f'sample_id={sampleid}, is_female={self.data[sampleid]["is_female"]}'
#            fig.suptitle(title, fontsize=title_size, y=title_y)
#
#        if draw_tumor_baf:
#            self.draw_baf_ax(
#                sampleid, 
#                axd['baf'], 
#                use_merged_segment=True,
#
#                draw_predicted=False,
#                draw_corrected=False,
#                draw_segmean=False,
#
#                n_xlabel=n_xlabel,
#
#                is_tumor=True,
#
#                mark_unfit_regions=False,
#
#                dot_kwargs=baf_dot_kwargs,
#                line_segmean_kwargs=baf_line_segmean_kwargs,
#                line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,
#
#                split_spines=split_spines,
#
#                ylabel=tumor_baf_ylabel,
#                ylabel_kwargs=tumor_baf_ylabel_kwargs,
#
#                modify_ax=axd_not_given,
#                merge_same_chroms=merge_same_chroms,
#            )
#        if draw_depthratio:
#            self.draw_depthratio_ax(
#                sampleid, 
#                axd['depthratio'], 
#                use_merged_segment=True,
#
#                draw_predicted=False,
#                draw_segmean=False,
#                draw_deviation=False,
#                draw_depthratio_peaks=False,
#
#                n_xlabel=n_xlabel,
#                dot_kwargs=depthratio_dot_kwargs,
#                line_segmean_kwargs=depthratio_line_segmean_kwargs,
#                ymax=depthratio_ymax,
#
#                split_spines=split_spines,
#
#                ylabel=depthratio_ylabel,
#                ylabel_kwargs=depthratio_ylabel_kwargs,
#
#                modify_ax=axd_not_given,
#                merge_same_chroms=merge_same_chroms,
#
#                use_upscaled=use_upscaled_depthratio,
#            )
#        if draw_normal_baf:
#            self.draw_baf_ax(
#                sampleid, 
#                axd['normal_baf'], 
#
#                use_merged_segment=True,
#
#                draw_predicted=False,
#                draw_corrected=False,
#                draw_segmean=False,
#
#                n_xlabel=n_xlabel,
#
#                is_tumor=False,
#
#                mark_unfit_regions=False,
#
#                dot_kwargs=baf_dot_kwargs,
#                line_segmean_kwargs=baf_line_segmean_kwargs,
#                line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,
#
#                split_spines=split_spines,
#
#                ylabel=normal_baf_ylabel,
#                ylabel_kwargs=normal_baf_ylabel_kwargs,
#
#                modify_ax=axd_not_given,
#                merge_same_chroms=merge_same_chroms,
#            )
#        if draw_normal_depth:
#            self.draw_depth_ax(
#                sampleid,
#                axd['normal_depth'],
#                n_xlabel=n_xlabel,
#                is_tumor=False,
#                is_rawdepth=is_rawdepth,
#                dot_kwargs=normal_depth_dot_kwargs,
#                ymax=depth_ymax,
#
#                split_spines=split_spines,
#
#                ylabel=normal_depth_ylabel,
#                ylabel_kwargs=normal_depth_ylabel_kwargs,
#
#                modify_ax=axd_not_given,
#                merge_same_chroms=merge_same_chroms,
#            )
#        if draw_tumor_depth:
#            self.draw_depth_ax(
#                sampleid,
#                axd['tumor_depth'],
#                n_xlabel=n_xlabel,
#                is_tumor=True,
#                is_rawdepth=is_rawdepth,
#                dot_kwargs=tumor_depth_dot_kwargs,
#                ymax=depth_ymax,
#
#                split_spines=split_spines,
#
#                ylabel=tumor_depth_ylabel,
#                ylabel_kwargs=tumor_depth_ylabel_kwargs,
#
#                modify_ax=axd_not_given,
#                merge_same_chroms=merge_same_chroms,
#            )
#        if draw_feature:
#            self.draw_feature_ax(
#                axd['feature'],
#                feature_df=feature_df,
#
#                n_xlabel=n_xlabel,
#
#                ylabel=feature_ylabel,
#                ylabel_kwargs=feature_ylabel_kwargs,
#
#                text_kwargs=feature_text_kwargs,
#                line_kwargs=feature_line_kwargs,
#
#                split_spines=split_spines,
#                merge_same_chroms=merge_same_chroms,
#
#                #feature_as_dot=feature_as_dot,
#                draw_label=draw_feature_label,
#            )
#
#        return fig, axd
#
#    ############
#    # plotting #
#    ############
#
#    def make_axd(
#        self, 
#        figsize=None, 
#        hspace=None, 
#        height_ratios=None,
#        draw_depthratio_hist=False, 
#        draw_solution=False,
#        draw_tumor_baf=False,
#        draw_depthratio=False,
#        draw_normal_baf=False,
#        draw_tumor_depth=False,
#        draw_normal_depth=False,
#        draw_feature=False,
#        draw_xlabel=False,
#        fig=None,
#    ):
#        # row names
#        row_names = list()
#        if draw_solution:
#            row_names.extend(['ccf', 'subclonal_CN', 'clonal_CN'])
#        if draw_tumor_baf:
#            row_names.append('baf')
#        if draw_depthratio:
#            row_names.append('depthratio')
#        if draw_normal_baf:
#            row_names.append('normal_baf')
#        if draw_normal_depth:
#            row_names.append('normal_depth')
#        if draw_tumor_depth:
#            row_names.append('tumor_depth')
#        if draw_feature:
#            row_names.append('feature')
#
#        # default figsize
#        if figsize is None:
#            figsize = (30, 5 * len(row_names))
#
#        # mosaic
#        if draw_depthratio_hist:
#            depthratio_idx = row_names.index('depthratio')
#            mosaic = [
#                (
#                    [name, 'empty_upper'] 
#                    if idx < depthratio_idx else
#                    (
#                        [name, 'empty_lower'] 
#                        if idx > depthratio_idx else
#                        [name, 'depthratio_hist']
#                    )
#                )
#                for idx, name in enumerate(row_names)
#            ]
#        else:
#            mosaic = [[name,] for name in row_names]
#
#        # gridspec_kw
#        if hspace is None:
#            hspace = (0.4 if draw_xlabel else 0.1)
#        gridspec_kw = dict(
#            hspace=hspace, 
#            height_ratios=height_ratios,
#        )
#        if draw_depthratio_hist:
#            gridspec_kw.update(dict(width_ratios=[1, 0.1], wspace=0.02))
#
#        # result
#        if fig is None:
#            fig, axd = plt.subplot_mosaic(
#                mosaic,
#                figsize=figsize,
#                gridspec_kw=gridspec_kw,
#            )
#        else:
#            axd = fig.subplot_mosaic(
#                mosaic,
#                figsize=figsize,
#                gridspec_kw=gridspec_kw,
#            )
#
#        if draw_depthratio_hist:
#            if 'empty_upper' in axd:
#                axd['empty_upper'].axis('off')
#            if 'empty_lower' in axd:
#                axd['empty_lower'].axis('off')
#
#        bottom_axkey = row_names[-1]
#
#        return fig, axd, bottom_axkey
#
#    @staticmethod
#    def get_yticklabel_size(yticks):
#        return min((200 /len(yticks)), 10)
#
##    def draw_grids(self, ax, ys, line_params=dict(), merge_same_chroms=True):
##        line_params = (
##            dict(color='black', linewidth=0.2, alpha=0.5)
##            | line_params
##        )
##        chroms, start0s, end0s = zip(
##            *self.genomeplotter.cconv.get_chrom_borders(
##                merge_same_chroms=merge_same_chroms,
##            )
##        )
##        for y in ys:
##            ax.hlines(
##                np.repeat(y, len(start0s)), start0s, end0s, 
##                **line_params,
##            )
#
#    #def draw_centromeres(self, ax):
#    #    self.genomeplotter.draw_centromeres(ax)
#
#    def draw_centromeres_type2(self, ax):
#        self.genomeplotter.draw_centromeres_type2(ax)
#
#    def make_plotdata_basic(self, sampleid, use_upscaled):
#        self.make_depthratio_plotdata(sampleid, use_upscaled)
#        self.make_tumor_baf_plotdata(sampleid)
#        self.make_segment_plotdata(sampleid)
#
#    def make_plotdata_fordepth(self, sampleid, use_upscaled=True, binsize=100000):
#        self.make_normal_baf_plotdata(sampleid)
#
#        LOGGER_INFO.info(f'Beginning tumor depth data processing')
#        self.make_depth_plotdata(
#            sampleid, is_tumor=True, use_upscaled=use_upscaled, binsize=binsize
#        )
#        LOGGER_INFO.info(f'Finished tumor depth data processing')
#
#        LOGGER_INFO.info(f'Beginning normal depth data processing')
#        self.make_depth_plotdata(
#            sampleid, is_tumor=False, use_upscaled=use_upscaled, binsize=binsize,
#        )
#        LOGGER_INFO.info(f'Finished normal depth data processing')
#
#    def make_depth_plotdata(self, sampleid, is_tumor, use_upscaled=True, binsize=1000):
#        # set params
#        datakey = ('tumor_depth' if is_tumor else 'normal_depth')
#        if use_upscaled:
#            datakey = datakey + '_upscaled'
#        plotdata_key = ('tumor_depth_plotdata' if is_tumor else 'normal_depth_plotdata')
#        #y_colname = ('mean_depth' if is_rawdepth else 'sequenza_style_norm_mean_depth')
#
#        # upscale raw depth df
#        relevant_chroms = [
#            x for x in self.genomeplotter.cconv.totalregion_df['Chromosome']
#            if not x.startswith('-')
#        ]
#        original_df = self.data[sampleid][datakey]
#        relevant_chroms_df = original_df.loc[
#            original_df['Chromosome'].isin(relevant_chroms), 
#            :
#        ]
#
##        if use_upscaled:
##            assert binsize is not None
##            input_df = cnvmisc.upsize_depth_df_bin(
##                relevant_chroms_df, 
##                size=binsize, 
##                refver=self.refver,
##            )
##        else:
##            input_df = relevant_chroms_df
#
#        # turn into plotdata
#        self.data[sampleid][plotdata_key] = self.genomeplotter.prepare_plot_data(
#            relevant_chroms_df
#        )
#
#    def make_depthratio_plotdata(self, sampleid, use_upscaled=True):
#        depthratio_df = (
#            self.data[sampleid]['depthratio_upscaled']
#            if use_upscaled else
#            self.data[sampleid]['depthratio']
#        )
#
#        relevant_chroms = [
#            x for x in self.genomeplotter.cconv.totalregion_df['Chromosome']
#            if not x.startswith('-')
#        ]
#        relevant_chroms_df = depthratio_df.loc[
#            depthratio_df['Chromosome'].isin(relevant_chroms), 
#            :
#        ]
#
#        self.data[sampleid]['depthratio_raw_plotdata'] = (
#            self.genomeplotter.prepare_plot_data(relevant_chroms_df)
#        )
#
#    def make_tumor_baf_plotdata(self, sampleid):
#        bafdf = self.data[sampleid]['original_baf']
#        tumor_baf_df = bafdf.loc[bafdf['baf_raw_tumor'] > 0, :]
#        self.data[sampleid]['baf_raw_plotdata'] = (
#            self.genomeplotter.prepare_plot_data(tumor_baf_df)
#        )
#
#    def make_normal_baf_plotdata(self, sampleid):
#        bafdf = self.data[sampleid]['original_baf']
#        normal_baf_df = bafdf.loc[bafdf['baf_raw_normal'] > 0, :]
#        self.data[sampleid]['normal_baf_raw_plotdata'] = (
#            self.genomeplotter.prepare_plot_data(normal_baf_df)
#        )
#
#    def make_segment_plotdata(self, sampleid):
#        self.data[sampleid]['merged_segment_plotdata'] = (
#            self.genomeplotter.prepare_plot_data(
#                self.data[sampleid]['merged_segment']
#            )
#        )
#
#    def make_plotdata_aftercp_wobaf(self, sampleid, use_saved_plotdata):
#        raw_plot_map = {
#            'depthratio_upscaled': 'depthratio_raw_plotdata',
#            'merged_segment': 'merged_segment_plotdata',
#        }
#        def helper(raw_key):
#            plot_key = raw_plot_map[raw_key]
#            if (
#                (not use_saved_plotdata) 
#                or (plot_key not in self.data[sampleid])
#            ):
#                self.data[sampleid][plot_key] = self.genomeplotter.prepare_plot_data(
#                    self.data[sampleid][raw_key]
#                )
#
#        helper('depthratio_upscaled')
#        helper('merged_segment')
#
#    def draw_feature_ax(
#        self, 
#        ax, 
#        feature_df,
#
#        n_xlabel=None,
#
#        ylabel=None,
#        ylabel_kwargs=dict(),
#
#        #feature_as_dot=False,
#        draw_label=True,
#
#        text_kwargs=dict(),
#        line_kwargs=dict(),
#
#        split_spines=True,
#        merge_same_chroms=True,
#
#        chromlabel_kwargs=dict(), 
#    ):
#        if ylabel is None:
#            ylabel = 'features'
#
#        ax.set_ylim(0, 1)
#        ax.set_ylabel(ylabel, **ylabel_kwargs)
#        ax.set_yticks([])
#        self.genomeplotter.draw_ax_common(
#            ax, 
#            n_xlabel=n_xlabel, 
#            split_spines=split_spines,
#            merge_same_chroms=merge_same_chroms,
#            chromlabel_kwargs=chromlabel_kwargs,
#        )
#
#        self.genomeplotter.draw_features(
#            ax,
#            df=feature_df,
#
#            draw_label=draw_label,
#            #feature_as_dot=feature_as_dot,
#
#            y_features=None,
#            y_labels=None,
#
#            text_kwargs=text_kwargs,
#            line_kwargs=line_kwargs,
#        )
#
#    def draw_depth_bundle(
#        self, 
#        sampleid, axd, n_xlabel, is_rawdepth, 
#        depth_dot_kwargs, depth_ymax,
#        baf_dot_kwargs, baf_line_segmean_kwargs, baf_line_corr_segmean_kwargs,
#        bottom_axkey,
#        use_upscaled,
#    ):
#        self.draw_depth_ax(
#            sampleid,
#            axd['tumor_depth'],
#            n_xlabel=(n_xlabel if 'tumor_depth' == bottom_axkey else None),
#            is_tumor=True,
#            is_rawdepth=is_rawdepth,
#            dot_kwargs=depth_dot_kwargs,
#            ymax=depth_ymax,
#
#            use_upscaled=use_upscaled,
#        )
#        self.draw_depth_ax(
#            sampleid,
#            axd['normal_depth'],
#            n_xlabel=(n_xlabel if 'normal_depth' == bottom_axkey else None),
#            is_tumor=False,
#            is_rawdepth=is_rawdepth,
#            dot_kwargs=depth_dot_kwargs,
#            ymax=depth_ymax,
#
#            use_upscaled=use_upscaled,
#        )
#        self.draw_baf_ax(
#            sampleid, 
#            axd['normal_baf'], 
#
#            use_merged_segment=True,
#
#            draw_predicted=False,
#            draw_corrected=False,
#            draw_segmean=False,
#
#            n_xlabel=(n_xlabel if 'normal_baf' == bottom_axkey else None),
#
#            is_tumor=False,
#
#            mark_unfit_regions=False,
#
#            dot_kwargs=baf_dot_kwargs,
#            line_segmean_kwargs=baf_line_segmean_kwargs,
#            line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,
#
#            make_segment_plotdata=False,
#        )
#
#    def draw_depth_ax(
#        self,
#        sampleid,
#        ax,
#        n_xlabel,
#        is_tumor,
#        is_rawdepth,
#        dot_kwargs=dict(),
#        ymax=None,
#        add_color=True,
#        split_spines=True,
#        ylabel=None,
#        ylabel_kwargs=dict(),
#
#        make_raw_plotdata=True,
#        use_upscaled=True,
#
#        modify_ax=True,
#        merge_same_chroms=True,
#
#        chromlabel_kwargs=dict(),
#    ):
#        # prepare plotdata
#        if make_raw_plotdata:
#            self.make_depth_plotdata(sampleid, is_tumor=is_tumor, use_upscaled=use_upscaled)
#
#        # set params
#        y_colname = (
#            'mean_depth' 
#            if is_rawdepth else 
#            'norm_mean_depth'
#        )
#        plotdata_key = ('tumor_depth_plotdata' if is_tumor else 'normal_depth_plotdata')
#        plotdata_df = self.data[sampleid][plotdata_key]
#
#        # calc dot alpha
#        n_dots = plotdata_df.shape[0]
#        default_alpha = calc_dot_alpha_depth(n_dots)
#
#        # handle kwargs
#        dot_kwargs = (
#            #{'color': 'black', 'markersize': 0.3, 'alpha': 0.05}
#            {'color': 'black', 'markersize': 0.3, 'alpha': default_alpha}
#            | dot_kwargs
#        )
#
#        # draw raw data
#        if add_color:
#            colors = np.repeat('blue', plotdata_df.shape[0])
#            colors[np.where(plotdata_df['excluded'])[0]] = 'red'
#            self.genomeplotter.draw_dots_scatter(
#                ax, 
#                df_plotdata=plotdata_df,
#                y_colname=y_colname, 
#                color_vals=colors,
#                plot_kwargs=dot_kwargs,
#            )
#        else:
#            self.genomeplotter.draw_dots(
#                ax, 
#                df_plotdata=plotdata_df,
#                y_colname=y_colname, 
#                plot_kwargs=dot_kwargs,
#            )
#
#        # set axes attributes
#        if modify_ax:
#            if ylabel is None:
#                ylabel_1 = ('tumor' if is_tumor else 'normal')
#                ylabel_2 = ('raw' if is_rawdepth else 'normalized')
#                ylabel = f'{ylabel_1} sample {ylabel_2} depth'
#            ax.set_ylabel(ylabel, **ylabel_kwargs)
#
#            if ymax is None:
#                ymax = np.nanmean(plotdata_df[y_colname]) * 2
#            ax.set_ylim(-ymax * 0.1, ymax)
#            roundnum = (1 if is_rawdepth else 2)
#            yticks = np.round(np.linspace(0, ymax, 10), roundnum)
#
#            ax.set_yticks(yticks)
#            ax.set_yticklabels(yticks)
#
#            self.genomeplotter.draw_ax_common(
#                ax, 
#                n_xlabel=n_xlabel, 
#                split_spines=split_spines, 
#                merge_same_chroms=merge_same_chroms,
#                chromlabel_kwargs=chromlabel_kwargs,
#            )
#
#    def draw_depthratio_ax(
#        self, 
#        sampleid, 
#        ax, 
#
#        use_merged_segment=True, 
#
#        cellularity=None,
#        ploidy=None,
#
#        draw_integerCN_lines=False,
#
#        draw_predicted=True,
#        draw_segmean=True,
#        draw_deviation=False,
#        draw_depthratio_peaks=False,
#
#        std_factor=1,
#
#        mark_unfit_regions=False,
#
#        n_xlabel=None,
#        dot_kwargs=dict(),
#        line_segmean_kwargs=dict(),
#        line_predict_kwargs=dict(),
#        line_predict_clonal_kwargs=dict(),
#        ymax=None,
#
#        split_spines=True,
#
#        ylabel=None,
#        ylabel_kwargs=dict(),
#
#        make_raw_plotdata=True,
#        make_segment_plotdata=True,
#        use_upscaled=True,
#
#        modify_ax=True,
#        merge_same_chroms=True,
#
#        chromlabel_kwargs=dict(),
#    ):
#        # determine if segment information is plotted
#        draw_segment_info = any([
#            draw_predicted,
#            draw_segmean,
#        ])
#
#        # prepare plotdata
#        if make_raw_plotdata:
#            self.make_depthratio_plotdata(sampleid, use_upscaled=use_upscaled)
#        if make_segment_plotdata:
#            self.make_segment_plotdata(sampleid)
#
#        # calc dot alpha
#        raw_plotdata = self.data[sampleid]['depthratio_raw_plotdata']
#        n_dots = raw_plotdata.shape[0]
#        default_alpha = calc_dot_alpha_depth(n_dots)
#
#        # handle kwargs
#        dot_kwargs = (
#            {'color': 'black', 'markersize': 0.3, 'alpha': default_alpha}
#            | dot_kwargs
#        )
#        line_segmean_kwargs = (
#            {'color': 'tab:blue', 'linewidth': 5, 'alpha': 0.6}
#            | line_segmean_kwargs
#        )
#        line_predict_kwargs = (
#            {'color': 'tab:red', 'linewidth': 1.5, 'alpha': 0.6}
#            | line_predict_kwargs
#        )
#        line_predict_clonal_kwargs = (
#            {'color': 'tab:orange', 'linewidth': 1.5, 'alpha': 0.6}
#            | line_predict_clonal_kwargs
#        )
#
#        # draw raw data
#        self.genomeplotter.draw_dots(
#            ax, 
#            y_colname='depthratio_raw', 
#            df_plotdata=raw_plotdata, 
#            plot_kwargs=dot_kwargs,
#        )
#
#        # draw segment information
#        if draw_segment_info:
#            if use_merged_segment:
#                segment_plotdata = self.data[sampleid]['merged_segment_plotdata']
#            else:
#                segment_plotdata = self.data[sampleid]['depthratio_segment_plotdata']
#
#            if draw_segmean:
#                self.genomeplotter.draw_hlines(
#                    ax, 
#                    df_plotdata=segment_plotdata,
#                    y_colname='depthratio_segment_mean', 
#                    plot_kwargs=line_segmean_kwargs,
#                )
#            if draw_predicted:
#                self.genomeplotter.draw_hlines(
#                    ax, 
#                    df_plotdata=segment_plotdata,
#                    y_colname='depthratio_predicted', 
#                    plot_kwargs=line_predict_kwargs,
#                )
#            if draw_depthratio_peaks:
#                segments_gr = pr.PyRanges(self.data[sampleid]['merged_segment'])
#                peak_xs = cnvmisc.get_depthratio_peaks(
#                    segments_gr.depthratio_segment_mean, 
#                    lengths=segments_gr.lengths(), 
#                    limits=(0, 2), 
#                    step=0.01, 
#                    peak_cutoff=1e8,
#                )
#            if draw_deviation:
#                self.genomeplotter.draw_boxes(
#                    ax,
#                    df_plotdata=segment_plotdata,
#                    ymins=(
#                        segment_plotdata['depthratio_segment_mean']
#                        - std_factor * segment_plotdata['depthratio_segment_std']
#                    ),
#                    ymaxs=(
#                        segment_plotdata['depthratio_segment_mean']
#                        + std_factor * segment_plotdata['depthratio_segment_std']
#                    ),
#                    colors='yellow',
#                    plot_kwargs=dict(alpha=0.4),
#                )
#
#            if mark_unfit_regions:
#                self.draw_unfit_region(ax, segment_plotdata)
#
#        # integer CN lines
#        if draw_integerCN_lines:
#            assert ((cellularity is not None) and (ploidy is not None))
#            self.draw_depthratio_ax_integer_CNs(ax, cellularity, ploidy, sampleid)
#
#        # modify axes
#        if modify_ax:
#            if ylabel is None:
#                ylabel = 'depth ratio'
#            ax.set_ylabel(ylabel, **ylabel_kwargs)
#
#            if ymax is None:
#                df = self.data[sampleid]['depthratio']
#                y_values = df['depthratio_raw'].loc[~df['excluded']].dropna()
#                ymax = y_values.quantile(0.999)
#
#                ax.set_ylim(-ymax * 0.1, ymax)
#                yticks = np.round(np.arange(0, ymax + 0.25, 0.25), 2)
#            else:
#                ax.set_ylim(-0.1, ymax)
#                yticks = np.round(np.linspace(0, ymax, 8), 2)
#
#            ax.set_yticks(yticks)
#            ax.set_yticklabels(yticks, size=self.get_yticklabel_size(yticks))
#
#            self.genomeplotter.draw_ax_common(
#                ax,
#                n_xlabel=n_xlabel, 
#                split_spines=split_spines,
#                merge_same_chroms=merge_same_chroms,
#                chromlabel_kwargs=chromlabel_kwargs,
#            )
#
#    def draw_depthratio_ax_integer_CNs(self, ax, cellularity, ploidy, sampleid):
#        plotregion_df_withCN = self.add_CNn_to_segment(
#            segment_df=self.genomeplotter.region_df,
#            mode=self.data[sampleid]['mode'],
#            refver=self.refver,
#            is_female=self.data[sampleid]['is_female'],
#            target_region=self.data[sampleid]['target_region'],
#        )
#        plotdata = self.genomeplotter.prepare_plot_data(plotregion_df_withCN)
#
#        integer_depthratios = cnvmisc.theoretical_depth_ratio(
#            CNt=np.arange(0, 10, 1)[:, np.newaxis],
#            cellularity=cellularity,
#            tumor_ploidy=ploidy,
#            CNn=plotregion_df_withCN['CNn'].to_numpy()[np.newaxis, :],
#            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
#        )  # ndim == 2
#        for ys in integer_depthratios:
#            plotdata['integer_depthratio'] = ys
#            self.genomeplotter.draw_hlines(
#                ax, 
#                y_colname='integer_depthratio',
#                df_plotdata=plotdata,
#                plot_kwargs=dict(linewidth=0.7, zorder=0, color='black', alpha=0.7),
#            )
#
#    def draw_unfit_region(self, ax, plotdata):
#        self.genomeplotter.draw_boxes(
#            ax,
#            df_plotdata=plotdata.loc[plotdata['polyploid_unfit'], :],
#            colors='yellow',
#            plot_kwargs=dict(alpha=0.2),
#        )
#        self.genomeplotter.draw_boxes(
#            ax,
#            df_plotdata=plotdata.loc[plotdata['polyploid_unfit_bafonly'], :],
#            colors='green',
#            plot_kwargs=dict(alpha=0.2),
#        )
#        self.genomeplotter.draw_boxes(
#            ax,
#            df_plotdata=plotdata.loc[plotdata['monoploid_unfit'], :],
#            colors='blue',
#            plot_kwargs=dict(alpha=0.2),
#        )
#
#    def draw_depthratio_hist_ax(
#        self,
#        sampleid,
#        ax,
#        use_merged_segment, 
#        depth_ylim,
#        #rm_haploid,
#        kde=True,
#        bw=None,
#        peak_threshold=None,
#        annotate_kwargs=dict(),
#        plot_kwargs=dict(),
#    ):
#        """Histogram only includes CNn == 2 positions"""
#        # handle kwargs
#        annotate_kwargs = (
#            dict(
#                va='center', ha='left', size=8,
#                arrowprops=dict(arrowstyle='-', color='black', linewidth=0.5),
#            )
#            | annotate_kwargs
#        )
#        plot_kwargs = (
#            #dict(linewidth=0.8)
#            dict(facecolor='tab:blue', alpha=1)
#            | plot_kwargs
#        )
#
#        if use_merged_segment:
#            segment_df = self.data[sampleid]['merged_segment']
#        else:
#            segment_df = self.data[sampleid]['depthratio_segment']
#
#        segment_df = segment_df.loc[segment_df['CNn'] == 2, :]
#
#        depthratio_list = segment_df['depthratio_segment_mean']
#        weights = (segment_df['End'] - segment_df.start0s)
#
#        # set ylim
#        ax.set_ylim(*depth_ylim)
#
#        if kde:
#            peak_values, peak_densities, density = cnvmisc.get_density_peaks(
#                depthratio_list, 
#                weights=weights, 
#                xs=None, 
#                threshold=peak_threshold, 
#                bw_method=bw, 
#                invert=False,
#            )
#            peak_depthratios = peak_values
#            ax_ymin, ax_ymax = ax.get_ylim()
#            fill_ys = np.linspace(ax_ymin, ax_ymax, 100)
#            fill_xs = density(fill_ys)
#        else:
#            histpeaks = cnvmisc.get_hist_peaks(
#                depthratio_list, 
#                weights=weights,
#                bins=np.arange(0, depthratio_list.max(), 0.01),
#                threshold=peak_threshold,
#            )
#            peak_depthratios = histpeaks['peak_values']
#            fill_ys = histpeaks['bin_midpoints']
#            fill_xs = histpeaks['hist']
#
#        # set xlim
#        ax.set_xlim(0, max(fill_xs))
#
#        # draw data
#        ax.fill_betweenx(y=fill_ys, x1=fill_xs, **plot_kwargs)
#        ytext_list = np.linspace(*ax.get_ylim(), len(peak_depthratios))
#        for y, ytext in zip(peak_depthratios, ytext_list):
#            ax.axhline(y, color='black', linewidth=0.5)
#            ax.annotate(
#                round(y, 3), 
#                (ax.get_xlim()[1], y),
#                xytext=(ax.get_xlim()[1] * 1.1, ytext),
#                annotation_clip=False,
#                **annotate_kwargs,
#            )
#        ax.set_yticks([])
#
#        return peak_depthratios
#
#    def draw_baf_ax(
#        self, 
#        sampleid, 
#        ax, 
#
#        use_merged_segment=True, 
#
#        draw_corrected=False,
#        draw_predicted=False,
#        draw_segmean=False,
#
#        n_xlabel=None,
#
#        is_tumor=True,
#
#        mark_unfit_regions=False,
#
#        dot_kwargs=dict(),
#        line_segmean_kwargs=dict(),
#        line_corr_segmean_kwargs=dict(),
#        line_predict_kwargs=dict(),
#        line_predict_clonal_kwargs=dict(),
#
#        split_spines=True,
#
#        ylabel=None,
#        ylabel_kwargs=dict(),
#
#        make_raw_plotdata=True,
#        make_segment_plotdata=True,
#
#        modify_ax=True,
#        merge_same_chroms=True,
#
#        chromlabel_kwargs=dict(),
#    ):
#        # prepare plotdata
#        if make_raw_plotdata:
#            if is_tumor:
#                self.make_tumor_baf_plotdata(sampleid)
#            else:
#                self.make_normal_baf_plotdata(sampleid)
#        if make_segment_plotdata:
#            self.make_segment_plotdata(sampleid)
#
#        # find plotdata
#        if is_tumor:
#            raw_plotdata = self.data[sampleid]['baf_raw_plotdata']
#        else:
#            raw_plotdata = self.data[sampleid]['normal_baf_raw_plotdata']
#
#        # calc default alpha
#        n_dots = raw_plotdata.shape[0]
#        default_alpha = calc_dot_alpha_baf(n_dots)
#
#        # handle kwargs
#        dot_kwargs = (
#            #{'color': 'black', 'markersize': 0.3, 'alpha': 0.01}
#            {'color': 'black', 'markersize': 0.3, 'alpha': default_alpha}
#            | dot_kwargs
#        )
#        line_segmean_kwargs = (
#            {'color': 'tab:blue', 'linewidth': 2, 'alpha': 0.4}
#            | line_segmean_kwargs
#        )
#        line_corr_segmean_kwargs = (
#            {'color': 'tab:green', 'linewidth': 2, 'alpha': 0.4}
#            | line_corr_segmean_kwargs
#        )
#        line_predict_kwargs = (
#            {'color': 'tab:red', 'linewidth': 2, 'alpha': 0.4}
#            | line_predict_kwargs
#        )
#        line_predict_clonal_kwargs = (
#            {'color': 'tab:orange', 'linewidth': 2, 'alpha': 0.4}
#            | line_predict_clonal_kwargs
#        )
#
#        # draw raw data
#        if is_tumor:
#            self.genomeplotter.draw_dots(
#                ax, 
#                df_plotdata=raw_plotdata, 
#                y_colname='baf_raw_tumor', 
#                plot_kwargs=dot_kwargs,
#            )
#        else:
#            self.genomeplotter.draw_dots(
#                ax, 
#                df_plotdata=raw_plotdata, 
#                y_colname='baf_raw_normal', 
#                plot_kwargs=dot_kwargs,
#            )
#
#        # draw segment information
#        draw_segment_info = any([
#            draw_corrected,
#            draw_predicted,
#            draw_segmean,
#        ])
#
#        if draw_segment_info:
#            if is_tumor:
#                segment_plotdata = self.data[sampleid]['merged_segment_plotdata']
#            else:
#                segment_plotdata = self.data[sampleid]['normal_baf_segment_plotdata']
#
#            if draw_segmean:
#                self.genomeplotter.draw_hlines(
#                    ax, 
#                    df_plotdata=segment_plotdata,
#                    y_colname='baf_segment_peak', 
#                    #y_colname='baf_segment_mean', 
#                    plot_kwargs=line_segmean_kwargs,
#                )
#            if draw_corrected:
#                self.genomeplotter.draw_hlines(
#                    ax, 
#                    df_plotdata=segment_plotdata,
#                    y_colname='corrected_baf_segment_mean', 
#                    plot_kwargs=line_corr_segmean_kwargs,
#                )
#
#            if draw_predicted:
#                self.genomeplotter.draw_hlines(
#                    ax, 
#                    df_plotdata=segment_plotdata,
#                    y_colname='baf_predicted', 
#                    plot_kwargs=line_predict_kwargs,
#                )
#
#            if mark_unfit_regions:
#                self.draw_unfit_region(ax, segment_plotdata)
#
#        # set axes attributes
#        if ylabel is None:
#            ylabel = 'tumor baf' if is_tumor else 'normal baf'
#        ax.set_ylabel(ylabel, **ylabel_kwargs)
#        ax.set_ylim(-0.6 * 0.1, 0.6)
#
#        yticks = np.round(np.arange(0, 0.6, 0.1), 1)
#        ax.set_yticks(yticks)
#        ax.set_yticklabels(yticks, size=self.get_yticklabel_size(yticks))
#
#        if modify_ax:
#            self.genomeplotter.draw_ax_common(
#                ax, 
#                n_xlabel=n_xlabel, 
#                split_spines=split_spines, 
#                merge_same_chroms=merge_same_chroms,
#                chromlabel_kwargs=chromlabel_kwargs,
#            )
#
#    def draw_CN_ax(
#        self, 
#        sampleid, 
#        ax,
#        n_xlabel=None,
#        line_A_kwargs=dict(),
#        line_B_kwargs=dict(),
#        line_CNt_kwargs=dict(),
#        ymax=None,
#        draw_CNt=True,
#        draw_A=True,
#        draw_B=True,
#
#        split_spines=True,
#
#        ylabel=None,
#        ylabel_kwargs=dict(),
#
#        merge_same_chroms=True,
#        chromlabel_kwargs=dict(),
#    ):
#        # handle kwargs
#        line_CNt_kwargs = (
#            {'color': 'black'}
#            | line_CNt_kwargs
#        )
#        line_A_kwargs = (
#            {'color': 'red'}
#            | line_A_kwargs
#        )
#        line_B_kwargs = (
#            {'color': 'blue'}
#            | line_B_kwargs
#        )
#
#        # draw data
#        if draw_CNt:
#            self.genomeplotter.draw_hlines(
#                ax, 
#                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
#                y_colname='clonal_CNt', 
#                offset=0.1,
#                plot_kwargs=line_CNt_kwargs,
#            )
#        if draw_A:
#            self.genomeplotter.draw_hlines(
#                ax, 
#                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
#                y_colname='clonal_A', 
#                offset=0,
#                plot_kwargs=line_A_kwargs,
#            )
#        if draw_B:
#            self.genomeplotter.draw_hlines(
#                ax, 
#                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
#                y_colname='clonal_B', 
#                offset=-0.1,
#                plot_kwargs=line_B_kwargs,
#            )
#
#        # set axes attributes
#        if ylabel is None:
#            ylabel = 'tumor clonal copy number'
#        ax.set_ylabel(ylabel, **ylabel_kwargs)
#
#        if ymax is None:
#            df = self.data[sampleid]['merged_segment']
#            ymax = df['clonal_CNt'].quantile(0.99)
#        else:
#            pass
#
#        max_ticknum = 15
#        step = np.ceil(ymax / 15).astype(int)
#        yticks = np.arange(0, ymax, step).astype(int)
#
#        ax.set_ylim(-0.5, ymax)
#        ax.set_yticks(yticks, minor=False)
#        #ax.set_yticks(np.arange(0, ymax, 1).astype(int), minor=True)
#        ax.set_yticklabels(yticks, size=self.get_yticklabel_size(yticks))
#
#        self.genomeplotter.draw_ax_common(
#            ax, 
#            n_xlabel=n_xlabel, 
#            split_spines=split_spines, 
#            merge_same_chroms=merge_same_chroms,
#            chromlabel_kwargs=chromlabel_kwargs,
#        )
#
#    def draw_subclonal_CN_ax(
#        self, 
#        sampleid, 
#        ax,
#        n_xlabel=None,
#        line_A_kwargs=dict(),
#        line_B_kwargs=dict(),
#        line_CNt_kwargs=dict(),
#        ymax=None,
#        draw_CNt=True,
#        draw_A=True,
#        draw_B=True,
#
#        split_spines=True,
#
#        ylabel=None,
#        ylabel_kwargs=dict(),
#
#        merge_same_chroms=True,
#        chromlabel_kwargs=dict(),
#    ):
#        # handle kwargs
#        line_CNt_kwargs = (
#            {'color': 'black'}
#            | line_CNt_kwargs
#        )
#        line_A_kwargs = (
#            {'color': 'red'}
#            | line_A_kwargs
#        )
#        line_B_kwargs = (
#            {'color': 'blue'}
#            | line_B_kwargs
#        )
#
#        # draw CNt
#        if draw_CNt:
#            self.genomeplotter.draw_hlines(
#                ax, 
#                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
#                y_colname='subclonal_CNt', 
#                offset=0.1,
#                plot_kwargs=line_CNt_kwargs,
#            )
#        if draw_A:
#            self.genomeplotter.draw_hlines(
#                ax, 
#                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
#                y_colname='subclonal_A', 
#                offset=0,
#                plot_kwargs=line_A_kwargs,
#            )
#        if draw_B:
#            self.genomeplotter.draw_hlines(
#                ax, 
#                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
#                y_colname='subclonal_B', 
#                offset=-0.1,
#                plot_kwargs=line_B_kwargs,
#            )
#
#        # set axes attributes
#        if ylabel is None:
#            ylabel = 'tumor subclonal copy number'
#        ax.set_ylabel(ylabel, **ylabel_kwargs)
#
#        if ymax is None:
#            df = self.data[sampleid]['merged_segment']
#            ymax = df['subclonal_CNt'].quantile(0.99)
#        else:
#            pass
#
#        max_ticknum = 15
#        step = np.ceil(ymax / 15).astype(int)
#        yticks = np.arange(0, ymax, step).astype(int)
#
#        ax.set_ylim(-0.5, ymax)
#        ax.set_yticks(yticks, minor=False)
#        #ax.set_yticks(np.arange(0, ymax, 1).astype(int), minor=True)
#        ax.set_yticklabels(yticks, size=self.get_yticklabel_size(yticks))
#
#        self.genomeplotter.draw_ax_common(
#            ax, 
#            n_xlabel=n_xlabel, 
#            split_spines=split_spines, 
#            merge_same_chroms=merge_same_chroms,
#            chromlabel_kwargs=chromlabel_kwargs,
#        )
#
#    def draw_ccf_ax(
#        self, 
#        sampleid, 
#        ax,
#        n_xlabel=None,
#        bar_kwargs=dict(),
#
#        split_spines=True,
#
#        ylabel=None,
#        ylabel_kwargs=dict(),
#
#        merge_same_chroms=True,
#        chromlabel_kwargs=dict(),
#
#        mark_clonal_region=False,
#    ):
#        bar_kwargs = (
#            dict()
#            | bar_kwargs
#        )
#
#        plotdata = self.data[sampleid]['merged_segment_plotdata']
#
#        self.genomeplotter.draw_bars(
#            ax, 
#            y_colname='ccf', 
#            df_plotdata=plotdata,
#            plot_kwargs=bar_kwargs,
#        )
#
#        if ylabel is None:
#            ylabel = 'subclonal fraction'
#        ax.set_ylabel(ylabel, **ylabel_kwargs)
#
#        ax.set_ylim(-0.05, 1.05)
#        yticks = np.round(np.arange(0, 1.1, 0.1), 1)
#        ax.set_yticks(yticks)
#        ax.set_yticklabels(yticks)
#
#        if mark_clonal_region:
#            self.genomeplotter.draw_boxes(
#                ax,
#                df_plotdata=plotdata.loc[plotdata['ccf'].isna(), :],
#                colors='green',
#                plot_kwargs=dict(alpha=0.1),
#            )
#
#        self.genomeplotter.draw_ax_common(
#            ax, 
#            n_xlabel=n_xlabel, 
#            split_spines=split_spines, 
#            merge_same_chroms=merge_same_chroms,
#            chromlabel_kwargs=chromlabel_kwargs,
#        )
#
#    def show_ccfs_main(
#        self,
#        ccfs,
#        lengths,
#        density,
#        peak_values,
#        draw_dots=False,
#    ):
#        fig, ax = plt.subplots(figsize=(30, 5))
#        ax.hist(
#            ccfs, 
#            range=(0, 1), 
#            bins=50, 
#            weights=lengths,
#            density=True,
#        )
#
#        density_xs = np.arange(0, 1, 0.01)
#        density_ys = density(density_xs)
#        ax.plot(density_xs, density_ys)
#
#        for x in peak_values:
#            ax.axvline(x, color='red')
#
#        if draw_dots:
#            ylim = ax.get_ylim()
#            dot_xs = ccfs
#            dot_ys = scipy.stats.uniform.rvs(loc=ylim[0], scale=(ylim[1] - ylim[0]), size=len(ccfs))
#            colors = fit.labels_
#            ax.scatter(dot_xs, dot_ys, c=colors, alpha=0.7, s=4)
#
#        return fig, ax
#
#    #############################
#    # helpers of add new sample #
#    #############################
#
#    def set_target_region(self, sampleid, mode, target_region_arg):
#        if mode == 'wgs':
#            if 'normal_depth' in self.data[sampleid]:
#                target_region_gr = self.find_germline_copyneutral_region(
#                    sampleid, factors=(0.8, 1.2),
#                )
#            else:
#                target_region_gr = refgenome.get_chromdict(self.refver).to_gr(
#                    assembled_only=True, as_gr=True,
#                )
#        elif mode == 'panel':
#            target_region_gr = cnvmisc.arg_into_gr(target_region_arg).drop()
#
#        # exclude y if female
#        if self.data[sampleid]['is_female']:
#            y_chrom = refgenome.get_chromdict(self.refver).XY_names[1]
#            target_region_gr = target_region_gr[
#                target_region_gr.Chromosome != y_chrom
#            ]
#
#        # merge
#        target_region_gr = target_region_gr.merge()
#
#        # result
#        self.data[sampleid]['target_region'] = target_region_gr
#
#    def find_germline_copyneutral_region(self, sampleid, factors=(0.3, 2)):
#        depth_df = self.data[sampleid]['normal_depth'].copy()
#
#        # make correction for male haploid chromosomes
#        if not self.data[sampleid]['is_female']:
#            XY = refgenome.get_chromdict(self.refver).XY_names
#            selector = depth_df['Chromosome'].isin(XY)
#            depth_df.loc[selector, 'mean_depth'] = depth_df.loc[selector, 'mean_depth'] * 2
#
#        # upscale depth bins for segmentation speed
#        upscaled_depth_df = cnvmisc.upsize_depth_df_bin(
#            depth_df, size=1000, refver=self.refver,
#        )
#        upscaled_depth_df['depth_raw'] = upscaled_depth_df['mean_depth']
#
#        # run segmentation
#        LOGGER_INFO.info('Running segmentation of normal bam depth in order to find copy-neutral region')
#        segdf, _ = rcopynumber.run_rcopynumber(
#            depth_df=upscaled_depth_df,
#            refver=self.refver, 
#            as_gr=False, 
#            winsorize=False, 
#            verbose=False,
#            remove_unassembled_contigs=True,
#        )
#        self.data[sampleid]['normal_depth_segment'] = segdf
#
#        # select copy-neutral depth segments
#        global_mean = np.average(
#            upscaled_depth_df['mean_depth'], 
#            weights=(upscaled_depth_df['End'] - upscaled_depth_df.start0s),
#        )
#        selector = segdf['depth_segment_mean'].between(
#            global_mean * factors[0], 
#            global_mean * factors[1], 
#        )
#        included_segments = segdf.loc[selector, :]
#
#        target_region_gr = pr.PyRanges(depth_df).drop().overlap(
#            pr.PyRanges(included_segments)
#        ).merge()
#
#        return target_region_gr
#
#    ###############################
#    # helpers of segment creation #
#    ###############################
#
##    def make_depth_segment(
##        self,
##        sampleid,
##        winsorize=False,
##        gamma=None,
##        kmin=None,
##        verbose=True,
##    ):
##        self.data[sampleid]['depthratio_segment'] = _make_depth_segment(
##            self.data[sampleid]['depthratio_upscaled'],
##            mode=self.data[sampleid]['mode'],
##            refver=self.refver,
##            winsorize=winsorize,
##            gamma=gamma,
##            kmin=kmin,
##            verbose=verbose,
##        )
##
##    def make_baf_segment(
##        self,
##        sampleid,
##        winsorize=False,
##        gamma=None,
##        kmin=None,
##        verbose=False,
##        bafcorrector=libbaf.load_bafcorrect_func(),
##    ):
##        self.data[sampleid]['baf_segment'] = _make_baf_segment(
##            baf_df=self.data[sampleid]['baf'],
##            target_region=self.data[sampleid]['target_region'],
##            mode=self.data[sampleid]['mode'],
##            refver=self.refver,
##            winsorize=winsorize,
##            gamma=gamma,
##            kmin=kmin,
##            verbose=verbose,
##            bafcorrector=bafcorrector,
##        )
#
#    def make_merged_segment(self, sampleid, depthratio_segment):
#        merged_segment = pyranges_helper.isec_union(
#            depthratio_segment,
#            self.data[sampleid]['baf_segment'],
#        )
#        merged_segment = pyranges_helper.join(
#            merged_segment, 
#            depthratio_segment,
#            how='left',
#            merge=None,
#            find_nearest=True,
#            sort=True,
#            refver=self.refver,
#        )
#        merged_segment = pyranges_helper.join(
#            merged_segment, 
#            self.data[sampleid]['baf_segment'],
#            how='left',
#            merge=None,
#            find_nearest=True,
#            sort=True,
#            refver=self.refver,
#        )
#
#        # reduce by merging adjacent segments with identical annotation values
#        merged_segment = self.deduplicate_merged_segment(merged_segment)
#
#        # add CNn
#        merged_segment = self.add_CNn_to_segment(
#            merged_segment,
#            self.data[sampleid]['mode'],
#            self.refver,
#            self.data[sampleid]['is_female'],
#            self.data[sampleid]['target_region'],
#        )
#
#        # add std and std/mean ratio
#        #self.add_depthratio_std_to_segment(sampleid)
#
#        # fit to target region
#        #merged_segment = self.fit_segment_to_targetregion(sampleid, merged_segment)
#
#        return merged_segment
#
#        #self.data[sampleid]['merged_segment'] = merged_segment
#
#    def fit_segment_to_targetregion(self, sampleid, segment_df):
#        segment_gr = cnvmisc.arg_into_gr(segment_df)
#        seg_subset = segment_gr.intersect(self.data[sampleid]['target_region'])
#        target_diff_seg = self.data[sampleid]['target_region'].subtract(segment_gr)
#        target_diff_seg = pyranges_helper.join(
#            target_diff_seg, 
#            segment_gr,
#            how='left',
#            merge=None,
#            find_nearest=True,
#            as_gr=True,
#        )
#        result = pr.concat([seg_subset, target_diff_seg]).df
#        result = cnvmisc.sort_genome_df(result, self.refver)
#
#        return result
#
#    def deduplicate_merged_segment(self, merged_segment):
#        merged_segment = cnvmisc.sort_genome_df(merged_segment, self.refver)
#
#        chromdict = refgenome.get_chromdict(self.refver)
#        annot_arr = np.concatenate(
#            [
#                merged_segment['Chromosome'].apply(chromdict.contigs.index).to_numpy()[:, np.newaxis],
#                merged_segment.loc[:, ['depthratio_segment_mean', 'baf_segment_mean']].to_numpy(),
#            ],
#            axis=1,
#        )
#        values, counts, groupkey = tools.array_grouper(annot_arr, omit_values=True)
#        groupby = merged_segment.groupby(groupkey)
#        dedup_segdf = groupby.first()
#        dedup_segdf['End'] = groupby['End'].last()
#
#        return dedup_segdf
#
##    def postprocess_segment(self, sampleid, cellularity, ploidy):
##        # add CNt and B
##        cpinfo = self.data[sampleid]['cpscores'][(cellularity, ploidy)]
##        self.data[sampleid]['merged_segment']['CNt'] = cpinfo['CNt_list']
##        self.data[sampleid]['merged_segment']['B'] = cpinfo['B_list']
##        self.data[sampleid]['merged_segment']['A'] = (
##            self.data[sampleid]['merged_segment']['CNt']
##            - self.data[sampleid]['merged_segment']['B']
##        )
##
##        # add theoreticals
##        self.data[sampleid]['merged_segment'] = cnvmisc.add_theoreticals_to_segment(
##            self.data[sampleid]['merged_segment'], 
##            cellularity=cellularity, 
##            tumor_ploidy=ploidy, 
##            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'], 
##            is_female=self.data[sampleid]['is_female'],
##        )
#
#    ####################
#    # solution finding #
#    ####################
#
#    def make_CN_solution_freeccf(
#        self,
#        sampleid,
#        cellularity,
#        ploidy,
#        depth_ratio_diff=None,
#        baf_diff=None,
#        min_N_CNt_candidates=5,
#        N_CNt_candidates_fraction=0.5,
#        limited_clonal=True,
#    ):
#        segdf = self.data[sampleid]['merged_segment']
#        (
#            clonal_solution, 
#            flags, 
#            freeccf_solution,
#            calculated_depth_ratio, 
#            calculated_baf,
#            average_CNt,
#        ) = cnvmisc.find_solution_before_ccfs(
#            depth_ratio=segdf['depthratio_segment_mean'], 
#            baf=segdf['corrected_baf_segment_mean'],
#            CNn=segdf['CNn'],
#            lengths=(segdf['End'] - segdf.start0s),
#            cellularity=cellularity,
#            tumor_ploidy=ploidy,
#            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
#            depth_ratio_diff=depth_ratio_diff,
#            baf_diff=baf_diff,
#            Bn=1,
#            min_N_CNt_candidates=min_N_CNt_candidates,
#            N_CNt_candidates_fraction=N_CNt_candidates_fraction,
#
#            limited_clonal=limited_clonal,
#        )
#
#        freeccf_result = {
#            #'fixed_ccfs': fixed_ccfs, 
#            #'ccf_plotdata': ccf_plotdata, 
#            'clonal_solution': clonal_solution, 
#            'flags': flags, 
#            'freeccf_solution': freeccf_solution,
#            'calculated_depth_ratio': calculated_depth_ratio, 
#            'calculated_baf': calculated_baf,
#        }
#        self.data[sampleid]['freeccf_result'] = freeccf_result
#        self.data[sampleid]['average_CNt'] = average_CNt
#
#    def select_fixed_ccfs(self, sampleid, bandwidth=0.1):
#        segdf = self.data[sampleid]['merged_segment']
#        lengths = (segdf['End'] - segdf.start0s).to_numpy()
#        fixed_ccfs, ccf_plotdata = cnvmisc.select_fixed_ccfs(
#            freeccf_solution=self.data[sampleid]['freeccf_result']['freeccf_solution'], 
#            lengths=lengths, 
#            flags=self.data[sampleid]['freeccf_result']['flags'], 
#            bandwidth=bandwidth,
#        )
#        self.data[sampleid]['fixed_ccfs'] = fixed_ccfs
#        self.data[sampleid]['ccf_plotdata'] = ccf_plotdata
#
#    def make_CN_solution_after_ccfs(
#        self,
#        sampleid,
#        cellularity,
#        ploidy,
#        min_N_CNt_candidates=5,
#        N_CNt_candidates_fraction=0.5,
#        CNt_diff_factor=0.1,
#        limited_clonal=True,
#    ):
#        segdf = self.data[sampleid]['merged_segment']
#        solution = cnvmisc.find_solution_after_ccfs(
#            depth_ratio=segdf['depthratio_segment_mean'].to_numpy(),
#            baf=segdf['corrected_baf_segment_mean'].to_numpy(),
#            CNn=segdf['CNn'].to_numpy(),
#            lengths=(segdf['End'] - segdf.start0s).to_numpy(),
#            cellularity=cellularity,
#            tumor_ploidy=ploidy,
#            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
#            average_CNt=self.data[sampleid]['average_CNt'],
#
#            fixed_ccfs=self.data[sampleid]['fixed_ccfs'], 
#            clonal_solution=self.data[sampleid]['freeccf_result']['clonal_solution'], 
#            flags=self.data[sampleid]['freeccf_result']['flags'],
#
#            Bn=1,
#            min_N_CNt_candidates=min_N_CNt_candidates,
#            N_CNt_candidates_fraction=N_CNt_candidates_fraction,
#
#            CNt_diff_factor=CNt_diff_factor,
#
#            limited_clonal=limited_clonal,
#        )
#
#        self.data[sampleid]['CN_solution'] = solution
#
##    def make_CN_solution_onestep(
##        self,
##        sampleid,
##        cellularity,
##        ploidy,
##        depth_ratio_diff=None,
##        baf_diff=0.05,
##        min_N_CNt_candidates=5,
##        N_CNt_candidates_fraction=0.5,
##        ccf_bw=0.1,
##    ):
##        self.make_CN_solution_freeccf(
##            sampleid,
##            cellularity,
##            ploidy,
##            depth_ratio_diff=depth_ratio_diff,
##            baf_diff=baf_diff,
##            min_N_CNt_candidates=min_N_CNt_candidates,
##            N_CNt_candidates_fraction=N_CNt_candidates_fraction,
##        )
##        self.select_fixed_ccfs(sampleid, bandwidth=ccf_bw)
##        self.make_CN_solution_after_ccfs(
##            sampleid,
##            cellularity,
##            ploidy,
##            min_N_CNt_candidates=min_N_CNt_candidates,
##            N_CNt_candidates_fraction=N_CNt_candidates_fraction,
##        )
#
#    def add_freeccf_solution_to_segment(self, sampleid):
#        #subclonal_theoretical_depthratio = self.data[sampleid]['freeccf_result']['calculated_depth_ratio']
#        #subclonal_theoretical_baf = self.data[sampleid]['freeccf_result']['calculated_baf']
#        cnvmisc.add_freeccf_solution_to_segment(
#            segment_df=self.data[sampleid]['merged_segment'], 
#            clonal_solution=self.data[sampleid]['freeccf_result']['clonal_solution'],
#            freeccf_subclonal_solution=self.data[sampleid]['freeccf_result']['freeccf_solution'], 
#            flags=self.data[sampleid]['freeccf_result']['flags'], 
#            unfit_region_depth_ratio=self.data[sampleid]['freeccf_result']['calculated_depth_ratio'], 
#            unfit_region_baf=self.data[sampleid]['freeccf_result']['calculated_baf'], 
#        )
#
#    def add_fixedccf_solution_to_segment(self, sampleid):
#        cnvmisc.add_fixedccf_solution_to_segment(
#            segment_df=self.data[sampleid]['merged_segment'], 
#            fixedccf_solution=self.data[sampleid]['CN_solution'],
#        )
#
#    def add_solution_to_plotdata(self, sampleid):
#        assert (
#            self.data[sampleid]['merged_segment'].loc[:, ['Chromosome', 'Start', 'End']]
#            == self.data[sampleid]['merged_segment_plotdata'].loc[:, ['Chromosome', 'Start', 'End']]
#        ).all(axis=None)
#
#        cnvmisc.add_fixedccf_solution_to_segment(
#            segment_df=self.data[sampleid]['merged_segment_plotdata'], 
#            fixedccf_solution=self.data[sampleid]['CN_solution'],
#        )
#
#    def add_CN_pred_to_segment(self, sampleid, cellularity, ploidy):
#        if 'cpscores' not in self.data[sampleid]:
#            self.data[sampleid]['cpscores'] = dict()
#
#        # add CNt and B
#        cpinfo = cnvmisc.calc_cp_score(
#            segment_df=self.data[sampleid]['merged_segment'], 
#            cellularity=cellularity, 
#            ploidy=ploidy, 
#            is_female=self.data[sampleid]['is_female'], 
#            CNt_weight=1, 
#            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
#        )
#        self.data[sampleid]['cpscores'][(cellularity, ploidy)] = cpinfo
#
#        self.data[sampleid]['merged_segment']['CNt'] = cpinfo['CNt_list']
#        self.data[sampleid]['merged_segment']['B'] = cpinfo['B_list']
#        self.data[sampleid]['merged_segment']['A'] = (
#            self.data[sampleid]['merged_segment']['CNt']
#            - self.data[sampleid]['merged_segment']['B']
#        )
#
#        # add theoreticals
#        self.data[sampleid]['merged_segment'] = cnvmisc.add_theoreticals_to_segment(
#            self.data[sampleid]['merged_segment'], 
#            cellularity=cellularity, 
#            tumor_ploidy=ploidy, 
#            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'], 
#            is_female=self.data[sampleid]['is_female'],
#        )
#
#    def add_CN_pred_to_segment_new(self, sampleid, cellularity, ploidy):
#        # calculate CNt values
#        depth_ratios = self.data[sampleid]['merged_segment']['depthratio_segment_mean'].to_numpy()
#        CNns = self.data[sampleid]['merged_segment']['CNn'].to_numpy()
#        clonal_CNts = cnvmisc.calc_clonal_CNt(
#            depth_ratio=depth_ratios,
#            CNn=CNns,
#
#            cellularity=cellularity,
#            tumor_ploidy=ploidy,
#            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
#        )
#
#        subclonal_CNts, ccfs = cnvmisc.calc_subclonal_CNt(
#            depth_ratio=depth_ratios,
#            clonal_CNt=clonal_CNts,
#            CNn=CNns,
#
#            cellularity=cellularity,
#            tumor_ploidy=ploidy,
#            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
#            tumor_avg_depth_ratio=1,
#            normal_avg_depth_ratio=1,
#            only_max_ccf=True,
#        )
#
#        # calculate theoretical depths
#        predicted_depth_ratios = cnvmisc.theoretical_depth_ratio_subclone(
#            clonal_CNt=clonal_CNts, 
#            subclonal_CNt=subclonal_CNts,
#            ccf=ccfs,
#            cellularity=cellularity, 
#            tumor_ploidy=ploidy, 
#            CNn=CNns, 
#            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
#            tumor_avg_depth_ratio=1, 
#            normal_avg_depth_ratio=1,
#        )
#
#        predicted_depth_ratios_clonal = cnvmisc.theoretical_depth_ratio(
#            CNt=clonal_CNts, 
#            cellularity=cellularity, 
#            tumor_ploidy=ploidy, 
#            CNn=CNns, 
#            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
#            tumor_avg_depth_ratio=1, 
#            normal_avg_depth_ratio=1,
#        )
#
#        # annotate segment dataframe
#        merged_segdf = self.data[sampleid]['merged_segment']
#
#        merged_segdf['clonal_CNt'] = clonal_CNts
#        merged_segdf['subclonal_CNt'] = subclonal_CNts
#        merged_segdf['ccf'] = ccfs
#        merged_segdf['depthratio_predicted'] = predicted_depth_ratios
#        merged_segdf['depthratio_predicted_clonal'] = predicted_depth_ratios_clonal
#
#        self.data[sampleid]['merged_segment'] = merged_segdf
#
#        # add standard deviation of depth ratios
#        self.add_depthratio_std_to_segment(sampleid)
#
#    def add_depthratio_std_to_segment(self, sampleid):
#        self.data[sampleid]['merged_segment'] = add_depthratio_std(
#            self.data[sampleid]['merged_segment'], 
#            self.data[sampleid]['depthratio_upscaled'],
#            self.refver,
#        )
#
#    def add_baf_std_to_segment(self, sampleid):
#        merged_segdf = self.data[sampleid]['merged_segment']
#        right_df = self.data[sampleid]['original_baf']
#        right_df = right_df.loc[
#            right_df['baf_raw_tumor'] > 0, 
#            ['Chromosome', 'Start', 'End', 'baf_raw_tumor'],
#        ]
#        joined_segdf = pyranges_helper.join(
#            merged_segdf,
#            right_df,
#            how='left', merge='mean', add_std=True, ddof=0,
#            sort=True, refver=self.refver,
#        )
#        joined_segdf.drop('baf_raw_tumor', axis=1, inplace=True)
#        joined_segdf.rename(
#            columns={'baf_raw_tumor_std': 'baf_segment_std'}, 
#            inplace=True,
#        )
#
#        self.data[sampleid]['merged_segment'] = joined_segdf
#
#    ####################
#    # cp value finding #
#    ####################
#
#    def get_valid_CNts(
#        self, sampleid, depthratio1, depthratio2, CNt1_range=range(0, 6), CNt2_maxdiff=6,
#    ):
#        return cnvmisc.get_valid_CNts_from_depthratios(
#            depthratio1, 
#            depthratio2, 
#            normal_mean_ploidy=self.data[sampleid]['normal_mean_ploidy'], 
#            CNn=2, 
#            tumor_avg_depth_ratio=1, 
#            normal_avg_depth_ratio=1, 
#            CNt1_range=CNt1_range, 
#            CNt2_maxdiff=CNt2_maxdiff,
#        )
#
#    def calc_cpscore(
#        self, 
#        sampleid, 
#        CNt_weight=cnvmisc.DEFAULT_CNT_WEIGHT,
#        nproc=1,
#    ):
#        cpscore_dict = cnvmisc.get_cp_score_dict(
#            self.data[sampleid]['merged_segment'], 
#            refver=self.refver, 
#            is_female=self.data[sampleid]['is_female'], 
#            target_region_gr=self.data[sampleid]['target_region'], 
#            CNt_weight=CNt_weight, 
#            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
#            nproc=nproc,
#        )
#        peak_values, dfs = cnvmisc.get_peak_info(cpscore_dict)
#
#        self.data[sampleid]['cpscores'] = cpscore_dict
#        self.data[sampleid]['peak_values'] = peak_values
#        self.data[sampleid]['peak_dfs'] = dfs
#
#    def show_peaks(self, sampleid, figsize=(20, 20), **kwargs):
#        fig, ax = cnvmisc.show_heatmap_peaks_new(
#            self.data[sampleid]['peak_dfs'], figsize=figsize, **kwargs,
#        )
#
#    #########
#    # depth #
#    #########
#
#    def load_bam_for_depth(self, sampleid, bam_path, sampletype):
#        assert sampletype in ('normal', 'tumor')
#        mosdepth_df = self._run_mosdepth(
#            bam_path,
#            binsize=self.default_binsize,
#            region_df=(
#                None 
#                if self.data[sampleid]['mode'] == 'wgs' else 
#                self.data[sampleid]['target_region']
#            ),
#        )
#
#        self.data.setdefault(sampleid, dict())
#        self.data[sampleid][f'{sampletype}_depth'] = mosdepth_df
#
#    def load_mosdepth_df(
#        self, sampleid, mosdepth_df, sampletype, use_default_gc=True, **kwargs,
#    ):
#        assert sampletype in ('normal', 'tumor')
#
#        self.data.setdefault(sampleid, dict())
#        depth_df, gcbin_average_depths = self._postprocess_mosdepth_df(
#            mosdepth_df, 
#            self.data[sampleid]['mode'],
#            sampletype, 
#            use_default_gc=use_default_gc, 
#            **kwargs,
#        )
#        self.data[sampleid][f'{sampletype}_depth'] = depth_df
#        self.data[sampleid][f'{sampletype}_gcdata'] = gcbin_average_depths
#
#    def _run_mosdepth(self, bam_path, binsize, region_df):
#        mosdepth_df, _ = libmosdepth.run_mosdepth(
#            bam_path, 
#            t=8, 
#            use_median=False, 
#            region_bed_path=None, 
#            region_gr=region_df, 
#            window_size=binsize, 
#            donot_subset_bam=True,
#            as_gr=False, 
#            load_perbase=False,
#        )
#
#        return mosdepth_df
#
#    def _postprocess_mosdepth_df(
#        self, 
#        mosdepth_df, 
#        mode,
#        sampletype, 
#        use_default_gc=True, 
#        **kwargs,
#    ):
#        assert sampletype in ('normal', 'tumor')
#
#        # set preset_cutoffs
#        #if mode == 'panel':
#            #preset_cutoffs = 'panel'
#        #elif mode == 'wgs':
#            #if sampletype == 'normal':
#                #preset_cutoffs = 'normal_wgs'
#            #elif sampletype == 'tumor':
#                #preset_cutoffs = 'wgs'
#        #kwargs['preset_cutoffs'] = preset_cutoffs
#
#        # set gc_df
#        if use_default_gc:
#            gc_df = libgcfraction.get_gc_df(
#                self.refver, self.default_binsize, coords_as_index=True,
#            )
#            kwargs['gc_df'] = gc_df
#
#        kwargs['as_gr'] = False
#        depth_df, gcbin_average_depths = cnvmisc.postprocess_depth_df(
#            mosdepth_df, 
#            **kwargs,
#        )
#        return depth_df, gcbin_average_depths
#
#    #######
#    # baf #
#    #######
#
#    def load_germline_vcf(
#        self, 
#        sampleid, 
#        vcf_path, 
#        vcf_sampleid_tumor,
#        vcf_sampleid_normal=None,
#        nproc=1, 
#        logging_lineno=50000,
#    ):
#        self.data.setdefault(sampleid, dict())
#
#        if vcf_sampleid_normal is None:
#            vcf_sampleids = [vcf_sampleid_tumor]
#        else:
#            vcf_sampleids = [vcf_sampleid_tumor, vcf_sampleid_normal]
#
#        # load vafdf
#        vaf_df = variantplus.get_vafdf(
#            vcf_path, 
#            sampleid=vcf_sampleids, 
#            nproc=nproc,
#        )
#        # rename vaf columns
#        vaf_df.rename(
#            columns={f'vaf_{vcf_sampleid_tumor}': 'vaf_raw_tumor'}, inplace=True,
#        )
#        if vcf_sampleid_normal is not None:
#            vaf_df.rename(
#                columns={f'vaf_{vcf_sampleid_normal}': 'vaf_raw_normal'}, inplace=True,
#            )
#
#        # postprocess
#        #vaf_df = vaf_df.loc[vaf_df['vaf_raw'].notna().to_numpy(), :]
#        vaf_df.reset_index(drop=True, inplace=True)
#        vaf_df['baf_raw_tumor'] = cnvmisc.get_bafs(vaf_df['vaf_raw_tumor'])
#        if 'vaf_raw_normal' in vaf_df:
#            vaf_df['baf_raw_normal'] = cnvmisc.get_bafs(vaf_df['vaf_raw_normal'])
#
#        selected_cols = [
#            'Chromosome', 
#            'Start', 
#            'End', 
#            'vaf_raw_tumor', 
#            'baf_raw_tumor', 
#        ]
#        if vcf_sampleid_normal is not None:
#            selected_cols.append('vaf_raw_normal')
#            selected_cols.append('baf_raw_normal')
#        vaf_df = vaf_df.loc[:, selected_cols]
#        self.data[sampleid]['original_baf'] = vaf_df
#
#    #################
#    # other helpers #
#    #################
#
#    def set_normal_mean_ploidy(self, sampleid):
#        self.data[sampleid]['normal_mean_ploidy'] = cnvmisc.get_normal_mean_ploidy(
#            self.refver, 
#            self.data[sampleid]['is_female'], 
#            self.data[sampleid]['target_region'],
#        )
#
#    @staticmethod
#    def add_CNn_to_segment(
#        segment_df,
#        mode,
#        refver,
#        is_female,
#        target_region=None,
#    ):
#        if mode == 'wgs':
#            segment_df = rcopynumber.add_CNn_to_wgs_segment_gr(
#                segment_df, refver, is_female, as_gr=False,
#            )
#        elif mode == 'panel':
#            assert target_region is not None
#            segment_df = rcopynumber.add_CNn_to_targetseq_segment_gr(
#                segment_df, target_region, refver, is_female, as_gr=False,
#            )
#        return segment_df
#
#
## segmentation functions for running with multiprocessing
#
#def handle_gamma_kmin_args(gamma, kmin, mode):
#    if gamma is None:
#        if mode == 'wgs':
#            gamma = 40
#        elif mode == 'panel':
#            gamma = 30
#    if kmin is None:
#        if mode == 'wgs':
#            kmin = 5
#        elif mode == 'panel':
#            kmin = 1
#
#    return gamma, kmin
#
#
#def _make_depth_segment(
#    depthratio_df,
#    mode,
#    refver,
#    winsorize=False,
#    gamma=None,
#    kmin=None,
#    verbose=False,
#):
#    gamma, kmin = handle_gamma_kmin_args(gamma, kmin, mode)
#
#    #if 'depthratio_raw' in depthratio_df.columns:
#    #    input_df = depthratio_df.rename(columns={'depthratio_raw': 'depth_raw'})
#    #else:
#    #    input_df = depthratio_df
#    input_df = depthratio_df.rename(columns={'depthratio_raw': 'depth_raw'})
#    segment_df, _ = rcopynumber.run_rcopynumber_unified(
#        depth_df=input_df,
#        refver=refver,
#
#        as_gr=False, 
#        winsorize=winsorize,
#        #compact=(mode == 'panel'), 
#        compact=False,
#        verbose=verbose,
#        gamma=gamma,
#        kmin=kmin,
#    )
#
#    segment_df.rename(columns={'depth_segment_mean': 'depthratio_segment_mean'}, inplace=True)
#
#    segment_df = add_depthratio_std(segment_df, depthratio_df, refver)
#
#    return segment_df
#
#
#def add_depthratio_std(segment_df, raw_df, refver):
#    # make left df
#    left_df = segment_df.drop(
#        segment_df.columns.intersection(
#            ['depthratio_segment_std', 'depthratio_segment_std_mean_ratio'], 
#        ),
#        axis=1, 
#        inplace=False,
#    )
#    # make right df
#    right_df = raw_df.loc[
#        :, ['Chromosome', 'Start', 'End', 'depthratio_raw']
#    ].copy()
#    # join
#    joined_segdf = pyranges_helper.join(
#        left_df,
#        right_df,
#        how='left', merge='mean', add_std=True, ddof=0,
#        sort=True, refver=refver,
#    )
#    joined_segdf.drop('depthratio_raw', axis=1, inplace=True)
#    joined_segdf.rename(
#        columns={'depthratio_raw_std': 'depthratio_segment_std'}, 
#        inplace=True,
#    )
#
#    # add std/mean ratio
#    joined_segdf['depthratio_segment_std_mean_ratio'] = (
#        joined_segdf['depthratio_segment_std']
#        / joined_segdf['depthratio_segment_mean']
#    ).to_numpy()
#
#    return joined_segdf
#
#
#def _make_baf_segment(
#    baf_df,
#    target_region,
#    mode,
#    refver,
#    winsorize=False,
#    gamma=None,
#    kmin=None,
#    verbose=False,
#    baf_cutoff=0.1,
#    #bafcorrector=libbaf.load_bafcorrect_func(),
#):
#    gamma, kmin = handle_gamma_kmin_args(gamma, kmin, mode)
#
#    targetovlp_baf_df = pr.PyRanges(baf_df).overlap(target_region).df
#    input_df = targetovlp_baf_df.loc[:, ['Chromosome', 'Start', 'End']]
#    input_df['depth_raw'] = targetovlp_baf_df['baf_raw_tumor']
#    input_df = input_df.loc[input_df['depth_raw'] > baf_cutoff, :]
#
#    segment_df, _ = rcopynumber.run_rcopynumber_unified(
#        depth_df=input_df,
#        refver=refver,
#
#        as_gr=False, 
#        winsorize=winsorize,
#        #compact=(mode == 'panel'), 
#        compact=False,
#        verbose=verbose,
#        gamma=gamma,
#        kmin=kmin,
#    )
#
#    segment_df.rename(columns={'depth_segment_mean': 'baf_segment_mean'}, inplace=True)
#
#    return segment_df
#
#
#def _make_segments_targetfunc_depth(shdict, shdict_key, **kwargs):
#    shdict[shdict_key] = _make_depth_segment(**kwargs)
#
#
#def _make_segments_targetfunc_baf(shdict, shdict_key, **kwargs):
#    shdict[shdict_key] = _make_baf_segment(**kwargs)
#
#
#def _make_segments_main(
#    depthratio_df,
#    mode,
#    refver,
#    winsorize,
#    depthratio_gamma,
#    depthratio_kmin,
#    baf_gamma,
#    baf_kmin,
#    verbose,
#
#    baf_df,
#    target_region,
#    baf_cutoff,
#    #bafcorrector,
#):
#    with multiprocessing.Manager() as manager:
#        shdict = manager.dict()
#        subp1 = multiprocessing.Process(
#            target=_make_segments_targetfunc_depth,
#            args=(
#                shdict,
#                'depth_segment',
#            ),
#            kwargs=dict(
#                depthratio_df=depthratio_df,
#                mode=mode,
#                refver=refver,
#                winsorize=winsorize,
#                gamma=depthratio_gamma,
#                kmin=depthratio_kmin,
#                verbose=verbose,
#            ),
#        )
#        subp2 = multiprocessing.Process(
#            target=_make_segments_targetfunc_baf,
#            args=(
#                shdict, 
#                'baf_segment',
#            ),
#            kwargs=dict(
#                baf_df=baf_df,
#                target_region=target_region,
#                mode=mode,
#                refver=refver,
#                winsorize=winsorize,
#                gamma=baf_gamma,
#                kmin=baf_kmin,
#                verbose=verbose,
#                baf_cutoff=baf_cutoff,
#            ),
#        )
#
#        subp1.start()
#        subp2.start()
#
#        subp1.join()
#        subp2.join()
#
#        # result
#        depth_segment = shdict['depth_segment']
#        baf_segment = shdict['baf_segment']
#
#    return depth_segment, baf_segment
#
#
################################################################
#
#
#def make_targetseq_cnvplot(data_df, draw_invalid_regions=False):
#    # sanity check
#    required_cols = {
#        'CNt', 'B', 
#        'baf_segment_mean', 'baf_predicted', 'baf_raw',
#        'depthratio_segment_mean', 'depthratio_predicted', 'depthratio_raw',
#    }
#    assert required_cols.issubset(data_df.columns)
#
#    # setup
#    if isinstance(data_df, pr.PyRanges):
#        data_df = data_df.df
#
#    gplotter = GenomePlotter(data_df)
#    fig, axd = plt.subplot_mosaic(
#        [
#            ['CN',], 
#            ['baf',], 
#            ['depth',],
#        ],
#        figsize=(30, 13),
#    )
#    for ax in axd.values():
#        gplotter.set_xlim(ax)
#
#    # CN
#    axd['CN'].set_ylabel('CN')
#
#    gplotter.draw_hlines(
#        axd['CN'], df=data_df, y_colname='CNt', offset=0.1,
#        plot_kwargs={'color': 'black'},
#    )
#    gplotter.draw_hlines(
#        axd['CN'], df=data_df, y_colname='B', 
#        plot_kwargs={'color': 'blue'},
#    )
#    gplotter.fit_spines_to_regions(axd['CN'])
#
#    # baf
#    axd['baf'].set_ylabel('baf')
#    axd['baf'].set_ylim(0, 0.6)
#
#    gplotter.draw_dots(
#        axd['baf'], df=data_df, y_colname='baf_raw', 
#        plot_kwargs={'color': 'gray', 'markersize': 1.5, 'alpha': 0.5},
#    )
#    gplotter.draw_hlines(
#        axd['baf'], df=data_df, y_colname='baf_segment_mean', 
#        plot_kwargs={'color': 'black', 'linewidth': 2, 'alpha': 0.3},
#    )
#    gplotter.draw_hlines(
#        axd['baf'], df=data_df, y_colname='baf_predicted', 
#        plot_kwargs={'color': 'green', 'linewidth': 2, 'alpha': 0.3},
#    )
#    gplotter.fit_spines_to_regions(axd['baf'])
#
#    # depthratio
#    axd['depth'].set_ylabel('depth ratio')
#    axd['depth'].set_ylim(
#        0, 
#        data_df['depthratio_raw'].max() * 1.1
#    )
#
#    gplotter.draw_dots(
#        axd['depth'], df=data_df, y_colname='depthratio_raw', 
#        plot_kwargs={'color': 'gray', 'markersize': 1.5, 'alpha': 0.5},
#    )
#    gplotter.draw_hlines(
#        axd['depth'], df=data_df, y_colname='depthratio_segment_mean', 
#        plot_kwargs={'color': 'black', 'linewidth': 2, 'alpha': 0.5},
#    )
#    gplotter.draw_hlines(
#        axd['depth'], df=data_df, y_colname='depthratio_predicted', 
#        plot_kwargs={'color': 'green', 'linewidth': 2, 'alpha': 0.5},
#    )
#
#    if draw_invalid_regions:
#        selector = np.logical_or(
#            data_df['depthratio_raw'].isna().to_numpy,
#            data_df['baf_raw'].isna().to_numpy,
#        )
#        gplotter.draw_boxes(
#            axd['depth'], df=data_df.loc[selector, :], 
#            plot_kwargs=dict(color='red', alpha=0.01)
#        )
#
#    gplotter.fit_spines_to_regions(axd['depth'])
#
#    return fig, axd
#
#
#def draw_targetseq_cnvplot_from_data_precp(
#    tumor_depth_df, 
#    normal_depth_df,
#    germline_df,
#
#    region_gr,
#    refver, 
#    is_female, 
#):
#    assert isinstance(tumor_depth_df, pd.DataFrame)
#    assert isinstance(normal_depth_df, pd.DataFrame)
#    assert isinstance(germline_df, pd.DataFrame)
#
#    if 'baf_raw' not in germline_df.columns:
#        germline_df['baf_raw'] = cnvmisc.get_bafs(germline_df['vaf'])
#    germline_gr = pr.PyRanges(germline_df)
#
#    # make depth ratio
#    depthratio_gr = cnvmisc.make_depth_ratio(tumor_depth_df, normal_depth_df, as_gr=True)
#    depthratio_df = depthratio_gr.df
#
#    # run R copynumber
#    segment_gr, depth_baf_gr = rcopynumber.run_rcopynumber_unified(
#        depthratio_df.rename(columns={'depth_ratio_sequenzastyle': 'depth_raw'}),
#        baf_df=germline_df, 
#        refver=refver, 
#        as_gr=True, 
#        compact=True, 
#        winsorize=False,
#        gamma=30, 
#        kmin=1,
#    )
#    segment_gr = segment_gr[[]]
#
#    # postprocess segment df
#    segment_gr = rcopynumber.add_CNn_to_targetseq_segment_gr(segment_gr, region_gr, refver=refver, is_female=is_female)
#    segment_gr = pyranges_helper.join(
#        segment_gr, 
#        depthratio_gr[['depth_ratio_sequenzastyle']], 
#        how='left', merge='mean', find_nearest=False, as_gr=True,
#    )
#    segment_gr = pyranges_helper.join(
#        segment_gr, 
#        germline_gr[['baf_raw']], 
#        how='left', merge='mean', find_nearest=False, as_gr=True,
#    )
#    segment_df = segment_gr.df
#    segment_df.rename(
#        {'depth_ratio_sequenzastyle': 'depthratio_segment_mean', 'baf_raw': 'baf_segment_mean'},
#        inplace=True,
#        axis=1,
#    )
#
#    return segment_df, depth_baf_gr
#
#
#def draw_targetseq_cnvplot_from_data_choosecp(
#    segment_df, 
#    depth_baf_gr,
#
#    region_gr,
#    refver,
#    is_female,
#
#    CNt_weight=cnvmisc.DEFAULT_CNT_WEIGHT,
#    nproc=None,
#):
#    cpscore_dict = cnvmisc.get_cp_score_dict(segment_df, refver=refver, is_female=is_female, target_region_gr=region_gr, CNt_weight=CNt_weight, nproc=nproc)
#    peak_values, dfs = cnvmisc.get_peak_info(cpscore_dict)
#
#    cnvmisc.show_heatmap_peaks_new(dfs, figsize=(20, 20))
#
#    return peak_values, dfs, cpscore_dict
#
#
#def draw_targetseq_cnvplot_from_data_postcp(
#    cellularity,
#    ploidy,
#    segment_df, 
#    depth_baf_gr,
#
#    region_gr,
#    refver,
#    is_female,
#
#    draw_invalid_regions=False,
#):
#    normal_mean_ploidy = cnvmisc.get_normal_mean_ploidy(refver, is_female, target_region_gr=region_gr)
#
#    # postprocess segment df
#    segment_df = cnvmisc.add_CNt_to_segment(
#        segment_df, 
#        cellularity=cellularity, 
#        tumor_ploidy=ploidy, 
#        is_female=is_female, 
#        normal_ploidy=normal_mean_ploidy,
#    )
#    segment_df = cnvmisc.add_theoreticals_to_segment(
#        segment_df, 
#        cellularity=cellularity, 
#        tumor_ploidy=ploidy, 
#        normal_ploidy=normal_mean_ploidy, 
#        is_female=is_female,
#    )
#
#    # make a single merged plot data df
#    depth_baf_gr.depthratio_raw = depth_baf_gr.depth_raw
#    depth_baf_gr = cnvmisc.annotate_region_with_segment(depth_baf_gr, pr.PyRanges(segment_df), as_gr=True)
#
#    # draw plot
#    fig, axd = make_targetseq_cnvplot(depth_baf_gr, draw_invalid_regions=draw_invalid_regions)
#
#    return fig, axd
#
#

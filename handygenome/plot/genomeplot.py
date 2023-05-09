import functools
import itertools
import multiprocessing
import operator
import pprint

import numpy as np
import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

import handygenome.common as common
import handygenome.pyranges_helper as pyranges_helper
import handygenome.cnv.misc as cnvmisc
import handygenome.cnv.rcopynumber as rcopynumber
import handygenome.cnv.mosdepth as libmosdepth
import handygenome.variant.variantplus as variantplus
import handygenome.cnv.gcfraction as libgcfraction
import handygenome.deco as deco
import handygenome.workflow as workflow
import handygenome.ucscdata as ucscdata
import handygenome.cnv.baf as libbaf


LOGGER_INFO = workflow.get_debugging_logger(verbose=False)
LOGGER_DEBUG = workflow.get_debugging_logger(verbose=True)
DOT_ALPHA_CONST = 0.01 / 3e6


def check_is_allregion(region_gr, refver):
    """checks if allregion_gr is included within region_gr"""
    region_gr = cnvmisc.arg_into_gr(region_gr)
    allregion_gr = common.DEFAULT_CHROMDICTS[refver].to_gr(assembled_only=True, as_gr=True)
    isec_gr = allregion_gr.intersect(region_gr)
    return (isec_gr[[]].sort().df == allregion_gr[[]].sort().df).all(axis=None)


class CoordConverter:
    def __init__(self, refver, df, region_gaps):
        """Args:
            df: pandas.DataFrame object with mandatory columns 'Chromosome', 
                'Start', 'End'. Each row represents a genomic interval, which 
                will be drawn separately in the plot. The rows need not be 
                sorted by chromosome or position. The order of rows of the 
                input DataFrame is respected.
        """
        self.refver = refver
        self.chromdict = common.DEFAULT_CHROMDICTS[refver]

        cnvmisc.genome_df_sanitycheck(df)
        df = cnvmisc.arg_into_df(df)
        df = self.handle_df(df, region_gaps)
        self.set_params(df)

    @staticmethod
    def handle_df(df, region_gaps):
        df = df.copy()

        if 'weight' in df.columns:
            df = df.loc[
                :, ['Chromosome', 'Start', 'End', 'weight']
            ].copy().reset_index(drop=True)
        else:
            df = df.loc[
                :, ['Chromosome', 'Start', 'End']
            ].copy().reset_index(drop=True)
            df['weight'] = 1

        #if region_gaps:
        if (region_gaps is not None) and (df.shape[0] > 1):
            total_gap_length = ((df['End'] - df['Start']) * df['weight']).sum() * region_gaps
            gap_weight = total_gap_length / (df.shape[0] - 1)

            tmplist = list()
            for idx, row in df.iterrows():
                chrom = str(-idx - 1)
                gap_row = pd.Series(
                    {
                        'Chromosome': chrom, 
                        'Start': 0, 
                        'End': 1, 
                        'weight': gap_weight,
                    }
                )

                tmplist.append(row)
                tmplist.append(gap_row)
            tmplist = tmplist[:-1]

            df = pd.DataFrame.from_records(tmplist)

        return df

    @staticmethod
    def handle_df_old(df, region_gaps):
        df = df.copy()

        if 'weight' in df.columns:
            df = df.loc[
                :, ['Chromosome', 'Start', 'End', 'weight']
            ].copy().reset_index(drop=True)
        else:
            df = df.loc[
                :, ['Chromosome', 'Start', 'End']
            ].copy().reset_index(drop=True)
            df['weight'] = 1

        #if region_gaps:
        if region_gaps is not None:
            total_length = (df['End'] - df['Start']).sum()
            gap_length = int(total_length * region_gaps)
            median_weight = np.median(df['weight']).astype(int)

            tmplist = list()
            for idx, row in df.iterrows():
                chrom = str(-idx - 1)
                if 'weight' in df.columns:
                    gap_row = pd.Series({'Chromosome': chrom, 'Start': 0, 'End': gap_length, 'weight': median_weight})
                else:
                    gap_row = pd.Series({'Chromosome': chrom, 'Start': 0, 'End': gap_length})

                tmplist.append(row)
                tmplist.append(gap_row)
            tmplist = tmplist[:-1]

            df = pd.DataFrame.from_records(tmplist)

        return df

    def set_params(self, df):
        """Set attributes:
            totalregion_df
            totalregion_gr
            plot_interval_start0s
            plot_interval_end0s
        """
        # sanity check
        tmp_gr = pr.PyRanges(df)
        if tmp_gr.length != tmp_gr.merge().length:
            raise Exception(f'Plot region dataframe must not have overlapping intervals.')

        # set totalregion_df and totalregion_gr
        totalregion_df = df
        totalregion_df['raw_region_length'] = (
            totalregion_df['End'] - totalregion_df['Start']
        ).array
        totalregion_df['plot_region_length'] = (
            totalregion_df['raw_region_length'] * totalregion_df['weight']
        ).array

        cumsum = totalregion_df['plot_region_length'].cumsum()
        cumsum_shift = cumsum.shift(1, fill_value=0)
        #totalregion_df['region_start_offset'] = cumsum_shift.array
        totalregion_df['plot_interval_start0s'] = cumsum_shift.array
        totalregion_df['plot_interval_end0s'] = cumsum.array

        self.totalregion_df = totalregion_df
        self.totalregion_gr = pr.PyRanges(self.totalregion_df)

        # set chromosome-wise params
        self.chromwise_params = dict()
        for chrom, subdf in self.totalregion_gr.items():
            self.chromwise_params[chrom] = {
                'start0': np.array(subdf['Start']),
                'end0': np.array(subdf['End']),
                'raw_region_length': np.array(subdf['raw_region_length']),
                'plot_region_length': np.array(subdf['plot_region_length']),
                'plot_region_start_offset': np.array(subdf['plot_interval_start0s']),
            }

    def iter_totalregion_df(self):
        chroms = self.totalregion_df['Chromosome']
        grouper = (chroms != chroms.shift(1)).cumsum()
        return iter(self.totalregion_df.groupby(grouper))

    @property
    def xlim(self):
        start0 = self.totalregion_df['plot_interval_start0s'].iloc[0]
        end0 = self.totalregion_df['plot_interval_end0s'].iloc[-1] - 1
        return (start0, end0)

    def genomic_to_plot(self, chrom, pos0_list):
        if chrom not in self.chromwise_params.keys():
            raise Exception(f'Input "chrom" argument is not included in the plotting region.')

        pos0_list = np.array(pos0_list)[:, np.newaxis]
        params = self.chromwise_params[chrom]

        contains = np.logical_and(
            (pos0_list >= params['start0']), (pos0_list < params['end0'])
        )
        pos0s_indexes, intv_indexes = np.where(contains)
            # np.ndarray composed of the indexes of the containing intervals
            # identical intervals can appear many times
        within_region_offsets = (
            params['plot_region_length'][intv_indexes]
            * (
                (pos0_list[pos0s_indexes, 0] - params['start0'][intv_indexes]) 
                / params['raw_region_length'][intv_indexes]
            )
        )
        return params['plot_region_start_offset'][intv_indexes] + within_region_offsets

    def genomic_to_plot_with_indexes(self, chrom, pos0_list, indexes):
        plot_coords = self.genomic_to_plot(chrom, pos0_list)
        return (indexes, plot_coords)

    def plot_to_genomic(self, plotcoord_list):
        plotcoord_list = np.array(plotcoord_list)
        xlim = self.xlim
        assert np.logical_and(
            plotcoord_list >= xlim[0],
            plotcoord_list <= xlim[1],
        ).all(), f'Input plot coordinates are out of plot limits'

        plotcoord_list_expand = plotcoord_list[:, np.newaxis]
        compare_result = np.logical_and(
            plotcoord_list_expand >= self.totalregion_df['plot_interval_start0s'].to_numpy(),
            plotcoord_list_expand < self.totalregion_df['plot_interval_end0s'].to_numpy(),
        )
        input_indexes, plotregion_indexes = np.where(compare_result)
        assert (plotcoord_list[input_indexes] == plotcoord_list).all()

        # results
        totalregion_subdf = self.totalregion_df.iloc[plotregion_indexes, :]
        result_chroms = totalregion_subdf['Chromosome'].to_numpy()
        offsets = (
            (plotcoord_list - totalregion_subdf['plot_interval_start0s'])
            / totalregion_subdf['plot_region_length']
        ) * (totalregion_subdf['End'] - totalregion_subdf['Start'])
        result_pos0s = np.rint(totalregion_subdf['Start'] + offsets).astype(int)
        return result_chroms, result_pos0s

    def plot_to_genomic_old(self, x):
        plot_intvlist = self.data.index.get_level_values('plot_interval')
        genome_intvlist = self.data.index.get_level_values('genome_interval')
        contains = plot_intvlist.contains(x)

        num_hit = contains.sum()
        if num_hit == 0:
            return None
        elif num_hit > 1:
            raise Exception(f'More than one intervals contains the input position.')

        idx = np.where(contains)[0][0]

        chrom = self.data.index.get_level_values('chromosome')[idx]

        plot_intv = plot_intvlist[idx]
        genome_intv = genome_intvlist[idx]
        regional_offset_fraction = (x - plot_intv.left) / plot_intv.length
        pos0 = int(np.rint(genome_intv.left + (genome_intv.length * regional_offset_fraction)))

        return chrom, pos0

    # Axes modification
    def get_chrom_borders(self):
        result = list()
        for key, subdf in self.iter_totalregion_df():
            chroms = set(subdf['Chromosome'])
            assert len(chroms) == 1
            chrom = chroms.pop()
            result.append(
                (
                    chrom, 
                    subdf['plot_interval_start0s'].iloc[0], 
                    subdf['plot_interval_end0s'].iloc[-1],
                )
            )
        return result


class GenomePlotter:
    def __init__(
        self, refver, 
        *, 
        region_df=None, 
        chroms=None, start0s=None, end0s=None, weights=None,
        region_gaps=None,
    ):
        if region_gaps is None:
            region_gaps = 0.1
            
        self.refver = refver
        self.region_gaps = region_gaps

        if region_df is None:
            if chroms is None:
                region_df = common.DEFAULT_CHROMDICTS[self.refver].to_gr(
                    assembled_only=True, as_gr=False,
                )
            else:
                region_df = self.make_new_region_df(
                    self.refver, chroms, start0s=start0s, end0s=end0s, weights=weights
                )

        self.region_df = region_df
        self.cconv = CoordConverter(refver=refver, df=region_df, region_gaps=region_gaps)

    @property
    def totalregion_df(self):
        return self.cconv.totalregion_df

    def set_xlim(self, ax):
        ax.set_xlim(*self.cconv.xlim)

    def draw_genomecoord_labels(self, ax, n=10):
        """Should be done after data drawings are finished"""
        xlim = self.cconv.xlim
        plotcoords = np.linspace(xlim[0], xlim[1], num=n, endpoint=True)
        chroms, pos0s = self.cconv.plot_to_genomic(plotcoords)

        chroms = [common.prefix_chr(x) for x in chroms]
        pos1s = pos0s + 1
        pos1_strings = common.shorten_int(pos1s)

        labels = [f'{x} : {y}' for x, y in zip(chroms, pos1_strings)]
        ax.set_xticks(plotcoords, labels=labels, minor=False, rotation=90)

    def draw_chrom_borders(
        self, ax, 
        draw_chrom_names=True,
        prefix_with_chr=True,
        text_kwargs=dict(), 
        line_kwargs=dict(),
    ):
        """Should be done after data drawings are finished"""
        # set plotting kwargs
        default_text_kwargs = dict(
            ha='center',
            va='bottom',
            size=8,
        )
        default_text_kwargs.update(text_kwargs)

        default_line_kwargs = dict(
            color='black', 
            linewidth=1,
        )
        default_line_kwargs.update(line_kwargs)

        # main
        chrom_borders = self.cconv.get_chrom_borders()

        # draw chromosome name texts
        if draw_chrom_names:
            for chrom, start0, end0 in chrom_borders:
                if chrom.startswith('-'):
                    continue
                   
                if prefix_with_chr:
                    if not chrom.startswith('chr'):
                        chrom = 'chr' + chrom

                ax.text(
                    (start0 + end0) / 2, 
                    ax.get_ylim()[1], 
                    chrom, 
                    **default_text_kwargs,
                )

        # draw chromosome region borderlines
        border_pos0s = set()
        xlim = self.cconv.xlim
        for chrom, start0, end0 in chrom_borders:
            if start0 != xlim[0]:
                border_pos0s.add(start0)
            if end0 != xlim[1]:
                border_pos0s.add(end0)

        for pos0 in border_pos0s:
            ax.axvline(pos0, **default_line_kwargs)

    def draw_hlines(
        self, ax, y_colname, *, 
        df=None, df_plotdata=None, offset=0,
        plot_kwargs=dict(),
    ):
        default_plot_kwargs = {}
        default_plot_kwargs.update(plot_kwargs)

        if df_plotdata is None:
            df_plotdata = self.prepare_plot_data(df)

        if df_plotdata is False:
            return

        ys, xmins, xmaxs = self._merge_adjacent_hlines_data(
            df_plotdata[y_colname].to_numpy() + offset,
            df_plotdata['plot_start0s'],
            df_plotdata['plot_end0s'],
        )

        ax.hlines(ys, xmins, xmaxs, **default_plot_kwargs)

    def draw_dots(
        self, ax, y_colname, *, 
        df=None, df_plotdata=None,
        plot_kwargs=dict(),
    ):
        # kwargs handling
        plot_kwargs = (
            {
                'color': 'black',
                'marker': 'o',
                'linestyle': '',
            } | plot_kwargs
        )

        # main
        if df_plotdata is None:
            df_plotdata = self.prepare_plot_data(df)

        if df_plotdata is False:
            return

        xs = (df_plotdata['plot_start0s'] + (df_plotdata['plot_end0s'] - 1)) / 2
        ys = df_plotdata[y_colname].to_numpy()
        ax.plot(xs, ys, **plot_kwargs)

    def draw_dots_scatter(
        self, ax, y_colname, *, 
        df=None, df_plotdata=None,
        color_colname=None,
        plot_kwargs=dict(),
    ):
        # kwargs handling
        plot_kwargs = (
            {
                'c': 'black',
                'marker': 'o',
            } | plot_kwargs
        )

        # main
        if df_plotdata is None:
            df_plotdata = self.prepare_plot_data(df)

        if df_plotdata is False:
            return

        if df_plotdata is not None:
            xs = (df_plotdata['plot_start0s'] + (df_plotdata['plot_end0s'] - 1)) / 2
            ys = df_plotdata[y_colname].to_numpy()
            if color_colname is not None:
                plot_kwargs['c'] = df_plotdata[color_colname]

            ax.scatter(xs, ys, **plot_kwargs)

    def draw_bgcolors(
        self, 
        ax, 
        df=None, 
        df_plotdata=None,
        ymins=None,
        ymaxs=None,
        colors='yellow',
        plot_kwargs=dict(),
    ):
        # setup plot_kwargs
        default_plot_kwargs = {
            'alpha': 0.1,
            'zorder': 0,
        }
        default_plot_kwargs.update(plot_kwargs)

        # main
        if df_plotdata is None:
            df_plotdata = self.prepare_plot_data(df)

        if df_plotdata is False:
            return

        # x
        xmaxs = df_plotdata['plot_end0s']
        xmins = df_plotdata['plot_start0s']
        widths = xmaxs - xmins

        # y
        ylims = ax.get_ylim()
        if ymins is None:
            ymins = np.repeat(ylims[0], len(widths))
        if ymaxs is None:
            ymaxs = np.repeat(ylims[1], len(widths))
        heights = ymaxs - ymins

        # color
        #colors = common.arg_into_list(colors)
        if np.isscalar(colors):
        #if len(colors) == 1:
            match_original = False
            default_plot_kwargs['color'] = colors
        else:
            match_original = True

        # final
        boxes = [
            Rectangle((xm, ym), width=w, height=h)
            for (xm, w, ym, h) in zip(xmins, widths, ymins, heights)
        ]
        if match_original:
            for col, box in zip(colors, boxes):
                box.set(color=col)

        ax.add_collection(
            PatchCollection(
                boxes, 
                match_original=match_original, 
                **default_plot_kwargs,
            )
        )

    def draw_ideogram(self, ax):
        cytoband_df = ucscdata.get_cytoband(self.refver, as_gr=False)
        colors = [ucscdata.CYTOBAND_COLORMAP[x] for x in cytoband_df['Stain']]
        self.draw_bgcolors(
            ax=ax, 
            df=cytoband_df, 
            plot_kwargs=dict(alpha=1),
            colors=colors,
        )

    def prepare_plot_data(self, df):
        gr = cnvmisc.arg_into_gr(df)
        isec_gr = gr.intersect(self.cconv.totalregion_gr).sort()
        if isec_gr.empty:
            return False

        result_start0s = list()
        result_end0s = list()
        ordered_chroms = [x[0] for x in itertools.groupby(isec_gr.Chromosome)]
        for chrom in ordered_chroms:
            subgr = isec_gr[chrom]
            result_start0s.extend(
                self.cconv.genomic_to_plot(chrom, subgr.Start)
            )
            result_end0s.extend(
                self.cconv.genomic_to_plot(chrom, subgr.End - 1) + 1
            )

        isec_gr.plot_start0s = result_start0s
        isec_gr.plot_end0s = result_end0s

        return isec_gr.df

    def prepare_plot_data_old(self, df, nproc=None):
        # create isec between total region and input data
        isec_gr, subgrs_bychrom = self._isec_trim_data_df(df)

        # Add "End_minus1" columns; "End" columns cannot be used for plot coordinate calculation
        for chrom, subgr in subgrs_bychrom.items():
            subgr.End_minus1 = subgr.End - 1

        xmins = self._get_ordered_plot_coords(subgrs_bychrom, 'Start', nproc=nproc)
        xmaxs_minus1 = self._get_ordered_plot_coords(subgrs_bychrom, 'End_minus1', nproc=nproc)
        xmaxs = xmaxs_minus1 + 1

        return {
            'isec_gr': isec_gr,
            'subgrs_bychrom': subgrs_bychrom,
            'xmins': xmins,
            'xmaxs': xmaxs,
        }

    ############################################

    @staticmethod
    def make_new_region_df(refver, chroms, start0s=None, end0s=None, weights=None):
        assert (start0s is None) == (end0s is None)
        # weights
        if weights is None:
            weights = 1

        weights = np.atleast_1d(weights)
        #if not all(str(x).isdecimal() for x in weights):
            #raise Exception(f'"weights" argument must be all integers')
        if not (weights > 0).all():
            raise Exception(f'All weights values must be greater than 0')

        # chroms
        chroms = np.atleast_1d(chroms)

        # start0s, end0s
        if start0s is None:
            chromdict = common.DEFAULT_CHROMDICTS[refver]
            start0s = np.repeat(0, len(chroms))
            end0s = np.fromiter((chromdict[x] for x in chroms), dtype=int)

        # broadcast
        chroms, start0s, end0s, weights = np.broadcast_arrays(
            chroms, start0s, end0s, weights,
        )

        return pd.DataFrame({
            'Chromosome': chroms,
            'Start': start0s,
            'End': end0s,
            'weight': weights,
        })

    @classmethod
    def _merge_adjacent_hlines_data(cls, ys, xmins, xmaxs):
        """Helper of draw_hlines"""
        ys = np.array(ys)
        xmins = np.array(xmins)
        xmaxs = np.array(xmaxs)

        flags = (ys[:-1] == ys[1:]) & (xmaxs[:-1] == xmins[1:])
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
        new_xmins = xmins[[x[0] for x in indexes]]
        new_xmaxs = xmaxs[[x[1] for x in indexes]]

        return new_ys, new_xmins, new_xmaxs

    def _isec_trim_data_df(self, df, asis=False):
        """helper of prepare_plot_data"""
        assert '_index' not in df.columns

        gr = cnvmisc.arg_into_gr(df)
        if asis:
            isec_gr = gr
        else:
            isec_gr = gr.intersect(self.cconv.totalregion_gr)

        isec_gr._index = list(range(isec_gr.df.shape[0]))

        subgrs_bychrom = dict()
        for chrom in isec_gr.Chromosome.unique():
            subgrs_bychrom[chrom] = isec_gr[chrom]

        return isec_gr, subgrs_bychrom

    def _get_ordered_plot_coords(self, subgrs_bychrom, pos0_colname, nproc=None):
        """helper of prepare_plot_data"""
        # Multiprocessing is slower than serial jobs!!
#        with multiprocessing.Pool(nproc) as pool:
#            result = pool.starmap(
#                self.genomic_to_plot_with_indexes, 
#                (
#                    (chrom, subdf[pos0_colname], subdf['_index']) 
#                    for chrom, subdf in subgrs_bychrom.items()
#                )
#            )

        result = list()
        for chrom, subgr in subgrs_bychrom.items():
            result.append(
                self.cconv.genomic_to_plot_with_indexes(
                    chrom, getattr(subgr, pos0_colname), subgr._index
                )
            )

        index_coord_pairs = itertools.chain.from_iterable(zip(*x) for x in result)
        plot_coords = np.fromiter(
            (x[1] for x in sorted(index_coord_pairs, key=operator.itemgetter(0))),
            dtype=np.int_,
        )

        return plot_coords

#    def _get_ordered_plot_coords_woidx(self, subgrs_bychrom, pos0_colname):
#        """helper of prepare_plot_data"""
#        return np.fromiter(
#            itertools.chain.from_iterable(
#                self.cconv.genomic_to_plot(chrom, getattr(subgr, pos0_colname))
#                for chrom, subgr in subgrs_bychrom.items()
#            ),
#            dtype=int,
#        )


class CNVPlotter:
    def __init__(self, refver, region_df=None):
        #super().__init__(refver=refver, region_df=region_df)
        self.refver = refver
        self.genomeplotter = GenomePlotter(refver, region_df=region_df)
        self.data = dict()
        self.default_binsize = 100

    def reset_genomeplotter(self, region_df=None):
        self.genomeplotter = GenomePlotter(self.refver, region_df=region_df)

    ##############
    # mainstream #
    ##############

    @deco.get_deco_num_set_differently(('normal_bam_path', 'normal_depth_path'), 1)
    @deco.get_deco_num_set_differently(('tumor_bam_path', 'tumor_depth_path'), 1)
    def add_sample_file(
        self, 
        sampleid, 
        is_female,
        tumor_vcf_path,
        tumor_vcf_sampleid,
        *,
        mode='wgs',
        target_region=None,

        normal_bam_path=None,
        normal_depth_path=None, 

        tumor_bam_path=None,
        tumor_depth_path=None, 

        vcfload_nproc=1,
        mosdepth_postprocess_kwargs=dict(),
    ):
        """Args:
            *_depth_path: mosdepth output file
            tumor_vcf_path: germline variant vcf
        """
        def helper(
            self, bam_path, depth_path, sampleid, sampletype, mosdepth_postprocess_kwargs,
        ):
            if bam_path is not None:
                self.load_bam_for_depth(
                    sampleid, 
                    bam_path=bam_path, 
                    sampletype=sampletype, 
                    use_default_gc=True, 
                    **mosdepth_postprocess_kwargs,
                )
            elif depth_path is not None:
                mosdepth_df = libmosdepth.load_mosdepth_output(depth_path, as_gr=False)
                self.load_mosdepth_df(
                    sampleid, 
                    mosdepth_df, 
                    sampletype=sampletype, 
                    use_default_gc=True, 
                    **mosdepth_postprocess_kwargs,
                )
        
        # sanity check
        assert mode in ('wgs', 'panel')
        if (mode == 'panel') and (target_region is None):
            raise Exception(f'"target_region" must be given when "mode" is "panel"')
        elif (mode == 'wgs') and (target_region is not None):
            raise Exception(f'"target_region" must not be given when "mode" is "wgs"')

        # main
        self.data[sampleid] = dict()
        self.data[sampleid]['is_female'] = is_female
        self.data[sampleid]['mode'] = mode

        # load normal
        LOGGER_INFO.info('Loading normal depth')
        helper(
            self, 
            normal_bam_path, 
            normal_depth_path, 
            sampleid, 
            'normal', 
            mosdepth_postprocess_kwargs=(
                dict(exclude_y=self.data[sampleid]['is_female'])
                | mosdepth_postprocess_kwargs
            ),
        )

        # set target_region
        chromdict = common.DEFAULT_CHROMDICTS[self.refver]
        all_assembled_chroms = chromdict.to_gr(assembled_only=True, as_gr=True)
        all_chroms = chromdict.to_gr(assembled_only=False, as_gr=True)
        if mode == 'wgs':
            # as non-nan region from normal_depth
            self.data[sampleid]['excluded_region'] = pr.PyRanges(
                self.data[sampleid]['normal_depth'].loc[
                    self.data[sampleid]['normal_depth']['excluded'], :
                ]
            ).drop()
            self.data[sampleid]['target_region'] = all_assembled_chroms.subtract(
                self.data[sampleid]['excluded_region']
            )
        elif mode == 'panel':
            self.data[sampleid]['target_region'] = cnvmisc.arg_into_gr(target_region).drop()
            if self.data[sampleid]['is_female']:
                self.data[sampleid]['target_region'] = (
                    self.data[sampleid]['target_region'].subtract(
                        all_assembled_chroms[chromdict.XY_names[1]]
                    )
                )

            self.data[sampleid]['excluded_region'] = all_chroms.subtract(
                self.data[sampleid]['target_region']
            )

        # normal mean ploidy
        self.set_normal_mean_ploidy(sampleid)

        # load tumor, using excluded region from normal depth
        LOGGER_INFO.info('Loading tumor depth')
        helper(
            self, 
            tumor_bam_path, 
            tumor_depth_path, 
            sampleid, 
            'tumor', 
            mosdepth_postprocess_kwargs=(
                dict(outlier_region=self.data[sampleid]['excluded_region'])
                | mosdepth_postprocess_kwargs
            ),
        )

        # sort tumor and normal depth dfs
        self.data[sampleid]['tumor_depth'] = cnvmisc.sort_genome_df(
            self.data[sampleid]['tumor_depth'], self.refver,
        )
        self.data[sampleid]['normal_depth'] = cnvmisc.sort_genome_df(
            self.data[sampleid]['normal_depth'], self.refver,
        )

        # load germline vcf
        LOGGER_INFO.info('Loading tumor germline vcf')
        self.load_germline_vcf(
            sampleid, tumor_vcf_path, tumor_vcf_sampleid, logging_lineno=50000,
            nproc=vcfload_nproc,
        )

    @deco.get_deco_num_set_differently(('normal_bam_path', 'normal_depth_path'), 1)
    @deco.get_deco_num_set_differently(('tumor_bam_path', 'tumor_depth_path'), 1)
    def add_sample_file_new(
        self, 
        sampleid, 
        is_female,
        tumor_vcf_path,
        vcf_sampleid_tumor,
        vcf_sampleid_normal=None,

        *,

        mode='wgs',
        target_region=None,

        normal_bam_path=None,
        normal_depth_path=None, 

        tumor_bam_path=None,
        tumor_depth_path=None, 

        vcfload_nproc=1,
        #mosdepth_postprocess_kwargs=dict(),
    ):
        """Args:
            *_depth_path: mosdepth output file
            tumor_vcf_path: germline variant vcf
            vcf_sampleids: tuple of (normal sample id, tumor sample id)
        """
        def helper(
            self, bam_path, depth_path, sampleid, sampletype,
        ):
            assert sampletype in ('tumor', 'normal')

            if bam_path is not None:
                mosdepth_df = libmosdepth.run_mosdepth(
                    bam_path, 
                    t=8, 
                    use_median=False, 
                    region_bed_path=None, 
                    region_gr=self.data[sampleid]['target_region'], 
                    window_size=self.default_binsize, 
                    donot_subset_bam=True,
                    as_gr=False, 
                    load_perbase=False,
                )
            elif depth_path is not None:
                mosdepth_df = libmosdepth.load_mosdepth_output(depth_path, as_gr=False)

            self.data.setdefault(sampleid, dict())
            self.data[sampleid][f'{sampletype}_depth'] = mosdepth_df
        
        # sanity check
        assert mode in ('wgs', 'panel')
        if (mode == 'panel') and (target_region is None):
            raise Exception(f'"target_region" must be given when "mode" is "panel"')
        elif (mode == 'wgs') and (target_region is not None):
            raise Exception(f'"target_region" must not be given when "mode" is "wgs"')

        if (mode == 'wgs') and (vcf_sampleid_normal is None):
            raise Exception(f'"vcf_sampleid_normal" must be given when "mode" is "wgs"')

        # main
        self.data[sampleid] = dict()
        self.data[sampleid]['is_female'] = is_female
        self.data[sampleid]['mode'] = mode

        # load normal
        LOGGER_INFO.info('Loading normal depth')
        helper(
            self, 
            normal_bam_path, 
            normal_depth_path, 
            sampleid, 
            'normal', 
        )

        # load tumor
        LOGGER_INFO.info('Loading tumor depth')
        helper(
            self, 
            tumor_bam_path, 
            tumor_depth_path, 
            sampleid, 
            'tumor', 
        )

        # set target_region
        self.set_target_region(sampleid, mode, target_region)

        # normal mean ploidy
        self.set_normal_mean_ploidy(sampleid)

        # load germline vcf
        LOGGER_INFO.info('Loading tumor germline vcf')
        self.load_germline_vcf(
            sampleid=sampleid, 
            vcf_path=tumor_vcf_path, 
            vcf_sampleid_tumor=vcf_sampleid_tumor,
            vcf_sampleid_normal=vcf_sampleid_normal,
            logging_lineno=50000,
            nproc=vcfload_nproc,
        )

    def postprocess_depth(self, sampleid):
        gc_df = libgcfraction.get_gc_df(
            self.refver, 
            self.default_binsize, 
            coords_as_index=True,
        )
        for key in ('normal', 'tumor'):
            output_depth_df, gcbin_average_depths = cnvmisc.postprocess_depth_df(
                self.data[sampleid][f'{key}_depth'],
                gc_df=gc_df,
                included_region=self.data[sampleid]['target_region'],
                as_gr=False,
            )
            self.data[sampleid][f'{key}_depth'] = output_depth_df
            self.data[sampleid][f'{key}_gcdata'] = gcbin_average_depths

    def set_depthratio(self, sampleid):
        depthratio_df = cnvmisc.make_depth_ratio(
            self.data[sampleid]['tumor_depth'], 
            self.data[sampleid]['normal_depth'],
            make_depthratio_mystyle=False,
            as_gr=False,
        )
        depthratio_df.rename(
            columns={'depth_ratio_sequenzastyle': 'depthratio_raw'}, inplace=True,
        )
        #assert (depthratio_df['normal_excluded'] == depthratio_df['tumor_excluded']).all()
        #depthratio_df['excluded'] = depthratio_df['normal_excluded']
        #depthratio_df.drop(['normal_excluded', 'tumor_excluded'], axis=1, inplace=True)

        self.data[sampleid]['depthratio'] = depthratio_df

    def upscale_preprocessing(self, input_df):
        result = input_df.copy()
        annot_cols = cnvmisc.get_genome_df_annotcols(input_df)
        result.loc[result['excluded'], annot_cols] = np.nan
        result.drop('excluded', axis=1, inplace=True)
        return result

    def upscale_depthratio(self, sampleid, binsize=1000):
        input_df = self.upscale_preprocessing(self.data[sampleid]['depthratio'])
        self.data[sampleid]['depthratio_upscaled'] = cnvmisc.upsize_depth_df_bin(
            input_df, 
            size=binsize, 
            refver=self.refver,
        )

    def upscale_depth(self, sampleid, binsize=10000):
        self.data[sampleid]['normal_depth_upscaled'] = cnvmisc.upsize_depth_df_bin(
            self.data[sampleid]['normal_depth'], 
            size=binsize, 
            refver=self.refver,
        )
        self.data[sampleid]['tumor_depth_upscaled'] = cnvmisc.upsize_depth_df_bin(
            self.data[sampleid]['tumor_depth'], 
            size=binsize, 
            refver=self.refver,
        )

#    def set_depth_baf_merge(
#        self, 
#        sampleid, 
#        nproc=1,
#    ):
#        depth_df = self.data[sampleid]['depthratio_upscaled'].loc[
#            :, ['Chromosome', 'Start', 'End', 'depth_raw']
#        ]
#        depth_df = cnvmisc.remove_unassembled_contigs(depth_df)
#
#        baf_df = self.data[sampleid]['tumor_baf'].loc[
#            :, ['Chromosome', 'Start', 'End', 'baf_raw']
#        ]
#        baf_df = cnvmisc.remove_unassembled_contigs(baf_df)
#
#        merged_df, merged_df_wona = rcopynumber.make_merged_df(
#            depth_df=depth_df, 
#            baf_df=baf_df, 
#            cytoband_df=None, 
#            refver=self.refver, 
#            logger=LOGGER_INFO, 
#            nproc=nproc,
#        )
#        self.data[sampleid]['depth_baf_merge'] = merged_df
#        self.data[sampleid]['depth_baf_merge_wona'] = merged_df_wona

    @staticmethod
    def handle_gamma_kmin_args(gamma, kmin, mode):
        if gamma is None:
            if mode == 'wgs':
                gamma = 40
            elif mode == 'panel':
                gamma = 30
        if kmin is None:
            if mode == 'wgs':
                kmin = 5
            elif mode == 'panel':
                kmin = 1

        return gamma, kmin

    def fit_segment_to_targetregion(self, sampleid, segment_df):
        segment_gr = cnvmisc.arg_into_gr(segment_df)
        seg_subset = segment_gr.intersect(self.data[sampleid]['target_region'])
        target_diff_seg = self.data[sampleid]['target_region'].subtract(segment_gr)
        target_diff_seg = pyranges_helper.join(
            target_diff_seg, 
            segment_gr,
            how='left',
            merge=None,
            find_nearest=True,
            as_gr=True,
        )
        result = pr.concat([seg_subset, target_diff_seg]).df
        result = cnvmisc.sort_genome_df(result, self.refver)

        return result

    def make_depth_segment(
        self,
        sampleid,
        winsorize=False,
        gamma=None,
        kmin=None,
        verbose=False,
    ):
        gamma, kmin = self.handle_gamma_kmin_args(gamma, kmin, self.data[sampleid]['mode'])

        if 'depthratio_raw' in self.data[sampleid]['depthratio_upscaled'].columns:
            input_df = self.data[sampleid]['depthratio_upscaled'].rename(
                columns={'depthratio_raw': 'depth_raw'}
            )
        else:
            input_df = self.data[sampleid]['depthratio_upscaled']

        segment_df, _ = rcopynumber.run_rcopynumber_unified(
            depth_df=input_df,
            refver=self.refver,

            as_gr=False, 
            winsorize=winsorize,
            compact=(True if self.data[sampleid]['mode'] == 'panel' else False), 
            verbose=verbose,
            gamma=gamma,
            kmin=kmin,
        )
        #segment_df.drop(columns='depth_segment_mean', inplace=True)
        segment_df.rename(columns={'depth_segment_mean': 'depthratio_segment_mean'}, inplace=True)
        self.data[sampleid]['depthratio_segment'] = segment_df

    def make_baf_segment(
        self,
        sampleid,
        winsorize=False,
        gamma=None,
        kmin=None,
        verbose=False,
        bafcorrector=libbaf.load_bafcorrect_func(),
    ):
        gamma, kmin = self.handle_gamma_kmin_args(gamma, kmin, self.data[sampleid]['mode'])

        filtered_baf_df = pr.PyRanges(self.data[sampleid]['baf']).overlap(
            self.data[sampleid]['target_region']
        ).df
        input_df = filtered_baf_df.loc[:, ['Chromosome', 'Start', 'End']]
        input_df['depth_raw'] = filtered_baf_df['baf_raw_tumor']

        segment_df, _ = rcopynumber.run_rcopynumber_unified(
            depth_df=input_df,
            refver=self.refver,

            as_gr=False, 
            winsorize=winsorize,
            compact=(self.data[sampleid]['mode'] == 'panel'), 
            verbose=verbose,
            gamma=gamma,
            kmin=kmin,
        )
        #segment_df.drop(columns='depth_segment_mean', inplace=True)
        segment_df.rename(columns={'depth_segment_mean': 'baf_segment_mean'}, inplace=True)
        segment_df['corrected_baf_segment_mean'] = bafcorrector(segment_df['baf_segment_mean'])

        self.data[sampleid]['baf_segment'] = segment_df

    def make_merged_segment(self, sampleid):
        merged_segment = pyranges_helper.isec_union(
            self.data[sampleid]['depthratio_segment'],
            self.data[sampleid]['baf_segment'],
        )
        merged_segment = pyranges_helper.join(
            merged_segment, 
            self.data[sampleid]['depthratio_segment'],
            how='left',
            merge=None,
            find_nearest=True,
            sort=True,
            refver=self.refver,
        )
        merged_segment = pyranges_helper.join(
            merged_segment, 
            self.data[sampleid]['baf_segment'],
            how='left',
            merge=None,
            find_nearest=True,
            sort=True,
            refver=self.refver,
        )

        # add CNn
        merged_segment = self._add_CNn_to_segment(
            merged_segment,
            self.data[sampleid]['mode'],
            self.refver,
            self.data[sampleid]['is_female'],
            self.data[sampleid]['target_region'],
        )

        # fit to target region
        #merged_segment = self.fit_segment_to_targetregion(sampleid, merged_segment)

        self.data[sampleid]['merged_segment'] = merged_segment

#    def make_segment(
#        self,
#        sampleid,
#        winsorize=False,
#        gamma=None,
#        kmin=None,
#        verbose=False,
#    ):
#        # set params
#        gamma, kmin = self.handle_gamma_kmin_args(gamma, kmin, self.data[sampleid]['mode'])
#
#        # run R-copynumber
#        if 'depthratio_raw' in self.data[sampleid]['depth_baf_merge'].columns:
#            self.data[sampleid]['depth_baf_merge'].rename(
#                columns={'depthratio_raw': 'depth_raw'}, inplace=True,
#            )
#
#        segment_df, depth_baf_df = rcopynumber.run_rcopynumber_unified(
#            merged_df=self.data[sampleid]['depth_baf_merge'], 
#            merged_df_wona=self.data[sampleid]['depth_baf_merge_wona'],
#
#            as_gr=False, 
#            winsorize=winsorize,
#            compact=(True if self.data[sampleid]['mode'] == 'panel' else False), 
#            verbose=verbose,
#            gamma=gamma,
#            kmin=kmin,
#        )
#
#        # modify segment_df
#        segment_df.rename(
#            columns={'depth_segment_mean': 'depthratio_segment_mean'}, inplace=True,
#        )
#
#        # add CNn
#        segment_df = self._add_CNn_to_segment(
#            segment_df,
#            self.data[sampleid]['mode'],
#            self.refver,
#            self.data[sampleid]['is_female'],
#            self.data[sampleid]['target_region'],
#        )
#
#        # recalc mean depthratio and baf with panel data
#        if self.data[sampleid]['mode'] == 'panel':
#            segment_df = segment_df.loc[:, ['Chromosome', 'Start', 'End', 'CNn']]
#
#            segment_df = pyranges_helper.join(
#                segment_df, 
#                self.data[sampleid]['depthratio'].loc[
#                    :, ['Chromosome', 'Start', 'End', 'depth_raw']
#                ], 
#                how='left', merge='mean', find_nearest=False, as_gr=False,
#            )
#            segment_df.rename(
#                columns={'depth_raw': 'depthratio_segment_mean'}, inplace=True,
#            )
#
#            segment_df = pyranges_helper.join(
#                segment_df, 
#                self.data[sampleid]['tumor_baf'].loc[
#                    :, ['Chromosome', 'Start', 'End', 'baf_raw']
#                ],
#                how='left', merge='mean', find_nearest=False, as_gr=True,
#            )
#            segment_df.rename(
#                columns={'baf_raw': 'baf_segment_mean'}, inplace=True,
#            )
#
#        # modify depth_baf_df
#        depth_baf_df.rename(
#            columns={'depth_raw': 'depthratio_raw'}, inplace=True,
#        )
#
#        # assign
#        self.data[sampleid]['segment'] = segment_df
#        self.data[sampleid]['depth_baf_merge'] = depth_baf_df

    def get_valid_CNts(
        self, sampleid, depthratio1, depthratio2, CNt1_range=range(0, 6), CNt2_maxdiff=6,
    ):
        return cnvmisc.get_valid_CNts_from_depthratios(
            depthratio1, 
            depthratio2, 
            normal_mean_ploidy=self.data[sampleid]['normal_mean_ploidy'], 
            CNn=2, 
            tumor_avg_depth_ratio=1, 
            normal_avg_depth_ratio=1, 
            CNt1_range=CNt1_range, 
            CNt2_maxdiff=CNt2_maxdiff,
        )

    def get_cp_from_twodata(
        self, sampleid, depthratio1, CNt1, depthratio2, CNt2,
    ):
        return cnvmisc.get_cp_from_twodata(
            depthratio1, 
            CNt1, 
            depthratio2, 
            CNt2, 
            normal_mean_ploidy=self.data[sampleid]['normal_mean_ploidy'], 
            CNn=2, 
            tumor_avg_depth_ratio=1, 
            normal_avg_depth_ratio=1,
        )

    def calc_cpscore(
        self, 
        sampleid, 
        CNt_weight=cnvmisc.DEFAULT_CNT_WEIGHT,
        nproc=1,
    ):
        cpscore_dict = cnvmisc.get_cp_score_dict(
            self.data[sampleid]['merged_segment'], 
            refver=self.refver, 
            is_female=self.data[sampleid]['is_female'], 
            target_region_gr=self.data[sampleid]['target_region'], 
            CNt_weight=CNt_weight, 
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
            nproc=nproc,
        )
        peak_values, dfs = cnvmisc.get_peak_info(cpscore_dict)

        self.data[sampleid]['cpscores'] = cpscore_dict
        self.data[sampleid]['peak_values'] = peak_values
        self.data[sampleid]['peak_dfs'] = dfs

    def show_peaks(self, sampleid, figsize=(20, 20), **kwargs):
        fig, ax = cnvmisc.show_heatmap_peaks_new(
            self.data[sampleid]['peak_dfs'], figsize=figsize, **kwargs,
        )

    def postprocess_segment(self, sampleid, cellularity, ploidy):
        # add CNt and B
        cpinfo = self.data[sampleid]['cpscores'][(cellularity, ploidy)]
        self.data[sampleid]['merged_segment']['CNt'] = cpinfo['CNt_list']
        self.data[sampleid]['merged_segment']['B'] = cpinfo['B_list']
        self.data[sampleid]['merged_segment']['A'] = (
            self.data[sampleid]['merged_segment']['CNt']
            - self.data[sampleid]['merged_segment']['B']
        )

        # add theoreticals
        self.data[sampleid]['merged_segment'] = cnvmisc.add_theoreticals_to_segment(
            self.data[sampleid]['merged_segment'], 
            cellularity=cellularity, 
            tumor_ploidy=ploidy, 
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'], 
            is_female=self.data[sampleid]['is_female'],
        )

    def add_CN_pred_to_segment(self, sampleid, cellularity, ploidy):
        if 'cpscores' not in self.data[sampleid]:
            self.data[sampleid]['cpscores'] = dict()

        # add CNt and B
        cpinfo = cnvmisc.calc_cp_score(
            segment_df=self.data[sampleid]['merged_segment'], 
            cellularity=cellularity, 
            ploidy=ploidy, 
            is_female=self.data[sampleid]['is_female'], 
            CNt_weight=1, 
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
        )
        self.data[sampleid]['cpscores'][(cellularity, ploidy)] = cpinfo

        self.data[sampleid]['merged_segment']['CNt'] = cpinfo['CNt_list']
        self.data[sampleid]['merged_segment']['B'] = cpinfo['B_list']
        self.data[sampleid]['merged_segment']['A'] = (
            self.data[sampleid]['merged_segment']['CNt']
            - self.data[sampleid]['merged_segment']['B']
        )

        # add theoreticals
        self.data[sampleid]['merged_segment'] = cnvmisc.add_theoreticals_to_segment(
            self.data[sampleid]['merged_segment'], 
            cellularity=cellularity, 
            tumor_ploidy=ploidy, 
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'], 
            is_female=self.data[sampleid]['is_female'],
        )

    def calculate_tumor_ploidy(self, sampleid):
        segment_df = self.data[sampleid]['merged_segment']
        if self.data[sampleid]['mode'] == 'wgs':
            weights = segment_df['End'] - segment_df['Start']
        else:
            stripped_segment_df = segment_df.loc[:, ['Chromosome', 'Start', 'End']].copy()
            all_indexes = list(range(stripped_segment_df.shape[0]))  # index of segments
            stripped_segment_df['index'] = all_indexes

            target_region = cnvmisc.arg_into_gr(self.data[sampleid]['target_region'])
            index_annotated_targetregion_df = cnvmisc.annotate_region_with_segment(  
                # each region is annotated with corresponding segment index
                target_region[[]],
                stripped_segment_df,
                as_gr=False,
            )

            index_annotated_targetregion_df['length'] = (
                index_annotated_targetregion_df['End']
                - index_annotated_targetregion_df['Start']
            )
            weights_dict = index_annotated_targetregion_df.loc[
                :, ['length', 'index']
            ].groupby('index').sum().to_dict()['length']
            weights = [
                (weights_dict[x] if x in weights_dict else 0)
                for x in all_indexes   
            ]

        return np.average(segment_df['CNt'], weights=weights)

    ############
    # plotting #
    ############

    def draw_spines_grids(self, ax, yticks):
        # remove existing grids and top/bottom spines
        #ax.grid(visible=False)
        ax.spines[['top', 'bottom']].set_visible(False)

        ylims = ax.get_ylim()
        chroms, start0s, end0s = zip(
            *[
                x for x in self.genomeplotter.cconv.get_chrom_borders() 
                if not x[0].startswith('-')
            ]
        )

        # spines
        ax.hlines(
            np.repeat(ylims[1], len(start0s)), start0s, end0s, color='black', linewidth=1,
        )
        ax.hlines(
            np.repeat(ylims[0], len(start0s)), start0s, end0s, color='black', linewidth=1.5,
        )
        # grids
        for y in yticks:
            ax.hlines(
                np.repeat(y, len(start0s)), start0s, end0s, 
                color='black', linewidth=0.5, alpha=0.5,
            )

    def draw_centromeres(self, ax):
        cytoband_gr = ucscdata.get_cytoband_gr(refver=self.refver, as_gr=True)

        self.genomeplotter.draw_bgcolors(
            ax, 
            df=cytoband_gr[cytoband_gr.Stain == 'acen'], 
            #df=cytoband_gr[cytoband_gr.Stain.isin(['acen', 'gvar', 'stalk'])], 
            colors='red',
            plot_kwargs=dict(alpha=0.3, linewidth=0),
        )

    def draw_centromeres_type2(self, ax):
        cytoband_gr = ucscdata.get_cytoband_gr(refver=self.refver, as_gr=True)

        mapping = {
            'acen': 'red', 
            'gvar': 'green', 
            'stalk': 'blue',
        }

        def helper(bandname):
            self.genomeplotter.draw_bgcolors(
                ax, 
                df=cytoband_gr[cytoband_gr.Stain == bandname], 
                colors=mapping[bandname],
                plot_kwargs=dict(alpha=0.3, linewidth=0),
            )

        helper('acen')
        helper('gvar')
        helper('stalk')

    def make_plotdata(self, sampleid):
        self.data[sampleid]['segment_plotdata'] = self.genomeplotter.prepare_plot_data(
            self.data[sampleid]['segment']
        )
        self.data[sampleid]['raw_plotdata'] = self.genomeplotter.prepare_plot_data(
            self.data[sampleid]['depth_baf_merge']
        )
        self.data[sampleid]['baf_plotdata'] = self.genomeplotter.prepare_plot_data(
            self.data[sampleid]['baf']
        )

    def make_plotdata_before_cp(self, sampleid, use_saved_plotdata, use_merged_segment):
        raw_plot_map = {
            'depthratio_upscaled': 'depthratio_raw_plotdata',
            'baf': 'baf_raw_plotdata',
            'merged_segment': 'merged_segment_plotdata',
            'depthratio_segment': 'depthratio_segment_plotdata',
            'baf_segment': 'baf_segment_plotdata',
        }
        def helper(raw_key):
            plot_key = raw_plot_map[raw_key]
            if (
                (not use_saved_plotdata) 
                or (plot_key not in self.data[sampleid])
            ):
                self.data[sampleid][plot_key] = self.genomeplotter.prepare_plot_data(
                    self.data[sampleid][raw_key]
                )

        helper('depthratio_upscaled')
        helper('baf')
        if use_merged_segment:
            helper('merged_segment')
        else:
            helper('depthratio_segment')
            helper('baf_segment')

    def make_plotdata_after_cp(self, sampleid, use_saved_plotdata):
        raw_plot_map = {
            'depthratio_upscaled': 'depthratio_raw_plotdata',
            'baf': 'baf_raw_plotdata',
            'merged_segment': 'merged_segment_plotdata',
        }
        def helper(raw_key):
            plot_key = raw_plot_map[raw_key]
            if (
                (not use_saved_plotdata) 
                or (plot_key not in self.data[sampleid])
            ):
                self.data[sampleid][plot_key] = self.genomeplotter.prepare_plot_data(
                    self.data[sampleid][raw_key]
                )

        helper('depthratio_upscaled')
        helper('baf')
        helper('merged_segment')

    def draw_ax_common(self, ax, n_xlabel):
        self.genomeplotter.set_xlim(ax)
        self.genomeplotter.draw_chrom_borders(ax)
        self.draw_spines_grids(ax, ax.get_yticks())
        self.draw_centromeres(ax)
        #self.draw_centromeres_type2(ax)
        if n_xlabel is not None:
            self.genomeplotter.draw_genomecoord_labels(ax, n=n_xlabel)
        else:
            ax.set_xticks([])

    def draw_depth_ax(
        self,
        sampleid,
        ax,
        n_xlabel,
        is_tumor,
        is_rawdepth,
        dot_kwargs,
        ymax=None,
        binsize=10000,
        add_color=True,
    ):
        # handle kwargs
        dot_kwargs = (
            {'color': 'black', 'markersize': 0.3, 'alpha': 0.05}
            | dot_kwargs
        )

        # set params
        datakey = ('tumor_depth' if is_tumor else 'normal_depth')
        y_colname = ('mean_depth' if is_rawdepth else 'sequenza_style_norm_mean_depth')

        # prepare data
        LOGGER_INFO.info(f'Beginning increasing depth data bin size')

        relevant_chroms = [
            x for x in self.genomeplotter.cconv.totalregion_df['Chromosome']
            if not x.startswith('-')
        ]
        original_df = self.data[sampleid][datakey]
        input_df = original_df.loc[
            original_df['Chromosome'].isin(relevant_chroms), :
        ]
        if binsize is None:
            upscaled_depth_df = input_df.copy()
        else:
            upscaled_depth_df = cnvmisc.upsize_depth_df_bin(
                input_df, 
                size=binsize, 
                refver=self.refver,
            )

        LOGGER_INFO.info(f'Finished increasing depth data bin size')

        # draw data
        if add_color:
            colors = np.repeat('blue', upscaled_depth_df.shape[0])
            colors[np.where(upscaled_depth_df['excluded'])[0]] = 'red'
            upscaled_depth_df['color'] = colors

            dot_kwargs['s'] = dot_kwargs['markersize']
            del dot_kwargs['markersize']
            del dot_kwargs['color']

            self.genomeplotter.draw_dots_scatter(
                ax, 
                df=upscaled_depth_df, 
                y_colname=y_colname, 
                color_colname='color',
                plot_kwargs=dot_kwargs,
            )
        else:
            self.genomeplotter.draw_dots(
                ax, 
                df=upscaled_depth_df, 
                y_colname=y_colname, 
                plot_kwargs=dot_kwargs,
            )

        # set axes attributes
        ylabel_1 = ('tumor' if is_tumor else 'normal')
        ylabel_2 = ('raw' if is_rawdepth else 'normalized')
        ylabel = f'{ylabel_1} sample {ylabel_2} depth'
        ax.set_ylabel(ylabel)

        if ymax is None:
            ymax = np.nanmean(self.data[sampleid][datakey][y_colname]) * 2
        ax.set_ylim(-ymax * 0.1, ymax)
        roundnum = (1 if is_rawdepth else 2)
        yticks = np.round(np.linspace(0, ymax, 10), roundnum)

        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)

        self.draw_ax_common(ax, n_xlabel)

    def draw_depthratio_ax(
        self, 
        sampleid, 
        ax, 
        use_merged_segment, 
        draw_predicted,
        n_xlabel,
        dot_kwargs,
        line_segmean_kwargs,
        line_predict_kwargs,
        ymax=None,
        draw_depthratio_peaks=True,
        depthratio_peaks_kwargs=dict(),
    ):
        # handle kwargs
        #n_dots = self.data[sampleid]['depthratio_raw_plotdata'].shape[0]
        #dot_alpha = n_dots * DOT_ALPHA_CONST

        dot_kwargs = (
            {'color': 'black', 'markersize': 0.3, 'alpha': 0.01}
            | dot_kwargs
        )
        line_segmean_kwargs = (
            {'color': 'tab:blue', 'linewidth': 2, 'alpha': 1}
            | line_segmean_kwargs
        )
        line_predict_kwargs = (
            {'color': 'tab:red', 'linewidth': 4, 'alpha': 0.4}
            | line_predict_kwargs
        )

        # draw data
        self.genomeplotter.draw_dots(
            ax, 
            df_plotdata=self.data[sampleid]['depthratio_raw_plotdata'], 
            y_colname='depthratio_raw', 
            plot_kwargs=dot_kwargs,
        )

        if use_merged_segment:
            segment_df_plotdata = self.data[sampleid]['merged_segment_plotdata']
        else:
            segment_df_plotdata = self.data[sampleid]['depthratio_segment_plotdata']

        self.genomeplotter.draw_hlines(
            ax, 
            df_plotdata=segment_df_plotdata,
            y_colname='depthratio_segment_mean', 
            plot_kwargs=line_segmean_kwargs,
        )

        if draw_predicted:
            if 'depthratio_predicted' in segment_df_plotdata.columns:
                self.genomeplotter.draw_hlines(
                    ax, 
                    df_plotdata=segment_df_plotdata,
                    y_colname='depthratio_predicted', 
                    plot_kwargs=line_predict_kwargs,
                )

        if draw_depthratio_peaks:
            segments_gr = pr.PyRanges(self.data[sampleid]['merged_segment'])
            peak_xs = cnvmisc.get_depthratio_peaks(
                segments_gr.depthratio_segment_mean, 
                lengths=segments_gr.lengths(), 
                limits=(0, 2), 
                step=0.01, 
                peak_cutoff=1e8,
            )

        # set axes attributes
        ax.set_ylabel('depth ratio')

        if ymax is None:
            ymax = np.nanmean(self.data[sampleid]['depthratio_upscaled']['depthratio_raw']) * 2
            ax.set_ylim(-ymax * 0.1, ymax)
            yticks = np.round(np.arange(0, ymax, 0.25), 2)
        else:
            ax.set_ylim(-0.1, ymax)
            yticks = np.round(np.linspace(0, ymax, 8), 2)

        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)

        self.draw_ax_common(ax, n_xlabel)

    def draw_depthratio_hist_ax(
        self,
        sampleid,
        ax,
        use_merged_segment, 
        depth_ylim,
        rm_haploid,
        peak_threshold=1,
        annotate_kwargs=dict(),
        plot_kwargs=dict(),
    ):
        # handle kwargs
        annotate_kwargs = (
            dict(
                va='center', ha='left', size=8,
                arrowprops=dict(arrowstyle='-', color='black', linewidth=0.5),
            )
            | annotate_kwargs
        )
        plot_kwargs = (
            #dict(linewidth=0.8)
            dict(facecolor='tab:blue', alpha=1)
            | plot_kwargs
        )

        if use_merged_segment:
            segment_df = self.data[sampleid]['merged_segment']
        else:
            segment_df = self.data[sampleid]['depthratio_segment']

        if rm_haploid:
            depthratio_list = segment_df.loc[
                segment_df['CNn'] == 2, :
            ]['depthratio_segment_mean']
        else:
            depthratio_list = segment_df['depthratio_segment_mean']
            
        weights = (segment_df['End'] - segment_df['Start'])
        peak_indexes, histbins, hist, bin_midpoints = cnvmisc.get_hist_peaks(
            depthratio_list, 
            weights=weights,
            bins=np.arange(0, depthratio_list.max(), 0.01),
            threshold=peak_threshold,
        )
        peak_depthratios = bin_midpoints[peak_indexes]

        # draw data
        #ax.plot(hist, bin_midpoints, **plot_kwargs)
        ax.fill_betweenx(bin_midpoints, hist, **plot_kwargs)
        ax.set_xlim(0, ax.get_xlim()[1])
        ax.set_ylim(*depth_ylim)
        ytext_list = np.linspace(*ax.get_ylim(), len(peak_depthratios))
        for y, ytext in zip(peak_depthratios, ytext_list):
            ax.axhline(y, color='black', linewidth=0.5)
            ax.annotate(
                round(y, 3), 
                (ax.get_xlim()[1], y),
                xytext=(ax.get_xlim()[1] * 1.1, ytext),
                annotation_clip=False,
                **annotate_kwargs,
            )
        ax.set_yticks([])

        return peak_depthratios

    def draw_baf_ax(
        self, 
        sampleid, 
        ax, 
        use_merged_segment, 
        draw_predicted,
        n_xlabel,
        dot_kwargs,
        line_segmean_kwargs,
        line_corr_segmean_kwargs,
        line_predict_kwargs,
    ):
        # handle kwargs
        dot_kwargs = (
            {'color': 'black', 'markersize': 0.3, 'alpha': 0.01}
            | dot_kwargs
        )
        line_segmean_kwargs = (
            {'color': 'tab:blue', 'linewidth': 2, 'alpha': 0.8}
            | line_segmean_kwargs
        )
        line_corr_segmean_kwargs = (
            {'color': 'tab:green', 'linewidth': 2, 'alpha': 0.8}
            | line_corr_segmean_kwargs
        )
        line_predict_kwargs = (
            {'color': 'tab:red', 'linewidth': 4, 'alpha': 0.4}
            | line_predict_kwargs
        )

        # draw data
        self.genomeplotter.draw_dots(
            ax, 
            df_plotdata=self.data[sampleid]['baf_raw_plotdata'], 
            y_colname='baf_raw_tumor', 
            plot_kwargs=dot_kwargs,
        )

        if use_merged_segment:
            segment_df_plotdata = self.data[sampleid]['merged_segment_plotdata']
        else:
            segment_df_plotdata = self.data[sampleid]['baf_segment_plotdata']

        self.genomeplotter.draw_hlines(
            ax, 
            df_plotdata=segment_df_plotdata,
            y_colname='baf_segment_mean', 
            plot_kwargs=line_segmean_kwargs,
        )
        #if 'corrected_baf_segment_mean' in segment_df_plotdata.columns:
        self.genomeplotter.draw_hlines(
            ax, 
            df_plotdata=segment_df_plotdata,
            y_colname='corrected_baf_segment_mean', 
            plot_kwargs=line_corr_segmean_kwargs,
        )

        if draw_predicted:
            if 'baf_predicted' in segment_df_plotdata.columns:
                self.genomeplotter.draw_hlines(
                    ax, 
                    df_plotdata=segment_df_plotdata,
                    y_colname='baf_predicted', 
                    plot_kwargs=line_predict_kwargs,
                )

        # set axes attributes
        ax.set_ylabel('baf')
        ax.set_ylim(-0.6 * 0.1, 0.6)

        yticks = np.round(np.arange(0, 0.6, 0.1), 1)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)

        self.draw_ax_common(ax, n_xlabel)

    def draw_CN_ax(
        self, 
        sampleid, 
        ax,
        n_xlabel,
        line_A_kwargs,
        line_B_kwargs,
        ymax=None,
    ):
        # handle kwargs
        line_A_kwargs = (
            {'color': 'red'}
            | line_A_kwargs
        )
        line_B_kwargs = (
            {'color': 'blue'}
            | line_B_kwargs
        )

        # draw data
        self.genomeplotter.draw_hlines(
            ax, 
            df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
            y_colname='A', 
            offset=0.1,
            plot_kwargs=line_A_kwargs,
        )
        self.genomeplotter.draw_hlines(
            ax, 
            df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
            y_colname='B', 
            offset=-0.1,
            plot_kwargs=line_B_kwargs,
        )

        # set axes attributes
        ax.set_ylabel('CN')

        if ymax is None:
            weights = (
                self.data[sampleid]['merged_segment']['End'] 
                - self.data[sampleid]['merged_segment']['Start']
            )
            CNts = self.data[sampleid]['merged_segment']['CNt']
            ymax = int(common.nanaverage(CNts, weights)) * 2

            ax.set_ylim(-0.5, ymax)
            yticks = np.arange(0, ymax, 1)
        else:
            ax.set_ylim(-0.5, ymax)
            #yticks = np.linspace(0, ymax, 8).astype(int)
            yticks = np.arange(0, ymax, 1).astype(int)

        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)

        self.draw_ax_common(ax, n_xlabel)

    #def plot_final_wgs(self, sampleid):
    #    return self.plot(sampleid, figsize=(30, 13), draw_invalid_regions=False)

    def plot_base_wrapper(
        self,
        basemethod, 

        region_chroms=None, 
        region_start0s=None, 
        region_end0s=None, 
        weights=None,
        region_gaps=None,

        **kwargs,
    ):
        if (
            (region_gaps is not None)
            or (region_chroms is not None)
        ):
            # make new genomeplotter
            if region_chroms is None:
                new_region_df = self.genomeplotter.region_df.copy()
            else:
                new_region_df = GenomePlotter.make_new_region_df(
                    self.refver, 
                    region_chroms, 
                    region_start0s, 
                    region_end0s, 
                    weights,
                )

            new_genomeplotter = GenomePlotter(
                self.refver, 
                region_df=new_region_df,
                region_gaps=region_gaps,
            )

            # switch and run
            old_genomeplotter = self.genomeplotter
            self.genomeplotter = new_genomeplotter

            fig, axd = basemethod(**kwargs)

            self.genomeplotter = old_genomeplotter
        else:
            fig, axd = basemethod(**kwargs)

        return fig, axd

    def plot_woCN_base(
        self, 
        sampleid, 
        figsize, 
        hspace,
        draw_invalid_regions, 
        use_saved_plotdata,
        use_merged_segment,

        draw_depthratio_hist,
        rm_haploid,
        depthratio_hist_threshold,

        draw_depth,
        is_rawdepth,
        depth_binsize,

        n_xlabel,
        depthratio_ymax,
        depth_ymax,
        
        # plot properties

        depthratio_dot_kwargs,
        depthratio_line_segmean_kwargs,
        depthratio_line_predict_kwargs,

        depthratio_hist_annotate_kwargs,
        depthratio_hist_plot_kwargs,

        depth_dot_kwargs,

        baf_dot_kwargs,
        baf_line_segmean_kwargs,
        baf_line_corr_segmean_kwargs,
        baf_line_predict_kwargs,
    ):
        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
        self.make_plotdata_before_cp(sampleid, use_saved_plotdata, use_merged_segment)
        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')

        gridspec_kw = dict(hspace=hspace)
        if draw_depthratio_hist:
            mosaic = [
                ['baf', 'empty'], 
                ['depthratio', 'depthratio_hist'],
            ]
            if draw_depth:
                mosaic.append(['depth_tumor', 'empty_lower'])
                mosaic.append(['depth_normal', 'empty_lower'])

            gridspec_kw.update(dict(width_ratios=[1, 0.1], wspace=0.02))
        else:
            mosaic = [
                ['baf'], 
                ['depthratio'],
            ]
            if draw_depth:
                mosaic.append(['depth_tumor'])
                mosaic.append(['depth_normal'])

        fig, axd = plt.subplot_mosaic(
            mosaic,
            figsize=figsize,
            gridspec_kw=gridspec_kw,
        )
        fig.suptitle(
            f'sample_id={sampleid}, is_female={self.data[sampleid]["is_female"]}',
            fontsize=20,
        )

        self.draw_baf_ax(
            sampleid, 
            axd['baf'], 
            use_merged_segment=use_merged_segment,
            draw_predicted=False,
            n_xlabel=n_xlabel,

            dot_kwargs=baf_dot_kwargs,
            line_segmean_kwargs=baf_line_segmean_kwargs,
            line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,
            line_predict_kwargs=baf_line_predict_kwargs,
        )
        self.draw_depthratio_ax(
            sampleid, 
            axd['depthratio'], 
            use_merged_segment=use_merged_segment,
            draw_predicted=False,
            n_xlabel=n_xlabel,
            dot_kwargs=depthratio_dot_kwargs,
            line_segmean_kwargs=depthratio_line_segmean_kwargs,
            line_predict_kwargs=depthratio_line_predict_kwargs,
            ymax=depthratio_ymax,
        )

        if draw_depthratio_hist:
            peak_depthratios = self.draw_depthratio_hist_ax(
                sampleid,
                axd['depthratio_hist'],
                use_merged_segment=use_merged_segment, 
                depth_ylim=axd['depthratio'].get_ylim(),
                rm_haploid=rm_haploid,
                peak_threshold=depthratio_hist_threshold,
                annotate_kwargs=depthratio_hist_annotate_kwargs,
                plot_kwargs=depthratio_hist_plot_kwargs,
            )

            for y in peak_depthratios:
                axd['depthratio'].axhline(y, color='orange', linewidth=1, alpha=0.6)

            axd['empty'].axis('off')

        if draw_depth:
            LOGGER_INFO.info(f'Beginning tumor depth data processing')
            self.draw_depth_ax(
                sampleid,
                axd['depth_tumor'],
                n_xlabel=n_xlabel,
                is_tumor=True,
                is_rawdepth=is_rawdepth,
                dot_kwargs=depth_dot_kwargs,
                ymax=depth_ymax,
                binsize=depth_binsize,
            )
            LOGGER_INFO.info(f'Finished tumor depth data processing')

            LOGGER_INFO.info(f'Beginning normal depth data processing')
            self.draw_depth_ax(
                sampleid,
                axd['depth_normal'],
                n_xlabel=n_xlabel,
                is_tumor=False,
                is_rawdepth=is_rawdepth,
                dot_kwargs=depth_dot_kwargs,
                ymax=depth_ymax,
                binsize=depth_binsize,
            )
            LOGGER_INFO.info(f'Finished normal depth data processing')

            if draw_depthratio_hist:
                axd['empty_lower'].axis('off')

        return fig, axd

    def plot_woCN(
        self, 
        sampleid, 
        figsize=None, 
        hspace=None,
        draw_invalid_regions=False,
        use_saved_plotdata=False,
        n_xlabel=None,
        depthratio_ymax=None,
        depth_ymax=None,
        use_merged_segment=True,

        draw_depthratio_hist=True,
        rm_haploid=True,
        depthratio_hist_threshold=3,

        draw_depth=False,
        is_rawdepth=True,
        depth_binsize=10000,

        depthratio_dot_kwargs=dict(),
        depthratio_line_segmean_kwargs=dict(),
        depthratio_line_predict_kwargs=dict(),

        depthratio_hist_annotate_kwargs=dict(),
        depthratio_hist_plot_kwargs=dict(),

        depth_dot_kwargs=dict(),

        baf_dot_kwargs=dict(),
        baf_line_segmean_kwargs=dict(),
        baf_line_corr_segmean_kwargs=dict(),
        baf_line_predict_kwargs=dict(),

        region_chroms=None, 
        region_start0s=None, 
        region_end0s=None, 
        weights=None,
        region_gaps=None,
    ):
        if figsize is None:
            if draw_depth:
                figsize = (30, 15)
            else:
                figsize = (30, 9)

        fig, axd = self.plot_base_wrapper(
            self.plot_woCN_base, 

            region_chroms=region_chroms, 
            region_start0s=region_start0s, 
            region_end0s=region_end0s, 
            weights=weights,
            region_gaps=region_gaps,

            sampleid=sampleid, 
            figsize=figsize, 
            hspace=hspace,
            draw_invalid_regions=draw_invalid_regions,
            use_saved_plotdata=use_saved_plotdata,
            n_xlabel=n_xlabel,
            depthratio_ymax=depthratio_ymax,
            depth_ymax=depth_ymax,
            use_merged_segment=use_merged_segment,

            draw_depthratio_hist=draw_depthratio_hist,
            rm_haploid=rm_haploid,
            depthratio_hist_threshold=depthratio_hist_threshold,

            draw_depth=draw_depth,
            is_rawdepth=is_rawdepth,
            depth_binsize=depth_binsize,

            depthratio_dot_kwargs=depthratio_dot_kwargs,
            depthratio_line_segmean_kwargs=depthratio_line_segmean_kwargs,
            depthratio_line_predict_kwargs=depthratio_line_predict_kwargs,
            depthratio_hist_annotate_kwargs=depthratio_hist_annotate_kwargs,
            depthratio_hist_plot_kwargs=depthratio_hist_plot_kwargs,
            depth_dot_kwargs=depth_dot_kwargs,
            baf_dot_kwargs=baf_dot_kwargs,
            baf_line_segmean_kwargs=baf_line_segmean_kwargs,
            baf_line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,
            baf_line_predict_kwargs=baf_line_predict_kwargs,
        )
        return fig, axd

    def plot_final_base(
        self, 
        sampleid, 
        cellularity,
        ploidy,

        figsize, 
        hspace,
        draw_invalid_regions,
        use_saved_plotdata,

        n_xlabel,
        depthratio_ymax,
        CN_ymax,

        depthratio_dot_kwargs,
        depthratio_line_segmean_kwargs,
        depthratio_line_predict_kwargs,

        baf_dot_kwargs,
        baf_line_segmean_kwargs,
        baf_line_corr_segmean_kwargs,
        baf_line_predict_kwargs,

        CN_line_A_kwargs,
        CN_line_B_kwargs,
    ):
        LOGGER_INFO.info(f'Beginning segment dataframe modification')
        self.add_CN_pred_to_segment(sampleid, cellularity, ploidy)
        LOGGER_INFO.info(f'Finished segment dataframe modification')

        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
        self.make_plotdata_after_cp(sampleid, use_saved_plotdata)
        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')

        fig, axd = plt.subplot_mosaic(
            [
                ['CN',], 
                ['baf',], 
                ['depth',],
            ],
            figsize=figsize,
            gridspec_kw=dict(hspace=hspace),
        )

        fig.suptitle(
            ', '.join([
                f'sample_id={sampleid}',
                f'is_female={self.data[sampleid]["is_female"]}',
                f'cellularity={round(cellularity, 2)}',
                f'ploidy={round(ploidy, 2)}',
            ]),
            fontsize=20,
        )

        self.draw_baf_ax(
            sampleid, 
            axd['baf'], 
            use_merged_segment=True,
            draw_predicted=True,
            n_xlabel=n_xlabel,

            dot_kwargs=baf_dot_kwargs,
            line_segmean_kwargs=baf_line_segmean_kwargs,
            line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,
            line_predict_kwargs=baf_line_predict_kwargs,
        )
        self.draw_depthratio_ax(
            sampleid, 
            axd['depth'], 
            use_merged_segment=True,
            draw_predicted=True,
            n_xlabel=n_xlabel,
            dot_kwargs=depthratio_dot_kwargs,
            line_segmean_kwargs=depthratio_line_segmean_kwargs,
            line_predict_kwargs=depthratio_line_predict_kwargs,
            ymax=depthratio_ymax,
        )

        self.draw_CN_ax(
            sampleid, 
            axd['CN'],
            n_xlabel=n_xlabel,
            line_A_kwargs=CN_line_A_kwargs,
            line_B_kwargs=CN_line_B_kwargs,
            ymax=CN_ymax,
        )

        return fig, axd

    def plot_final(
        self, 
        sampleid, 
        cellularity,
        ploidy,

        figsize=(30, 13.5), 
        hspace=None,
        draw_invalid_regions=False,
        use_saved_plotdata=False,
        n_xlabel=None,
        depthratio_ymax=None,
        CN_ymax=None,

        depthratio_dot_kwargs=dict(),
        depthratio_line_segmean_kwargs=dict(),
        depthratio_line_predict_kwargs=dict(),
        baf_dot_kwargs=dict(),
        baf_line_segmean_kwargs=dict(),
        baf_line_corr_segmean_kwargs=dict(),
        baf_line_predict_kwargs=dict(),
        CN_line_A_kwargs=dict(),
        CN_line_B_kwargs=dict(),

        region_chroms=None, 
        region_start0s=None, 
        region_end0s=None, 
        weights=None,
        region_gaps=None,
    ):
        fig, axd = self.plot_base_wrapper(
            self.plot_final_base, 

            region_chroms=region_chroms, 
            region_start0s=region_start0s, 
            region_end0s=region_end0s, 
            weights=weights,
            region_gaps=region_gaps,

            sampleid=sampleid, 
            cellularity=cellularity,
            ploidy=ploidy,
            figsize=figsize, 
            hspace=hspace,
            draw_invalid_regions=draw_invalid_regions,
            use_saved_plotdata=use_saved_plotdata,
            n_xlabel=n_xlabel,
            depthratio_ymax=depthratio_ymax,
            CN_ymax=CN_ymax,

            depthratio_dot_kwargs=depthratio_dot_kwargs,
            depthratio_line_segmean_kwargs=depthratio_line_segmean_kwargs,
            depthratio_line_predict_kwargs=depthratio_line_predict_kwargs,
            baf_dot_kwargs=baf_dot_kwargs,
            baf_line_segmean_kwargs=baf_line_segmean_kwargs,
            baf_line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,
            baf_line_predict_kwargs=baf_line_predict_kwargs,
            CN_line_A_kwargs=CN_line_A_kwargs,
            CN_line_B_kwargs=CN_line_B_kwargs,
        )

        return fig, axd

#    def plot(self, sampleid, figsize=(30, 13), draw_invalid_regions=False):
#        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
#        self.make_plotdata(sampleid)
#        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')
#
#        fig, axd = plt.subplot_mosaic(
#            [
#                ['CN',], 
#                ['baf',], 
#                ['depth',],
#            ],
#            figsize=figsize,
#        )
#        for ax in axd.values():
#            self.genomeplotter.set_xlim(ax)
#
#        # CN
#        self.draw_CN_ax(sampleid, axd['CN'])
#
#        # baf
#        self.draw_baf_ax(sampleid, axd['baf'])
#
#        # depthratio
#        self.draw_depthratio_ax(sampleid, axd['depth'])
#
#        # invalid regions
#        if draw_invalid_regions:
#            selector = np.logical_or(
#                self.data[sampleid]['depth_baf_merge']['depthratio_raw'].isna().to_numpy,
#                self.data[sampleid]['depth_baf_merge']['baf_raw'].isna().to_numpy,
#            )
#            self.genomeplotter.draw_bgcolors(
#                axd['depth'], 
#                df=self.data[sampleid]['depth_baf_merge'].loc[selector, :], 
#                colors='red',
#                plot_kwargs=dict(alpha=0.01),
#            )
#
#        # draw centromeric regions
#        self.draw_centromeres(axd.values())
#
#        # return
#        return fig, axd

    #############################
    # helpers of add new sample #
    #############################

    def set_target_region(self, sampleid, mode, target_region_arg):
        if mode == 'wgs':
            target_region_gr = self.find_germline_copyneutral_region(
                sampleid, factors=(0.8, 1.2),
            )
        elif mode == 'panel':
            target_region_gr = cnvmisc.arg_into_gr(target_region_arg).drop()

        # exclude y if female
        if self.data[sampleid]['is_female']:
            y_chrom = common.DEFAULT_CHROMDICTS[self.refver].XY_names[1]
            target_region_gr = target_region_gr[
                target_region_gr.Chromosome != y_chrom
            ]

        # merge
        target_region_gr = target_region_gr.merge()

        # result
        self.data[sampleid]['target_region'] = target_region_gr

    def find_germline_copyneutral_region(self, sampleid, factors=(0.8, 1.2)):
        depth_df = self.data[sampleid]['normal_depth']

        # upscale depth bins for segmentation speed
        upscaled_depth_df = cnvmisc.upsize_depth_df_bin(
            depth_df, size=1000, refver=self.refver,
        )

        # run segmentation
        segdf, _ = rcopynumber.run_rcopynumber(
            depth_df=upscaled_depth_df.rename(columns={'mean_depth': 'depth_raw'}), 
            refver=self.refver, 
            as_gr=False, 
            winsorize=False, 
            verbose=True,
            remove_unassembled_contigs=True,
        )
        self.data[sampleid]['normal_depth_segment'] = segdf

        # select copy-neutral depth segments
        global_mean = np.average(
            upscaled_depth_df['mean_depth'], 
            weights=(upscaled_depth_df['End'] - upscaled_depth_df['Start']),
        )
        selector = segdf['depth_segment_mean'].between(
            global_mean * factors[0], 
            global_mean * factors[1], 
        )
        included_segments = segdf.loc[selector, :]

        target_region_gr = pr.PyRanges(depth_df).drop().overlap(
            pr.PyRanges(included_segments)
        )

        return target_region_gr

    #########
    # depth #
    #########

    def load_bam_for_depth(self, sampleid, bam_path, sampletype):
        assert sampletype in ('normal', 'tumor')
        mosdepth_df = self._run_mosdepth(
            bam_path,
            binsize=self.default_binsize,
            region_df=(
                None 
                if self.data[sampleid]['mode'] == 'wgs' else 
                self.data[sampleid]['target_region']
            ),
        )

        self.data.setdefault(sampleid, dict())
        self.data[sampleid][f'{sampletype}_depth'] = mosdepth_df

    def load_mosdepth_df(
        self, sampleid, mosdepth_df, sampletype, use_default_gc=True, **kwargs,
    ):
        assert sampletype in ('normal', 'tumor')

        self.data.setdefault(sampleid, dict())
        depth_df, gcbin_average_depths = self._postprocess_mosdepth_df(
            mosdepth_df, 
            self.data[sampleid]['mode'],
            sampletype, 
            use_default_gc=use_default_gc, 
            **kwargs,
        )
        self.data[sampleid][f'{sampletype}_depth'] = depth_df
        self.data[sampleid][f'{sampletype}_gcdata'] = gcbin_average_depths

    def _run_mosdepth(self, bam_path, binsize, region_df):
        mosdepth_df, _ = libmosdepth.run_mosdepth(
            bam_path, 
            t=8, 
            use_median=False, 
            region_bed_path=None, 
            region_gr=region_df, 
            window_size=binsize, 
            donot_subset_bam=True,
            as_gr=False, 
            load_perbase=False,
        )

        return mosdepth_df

    def _postprocess_mosdepth_df(
        self, 
        mosdepth_df, 
        mode,
        sampletype, 
        use_default_gc=True, 
        **kwargs,
    ):
        assert sampletype in ('normal', 'tumor')

        # set preset_cutoffs
        #if mode == 'panel':
            #preset_cutoffs = 'panel'
        #elif mode == 'wgs':
            #if sampletype == 'normal':
                #preset_cutoffs = 'normal_wgs'
            #elif sampletype == 'tumor':
                #preset_cutoffs = 'wgs'
        #kwargs['preset_cutoffs'] = preset_cutoffs

        # set gc_df
        if use_default_gc:
            gc_df = libgcfraction.get_gc_df(
                self.refver, self.default_binsize, coords_as_index=True,
            )
            kwargs['gc_df'] = gc_df

        kwargs['as_gr'] = False
        depth_df, gcbin_average_depths = cnvmisc.postprocess_depth_df(
            mosdepth_df, 
            **kwargs,
        )
        return depth_df, gcbin_average_depths

    #######
    # baf #
    #######

    def load_germline_vcf(
        self, 
        sampleid, 
        vcf_path, 
        vcf_sampleid_tumor,
        vcf_sampleid_normal=None,
        nproc=1, 
        logging_lineno=50000,
    ):
        self.data.setdefault(sampleid, dict())

        if vcf_sampleid_normal is None:
            vcf_sampleids = [vcf_sampleid_tumor]
        else:
            vcf_sampleids = [vcf_sampleid_tumor, vcf_sampleid_normal]

        # load vafdf
        vaf_df = variantplus.get_vafdf(
            vcf_path, 
            sampleid=vcf_sampleids, 
            nproc=nproc,
        )
        # rename vaf columns
        vaf_df.rename(
            columns={f'vaf_{vcf_sampleid_tumor}': 'vaf_raw_tumor'}, inplace=True,
        )
        if vcf_sampleid_normal is not None:
            vaf_df.rename(
                columns={f'vaf_{vcf_sampleid_normal}': 'vaf_raw_normal'}, inplace=True,
            )

        # postprocess
        #vaf_df = vaf_df.loc[vaf_df['vaf_raw'].notna().to_numpy(), :]
        vaf_df.reset_index(drop=True, inplace=True)
        vaf_df['baf_raw_tumor'] = cnvmisc.get_bafs(vaf_df['vaf_raw_tumor'])
        if 'vaf_raw_normal' in vaf_df:
            vaf_df['baf_raw_normal'] = cnvmisc.get_bafs(vaf_df['vaf_raw_normal'])

        vaf_df = vaf_df.loc[
            :, [
                'Chromosome', 'Start', 'End', 
                'vaf_raw_tumor', 
                'baf_raw_tumor', 
                'vaf_raw_normal', 
                'baf_raw_normal', 
            ]
        ]

        self.data[sampleid]['baf'] = vaf_df

    #################
    # other helpers #
    #################

    def set_normal_mean_ploidy(self, sampleid):
        self.data[sampleid]['normal_mean_ploidy'] = cnvmisc.get_normal_mean_ploidy(
            self.refver, 
            self.data[sampleid]['is_female'], 
            self.data[sampleid]['target_region'],
        )

    @staticmethod
    def _add_CNn_to_segment(
        segment_df,
        mode,
        refver,
        is_female,
        target_region=None,
    ):
        if mode == 'wgs':
            segment_df = rcopynumber.add_CNn_to_wgs_segment_gr(
                segment_df, refver, is_female, as_gr=False,
            )
        elif mode == 'panel':
            assert target_region is not None
            segment_df = rcopynumber.add_CNn_to_targetseq_segment_gr(
                segment_df, target_region, refver, is_female, as_gr=False,
            )
        return segment_df


###############################################################


def make_targetseq_cnvplot(data_df, draw_invalid_regions=False):
    # sanity check
    required_cols = {
        'CNt', 'B', 
        'baf_segment_mean', 'baf_predicted', 'baf_raw',
        'depthratio_segment_mean', 'depthratio_predicted', 'depthratio_raw',
    }
    assert required_cols.issubset(data_df.columns)

    # setup
    if isinstance(data_df, pr.PyRanges):
        data_df = data_df.df

    gplotter = GenomePlotter(data_df)
    fig, axd = plt.subplot_mosaic(
        [
            ['CN',], 
            ['baf',], 
            ['depth',],
        ],
        figsize=(30, 13),
    )
    for ax in axd.values():
        gplotter.set_xlim(ax)

    # CN
    axd['CN'].set_ylabel('CN')

    gplotter.draw_hlines(
        axd['CN'], df=data_df, y_colname='CNt', offset=0.1,
        plot_kwargs={'color': 'black'},
    )
    gplotter.draw_hlines(
        axd['CN'], df=data_df, y_colname='B', 
        plot_kwargs={'color': 'blue'},
    )
    gplotter.draw_chrom_borders(axd['CN'])

    # baf
    axd['baf'].set_ylabel('baf')
    axd['baf'].set_ylim(0, 0.6)

    gplotter.draw_dots(
        axd['baf'], df=data_df, y_colname='baf_raw', 
        plot_kwargs={'color': 'gray', 'markersize': 1.5, 'alpha': 0.5},
    )
    gplotter.draw_hlines(
        axd['baf'], df=data_df, y_colname='baf_segment_mean', 
        plot_kwargs={'color': 'black', 'linewidth': 2, 'alpha': 0.3},
    )
    gplotter.draw_hlines(
        axd['baf'], df=data_df, y_colname='baf_predicted', 
        plot_kwargs={'color': 'green', 'linewidth': 2, 'alpha': 0.3},
    )
    gplotter.draw_chrom_borders(axd['baf'])

    # depthratio
    axd['depth'].set_ylabel('depth ratio')
    axd['depth'].set_ylim(
        0, 
        data_df['depthratio_raw'].max() * 1.1
    )

    gplotter.draw_dots(
        axd['depth'], df=data_df, y_colname='depthratio_raw', 
        plot_kwargs={'color': 'gray', 'markersize': 1.5, 'alpha': 0.5},
    )
    gplotter.draw_hlines(
        axd['depth'], df=data_df, y_colname='depthratio_segment_mean', 
        plot_kwargs={'color': 'black', 'linewidth': 2, 'alpha': 0.5},
    )
    gplotter.draw_hlines(
        axd['depth'], df=data_df, y_colname='depthratio_predicted', 
        plot_kwargs={'color': 'green', 'linewidth': 2, 'alpha': 0.5},
    )

    if draw_invalid_regions:
        selector = np.logical_or(
            data_df['depthratio_raw'].isna().to_numpy,
            data_df['baf_raw'].isna().to_numpy,
        )
        gplotter.draw_bgcolors(
            axd['depth'], df=data_df.loc[selector, :], 
            plot_kwargs=dict(color='red', alpha=0.01)
        )

    gplotter.draw_chrom_borders(axd['depth'])

    return fig, axd


def draw_targetseq_cnvplot_from_data_precp(
    tumor_depth_df, 
    normal_depth_df,
    germline_df,

    region_gr,
    refver, 
    is_female, 
):
    assert isinstance(tumor_depth_df, pd.DataFrame)
    assert isinstance(normal_depth_df, pd.DataFrame)
    assert isinstance(germline_df, pd.DataFrame)

    if 'baf_raw' not in germline_df.columns:
        germline_df['baf_raw'] = cnvmisc.get_bafs(germline_df['vaf'])
    germline_gr = pr.PyRanges(germline_df)

    # make depth ratio
    depthratio_gr = cnvmisc.make_depth_ratio(tumor_depth_df, normal_depth_df, as_gr=True)
    depthratio_df = depthratio_gr.df

    # run R copynumber
    segment_gr, depth_baf_gr = rcopynumber.run_rcopynumber_unified(
        depthratio_df.rename(columns={'depth_ratio_sequenzastyle': 'depth_raw'}),
        baf_df=germline_df, 
        refver=refver, 
        as_gr=True, 
        compact=True, 
        winsorize=False,
        gamma=30, 
        kmin=1,
    )
    segment_gr = segment_gr[[]]

    # postprocess segment df
    segment_gr = rcopynumber.add_CNn_to_targetseq_segment_gr(segment_gr, region_gr, refver=refver, is_female=is_female)
    segment_gr = pyranges_helper.join(
        segment_gr, 
        depthratio_gr[['depth_ratio_sequenzastyle']], 
        how='left', merge='mean', find_nearest=False, as_gr=True,
    )
    segment_gr = pyranges_helper.join(
        segment_gr, 
        germline_gr[['baf_raw']], 
        how='left', merge='mean', find_nearest=False, as_gr=True,
    )
    segment_df = segment_gr.df
    segment_df.rename(
        {'depth_ratio_sequenzastyle': 'depthratio_segment_mean', 'baf_raw': 'baf_segment_mean'},
        inplace=True,
        axis=1,
    )

    return segment_df, depth_baf_gr


def draw_targetseq_cnvplot_from_data_choosecp(
    segment_df, 
    depth_baf_gr,

    region_gr,
    refver,
    is_female,

    CNt_weight=cnvmisc.DEFAULT_CNT_WEIGHT,
    nproc=None,
):
    cpscore_dict = cnvmisc.get_cp_score_dict(segment_df, refver=refver, is_female=is_female, target_region_gr=region_gr, CNt_weight=CNt_weight, nproc=nproc)
    peak_values, dfs = cnvmisc.get_peak_info(cpscore_dict)

    cnvmisc.show_heatmap_peaks_new(dfs, figsize=(20, 20))

    return peak_values, dfs, cpscore_dict


def draw_targetseq_cnvplot_from_data_postcp(
    cellularity,
    ploidy,
    segment_df, 
    depth_baf_gr,

    region_gr,
    refver,
    is_female,

    draw_invalid_regions=False,
):
    normal_mean_ploidy = cnvmisc.get_normal_mean_ploidy(refver, is_female, target_region_gr=region_gr)

    # postprocess segment df
    segment_df = cnvmisc.add_CNt_to_segment(
        segment_df, 
        cellularity=cellularity, 
        tumor_ploidy=ploidy, 
        is_female=is_female, 
        normal_ploidy=normal_mean_ploidy,
    )
    segment_df = cnvmisc.add_theoreticals_to_segment(
        segment_df, 
        cellularity=cellularity, 
        tumor_ploidy=ploidy, 
        normal_ploidy=normal_mean_ploidy, 
        is_female=is_female,
    )

    # make a single merged plot data df
    depth_baf_gr.depthratio_raw = depth_baf_gr.depth_raw
    depth_baf_gr = cnvmisc.annotate_region_with_segment(depth_baf_gr, pr.PyRanges(segment_df), as_gr=True)

    # draw plot
    fig, axd = make_targetseq_cnvplot(depth_baf_gr, draw_invalid_regions=draw_invalid_regions)

    return fig, axd



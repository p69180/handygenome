import functools
import itertools
import multiprocessing
import operator

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


LOGGER_INFO = workflow.get_debugging_logger(verbose=False)
LOGGER_DEBUG = workflow.get_debugging_logger(verbose=True)


class CoordConverter:
    def __init__(self, refver, df=None, region_gaps=True):
        """Args:
            df: pandas.DataFrame object with mandatory columns 'Chromosome', 
                'Start', 'End'. Each row represents a genomic interval, which 
                will be drawn separately in the plot. The rows need not be 
                sorted by chromosome or position. The order of rows of the 
                input DataFrame is respected.
        """
        self.refver = refver
        self.chromdict = common.DEFAULT_CHROMDICTS[refver]
        if df is None:
            df = self.chromdict.to_gr(assembled_only=True, as_gr=False)
        else:
            cnvmisc.genome_df_sanitycheck(df)
            df = cnvmisc.arg_into_df(df)

        df = self.handle_df(df, region_gaps)
        self.set_params(df)
        #self.is_allregion = (df is None)

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

        if region_gaps:
            total_length = (df['End'] - df['Start']).sum()
            gap_length = int(total_length / 100)
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

    def check_is_allregion(self):
        allregion_gr = self.chromdict.to_gr(assembled_only=True, as_gr=True)
        isec_gr = allregion_gr.intersect(self.totalregion_gr)
        return (isec_gr[[]].sort().df == allregion_gr[[]].sort().df).all(axis=None)

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

#    def genomic_to_plot_new(self, chrom, pos0_list):
#        if chrom not in self.chromwise_params.keys():
#            raise Exception(f'Input "chrom" argument is not included in the plotting region.')
#
#        pos0_list = np.array(pos0_list)[:, np.newaxis]
#        params = self.chromwise_params[chrom]
#
#        contains = np.logical_and(
#            (pos0_list >= params['start0']), (pos0_list < params['end0'])
#        )
#        pos0s_indexes, intv_indexes = np.where(contains)
#            # np.ndarray composed of the indexes of the containing intervals
#            # identical intervals can appear many times
#        within_region_offsets = (
#            params['plot_region_length'][intv_indexes]
#            * (
#                (pos0_list[pos0s_indexes, 0] - params['start0'][intv_indexes]) 
#                / params['raw_region_length'][intv_indexes]
#            )
#        )
#        plot_coords = params['plot_region_start_offset'][intv_indexes] + within_region_offsets
#        return plot_coords, pos0s_indexes, intv_indexes

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
        result_pos0s = (totalregion_subdf['Start'] + offsets).to_numpy().astype(int)
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
        region_df=None, chroms=None, start0s=None, end0s=None, weights=1,
        region_gaps=True,
    ):
        self.refver = refver

        if region_df is None:
            if chroms is None:
                region_df = None
            else:
                region_df = self.make_new_region_df(
                    self.refver, chroms, start0s=start0s, end0s=end0s, weights=weights
                )

        self.cconv = CoordConverter(refver=refver, df=region_df, region_gaps=region_gaps)

    def check_is_allregion(self):
        return self.cconv.check_is_allregion()

    #@property
    #def refver(self):
    #    return self.cconv.refver

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
        default_plot_kwargs = {
            'color': 'black',
            'marker': 'o',
            'linestyle': '',
        }
        default_plot_kwargs.update(plot_kwargs)

        if df_plotdata is None:
            df_plotdata = self.prepare_plot_data(df)

        xs = (df_plotdata['plot_start0s'] + (df_plotdata['plot_end0s'] - 1)) / 2
        ys = df_plotdata[y_colname].to_numpy()
        ax.plot(xs, ys, **default_plot_kwargs)

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
        colors = common.arg_into_list(colors)
        if len(colors) == 1:
            match_original = False
            default_plot_kwargs['color'] = colors[0]
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
        ordered_chroms = [x[0] for x in itertools.groupby(isec_gr.Chromosome)]

        result_start0s = list()
        result_end0s = list()
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
    def make_new_region_df(refver, chroms, start0s=None, end0s=None, weights=1):
        assert (start0s is None) == (end0s is None)
        # weights
        weights = common.arg_into_list(weights)
        if not all(str(x).isdecimal() for x in weights):
            raise Exception(f'"weights" argument must be all intergers')
        # chroms
        chroms = common.arg_into_list(chroms)
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


class CNVPlotter(GenomePlotter):
    def __init__(self, refver, region_df=None):
        super().__init__(refver=refver, region_df=region_df)
        self.data = dict()
        self.default_binsize = 100

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

        # main
        self.data[sampleid] = dict()
        self.data[sampleid]['is_female'] = is_female
        self.data[sampleid]['mode'] = mode
        self.data[sampleid]['target_region'] = target_region
        self.set_normal_mean_ploidy(sampleid)

        LOGGER_INFO.info('Loading normal depth')
        helper(self, normal_bam_path, normal_depth_path, sampleid, 'normal', mosdepth_postprocess_kwargs)
        LOGGER_INFO.info('Loading tumor depth')
        helper(self, tumor_bam_path, tumor_depth_path, sampleid, 'tumor', mosdepth_postprocess_kwargs)
        LOGGER_INFO.info('Loading tumor germline vcf')
        self.load_tumor_germline_vcf(
            sampleid, tumor_vcf_path, tumor_vcf_sampleid, logging_lineno=50000,
            nproc=vcfload_nproc,
        )

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
        self.data[sampleid]['depthratio'] = depthratio_df

    def upscale_depthratio(self, sampleid, binsize=10000):
        self.data[sampleid]['depthratio_upscaled'] = cnvmisc.upsize_depth_df_bin(
            self.data[sampleid]['depthratio'], 
            size=binsize, 
            refver=self.refver,
        )

    def set_depth_baf_merge(
        self, 
        sampleid, 
        nproc=1,
    ):
        depth_df = self.data[sampleid]['depthratio_upscaled'].loc[
            :, ['Chromosome', 'Start', 'End', 'depth_raw']
        ]
        depth_df = cnvmisc.remove_unassembled_contigs(depth_df)

        baf_df = self.data[sampleid]['tumor_baf'].loc[
            :, ['Chromosome', 'Start', 'End', 'baf_raw']
        ]
        baf_df = cnvmisc.remove_unassembled_contigs(baf_df)

        merged_df, merged_df_wona = rcopynumber.make_merged_df(
            depth_df=depth_df, 
            baf_df=baf_df, 
            cytoband_df=None, 
            refver=self.refver, 
            logger=LOGGER_INFO, 
            nproc=nproc,
        )
        self.data[sampleid]['depth_baf_merge'] = merged_df
        self.data[sampleid]['depth_baf_merge_wona'] = merged_df_wona

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
        bafcorrector=None,
    ):
        gamma, kmin = self.handle_gamma_kmin_args(gamma, kmin, self.data[sampleid]['mode'])

        input_df = self.data[sampleid]['tumor_baf'].loc[:, ['Chromosome', 'Start', 'End']]
        input_df['depth_raw'] = self.data[sampleid]['tumor_baf']['baf_raw']

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
        segment_df.rename(columns={'depth_segment_mean': 'baf_segment_mean'}, inplace=True)
        if bafcorrector is not None:
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

        self.data[sampleid]['merged_segment'] = merged_segment

    def make_segment(
        self,
        sampleid,
        winsorize=False,
        gamma=None,
        kmin=None,
        verbose=False,
    ):
        # set params
        gamma, kmin = self.handle_gamma_kmin_args(gamma, kmin, self.data[sampleid]['mode'])

        # run R-copynumber
        if 'depthratio_raw' in self.data[sampleid]['depth_baf_merge'].columns:
            self.data[sampleid]['depth_baf_merge'].rename(
                columns={'depthratio_raw': 'depth_raw'}, inplace=True,
            )

        segment_df, depth_baf_df = rcopynumber.run_rcopynumber_unified(
            merged_df=self.data[sampleid]['depth_baf_merge'], 
            merged_df_wona=self.data[sampleid]['depth_baf_merge_wona'],

            as_gr=False, 
            winsorize=winsorize,
            compact=(True if self.data[sampleid]['mode'] == 'panel' else False), 
            verbose=verbose,
            gamma=gamma,
            kmin=kmin,
        )

        # modify segment_df
        segment_df.rename(
            columns={'depth_segment_mean': 'depthratio_segment_mean'}, inplace=True,
        )

        # add CNn
        segment_df = self._add_CNn_to_segment(
            segment_df,
            self.data[sampleid]['mode'],
            self.refver,
            self.data[sampleid]['is_female'],
            self.data[sampleid]['target_region'],
        )

        # recalc mean depthratio and baf with panel data
        if self.data[sampleid]['mode'] == 'panel':
            segment_df = segment_df.loc[:, ['Chromosome', 'Start', 'End', 'CNn']]

            segment_df = pyranges_helper.join(
                segment_df, 
                self.data[sampleid]['depthratio'].loc[
                    :, ['Chromosome', 'Start', 'End', 'depth_raw']
                ], 
                how='left', merge='mean', find_nearest=False, as_gr=False,
            )
            segment_df.rename(
                columns={'depth_raw': 'depthratio_segment_mean'}, inplace=True,
            )

            segment_df = pyranges_helper.join(
                segment_df, 
                self.data[sampleid]['tumor_baf'].loc[
                    :, ['Chromosome', 'Start', 'End', 'baf_raw']
                ],
                how='left', merge='mean', find_nearest=False, as_gr=True,
            )
            segment_df.rename(
                columns={'baf_raw': 'baf_segment_mean'}, inplace=True,
            )

        # modify depth_baf_df
        depth_baf_df.rename(
            columns={'depth_raw': 'depthratio_raw'}, inplace=True,
        )

        # assign
        self.data[sampleid]['segment'] = segment_df
        self.data[sampleid]['depth_baf_merge'] = depth_baf_df

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

    ##################
    # final plotting #
    ##################

    def draw_centromeres(self, axlist):
        cytoband_gr = ucscdata.get_cytoband_gr(refver=self.refver, as_gr=True)
        cytoband_gr = cytoband_gr[cytoband_gr.Stain == 'acen']

        for ax in axlist:
            self.draw_bgcolors(
                ax, 
                df=cytoband_gr, 
                plot_kwargs=dict(color='green', alpha=0.5, linewidth=0),
            )

    def make_plotdata(self, sampleid):
        self.data[sampleid]['segment_plotdata'] = self.prepare_plot_data(
            self.data[sampleid]['segment']
        )
        self.data[sampleid]['raw_plotdata'] = self.prepare_plot_data(
            self.data[sampleid]['depth_baf_merge']
        )
        self.data[sampleid]['baf_plotdata'] = self.prepare_plot_data(
            self.data[sampleid]['tumor_baf']
        )

    def make_plotdata_before_cp(self, sampleid):
        self.data[sampleid]['depthratio_raw_plotdata'] = self.prepare_plot_data(
            self.data[sampleid]['depthratio_upscaled']
        )
        self.data[sampleid]['depthratio_segment_plotdata'] = self.prepare_plot_data(
            self.data[sampleid]['depthratio_segment']
        )

        self.data[sampleid]['baf_raw_plotdata'] = self.prepare_plot_data(
            self.data[sampleid]['tumor_baf']
        )
        self.data[sampleid]['baf_segment_plotdata'] = self.prepare_plot_data(
            self.data[sampleid]['baf_segment']
        )

    def make_plotdata_after_cp(self, sampleid):
        self.data[sampleid]['depthratio_raw_plotdata'] = self.prepare_plot_data(
            self.data[sampleid]['depthratio_upscaled']
        )

        self.data[sampleid]['baf_raw_plotdata'] = self.prepare_plot_data(
            self.data[sampleid]['tumor_baf']
        )

        self.data[sampleid]['merged_segment_plotdata'] = self.prepare_plot_data(
            self.data[sampleid]['merged_segment']
        )

    def draw_depthratio_ax(self, sampleid, ax, use_merged_segment=True, draw_predicted=True):
        ax.set_ylabel('depth ratio')
        ax.set_ylim(
            0, 
            np.nanmean(self.data[sampleid]['depthratio_upscaled']['depthratio_raw']) * 2,
        )

        self.draw_dots(
            ax, 
            df_plotdata=self.data[sampleid]['depthratio_raw_plotdata'], 
            y_colname='depthratio_raw', 
            plot_kwargs={'color': 'black', 'markersize': 0.3, 'alpha': 0.01},
        )

        if use_merged_segment:
            segment_df_plotdata = self.data[sampleid]['merged_segment_plotdata']
        else:
            segment_df_plotdata = self.data[sampleid]['depthratio_segment_plotdata']

        self.draw_hlines(
            ax, 
            df_plotdata=segment_df_plotdata,
            y_colname='depthratio_segment_mean', 
            plot_kwargs={'color': 'tab:blue', 'linewidth': 2, 'alpha': 0.5},
        )

        if draw_predicted:
            if 'depthratio_predicted' in segment_df_plotdata.columns:
                self.draw_hlines(
                    ax, 
                    df_plotdata=segment_df_plotdata,
                    y_colname='depthratio_predicted', 
                    plot_kwargs={'color': 'red', 'linewidth': 2, 'alpha': 0.5},
                )

        self.draw_chrom_borders(ax)
        self.draw_genomecoord_labels(ax, n=10)

    def draw_baf_ax(self, sampleid, ax, use_merged_segment=True, draw_predicted=True):
        ax.set_ylabel('baf')
        ax.set_ylim(0, 1)

        self.draw_dots(
            ax, 
            df_plotdata=self.data[sampleid]['baf_raw_plotdata'], 
            y_colname='baf_raw', 
            plot_kwargs={'color': 'black', 'markersize': 0.3, 'alpha': 0.01},
        )

        if use_merged_segment:
            segment_df_plotdata = self.data[sampleid]['merged_segment_plotdata']
        else:
            segment_df_plotdata = self.data[sampleid]['baf_segment_plotdata']

        self.draw_hlines(
            ax, 
            df_plotdata=segment_df_plotdata,
            y_colname='baf_segment_mean', 
            plot_kwargs={'color': 'tab:blue', 'linewidth': 2, 'alpha': 0.8},
        )
        if 'corrected_baf_segment_mean' in segment_df_plotdata.columns:
            self.draw_hlines(
                ax, 
                df_plotdata=segment_df_plotdata,
                y_colname='corrected_baf_segment_mean', 
                plot_kwargs={'color': 'green', 'linewidth': 2, 'alpha': 0.8},
            )

        if draw_predicted:
            if 'baf_predicted' in segment_df_plotdata.columns:
                self.draw_hlines(
                    ax, 
                    df_plotdata=segment_df_plotdata,
                    y_colname='baf_predicted', 
                    plot_kwargs={'color': 'red', 'linewidth': 2, 'alpha': 0.8},
                )

        self.draw_chrom_borders(ax)
        self.draw_genomecoord_labels(ax, n=10)

    def draw_CN_ax(self, sampleid, ax):
        ax.set_ylabel('CN')
        weights = self.data[sampleid]['merged_segment']['End'] - self.data[sampleid]['merged_segment']['Start']
        CNts = self.data[sampleid]['merged_segment']['CNt']
        ax.set_ylim(
            -1, 
            int(common.nanaverage(CNts, weights) * 5),
        )

        self.draw_hlines(
            ax, 
            df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
            y_colname='A', 
            offset=0.2,
            plot_kwargs={'color': 'red'},
        )
        self.draw_hlines(
            ax, 
            df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
            y_colname='B', 
            offset=0,
            plot_kwargs={'color': 'blue'},
        )
        self.draw_chrom_borders(ax)
        self.draw_genomecoord_labels(ax, n=10)

    def plot_final_wgs(self, sampleid):
        return self.plot(sampleid, figsize=(30, 13), draw_invalid_regions=False)

    def plot_woCN_base(self, sampleid, figsize, draw_invalid_regions):
        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
        self.make_plotdata_before_cp(sampleid)
        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')

        fig, axd = plt.subplot_mosaic(
            [
                ['baf',], 
                ['depth',],
            ],
            figsize=figsize,
            gridspec_kw=dict(hspace=0.7),
        )
        for ax in axd.values():
            self.set_xlim(ax)

        fig.suptitle(
            f'sample_id={sampleid}, is_female={self.data[sampleid]["is_female"]}',
            fontsize=20,
        )

        self.draw_baf_ax(sampleid, axd['baf'], draw_predicted=False)
        self.draw_depthratio_ax(sampleid, axd['depth'], draw_predicted=False)
        self.draw_centromeres(axd.values())

        return fig, axd

    def plot_woCN(
        self, sampleid, figsize=(30, 9), draw_invalid_regions=False,
        region_chroms=None, region_start0s=None, region_end0s=None, weights=1,
    ):
        if region_chroms is not None:
            new_region_df = self.make_new_region_df(
                self.refver, region_chroms, region_start0s, region_end0s, weights
            )

            self.old_cconv = self.cconv
            self.cconv = CoordConverter(refver=self.refver, df=new_region_df)

            fig, axd = self.plot_woCN_base(sampleid, figsize=figsize, draw_invalid_regions=draw_invalid_regions)

            self.cconv = self.old_cconv
        else:
            fig, axd = self.plot_woCN_base(sampleid, figsize=figsize, draw_invalid_regions=draw_invalid_regions)

        return fig, axd

    def plot_final_base(self, sampleid, figsize, draw_invalid_regions):
        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
        self.make_plotdata_after_cp(sampleid)
        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')

        fig, axd = plt.subplot_mosaic(
            [
                ['CN',], 
                ['baf',], 
                ['depth',],
            ],
            figsize=figsize,
            gridspec_kw=dict(hspace=0.7),
        )
        for ax in axd.values():
            self.set_xlim(ax)

        fig.suptitle(
            f'sample_id={sampleid}, is_female={self.data[sampleid]["is_female"]}',
            fontsize=20,
        )

        self.draw_baf_ax(sampleid, axd['baf'])
        self.draw_depthratio_ax(sampleid, axd['depth'])
        self.draw_CN_ax(sampleid, axd['CN'])
        self.draw_centromeres(axd.values())

        return fig, axd

    def plot_final(
        self, sampleid, figsize=(30, 12), draw_invalid_regions=False,
        region_chroms=None, region_start0s=None, region_end0s=None, weights=1,
    ):
        if region_chroms is not None:
            new_region_df = self.make_new_region_df(
                self.refver, region_chroms, region_start0s, region_end0s, weights
            )

            self.old_cconv = self.cconv
            self.cconv = CoordConverter(refver=self.refver, df=new_region_df)

            fig, axd = self.plot_final_base(sampleid, figsize=figsize, draw_invalid_regions=draw_invalid_regions)

            self.cconv = self.old_cconv
        else:
            fig, axd = self.plot_final_base(sampleid, figsize=figsize, draw_invalid_regions=draw_invalid_regions)

        return fig, axd

    def plot(self, sampleid, figsize=(30, 13), draw_invalid_regions=False):
        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
        self.make_plotdata(sampleid)
        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')

        fig, axd = plt.subplot_mosaic(
            [
                ['CN',], 
                ['baf',], 
                ['depth',],
            ],
            figsize=figsize,
        )
        for ax in axd.values():
            self.set_xlim(ax)

        # CN
        self.draw_CN_ax(sampleid, axd['CN'])

        # baf
        self.draw_baf_ax(sampleid, axd['baf'])

        # depthratio
        self.draw_depthratio_ax(sampleid, axd['depth'])

        # invalid regions
        if draw_invalid_regions:
            selector = np.logical_or(
                self.data[sampleid]['depth_baf_merge']['depthratio_raw'].isna().to_numpy,
                self.data[sampleid]['depth_baf_merge']['baf_raw'].isna().to_numpy,
            )
            self.draw_bgcolors(
                axd['depth'], 
                df=self.data[sampleid]['depth_baf_merge'].loc[selector, :], 
                plot_kwargs=dict(color='red', alpha=0.01),
            )

        # draw centromeric regions
        self.draw_centromeres(axd.values())

        # return
        return fig, axd

    #########
    # depth #
    #########

    def load_bam_for_depth(self, sampleid, bam_path, sampletype, use_default_gc=True, **kwargs):
        assert sampletype in ('normal', 'tumor')
        mosdepth_df = self._run_mosdepth(bam_path)
        self.load_mosdepth_df(
            sampleid, mosdepth_df, sampletype, use_default_gc=use_default_gc, **kwargs,
        )

    def load_mosdepth_df(self, sampleid, mosdepth_df, sampletype, use_default_gc=True, **kwargs):
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

    def _run_mosdepth(self, bam_path, binsize=None, region_df=None):
        if region_df is None:
            region_df = self.totalregion_df
        if binsize is None:
            binsize = self.default_binsize

        mosdepth_df, _ = libmosdepth.run_mosdepth(
            bam_path, 
            t=8, 
            use_median=False, 
            region_bed_path=None, 
            region_gr=region_df, 
            window_size=binsize, 
            donot_subset_bam=self.check_is_allregion(),
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
        if mode == 'panel':
            preset_cutoffs = 'panel'
        elif mode == 'wgs':
            if sampletype == 'normal':
                preset_cutoffs = 'normal_wgs'
            elif sampletype == 'tumor':
                preset_cutoffs = 'wgs'
        kwargs['preset_cutoffs'] = preset_cutoffs

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

    def load_tumor_germline_vcf(
        self, sampleid, vcf_path, vcf_sampleid, nproc=1, logging_lineno=50000
    ):
        self.data.setdefault(sampleid, dict())
        vaf_df = variantplus.get_vafdf(
            vcf_path, 
            sampleid=vcf_sampleid, 
            nproc=nproc,
        )
        vaf_df.rename(columns={f'vaf_{vcf_sampleid}': 'vaf_raw'}, inplace=True)
        vaf_df = vaf_df.loc[vaf_df['vaf_raw'].notna().to_numpy(), :]
        vaf_df.reset_index(drop=True, inplace=True)
        vaf_df['baf_raw'] = cnvmisc.get_bafs(vaf_df['vaf_raw'])
        vaf_df = vaf_df.loc[:, ['Chromosome', 'Start', 'End', 'vaf_raw', 'baf_raw']]

        self.data[sampleid]['tumor_baf'] = vaf_df

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
                segment_df, refver, is_female,
            )
        elif mode == 'panel':
            assert target_region is not None
            segment_df = rcopynumber.add_CNn_to_targetseq_segment_gr(
                segment_df, target_region, refver, is_female,
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



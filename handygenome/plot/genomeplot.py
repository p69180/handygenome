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


class CoordConverter:
    def __init__(self, df):
        """Args:
            df: pandas.DataFrame object with mandatory columns 'Chromosome', 
                'Start', 'End'. Each row represents a genomic interval, which 
                will be drawn separately in the plot. The rows need not be 
                sorted by chromosome or position. The order of rows of the 
                input DataFrame is respected.
        """
        # sanity check
        assert (
            isinstance(df, pr.PyRanges)
            or (
                isinstance(df, pd.DataFrame)
                and {'Chromosome', 'Start', 'End'}.issubset(df.columns)
            )
        )

        if isinstance(df, pr.PyRanges):
            df = df.df
        self.set_params(df)

    def set_params(self, df):
        """Set attributes:
            totalregion_df
            totalregion_gr
            plot_interval_start0s
            plot_interval_end0s
        """
        # set totalregion_df and totalregion_gr
        if 'weight' in df.columns:
            totalregion_df = df.loc[
                :, ['Chromosome', 'Start', 'End', 'weight']
            ].copy().reset_index(drop=True)
        else:
            totalregion_df = df.loc[
                :, ['Chromosome', 'Start', 'End']
            ].copy().reset_index(drop=True)
            totalregion_df['weight'] = 1

        totalregion_df['raw_region_length'] = (
            totalregion_df['End'] - totalregion_df['Start']
        ).array
        totalregion_df['plot_region_length'] = (
            totalregion_df['raw_region_length'] * totalregion_df['weight']
        ).array

        cumsum = totalregion_df['plot_region_length'].cumsum()
        cumsum_shift = cumsum.shift(1, fill_value=0)
        totalregion_df['region_start_offset'] = cumsum_shift.array
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
                'region_start_offset': np.array(subdf['region_start_offset']),
            }

    def iter_totalregion_df(self):
        chroms = self.totalregion_df['Chromosome']
        grouper = (chroms != chroms.shift(1)).cumsum()
        return iter(self.totalregion_df.groupby(grouper))

    @property
    def xlim(self):
        start0 = self.totalregion_df['plot_interval_start0s'].iloc[0]
        end0 = self.totalregion_df['plot_interval_end0s'].iloc[-1]
        return (start0, end0)

    def genomic_to_plot(self, chrom, pos0_list):
        if chrom not in self.chromwise_params.keys():
            raise Exception(f'Input "chrom" argument is not included in the plotting region.')

        pos0s = np.array(pos0_list)
        pos0s_expand = pos0s[:, np.newaxis]
        params = self.chromwise_params[chrom]
        contains = np.logical_and(
            (pos0s_expand >= params['start0']), (pos0s_expand < params['end0'])
        )
        pos0s_indexes, intv_indexes = np.where(contains)
            # np.ndarray composed of the indexes of the containing intervals
            # identical intervals can appear many times
        within_region_offsets = (
            params['plot_region_length'][intv_indexes]
            * (
                (pos0s[pos0s_indexes] - params['start0'][intv_indexes]) 
                / params['raw_region_length'][intv_indexes]
            )
        )
        return params['region_start_offset'][intv_indexes] + within_region_offsets

    def genomic_to_plot_with_indexes(self, chrom, pos0_list, indexes):
        plot_coords = self.genomic_to_plot(chrom, pos0_list)
        return (indexes, plot_coords)

    def plot_to_genomic(self, x):
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
    def __init__(self, region_df):
        self.cconv = CoordConverter(region_df)

    def set_xlim(self, ax):
        ax.set_xlim(*self.cconv.xlim)

    def draw_chrom_borders(
        self, ax, 
        text_kwargs=dict(), 
        line_kwargs=dict(),
    ):
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
        for chrom, start0, end0 in chrom_borders:
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
        df=None, df_plotdata=None, offset=0, nproc=None, 
        plot_kwargs=dict(),
    ):
        default_plot_kwargs = {}
        default_plot_kwargs.update(plot_kwargs)

        if df_plotdata is None:
            df_plotdata = self.prepare_plot_data(df, nproc=nproc)
            
        ys, xmins, xmaxs = self._merge_adjacent_hlines_data(
            getattr(df_plotdata['isec_gr'], y_colname).to_numpy() + offset,
            df_plotdata['xmins'],
            df_plotdata['xmaxs'],
        )
        ax.hlines(ys, xmins, xmaxs, **default_plot_kwargs)

    def draw_dots(
        self, ax, y_colname, *, 
        df=None, df_plotdata=None, nproc=None, 
        plot_kwargs=dict(),
    ):
        default_plot_kwargs = {
            'color': 'black',
            'marker': 'o',
            'linestyle': '',
        }
        default_plot_kwargs.update(plot_kwargs)

        if df_plotdata is None:
            df_plotdata = self.prepare_plot_data(df, nproc=nproc)

        xs = (df_plotdata['xmins'] + (df_plotdata['xmaxs'] - 1)) / 2
        ys = getattr(df_plotdata['isec_gr'], y_colname).to_numpy()
        ax.plot(xs, ys, **default_plot_kwargs)

    def draw_bgcolors(
        self, ax, *, 
        df=None, df_plotdata=None, nproc=None,
        plot_kwargs=dict(),
    ):
        default_plot_kwargs = {
            'color': 'yellow',
            'alpha': 0.1,
            'zorder': 0,
        }
        default_plot_kwargs.update(plot_kwargs)

        if df_plotdata is None:
            df_plotdata = self.prepare_plot_data(df, nproc=nproc)

        ylims = ax.get_ylim()
        ymin = ylims[0]
        height = ylims[1] - ylims[0]
        xmaxs = df_plotdata['xmaxs']
        xmins = df_plotdata['xmins']
        widths = xmaxs - xmins

        boxes = [
            Rectangle((xm, ymin), width=w, height=height)
            for (xm, w) in zip(xmins, widths)
        ]
        ax.add_collection(PatchCollection(boxes, **default_plot_kwargs))

    def prepare_plot_data(self, df, nproc=None):
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

    @classmethod
    def _merge_adjacent_hlines_data(cls, ys, xmins, xmaxs):
        """Helper of draw_hlines"""
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

    def _isec_trim_data_df(self, df):
        """helper of prepare_plot_data"""
        assert '_index' not in df.columns

        if isinstance(df, pd.DataFrame):
            gr = pr.PyRanges(df)
            nrow = df.shape[0]
        elif isinstance(df, pr.PyRanges):
            gr = df
            nrow = gr.df.shape[0]

        isec_gr = gr.intersect(self.cconv.totalregion_gr)
        isec_gr._index = list(range(nrow))

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


class CNVPlotter(GenomePlotter):
    def __init__(self, region_df):
        super().__init__(region_df)
        self.data = dict()

    def df_arg_sanitycheck(self, arg):
        if arg is not None:
            assert isinstance(arg, pd.DataFrame)

    def add_sample(
        self, 
        sampleid, 
        normal_depth_df, 
        tumor_depth_df=None, 
        normal_baf_df=None,
        tumor_baf_df=None,
    ):
        for arg in (
            normal_depth_df, 
            tumor_depth_df,
            normal_baf_df,
            tumor_baf_df,
        ):
            self.df_arg_sanitycheck(arg)

        self.data[sampleid] = {
            'normal_depth': normal_depth_df,
            'tumor_depth': tumor_depth_df,
            'normal_baf': normal_baf_df,
            'tumor_baf': tumor_baf_df,
        }

    def make_depthratio(self, sampleid):
        self.data[sampleid]['depthratio'] = cnvmisc.make_depth_ratio_df(
            self.data[sampleid]['tumor_depth'], 
            self.data[sampleid]['normal_depth'], 
            as_gr=False,
        )

    


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
    depthratio_gr = cnvmisc.make_depth_ratio_df(tumor_depth_df, normal_depth_df, as_gr=True)
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
    segment_gr = pyranges_helper.join_new(
        segment_gr, 
        depthratio_gr[['depth_ratio_sequenzastyle']], 
        how='left', merge='mean', find_nearest=False, as_gr=True,
    )
    segment_gr = pyranges_helper.join_new(
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

    CNt_weight=5,
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



import functools
import itertools
import multiprocessing
import operator
import pprint
import inspect
import pickle

import numpy as np
import pandas as pd
import pyranges as pr
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn.cluster
import scipy
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

#####################

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


####################


def check_is_allregion(region_gr, refver):
    """checks if allregion_gr is included within region_gr"""
    region_gr = cnvmisc.arg_into_gr(region_gr)
    allregion_gr = common.DEFAULT_CHROMDICTS[refver].to_gr(assembled_only=True, as_gr=True)
    isec_gr = allregion_gr.intersect(region_gr)
    return (isec_gr[[]].sort().df == allregion_gr[[]].sort().df).all(axis=None)


# region argument handler

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


def handle_region_args(
    refver, 
    region_df=None, 
    chroms=None, start0s=None, end0s=None, weights=None,
    region_gaps=None,
):
    if region_gaps is None:
        region_gaps = 0.1

    if region_df is None:
        if chroms is None:
            region_df = common.DEFAULT_CHROMDICTS[refver].to_gr(
                assembled_only=True, as_gr=False,
            )
        else:
            region_df = make_new_region_df(
                refver, chroms, start0s=start0s, end0s=end0s, weights=weights
            )

    return region_df, region_gaps


# decorator for making plot method

PLOTTER_DECORATOR_REGION_ARG_NAMES = (
    'region_chroms',
    'region_start0s',
    'region_end0s',
    'weights',
    'region_gaps',
)


def plotter_decorator(func):
    sig = inspect.signature(func)
    new_sig = sig.replace(
        parameters=itertools.chain(
            sig.parameters.values(),
            (
                inspect.Parameter(
                    x, inspect.Parameter.POSITIONAL_OR_KEYWORD, default=None,
                ) for x in PLOTTER_DECORATOR_REGION_ARG_NAMES
            ),
        )
    )

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        ba = new_sig.bind(*args, **kwargs)
        ba.apply_defaults()
        argdict = ba.arguments

        fig, axd = plotter_decorator_core(
            basemethod=func, 
            **argdict,
        )

        return fig, axd

    return wrapper


def plotter_decorator_core(
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
        self = kwargs['self']

        # make new genomeplotter
        if region_chroms is None:
            new_region_df = self.genomeplotter.region_df.copy()
        else:
            new_region_df = make_new_region_df(
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


class CoordConverter:
    def __init__(self, refver, df, region_gaps):
        """Args:
            df: pandas.DataFrame object with mandatory columns 'Chromosome', 
                'Start', 'End'. Each row represents a genomic interval, which 
                will be drawn separately in the plot. The rows need not be 
                sorted by chromosome or position. The order of rows of the 
                input DataFrame is respected.
            region_gaps: fraction of sum of gap lengths over plotting region lengths
        """
        assert region_gaps is not None

        self.refver = refver
        self.chromdict = common.DEFAULT_CHROMDICTS[refver]
        self.region_gaps = region_gaps

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
        if (region_gaps != 0) and (df.shape[0] > 1):
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

        # set dfs without gap regions
        gap_indexes = np.char.startswith(
            totalregion_df['Chromosome'].to_numpy().astype(str), '-',
        )
        self.totalregion_df_wogap = totalregion_df.loc[~gap_indexes, :]
        self.totalregion_gr_wogap = pr.PyRanges(self.totalregion_df_wogap)

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

    def iter_totalregion_df(self, merge_same_chroms=True):
        totalregion_df_subset = self.totalregion_df.loc[
            ~np.char.startswith(
                self.totalregion_df['Chromosome'].to_numpy().astype(str),
                '-',
            ),
            :
        ]
        if merge_same_chroms:
            chroms = totalregion_df_subset['Chromosome']
            grouper = (chroms != chroms.shift(1)).cumsum()
        else:
            grouper = np.arange(totalregion_df_subset.shape[0])
        return iter(totalregion_df_subset.groupby(grouper))

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
    def get_chrom_borders(self, merge_same_chroms=True):
        """Chromosome names starting with "-1" are omitted"""
        result = list()
        for key, subdf in self.iter_totalregion_df(merge_same_chroms=merge_same_chroms):
            chroms = set(subdf['Chromosome'])
            assert len(chroms) == 1
            chrom = chroms.pop()
            if chrom.startswith('-'):
                continue

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
        region_df, region_gaps = handle_region_args(
            refver, 
            region_df=region_df, 
            chroms=chroms, 
            start0s=start0s, 
            end0s=end0s, 
            weights=weights,
            region_gaps=region_gaps,
        )
            
        self.refver = refver
        self.region_gaps = region_gaps
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
        #pos1_strings = common.shorten_int(pos1s)
        pos1_strings = [f'{x:,}' for x in pos1s]

        labels = [f'{x} : {y}' for x, y in zip(chroms, pos1_strings)]
        ax.set_xticks(plotcoords, labels=labels, minor=False, rotation=90)

    def fit_spines_to_regions(
        self, ax, ylims,
        draw_chromlabel=True,
        prefix_with_chr=True,
        chromlabel_offset=0.01,
        chromlabel_kwargs=dict(), 
        line_kwargs=dict(),
        merge_same_chroms=True,
    ):
        """Should be done after data drawings are finished"""
        # set plotting kwargs
        chromlabel_kwargs = (
            dict(ha='center', va='bottom', size=12)
            | chromlabel_kwargs
        )
        line_kwargs = (
            dict(color='black', linewidth=1)
            | line_kwargs
        )

        # main
        chrom_borders = self.cconv.get_chrom_borders(
            merge_same_chroms=merge_same_chroms,
        )

        # horizontal spines
        ax.spines[['top', 'bottom']].set_visible(False)
        _, start0s, end0s = zip(*chrom_borders)
        ax.hlines(
            np.repeat(ylims[1], len(start0s)), start0s, end0s, 
            color='black', linewidth=1,
        )
        ax.hlines(
            np.repeat(ylims[0], len(start0s)), start0s, end0s, 
            color='black', linewidth=1.5,
        )

        # draw vertical spines
        ax.spines[['left', 'right']].set_visible(False)
        ax.vlines(ax.get_xlim(), ylims[0], ylims[1], **line_kwargs)

        # draw chromosome region borderlines
        border_pos0s = set()
        xlim = self.cconv.xlim
        for _, start0, end0 in chrom_borders:
            if start0 != xlim[0]:
                border_pos0s.add(start0)
            if end0 != xlim[1]:
                border_pos0s.add(end0)

        ax.vlines(tuple(border_pos0s), ylims[0], ylims[1], **line_kwargs)

        # draw chromosome name texts
        if draw_chromlabel:
            chromlabel_y = ylims[1] + chromlabel_offset * (ylims[1] - ylims[0])
            for chrom, start0, end0 in chrom_borders:
                if chrom.startswith('-'):
                    continue
                   
                if prefix_with_chr:
                    if not chrom.startswith('chr'):
                        chrom = 'chr' + chrom

                ax.text(
                    (start0 + end0) / 2, 
                    chromlabel_y, 
                    chrom, 
                    **chromlabel_kwargs,
                )

    def draw_hlines(
        self, 
        ax, 
        *, 
        y_colname=None, 
        y_val=None,
        df=None, 
        df_plotdata=None, 
        offset=0,
        plot_kwargs=dict(),
    ):
        plot_kwargs = (
            dict()
            | plot_kwargs
        )

        if df_plotdata is None:
            df_plotdata = self.prepare_plot_data(df)

        if df_plotdata is False:
            return

        # draw data
        if y_colname is None:
            ys = np.repeat(y_val, df_plotdata.shape[0]) + offset
        else:
            ys = df_plotdata[y_colname].to_numpy() + offset

        if df_plotdata.shape[0] > 1:
            ys, xmins, xmaxs = self._merge_adjacent_data_new(
                genome_xmins=df_plotdata['Start'], 
                genome_xmaxs=df_plotdata['End'], 
                plot_xmins=df_plotdata['plot_start0s'], 
                plot_xmaxs=df_plotdata['plot_end0s'], 
                ys=ys,
            )
        else:
            row = df_plotdata.iloc[0, :]
            xmins = row['plot_start0s']
            xmaxs = row['plot_end0s']

        ax.hlines(ys, xmins, xmaxs, **plot_kwargs)

    def draw_dots(
        self, ax, 
        *, 
        y_colname=None, 
        y_val=None,
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
        if y_colname is None:
            ys = np.repeat(y_val, len(xs))
        else:
            ys = df_plotdata[y_colname].to_numpy()

        ax.plot(xs, ys, **plot_kwargs)

    def draw_dots_scatter(
        self, ax, y_colname, *, 
        df=None, df_plotdata=None,
        color_colname=None,
        color_vals=None,
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
            elif color_vals is not None:
                plot_kwargs['c'] = color_vals

            if 'color' in plot_kwargs:
                del plot_kwargs['color']
            if 'markersize' in plot_kwargs:
                plot_kwargs['s'] = plot_kwargs['markersize']
                del plot_kwargs['markersize']

            ax.scatter(xs, ys, **plot_kwargs)

    def draw_bars(
        self, ax, y_colname, *, 
        df=None, df_plotdata=None,
        bottom_colname=None,
        plot_kwargs=dict(),
    ):
        plot_kwargs = (
            dict(alpha=0.5, color='tab:blue')
            | plot_kwargs
        )

        if df_plotdata is None:
            df_plotdata = self.prepare_plot_data(df)

        if df_plotdata is False:
            return

        xs = (df_plotdata['plot_end0s'] + df_plotdata['plot_start0s']) / 2
        widths = df_plotdata['plot_end0s'] - df_plotdata['plot_start0s']
        heights = df_plotdata[y_colname]
        if bottom_colname is None:
            bottoms = 0 
        else:
            bottoms = df_plotdata[bottom_colname]

        ax.bar(xs, height=heights, width=widths, bottom=bottoms, **plot_kwargs)

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
        #xmaxs = df_plotdata['plot_end0s']
        #xmins = df_plotdata['plot_start0s']
        _, xmins, xmaxs = self._merge_adjacent_data_new(
            genome_xmins=df_plotdata['Start'], 
            genome_xmaxs=df_plotdata['End'], 
            plot_xmins=df_plotdata['plot_start0s'], 
            plot_xmaxs=df_plotdata['plot_end0s'], 
            ys=None,
        )
        #_, xmins, xmaxs = self._merge_adjacent_data(
        #    xmins=xmins, xmaxs=xmaxs, ys=None,
        #)
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

    def draw_centromeres(self, ax, ymins=None, ymaxs=None):
        cytoband_gr = ucscdata.get_cytoband_gr(refver=self.refver, as_gr=True)
        self.draw_bgcolors(
            ax, 
            df=cytoband_gr[cytoband_gr.Stain == 'acen'], 
            #df=cytoband_gr[cytoband_gr.Stain.isin(['acen', 'gvar', 'stalk'])], 
            ymins=ymins,
            ymaxs=ymaxs,
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
            self.draw_bgcolors(
                ax, 
                df=cytoband_gr[cytoband_gr.Stain == bandname], 
                colors=mapping[bandname],
                plot_kwargs=dict(alpha=0.3, linewidth=0),
            )

        helper('acen')
        helper('gvar')
        helper('stalk')

    def draw_grids(self, ax, ys, line_params=dict(), merge_same_chroms=True):
        line_params = (
            dict(color='black', linewidth=0.2, alpha=0.5)
            | line_params
        )
        chroms, start0s, end0s = zip(
            *self.cconv.get_chrom_borders(
                merge_same_chroms=merge_same_chroms,
            )
        )
        for y in ys:
            ax.hlines(
                np.repeat(y, len(start0s)), start0s, end0s, 
                **line_params,
            )

    def draw_ax_common(
        self, 
        ax, 
        n_xlabel=None,
        split_spines=True,
        merge_same_chroms=True,
        chromlabel_kwargs=dict(), 
        draw_chromlabel=True,
    ):
        self.set_xlim(ax)
        ylims = ax.get_ylim()

        if split_spines:
            self.fit_spines_to_regions(
                ax,
                ylims=ylims,
                draw_chromlabel=draw_chromlabel,
                prefix_with_chr=True,
                merge_same_chroms=merge_same_chroms,
                chromlabel_kwargs=chromlabel_kwargs,
            )

        yticks = [
            x for x in ax.get_yticks()
            if (x > ylims[0]) and (x < ylims[1])
        ]
        self.draw_grids(
            ax, 
            ys=yticks, 
            line_params=dict(), 
            merge_same_chroms=merge_same_chroms,
        )

        self.draw_centromeres(ax, ymins=ylims[0], ymaxs=ylims[1])
        #self.draw_centromeres_type2(ax)
        if n_xlabel is not None:
            self.draw_genomecoord_labels(ax, n=n_xlabel)
        else:
            ax.set_xticks([])

    def draw_features(
        self, 
        ax,
        df=None,
        df_plotdata=None,

        y_features=None,
        y_labels=None,
        #feature_as_dot=False,
        draw_label=True,

        text_kwargs=dict(),
        line_kwargs=dict(),
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

        # make plotdata
        if df_plotdata is None:
            df_plotdata = self.prepare_plot_data(df)
        if df_plotdata is False:
            return


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
            y_val=y_features,
            df_plotdata=df_plotdata, 
            offset=0,
            plot_kwargs=line_kwargs,
        )

        # draw texts
        if draw_label:
            assert 'name' in df_plotdata.columns
            for name, subdf in df_plotdata.groupby('name'):
                subdf = cnvmisc.sort_genome_df(subdf, refver=self.refver)
                if subdf.shape[0] == 1:
                    xmins = subdf['plot_start0s']
                    xmaxs = subdf['plot_end0s']
                else:
                    _, xmins, xmaxs = self._merge_adjacent_data_new(
                        genome_xmins=subdf['Start'], 
                        genome_xmaxs=subdf['End'], 
                        plot_xmins=subdf['plot_start0s'], 
                        plot_xmaxs=subdf['plot_end0s'], 
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

    def prepare_plot_data(self, df):
        gr = cnvmisc.arg_into_gr(df)
        isec_gr = gr.intersect(self.cconv.totalregion_gr_wogap).sort()
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

    @classmethod
    def _merge_adjacent_data(cls, xmins, xmaxs, ys=None):
        """Helper of draw_hlines"""
        if ys is None:
            ys = np.repeat(1, len(xmins))

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
    def __init__(
        self, 
        refver, 
        region_df=None, 
        chroms=None, start0s=None, end0s=None, weights=None,
        region_gaps=None,
    ):
        #super().__init__(refver=refver, region_df=region_df)
        self.refver = refver

        region_df, region_gaps = handle_region_args(
            refver, 
            region_df=region_df, 
            chroms=chroms, 
            start0s=start0s, 
            end0s=end0s, 
            weights=weights,
            region_gaps=region_gaps,
        )

        self.genomeplotter = GenomePlotter(
            refver, region_df=region_df, region_gaps=region_gaps,
        )
        self.data = dict()
        self.default_binsize = 100

    def reset_genomeplotter(
        self, 
        region_df=None, 
        chroms=None, start0s=None, end0s=None, weights=None,
        region_gaps=None,
    ):
        region_df, region_gaps = handle_region_args(
            refver=self.refver, 
            region_df=region_df, 
            chroms=chroms, 
            start0s=start0s, 
            end0s=end0s, 
            weights=weights,
            region_gaps=region_gaps,
        )
        self.genomeplotter = GenomePlotter(
            self.refver, region_df=region_df, region_gaps=region_gaps,
        )

    def save_data(self, sampleid, outfile_path):
        with open(outfile_path, 'wb') as outfile:
            pickle.dump(self.data[sampleid], outfile)

    def load_data(self, sampleid, infile_path):
        with open(infile_path, 'rb') as infile:
            self.data[sampleid] = pickle.load(infile)

    ##############
    # mainstream #
    ##############

    @deco.get_deco_num_set_differently(
        ('normal_bam_path', 'normal_depth_path', 'normal_depth_df'), 2, 'lt',
    )
    @deco.get_deco_num_set_differently(
        ('tumor_bam_path', 'tumor_depth_path', 'tumor_depth_df'), 2, 'lt',
    )
    @deco.get_deco_num_set_differently(
        ('germline_vcf_path', 'germline_vafdf_path'), 2, 'lt',
    )
    def add_sample_file_new(
        self, 
        sampleid, 
        is_female,

        germline_vcf_path=None,
        germline_vafdf_path=None,
        vcf_sampleid_tumor=None,
        vcf_sampleid_normal=None,

        *,

        mode='wgs',
        target_region=None,

        normal_bam_path=None,
        normal_depth_path=None, 
        normal_depth_df=None, 

        tumor_bam_path=None,
        tumor_depth_path=None, 
        tumor_depth_df=None, 

        vcfload_nproc=1,
        #mosdepth_postprocess_kwargs=dict(),

        verbose=True,
    ):
        """Args:
            *_depth_path: mosdepth output file
            germline_vcf_path: germline variant vcf
            vcf_sampleids: tuple of (normal sample id, tumor sample id)
        """
        def depth_loading_helper(
            self, 
            bam_path, 
            depth_path, 
            depth_df, 
            sampleid, 
            sampletype,
        ):
            assert sampletype in ('tumor', 'normal')

            if bam_path is not None:
                LOGGER_INFO.info(f'Loading {sampletype} depth - running mosdepth')
                depth_df = libmosdepth.run_mosdepth(
                    bam_path, 
                    t=8, 
                    use_median=False, 
                    region_bed_path=None, 
                    region_gr=(
                        None
                        if self.data[sampleid]['mode'] == 'wgs' else
                        self.data[sampleid]['target_region']
                    ), 
                    window_size=self.default_binsize, 
                    donot_subset_bam=True,
                    as_gr=False, 
                    load_perbase=False,
                )
            elif depth_path is not None:
                LOGGER_INFO.info(f'Loading {sampletype} depth - reading mosdepth output file')
                depth_df = libmosdepth.load_mosdepth_output(
                    depth_path, depth_colname='mean_depth', as_gr=False,
                )
            elif depth_df is not None:
                LOGGER_INFO.info(f'Loading {sampletype} depth - using the given depth dataframe')
                assert isinstance(depth_df, pd.DataFrame)
                assert set(depth_df.columns) == {'Chromosome', 'Start', 'End', 'mean_depth'}

            if depth_df is not None:
                self.data[sampleid][f'{sampletype}_depth'] = depth_df

        def sanity_check(
            mode, target_region, vcf_sampleid_tumor, vcf_sampleid_normal, germline_vcf_path,
        ):
            # target_region
            assert mode in ('wgs', 'panel')
            if (mode == 'panel') and (target_region is None):
                raise Exception(f'"target_region" must be given when "mode" is "panel"')
            elif (mode == 'wgs') and (target_region is not None):
                raise Exception(f'"target_region" must not be given when "mode" is "wgs"')

            # germline VCF file arguments
            if germline_vcf_path is not None:
                if vcf_sampleid_tumor is None:
                    raise Exception(f'When "germline_vcf_path" is used, "vcf_sampleid_tumor" must be given.')
                if (mode == 'wgs') and (vcf_sampleid_normal is None):
                    raise Exception(f'"vcf_sampleid_normal" must be given when "mode" is "wgs"')

        # main
        sanity_check(mode, target_region, vcf_sampleid_tumor, vcf_sampleid_normal, germline_vcf_path)
        self.data[sampleid] = dict()
        self.data[sampleid]['is_female'] = is_female
        self.data[sampleid]['mode'] = mode

        # load normal
        depth_loading_helper(
            self, 
            normal_bam_path, 
            normal_depth_path, 
            normal_depth_df,
            sampleid, 
            'normal', 
        )

        # load tumor
        depth_loading_helper(
            self, 
            tumor_bam_path, 
            tumor_depth_path, 
            tumor_depth_df,
            sampleid, 
            'tumor', 
        )

        # set target_region
        self.set_target_region(sampleid, mode, target_region)

        # normal mean ploidy
        self.set_normal_mean_ploidy(sampleid)

        # load germline vcf
        if germline_vcf_path is not None:
            LOGGER_INFO.info('Loading tumor germline vcf')
            self.load_germline_vcf(
                sampleid=sampleid, 
                vcf_path=germline_vcf_path, 
                vcf_sampleid_tumor=vcf_sampleid_tumor,
                vcf_sampleid_normal=vcf_sampleid_normal,
                logging_lineno=50000,
                nproc=vcfload_nproc,
            )
        elif germline_vafdf_path is not None:
            LOGGER_INFO.info('Loading tumor germline vaf dataframe')
            self.load_germline_vafdf(sampleid, germline_vafdf_path)

        #self.postprocess_bafdf(sampleid)

        # postprocess depths
        self.postprocess_depth(sampleid, verbose=verbose)

    def load_germline_vafdf(self, sampleid, vafdf_path):
        self.data[sampleid]['original_baf'] = pd.read_csv(
            vafdf_path,
            sep='\t',
            dtype={
                'Chromosome': 'string',  
                'Start': int,   
                'End': int,     
                'vaf_raw_tumor': float,   
                'baf_raw_tumor': float,   
                'vaf_raw_normal': float,  
                'baf_raw_normal': float,
            },
        )

    def postprocess_bafdf(self, sampleid):
        modified_baf = self.data[sampleid]['original_baf'].copy()
        #self.data[sampleid]['baf'] = modified_baf.loc[modified_baf['baf_raw_tumor'] > 0, :]

    def add_bafpeak_to_segment(self, sampleid, bw=1):
        # join segment and raw bafs
        left = self.data[sampleid]['baf_segment']

        right = self.data[sampleid]['original_baf']
        right = right.loc[
            right['baf_raw_tumor'] > 0, 
            ['Chromosome', 'Start', 'End', 'baf_raw_tumor'],
        ]

        joined = pyranges_helper.join(
            left, right, how='left', merge=None, sort=True, refver=self.refver,
        )
        # find peaks
        groupkey = cnvmisc.genome_df_groupkey(joined, refver=self.refver)
        peaks = list()
        for k, v in joined.groupby(groupkey)[['baf_raw_tumor']]:
            baf_values = v['baf_raw_tumor'].to_numpy()
            peaks.append(
                libbaf.infer_baf_density(baf_values, bw=bw, rmzero=False)
            )
        # assign values
        assert len(peaks) == left.shape[0], (
            f'The number of groupby groups and segment df row number are different'
        )
        left['baf_segment_peak'] = peaks
        self.data[sampleid]['baf_segment'] = left

    def postprocess_depth(self, sampleid, verbose=False):
        logger = (LOGGER_DEBUG if verbose else LOGGER_INFO)

        # get GC df
        if self.data[sampleid]['mode'] == 'wgs':
            logger.info(f'Getting gc fraction dataframe')
            gc_df = libgcfraction.get_gc_df(
                self.refver, 
                self.default_binsize, 
                coords_as_index=True,
            )
        elif self.data[sampleid]['mode'] == 'panel':
            gc_df = None

        # main
        for key in ('normal', 'tumor'):
            if f'{key}_depth' not in self.data[sampleid].keys():
                continue

            logger.info(f'Beginning postprocess of {key} depth')
            refver_arg = (
                self.refver 
                if self.data[sampleid]['mode'] == 'panel' else 
                None
            )
            output_depth_df, gcbin_average_depths = cnvmisc.postprocess_depth_df(
                self.data[sampleid][f'{key}_depth'],
                refver=refver_arg,
                gc_df=gc_df,
                included_region=self.data[sampleid]['target_region'],
                as_gr=False,
                verbose=verbose,
                add_norm_depth=True,
            )

            self.data[sampleid][f'{key}_depth'] = output_depth_df
            self.data[sampleid][f'{key}_gcdata'] = gcbin_average_depths

    def set_depthratio(self, sampleid):
        depthratio_df = cnvmisc.make_depth_ratio(
            self.data[sampleid]['tumor_depth'], 
            self.data[sampleid]['normal_depth'],
            make_depthratio_mystyle=False,
            make_depthratio_plain=True,
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

    def upscale_depth(self, sampleid, binsize=1000, do_normal=True, do_tumor=True):
        if do_normal:
            self.data[sampleid]['normal_depth_upscaled'] = cnvmisc.upsize_depth_df_bin(
                self.data[sampleid]['normal_depth'], 
                size=binsize, 
                refver=self.refver,
            )

        if do_tumor:
            self.data[sampleid]['tumor_depth_upscaled'] = cnvmisc.upsize_depth_df_bin(
                self.data[sampleid]['tumor_depth'], 
                size=binsize, 
                refver=self.refver,
            )

    def make_segments(
        self,
        sampleid,
        winsorize=False,
        depthratio_gamma=None,
        depthratio_kmin=None,
        baf_gamma=100,
        baf_kmin=None,
        verbose=False,
        segment_baf_cutoff=0.1,

        bafcorrection_cutoff=None,
        bafcorrector=None,

        bw=1,
    ):
        depthratio_df = (
            self.data[sampleid]['depthratio_upscaled']
            if 'depthratio_upscaled' in self.data[sampleid] else
            self.data[sampleid]['depthratio']
        )
        depth_segment, baf_segment = _make_segments_main(
            depthratio_df=depthratio_df,
            mode=self.data[sampleid]['mode'],
            refver=self.refver,
            winsorize=winsorize,
            depthratio_gamma=depthratio_gamma,
            depthratio_kmin=depthratio_kmin,
            baf_gamma=baf_gamma,
            baf_kmin=baf_kmin,
            verbose=verbose,

            baf_df=self.data[sampleid]['original_baf'],
            target_region=self.data[sampleid]['target_region'],
            baf_cutoff=segment_baf_cutoff,
        )
        self.data[sampleid]['depthratio_segment'] = depth_segment
        self.data[sampleid]['baf_segment'] = baf_segment
        #if self.data[sampleid]['mode'] == 'wgs':
        self.add_bafpeak_to_segment(sampleid, bw=bw)
        self.make_merged_segment(sampleid)
        self.add_corrected_baf(
            sampleid, 
            round_cutoff=bafcorrection_cutoff,
            bafcorrector=bafcorrector,
        )

    def add_corrected_baf(
        self, 
        sampleid, 
        round_cutoff=None,
        bafcorrector=None,
    ):
        if bafcorrector is None:
            bafcorrector = libbaf.load_bafcorrect_func(x_cutoff=round_cutoff)

        for df in (
            self.data[sampleid]['baf_segment'], 
            self.data[sampleid]['merged_segment'], 
        ):
            df['corrected_baf_segment_mean'] = bafcorrector(df['baf_segment_peak'].to_numpy())

    def get_cp_from_twodata(
        self, sampleid, depthratio1, CNt1, depthratio2, CNt2,
    ):
        return cnvmisc.get_cp_from_twodata(
            depthratio1, 
            CNt1, 
            depthratio2, 
            CNt2, 
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'], 
            CNn=2, 
            tumor_avg_depth_ratio=1, 
            normal_avg_depth_ratio=1,
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

    @plotter_decorator
    def plot_beforecp(
        self, 
        sampleid, 
        figsize=None, 
        hspace=None,
        draw_invalid_regions=False, 
        use_saved_plotdata=False,
        use_merged_segment=True,

        draw_depthratio_hist=True,
        rm_haploid_from_hist=True,
        depthratio_hist_threshold=3,

        draw_depth=False,
        is_rawdepth=True,
        depth_binsize=10000,
        depth_ymax=None,

        n_xlabel=None,
        depthratio_ymax=None,
        
        # plot properties

        depthratio_dot_kwargs=dict(),
        depthratio_line_segmean_kwargs=dict(),
        #depthratio_line_predict_kwargs=dict(),

        depthratio_hist_annotate_kwargs=dict(),
        depthratio_hist_plot_kwargs=dict(),

        depth_dot_kwargs=dict(),

        baf_dot_kwargs=dict(),
        baf_line_segmean_kwargs=dict(),
        baf_line_corr_segmean_kwargs=dict(),
        #baf_line_predict_kwargs=dict(),
    ):
        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
        self.make_plotdata_basic(sampleid)
        if draw_depth:
            self.make_plotdata_fordepth(sampleid, binsize=depth_binsize)
        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')

        fig, axd = self.make_axd(
            figsize=figsize, 
            hspace=hspace, 
            draw_depthratio_hist=draw_depthratio_hist, 
            draw_solution=False,
            draw_tumor_baf=True,
            draw_depthratio=True,
            draw_normal_baf=draw_depth,
            draw_tumor_depth=draw_depth,
            draw_normal_depth=draw_depth,
        )

        fig.suptitle(
            f'sample_id={sampleid}, is_female={self.data[sampleid]["is_female"]}',
            fontsize=20,
        )

        self.draw_baf_ax(
            sampleid, 
            axd['baf'], 
            use_merged_segment=True,

            draw_predicted=False,
            draw_corrected=True,
            draw_segmean=True,
            n_xlabel=n_xlabel,

            is_tumor=True,

            mark_unfit_regions=False,

            dot_kwargs=baf_dot_kwargs,
            line_segmean_kwargs=baf_line_segmean_kwargs,
            line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,

            make_raw_plotdata=False,
            make_segment_plotdata=False,
        )
        self.draw_depthratio_ax(
            sampleid, 
            axd['depthratio'], 
            use_merged_segment=True,

            draw_predicted=False,
            draw_segmean=True,
            draw_deviation=False,
            draw_depthratio_peaks=False,

            n_xlabel=n_xlabel,
            dot_kwargs=depthratio_dot_kwargs,
            line_segmean_kwargs=depthratio_line_segmean_kwargs,
            ymax=depthratio_ymax,

            make_raw_plotdata=False,
            make_segment_plotdata=False,
        )

        if draw_depthratio_hist:
            peak_depthratios = self.draw_depthratio_hist_ax(
                sampleid,
                axd['depthratio_hist'],
                use_merged_segment=True, 
                depth_ylim=axd['depthratio'].get_ylim(),
                rm_haploid=rm_haploid_from_hist,
                peak_threshold=depthratio_hist_threshold,
                annotate_kwargs=depthratio_hist_annotate_kwargs,
                plot_kwargs=depthratio_hist_plot_kwargs,
            )

            for y in peak_depthratios:
                axd['depthratio'].axhline(y, color='orange', linewidth=1, alpha=0.6)

        if draw_depth:
            self.draw_depth_bundle(
                sampleid, axd, n_xlabel, is_rawdepth, 
                depth_dot_kwargs, depth_ymax, 
                baf_dot_kwargs, baf_line_segmean_kwargs, baf_line_corr_segmean_kwargs,
            )

        return fig, axd

    @plotter_decorator
    def plot_aftercp_freeccf(
        self, 
        sampleid, 
        cellularity,
        ploidy,

        figsize=(30, 16), 
        hspace=None,

        n_xlabel=None,
        depthratio_ymax=None,
        CN_ymax=None,
        subCN_ymax=None,

        depthratio_std_factor=1,

        draw_depth=False,
        is_rawdepth=True,
        depth_binsize=10000,
        depth_ymax=None,

        depthratio_dot_kwargs=dict(),
        depthratio_line_segmean_kwargs=dict(),
        depthratio_line_predict_kwargs=dict(),
        depthratio_line_predict_clonal_kwargs=dict(),

        baf_dot_kwargs=dict(),
        baf_line_segmean_kwargs=dict(),
        baf_line_corr_segmean_kwargs=dict(),
        baf_line_predict_kwargs=dict(),
        baf_line_predict_clonal_kwargs=dict(),

        CN_line_CNt_kwargs=dict(),
        subCN_line_CNt_kwargs=dict(),
        ccf_bar_kwargs=dict(),

        depth_ratio_diff=None,
        baf_diff=0.05,

        Bn=1,

        min_N_CNt_candidates=5,
        N_CNt_candidates_fraction=0.5,

        limited_clonal=True,
    ):
        LOGGER_INFO.info(f'Beginning calculation of subclonal solution')
        self.make_CN_solution_freeccf(
            sampleid,
            cellularity,
            ploidy,
            depth_ratio_diff=depth_ratio_diff,
            baf_diff=baf_diff,
            min_N_CNt_candidates=min_N_CNt_candidates,
            N_CNt_candidates_fraction=N_CNt_candidates_fraction,
            limited_clonal=limited_clonal,
        )
        self.add_freeccf_solution_to_segment(sampleid)

        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
        self.make_plotdata_basic(sampleid)
        if draw_depth:
            self.make_plotdata_fordepth(sampleid, binsize=depth_binsize)
        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')

        fig, axd = self.make_axd(
            figsize=figsize, 
            hspace=hspace, 
            draw_depthratio_hist=False, 
            draw_solution=True,
            draw_tumor_baf=True,
            draw_depthratio=True,
            draw_normal_baf=draw_depth,
            draw_tumor_depth=draw_depth,
            draw_normal_depth=draw_depth,
        )

        fig.suptitle(
            ', '.join([
                f'sample_id={sampleid}',
                f'is_female={self.data[sampleid]["is_female"]}',
                f'cellularity={round(cellularity, 3)}',
                f'ploidy={round(ploidy, 3)}',
            ]),
            fontsize=20,
        )

        # depth
        self.draw_depthratio_ax(
            sampleid, 
            axd['depthratio'], 
            use_merged_segment=True,

            draw_predicted=True,
            draw_segmean=True,
            draw_deviation=False,
            draw_depthratio_peaks=False,

            mark_unfit_regions=True,

            cellularity=cellularity,
            ploidy=ploidy,
            draw_integerCN_lines=True,

            std_factor=depthratio_std_factor,
            n_xlabel=n_xlabel,
            dot_kwargs=depthratio_dot_kwargs,
            line_segmean_kwargs=depthratio_line_segmean_kwargs,
            line_predict_kwargs=depthratio_line_predict_kwargs,
            line_predict_clonal_kwargs=depthratio_line_predict_clonal_kwargs,
            ymax=depthratio_ymax,

            make_raw_plotdata=False,
            make_segment_plotdata=False,
        )
        if draw_depth:
            self.draw_depth_bundle(
                sampleid, axd, n_xlabel, is_rawdepth, 
                depth_dot_kwargs, depth_ymax, 
                baf_dot_kwargs, baf_line_segmean_kwargs, baf_line_corr_segmean_kwargs,
            )

        # baf
        self.draw_baf_ax(
            sampleid, 
            axd['baf'], 
            use_merged_segment=True,
            draw_predicted=True,
            draw_corrected=True,
            draw_segmean=True,
            n_xlabel=n_xlabel,

            is_tumor=True,
            mark_unfit_regions=True,

            dot_kwargs=baf_dot_kwargs,
            line_segmean_kwargs=baf_line_segmean_kwargs,
            line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,
            line_predict_kwargs=baf_line_predict_kwargs,
            line_predict_clonal_kwargs=baf_line_predict_clonal_kwargs,

            make_raw_plotdata=False,
            make_segment_plotdata=False,
        )

        # clonal CN
        self.draw_CN_ax(
            sampleid, 
            axd['clonal_CN'],
            n_xlabel=n_xlabel,
            line_CNt_kwargs=CN_line_CNt_kwargs,
            ymax=CN_ymax,
            draw_CNt=True,
            draw_A=False,
            draw_B=True,
        )
        # subclonal CN
        self.draw_subclonal_CN_ax(
            sampleid, 
            axd['subclonal_CN'],
            n_xlabel=n_xlabel,
            line_CNt_kwargs=subCN_line_CNt_kwargs,
            ymax=subCN_ymax,
            draw_CNt=True,
            draw_A=False,
            draw_B=True,
        )
        # ccf
        self.draw_ccf_ax(
            sampleid,
            axd['ccf'],
            n_xlabel=n_xlabel,
            bar_kwargs=ccf_bar_kwargs,
        )

        return fig, axd

    def show_ccfs(self, sampleid, bandwidth=0.1):
        self.select_fixed_ccfs(sampleid, bandwidth=bandwidth)
        ccf_plotdata = self.data[sampleid]['ccf_plotdata']
        self.show_ccfs_main(
            ccfs=ccf_plotdata['ccfs'],
            lengths=ccf_plotdata['lengths'],
            density=ccf_plotdata['density'],
            peak_values=ccf_plotdata['peak_values'],
        )

    @plotter_decorator
    def plot_aftercp_fixedccf(
        self, 
        sampleid, 
        cellularity,
        ploidy,

        figsize=(30, 16), 
        hspace=None,

        n_xlabel=None,
        depthratio_ymax=None,
        CN_ymax=None,
        subCN_ymax=None,

        depthratio_std_factor=1,

        draw_depth=False,
        is_rawdepth=True,
        depth_binsize=10000,
        depth_ymax=None,

        depthratio_dot_kwargs=dict(),
        depthratio_line_segmean_kwargs=dict(),
        depthratio_line_predict_kwargs=dict(),
        depthratio_line_predict_clonal_kwargs=dict(),

        baf_dot_kwargs=dict(),
        baf_line_segmean_kwargs=dict(),
        baf_line_corr_segmean_kwargs=dict(),
        baf_line_predict_kwargs=dict(),
        baf_line_predict_clonal_kwargs=dict(),

        CN_line_CNt_kwargs=dict(),
        subCN_line_CNt_kwargs=dict(),
        ccf_bar_kwargs=dict(),

        depth_ratio_diff=None,
        baf_diff=0.05,

        min_N_CNt_candidates=5,
        N_CNt_candidates_fraction=0.5,

        ccf_bw=0.1,

        update_plotdata=False,
        CNt_diff_factor=0.1,

        mark_unfit_regions=False,

        limited_clonal=True,
    ):
        LOGGER_INFO.info(f'Beginning calculation of subclonal solution')
        self.make_CN_solution_after_ccfs(
            sampleid,
            cellularity,
            ploidy,
            min_N_CNt_candidates=min_N_CNt_candidates,
            N_CNt_candidates_fraction=N_CNt_candidates_fraction,
            CNt_diff_factor=CNt_diff_factor,
            limited_clonal=limited_clonal,
        )
        self.add_fixedccf_solution_to_segment(sampleid)
        LOGGER_INFO.info(f'Finished calculation of subclonal solution')

        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
        if update_plotdata:
            self.add_solution_to_plotdata(sampleid)
        else:
            self.make_plotdata_basic(sampleid)
            if draw_depth:
                self.make_plotdata_fordepth(sampleid, binsize=depth_binsize)
        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')

        fig, axd = self.make_axd(
            figsize=figsize, 
            hspace=hspace, 
            draw_depthratio_hist=False, 
            draw_solution=True,
            draw_tumor_baf=True,
            draw_depthratio=True,
            draw_normal_baf=draw_depth,
            draw_tumor_depth=draw_depth,
            draw_normal_depth=draw_depth,
        )

        fig.suptitle(
            ', '.join([
                f'sample_id={sampleid}',
                f'is_female={self.data[sampleid]["is_female"]}',
                f'cellularity={round(cellularity, 3)}',
                f'ploidy={round(ploidy, 3)}',
            ]),
            fontsize=20,
        )

        # depth
        self.draw_depthratio_ax(
            sampleid, 
            axd['depthratio'], 
            use_merged_segment=True,

            draw_predicted=True,
            draw_segmean=True,
            draw_deviation=False,
            draw_depthratio_peaks=False,

            mark_unfit_regions=mark_unfit_regions,

            cellularity=cellularity,
            ploidy=ploidy,
            draw_integerCN_lines=True,

            std_factor=depthratio_std_factor,
            n_xlabel=n_xlabel,
            dot_kwargs=depthratio_dot_kwargs,
            line_segmean_kwargs=depthratio_line_segmean_kwargs,
            line_predict_kwargs=depthratio_line_predict_kwargs,
            line_predict_clonal_kwargs=depthratio_line_predict_clonal_kwargs,
            ymax=depthratio_ymax,

            make_raw_plotdata=False,
            make_segment_plotdata=False,
        )
        if draw_depth:
            self.draw_depth_bundle(
                sampleid, axd, n_xlabel, is_rawdepth, 
                depth_dot_kwargs, depth_ymax, 
                baf_dot_kwargs, baf_line_segmean_kwargs, baf_line_corr_segmean_kwargs,
            )

        # baf
        self.draw_baf_ax(
            sampleid, 
            axd['baf'], 
            use_merged_segment=True,
            draw_predicted=True,
            draw_corrected=True,
            draw_segmean=True,
            n_xlabel=n_xlabel,

            is_tumor=True,
            mark_unfit_regions=mark_unfit_regions,

            dot_kwargs=baf_dot_kwargs,
            line_segmean_kwargs=baf_line_segmean_kwargs,
            line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,
            line_predict_kwargs=baf_line_predict_kwargs,
            line_predict_clonal_kwargs=baf_line_predict_clonal_kwargs,

            make_raw_plotdata=False,
            make_segment_plotdata=False,
        )

        # clonal CN
        self.draw_CN_ax(
            sampleid, 
            axd['clonal_CN'],
            n_xlabel=n_xlabel,
            line_CNt_kwargs=CN_line_CNt_kwargs,
            ymax=CN_ymax,
            draw_CNt=True,
            draw_A=False,
            draw_B=True,
        )
        # subclonal CN
        self.draw_subclonal_CN_ax(
            sampleid, 
            axd['subclonal_CN'],
            n_xlabel=n_xlabel,
            line_CNt_kwargs=subCN_line_CNt_kwargs,
            ymax=subCN_ymax,
            draw_CNt=True,
            draw_A=False,
            draw_B=True,
        )
        # ccf
        self.draw_ccf_ax(
            sampleid,
            axd['ccf'],
            n_xlabel=n_xlabel,
            bar_kwargs=ccf_bar_kwargs,
        )
        fixed_ccfs = self.data[sampleid]['fixed_ccfs']
        for y in fixed_ccfs:
            axd['ccf'].axhline(y, color='red', linewidth=1)

        return fig, axd

    @plotter_decorator
    def plot_custom(
        self, 
        sampleid, 

        fig=None,
        axd=None,
        figsize=(30, 16), 
        hspace=None,
        height_ratios=None,
        title=None,
        title_size=20,
        title_y=0.95,

        draw_depthratio_hist=False, 
        draw_tumor_baf=False,
        draw_depthratio=False,
        draw_normal_baf=False,
        draw_normal_depth=False,
        draw_tumor_depth=False,
        draw_feature=False,

        use_plain_depthratio=False,

        tumor_baf_ylabel=None,
        depthratio_ylabel=None,
        normal_baf_ylabel=None,
        normal_depth_ylabel=None,
        tumor_depth_ylabel=None,
        feature_ylabel=None,
        tumor_baf_ylabel_kwargs=dict(),
        depthratio_ylabel_kwargs=dict(),
        normal_baf_ylabel_kwargs=dict(),
        normal_depth_ylabel_kwargs=dict(),
        tumor_depth_ylabel_kwargs=dict(),
        feature_ylabel_kwargs=dict(),

        n_xlabel=None,
        depthratio_ymax=None,
        CN_ymax=None,
        subCN_ymax=None,

        is_rawdepth=True,
        depth_binsize=1000,
        depth_ymax=None,
        use_upscaled_depthratio=True,
        use_upscaled_tumor_depth=True,
        use_upscaled_normal_depth=True,

        depthratio_dot_kwargs=dict(),
        depthratio_line_segmean_kwargs=dict(),
        depthratio_line_predict_kwargs=dict(),
        depthratio_line_predict_clonal_kwargs=dict(),

        baf_dot_kwargs=dict(),
        baf_line_segmean_kwargs=dict(),
        baf_line_corr_segmean_kwargs=dict(),
        baf_line_predict_kwargs=dict(),
        baf_line_predict_clonal_kwargs=dict(),

        feature_text_kwargs=dict(),
        feature_line_kwargs=dict(),

        CN_line_CNt_kwargs=dict(),
        subCN_line_CNt_kwargs=dict(),
        ccf_bar_kwargs=dict(),

        normal_depth_dot_kwargs=dict(),
        tumor_depth_dot_kwargs=dict(),

        feature_df=None,
        draw_feature_label=True,
        #feature_as_dot=False,

        #depth_ratio_diff=None,
        #baf_diff=0.05,

        #min_N_CNt_candidates=5,
        #N_CNt_candidates_fraction=0.5,

        #ccf_bw=0.1,

        #update_plotdata=False,
        #CNt_diff_factor=0.1,

        #mark_unfit_regions=False,

        split_spines=True,
        merge_same_chroms=True,
    ):
        # make plotdata
        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
        if draw_tumor_baf:
            self.make_tumor_baf_plotdata(sampleid)
        if draw_depthratio:
            self.make_depthratio_plotdata(sampleid, use_upscaled=use_upscaled_depthratio)
        if draw_normal_baf:
            self.make_normal_baf_plotdata(sampleid)
        if draw_normal_depth:
            self.make_depth_plotdata(
                sampleid, 
                is_tumor=False, 
                use_upscaled=use_upscaled_normal_depth, 
                binsize=depth_binsize,
            )
        if draw_tumor_depth:
            self.make_depth_plotdata(
                sampleid, 
                is_tumor=True, 
                use_upscaled=use_upscaled_tumor_depth, 
                binsize=depth_binsize,
            )
        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')

        # main
        assert not (
            (fig is None)
            and (axd is not None)
        )
        fig_not_given = (fig is None)
        axd_not_given = (axd is None)
        if axd is None:
            fig, axd = self.make_axd(
                figsize=figsize, 
                hspace=hspace, 
                height_ratios=height_ratios,
                draw_depthratio_hist=draw_depthratio_hist, 
                draw_solution=False,
                draw_tumor_baf=draw_tumor_baf,
                draw_depthratio=draw_depthratio,
                draw_normal_baf=draw_normal_baf,
                draw_tumor_depth=draw_tumor_depth,
                draw_normal_depth=draw_normal_depth,
                draw_feature=draw_feature,
                draw_xlabel=(n_xlabel is not None),
                fig=fig,
            )
            
        if fig_not_given:
            if title is None:
                title = f'sample_id={sampleid}, is_female={self.data[sampleid]["is_female"]}'
            fig.suptitle(title, fontsize=title_size, y=title_y)

        if draw_tumor_baf:
            self.draw_baf_ax(
                sampleid, 
                axd['baf'], 
                use_merged_segment=True,

                draw_predicted=False,
                draw_corrected=False,
                draw_segmean=False,

                n_xlabel=n_xlabel,

                is_tumor=True,

                mark_unfit_regions=False,

                dot_kwargs=baf_dot_kwargs,
                line_segmean_kwargs=baf_line_segmean_kwargs,
                line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,

                split_spines=split_spines,

                ylabel=tumor_baf_ylabel,
                ylabel_kwargs=tumor_baf_ylabel_kwargs,

                modify_ax=axd_not_given,
                merge_same_chroms=merge_same_chroms,

                make_raw_plotdata=False,
                make_segment_plotdata=False,
            )
        if draw_depthratio:
            self.draw_depthratio_ax(
                sampleid, 
                axd['depthratio'], 
                use_merged_segment=True,

                use_plain_depthratio=use_plain_depthratio,

                draw_predicted=False,
                draw_segmean=False,
                draw_deviation=False,
                draw_depthratio_peaks=False,

                n_xlabel=n_xlabel,
                dot_kwargs=depthratio_dot_kwargs,
                line_segmean_kwargs=depthratio_line_segmean_kwargs,
                ymax=depthratio_ymax,

                split_spines=split_spines,

                ylabel=depthratio_ylabel,
                ylabel_kwargs=depthratio_ylabel_kwargs,

                modify_ax=axd_not_given,
                merge_same_chroms=merge_same_chroms,

                make_raw_plotdata=False,
                make_segment_plotdata=False,
            )
        if draw_normal_baf:
            self.draw_baf_ax(
                sampleid, 
                axd['normal_baf'], 

                use_merged_segment=True,

                draw_predicted=False,
                draw_corrected=False,
                draw_segmean=False,

                n_xlabel=n_xlabel,

                is_tumor=False,

                mark_unfit_regions=False,

                dot_kwargs=baf_dot_kwargs,
                line_segmean_kwargs=baf_line_segmean_kwargs,
                line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,

                split_spines=split_spines,

                ylabel=normal_baf_ylabel,
                ylabel_kwargs=normal_baf_ylabel_kwargs,

                modify_ax=axd_not_given,
                merge_same_chroms=merge_same_chroms,

                make_raw_plotdata=False,
                make_segment_plotdata=False,
            )
        if draw_normal_depth:
            self.draw_depth_ax(
                sampleid,
                axd['normal_depth'],
                n_xlabel=n_xlabel,
                is_tumor=False,
                is_rawdepth=is_rawdepth,
                dot_kwargs=normal_depth_dot_kwargs,
                ymax=depth_ymax,

                split_spines=split_spines,

                ylabel=normal_depth_ylabel,
                ylabel_kwargs=normal_depth_ylabel_kwargs,

                modify_ax=axd_not_given,
                merge_same_chroms=merge_same_chroms,

                make_raw_plotdata=False,
            )
        if draw_tumor_depth:
            self.draw_depth_ax(
                sampleid,
                axd['tumor_depth'],
                n_xlabel=n_xlabel,
                is_tumor=True,
                is_rawdepth=is_rawdepth,
                dot_kwargs=tumor_depth_dot_kwargs,
                ymax=depth_ymax,

                split_spines=split_spines,

                ylabel=tumor_depth_ylabel,
                ylabel_kwargs=tumor_depth_ylabel_kwargs,

                modify_ax=axd_not_given,
                merge_same_chroms=merge_same_chroms,

                make_raw_plotdata=False,
            )
        if draw_feature:
            self.draw_feature_ax(
                axd['feature'],
                feature_df=feature_df,

                n_xlabel=n_xlabel,

                ylabel=feature_ylabel,
                ylabel_kwargs=feature_ylabel_kwargs,

                text_kwargs=feature_text_kwargs,
                line_kwargs=feature_line_kwargs,

                split_spines=split_spines,
                merge_same_chroms=merge_same_chroms,

                #feature_as_dot=feature_as_dot,
                draw_label=draw_feature_label,
            )

        return fig, axd

    ############
    # plotting #
    ############

    def make_axd(
        self, 
        figsize=None, 
        hspace=None, 
        height_ratios=None,
        draw_depthratio_hist=False, 
        draw_solution=False,
        draw_tumor_baf=False,
        draw_depthratio=False,
        draw_normal_baf=False,
        draw_tumor_depth=False,
        draw_normal_depth=False,
        draw_feature=False,
        draw_xlabel=False,
        fig=None,
    ):
        # row names
        row_names = list()
        if draw_solution:
            row_names.extend(['ccf', 'subclonal_CN', 'clonal_CN'])
        if draw_tumor_baf:
            row_names.append('baf')
        if draw_depthratio:
            row_names.append('depthratio')
        if draw_normal_baf:
            row_names.append('normal_baf')
        if draw_normal_depth:
            row_names.append('normal_depth')
        if draw_tumor_depth:
            row_names.append('tumor_depth')
        if draw_feature:
            row_names.append('feature')

        if figsize is None:
            figsize = (30, 5 * len(row_names))

        # mosaic
        if draw_depthratio_hist:
            depthratio_idx = row_names.index('depthratio')
            mosaic = [
                (
                    [name, 'empty_upper'] 
                    if idx < depthratio_idx else
                    (
                        [name, 'empty_lower'] 
                        if idx > depthratio_idx else
                        [name, 'depthratio_hist']
                    )
                )
                for idx, name in enumerate(row_names)
            ]
        else:
            mosaic = [[name,] for name in row_names]

        # gridspec_kw
        if hspace is None:
            hspace = (0.4 if draw_xlabel else 0.1)
        gridspec_kw = dict(
            hspace=hspace, 
            height_ratios=height_ratios,
        )
        if draw_depthratio_hist:
            gridspec_kw.update(dict(width_ratios=[1, 0.1], wspace=0.02))

        # result
        if fig is None:
            fig, axd = plt.subplot_mosaic(
                mosaic,
                figsize=figsize,
                gridspec_kw=gridspec_kw,
            )
        else:
            axd = fig.subplot_mosaic(
                mosaic,
                figsize=figsize,
                gridspec_kw=gridspec_kw,
            )

        if draw_depthratio_hist:
            if 'empty_upper' in axd:
                axd['empty_upper'].axis('off')
            if 'empty_lower' in axd:
                axd['empty_lower'].axis('off')

        return fig, axd

    @staticmethod
    def get_yticklabel_size(yticks):
        return min((200 /len(yticks)), 10)

#    def draw_grids(self, ax, ys, line_params=dict(), merge_same_chroms=True):
#        line_params = (
#            dict(color='black', linewidth=0.2, alpha=0.5)
#            | line_params
#        )
#        chroms, start0s, end0s = zip(
#            *self.genomeplotter.cconv.get_chrom_borders(
#                merge_same_chroms=merge_same_chroms,
#            )
#        )
#        for y in ys:
#            ax.hlines(
#                np.repeat(y, len(start0s)), start0s, end0s, 
#                **line_params,
#            )

    #def draw_centromeres(self, ax):
    #    self.genomeplotter.draw_centromeres(ax)

    def draw_centromeres_type2(self, ax):
        self.genomeplotter.draw_centromeres_type2(ax)

    def make_plotdata_basic(self, sampleid):
        self.make_depthratio_plotdata(sampleid)
        self.make_tumor_baf_plotdata(sampleid)
        self.make_segment_plotdata(sampleid)

    def make_plotdata_fordepth(self, sampleid, use_upscaled=True, binsize=100000):
        self.make_normal_baf_plotdata(sampleid)

        LOGGER_INFO.info(f'Beginning tumor depth data processing')
        self.make_depth_plotdata(
            sampleid, is_tumor=True, use_upscaled=use_upscaled, binsize=binsize
        )
        LOGGER_INFO.info(f'Finished tumor depth data processing')

        LOGGER_INFO.info(f'Beginning normal depth data processing')
        self.make_depth_plotdata(
            sampleid, is_tumor=False, use_upscaled=use_upscaled, binsize=binsize,
        )
        LOGGER_INFO.info(f'Finished normal depth data processing')

    def make_depth_plotdata(self, sampleid, is_tumor, use_upscaled=True, binsize=1000):
        # set params
        datakey = ('tumor_depth' if is_tumor else 'normal_depth')
        plotdata_key = ('tumor_depth_plotdata' if is_tumor else 'normal_depth_plotdata')
        #y_colname = ('mean_depth' if is_rawdepth else 'sequenza_style_norm_mean_depth')

        # upscale raw depth df
        relevant_chroms = [
            x for x in self.genomeplotter.cconv.totalregion_df['Chromosome']
            if not x.startswith('-')
        ]
        original_df = self.data[sampleid][datakey]
        relevant_chroms_df = original_df.loc[
            original_df['Chromosome'].isin(relevant_chroms), 
            :
        ]

        if use_upscaled:
            assert binsize is not None
            input_df = cnvmisc.upsize_depth_df_bin(
                relevant_chroms_df, 
                size=binsize, 
                refver=self.refver,
            )
        else:
            input_df = relevant_chroms_df

        # turn into plotdata
        self.data[sampleid][plotdata_key] = self.genomeplotter.prepare_plot_data(input_df)

    def make_depthratio_plotdata(self, sampleid, use_upscaled=True):
        depthratio_df = (
            self.data[sampleid]['depthratio_upscaled']
            if use_upscaled else
            self.data[sampleid]['depthratio']
        )
        self.data[sampleid]['depthratio_raw_plotdata'] = (
            self.genomeplotter.prepare_plot_data(depthratio_df)
        )

    def make_tumor_baf_plotdata(self, sampleid):
        bafdf = self.data[sampleid]['original_baf']
        tumor_baf_df = bafdf.loc[bafdf['baf_raw_tumor'] > 0, :]
        self.data[sampleid]['baf_raw_plotdata'] = (
            self.genomeplotter.prepare_plot_data(tumor_baf_df)
        )

    def make_normal_baf_plotdata(self, sampleid):
        bafdf = self.data[sampleid]['original_baf']
        normal_baf_df = bafdf.loc[bafdf['baf_raw_normal'] > 0, :]
        self.data[sampleid]['normal_baf_raw_plotdata'] = (
            self.genomeplotter.prepare_plot_data(normal_baf_df)
        )

    def make_segment_plotdata(self, sampleid):
        self.data[sampleid]['merged_segment_plotdata'] = (
            self.genomeplotter.prepare_plot_data(
                self.data[sampleid]['merged_segment']
            )
        )

    def make_plotdata_aftercp_wobaf(self, sampleid, use_saved_plotdata):
        raw_plot_map = {
            'depthratio_upscaled': 'depthratio_raw_plotdata',
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
        helper('merged_segment')

#    def make_plotdata_after_cp(self, sampleid, use_saved_plotdata):
#        raw_plot_map = {
#            'depthratio_upscaled': 'depthratio_raw_plotdata',
#            'baf': 'baf_raw_plotdata',
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
#        helper('baf')
#        helper('merged_segment')

    #def draw_ax_common(
    #    self, ax, **kwargs, 
    #):
    #    self.genomeplotter.draw_ax_common(
    #        ax=ax, **kwargs,
    #    )

#        self.genomeplotter.set_xlim(ax)
#
#        if split_spines:
#            self.genomeplotter.fit_spines_to_regions(
#                ax,
#                draw_chrom_names=True,
#                prefix_with_chr=True,
#                merge_same_chroms=merge_same_chroms,
#            )
#
#        self.draw_grids(ax, ax.get_yticks(), line_params=dict(), merge_same_chroms=merge_same_chroms)
#
#        self.draw_centromeres(ax)
#        #self.draw_centromeres_type2(ax)
#        if n_xlabel is not None:
#            self.genomeplotter.draw_genomecoord_labels(ax, n=n_xlabel)
#        else:
#            ax.set_xticks([])

    def draw_feature_ax(
        self, 
        ax, 
        feature_df,

        n_xlabel=None,

        ylabel=None,
        ylabel_kwargs=dict(),

        #feature_as_dot=False,
        draw_label=True,

        text_kwargs=dict(),
        line_kwargs=dict(),

        split_spines=True,
        merge_same_chroms=True,

        chromlabel_kwargs=dict(), 
    ):
        if ylabel is None:
            ylabel = 'features'

        ax.set_ylim(0, 1)
        ax.set_ylabel(ylabel, **ylabel_kwargs)
        ax.set_yticks([])
        self.draw_ax_common(
            ax, 
            n_xlabel=n_xlabel, 
            split_spines=split_spines,
            merge_same_chroms=merge_same_chroms,
            chromlabel_kwargs=chromlabel_kwargs,
        )

        self.genomeplotter.draw_features(
            ax,
            df=feature_df,

            draw_label=draw_label,
            #feature_as_dot=feature_as_dot,

            y_features=None,
            y_labels=None,

            text_kwargs=text_kwargs,
            line_kwargs=line_kwargs,
        )

    def draw_depth_bundle(
        self, 
        sampleid, axd, n_xlabel, is_rawdepth, 
        depth_dot_kwargs, depth_ymax,
        baf_dot_kwargs, baf_line_segmean_kwargs, baf_line_corr_segmean_kwargs,
    ):
        self.draw_depth_ax(
            sampleid,
            axd['tumor_depth'],
            n_xlabel=n_xlabel,
            is_tumor=True,
            is_rawdepth=is_rawdepth,
            dot_kwargs=depth_dot_kwargs,
            ymax=depth_ymax,

            make_raw_plotdata=False,
        )
        self.draw_depth_ax(
            sampleid,
            axd['normal_depth'],
            n_xlabel=n_xlabel,
            is_tumor=False,
            is_rawdepth=is_rawdepth,
            dot_kwargs=depth_dot_kwargs,
            ymax=depth_ymax,

            make_raw_plotdata=False,
        )
        self.draw_baf_ax(
            sampleid, 
            axd['normal_baf'], 

            use_merged_segment=True,

            draw_predicted=False,
            draw_corrected=False,
            draw_segmean=False,

            n_xlabel=n_xlabel,

            is_tumor=False,

            mark_unfit_regions=False,

            dot_kwargs=baf_dot_kwargs,
            line_segmean_kwargs=baf_line_segmean_kwargs,
            line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,

            make_raw_plotdata=False,
            make_segment_plotdata=False,
        )

    def draw_depth_ax(
        self,
        sampleid,
        ax,
        n_xlabel,
        is_tumor,
        is_rawdepth,
        dot_kwargs=dict(),
        ymax=None,
        add_color=True,
        split_spines=True,
        ylabel=None,
        ylabel_kwargs=dict(),

        make_raw_plotdata=True,
        use_upscaled=True,

        modify_ax=True,
        merge_same_chroms=True,

        chromlabel_kwargs=dict(),
    ):
        # prepare plotdata
        if make_raw_plotdata:
            self.make_depth_plotdata(sampleid, is_tumor=is_tumor, use_upscaled=use_upscaled)

        # set params
        y_colname = ('mean_depth' if is_rawdepth else 'sequenza_style_norm_mean_depth')
        plotdata_key = ('tumor_depth_plotdata' if is_tumor else 'normal_depth_plotdata')
        plotdata_df = self.data[sampleid][plotdata_key]

        # calc dot alpha
        n_dots = plotdata_df.shape[0]
        default_alpha = calc_dot_alpha_depth(n_dots)

        # handle kwargs
        dot_kwargs = (
            #{'color': 'black', 'markersize': 0.3, 'alpha': 0.05}
            {'color': 'black', 'markersize': 0.3, 'alpha': default_alpha}
            | dot_kwargs
        )

        # draw raw data
        if add_color:
            colors = np.repeat('blue', plotdata_df.shape[0])
            colors[np.where(plotdata_df['excluded'])[0]] = 'red'
            self.genomeplotter.draw_dots_scatter(
                ax, 
                df_plotdata=plotdata_df,
                y_colname=y_colname, 
                color_vals=colors,
                plot_kwargs=dot_kwargs,
            )
        else:
            self.genomeplotter.draw_dots(
                ax, 
                df_plotdata=plotdata_df,
                y_colname=y_colname, 
                plot_kwargs=dot_kwargs,
            )

        # set axes attributes
        if modify_ax:
            if ylabel is None:
                ylabel_1 = ('tumor' if is_tumor else 'normal')
                ylabel_2 = ('raw' if is_rawdepth else 'normalized')
                ylabel = f'{ylabel_1} sample {ylabel_2} depth'
            ax.set_ylabel(ylabel, **ylabel_kwargs)

            if ymax is None:
                ymax = np.nanmean(plotdata_df[y_colname]) * 2
            ax.set_ylim(-ymax * 0.1, ymax)
            roundnum = (1 if is_rawdepth else 2)
            yticks = np.round(np.linspace(0, ymax, 10), roundnum)

            ax.set_yticks(yticks)
            ax.set_yticklabels(yticks)

            self.draw_ax_common(
                ax, 
                n_xlabel=n_xlabel, 
                split_spines=split_spines, 
                merge_same_chroms=merge_same_chroms,
                chromlabel_kwargs=chromlabel_kwargs,
            )

    def draw_depthratio_ax(
        self, 
        sampleid, 
        ax, 

        use_plain_depthratio=False,
        use_merged_segment=True, 

        cellularity=None,
        ploidy=None,

        draw_integerCN_lines=False,

        draw_predicted=True,
        draw_segmean=True,
        draw_deviation=False,
        draw_depthratio_peaks=False,

        std_factor=1,

        mark_unfit_regions=False,

        n_xlabel=None,
        dot_kwargs=dict(),
        line_segmean_kwargs=dict(),
        line_predict_kwargs=dict(),
        line_predict_clonal_kwargs=dict(),
        ymax=None,

        split_spines=True,

        ylabel=None,
        ylabel_kwargs=dict(),

        make_raw_plotdata=True,
        make_segment_plotdata=True,
        use_upscaled=True,

        modify_ax=True,
        merge_same_chroms=True,

        chromlabel_kwargs=dict(),
    ):
        # determine if segment information is plotted
        draw_segment_info = any([
            draw_predicted,
            draw_segmean,
        ])

        # prepare plotdata
        if make_raw_plotdata:
            self.make_depthratio_plotdata(sampleid, use_upscaled=use_upscaled)
        if make_segment_plotdata:
            self.make_segment_plotdata(sampleid)

        # calc dot alpha
        raw_plotdata = self.data[sampleid]['depthratio_raw_plotdata']
        n_dots = raw_plotdata.shape[0]
        default_alpha = calc_dot_alpha_depth(n_dots)

        # handle kwargs
        dot_kwargs = (
            {'color': 'black', 'markersize': 0.3, 'alpha': default_alpha}
            | dot_kwargs
        )
        line_segmean_kwargs = (
            {'color': 'tab:blue', 'linewidth': 5, 'alpha': 0.6}
            | line_segmean_kwargs
        )
        line_predict_kwargs = (
            {'color': 'tab:red', 'linewidth': 1.5, 'alpha': 0.6}
            | line_predict_kwargs
        )
        line_predict_clonal_kwargs = (
            {'color': 'tab:orange', 'linewidth': 1.5, 'alpha': 0.6}
            | line_predict_clonal_kwargs
        )

        # draw raw data
        if use_plain_depthratio:
            self.genomeplotter.draw_dots(
                ax, 
                y_colname='depth_ratio_plain', 
                df_plotdata=raw_plotdata, 
                plot_kwargs=dot_kwargs,
            )
        else:
            self.genomeplotter.draw_dots(
                ax, 
                y_colname='depthratio_raw', 
                df_plotdata=raw_plotdata, 
                plot_kwargs=dot_kwargs,
            )

        # draw segment information
        if draw_segment_info:
            if use_merged_segment:
                segment_plotdata = self.data[sampleid]['merged_segment_plotdata']
            else:
                segment_plotdata = self.data[sampleid]['depthratio_segment_plotdata']

            if draw_segmean:
                self.genomeplotter.draw_hlines(
                    ax, 
                    df_plotdata=segment_plotdata,
                    y_colname='depthratio_segment_mean', 
                    plot_kwargs=line_segmean_kwargs,
                )
            if draw_predicted:
                self.genomeplotter.draw_hlines(
                    ax, 
                    df_plotdata=segment_plotdata,
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
            if draw_deviation:
                self.genomeplotter.draw_bgcolors(
                    ax,
                    df_plotdata=segment_plotdata,
                    ymins=(
                        segment_plotdata['depthratio_segment_mean']
                        - std_factor * segment_plotdata['depthratio_segment_std']
                    ),
                    ymaxs=(
                        segment_plotdata['depthratio_segment_mean']
                        + std_factor * segment_plotdata['depthratio_segment_std']
                    ),
                    colors='yellow',
                    plot_kwargs=dict(alpha=0.4),
                )

            if mark_unfit_regions:
                self.draw_unfit_region(ax, segment_plotdata)

        # integer CN lines
        if draw_integerCN_lines:
            assert ((cellularity is not None) and (ploidy is not None))
            self.draw_depthratio_ax_integer_CNs(ax, cellularity, ploidy, sampleid)

        # modify axes
        if modify_ax:
            if ylabel is None:
                ylabel = 'depth ratio'
            ax.set_ylabel(ylabel, **ylabel_kwargs)

            if ymax is None:
                df = self.data[sampleid]['depthratio']
                y_values = df['depthratio_raw'].loc[~df['excluded']].dropna()
                ymax = y_values.quantile(0.999)

                ax.set_ylim(-ymax * 0.1, ymax)
                yticks = np.round(np.arange(0, ymax + 0.25, 0.25), 2)
            else:
                ax.set_ylim(-0.1, ymax)
                yticks = np.round(np.linspace(0, ymax, 8), 2)

            ax.set_yticks(yticks)
            ax.set_yticklabels(yticks, size=self.get_yticklabel_size(yticks))

            self.draw_ax_common(
                ax,
                n_xlabel=n_xlabel, 
                split_spines=split_spines,
                merge_same_chroms=merge_same_chroms,
                chromlabel_kwargs=chromlabel_kwargs,
            )

    def draw_depthratio_ax_integer_CNs(self, ax, cellularity, ploidy, sampleid):
        plotregion_df_withCN = self.add_CNn_to_segment(
            segment_df=self.genomeplotter.region_df,
            mode=self.data[sampleid]['mode'],
            refver=self.refver,
            is_female=self.data[sampleid]['is_female'],
            target_region=self.data[sampleid]['target_region'],
        )
        plotdata = self.genomeplotter.prepare_plot_data(plotregion_df_withCN)

        integer_depthratios = cnvmisc.theoretical_depth_ratio(
            CNt=np.arange(0, 10, 1)[:, np.newaxis],
            cellularity=cellularity,
            tumor_ploidy=ploidy,
            CNn=plotregion_df_withCN['CNn'].to_numpy()[np.newaxis, :],
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
        )  # ndim == 2
        for ys in integer_depthratios:
            plotdata['integer_depthratio'] = ys
            self.genomeplotter.draw_hlines(
                ax, 'integer_depthratio',
                df_plotdata=plotdata,
                plot_kwargs=dict(linewidth=0.7, zorder=0, color='black', alpha=0.7),
            )

    def draw_unfit_region(self, ax, plotdata):
        self.genomeplotter.draw_bgcolors(
            ax,
            df_plotdata=plotdata.loc[plotdata['polyploid_unfit'], :],
            colors='yellow',
            plot_kwargs=dict(alpha=0.2),
        )
        self.genomeplotter.draw_bgcolors(
            ax,
            df_plotdata=plotdata.loc[plotdata['polyploid_unfit_bafonly'], :],
            colors='green',
            plot_kwargs=dict(alpha=0.2),
        )
        self.genomeplotter.draw_bgcolors(
            ax,
            df_plotdata=plotdata.loc[plotdata['monoploid_unfit'], :],
            colors='blue',
            plot_kwargs=dict(alpha=0.2),
        )

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
        histpeaks = cnvmisc.get_hist_peaks(
            depthratio_list, 
            weights=weights,
            bins=np.arange(0, depthratio_list.max(), 0.01),
            threshold=peak_threshold,
        )
        peak_depthratios = histpeaks['peak_values']

        # draw data
        ax.fill_betweenx(histpeaks['bin_midpoints'], histpeaks['hist'], **plot_kwargs)
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

        use_merged_segment=True, 

        draw_corrected=False,
        draw_predicted=False,
        draw_segmean=False,

        n_xlabel=None,

        is_tumor=True,

        mark_unfit_regions=False,

        dot_kwargs=dict(),
        line_segmean_kwargs=dict(),
        line_corr_segmean_kwargs=dict(),
        line_predict_kwargs=dict(),
        line_predict_clonal_kwargs=dict(),

        split_spines=True,

        ylabel=None,
        ylabel_kwargs=dict(),

        make_raw_plotdata=True,
        make_segment_plotdata=True,

        modify_ax=True,
        merge_same_chroms=True,

        chromlabel_kwargs=dict(),
    ):
        # prepare plotdata
        if make_raw_plotdata:
            if is_tumor:
                self.make_tumor_baf_plotdata(sampleid)
            else:
                self.make_normal_baf_plotdata(sampleid)
        if make_segment_plotdata:
            self.make_segment_plotdata(sampleid)

        # find plotdata
        if is_tumor:
            raw_plotdata = self.data[sampleid]['baf_raw_plotdata']
        else:
            raw_plotdata = self.data[sampleid]['normal_baf_raw_plotdata']

        # calc default alpha
        n_dots = raw_plotdata.shape[0]
        default_alpha = calc_dot_alpha_baf(n_dots)

        # handle kwargs
        dot_kwargs = (
            #{'color': 'black', 'markersize': 0.3, 'alpha': 0.01}
            {'color': 'black', 'markersize': 0.3, 'alpha': default_alpha}
            | dot_kwargs
        )
        line_segmean_kwargs = (
            {'color': 'tab:blue', 'linewidth': 2, 'alpha': 0.4}
            | line_segmean_kwargs
        )
        line_corr_segmean_kwargs = (
            {'color': 'tab:green', 'linewidth': 2, 'alpha': 0.4}
            | line_corr_segmean_kwargs
        )
        line_predict_kwargs = (
            {'color': 'tab:red', 'linewidth': 2, 'alpha': 0.4}
            | line_predict_kwargs
        )
        line_predict_clonal_kwargs = (
            {'color': 'tab:orange', 'linewidth': 2, 'alpha': 0.4}
            | line_predict_clonal_kwargs
        )

        # draw raw data
        if is_tumor:
            self.genomeplotter.draw_dots(
                ax, 
                df_plotdata=raw_plotdata, 
                y_colname='baf_raw_tumor', 
                plot_kwargs=dot_kwargs,
            )
        else:
            self.genomeplotter.draw_dots(
                ax, 
                df_plotdata=raw_plotdata, 
                y_colname='baf_raw_normal', 
                plot_kwargs=dot_kwargs,
            )

        # draw segment information
        draw_segment_info = any([
            draw_corrected,
            draw_predicted,
            draw_segmean,
        ])

        if draw_segment_info:
            if is_tumor:
                segment_plotdata = self.data[sampleid]['merged_segment_plotdata']
            else:
                segment_plotdata = self.data[sampleid]['normal_baf_segment_plotdata']

            if draw_segmean:
                self.genomeplotter.draw_hlines(
                    ax, 
                    df_plotdata=segment_plotdata,
                    y_colname='baf_segment_peak', 
                    #y_colname='baf_segment_mean', 
                    plot_kwargs=line_segmean_kwargs,
                )
            if draw_corrected:
                self.genomeplotter.draw_hlines(
                    ax, 
                    df_plotdata=segment_plotdata,
                    y_colname='corrected_baf_segment_mean', 
                    plot_kwargs=line_corr_segmean_kwargs,
                )

            if draw_predicted:
                self.genomeplotter.draw_hlines(
                    ax, 
                    df_plotdata=segment_plotdata,
                    y_colname='baf_predicted', 
                    plot_kwargs=line_predict_kwargs,
                )

            if mark_unfit_regions:
                self.draw_unfit_region(ax, segment_plotdata)

        # set axes attributes
        if ylabel is None:
            ylabel = 'tumor baf' if is_tumor else 'normal baf'
        ax.set_ylabel(ylabel, **ylabel_kwargs)
        ax.set_ylim(-0.6 * 0.1, 0.6)

        yticks = np.round(np.arange(0, 0.6, 0.1), 1)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks, size=self.get_yticklabel_size(yticks))

        if modify_ax:
            self.draw_ax_common(
                ax, 
                n_xlabel=n_xlabel, 
                split_spines=split_spines, 
                merge_same_chroms=merge_same_chroms,
                chromlabel_kwargs=chromlabel_kwargs,
            )

    def draw_CN_ax(
        self, 
        sampleid, 
        ax,
        n_xlabel=None,
        line_A_kwargs=dict(),
        line_B_kwargs=dict(),
        line_CNt_kwargs=dict(),
        ymax=None,
        draw_CNt=True,
        draw_A=True,
        draw_B=True,

        split_spines=True,

        ylabel=None,
        ylabel_kwargs=dict(),

        merge_same_chroms=True,
        chromlabel_kwargs=dict(),
    ):
        # handle kwargs
        line_CNt_kwargs = (
            {'color': 'black'}
            | line_CNt_kwargs
        )
        line_A_kwargs = (
            {'color': 'red'}
            | line_A_kwargs
        )
        line_B_kwargs = (
            {'color': 'blue'}
            | line_B_kwargs
        )

        # draw data
        if draw_CNt:
            self.genomeplotter.draw_hlines(
                ax, 
                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
                y_colname='clonal_CNt', 
                offset=0.1,
                plot_kwargs=line_CNt_kwargs,
            )
        if draw_A:
            self.genomeplotter.draw_hlines(
                ax, 
                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
                y_colname='clonal_A', 
                offset=0,
                plot_kwargs=line_A_kwargs,
            )
        if draw_B:
            self.genomeplotter.draw_hlines(
                ax, 
                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
                y_colname='clonal_B', 
                offset=-0.1,
                plot_kwargs=line_B_kwargs,
            )

        # set axes attributes
        if ylabel is None:
            ylabel = 'tumor clonal copy number'
        ax.set_ylabel(ylabel, **ylabel_kwargs)

        if ymax is None:
            #weights = (
            #    self.data[sampleid]['merged_segment']['End'] 
            #    - self.data[sampleid]['merged_segment']['Start']
            #)
            #CNts = self.data[sampleid]['merged_segment']['clonal_CNt']
            #ymax = int(common.nanaverage(CNts, weights)) * 2

            df = self.data[sampleid]['merged_segment']
            ymax = df['clonal_CNt'].quantile(0.99)
        else:
            pass

        max_ticknum = 15
        step = np.ceil(ymax / 15).astype(int)
        yticks = np.arange(0, ymax, step).astype(int)

        ax.set_ylim(-0.5, ymax)
        ax.set_yticks(yticks, minor=False)
        #ax.set_yticks(np.arange(0, ymax, 1).astype(int), minor=True)
        ax.set_yticklabels(yticks, size=self.get_yticklabel_size(yticks))

        self.draw_ax_common(
            ax, 
            n_xlabel=n_xlabel, 
            split_spines=split_spines, 
            merge_same_chroms=merge_same_chroms,
            chromlabel_kwargs=chromlabel_kwargs,
        )

    def draw_subclonal_CN_ax(
        self, 
        sampleid, 
        ax,
        n_xlabel=None,
        line_A_kwargs=dict(),
        line_B_kwargs=dict(),
        line_CNt_kwargs=dict(),
        ymax=None,
        draw_CNt=True,
        draw_A=True,
        draw_B=True,

        split_spines=True,

        ylabel=None,
        ylabel_kwargs=dict(),

        merge_same_chroms=True,
        chromlabel_kwargs=dict(),
    ):
        # handle kwargs
        line_CNt_kwargs = (
            {'color': 'black'}
            | line_CNt_kwargs
        )
        line_A_kwargs = (
            {'color': 'red'}
            | line_A_kwargs
        )
        line_B_kwargs = (
            {'color': 'blue'}
            | line_B_kwargs
        )

        # draw CNt
        if draw_CNt:
            self.genomeplotter.draw_hlines(
                ax, 
                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
                y_colname='subclonal_CNt', 
                offset=0.1,
                plot_kwargs=line_CNt_kwargs,
            )
        if draw_A:
            self.genomeplotter.draw_hlines(
                ax, 
                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
                y_colname='subclonal_A', 
                offset=0,
                plot_kwargs=line_A_kwargs,
            )
        if draw_B:
            self.genomeplotter.draw_hlines(
                ax, 
                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
                y_colname='subclonal_B', 
                offset=-0.1,
                plot_kwargs=line_B_kwargs,
            )

        # set axes attributes
        if ylabel is None:
            ylabel = 'tumor subclonal copy number'
        ax.set_ylabel(ylabel, **ylabel_kwargs)

        if ymax is None:
            #weights = (
            #    self.data[sampleid]['merged_segment']['End'] 
            #    - self.data[sampleid]['merged_segment']['Start']
            #)
            #CNts = self.data[sampleid]['merged_segment']['subclonal_CNt']
            #ymax = int(common.nanaverage(CNts, weights)) * 2

            df = self.data[sampleid]['merged_segment']
            ymax = df['subclonal_CNt'].quantile(0.99)
        else:
            pass

        max_ticknum = 15
        step = np.ceil(ymax / 15).astype(int)
        yticks = np.arange(0, ymax, step).astype(int)

        ax.set_ylim(-0.5, ymax)
        ax.set_yticks(yticks, minor=False)
        #ax.set_yticks(np.arange(0, ymax, 1).astype(int), minor=True)
        ax.set_yticklabels(yticks, size=self.get_yticklabel_size(yticks))

        self.draw_ax_common(
            ax, 
            n_xlabel=n_xlabel, 
            split_spines=split_spines, 
            merge_same_chroms=merge_same_chroms,
            chromlabel_kwargs=chromlabel_kwargs,
        )

    def draw_ccf_ax(
        self, 
        sampleid, 
        ax,
        n_xlabel=None,
        bar_kwargs=dict(),

        split_spines=True,

        ylabel=None,
        ylabel_kwargs=dict(),

        merge_same_chroms=True,
        chromlabel_kwargs=dict(),
    ):
        bar_kwargs = (
            dict()
            | bar_kwargs
        )

        plotdata = self.data[sampleid]['merged_segment_plotdata']

        self.genomeplotter.draw_bars(
            ax, 
            y_colname='ccf', 
            df_plotdata=plotdata,
            plot_kwargs=bar_kwargs,
        )

        if ylabel is None:
            ylabel = 'subclonal fraction'
        ax.set_ylabel(ylabel, **ylabel_kwargs)

        ax.set_ylim(-0.05, 1.05)
        yticks = np.round(np.arange(0, 1.1, 0.1), 1)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)

        self.genomeplotter.draw_bgcolors(
            ax,
            df_plotdata=plotdata.loc[plotdata['ccf'].isna(), :],
            colors='green',
            plot_kwargs=dict(alpha=0.1),
        )

        self.draw_ax_common(
            ax, 
            n_xlabel=n_xlabel, 
            split_spines=split_spines, 
            merge_same_chroms=merge_same_chroms,
            chromlabel_kwargs=chromlabel_kwargs,
        )

    def show_ccfs_main(
        self,
        ccfs,
        lengths,
        density,
        peak_values,
        draw_dots=False,
    ):
        fig, ax = plt.subplots(figsize=(30, 5))
        ax.hist(
            ccfs, 
            range=(0, 1), 
            bins=50, 
            weights=lengths,
            density=True,
        )

        density_xs = np.arange(0, 1, 0.01)
        density_ys = density(density_xs)
        ax.plot(density_xs, density_ys)

        for x in peak_values:
            ax.axvline(x, color='red')

        if draw_dots:
            ylim = ax.get_ylim()
            dot_xs = ccfs
            dot_ys = scipy.stats.uniform.rvs(loc=ylim[0], scale=(ylim[1] - ylim[0]), size=len(ccfs))
            colors = fit.labels_
            ax.scatter(dot_xs, dot_ys, c=colors, alpha=0.7, s=4)

        return fig, ax

#    def plot_base_wrapper(
#        self,
#        basemethod, 
#
#        region_chroms=None, 
#        region_start0s=None, 
#        region_end0s=None, 
#        weights=None,
#        region_gaps=None,
#
#        **kwargs,
#    ):
#        if (
#            (region_gaps is not None)
#            or (region_chroms is not None)
#        ):
#            # make new genomeplotter
#            if region_chroms is None:
#                new_region_df = self.genomeplotter.region_df.copy()
#            else:
#                new_region_df = GenomePlotter.make_new_region_df(
#                    self.refver, 
#                    region_chroms, 
#                    region_start0s, 
#                    region_end0s, 
#                    weights,
#                )
#
#            new_genomeplotter = GenomePlotter(
#                self.refver, 
#                region_df=new_region_df,
#                region_gaps=region_gaps,
#            )
#
#            # switch and run
#            old_genomeplotter = self.genomeplotter
#            self.genomeplotter = new_genomeplotter
#
#            fig, axd = basemethod(**kwargs)
#
#            self.genomeplotter = old_genomeplotter
#        else:
#            fig, axd = basemethod(**kwargs)
#
#        return fig, axd

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
            if 'normal_depth' in self.data[sampleid]:
                target_region_gr = self.find_germline_copyneutral_region(
                    sampleid, factors=(0.8, 1.2),
                )
            else:
                target_region_gr = common.DEFAULT_CHROMDICTS[self.refver].to_gr(
                    assembled_only=True, as_gr=True,
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

    def find_germline_copyneutral_region(self, sampleid, factors=(0.3, 2)):
        depth_df = self.data[sampleid]['normal_depth']

        # upscale depth bins for segmentation speed
        upscaled_depth_df = cnvmisc.upsize_depth_df_bin(
            depth_df, size=1000, refver=self.refver,
        )

        # run segmentation
        LOGGER_INFO.info('Running segmentation of normal bam depth in order to find copy-neutral region')
        segdf, _ = rcopynumber.run_rcopynumber(
            depth_df=upscaled_depth_df.rename(columns={'mean_depth': 'depth_raw'}), 
            refver=self.refver, 
            as_gr=False, 
            winsorize=False, 
            verbose=False,
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

    ###############################
    # helpers of segment creation #
    ###############################

#    def make_depth_segment(
#        self,
#        sampleid,
#        winsorize=False,
#        gamma=None,
#        kmin=None,
#        verbose=True,
#    ):
#        self.data[sampleid]['depthratio_segment'] = _make_depth_segment(
#            self.data[sampleid]['depthratio_upscaled'],
#            mode=self.data[sampleid]['mode'],
#            refver=self.refver,
#            winsorize=winsorize,
#            gamma=gamma,
#            kmin=kmin,
#            verbose=verbose,
#        )
#
#    def make_baf_segment(
#        self,
#        sampleid,
#        winsorize=False,
#        gamma=None,
#        kmin=None,
#        verbose=False,
#        bafcorrector=libbaf.load_bafcorrect_func(),
#    ):
#        self.data[sampleid]['baf_segment'] = _make_baf_segment(
#            baf_df=self.data[sampleid]['baf'],
#            target_region=self.data[sampleid]['target_region'],
#            mode=self.data[sampleid]['mode'],
#            refver=self.refver,
#            winsorize=winsorize,
#            gamma=gamma,
#            kmin=kmin,
#            verbose=verbose,
#            bafcorrector=bafcorrector,
#        )

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

        # reduce by merging adjacent segments with identical annotation values
        merged_segment = self.deduplicate_merged_segment(merged_segment)

        # add CNn
        merged_segment = self.add_CNn_to_segment(
            merged_segment,
            self.data[sampleid]['mode'],
            self.refver,
            self.data[sampleid]['is_female'],
            self.data[sampleid]['target_region'],
        )

        # add std and std/mean ratio
        #self.add_depthratio_std_to_segment(sampleid)

        # fit to target region
        #merged_segment = self.fit_segment_to_targetregion(sampleid, merged_segment)

        self.data[sampleid]['merged_segment'] = merged_segment

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

    def deduplicate_merged_segment(self, merged_segment):
        merged_segment = cnvmisc.sort_genome_df(merged_segment, self.refver)

        chromdict = common.DEFAULT_CHROMDICTS[self.refver]
        annot_arr = np.concatenate(
            [
                merged_segment['Chromosome'].apply(chromdict.contigs.index).to_numpy()[:, np.newaxis],
                merged_segment.loc[:, ['depthratio_segment_mean', 'baf_segment_mean']].to_numpy(),
            ],
            axis=1,
        )
        values, counts, groupkey = common.array_grouper(annot_arr, omit_values=True)
        groupby = merged_segment.groupby(groupkey)
        dedup_segdf = groupby.first()
        dedup_segdf['End'] = groupby['End'].last()

        return dedup_segdf

#    def postprocess_segment(self, sampleid, cellularity, ploidy):
#        # add CNt and B
#        cpinfo = self.data[sampleid]['cpscores'][(cellularity, ploidy)]
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

    ####################
    # solution finding #
    ####################

    def make_CN_solution_freeccf(
        self,
        sampleid,
        cellularity,
        ploidy,
        depth_ratio_diff=None,
        baf_diff=None,
        min_N_CNt_candidates=5,
        N_CNt_candidates_fraction=0.5,
        limited_clonal=True,
    ):
        segdf = self.data[sampleid]['merged_segment']
        (
            clonal_solution, 
            flags, 
            freeccf_solution,
            calculated_depth_ratio, 
            calculated_baf,
            average_CNt,
        ) = cnvmisc.find_solution_before_ccfs(
            depth_ratio=segdf['depthratio_segment_mean'], 
            baf=segdf['corrected_baf_segment_mean'],
            CNn=segdf['CNn'],
            lengths=(segdf['End'] - segdf['Start']),
            cellularity=cellularity,
            tumor_ploidy=ploidy,
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
            depth_ratio_diff=depth_ratio_diff,
            baf_diff=baf_diff,
            Bn=1,
            min_N_CNt_candidates=min_N_CNt_candidates,
            N_CNt_candidates_fraction=N_CNt_candidates_fraction,

            limited_clonal=limited_clonal,
        )

        freeccf_result = {
            #'fixed_ccfs': fixed_ccfs, 
            #'ccf_plotdata': ccf_plotdata, 
            'clonal_solution': clonal_solution, 
            'flags': flags, 
            'freeccf_solution': freeccf_solution,
            'calculated_depth_ratio': calculated_depth_ratio, 
            'calculated_baf': calculated_baf,
        }
        self.data[sampleid]['freeccf_result'] = freeccf_result
        self.data[sampleid]['average_CNt'] = average_CNt

    def select_fixed_ccfs(self, sampleid, bandwidth=0.1):
        segdf = self.data[sampleid]['merged_segment']
        lengths = (segdf['End'] - segdf['Start']).to_numpy()
        fixed_ccfs, ccf_plotdata = cnvmisc.select_fixed_ccfs(
            freeccf_solution=self.data[sampleid]['freeccf_result']['freeccf_solution'], 
            lengths=lengths, 
            flags=self.data[sampleid]['freeccf_result']['flags'], 
            bandwidth=bandwidth,
        )
        self.data[sampleid]['fixed_ccfs'] = fixed_ccfs
        self.data[sampleid]['ccf_plotdata'] = ccf_plotdata

    def make_CN_solution_after_ccfs(
        self,
        sampleid,
        cellularity,
        ploidy,
        min_N_CNt_candidates=5,
        N_CNt_candidates_fraction=0.5,
        CNt_diff_factor=0.1,
        limited_clonal=True,
    ):
        segdf = self.data[sampleid]['merged_segment']
        solution = cnvmisc.find_solution_after_ccfs(
            depth_ratio=segdf['depthratio_segment_mean'].to_numpy(),
            baf=segdf['corrected_baf_segment_mean'].to_numpy(),
            CNn=segdf['CNn'].to_numpy(),
            lengths=(segdf['End'] - segdf['Start']).to_numpy(),
            cellularity=cellularity,
            tumor_ploidy=ploidy,
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
            average_CNt=self.data[sampleid]['average_CNt'],

            fixed_ccfs=self.data[sampleid]['fixed_ccfs'], 
            clonal_solution=self.data[sampleid]['freeccf_result']['clonal_solution'], 
            flags=self.data[sampleid]['freeccf_result']['flags'],

            Bn=1,
            min_N_CNt_candidates=min_N_CNt_candidates,
            N_CNt_candidates_fraction=N_CNt_candidates_fraction,

            CNt_diff_factor=CNt_diff_factor,

            limited_clonal=limited_clonal,
        )

        self.data[sampleid]['CN_solution'] = solution

#    def make_CN_solution_onestep(
#        self,
#        sampleid,
#        cellularity,
#        ploidy,
#        depth_ratio_diff=None,
#        baf_diff=0.05,
#        min_N_CNt_candidates=5,
#        N_CNt_candidates_fraction=0.5,
#        ccf_bw=0.1,
#    ):
#        self.make_CN_solution_freeccf(
#            sampleid,
#            cellularity,
#            ploidy,
#            depth_ratio_diff=depth_ratio_diff,
#            baf_diff=baf_diff,
#            min_N_CNt_candidates=min_N_CNt_candidates,
#            N_CNt_candidates_fraction=N_CNt_candidates_fraction,
#        )
#        self.select_fixed_ccfs(sampleid, bandwidth=ccf_bw)
#        self.make_CN_solution_after_ccfs(
#            sampleid,
#            cellularity,
#            ploidy,
#            min_N_CNt_candidates=min_N_CNt_candidates,
#            N_CNt_candidates_fraction=N_CNt_candidates_fraction,
#        )

    def add_freeccf_solution_to_segment(self, sampleid):
        #subclonal_theoretical_depthratio = self.data[sampleid]['freeccf_result']['calculated_depth_ratio']
        #subclonal_theoretical_baf = self.data[sampleid]['freeccf_result']['calculated_baf']
        cnvmisc.add_freeccf_solution_to_segment(
            segment_df=self.data[sampleid]['merged_segment'], 
            clonal_solution=self.data[sampleid]['freeccf_result']['clonal_solution'],
            freeccf_subclonal_solution=self.data[sampleid]['freeccf_result']['freeccf_solution'], 
            flags=self.data[sampleid]['freeccf_result']['flags'], 
        )

    def add_fixedccf_solution_to_segment(self, sampleid):
        cnvmisc.add_fixedccf_solution_to_segment(
            segment_df=self.data[sampleid]['merged_segment'], 
            fixedccf_solution=self.data[sampleid]['CN_solution'],
        )

    def add_solution_to_plotdata(self, sampleid):
        assert (
            self.data[sampleid]['merged_segment'].loc[:, ['Chromosome', 'Start', 'End']]
            == self.data[sampleid]['merged_segment_plotdata'].loc[:, ['Chromosome', 'Start', 'End']]
        ).all(axis=None)

        cnvmisc.add_fixedccf_solution_to_segment(
            segment_df=self.data[sampleid]['merged_segment_plotdata'], 
            fixedccf_solution=self.data[sampleid]['CN_solution'],
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

    def add_CN_pred_to_segment_new(self, sampleid, cellularity, ploidy):
        # calculate CNt values
        depth_ratios = self.data[sampleid]['merged_segment']['depthratio_segment_mean'].to_numpy()
        CNns = self.data[sampleid]['merged_segment']['CNn'].to_numpy()
        clonal_CNts = cnvmisc.calc_clonal_CNt(
            depth_ratio=depth_ratios,
            CNn=CNns,

            cellularity=cellularity,
            tumor_ploidy=ploidy,
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
        )

        subclonal_CNts, ccfs = cnvmisc.calc_subclonal_CNt(
            depth_ratio=depth_ratios,
            clonal_CNt=clonal_CNts,
            CNn=CNns,

            cellularity=cellularity,
            tumor_ploidy=ploidy,
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
            tumor_avg_depth_ratio=1,
            normal_avg_depth_ratio=1,
            only_max_ccf=True,
        )

        # calculate theoretical depths
        predicted_depth_ratios = cnvmisc.theoretical_depth_ratio_subclone(
            clonal_CNt=clonal_CNts, 
            subclonal_CNt=subclonal_CNts,
            ccf=ccfs,
            cellularity=cellularity, 
            tumor_ploidy=ploidy, 
            CNn=CNns, 
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
            tumor_avg_depth_ratio=1, 
            normal_avg_depth_ratio=1,
        )

        predicted_depth_ratios_clonal = cnvmisc.theoretical_depth_ratio(
            CNt=clonal_CNts, 
            cellularity=cellularity, 
            tumor_ploidy=ploidy, 
            CNn=CNns, 
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
            tumor_avg_depth_ratio=1, 
            normal_avg_depth_ratio=1,
        )

        # annotate segment dataframe
        merged_segdf = self.data[sampleid]['merged_segment']

        merged_segdf['clonal_CNt'] = clonal_CNts
        merged_segdf['subclonal_CNt'] = subclonal_CNts
        merged_segdf['ccf'] = ccfs
        merged_segdf['depthratio_predicted'] = predicted_depth_ratios
        merged_segdf['depthratio_predicted_clonal'] = predicted_depth_ratios_clonal

        self.data[sampleid]['merged_segment'] = merged_segdf

        # add standard deviation of depth ratios
        self.add_depthratio_std_to_segment(sampleid)

    def add_depthratio_std_to_segment(self, sampleid):
        self.data[sampleid]['merged_segment'] = add_depthratio_std(
            self.data[sampleid]['merged_segment'], 
            self.data[sampleid]['depthratio_upscaled'],
            self.refver,
        )

    def add_baf_std_to_segment(self, sampleid):
        merged_segdf = self.data[sampleid]['merged_segment']
        right_df = self.data[sampleid]['original_baf']
        right_df = right_df.loc[
            right_df['baf_raw_tumor'] > 0, 
            ['Chromosome', 'Start', 'End', 'baf_raw_tumor'],
        ]
        joined_segdf = pyranges_helper.join(
            merged_segdf,
            right_df,
            how='left', merge='mean', add_std=True, ddof=0,
            sort=True, refver=self.refver,
        )
        joined_segdf.drop('baf_raw_tumor', axis=1, inplace=True)
        joined_segdf.rename(
            columns={'baf_raw_tumor_std': 'baf_segment_std'}, 
            inplace=True,
        )

        self.data[sampleid]['merged_segment'] = joined_segdf

    ####################
    # cp value finding #
    ####################

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

        selected_cols = [
            'Chromosome', 
            'Start', 
            'End', 
            'vaf_raw_tumor', 
            'baf_raw_tumor', 
        ]
        if vcf_sampleid_normal is not None:
            selected_cols.append('vaf_raw_normal')
            selected_cols.append('baf_raw_normal')
        vaf_df = vaf_df.loc[:, selected_cols]
        self.data[sampleid]['original_baf'] = vaf_df

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
    def add_CNn_to_segment(
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


# segmentation functions for running with multiprocessing

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


def _make_depth_segment(
    depthratio_df,
    mode,
    refver,
    winsorize=False,
    gamma=None,
    kmin=None,
    verbose=False,
):
    gamma, kmin = handle_gamma_kmin_args(gamma, kmin, mode)

    if 'depthratio_raw' in depthratio_df.columns:
        input_df = depthratio_df.rename(columns={'depthratio_raw': 'depth_raw'})
    else:
        input_df = depthratio_df

    segment_df, _ = rcopynumber.run_rcopynumber_unified(
        depth_df=input_df,
        refver=refver,

        as_gr=False, 
        winsorize=winsorize,
        compact=(mode == 'panel'), 
        verbose=verbose,
        gamma=gamma,
        kmin=kmin,
    )

    segment_df.rename(columns={'depth_segment_mean': 'depthratio_segment_mean'}, inplace=True)

    segment_df = add_depthratio_std(segment_df, depthratio_df, refver)

    return segment_df


def add_depthratio_std(segment_df, raw_df, refver):
    left_df = segment_df.drop(
        segment_df.columns.intersection(
            ['depthratio_segment_std', 'depthratio_segment_std_mean_ratio'], 
        ),
        axis=1, 
        inplace=False,
    )
    right_df = raw_df.loc[
        :, ['Chromosome', 'Start', 'End', 'depthratio_raw']
    ]
    joined_segdf = pyranges_helper.join(
        left_df,
        right_df,
        how='left', merge='mean', add_std=True, ddof=0,
        sort=True, refver=refver,
    )
    joined_segdf.drop('depthratio_raw', axis=1, inplace=True)
    joined_segdf.rename(
        columns={'depthratio_raw_std': 'depthratio_segment_std'}, 
        inplace=True,
    )

    # add std/mean ratio
    joined_segdf['depthratio_segment_std_mean_ratio'] = (
        joined_segdf['depthratio_segment_std']
        / joined_segdf['depthratio_segment_mean']
    ).to_numpy()

    return joined_segdf


def _make_baf_segment(
    baf_df,
    target_region,
    mode,
    refver,
    winsorize=False,
    gamma=None,
    kmin=None,
    verbose=False,
    baf_cutoff=0.1,
    #bafcorrector=libbaf.load_bafcorrect_func(),
):
    gamma, kmin = handle_gamma_kmin_args(gamma, kmin, mode)

    targetovlp_baf_df = pr.PyRanges(baf_df).overlap(target_region).df
    input_df = targetovlp_baf_df.loc[:, ['Chromosome', 'Start', 'End']]
    input_df['depth_raw'] = targetovlp_baf_df['baf_raw_tumor']
    input_df = input_df.loc[input_df['depth_raw'] > baf_cutoff, :]

    segment_df, _ = rcopynumber.run_rcopynumber_unified(
        depth_df=input_df,
        refver=refver,

        as_gr=False, 
        winsorize=winsorize,
        compact=(mode == 'panel'), 
        verbose=verbose,
        gamma=gamma,
        kmin=kmin,
    )

    segment_df.rename(columns={'depth_segment_mean': 'baf_segment_mean'}, inplace=True)

    return segment_df


def _make_segments_targetfunc_depth(shdict, *args, **kwargs):
    shdict['depth_segment'] = _make_depth_segment(*args, **kwargs)


def _make_segments_targetfunc_baf(shdict, *args, **kwargs):
    shdict['baf_segment'] = _make_baf_segment(*args, **kwargs)


def _make_segments_main(
    depthratio_df,
    mode,
    refver,
    winsorize,
    depthratio_gamma,
    depthratio_kmin,
    baf_gamma,
    baf_kmin,
    verbose,

    baf_df,
    target_region,
    baf_cutoff,
    #bafcorrector,
):
    with multiprocessing.Manager() as manager:
        shdict = manager.dict()
        subp1 = multiprocessing.Process(
            target=_make_segments_targetfunc_depth,
            args=(
                shdict, 
                depthratio_df,
                mode,
                refver,
                winsorize,
                depthratio_gamma,
                depthratio_kmin,
                verbose,
            ),
        )
        subp2 = multiprocessing.Process(
            target=_make_segments_targetfunc_baf,
            args=(
                shdict, 
                baf_df,
                target_region,
                mode,
                refver,
                winsorize,
                baf_gamma,
                baf_kmin,
                verbose,
                baf_cutoff,
                #bafcorrector,
            ),
        )

        subp1.start()
        subp2.start()
        subp1.join()
        subp2.join()

        depth_segment, baf_segment = shdict['depth_segment'], shdict['baf_segment']

    return depth_segment, baf_segment


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
    gplotter.fit_spines_to_regions(axd['CN'])

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
    gplotter.fit_spines_to_regions(axd['baf'])

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

    gplotter.fit_spines_to_regions(axd['depth'])

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



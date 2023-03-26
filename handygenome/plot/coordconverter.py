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
#import handygenome.cnv.mosdepth as libmosdepth


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

#    def isec_trim_data_df(self, df):
#        assert '_index' not in df.columns
#
#        isec_gr = pr.PyRanges(df).intersect(self.totalregion_gr)
#        isec_gr._index = list(range(isec_gr.df.shape[0]))
#
#        subgrs_bychrom = dict()
#        for chrom in isec_gr.Chromosome.unique():
#            subgrs_bychrom[chrom] = isec_gr[chrom]
#
#        return isec_gr, subgrs_bychrom

#    def get_ordered_plot_coords(self, subgrs_bychrom, pos0_colname, nproc=None):
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
#                self.genomic_to_plot_with_indexes(chrom, getattr(subgr, pos0_colname), subgr._index)
#            )
#
#        index_coord_pairs = itertools.chain.from_iterable(zip(*x) for x in result)
#        plot_coords = np.fromiter(
#            (x[1] for x in sorted(index_coord_pairs, key=operator.itemgetter(0))),
#            dtype=np.int_,
#        )
#
#        return plot_coords

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
    #def set_xlim(self, ax):
    #    ax.set_xlim(*self.xlim)

    def get_chrom_borders(self):
        result = list()
        for key, subdf in self.iter_totalregion_df():
            chroms = set(subdf['Chromosome'])
            assert len(chroms) == 1
            chrom = chroms.pop()
            result.append(
                (chrom, subdf['plot_interval_start0s'].iloc[0], subdf['plot_interval_end0s'].iloc[-1])
            )
        return result

#    def draw_chrom_borders(self, ax, color='black', linewidth=1, fontsize=8, **kwargs):
#        border_pos0s = set()
#        xlim = self.xlim
#        chrom_borders = self.get_chrom_borders()
#        for chrom, start0, end0 in chrom_borders:
#            ax.text(
#                (start0 + end0) / 2, 
#                ax.get_ylim()[1], 
#                chrom, 
#                ha='center', va='bottom',
#                size=fontsize,
#            )
#
#        # draw chromosome region borderlines
#        for chrom, start0, end0 in chrom_borders:
#            if start0 != xlim[0]:
#                border_pos0s.add(start0)
#            if end0 != xlim[1]:
#                border_pos0s.add(end0)
#        for pos0 in border_pos0s:
#            ax.axvline(pos0, color=color, linewidth=linewidth, **kwargs)

#    def prepare_plot_data(self, df, nproc=None):
#        # create isec between total region and input data
#        isec_gr, subgrs_bychrom = self.isec_trim_data_df(df)
#        # Add "End_minus1" columns; "End" columns cannot be used for plot coordinate calculation
#        for chrom, subgr in subgrs_bychrom.items():
#            subgr.End_minus1 = subgr.End - 1
#
#        xmins = self.get_ordered_plot_coords(subgrs_bychrom, 'Start', nproc=nproc)
#        xmaxs_minus1 = self.get_ordered_plot_coords(subgrs_bychrom, 'End_minus1', nproc=nproc)
#        xmaxs = xmaxs_minus1 + 1
#
#        return {
#            'isec_gr': isec_gr,
#            'subgrs_bychrom': subgrs_bychrom,
#            'xmins': xmins,
#            'xmaxs': xmaxs,
#        }

#    def draw_hlines(
#        self, ax, y_colname, *, 
#        df=None, df_plotdata=None, offset=0, nproc=None, 
#        plot_kwargs={},
#    ):
#        default_plot_kwargs = {}
#        default_plot_kwargs.update(plot_kwargs)
#
#        if df_plotdata is None:
#            df_plotdata = self.prepare_plot_data(df, nproc=nproc)
#            
#        ys, xmins, xmaxs = self.merge_adjacent_hlines_data(
#            getattr(df_plotdata['isec_gr'], y_colname).to_numpy() + offset,
#            df_plotdata['xmins'],
#            df_plotdata['xmaxs'],
#        )
#        ax.hlines(ys, xmins, xmaxs, **default_plot_kwargs)
#
#    @classmethod
#    def merge_adjacent_hlines_data(cls, ys, xmins, xmaxs):
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
#
#    def draw_dots(
#        self, ax, y_colname, *, 
#        df=None, df_plotdata=None, nproc=None, 
#        plot_kwargs={},
#    ):
#        default_plot_kwargs = {
#            'color': 'black',
#            'marker': 'o',
#            'linestyle': '',
#        }
#        default_plot_kwargs.update(plot_kwargs)
#
#        if df_plotdata is None:
#            df_plotdata = self.prepare_plot_data(df, nproc=nproc)
#
#        xs = (df_plotdata['xmins'] + (df_plotdata['xmaxs'] - 1)) / 2
#        ys = getattr(df_plotdata['isec_gr'], y_colname).to_numpy()
#        ax.plot(xs, ys, **default_plot_kwargs)
#
#    def draw_bgcolors(
#        self, ax, *, 
#        df=None, df_plotdata=None, nproc=None,
#        plot_kwargs={},
#    ):
#        default_plot_kwargs = {
#            'color': 'yellow',
#            'alpha': 0.1,
#            'zorder': 0,
#        }
#        default_plot_kwargs.update(plot_kwargs)
#
#        if df_plotdata is None:
#            df_plotdata = self.prepare_plot_data(df, nproc=nproc)
#
#        ylims = ax.get_ylim()
#        ymin = ylims[0]
#        height = ylims[1] - ylims[0]
#        xmaxs = df_plotdata['xmaxs']
#        xmins = df_plotdata['xmins']
#        widths = xmaxs - xmins
#
#        boxes = [
#            Rectangle((xm, ymin), width=w, height=height)
#            for (xm, w) in zip(xmins, widths)
#        ]
#        ax.add_collection(PatchCollection(boxes, **default_plot_kwargs))


#class CoordConverter_old:
#    def __init__(self, gr):
#        assert isinstance(gr, pr.PyRanges), f'"gr" must be a pyranges.PyRanges object.'
#        self.totalregion_gr = gr.copy()
#        if 'weight' not in self.totalregion_gr.columns:
#            self.totalregion_gr.weight = 1
#        self.chromosomes = set(self.totalregion_gr.keys())
#        self._set_data()
#
#    def _set_data(self):
#        genome_intervals = pd.IntervalIndex.from_arrays(
#            left=list(self.totalregion_gr.Start), right=list(self.totalregion_gr.End), closed='left',
#        )
#        raw_region_lengths = self.totalregion_gr.lengths().reset_index(drop=True)
#        plot_region_lengths = raw_region_lengths * self.totalregion_gr.weight.reset_index(drop=True)
#        cumsum = plot_region_lengths.cumsum()
#        global_offsets = cumsum.shift(1, fill_value=0)
#        plot_intervals = pd.IntervalIndex.from_breaks(
#            ([0] + list(cumsum)), closed='left',
#        )
#
#        self.data = pd.DataFrame(
#            data={
#                'raw_region_lengths': raw_region_lengths.array,
#                'plot_region_lengths': plot_region_lengths.array,
#                'global_offsets': global_offsets.array,
#            },
#            index=pd.MultiIndex.from_arrays(
#                [self.totalregion_gr.Chromosome.array, genome_intervals, plot_intervals],
#                names=('chromosome', 'genome_interval', 'plot_interval'),
#            )
#        )
#
#        #self.genome_intervals = genome_intervals
#        #self.plot_intervals = plot_intervals
#
#    @property
#    def xlim(self):
#        start0 = self.data.index.values[0][2].left
#        end0 = self.data.index.values[-1][2].right
#        return (start0, end0)
#
#    def genomic_to_plot_vectorized(self, chrom, pos0_list):
#        if chrom not in self.totalregion_gr.keys():
#            raise Exception(f'Input "chrom" argument is not included in the plotting region.')
#
#
#        subdata = self.data.loc[chrom, :]
#        genome_intvlist = subdata.index.get_level_values('genome_interval')
#        contains = genome_intvlist.contains(pos0)
#
#        num_hit = contains.sum()
#        if num_hit == 0:
#            return None
#        elif num_hit > 1:
#            raise Exception(f'More than one intervals contains the input position.')
#
#        idx = np.where(contains)[0][0]
#        genome_intv = genome_intvlist[idx]
#        global_offset = subdata.iloc[idx, :].loc['global_offsets']
#        raw_region_length = subdata.iloc[idx, :].loc['raw_region_lengths']
#        plot_region_length = subdata.iloc[idx, :].loc['plot_region_lengths']
#        regional_offset = plot_region_length * ((pos0 - genome_intv.left) / raw_region_length)
#
#        return global_offset + regional_offset
#
#    def genomic_to_plot(self, chrom, pos0):
#        if chrom not in self.chromosomes:
#            return None
#
#        subdata = self.data.loc[chrom, :]
#        genome_intvlist = subdata.index.get_level_values('genome_interval')
#        contains = genome_intvlist.contains(pos0)
#
#        num_hit = contains.sum()
#        if num_hit == 0:
#            return None
#        elif num_hit > 1:
#            raise Exception(f'More than one intervals contains the input position.')
#
#        idx = np.where(contains)[0][0]
#        genome_intv = genome_intvlist[idx]
#        global_offset = subdata.iloc[idx, :].loc['global_offsets']
#        raw_region_length = subdata.iloc[idx, :].loc['raw_region_lengths']
#        plot_region_length = subdata.iloc[idx, :].loc['plot_region_lengths']
#        regional_offset = plot_region_length * ((pos0 - genome_intv.left) / raw_region_length)
#
#        return global_offset + regional_offset
#
#    def plot_to_genomic(self, x):
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
#
#    # Axes modification
#    def set_xlim(self, ax):
#        ax.set_xlim(*self.xlim)
#
#    def get_chrom_borders(self):
#        result = dict()
#        #for chrom in self.data.index.levels[0]:
#        for chrom in set(self.data.index.get_level_values('chromosome')):
#            subdf = self.data.xs(chrom, level='chromosome')
#            intvindex = subdf.index.get_level_values('plot_interval')
#            result[chrom] = (intvindex[0].left, intvindex[-1].right)
#        return result
#
#    def draw_chrom_borders(self, ax, color='black', linewidth=1, fontsize=8, **kwargs):
#        borders_set = set()
#        xlim = self.xlim
#        chrom_borders = self.get_chrom_borders()
#        # draw chromosome names
#        for chrom, (start0, end0) in chrom_borders.items():
#            ax.text(
#                (start0 + end0) / 2, 
#                ax.get_ylim()[1], 
#                chrom, 
#                ha='center', va='bottom',
#                size=fontsize,
#            )
#
#        # draw chromosome region borderlines
#        for chrom, (start0, end0) in chrom_borders.items():
#            if start0 != xlim[0]:
#                borders_set.add(start0)
#            if end0 != xlim[1]:
#                borders_set.add(end0)
#        for pos0 in borders_set:
#            ax.axvline(pos0, color=color, linewidth=linewidth, **kwargs)
#
#    def draw_hlines(self, ax, df, y_colname, offset=0, nproc=None, **kwargs):
#        joined_gr = pyranges_helper.join(self.totalregion_gr, pr.PyRanges(df), how='inner', merge='first', find_nearest=False, as_gr=True)
#        joined_df = joined_gr.df
#
#        ys = joined_df.loc[:, y_colname].array + offset
#
#        with multiprocessing.Pool(nproc) as pool:
#            xmins = pool.starmap(
#                self.genomic_to_plot, 
#                ((row['Chromosome'], row['Start']) for idx, row in joined_df.iterrows())
#            )
#            xmaxs = pool.starmap(
#                self.genomic_to_plot, 
#                ((row['Chromosome'], row['End']) for idx, row in joined_df.iterrows())
#            )
#
#        ax.hlines(ys, xmins, xmaxs, **kwargs)
#
#    def draw_dots(self, ax, df, y_colname, nproc=None, color='black', marker='o', linestyle='', **kwargs):
#        joined_gr = pyranges_helper.join(self.totalregion_gr, pr.PyRanges(df), how='inner', merge='first', find_nearest=False, as_gr=True)
#        joined_df = joined_gr.df
#
#        ys = joined_df.loc[:, y_colname].array
#        with multiprocessing.Pool(nproc) as pool:
#            xs = pool.starmap(
#                self.genomic_to_plot, 
#                ((row['Chromosome'], (row['Start'] + row['End']) / 2) for idx, row in joined_df.iterrows())
#            )
#
#        ax.plot(xs, ys, color=color, marker=marker, linestyle=linestyle, **kwargs)
#
#    def draw_bgcolors(self, ax, df, nproc=None, **kwargs):
#        joined_gr = pyranges_helper.join(self.totalregion_gr, pr.PyRanges(df), how='inner', merge='first', find_nearest=False, as_gr=True)
#        joined_df = joined_gr.df
#
#        with multiprocessing.Pool(nproc) as pool:
#            xmins = pool.starmap(
#                self.genomic_to_plot, 
#                ((row['Chromosome'], row['Start']) for idx, row in joined_df.iterrows())
#            )
#            xmaxs = pool.starmap(
#                self.genomic_to_plot, 
#                ((row['Chromosome'], row['End']) for idx, row in joined_df.iterrows())
#            )
#
#        ylims = ax.get_ylim()
#        ymin = ylims[0]
#        height = ylims[1] - ylims[0]
#        widths = np.array(xmaxs) - np.array(xmins)
#
#        boxes = [
#            Rectangle((xm, ymin), width=w, height=height)
#            for (xm, w) in zip(xmins, widths)
#        ]
#        ax.add_collection(PatchCollection(boxes, **kwargs))
#
#
## not used because slower than CoordConverter_old
#class CoordConverter2:
#    def __init__(self, gr):
#        assert isinstance(gr, pr.PyRanges), f'"gr" must be a pyranges.PyRanges object.'
#        #assert {'weight'}.issubset(set(gr.columns))
#        assert 'weight' in gr.columns, f'"gr" must include a column named "weight"'
#        self.gr = gr.copy()
#        self.set_data()
#
#    def set_data(self):
#        genome_intervals = pd.IntervalIndex.from_arrays(
#            left=list(self.gr.Start), right=list(self.gr.End), closed='left',
#        )
#        raw_region_lengths = self.gr.lengths().reset_index(drop=True)
#        plot_region_lengths = raw_region_lengths * self.gr.weight.reset_index(drop=True)
#        cumsum = plot_region_lengths.cumsum()
#        global_offsets = cumsum.shift(1, fill_value=0)
#        plot_intervals = pd.IntervalIndex.from_breaks(
#            ([0] + list(cumsum)), closed='left',
#        )
#
#        self.gr.raw_region_lengths = list(raw_region_lengths)
#        self.gr.plot_region_lengths = list(plot_region_lengths)
#        self.gr.global_offsets = list(global_offsets)
#
#        self.chromosomes = self.gr.Chromosome
#
#        self.genome_intervals = dict()
#        self.plot_intervals = dict()
#        for chrom in self.gr.keys():
#            self.genome_intervals[chrom] = genome_intervals[self.chromosomes == chrom]
#            self.plot_intervals[chrom] = plot_intervals[self.chromosomes == chrom]
#
#    def genomic_to_plot(self, chrom, pos0):
#        if chrom not in self.gr.keys():
#            return None
#
#        genome_intvlist = self.genome_intervals[chrom]
#        contains = genome_intvlist.contains(pos0)
#
#        num_hit = contains.sum()
#        if num_hit == 0:
#            return None
#        elif num_hit > 1:
#            raise Exception(f'More than one intervals contains the input position.')
#
#        idx = np.where(contains)[0][0]
#        genome_intv = genome_intvlist[idx]
#
#        subgr = self.gr[chrom]
#        global_offset = subgr.global_offsets.iloc[idx]
#        raw_region_length = subgr.raw_region_lengths.iloc[idx]
#        plot_region_length = subgr.plot_region_lengths.iloc[idx]
#        regional_offset = plot_region_length * ((pos0 - genome_intv.left) / raw_region_length)
#
#        return global_offset + regional_offset
#
#    def plot_to_genomic(self, x):
#        hit = False
#        for chrom, plot_intvlist in self.plot_intervals.items():
#            contains = plot_intvlist.contains(x)
#            if contains.any():
#                hit = True
#                break
#        if not hit:
#            return None
#
#        num_hit = contains.sum()
#        if num_hit == 0:
#            return None
#        elif num_hit > 1:
#            raise Exception(f'More than one intervals contains the input position.')
#
#        idx = np.where(contains)[0][0]
#        plot_intv = plot_intvlist[idx]
#        genome_intv = self.genome_intervals[chrom][idx]
#
#        regional_offset_fraction = (x - plot_intv.left) / plot_intv.length
#        pos0 = int(np.rint(genome_intv.left + (genome_intv.length * regional_offset_fraction)))
#
#        return chrom, pos0



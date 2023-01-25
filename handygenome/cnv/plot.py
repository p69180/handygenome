import itertools

import numpy as np
import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt
import seaborn as sns

import handygenome.common as common
import handygenome.cnv.mosdepth as libmosdepth


class CoordConverter:
    def __init__(self, gr):
        assert isinstance(gr, pr.PyRanges), f'"gr" must be a pyranges.PyRanges object.'
        #assert {'weight'}.issubset(set(gr.columns))
        assert 'weight' in gr.columns, f'"gr" must include a column named "weight"'
        self.raw_gr = gr
        self.chromosomes = set(self.raw_gr.keys())
        self.set_data()

    def set_data(self):
        genome_intervals = pd.IntervalIndex.from_arrays(
            left=list(self.raw_gr.Start), right=list(self.raw_gr.End), closed='left',
        )
        raw_region_lengths = self.raw_gr.lengths().reset_index(drop=True)
        plot_region_lengths = raw_region_lengths * self.raw_gr.weight.reset_index(drop=True)
        cumsum = plot_region_lengths.cumsum()
        global_offsets = cumsum.shift(1, fill_value=0)
        plot_intervals = pd.IntervalIndex.from_breaks(
            ([0] + list(cumsum)), closed='left',
        )

        self.data = pd.DataFrame(
            data={
                'raw_region_lengths': raw_region_lengths.array,
                'plot_region_lengths': plot_region_lengths.array,
                'global_offsets': global_offsets.array,
            },
            index=pd.MultiIndex.from_arrays(
                [self.raw_gr.Chromosome.array, genome_intervals, plot_intervals],
                names=('chromosome', 'genome_interval', 'plot_interval'),
            )
        )

        #self.genome_intervals = genome_intervals
        #self.plot_intervals = plot_intervals

    def genomic_to_plot(self, chrom, pos0):
        if chrom not in self.chromosomes:
            return None

        subdata = self.data.loc[chrom, :]
        genome_intvlist = subdata.index.get_level_values('genome_interval')
        contains = genome_intvlist.contains(pos0)

        num_hit = contains.sum()
        if num_hit == 0:
            return None
        elif num_hit > 1:
            raise Exception(f'More than one intervals contains the input position.')

        idx = np.where(contains)[0][0]
        genome_intv = genome_intvlist[idx]
        global_offset = subdata.iloc[idx, :].loc['global_offsets']
        raw_region_length = subdata.iloc[idx, :].loc['raw_region_lengths']
        plot_region_length = subdata.iloc[idx, :].loc['plot_region_lengths']
        regional_offset = plot_region_length * ((pos0 - genome_intv.left) / raw_region_length)

        return global_offset + regional_offset

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



# not used because slower than CoordConverter
class CoordConverter2:
    def __init__(self, gr):
        assert isinstance(gr, pr.PyRanges), f'"gr" must be a pyranges.PyRanges object.'
        #assert {'weight'}.issubset(set(gr.columns))
        assert 'weight' in gr.columns, f'"gr" must include a column named "weight"'
        self.gr = gr.copy()
        self.set_data()

    def set_data(self):
        genome_intervals = pd.IntervalIndex.from_arrays(
            left=list(self.gr.Start), right=list(self.gr.End), closed='left',
        )
        raw_region_lengths = self.gr.lengths().reset_index(drop=True)
        plot_region_lengths = raw_region_lengths * self.gr.weight.reset_index(drop=True)
        cumsum = plot_region_lengths.cumsum()
        global_offsets = cumsum.shift(1, fill_value=0)
        plot_intervals = pd.IntervalIndex.from_breaks(
            ([0] + list(cumsum)), closed='left',
        )

        self.gr.raw_region_lengths = list(raw_region_lengths)
        self.gr.plot_region_lengths = list(plot_region_lengths)
        self.gr.global_offsets = list(global_offsets)

        self.chromosomes = self.gr.Chromosome

        self.genome_intervals = dict()
        self.plot_intervals = dict()
        for chrom in self.gr.keys():
            self.genome_intervals[chrom] = genome_intervals[self.chromosomes == chrom]
            self.plot_intervals[chrom] = plot_intervals[self.chromosomes == chrom]

    def genomic_to_plot(self, chrom, pos0):
        if chrom not in self.gr.keys():
            return None

        genome_intvlist = self.genome_intervals[chrom]
        contains = genome_intvlist.contains(pos0)

        num_hit = contains.sum()
        if num_hit == 0:
            return None
        elif num_hit > 1:
            raise Exception(f'More than one intervals contains the input position.')

        idx = np.where(contains)[0][0]
        genome_intv = genome_intvlist[idx]

        subgr = self.gr[chrom]
        global_offset = subgr.global_offsets.iloc[idx]
        raw_region_length = subgr.raw_region_lengths.iloc[idx]
        plot_region_length = subgr.plot_region_lengths.iloc[idx]
        regional_offset = plot_region_length * ((pos0 - genome_intv.left) / raw_region_length)

        return global_offset + regional_offset

    def plot_to_genomic(self, x):
        hit = False
        for chrom, plot_intvlist in self.plot_intervals.items():
            contains = plot_intvlist.contains(x)
            if contains.any():
                hit = True
                break
        if not hit:
            return None

        num_hit = contains.sum()
        if num_hit == 0:
            return None
        elif num_hit > 1:
            raise Exception(f'More than one intervals contains the input position.')

        idx = np.where(contains)[0][0]
        plot_intv = plot_intvlist[idx]
        genome_intv = self.genome_intervals[chrom][idx]

        regional_offset_fraction = (x - plot_intv.left) / plot_intv.length
        pos0 = int(np.rint(genome_intv.left + (genome_intv.length * regional_offset_fraction)))

        return chrom, pos0



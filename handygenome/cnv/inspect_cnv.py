import gzip
import subprocess

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyranges as pr

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
cnvmisc = importlib.import_module('.'.join([top_package_name, 'cnv', 'misc']))
sequenza_handler = importlib.import_module('.'.join([top_package_name, 'cnv', 'sequenza_handler']))


def draw_CNVcall(extract_df, seg_df, tumor_depth_df, normal_depth_df, cellularity, ploidy, is_female, chromdict, sampleid, figsize=(25, 12), CNAB_offset=0.1):
    # basic params
    male_haploid_chroms = ('X', 'Y', 'chrX', 'chrY')
    male_only_chroms = ('Y', 'chrY')

    if is_female:
        irrelevant_chroms = male_only_chroms
        haploid_chroms = tuple()
    else:
        irrelevant_chroms = tuple()
        haploid_chroms = male_haploid_chroms

    # select assembled chromosomes
    from collections import OrderedDict
    chromdict_assembled = OrderedDict({key: val for key, val in chromdict.items() if common.RE_PATS['assembled_chromosome'].fullmatch(key)})
    contigs = list(chromdict_assembled.keys())
    lengths = list(chromdict_assembled.values())

    # make a copy of segments dataframe and modify it
    seg_df = seg_df.copy()
    add_expected_values(seg_df, cellularity, ploidy, irrelevant_chroms, haploid_chroms)
    add_ax_xmin_xmax(seg_df, chromdict_assembled)
    
    # init figure and axes
    fig = plt.figure(figsize=figsize)
    fig.suptitle(
        ', '.join([
            f'SampleID = {sampleid}',
            f'Cellularity = {cellularity}',
            f'Ploidy = {ploidy}',
            f'Gender = {"F" if is_female else "M"}',
        ])
    )
    axs = fig.subplots(
        6, len(chromdict_assembled),
        width_ratios=lengths,
        gridspec_kw={'wspace': 0, 'hspace': 0.05},
    )
    row_idx_normaldepth = -1
    row_idx_tumordepth = -2
    row_idx_depthr = -3
    row_idx_baf = -4
    row_idx_CNt = -5
    row_idx_CNAB = -6

    # main drawing processes
    basic_setups(axs, contigs, lengths)

    draw_depth(axs, row_idx_normaldepth, normal_depth_df, contigs, 'Normal depth')
    draw_depth(axs, row_idx_tumordepth, tumor_depth_df, contigs, 'Tumor depth')
    draw_depth_ratios(
        axs, row_idx_depthr, extract_df, cellularity, ploidy, contigs,
        irrelevant_chroms, haploid_chroms,
    )
    draw_baf(axs, row_idx_baf, extract_df, contigs)
    draw_CNt(axs, row_idx_CNt, seg_df, contigs, lengths)
    draw_CNAB(axs, row_idx_CNAB, seg_df, contigs, lengths, CNAB_offset=CNAB_offset)
    draw_expected_values(axs, row_idx_depthr, row_idx_baf, seg_df, contigs, alpha=0.5)

    return fig


def add_expected_values(seg_df, cellularity, ploidy, irrelevant_chroms, haploid_chroms):
    expected_depthrs = list()
    expected_bafs = list()
    for _, row in seg_df.iterrows():
        # depth ratio
        if row.Chromosome in irrelevant_chroms:
            depthr = np.nan
        elif row.Chromosome in haploid_chroms:
            depthr = cnvmisc.theoretical_depth_ratio(row.CNt, cellularity, ploidy, CNn=1)
        else:
            depthr = cnvmisc.theoretical_depth_ratio(row.CNt, cellularity, ploidy, CNn=2)
        expected_depthrs.append(depthr)
        # baf
        if np.isnan(row.B):
            baf = np.nan
        else:
            baf = cnvmisc.theoretical_baf(row.CNt, row.B, cellularity)
        expected_bafs.append(baf)

    seg_df.insert(seg_df.shape[1], 'expected_depthr', expected_depthrs)
    seg_df.insert(seg_df.shape[1], 'expected_baf', expected_bafs)


def add_ax_xmin_xmax(seg_df, chromdict_assembled):
    chromlens = seg_df.Chromosome.apply(lambda x: chromdict_assembled[x])
    seg_df.insert(seg_df.shape[1], 'xmin', seg_df.Start / chromlens)
    seg_df.insert(seg_df.shape[1], 'xmax', seg_df.End / chromlens)


def basic_setups(axs, contigs, lengths):
    # ticks and grids
    for idx, ax in np.ndenumerate(axs):
        chromlen = lengths[idx[1]]
        ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
        ax.set(xlim=(0, chromlen))
        ax.set_xticks(np.arange(0, chromlen, int(1e7)))
        ax.grid(axis='x', which='major', zorder=-1, alpha=0.4)

    # chromosome name labels
    for idx, ax in enumerate(axs[-1, :]):
        ax.set_xlabel(contigs[idx])


def draw_depth(axs, row_idx, depth_df, contigs, ylabel, alpha=0.2):
    ylims = (0, 2 * int(depth_df['Depth'].quantile(0.95)))
    yticks = np.arange(ylims[0], ylims[1], 10)

    axs[row_idx, 0].set_ylabel(ylabel)
    axs[row_idx, 0].set_yticks(yticks)
    axs[row_idx, 0].tick_params(labelleft=True, left=True)

    for idx, ax in enumerate(axs[row_idx, :]):
        ax.set(ylim=ylims)
        df_subset = depth_df.loc[depth_df.Chromosome == contigs[idx], :]
        if df_subset.shape[0] > 0:
            x = (df_subset.Start + df_subset.End) / 2
            y = df_subset['Depth']
            ax.scatter(x, y, s=2, c='black', alpha=alpha)


def draw_depth_ratios(
    axs, row_idx, extract_df, cellularity, ploidy, contigs,
    irrelevant_chroms, haploid_chroms,
    alpha_CNlines=0.6, alpha_dot=0.5,
):
    depthr_ylims = (0, np.ceil(extract_df.ratio_mean.quantile(0.95)) + 1)
    yticks = np.arange(depthr_ylims[0], depthr_ylims[1], 0.5)

    repr_depthrs = [
        cnvmisc.theoretical_depth_ratio(CNt, cellularity, ploidy, CNn=2)
        for CNt in range(7)
    ]
    repr_depthrs_haploid = [
        cnvmisc.theoretical_depth_ratio(CNt, cellularity, ploidy, CNn=1)
        for CNt in range(7)
    ]
    #male_haploid_chroms = ('X', 'Y', 'chrX', 'chrY')
    #male_only_chroms = ('Y', 'chrY')

    axs[row_idx, 0].set_ylabel('depth ratio')
    axs[row_idx, 0].set_yticks(yticks)
    axs[row_idx, 0].tick_params(labelleft=True, left=True)

    for idx, ax in enumerate(axs[row_idx, :]):
        # set ylim
        ax.set(ylim=depthr_ylims)
        
        chrom = contigs[idx]
        # draw representative depth ratio hlines
        if not (chrom in irrelevant_chroms):
            if chrom in haploid_chroms:
                repr_ys = repr_depthrs_haploid
            else:
                repr_ys = repr_depthrs

            for y in repr_ys:
                ax.axhline(y, linestyle='--', color='black', linewidth=1, zorder=-1, alpha=alpha_CNlines)
        # draw depth scatterplot
        df_subset = extract_df.loc[extract_df.Chromosome == chrom, :]
        if df_subset.shape[0] > 0:
            x = (df_subset.Start + df_subset.End) / 2
            y = df_subset.ratio_mean
            ax.scatter(x, y, s=2, c='black', alpha=alpha_dot)
            
            
def draw_baf(axs, row_idx, extract_df, contigs, alpha_tickgrid=0.4, alpha_dot=0.5):
    ylims = (0, 0.6)
    baf_yticks = np.arange(ylims[0], ylims[1], 0.1)

    axs[row_idx, 0].set_ylabel('BAF')
    axs[row_idx, 0].set_yticks(baf_yticks)
    axs[row_idx, 0].tick_params(labelleft=True, left=True)

    for idx, ax in enumerate(axs[row_idx, :]):
        ax.set(ylim=ylims)
        # ax.set(ylim=(0, 1.5))
        for y in baf_yticks:
            ax.axhline(y, linestyle='-', color='black', linewidth=1, zorder=-1, alpha=alpha_tickgrid)

        df_subset = extract_df.loc[extract_df.Chromosome == contigs[idx], :]
        if df_subset.shape[0] > 0:
            x = (df_subset.Start + df_subset.End) / 2
            y = df_subset.baf_mean
            ax.scatter(x, y, s=2, c='black', alpha=alpha_dot)
            
            
def draw_CNt(axs, row_idx, seg_df, contigs, lengths, alpha_tickgrid=0.4):
    CNt_ylims = (0, int(seg_df.CNt.max()))
    CNt_yticks = np.arange(CNt_ylims[0], CNt_ylims[1], 2)

    axs[row_idx, 0].set_ylabel('CNt')
    axs[row_idx, 0].set_yticks(CNt_yticks)
    axs[row_idx, 0].tick_params(labelleft=True, left=True)

    for idx, ax in enumerate(axs[row_idx, :]):
        ax.set(ylim=CNt_ylims)
        # ax.set(ylim=(0, 1.5))
        for y in CNt_yticks:
            ax.axhline(y, linestyle='-', color='black', linewidth=1, zorder=-1, alpha=alpha_tickgrid)

        df_subset = seg_df.loc[seg_df.Chromosome == contigs[idx], :]
        if df_subset.shape[0] > 0:
            for _, row in df_subset.iterrows():
                ax.axhline(row['CNt'], xmin=row['xmin'], xmax=row['xmax'], color='black')
                
                
def draw_CNAB(axs, row_idx, seg_df, contigs, lengths, CNAB_offset=0.1, alpha_tickgrid=0.4):
    CNAB_ylims = (0, int(max(seg_df.A.max(), seg_df.B.max())))
    CNAB_yticks = np.arange(CNAB_ylims[0], CNAB_ylims[1], 2)

    axs[row_idx, 0].set_ylabel('CNA, CNB')
    axs[row_idx, 0].set_yticks(CNAB_yticks)
    axs[row_idx, 0].tick_params(labelleft=True, left=True)

    for idx, ax in enumerate(axs[row_idx, :]):
        ax.set(ylim=CNAB_ylims)
        # ax.set(ylim=(0, 1.5))
        for y in CNAB_yticks:
            ax.axhline(y, linestyle='-', color='black', linewidth=1, zorder=-1, alpha=alpha_tickgrid)

        df_subset = seg_df.loc[seg_df.Chromosome == contigs[idx], :]
        if df_subset.shape[0] > 0:
            for _, row in df_subset.iterrows():
                ax.axhline(row['A'] + (2 * CNAB_offset), xmin=row['xmin'], xmax=row['xmax'], color='red')
                ax.axhline(row['B'] + CNAB_offset, xmin=row['xmin'], xmax=row['xmax'], color='blue')


def draw_expected_values(axs, row_idx_depthr, row_idx_baf, seg_df, contigs, alpha=0.5):
    for axs_col_idx in range(axs.shape[1]):
        ax_depthr = axs[row_idx_depthr, axs_col_idx]
        ax_baf = axs[row_idx_baf, axs_col_idx]
        df_subset = seg_df.loc[seg_df.Chromosome == contigs[axs_col_idx], :]
        if df_subset.shape[0] > 0:
            for _, row in df_subset.iterrows():
                ax_depthr.axhline(row['expected_depthr'], xmin=row['xmin'], xmax=row['xmax'], color='green', alpha=alpha)
                if not np.isnan(row['expected_baf']):
                    ax_baf.axhline(row['expected_baf'], xmin=row['xmin'], xmax=row['xmax'], color='green', alpha=alpha)



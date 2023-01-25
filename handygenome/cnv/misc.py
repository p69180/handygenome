import collections
import random
import itertools
import array

import Bio.SeqUtils
import numpy as np
import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
#from scipy.stats.mstats import winsorize

import handygenome.common as common
import handygenome.cnv.mosdepth as libmosdepth


MALE_HAPLOID_CHROMS = ('X', 'Y', 'chrX', 'chrY')

CCFResult = collections.namedtuple('CCFResult', ('CNm', 'ccf', 'mutated_allele'))


def theoretical_depth_ratio(CNt, cellularity, ploidy, CNn=2, normal_ploidy=2, avg_depth_ratio=1):
    cellularity_term = ((CNt / CNn) * cellularity) + (1 - cellularity)
    ploidy_term = ((ploidy / normal_ploidy) * cellularity) + (1 - cellularity)
    return (avg_depth_ratio * cellularity_term) / ploidy_term


def theoretical_baf(CNt, CN_B, cellularity, CNn=2):
    if CNn <= 1:
        return None
    else:
        return (
            ((CN_B * cellularity) + (1 - cellularity)) /
            ((CNt * cellularity) + CNn * (1 - cellularity))
        )


def get_mean_ploidy(seg_df):
    seg_lengths = seg_df.End - seg_df.Start
    numerator = (seg_df.CNt * seg_lengths).sum()
    denominator = seg_lengths.sum()
    return numerator / denominator


def theoretical_somatic_vaf(CNt, CNm, cellularity, ccf, CNn=2):
    """Assumes no overlapping subclonal CNV.
    Args:
        CNm: Mutated copy number
        CNn: Germline copy number (should be 1 in haploid area such as male sex chromosome)
        ccf: Cancer cell fraction
    """
    numerator = CNm * cellularity * ccf
    denominator = CNt * cellularity + CNn * (1 - cellularity)
    return numerator / denominator


def ccf_from_params(vaf, CNt, CNm, cellularity, CNn):
    """Reverse of 'theoretical_somatic_vaf'"""
    numerator = vaf * (CNt * cellularity + CNn * (1 - cellularity))
    denominator = CNm * cellularity
    return numerator / denominator


def CNm_likelihood_difference(observed_vaf, CNt, CNm, cellularity, CNn):
    expected_vaf = theoretical_somatic_vaf(
        CNt=CNt, CNm=CNm, cellularity=cellularity, ccf=1, CNn=CNn,
    )
    return -abs(expected_vaf - observed_vaf)


def CNm_likelihood_binom(total_read, var_read, CNt, CNm, cellularity, CNn):
    expected_vaf = theoretical_somatic_vaf(
        CNt=CNt, CNm=CNm, cellularity=cellularity, ccf=1, CNn=CNn,
    )
    return scipy.stats.binom.pmf(var_read, total_read, expected_vaf)


def get_ccf_CNm(vaf, cellularity, CNt, CNB=None, likelihood='diff', total_read=None, var_read=None, CNn=2):
    """Notes:
        - Assumes no overlapping subclonal CNV
        - If subclonal, mutation may reside in only one copy. CNm must be 1.
        - If clonal, CNm can be any value between 1 and max(CNA, CNB).
    Args:
        CNn: Should be 1 in haploid area such as male sex chromosome
        CNB: None when haploid
    """
    assert likelihood in ('diff', 'binom')

    # NA cases
    if (
        CNt == 0 or
        np.isnan(vaf)
    ):
        return CCFResult(None, None, None)

    # get CNA
    if CNB is None:
        CNA = None
    else:
        CNA = CNt - CNB
    # max_subclonal_vaf results when CNm == 1 and ccf == 1
    max_subclonal_vaf = theoretical_somatic_vaf(
        CNt=CNt, CNm=1, cellularity=cellularity, ccf=1, CNn=CNn,
    )

    if vaf >= max_subclonal_vaf:
        # Assumes clonal mutation
        ccf = 1
        max_CNm = CNt if (CNB is None) else max(CNA, CNB)
        CNm_candidates = list(range(1, max_CNm + 1))
        if likelihood == 'diff':
            CNm = max(
                CNm_candidates, 
                key=(lambda x: CNm_likelihood_difference(
                    observed_vaf=vaf, CNt=CNt, CNm=x, cellularity=cellularity, CNn=CNn,
                )),
            )
        elif likelihood == 'binom':
            CNm = max(
                CNm_candidates, 
                key=(lambda x: CNm_likelihood_binom(
                    total_read=total_read, 
                    var_read=var_read, 
                    CNt=CNt, 
                    CNm=x, 
                    cellularity=cellularity, 
                    CNn=CNn,
                )),
            )
    else:
        # Assumes subclonal mutation
        CNm = 1
        ccf = ccf_from_params(vaf=vaf, CNt=CNt, CNm=CNm, cellularity=cellularity, CNn=CNn)

    # get which allele is mutated
    if CNB is None:
        mutated_allele = None
    else:
        if CNA == CNB:
            mutated_allele = None
        else:
            if CNm > min(CNA, CNB):
                if CNA > CNB:
                    mutated_allele = 'A'
                else:
                    mutated_allele = 'B'
            else:
                mutated_allele = None

    return CCFResult(CNm, ccf, mutated_allele)


def check_haploid(is_female, chrom):
    return (not is_female) and (chrom in MALE_HAPLOID_CHROMS)



###############################################

# gc & depth processsing

def add_gc_to_gr(gr, fasta, window=None):
    """Args:
        gr: pyranges.PyRanges

    Changes in-place
    """
    gr.GC = make_gc(gr.Chromosome, gr.Start, gr.End, fasta, window=window)
        

def add_gc_to_df(df, fasta, window=None):
    """Args:
        df: pandas.DataFrame

    Changes in-place
    """
    df.insert(
        df.shape[1], 'GC', make_gc(df['Chromosome'], df['Start'], df['End'], fasta, window=window)
    )
        

def make_gc(chroms, start0s, end0s, fasta, window=None):
    if window is None:
        def get_fetchargs(chrom, start0, end0):
            return chrom, start0, end0
    else:
        chromlens = dict(zip(fasta.references, fasta.lengths))

        pad_left = int(window / 2)
        if window % 2 == 0:
            pad_right = pad_left
        else:
            pad_right = pad_left + 1

        def get_fetchargs(chrom, start0, end0):
            if end0 - start0 >= window:
                return chrom, start0, end0
            else:
                mid = int((start0 + end0) / 2)
                new_start0 = max(0, mid - pad_left)
                new_end0 = min(chromlens[chrom], mid + pad_right)
                return chrom, new_start0, new_end0

    return array.array(
        'f',
        (
            Bio.SeqUtils.gc_fraction(
                fasta.fetch(*get_fetchargs(chrom, start0, end0))
            )
            for chrom, start0, end0 in zip(chroms, start0s, end0s)
        )
    )


def get_gc_breaks(n_bin):
    breaks = np.linspace(0, 1, (n_bin + 1), endpoint=True)
    breaks[0] = -0.001
    return breaks


# not used
def get_gc_depth_data(depth_df, n_bin):
    assert 'GC' in depth_df.columns

    breaks = get_gc_breaks(n_bin)
    cutresult = pd.cut(depth_df['GC'], bins=breaks)

    gcbin_average_depths = dict()
    for intv, subdf in depth_df.groupby(cutresult, axis=0):
        # since "cutresult" is a categorical series, iteration over groupby object yields all bins including empty ones
        if subdf.shape[0] == 0:
            gcbin_average_depths[intv] = np.nan
        else:
            gcbin_average_depths[intv] = np.average(subdf['mean_depth'], weights=(subdf['End'] - subdf['Start']))
                # weighted average depth for the given GC bin

    gcbin_average_depths = pd.Series(gcbin_average_depths)
    gcbin_norm_average_depths = gcbin_average_depths / gcbin_average_depths.mean(skipna=True)

    return gcbin_average_depths, gcbin_norm_average_depths, cutresult


def get_gcbin_data(depth_df, gc_breaks):
    """This does not apply weights when calculating averages. Suitable for even-sized bins"""
    # cut by gc bins
    cutresult = pd.cut(depth_df['GC'], bins=gc_breaks)
    # make results
    gcbin_depth_data = {intv: list() for intv in cutresult.dtype.categories}
    for intv, depth in zip(cutresult, depth_df['mean_depth']):
        gcbin_depth_data[intv].append(depth)

    gcbin_average_depths = {intv: np.mean(data) for intv, data in gcbin_depth_data.items()}
    gcbin_average_depths = pd.Series(gcbin_average_depths)
    gcbin_norm_average_depths = gcbin_average_depths / gcbin_average_depths.mean(skipna=True)

    return gcbin_depth_data, gcbin_average_depths, gcbin_norm_average_depths, cutresult


def postprocess_depth_df(
    depth_df, fasta, 

    nan_lower_cutoff=None,
    nan_upper_cutoff=None,

    gcdata_trim_limits=(0.05, 0.05),
    gcdata_lower_cutoff=None,
    gcdata_upper_cutoff=None,

    n_gcbin=100, 
    as_gr=True, 
    gc_window=None,
):
    """Adds 4 columns:
        - GC
        - norm_mean_depth
        - gc_corrected_mean_depth
        - sequenza_style_norm_mean_depth
    """
    #depth_df = depth_df.copy()

    # replace outliers with np.nan
    if nan_lower_cutoff is not None:
        depth_df.loc[depth_df['mean_depth'] < nan_lower_cutoff, ['mean_depth']] = np.nan
    if nan_upper_cutoff is not None:
        depth_df.loc[depth_df['mean_depth'] >= nan_upper_cutoff, ['mean_depth']] = np.nan

    # GC
    add_gc_to_df(depth_df, fasta, window=gc_window)

    # norm_mean_depth
    depth_df_wo_nan = depth_df.loc[~np.isnan(depth_df['mean_depth']), :]
    global_mean_depth = np.average(depth_df_wo_nan['mean_depth'], weights=(depth_df_wo_nan['End'] - depth_df_wo_nan['Start']))
    depth_df.insert(depth_df.shape[1], 'norm_mean_depth', (depth_df['mean_depth'] / global_mean_depth).array)

    # get gc depth data
    gc_breaks = get_gc_breaks(n_gcbin)
    trimmed_depth_df = remove_outliers(depth_df, trim_limits=gcdata_trim_limits, lower_cutoff=gcdata_lower_cutoff, upper_cutoff=gcdata_upper_cutoff)
    gcbin_depth_data, gcbin_average_depths, gcbin_norm_average_depths, cutresult = get_gcbin_data(trimmed_depth_df, gc_breaks)
    #gcbin_average_depths, gcbin_norm_average_depths, cutresult = get_gc_depth_data(depth_df, n_bin=n_gcbin)

    # gc_corrected_mean_depth
    cutresult = pd.cut(depth_df['GC'], bins=gc_breaks)
    gcbin_norm_average_depths_selected = gcbin_norm_average_depths.loc[cutresult]
    depth_df.insert(depth_df.shape[1], 'gc_corrected_mean_depth', (depth_df['mean_depth'].array / gcbin_norm_average_depths_selected.array))

    # sequenza_style_norm_mean_depth
    gcbin_average_depths_selected = gcbin_average_depths.loc[cutresult]
    depth_df.insert(depth_df.shape[1], 'sequenza_style_norm_mean_depth', (depth_df['mean_depth'].array / gcbin_average_depths_selected.array))

    if as_gr:
        return pr.PyRanges(depth_df, int64=True)
    else:
        return depth_df


def arghandler_ylims(ydata, ylims):
    if ylims is None:
        return None
    else:
        return (
            min(ylims[0], min(ydata)), 
            max(ylims[1], max(ydata)),
        )


def remove_outliers(depth_df, trim_limits=(0.05, 0.05), lower_cutoff=None, upper_cutoff=None):
    # remove nan
    #selector = np.repeat(True, depth_df.shape[0])
    selector = ~np.isnan(depth_df['mean_depth'])

    # trim by absolute cutoff
    if lower_cutoff is not None:
        selector = np.logical_and(selector, depth_df['mean_depth'] >= lower_cutoff)
    if upper_cutoff is not None:
        selector = np.logical_and(selector, depth_df['mean_depth'] <= upper_cutoff)

    depth_df = depth_df.loc[selector, :]

    # trim by proportion
    trimmed_depth = scipy.stats.mstats.trim(
        depth_df['mean_depth'], limits=trim_limits, relative=True,
    )
    depth_df = depth_df.loc[~trimmed_depth.mask, :]

    return depth_df


def plot_gc_depth(
    fasta, bam_path, 
    gc_window=50, 
    xwidth=1000, 
    chrom=None,
    start0=None,
    end0=None,
    return_gr=True, 
    depth_ylims=None, 
    gc_ylims=None,
    alpha=0.6,
    linewidth=1,
):
    pad_left = int(gc_window / 2)
    if gc_window % 2 == 0:
        pad_right = pad_left
    else:
        pad_right = pad_left + 1    

    if any(x is None for x in [chrom, start0, end0]):
        chrom = random.choice(fasta.references[:24])
        chromlen = fasta.lengths[fasta.references.index(chrom)]
        start0 = random.choice(range(pad_left, (chromlen - xwidth - pad_right) + 1))
        end0 = start0 + xwidth

    fetch_start0 = start0 - pad_left
    fetch_end0 = end0 + pad_right
    
    ### GC ###
    seq = fasta.fetch(chrom, fetch_start0, fetch_end0)
    x = list()
    y_gc = list()
    for idx in range(0, end0 - start0):
        subseq = seq[idx:(idx + gc_window)]
        gcfrac = Bio.SeqUtils.gc_fraction(subseq)
        x.append(start0 + idx)
        y_gc.append(gcfrac)

    # modify ylims
    gc_ylims = arghandler_ylims(y_gc, gc_ylims)
    
    # plot
    fig, ax_gc = plt.subplots(figsize=(20, 5))
    if gc_ylims is not None:
        ax_gc.set_ylim(*gc_ylims)
    ax_gc.set_title(f'{chrom}:{start0 + 1:,}-{end0:,}')
    ax_gc.plot(
        x, y_gc, 
        color='green', marker=None, linestyle='-', markersize=1, linewidth=linewidth, alpha=alpha
    )
    
    ### depth ###
    region_gr = pr.from_dict({'Chromosome': [chrom], 'Start': [start0], 'End': [end0]})
    gr, gr_perbase = libmosdepth.run_mosdepth(
        bam_path, 
        use_median=False,
        region_gr=region_gr, 
        as_gr=True, 
        load_perbase=True,
    )
    y_depth = list(gr_perbase.depth)

    # modify ylims
    depth_ylims = arghandler_ylims(y_depth, depth_ylims)

    # plot
    ax_depth = ax_gc.twinx()
    if depth_ylims is not None:
        ax_depth.set_ylim(*depth_ylims)
    ax_depth.plot(
        x, y_depth, 
        color='black', marker=None, linestyle='-', markersize=1, linewidth=linewidth, alpha=alpha
    )

    # return
    if return_gr:
        gr_perbase.GC = y_gc
        return gr_perbase


def plot_gc_distribution(
    depth_df, 

    n_gcbin=100,

    trim_limits=(0.05, 0.05),
    lower_cutoff=None,
    upper_cutoff=None,

    histogram_mode='2d',
    scatter_ylims=None,
    heatmap_ylims=None, 
    heatmap_vmax=None,
):
    def make_average_scatter(gcbin_average_depths, ylims=None):
        x = list()
        y = list()
        for key, val in gcbin_average_depths.items():
            x.append(key.mid)
            y.append(np.mean(val))

        ylims = arghandler_ylims(y, ylims)

        fig, ax = plt.subplots(figsize=(10, 5))
        if ylims is not None:
            ax.set_ylim(*ylims)
        ax.set_xlabel('GC fraction')
        ax.set_ylabel('mean depth')
        ax.plot(x, y, 'ko', markersize=2)
        return fig

    def make_heatmap(gcbin_depth_data, ylims=None, vmax=None):
        ydata = list(itertools.chain.from_iterable(gcbin_depth_data.values()))
        if ylims is None:
            ylims = (min(ydata), max(ydata)) 
        else:
            ylims = arghandler_ylims(ydata, ylims)

        # set breaks
        breaks = np.linspace(ylims[0], ylims[1], 100, endpoint=True)

        # get histogram
        data = dict()
        for gc_interval, val in gcbin_depth_data.items():
            data[gc_interval] = np.histogram(val, bins=breaks, density=True)[0]
        df_plot = pd.DataFrame.from_dict(data)
        df_plot.index = pd.IntervalIndex.from_breaks(np.round(breaks, 1))

        # draw
        fig, ax = plt.subplots(figsize=(15, 10))
        sns.heatmap(df_plot, ax=ax, vmax=vmax)
        ax.set_xlabel('GC fraction')
        ax.set_ylabel('mean depth')
        return fig

    def make_heatmap_hist2d(depth_df, n_gcbin, ylims=None, vmax=None):
        # set breaks
        if ylims is None:
            ylims = (min(depth_df['mean_depth']), max(depth_df['mean_depth'])) 
        else:
            ylims = arghandler_ylims(depth_df['mean_depth'], ylims)
        depth_breaks = np.linspace(ylims[0], ylims[1], 100, endpoint=True)
        #gc_breaks = np.linspace(0, 1, 101, endpoint=True)
        gc_breaks = get_gc_breaks(n_gcbin)

        # get histogram
        H, xedges, yedges = np.histogram2d(
            x=depth_df['mean_depth'], y=depth_df['GC'], bins=(depth_breaks, gc_breaks), density=True,
        )
        plot_df = pd.DataFrame(
            H, 
            index=pd.IntervalIndex.from_breaks(np.round(xedges, 1)), 
            columns=pd.IntervalIndex.from_breaks(np.round(yedges, 2)),
        )
      
        # draw
        fig, ax = plt.subplots(figsize=(15, 10))
        sns.heatmap(plot_df, ax=ax, vmax=vmax)
        ax.set_xlabel('GC fraction')
        ax.set_ylabel('mean depth')
        return fig

    # main
    assert histogram_mode in ('1d', '2d')

    depth_df = remove_outliers(depth_df, trim_limits=trim_limits, lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff)
    gc_breaks = get_gc_breaks(n_gcbin)
    gcbin_depth_data, gcbin_average_depths, gcbin_norm_average_depths, cutresult = get_gcbin_data(depth_df, gc_breaks)
    fig_avg_scatter = make_average_scatter(gcbin_average_depths, ylims=scatter_ylims)
    if histogram_mode == '2d':
        fig_heatmap = make_heatmap_hist2d(depth_df, n_gcbin, ylims=heatmap_ylims, vmax=heatmap_vmax)
    elif histogram_mode == '1d':
        fig_heatmap = make_heatmap(gcbin_depth_data, ylims=heatmap_ylims, vmax=heatmap_vmax)

    return fig_avg_scatter, fig_heatmap



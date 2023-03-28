import os
import collections
import random
import itertools
import functools
import multiprocessing
import operator
import warnings
import inspect

import Bio.SeqUtils
import numpy as np
import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
#from scipy.stats.mstats import winsorize

import handygenome.common as common
import handygenome.pyranges_helper as pyranges_helper
import handygenome.cnv.mosdepth as libmosdepth
import handygenome.cnv.gcfraction as gcfraction
import handygenome.assemblyspec as libassemblyspec


MALE_HAPLOID_CHROMS = ('X', 'Y', 'chrX', 'chrY')
DEFAULT_CNT_WEIGHT = 5


class CCFInfo(
    collections.namedtuple('CCFInfo', ('CNm', 'ccf', 'mutated_allele'))
):
    pass


class CPPair(
    collections.namedtuple('CPPair', ('cellularity', 'ploidy'))
):
    pass


####################
# argument handler #
####################

def check_having_coord_cols(arg):
    if not (
        {'Chromosome', 'Start', 'End'}.issubset(arg.columns)
        or arg.index.names == ['Chromosome', 'Start', 'End']
    ):
        raise Exception(f'Input dataframe must have columns "Chromosome", "Start", and "End".')


def arg_into_df(arg):
    if arg is None:
        return None
    elif isinstance(arg, pd.DataFrame):
        check_having_coord_cols(arg)
        return arg
    elif isinstance(arg, pr.PyRanges):
        return arg.df
    else:
        raise Exception(f'Argument must be either None, pd.DataFrame, or pr.PyRanges')


def arg_into_gr(arg):
    if arg is None:
        return None
    elif isinstance(arg, pd.DataFrame):
        check_having_coord_cols(arg)
        return pr.PyRanges(arg)
    elif isinstance(arg, pr.PyRanges):
        return arg
    else:
        raise Exception(f'Argument must be either None, pd.DataFrame, or pr.PyRanges')


################
# peak finding #
################

def get_1d_peaks(row, invert=False):
    #assert not hasattr(row, '__next__'), f'"row" must not be an iterator'

    row = np.array(row)
    if invert:
        row = -1 * row

    peaks_result, _ = scipy.signal.find_peaks(row)
    if len(peaks_result) == 0:  # uniformly increasing
        peak_indexes = list()
        if len(row) > 1:
            if row[0] > row[1]:
                peak_indexes.append(0)
            if row[-1] > row[-2]:
                peak_indexes.append(len(row) - 1)
        elif len(row) == 1:
            peak_indexes.append(0)
            
#        peak_indexes = sorted(
#            (
#                x[0] for x in
#                common.multi_max(enumerate(row), key=operator.itemgetter(1))
#            )
#        )
    else:
        diff = np.diff(row)

        peak_indexes = list()
        for within_row_idx in peaks_result:
            peak_indexes.append(within_row_idx)
            # leftward
            current_idx = within_row_idx
            while True:
                current_idx -= 1
                if diff[current_idx] != 0:
                    break
                else:
                    peak_indexes.append(current_idx)

                if current_idx == 0:
                    break
            # rightward
            current_idx = within_row_idx - 1
            while True:
                current_idx += 1
                if diff[current_idx] != 0:
                    break
                else:
                    peak_indexes.append(current_idx + 1)

                if current_idx == len(diff) - 1:
                    break

        peak_indexes.sort()

    return peak_indexes


def get_df_diagonals(df, upslope=True):
    """Sorting of diagonal rows: top-down
    Sorting of values in a row: left-right
    """
    def get_start_coord_upslope(n, nrow):
        if n < nrow:
            return (n, 0)
        else:
            return (nrow - 1, n - (nrow - 1))
        
    def get_diags_upslope(df, N, nrow, ncol):
        result = list()
        for n in range(N):
            pairs = list() 
            coord = get_start_coord_upslope(n, nrow)
            while True:
                if (coord[0] < 0) or (coord[1] >= ncol):
                    break

                val = df.iloc[coord[0], coord[1]]
                pairs.append((coord, val))
                coord = (coord[0] - 1, coord[1] + 1)
            result.append(pairs)
        return result

    def get_start_coord_downslope(n, ncol):
        if n < ncol:
            return (0, ncol - 1 - n)
        else:
            return (n - (ncol - 1), 0)
        
    def get_diags_downslope(df, N, nrow, ncol):
        result = list()
        for n in range(N):
            pairs = list() 
            coord = get_start_coord_downslope(n, ncol)
            while True:
                if (coord[0] >= nrow) or (coord[1] >= ncol):
                    break

                val = df.iloc[coord[0], coord[1]]
                pairs.append((coord, val))
                coord = (coord[0] + 1, coord[1] + 1)
            result.append(pairs)
        return result

    N = sum(df.shape) - 1
    nrow = df.shape[0]
    ncol = df.shape[1]

    if upslope:
        return get_diags_upslope(df, N, nrow, ncol)
    else:
        return get_diags_downslope(df, N, nrow, ncol)


def get_peak_coords_from_diagonals(diagonals, invert=False):
    peak_coords = list()
    for diag in diagonals:
        peak_indexes = get_1d_peaks([x[1] for x in diag], invert=invert)
        peak_coords.extend(diag[idx][0] for idx in peak_indexes)
    return peak_coords


def get_df_diagonal_peaks(df, invert=False):
    upslope_diags = get_df_diagonals(df, upslope=True)
    upslope_peak_coords = get_peak_coords_from_diagonals(upslope_diags, invert=invert)

    downslope_diags = get_df_diagonals(df, upslope=False)
    downslope_peak_coords = get_peak_coords_from_diagonals(downslope_diags, invert=invert)
        
    return upslope_peak_coords, downslope_peak_coords


def find_df_peaks(df, invert=False):
    row_peaks = df.apply(lambda x: get_1d_peaks(x, invert=invert), axis=1)
    row_peak_indexes = list()
    for row_idx, peaks in enumerate(row_peaks):
        row_peak_indexes.extend((row_idx, col_idx) for col_idx in peaks)

    col_peaks = df.apply(lambda x: get_1d_peaks(x, invert=invert), axis=0)
    col_peak_indexes = list()
    for col_idx, peaks in enumerate(col_peaks):
        col_peak_indexes.extend((row_idx, col_idx) for row_idx in peaks)

    return list(set(row_peak_indexes).intersection(set(col_peak_indexes)))


def find_df_peaks_4directions(df, invert=False):
    row_peaks = df.apply(lambda x: get_1d_peaks(x, invert=invert), axis=1)
    row_peak_indexes = list()
    for row_idx, peaks in enumerate(row_peaks):
        row_peak_indexes.extend((row_idx, col_idx) for col_idx in peaks)

    col_peaks = df.apply(lambda x: get_1d_peaks(x, invert=invert), axis=0)
    col_peak_indexes = list()
    for col_idx, peaks in enumerate(col_peaks):
        col_peak_indexes.extend((row_idx, col_idx) for row_idx in peaks)

    upslope_peak_coords, downslope_peak_coords = get_df_diagonal_peaks(df, invert=invert)

    return set.intersection(
        set(row_peak_indexes), set(col_peak_indexes),
        set(upslope_peak_coords), set(downslope_peak_coords),
    )


def find_df_peak_cpvalues(df, invert=False):
    assert df.index.name == 'cellularity'
    assert df.columns.name == 'ploidy'

    #peaks = find_df_peaks(df, invert=invert)
    peaks = find_df_peaks_4directions(df, invert=invert)
    cpvalues = [CPPair(df.index[x], df.columns[y]) for (x, y) in peaks]
    return sorted(cpvalues)


def get_peak_info(cp_score_dict, key='segfit', invert=True):
    dfs = make_cpscore_dfs(cp_score_dict)
    peaks_cpvalues = find_df_peak_cpvalues(dfs[key], invert=invert)
    peak_values = list()
    for c, p in peaks_cpvalues:
        data = {
            k: v for k, v in cp_score_dict[(c, p)].items()
            if k != 'CNt_list'
        }
        data['cellularity'] = c
        data['ploidy'] = p
        peak_values.append(data)
    #peak_values = sorted(peak_values, key=(lambda x: x['segfit_score']))

    return peak_values, dfs


###############
# CNn related #
###############

@functools.cache
def get_CNn_gr(refver, is_female):
    """Only accepts PAR-available reference versions"""
    chromdict = common.DEFAULT_CHROMDICTS[refver]
    autosomal_chroms = [
        x for x in chromdict.contigs 
        if (
            (common.RE_PATS['assembled_chromosome'].fullmatch(x) is not None)
            and (x not in MALE_HAPLOID_CHROMS)
        )
    ]  # 1, 2, ...,
    X_chrom, Y_chrom = chromdict.XY_names
    par_gr = libassemblyspec.get_par_gr(refver)

    # autosomal
    all_chrom_gr = chromdict.to_gr()
    autosomal_gr = all_chrom_gr[all_chrom_gr.Chromosome.isin(autosomal_chroms)]
    autosomal_gr.CNn = 2
    # sex
    if is_female:
        sex_gr = all_chrom_gr[all_chrom_gr.Chromosome.isin([X_chrom, Y_chrom])]
        sex_gr.CNn = 2
    else:
        all_sex_gr = all_chrom_gr[all_chrom_gr.Chromosome.isin([X_chrom, Y_chrom])]
        nonpar_sex_gr = all_sex_gr.subtract(par_gr)
        nonpar_sex_gr.CNn = 1
        par_sex_gr = par_gr.copy()
        par_sex_gr.CNn = [(2 if x == X_chrom else 0) for x in par_sex_gr.Chromosome]
        sex_gr = pr.concat([nonpar_sex_gr, par_sex_gr]).sort()
    # result
    return pr.concat([autosomal_gr, sex_gr]).sort()


def add_CNn_to_gr(gr, refver, is_female):
    return pyranges_helper.join(
        gr, get_CNn_gr(refver, is_female), how='left', merge='first',
        as_gr=True,
    )


def get_normal_mean_ploidy(refver, is_female, target_region_gr=None):
    """Only supports hg19 or hg38"""
    if is_female:
        return 2
    else:
        # sanity check
        if target_region_gr is not None:
            assert 'CNn' not in target_region_gr.columns

        chromdict = common.DEFAULT_CHROMDICTS[refver]
        if target_region_gr is None:
            N_region_gr = libassemblyspec.get_N_region_gr(refver)
            target_region_gr = chromdict.to_gr().subtract(N_region_gr)

        CNn_gr = get_CNn_gr(refver, is_female)
        target_region_gr = pyranges_helper.join(
            target_region_gr, CNn_gr, how='left', merge='longest', as_gr=True,
        )

        indexer = target_region_gr.CNn.notna()
        CNn_list = target_region_gr.CNn.array[indexer]
        weights = target_region_gr.lengths().array[indexer]
        return np.average(CNn_list, weights=weights)


################
# miscellanous #
################

def check_haploid(is_female, chrom):
    return (not is_female) and (chrom in MALE_HAPLOID_CHROMS)


def get_quantile_df(df):
    return df.applymap(lambda x: (df < x).to_numpy().mean())


def get_bafs(vaf_list):
    bafs = np.array(vaf_list)
    gt_half_indexes = bafs > 0.5
    bafs[gt_half_indexes] = 1 - bafs[gt_half_indexes]
    return bafs


#################################
# theoretical value calculation # 
#################################

def theoretical_depth_ratio(
    CNt, cellularity, tumor_mean_ploidy, CNn, normal_mean_ploidy, 
    tumor_avg_depth_ratio=1, normal_avg_depth_ratio=1,
):
    """Args:
        CNn, normal_ploidy must be nonzero
    Can be vectorized
    """
    tumor_upper_term = (CNt * cellularity) + (CNn * (1 - cellularity))
    tumor_lower_term = (tumor_mean_ploidy * cellularity) + (normal_mean_ploidy * (1 - cellularity))
    tumor_norm_depth = tumor_avg_depth_ratio * (tumor_upper_term / tumor_lower_term)

    #normal_upper_term = CNn
    #normal_lower_term = normal_mean_ploidy
    #normal_norm_depth = normal_avg_depth_ratio * (normal_upper_term / normal_lower_term)
    normal_norm_depth = normal_avg_depth_ratio * (CNn / normal_mean_ploidy)

    return tumor_norm_depth / normal_norm_depth


def inverse_theoretical_depth_ratio(
    depth_ratio, cellularity, tumor_mean_ploidy, CNn, normal_mean_ploidy, 
    tumor_avg_depth_ratio=1, normal_avg_depth_ratio=1,
):
    """Args:
        cellularity must be nonzero
    Returns:
        CNt estimate

    """
    '''
    <Induction>
    tumor_upper_term = (CNt * cellularity) + (CNn * (1 - cellularity))
    tumor_lower_term = (tumor_mean_ploidy * cellularity) + (normal_mean_ploidy * (1 - cellularity))
    tumor_norm_depth = tumor_avg_depth_ratio * (tumor_upper_term / tumor_lower_term)
    normal_norm_depth = normal_avg_depth_ratio * (CNn / normal_mean_ploidy)
    depth_ratio = tumor_norm_depth / normal_norm_depth
    ---------------------------------
    tumor_norm_depth = depth_ratio * normal_norm_depth
    tumor_avg_depth_ratio * (tumor_upper_term / tumor_lower_term) = depth_ratio * normal_norm_depth
    tumor_upper_term = (depth_ratio * normal_norm_depth * tumor_lower_term) / tumor_avg_depth_ratio
    (CNt * cellularity) + (CNn * (1 - cellularity)) = (depth_ratio * normal_norm_depth * tumor_lower_term) / tumor_avg_depth_ratio
    CNt * cellularity = ((depth_ratio * normal_norm_depth * tumor_lower_term) / tumor_avg_depth_ratio) - (CNn * (1 - cellularity))
    CNt = (
        ((depth_ratio * normal_norm_depth * tumor_lower_term) / tumor_avg_depth_ratio) - (CNn * (1 - cellularity))
    ) / cellularity
    '''
    normal_norm_depth = normal_avg_depth_ratio * (CNn / normal_mean_ploidy)
    tumor_lower_term = (tumor_mean_ploidy * cellularity) + (normal_mean_ploidy * (1 - cellularity))
    CNt = (
        ((depth_ratio * normal_norm_depth * tumor_lower_term) / tumor_avg_depth_ratio) - (CNn * (1 - cellularity))
    ) / cellularity
    return CNt


@functools.cache
def delta_depth_ratio(
    cellularity, tumor_mean_ploidy, CNn, normal_mean_ploidy, 
    tumor_avg_depth_ratio=1, normal_avg_depth_ratio=1,
):
    """Returns:
        Difference of depth ratio caused by 1 CNt change
    """
    '''
    <Induction>
    tumor_upper_term = (CNt * cellularity) + (CNn * (1 - cellularity))
    tumor_lower_term = (tumor_mean_ploidy * cellularity) + (normal_mean_ploidy * (1 - cellularity))
    tumor_norm_depth = tumor_avg_depth_ratio * (tumor_upper_term / tumor_lower_term)
    normal_norm_depth = normal_avg_depth_ratio * (CNn / normal_mean_ploidy)
    depth_ratio = tumor_norm_depth / normal_norm_depth
    ------------------------------------------
    depth_ratio = (tumor_avg_depth_ratio * (tumor_upper_term / tumor_lower_term)) / normal_norm_depth
    depth_ratio = tumor_upper_term * (tumor_avg_depth_ratio / (tumor_lower_term * normal_norm_depth))
    depth_ratio = (
        ((CNt * cellularity) + (CNn * (1 - cellularity))) 
        * 
        (tumor_avg_depth_ratio / (tumor_lower_term * normal_norm_depth))
    )
    delta_depth_ratio = (
        ((CNt_1 * cellularity) - (CNt_2 * cellularity)) 
        * 
        (tumor_avg_depth_ratio / (tumor_lower_term * normal_norm_depth))
    )
    delta_depth_ratio = (
        cellularity
        * 
        (tumor_avg_depth_ratio / (tumor_lower_term * normal_norm_depth))
    )
    delta_depth_ratio = (cellularity * tumor_avg_depth_ratio) / (tumor_lower_term * normal_norm_depth)
    '''
    tumor_lower_term = (tumor_mean_ploidy * cellularity) + (normal_mean_ploidy * (1 - cellularity))
    normal_norm_depth = normal_avg_depth_ratio * (CNn / normal_mean_ploidy)
    delta_depth_ratio = (cellularity * tumor_avg_depth_ratio) / (tumor_lower_term * normal_norm_depth)
    return delta_depth_ratio


def theoretical_depth_ratio_sequenza(CNt, cellularity, ploidy, CNn=2, normal_ploidy=2, avg_depth_ratio=1):
    """Args:
        CNn, normal_ploidy must be nonzero
    """
    cellularity_term = ((CNt / CNn) * cellularity) + (1 - cellularity)
    ploidy_term = ((ploidy / normal_ploidy) * cellularity) + (1 - cellularity)
    return avg_depth_ratio * (cellularity_term / ploidy_term)


def inverse_theoretical_depth_ratio_sequenza(depth_ratio, cellularity, ploidy, CNn=2, normal_ploidy=2, avg_depth_ratio=1):
    """Args:
        cellularity must be nonzero
    Returns:
        CNt estimate
    """
    ploidy_term = ((ploidy / normal_ploidy) * cellularity) + (1 - cellularity)
    # depth_ratio * ploidy_term == avg_depth_ratio * cellularity_term
    return (((depth_ratio * ploidy_term) / avg_depth_ratio) - (1 - cellularity)) * (CNn / cellularity)


@functools.cache
def delta_depth_ratio_sequenza(cellularity, ploidy, CNn=2, normal_ploidy=2, avg_depth_ratio=1):
    ploidy_term = ((ploidy / normal_ploidy) * cellularity) + (1 - cellularity)
    return (avg_depth_ratio / ploidy_term) * (cellularity / CNn)


def theoretical_baf(CNt, B, cellularity, CNn):
    if CNn <= 1:
        return np.nan
    else:
        numerator = (B * cellularity) + (1 - cellularity)
        denominator = (CNt * cellularity) + (CNn * (1 - cellularity))
        return numerator / denominator


def theoretical_baf_vectorized(CNt, B, cellularity, CNn):
    numerator = (B * cellularity) + (1 - cellularity)
    denominator = (CNt * cellularity) + (CNn * (1 - cellularity))
    notna_values = numerator / denominator

    result = pd.Series(np.repeat(np.nan, len(CNn)))
    result.where(CNn <= 1, notna_values, inplace=True)
    return result


def inverse_theoretical_baf(baf, cellularity, CNt, CNn):
    """Args:
        cellularity must be nonzero
    Returns:
        B estimate
    """
    '''
    numerator = (B * cellularity) + (1 - cellularity)
    denominator = (CNt * cellularity) + (CNn * (1 - cellularity))
    baf = numerator / denominator

    baf = ((B * cellularity) + (1 - cellularity)) / denominator
    (B * cellularity) + (1 - cellularity) = baf * denominator
    B * cellularity = (baf * denominator) - (1 - cellularity)
    B = ((baf * denominator) - (1 - cellularity)) / cellularity
    '''
    denominator = (CNt * cellularity) + (CNn * (1 - cellularity))
    B = ((baf * denominator) - (1 - cellularity)) / cellularity
    return B


#############################################
# find an optimal copy number from cp value # 
#############################################


def get_B(CNt, cellularity, baf, CNn):
    B_estimate = inverse_theoretical_baf(baf, cellularity, CNt, CNn)
    if B_estimate <= 0:
        B = 0
        theo_baf = theoretical_baf(CNt, B, cellularity, CNn)
        if np.isnan(theo_baf):
            B = np.nan
            diff = np.nan
        else:
            diff = abs(baf - theo_baf)
    else:
        B_candidate_upper = int(np.ceil(B_estimate))
        B_candidate_lower = int(np.floor(B_estimate))

        theo_baf_upper = theoretical_baf(CNt, B_candidate_upper, cellularity, CNn)
        theo_baf_lower = theoretical_baf(CNt, B_candidate_lower, cellularity, CNn)

        if np.isnan(theo_baf_upper) and np.isnan(theo_baf_lower):
            B = np.nan
            diff = np.nan
        elif (not np.isnan(theo_baf_upper)) and np.isnan(theo_baf_lower):
            B = B_candidate_upper
            diff = abs(baf - theo_baf_upper)
        elif np.isnan(theo_baf_upper) and (not np.isnan(theo_baf_lower)):
            B = B_candidate_lower
            diff = abs(baf - theo_baf_lower)
        else:
            diff_upper = abs(theo_baf_upper - baf)
            diff_lower = abs(theo_baf_lower - baf)
            if diff_upper <= diff_lower:
                B = B_candidate_upper
                diff = diff_upper
            else:
                B = B_candidate_lower
                diff = diff_lower

    return B, diff


def get_CN_from_cp(cellularity, ploidy, depth_ratio, baf, CNt_weight, CNn, normal_ploidy):
    def save_cache(CNt_candidate, B_cache, segfit_score_cache, baf_diff_cache):
        B, baf_diff = get_B(CNt_candidate, cellularity, baf, CNn)
        if np.isnan(B):
            # drop this CNt candidate if B cannot be calculated
            return

        ratio_diff = abs(
            theoretical_depth_ratio(CNt_candidate, cellularity, ploidy, CNn, normal_ploidy) 
            - depth_ratio
        )
        B_cache[CNt_candidate] = B
        baf_diff_cache[CNt_candidate] = baf_diff
        segfit_score_cache[CNt_candidate] = ratio_diff * CNt_weight + baf_diff

    # get initial CNt candidates
    CNt_estimate = inverse_theoretical_depth_ratio(
        depth_ratio, cellularity, ploidy, CNn, normal_ploidy,
    )
    if CNt_estimate <= 0:
        initial_candidate_upper = 0
        initial_candidate_lower = None
    else:
        initial_candidate_upper = int(np.ceil(CNt_estimate))
        initial_candidate_lower = int(np.floor(CNt_estimate))

    # set caches
    B_cache = dict()
    segfit_score_cache = dict()
    baf_diff_cache = dict()

    delta_ratio = delta_depth_ratio(cellularity, ploidy, CNn, normal_ploidy) * CNt_weight

    # upper
    save_cache(initial_candidate_upper, B_cache, segfit_score_cache, baf_diff_cache)
    if initial_candidate_upper in baf_diff_cache:
        n_further_candidates = int(baf_diff_cache[initial_candidate_upper] / delta_ratio)
        for offset in range(n_further_candidates):
            CNt_candidate = initial_candidate_upper + (offset + 1)
            save_cache(CNt_candidate, B_cache, segfit_score_cache, baf_diff_cache)
    # lower
    if initial_candidate_lower is not None:
        save_cache(initial_candidate_lower, B_cache, segfit_score_cache, baf_diff_cache)
        if initial_candidate_lower in baf_diff_cache:
            try:
                n_further_candidates = int(baf_diff_cache[initial_candidate_lower] / delta_ratio)
            except:
                print(cellularity, ploidy, depth_ratio, baf)
                print(baf, baf_diff_cache[initial_candidate_lower], initial_candidate_lower, delta_ratio)
                raise

            for offset in range(n_further_candidates):
                CNt_candidate = initial_candidate_lower - (offset + 1)
                if CNt_candidate < 0:
                    break
                save_cache(CNt_candidate, B_cache, segfit_score_cache, baf_diff_cache)

    # result
    if len(segfit_score_cache) == 0:
        return np.nan, np.nan, np.nan
    else:
        CNt, segfit_score = min(segfit_score_cache.items(), key=operator.itemgetter(1))
        B = B_cache[CNt]
        #selected_score = segfit_score_cache[CNt]

        return CNt, B, segfit_score


def get_CN_from_cp_wo_baf(cellularity, ploidy, depth_ratio, CNt_weight, CNn, normal_ploidy):
    CNt_estimate = inverse_theoretical_depth_ratio(
        depth_ratio, cellularity, ploidy, CNn, normal_ploidy,
    )
    if CNt_estimate <= 0:
        CNt = 0
        diff = abs(
            theoretical_depth_ratio(CNt, cellularity, ploidy, CNn, normal_ploidy) 
            - depth_ratio
        )
        segfit_score = diff * CNt_weight
    else:
        upper_candidate = int(np.ceil(CNt_estimate))
        lower_candidate = int(np.floor(CNt_estimate))

        upper_diff = abs(
            theoretical_depth_ratio(upper_candidate, cellularity, ploidy, CNn, normal_ploidy) 
            - depth_ratio
        )
        lower_diff = abs(
            theoretical_depth_ratio(lower_candidate, cellularity, ploidy, CNn, normal_ploidy) 
            - depth_ratio
        )
        if upper_diff <= lower_diff:
            CNt = upper_candidate
            segfit_score = upper_diff * CNt_weight
        else:
            CNt = lower_candidate
            segfit_score = lower_diff * CNt_weight

    return CNt, np.nan, segfit_score 


#########################
# find optimal cp value #
#########################


def calc_cp_score(
    segment_df, cellularity, ploidy, is_female, CNt_weight, normal_ploidy,
):
    def applied_func(row):
        if check_haploid(is_female, row.iloc[0]) or np.isnan(row.iloc[4]):
            CNt, B, segfit_score = get_CN_from_cp_wo_baf(cellularity, ploidy, row.iloc[3], CNt_weight, row.iloc[5], normal_ploidy)
        else:
            CNt, B, segfit_score = get_CN_from_cp(cellularity, ploidy, row.iloc[3], row.iloc[4], CNt_weight, row.iloc[5], normal_ploidy)

        return CNt, B, segfit_score

    assert {'depthratio_segment_mean', 'baf_segment_mean', 'CNn'}.issubset(segment_df.columns)
    assert segment_df['depthratio_segment_mean'].notna().all()

    # set column order
    segment_df = segment_df.loc[
        :, 
        ['Chromosome', 'Start', 'End', 'depthratio_segment_mean', 'baf_segment_mean', 'CNn']
    ]

    apply_result = segment_df.apply(applied_func, axis=1)
    segfit_score_sum = sum(x[2] for x in apply_result if not np.isnan(x[2]))
    CNt_list = np.fromiter((x[0] for x in apply_result), dtype=np.float_)
    B_list = np.fromiter((x[1] for x in apply_result), dtype=np.float_)

    return {'segfit_score': segfit_score_sum, 'CNt_list': CNt_list, 'B_list': B_list}


def get_cp_score_dict(
    segment_df, refver, is_female, target_region_gr, 
    CNt_weight=DEFAULT_CNT_WEIGHT, nproc=None,
):
    def get_calculated_tumor_ploidies_new(segment_df, scorelist, target_region_gr):
        stripped_segment_df = segment_df.loc[:, ['Chromosome', 'Start', 'End']].copy()
        all_indexes = list(range(stripped_segment_df.shape[0]))  # means index of segments
        stripped_segment_df.insert(3, 'index', all_indexes)
        index_annotated_targetregion_gr = annotate_region_with_segment(  # each region is annotated with corresponding segment index
            target_region_gr[['Chromosome', 'Start', 'End']], 
            pr.PyRanges(stripped_segment_df),
        )
        # make weights of CNt values
        index_annotated_targetregion_gr.length = index_annotated_targetregion_gr.lengths()
        weights_dict = index_annotated_targetregion_gr.df.loc[:, ['length', 'index']].groupby('index').sum().to_dict()['length']
        weights = [
            (weights_dict[x] if x in weights_dict else 0)
            for x in all_indexes   
        ]

        CNt_values = np.array([x['CNt_list'] for x in scorelist])
        calculated_tumor_ploidies = np.average(CNt_values, axis=1, weights=weights)

        return calculated_tumor_ploidies

    def get_calculated_tumor_ploidies(segment_df, scorelist, target_region_gr):
        stripped_segment_df = segment_df.loc[:, ['Chromosome', 'Start', 'End']]
        CNt_col_list = [
            pd.Series(x['CNt_list'], name=f'CNt_{idx}') 
            for idx, x in enumerate(scorelist)
        ]
        CNt_annotated_segment_df = pd.concat(
            ([stripped_segment_df] + CNt_col_list), axis=1
        )
        CNt_annotated_segment_gr = pr.PyRanges(CNt_annotated_segment_df)
        CNt_annotated_targetregion_gr = annotate_region_with_segment(
            target_region_gr[['Chromosome', 'Start', 'End']], CNt_annotated_segment_gr
        )
        calculated_tumor_ploidies = np.average(
            CNt_annotated_targetregion_gr.df.iloc[:, 3:],
            axis=0,
            weights=CNt_annotated_targetregion_gr.lengths(),
        )
        return calculated_tumor_ploidies

    # get calculated CNt and segment fit cores
    normal_ploidy = get_normal_mean_ploidy(refver, is_female, target_region_gr)
    c_candidates = np.round(np.arange(0.01, 1.00, 0.01), 2)  # 0 and 1 are excluded
    p_candidates = np.round(np.arange(0.5, 7.1, 0.1), 1)  # 0.5 to 7.0
    cp_pairs = tuple(CPPair(x, y) for (x, y) in itertools.product(c_candidates, p_candidates))

    with multiprocessing.Pool(nproc) as pool:
        scorelist = pool.starmap(
            calc_cp_score, 
            (   
                (segment_df, x.cellularity, x.ploidy, is_female, CNt_weight, normal_ploidy) 
                for x in cp_pairs
            ),
        )

    # get calculated tumor mean ploidy for each cp pair
    calculated_tumor_ploidies = get_calculated_tumor_ploidies_new(segment_df, scorelist, target_region_gr)

    # calculate ploidy fitting scores
    for cpppair, calc_ploidy, dic in zip(cp_pairs, calculated_tumor_ploidies, scorelist):
        dic['calculated_tumor_ploidy'] = calc_ploidy
        dic['ploidy_diff'] = calc_ploidy - cpppair.ploidy
        dic['ploidy_diff_ratio'] = dic['ploidy_diff'] / cpppair.ploidy

    cp_scores = dict(zip(cp_pairs, scorelist))
    return cp_scores


def make_cpscore_dfs(cpscore_dict):
    clist = sorted(set(x.cellularity for x in cpscore_dict.keys()))
    plist = sorted(set(x.ploidy for x in cpscore_dict.keys()))
    data = {
        'segfit': list(), 
        'ploidy_diff': list(), 
        'ploidy_diff_abs': list(),
        'ploidy_diff_ratio': list(),
        'ploidy_diff_ratio_abs': list(),
    }
    for c in clist:
        values = [
            val for key, val in 
            sorted(
                ((k, v) for (k, v) in cpscore_dict.items() if k.cellularity == c),
                key=(lambda y: y[0].ploidy),
            )
        ]
        data['segfit'].append([x['segfit_score'] for x in values])
        data['ploidy_diff'].append([x['ploidy_diff'] for x in values])
        data['ploidy_diff_abs'].append([abs(x['ploidy_diff']) for x in values])
        data['ploidy_diff_ratio'].append([x['ploidy_diff_ratio'] for x in values])
        data['ploidy_diff_ratio_abs'].append([abs(x['ploidy_diff_ratio']) for x in values])

    dfs = dict()
    for key, val in data.items():
        df = pd.DataFrame.from_records(val, index=clist, columns=plist)
        df.index.name = 'cellularity'
        df.columns.name = 'ploidy'
        dfs[key] = df

    return dfs


def show_heatmap_peaks(score_df, invert=True, quantile=True, figsize=None):
    #peaks = find_df_peaks(score_df, invert=invert)
    peaks = find_df_peaks_4directions(score_df, invert=invert)
    ys, xs = zip(*peaks)

    fig, ax = plt.subplots(figsize=figsize)
    if quantile:
        ax = sns.heatmap(get_quantile_df(score_df), ax=ax)
    else:
        ax = sns.heatmap(score_df, ax=ax)
    ax.plot(xs, ys, linestyle='', marker='o', markersize=1.5)

    return fig, ax


def show_heatmap_peaks_new(dfs, invert=True, quantile=True, figsize=None, fourway=True):
    if fourway:
        peaks = find_df_peaks_4directions(dfs['segfit'], invert=invert)
    else:
        peaks = find_df_peaks(dfs['segfit'], invert=invert)

    ys, xs = zip(*peaks)
    ploidy_ratios = [dfs['ploidy_diff_ratio'].iloc[y, x] for y, x in zip(ys, xs)]
    segscores = [dfs['segfit'].iloc[y, x] for y, x in zip(ys, xs)]

    fig, ax = plt.subplots(figsize=figsize)
    if quantile:
        ax = sns.heatmap(get_quantile_df(dfs['segfit']), ax=ax)
    else:
        ax = sns.heatmap(dfs['segfit'], ax=ax)
    ax.plot(xs, ys, linestyle='', marker='o', markersize=3)
    for x, y, pratio, segscore in zip(xs, ys, ploidy_ratios, segscores):
        cellularity = dfs['segfit'].index[y]
        ploidy = dfs['segfit'].columns[x]
        text = '\n'.join([
            f'cellularity: {cellularity}',
            f'ploidy: {ploidy}',
            f'segfit score: {round(segscore, 6)}',
            f'ploidy diff ratio: {round(pratio, 6)}',
        ])
        ax.text(
            x + 0.5, y - 0.5, text, 
            fontsize=6,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
        )

    return fig, ax


##################################
# segment dataframe modification #
##################################


def annotate_region_with_segment(target_region_gr, segment_gr, as_gr=True):
    """Add values of segment_gr, such as CNt or B, to target_region_gr"""
    return pyranges_helper.join(
        target_region_gr, segment_gr, how='left', merge='first', 
        find_nearest=True, as_gr=as_gr,
    )


def add_CNt_to_segment(
    segment_df, cellularity, tumor_ploidy, normal_ploidy, is_female, 
    CNt_weight=DEFAULT_CNT_WEIGHT,
):
    assert isinstance(segment_df, pd.DataFrame)

    data = calc_cp_score(
        segment_df, cellularity, tumor_ploidy, is_female, CNt_weight, normal_ploidy,
    )
    result = segment_df.copy()
    result['CNt'] = data['CNt_list']
    result['B'] = data['B_list']
    return result


def add_theoreticals_to_segment(
    segment_df, cellularity, tumor_ploidy, normal_ploidy, is_female,
):
    assert isinstance(segment_df, pd.DataFrame)
    assert {'CNn', 'CNt', 'B'}.issubset(segment_df.columns)

    def depthr_getter(row):
        return theoretical_depth_ratio(
            row['CNt'], cellularity, tumor_ploidy, row['CNn'], normal_ploidy, 
            tumor_avg_depth_ratio=1, normal_avg_depth_ratio=1,
        )

    def baf_getter(row):
        if np.isnan(row['B']):
            return np.nan
        else:
            return theoretical_baf(row['CNt'], row['B'], cellularity, row['CNn'])

    theo_depthr_list = segment_df.apply(depthr_getter, axis=1)
    theo_baf_list = segment_df.apply(baf_getter, axis=1)

    return segment_df.assign(
        **{'depthratio_predicted': theo_depthr_list, 'baf_predicted': theo_baf_list}
    )
    

###################################


def _get_mean_ploidy_from_values(start0_list, end0_list, CNt_list):
    lengths = pd.Series(end0_list) - pd.Series(start0_list)
    return np.average(CNt_list, weights=lengths)


def get_mean_ploidy(seg_df):
    seg_lengths = seg_df.End - seg_df.Start
    numerator = (seg_df.CNt * seg_lengths).sum()
    denominator = seg_lengths.sum()
    return numerator / denominator


def get_tumor_mean_ploidy(CNt_annotated_gr):
    return np.average(CNt_annotated_gr.CNt, weights=CNt_annotated_gr.lengths())


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
        return CCFInfo(None, None, None)

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

    return CCFInfo(CNm, ccf, mutated_allele)


#################################
# depth df processing functions #
#################################

def get_gc_breaks(n_bin):
    breaks = np.linspace(0, 1, (n_bin + 1), endpoint=True)
    breaks[0] = -0.001
    intervals = pd.IntervalIndex.from_breaks(breaks, closed='right')

    return breaks, intervals


def handle_outliers(depth_df, trim_limits=None, lower_cutoff=None, upper_cutoff=None):
    """Args:
        trim_limits: Tuple with length 2 (lower proportion, upper proportion).
            Used as an argument to "scipy.stats.mstats.trim".
            Example: (0.05, 0.05)
    Returns:
        Tuple (depth_df, selector)
            depth_df: The input dataframe, with outlier depth values replaced with np.nan
            selector: boolean array indicating non-outlier rows
        Outlier handling processes:
            1) Removes/masks oringinal nan values
            2) Removes/masks based on relative proportions
            3) Removes/masks based on absolute cutoff values
    """
    assert 'mean_depth' in depth_df.columns

    # remove nan
    selector = depth_df['mean_depth'].notna()

    # trim by proportion
    if trim_limits is not None:
        trimmed_depth = scipy.stats.mstats.trim(
            depth_df['mean_depth'], limits=trim_limits, relative=True,
        )
        selector = np.logical_and(selector, ~trimmed_depth.mask)

    # trim by absolute cutoff
    if lower_cutoff is not None:
        selector = np.logical_and(selector, depth_df['mean_depth'] >= lower_cutoff)
    if upper_cutoff is not None:
        selector = np.logical_and(selector, depth_df['mean_depth'] <= upper_cutoff)

    # result
    selector = selector.to_numpy()
    depth_df.loc[~selector, 'mean_depth'] = np.nan

    return depth_df, selector


def get_depth_df_from_file(chroms, start0s, end0s, depth_file_path):
    preset_depth_df = pd.read_csv(
        depth_file_path, 
        sep='\t', 
        names=('Chromosome', 'Start', 'End', 'depth'), 
        dtype={'Chromosome': str, 'Start': int, 'End': int, 'depth': float},
    )
    preset_depth_gr = pr.PyRanges(preset_depth_df)
    result_gr = pr.PyRanges(chromosomes=chroms, starts=start0s, ends=end0s)
    result_gr = pyranges_helper.join(result_gr, preset_depth_gr, how='left', merge='weighted_mean', as_gr=True)
    return result_gr


# NOT USED
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


def get_gcbin_data(depth_df, gc_breaks, make_raw_data=False):
    """This does not apply weights when calculating averages. Suitable for even-sized bins"""
    # cut by gc bins
    cutresult = pd.cut(depth_df['GC'], bins=gc_breaks, include_lowest=False)
    # make results
    if make_raw_data:
        gcbin_depth_data = {intv: list() for intv in cutresult.dtype.categories}
        for intv, depth in zip(cutresult, depth_df['mean_depth']):
            gcbin_depth_data[intv].append(depth)

        gcbin_average_depths = {
            intv: (np.nan if len(data) == 0 else np.mean(data))
            for intv, data in gcbin_depth_data.items()
        }
    else:
        gcbin_depth_data = None
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)

            gcbin_average_depths = dict()
            for intv in cutresult.dtype.categories:
                gcbin_average_depths[intv] = np.nanmean(
                    depth_df['mean_depth'].loc[(cutresult == intv).array]
                )

    gcbin_average_depths = pd.Series(gcbin_average_depths)
    gcbin_norm_average_depths = gcbin_average_depths / gcbin_average_depths.mean(skipna=True)

    return gcbin_depth_data, gcbin_average_depths, gcbin_norm_average_depths, cutresult


def get_normal_wgs_depth_cutoffs(depth_df):
    mean = depth_df['mean_depth'].mean()


@common.get_deco_num_set_differently(('fasta', 'refver', 'gc_df'), 1)
def postprocess_depth_df(
    depth_df, 
    *,
    refver=None,
    fasta=None, 
    gc_df=None,

    gc_window=None,

    preset_cutoffs=None,

    lower_cutoff=None,
    upper_cutoff=None,
    trim_limits=None,

    #nan_lower_cutoff=None,
    #nan_upper_cutoff=None,

    #gcdata_trim_limits=(0.05, 0.05),
    #gcdata_lower_cutoff=None,
    #gcdata_upper_cutoff=None,

    n_gcbin=100, 
    as_gr=True, 

    verbose=False,
):
    """Args:
        fasta, refver: Only used for gc fraction calculation. Only one of
            "fasta" or "refver" must be set. Using "refver" is recommended.

    Returns:
        postprocessed df & gcbin average depth dict
        Adds 4 columns:
            - GC
            - norm_mean_depth
            - gc_corrected_mean_depth
            - sequenza_style_norm_mean_depth
    """
    if verbose:
        def printlog(msg):
            funcname = inspect.stack()[0].function
            common.print_timestamp(f'{funcname}: {msg}')
    else:
        def printlog(msg):
            pass

    # sanity check
    depth_df = arg_into_df(depth_df)
    assert 'mean_depth' in depth_df.columns, f'"depth_df" must include a column named "mean_depth"'
    assert preset_cutoffs in ('wgs', 'normal_wgs', 'panel', None)

    # set outlier cutoffs
    if preset_cutoffs == 'wgs':
        printlog(f'Running "postprocess_depth_df" function with preset cutoff mode "wgs"')
        lower_cutoff = 0
        upper_cutoff = 2000
    elif preset_cutoffs == 'panel':
        printlog(f'Running "postprocess_depth_df" function with preset cutoff mode "panel"')
        lower_cutoff = 50
        upper_cutoff = np.inf

    # add GC fractions
    printlog(f'Adding GC fraction values')
    if gc_df is None:
        gcfraction.add_gc_calculating(depth_df, refver=refver, fasta=fasta, window=gc_window)
    else:
        gc_df = arg_into_df(gc_df)
        assert 'GC' in gc_df.columns
        assert gc_df.index.names == ['Chromosome', 'Start', 'End']
            # gc_df must have a MultiIndex ['Chromosome', 'Start', 'End']

        depth_df = depth_df.join(gc_df, on=['Chromosome', 'Start', 'End'], how='left')

    # replace outliers with nan
    depth_df, selector = handle_outliers(
        depth_df, trim_limits=trim_limits, lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff
    )

    # norm_mean_depth
    printlog(f'Getting normalized depths')
    global_mean_depth = common.nanaverage(
        depth_df['mean_depth'].to_numpy(), 
        weights=(depth_df['End'] - depth_df['Start']).to_numpy(),
    )
    depth_df['norm_mean_depth'] = (depth_df['mean_depth'] / global_mean_depth).array

    # get gc depth data
    printlog(f'Getting depth data by gc value ranges')
    gc_breaks, gc_intervals = get_gc_breaks(n_gcbin)
    (
        gcbin_depth_data, 
        gcbin_average_depths, 
        gcbin_norm_average_depths, 
        cutresult 
    ) = get_gcbin_data(depth_df, gc_breaks, make_raw_data=False)
    cutresult_idx = common.get_indexes_of_array(cutresult, gcbin_norm_average_depths.index)

    # gc_corrected_mean_depth
    printlog(f'Getting GC-corrected depths')
    gcbin_norm_average_depths_selected = gcbin_norm_average_depths.iloc[cutresult_idx]
    depth_df['gc_corrected_mean_depth'] = (
        depth_df['mean_depth'].array / gcbin_norm_average_depths_selected.array
    )

    # sequenza_style_norm_mean_depth
    printlog(f'Getting GC-corrected normalized depths (sequenza style)')
    gcbin_average_depths_selected = gcbin_average_depths.iloc[cutresult_idx]
    depth_df['sequenza_style_norm_mean_depth'] = (
        depth_df['mean_depth'].array / gcbin_average_depths_selected.array
    )

    # result
    if as_gr:
        depth_df = pr.PyRanges(depth_df, int64=False)

    printlog(f'Finished')

    return depth_df, gcbin_average_depths


def get_processed_depth_df(
    bam_path, fasta, region_gr, 
    gc_window=None, 
    preset_cutoffs='wgs',
    lower_cutoff=None,
    upper_cutoff=None,
    trim_limits=None,
    donot_subset_bam=False, 
    as_gr=False,
):
    # run mosdepth
    depth_df, _ = libmosdepth.run_mosdepth(
        bam_path,
        t=8, 
        region_gr=region_gr, 
        donot_subset_bam=donot_subset_bam, 
        as_gr=False, 
        load_perbase=False,
    )

    # postprocess
    depth_gr, gcbin_average_depths = postprocess_depth_df(
        depth_df, 

        fasta=fasta, 
        gc_window=gc_window,

        preset_cutoffs=preset_cutoffs,
        lower_cutoff=lower_cutoff,
        upper_cutoff=upper_cutoff,
        trim_limits=trim_limits,

        n_gcbin=100,
        as_gr=as_gr, 
    )

    return depth_gr, gcbin_average_depths


def make_depth_ratio_df(tumor_depth_df, normal_depth_df, as_gr=False):
    def assign_ratio(result, new_colname, numerator, denominator):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)

            arr = numerator / denominator
            arr[np.isinf(arr)] = np.nan
            result[new_colname] = arr

    # sanity check
    tumor_depth_df = arg_into_df(tumor_depth_df)
    normal_depth_df = arg_into_df(normal_depth_df)
    required_cols = {'mean_depth', 'GC', 'norm_mean_depth', 'gc_corrected_mean_depth', 'sequenza_style_norm_mean_depth'}
    assert all(
        required_cols.issubset(df.columns)
        for df in (tumor_depth_df, normal_depth_df)
    )

    compared_columns = ['Chromosome', 'Start', 'End', 'GC']
    tumor_common_cols = tumor_depth_df.loc[:, compared_columns]
    normal_common_cols = normal_depth_df.loc[:, compared_columns]
    if not (normal_common_cols == tumor_common_cols).all(axis=None):
        raise Exception(
            f'The columns {compared_columns} of the two input dataframes must be identical.'
        )

    # main
    result = tumor_common_cols
    assign_ratio(
        result, 
        'depth_ratio_sequenzastyle', 
        tumor_depth_df['sequenza_style_norm_mean_depth'].values, 
        normal_depth_df['sequenza_style_norm_mean_depth'].values,
    )
    assign_ratio(
        result, 
        'depth_ratio_mystyle', 
        tumor_depth_df['gc_corrected_mean_depth'].values, 
        normal_depth_df['gc_corrected_mean_depth'].values,
    )

    for colname in (
        'mean_depth', 
        'norm_mean_depth', 
        'gc_corrected_mean_depth', 
        'sequenza_style_norm_mean_depth',
    ):
        result[f'tumor_{colname}'] = tumor_depth_df.loc[:, colname].array
        result[f'normal_{colname}'] = normal_depth_df.loc[:, colname].array

    if as_gr:
        return pr.PyRanges(result)
    else:
        return result


def merge_normal_tumor_bafs(tumor_baf_gr, normal_baf_gr, region_gr, as_gr=False):
    # convert df types
    #tumor_baf_gr = arg_into_gr(tumor_baf_gr)
    #normal_baf_gr = arg_into_gr(normal_baf_gr)
    #region_gr = arg_into_gr(region_gr)

    # sanity check
    assert 'baf' in tumor_baf_gr.columns
    assert 'baf' in normal_baf_gr.columns

    # remove annotation columns
    result = region_gr[[]]  

    # tumor
    result = pyranges_helper.join_newest(
        result, tumor_baf_gr,
        how='left', find_nearest=False, merge='mean',
        as_gr=False,
    )
    result.rename(columns={'baf': 'tumor_baf'}, inplace=True)
    result = pr.PyRanges(result)

    # normal
    result = pyranges_helper.join_newest(
        result, normal_baf_gr,
        how='left', find_nearest=False, merge='mean',
        as_gr=False,
    )
    result.rename(columns={'baf': 'normal_baf'}, inplace=True)

    # return
    if as_gr:
        return pr.PyRanges(result)
    else:
        return result


def handle_merged_baf_df(merged_baf_df):
    """Find regions with normal sample baf deviating from 0.5"""
    pass

    


#########################################
# plotting GCfrac - depth relationships #
#########################################

def arghandler_ylims(ydata, ylims):
    if ylims is None:
        return None
    else:
        return (
            min(ylims[0], min(ydata)), 
            max(ylims[1], max(ydata)),
        )


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
    """Intended to be used with short regions (~10kb)
    GC fraction value is calculated for every single base, with a window "gc_window" wide.
    """
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
        gc_breaks, gc_intervals = get_gc_breaks(n_gcbin)

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

    depth_df, selector = handle_outliers(depth_df, trim_limits=trim_limits, lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff)
    trimmed_depth_df = depth_df.iloc[selector, :]
    gc_breaks, gc_intervals = get_gc_breaks(n_gcbin)
    gcbin_depth_data, gcbin_average_depths, gcbin_norm_average_depths, cutresult = get_gcbin_data(trimmed_depth_df, gc_breaks, make_raw_data=True)
    fig_avg_scatter = make_average_scatter(gcbin_average_depths, ylims=scatter_ylims)
    if histogram_mode == '2d':
        fig_heatmap = make_heatmap_hist2d(trimmed_depth_df, n_gcbin, ylims=heatmap_ylims, vmax=heatmap_vmax)
    elif histogram_mode == '1d':
        fig_heatmap = make_heatmap(gcbin_depth_data, ylims=heatmap_ylims, vmax=heatmap_vmax)

    return fig_avg_scatter, fig_heatmap



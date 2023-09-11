import os
import collections
import random
import itertools
import functools
import multiprocessing
import operator

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
import handygenome.cnv.sequenza_handler as sequenza_handler
import handygenome.cnv.gcfraction as gcfraction
import handygenome.assemblyspec as libassemblyspec


MALE_HAPLOID_CHROMS = ('X', 'Y', 'chrX', 'chrY')
DEFAULT_CNT_WEIGHT = 5


class CCFInfo(
    collections.namedtuple('CCFInfo', ('CNm', 'ccf', 'mutated_allele'))
):
    pass


def get_normal_mean_ploidy(refver, is_female):
    """Only supports hg19 or hg38"""
    if is_female:
        return 2
    else:
        chromdict = common.DEFAULT_CHROMDICTS[refver]
        N_region_gr = libassemblyspec.get_N_region_gr(refver)
        par_gr = libassemblyspec.get_par_gr(refver)

        autosomal_chroms = [
            x for x in chromdict.contigs 
            if (
                (common.RE_PATS['assembled_chromosome'] is not None)
                and (x not in MALE_HAPLOID_CHROMS)
            )
        ]  # 1, 2, ...,
        X_chrom, Y_chrom = chromdict.XY_names
            # X, Y or chrX, chrY

        CNs = list()
        lengths = list()
        # autosomes
        for chrom in autosomal_chroms:
            CNs.append(2)
            N_length = N_region_gr[N_region_gr.Chromosome == chrom].lengths().sum()
            lengths.append(chromdict[chrom] - N_length)
        # X - PAR
        CNs.append(2)
        PAR_X_length = par_gr[par_gr.Chromosome == X_chrom].lengths().sum()
        lengths.append(PAR_X_length)
        # X - non-PAR
        CNs.append(1)
        X_N_length = N_region_gr[N_region_gr.Chromosome == X_chrom].lengths().sum()
        lengths.append(chromdict[X_chrom] - X_N_length - PAR_X_length)
        # Y
        CNs.append(1)
        PAR_Y_length = par_gr[par_gr.Chromosome == Y_chrom].lengths().sum()
        Y_N_length = N_region_gr[N_region_gr.Chromosome == Y_chrom].lengths().sum()
        lengths.append(chromdict[Y_chrom] - Y_N_length - PAR_Y_length)

        return np.average(CNs, weights=lengths)


def check_haploid(is_female, chrom):
    return (not is_female) and (chrom in MALE_HAPLOID_CHROMS)


def get_CNn_list(refver, is_female, chroms, start0s, end0s):
    """Assumes that N-masked regions or Y-PAR is not given as arguments"""
    chromdict = common.DEFAULT_CHROMDICTS[refver]
    X_chrom, Y_chrom = chromdict.XY_names
    if refver in libassemblyspec.PARS:
        pars = libassemblyspec.PARS[refver]
        par1x = pars['PAR1_X']
        par2x = pars['PAR2_X']
        par_unavailable = False
    else:
        par_unavailable = True

    is_haploid_dict = {
        chrom: check_haploid(is_female, chrom)
        for chrom in set(chroms)
    }

    # main
    result = list()
    for chrom, start0, end0 in zip(chroms, start0s, end0s):
        if is_haploid_dict[chrom]:
            if chrom == X_chrom:
                if par_unavailable:  # Assumes CNn is globally 1 if PAR is not available
                    CNn = 1
                else:
                    if any(
                        common.check_overlaps(start0, end0, par[1], par[2])
                        for par in (par1x, par2x)
                    ):  # any overlap with PAR is regarded as CNn == 2
                        CNn = 2
                    else:
                        CNn = 1
            elif chrom == Y_chrom:
                CNn = 1
            else:
                raise Exception(f'Non-XY haploid chromosome')
        else:
            CNn = 2

        result.append(CNn)

    return result


##############################################


def theoretical_depth_ratio(
    CNt, cellularity, tumor_mean_ploidy, CNn, normal_mean_ploidy, 
    tumor_avg_depth_ratio=1, normal_avg_depth_ratio=1,
):
    """Args:
        CNn, normal_ploidy must be nonzero
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
        return None
    else:
        numerator = (B * cellularity) + (1 - cellularity)
        denominator = (CNt * cellularity) + (CNn * (1 - cellularity))
        return numerator / denominator


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


def get_B(CNt, cellularity, baf):
    B_estimate = inverse_theoretical_baf(baf, cellularity, CNt)
    if B_estimate <= 0:
        B = 0
        theo_baf = theoretical_baf(CNt, B, cellularity)
        if theo_baf is None:
            B = None
            diff = None
        else:
            diff = abs(baf - theo_baf)
    else:
        B_candidate_upper = int(np.ceil(B_estimate))
        B_candidate_lower = int(np.floor(B_estimate))

        theo_baf_upper = theoretical_baf(CNt, B_candidate_upper, cellularity)
        theo_baf_lower = theoretical_baf(CNt, B_candidate_lower, cellularity)

        if (theo_baf_upper is None) and (theo_baf_lower is None):
            B = None
            diff = None
        elif (theo_baf_upper is not None) and (theo_baf_lower is None):
            B = B_candidate_upper
            diff = abs(baf - theo_baf_upper)
        elif (theo_baf_upper is None) and (theo_baf_lower is not None):
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


def get_CN_from_cp(cellularity, ploidy, depth_ratio, baf, CNt_weight):
    def save_cache(CNt_candidate, B_cache, segfit_score_cache, baf_diff_cache):
        B, baf_diff = get_B(CNt_candidate, cellularity, baf)
        if B is None:
            # drop this CNt candidate if B cannot be calculated
            return

        ratio_diff = abs(theoretical_depth_ratio(CNt_candidate, cellularity, ploidy) - depth_ratio)
        B_cache[CNt_candidate] = B
        baf_diff_cache[CNt_candidate] = baf_diff
        segfit_score_cache[CNt_candidate] = ratio_diff * CNt_weight + baf_diff

    # get initial CNt candidates
    CNt_estimate = inverse_theoretical_depth_ratio(depth_ratio, cellularity, ploidy)
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

    delta_ratio = delta_depth_ratio(cellularity, ploidy) * CNt_weight

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
        return None, None, None
    else:
        CNt, segfit_score = min(segfit_score_cache.items(), key=operator.itemgetter(1))
        B = B_cache[CNt]
        #selected_score = segfit_score_cache[CNt]

        return CNt, B, segfit_score


def get_CN_from_cp_wo_baf(cellularity, ploidy, depth_ratio, CNt_weight):
    CNt_estimate = inverse_theoretical_depth_ratio(depth_ratio, cellularity, ploidy)
    if CNt_estimate <= 0:
        CNt = 0
        diff = abs(theoretical_depth_ratio(CNt, cellularity, ploidy) - depth_ratio)
        segfit_score = diff * CNt_weight
    else:
        upper_candidate = int(np.ceil(CNt_estimate))
        lower_candidate = int(np.floor(CNt_estimate))

        upper_diff = abs(theoretical_depth_ratio(upper_candidate, cellularity, ploidy) - depth_ratio)
        lower_diff = abs(theoretical_depth_ratio(lower_candidate, cellularity, ploidy) - depth_ratio)
        if upper_diff <= lower_diff:
            CNt = upper_candidate
            segfit_score = upper_diff * CNt_weight
        else:
            CNt = lower_candidate
            segfit_score = lower_diff * CNt_weight

    return CNt, None, segfit_score  # B is None


def calc_cp_score(segment_df, cellularity, ploidy, is_female, CNt_weight):
    segfit_score_sum = 0
    start0_list = list()
    end0_list = list()
    CNt_list = list()
    for idx, row in segment_df.iterrows():
        if np.isnan(row['depth_mean']):
            continue

        if check_haploid(is_female, row['Chromosome']) or np.isnan(row['baf_mean']):
            CNt, B, segfit_score = get_CN_from_cp_wo_baf(cellularity, ploidy, row['depth_mean'], CNt_weight)
        else:
            CNt, B, segfit_score = get_CN_from_cp(cellularity, ploidy, row['depth_mean'], row['baf_mean'], CNt_weight)

        if segfit_score is not None:
            #segfit_score_sum += segfit_score
            segfit_score_sum -= segfit_score  # to make greater score better

        start0_list.append(row['Start'])
        end0_list.append(row['End'])
        CNt_list.append(CNt)

    seg_mean_ploidy = _get_mean_ploidy_from_values(start0_list, end0_list, CNt_list)
    ploidyfit_diff = abs(seg_mean_ploidy - ploidy)

    return {'segfit_score': segfit_score_sum, 'ploidyfit_diff': ploidyfit_diff}


def get_cp_score_dict(segment_df, is_female, CNt_weight=DEFAULT_CNT_WEIGHT, nproc=None):
    c_candidates = np.round(np.arange(0.01, 1.00, 0.01), 2)  # 0 and 1 are excluded
    p_candidates = np.round(np.arange(0.5, 7.1, 0.1), 1)  # 0.5 to 7.0
    cp_pairs = tuple(itertools.product(c_candidates, p_candidates))
    with multiprocessing.Pool(nproc) as pool:
        scorelist = pool.starmap(
            calc_cp_score, 
            ((segment_df, x[0], x[1], is_female, CNt_weight) for x in cp_pairs),
        )
    cp_scores = dict(zip(cp_pairs, scorelist))
    return cp_scores


def get_cpscore_peaks(cp_scores):
    c_list = sorted(set((x[0] for x in cp_scores.keys())))
    p_list = sorted(set((x[1] for x in cp_scores.keys())))
    cpscore_df_source = dict()
    for c in c_list:
        cpscore_df_source[c] = [cp_scores[(c, p)] for p in p_list]
    cpscore_df = pd.DataFrame.from_dict(cpscore_df_source)
    cpscore_df.index = p_list
    cpscore_df = -1 * cpscore_df  # to be run with get_fitresult_peaks function

    peaks = sequenza_handler.get_fitresult_peaks(cpscore_df)
    for x in peaks:
        x['score'] = -1 * x['lpp'] 
        del x['lpp']

    return peaks


def add_CN_to_segment(segment_df, cellularity, ploidy, is_female, CNt_weight=DEFAULT_CNT_WEIGHT):
    CNt_list = list()
    B_list = list()
    score_list = list()
    for idx, row in segment_df.iterrows():
        if np.isnan(row['depth_mean']):
            CNt, B, score = None, None, None
        else:
            if check_haploid(is_female, row['Chromosome']) or np.isnan(row['baf_mean']):
                CNt, B, score = get_CN_from_cp_wo_baf(cellularity, ploidy, row['depth_mean'], CNt_weight)
            else:
                CNt, B, score = get_CN_from_cp(cellularity, ploidy, row['depth_mean'], row['baf_mean'], CNt_weight)

        if CNt is None:
            CNt = np.nan
        if B is None:
            B = np.nan
        if score is None:
            score = np.nan

        CNt_list.append(CNt)
        B_list.append(B)
        score_list.append(score)

    return segment_df.assign(**{'CNt': CNt_list, 'B': B_list, 'score': score_list})


def add_theoreticals_to_segment(segment_df, cellularity, ploidy, is_female):
    theo_depthratio_list = list()
    theo_baf_list = list()
    for idx, row in segment_df.iterrows():
        if np.isnan(row['CNt']):
            theo_depthratio = np.nan
            theo_baf = np.nan
        else:
            if check_haploid(is_female, row['Chromosome']):
                CNn = 1
            else:
                CNn = 2

            # depthratio
            theo_depthratio = theoretical_depth_ratio(row['CNt'], cellularity, ploidy, CNn=CNn)
            # baf
            if np.isnan(row['B']):
                theo_baf = np.nan
            else:
                theo_baf = theoretical_baf(row['CNt'], row['B'], cellularity, CNn=CNn)
                if theo_baf is None:
                    theo_baf = np.nan

        theo_depthratio_list.append(theo_depthratio)
        theo_baf_list.append(theo_baf)

    return segment_df.assign(
        **{'predicted_depth_ratio': theo_depthratio_list, 'predicted_baf': theo_baf_list}
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
    return breaks


def remove_outliers(depth_df, trim_limits=(0.05, 0.05), lower_cutoff=None, upper_cutoff=None):
    """Absolute cutoff-based trimming is done first, then proportional trimming is later
    """
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


def get_gcbin_data(depth_df, gc_breaks):
    """This does not apply weights when calculating averages. Suitable for even-sized bins"""
    # cut by gc bins
    cutresult = pd.cut(depth_df['GC'], bins=gc_breaks)
    # make results
    gcbin_depth_data = {intv: list() for intv in cutresult.dtype.categories}
    for intv, depth in zip(cutresult, depth_df['mean_depth']):
        gcbin_depth_data[intv].append(depth)

    gcbin_average_depths = {
        intv: (np.nan if len(data) == 0 else np.mean(data))
        for intv, data in gcbin_depth_data.items()
    }
    gcbin_average_depths = pd.Series(gcbin_average_depths)
    gcbin_norm_average_depths = gcbin_average_depths / gcbin_average_depths.mean(skipna=True)

    return gcbin_depth_data, gcbin_average_depths, gcbin_norm_average_depths, cutresult


def postprocess_depth_df(
    depth_df, 
    *,
    fasta=None, 
    refver=None,
    binsize=None,
    gc_window=None,

    nan_lower_cutoff=None,
    nan_upper_cutoff=None,

    gcdata_trim_limits=(0.05, 0.05),
    gcdata_lower_cutoff=None,
    gcdata_upper_cutoff=None,

    n_gcbin=100, 
    as_gr=True, 
):
    """Args:
        fasta, refver, binsize: Only used for gc fraction calculation.
            Allowed argument usage: "fasta" only OR "refver" and "binsize" only

    Returns:
        postprocessed df & gcbin average depth dict
        Adds 4 columns:
            - GC
            - norm_mean_depth
            - gc_corrected_mean_depth
            - sequenza_style_norm_mean_depth
    """
    # sanity check
    if not (
        (
            (fasta is not None )
            and (refver is None)
            and (binsize is None)
        )
        or (
            (fasta is None )
            and (refver is not None)
            and (binsize is not None)
        )
    ):
        raise Exception(
            f'Allowed argument usage for GC fraction annotation: '
            f'1) use "fasta", do NOT use "refver" and "binsize" ; '
            f'2) use "refver" and "binsize", do NOT use "fasta"'
        )

    assert 'mean_depth' in depth_df.columns, f'"depth_df" must include a column named "mean_depth"'

    # replace outliers with np.nan
    if nan_lower_cutoff is not None:
        depth_df.loc[depth_df['mean_depth'] < nan_lower_cutoff, ['mean_depth']] = np.nan
    if nan_upper_cutoff is not None:
        depth_df.loc[depth_df['mean_depth'] >= nan_upper_cutoff, ['mean_depth']] = np.nan

    # GC
    gcfraction.add_gc_calculating(depth_df, fasta, window=gc_window)

    # norm_mean_depth
    depth_df_wo_nan = depth_df.loc[~np.isnan(depth_df['mean_depth']), :]
    global_mean_depth = np.average(depth_df_wo_nan['mean_depth'], weights=(depth_df_wo_nan['End'] - depth_df_wo_nan['Start']))
    depth_df.insert(depth_df.shape[1], 'norm_mean_depth', (depth_df['mean_depth'] / global_mean_depth).array)

    # get gc depth data
    gc_breaks = get_gc_breaks(n_gcbin)
    trimmed_depth_df = remove_outliers(depth_df, trim_limits=gcdata_trim_limits, lower_cutoff=gcdata_lower_cutoff, upper_cutoff=gcdata_upper_cutoff)
    gcbin_depth_data, gcbin_average_depths, gcbin_norm_average_depths, cutresult = get_gcbin_data(trimmed_depth_df, gc_breaks)

    # gc_corrected_mean_depth
    cutresult = pd.cut(depth_df['GC'], bins=gc_breaks)
    gcbin_norm_average_depths_selected = gcbin_norm_average_depths.loc[cutresult]
    depth_df.insert(depth_df.shape[1], 'gc_corrected_mean_depth', (depth_df['mean_depth'].array / gcbin_norm_average_depths_selected.array))

    # sequenza_style_norm_mean_depth
    gcbin_average_depths_selected = gcbin_average_depths.loc[cutresult]
    depth_df.insert(depth_df.shape[1], 'sequenza_style_norm_mean_depth', (depth_df['mean_depth'].array / gcbin_average_depths_selected.array))

    if as_gr:
        return pr.PyRanges(depth_df, int64=False), gcbin_average_depths
    else:
        return depth_df, gcbin_average_depths


def get_processed_depth_df(bam_path, fasta, region_gr, gc_window=None, outlier_cutoffs='wgs', donot_subset_bam=False, as_gr=False):
    # set params
    assert isinstance(outlier_cutoffs, (tuple, list)) or outlier_cutoffs in ('wgs', 'panel')

    if outlier_cutoffs == 'wgs':
        lower_cutoff = 0
        upper_cutoff = 2000
    elif outlier_cutoffs == 'panel':
        lower_cutoff = 50
        upper_cutoff = np.inf
    else:
        lower_cutoff, upper_cutoff = outlier_cutoffs

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

        gcdata_trim_limits=(0.05, 0.05),
        gcdata_lower_cutoff=lower_cutoff,
        gcdata_upper_cutoff=upper_cutoff,
        nan_lower_cutoff=lower_cutoff,
        nan_upper_cutoff=upper_cutoff,

        n_gcbin=100,
        as_gr=as_gr, 
    )

    return depth_gr, gcbin_average_depths


def make_depth_ratio_df(tumor_depth_df, normal_depth_df, as_gr=False):
    # sanity check
    assert isinstance(tumor_depth_df, pd.DataFrame) and isinstance(normal_depth_df, pd.DataFrame)
    assert all(
        {'mean_depth', 'GC', 'norm_mean_depth', 'gc_corrected_mean_depth', 'sequenza_style_norm_mean_depth'}.issubset(df.columns)
        for df in (tumor_depth_df, normal_depth_df)
    )

    compared_columns = ['Chromosome', 'Start', 'End', 'GC']
    if not (tumor_depth_df.loc[:, compared_columns] == normal_depth_df.loc[:, compared_columns]).to_numpy().all():
        raise Exception(f'The columns {compared_columns} of the two input dataframes must be identical.')

    # main
    result = tumor_depth_df.loc[:, compared_columns].copy()
    result['depth_ratio_sequenzastyle'] = tumor_depth_df['sequenza_style_norm_mean_depth'] / normal_depth_df['sequenza_style_norm_mean_depth']
    result['depth_ratio_mystyle'] = tumor_depth_df['gc_corrected_mean_depth'] / normal_depth_df['gc_corrected_mean_depth']
    for colname in ('mean_depth', 'norm_mean_depth', 'gc_corrected_mean_depth', 'sequenza_style_norm_mean_depth'):
        result[f'tumor_{colname}'] = tumor_depth_df.loc[:, colname]
        result[f'normal_{colname}'] = normal_depth_df.loc[:, colname]

    if as_gr:
        return pr.PyRanges(result)
    else:
        return result


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



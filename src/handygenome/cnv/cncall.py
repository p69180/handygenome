import inspect
import os
import collections
import functools

import numpy as np
import pandas as pd

import handygenome.tools as tools
import handygenome.deco as deco
import handygenome.refgenome.refgenome as refgenome
from handygenome.refgenome.refgenome import NoXYError
import handygenome.refgenome.refverfile as refverfile
import handygenome.genomedf.genomedf as libgenomedf
from handygenome.genomedf.genomedf import GenomeDataFrame as GDF


###############
# terminology #
###############

'''
CNg: total copy number in germline
CNt: total copy number of the target (tumor) sample
Bg: B allele copy number in germline
Bt: B allele copy number of the target (tumor) sample
'''


#############
# constants #
#############

DEFAULT_CNG_COLNAME = 'default_CNg'
DEFAULT_BG_COLNAME = 'default_Bg'
MALE_HAPLOID_CHROMS = ('X', 'Y', 'chrX', 'chrY')
DEFAULT_CNT_WEIGHT = 5


###############
# CCF classes #
###############

class CCFInfo(
    collections.namedtuple('CCFInfo', ('CNm', 'ccf', 'mutated_allele'))
):
    pass


class CPPair(
    collections.namedtuple('CPPair', ('cellularity', 'ploidy'))
):
    pass


##############
# exceptions #
##############

class DefaultCNgUnavailableError(Exception):
    pass


##############
# decorators #
##############

def deco_cellularity_sanitycheck(func):
    sig = inspect.signature(func)
    if 'cellularity' not in sig.parameters.keys():
        raise Exception(f'"cellularity" is not included in the paramter list.')
    
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        ba = sig.bind(*args, **kwargs)
        ba.apply_defaults()
        if ba.arguments['cellularity'] == 0:
            raise Exception(f'"cellularity" must not be 0.')

        return func(*ba.args, **ba.kwargs)

    return wrapper


################
# miscellanous #
################

def check_haploid(is_female, chrom):
    return (
        (not is_female) 
        and refgenome.check_xy_chrom(chrom)
    )


############################
# germline sample solution #
############################

def get_default_CNg_Bg_diploid(refver, is_female):
    """Assumptions
        1) X, Y is sex chromosome
        2) Male has one copy of each sex chromosomes
        3) pseudoautosomal region is available
    """
    chromdict = refgenome.get_chromdict(refver)
    autosomal_chroms = [
        x for x in chromdict.contigs 
        if (
            refgenome.check_assembled_chrom(x)
            and (not refgenome.check_xy_chrom(x))
        )
    ]  # 1, 2, ...,
    try:
        X_chrom, Y_chrom = chromdict.XY_names
    except NoXYError:
        raise DefaultCNgUnavailableError(f'Default germline CN information is not available for this reference version.')

    par_gdf = refverfile.get_par_gdf(refver)

    # autosomal
    all_chrom_gdf = GDF.all_regions(refver, assembled_only=True)
    autosomal_gdf = all_chrom_gdf.subset_chroms(autosomal_chroms)
    autosomal_gdf[DEFAULT_CNG_COLNAME] = float(2)
    autosomal_gdf[DEFAULT_BG_COLNAME] = float(1)

    # sex
    if is_female:
        sex_gdf = all_chrom_gdf.subset_chroms([X_chrom, Y_chrom])
        is_X = (sex_gdf.chroms == X_chrom)
        sex_gdf.loc[is_X, DEFAULT_CNG_COLNAME] = float(2)
        sex_gdf.loc[is_X, DEFAULT_BG_COLNAME] = float(1)
        sex_gdf.loc[~is_X, DEFAULT_CNG_COLNAME] = float(0)
        sex_gdf.loc[~is_X, DEFAULT_BG_COLNAME] = np.nan
    else:
        all_sex_gdf = all_chrom_gdf.subset_chroms([X_chrom, Y_chrom])
        nonpar_sex_gdf = all_sex_gdf.subtract(par_gdf)
        nonpar_sex_gdf[DEFAULT_CNG_COLNAME] = float(1)
        nonpar_sex_gdf[DEFAULT_BG_COLNAME] = np.nan

        par_sex_gdf = par_gdf.copy()
        is_X = (par_sex_gdf['Chromosome'] == X_chrom)
        par_sex_gdf.loc[is_X, DEFAULT_CNG_COLNAME] = float(2)
        par_sex_gdf.loc[is_X, DEFAULT_BG_COLNAME] = float(1)
        par_sex_gdf.loc[~is_X, DEFAULT_CNG_COLNAME] = float(0)
        par_sex_gdf.loc[~is_X, DEFAULT_BG_COLNAME] = np.nan

        sex_gdf = GDF.concat([nonpar_sex_gdf, par_sex_gdf])

    # result
    result = GDF.concat([autosomal_gdf, sex_gdf])
    result.sort()
    return result


def clonal_solution_from_germline(cnv_seg_gdf, is_female, ploidy):
    # add default CNg and Bg
    cnv_seg_gdf_copy = cnv_seg_gdf.copy()

    # find clonal CNt
    onecopy_depth = cnv_seg_gdf_copy.find_onecopy_depth(is_female, ploidy)
    K = get_K_from_onecopy_depth(onecopy_depth)
    clonal_CNt = find_clonal_CNt(
        corrected_depth=cnv_seg_gdf_copy.norm_depth_mean, 
        cellularity=1, 
        K=K,
        CNg=None, 
    )

    # find clonal B
    B_isnan_selector = np.isnan(cnv_seg_gdf_copy.get_default_Bg(is_female))
    B_notnan_selector = ~B_isnan_selector

    clonal_Bt = np.empty(cnv_seg_gdf_copy.nrow, dtype=float)
    clonal_Bt[B_isnan_selector] = np.nan

    relevant_bafs = cnv_seg_gdf_copy.get_corrected_baf()[B_notnan_selector]
    assert not np.isnan(relevant_bafs).any()
    clonal_Bt[B_notnan_selector] = find_clonal_Bt(
        baf=relevant_bafs, 
        CNt=clonal_CNt[B_notnan_selector], 
        cellularity=1, 
        CNg=None, 
        Bg=None,
    )

    return clonal_CNt, clonal_Bt


def solution_from_targetsample():
    average_CNt = depth_seg_gdf.norm_depth / onecopy_depth
    clonal_CNt = np.rint(average_CNt)

    # 4. determine regions where clonal solutions are unsatisfactory
    depth_diffs = np.abs((onecopy_depth * clonal_CNt) - depth_seg_gdf.norm_depth)
    diff_thresholds = (onecopy_depth * threshold_factor) * average_CNt
        # Where "average_CNt" is 2, "diff_thresholds" value becomes:
        # 2 * onecopy_depth * threshold_factor
    unfit_selector = depth_diffs > diff_thresholds
    unfit_selector_idxs = np.nonzero(unfit_selector)[0]

    # 5. find subclonal solutions
    unfit_average_CNt = average_CNt[unfit_selector_idxs]
    unfit_clonal_CNt = depth_seg_gdf['default_CNg'][unfit_selector_idxs]

    unfit_clonal_CNt, unfit_subclonal_CNt = get_CN_candidates_type2(unfit_average_CNt, clonal_CNt=unfit_clonal_CNt)

    return depth_seg_gdf


@deco_cellularity_sanitycheck
def find_clonal_CNt(corrected_depth, cellularity, K, CNg=None):
    """Assume 'corrected_depth' is corrected so that position-dependent bias has been removed.
    It can be formulated like:
        'corrected_depth' = K * (CNg * (1 - cellularity) + CNt * cellularity)

    derivation:
        corrected_depth = K * (CNg * (1 - cellularity) + CNt * cellularity)
        corrected_depth / K = CNg * (1 - cellularity) + CNt * cellularity
        (corrected_depth / K) - (CNg * (1 - cellularity)) = CNt * cellularity
        CNt = ( (corrected_depth / K) - (CNg * (1 - cellularity)) ) / cellularity
    """
    if (cellularity != 1) and (CNg is None):
        raise Exception(f'When cellularity is not 1, "CNg" must be given.')

    if cellularity == 1:
        CNt = (corrected_depth / K)
    else:
        CNt = (
            (corrected_depth / K) 
            - (CNg * (1 - cellularity))
        ) / cellularity

    CNt = np.rint(CNt)
    return CNt


@deco_cellularity_sanitycheck
def find_clonal_Bt(baf, CNt, cellularity, CNg=None, Bg=None):
    """
    derivation:
        baf = (
            (Bt * cellularity + Bg * (1 - cellularity))
            / (CNt * cellularity + CNg * (1 - cellularity))
        )
        Bt * cellularity + Bg * (1 - cellularity) = (
            baf * (CNt * cellularity + CNg * (1 - cellularity))
        )
        Bt * cellularity = (
            (baf * (CNt * cellularity + CNg * (1 - cellularity)))
            - (Bg * (1 - cellularity))
        )
        Bt = (
            (baf * (CNt * cellularity + CNg * (1 - cellularity)))
            - (Bg * (1 - cellularity))
        ) / cellularity
        Bt = (
            (baf * (CNt + CNg * (1 / cellularity - 1)))
            - (Bg * (1 / cellularity - 1))
        )
    """

    if (
        (cellularity != 1) 
        and (
            (CNg is None)
            or (Bg is None)
        )
    ):
        raise Exception(f'When cellularity is not 1, "CNg" and "Bg" must be given.')

    c = (1 / cellularity) - 1
    if c == 0:
        Bt = baf * CNt 
    else:
        Bt = (baf * (CNt + CNg * c)) - (Bg * c)

    Bt = np.rint(Bt)
    Bt = np.clip(Bt, 0, np.floor(CNt / 2))
    return Bt


def get_K_from_onecopy_depth(onecopy_depth):
    """Only make sense with germline wgs sample

    derivation:
        corrected_depths = K * (CNg * (1 - cellularity) + CNt * cellularity)
        corrected_depths = K * (CNt)  (since cellularity == 1)
        K = corrected_depths / CNt
        K = onecopy_depth / 1
    """
    return onecopy_depth


#####################################
# helpers - CN candidate generators #
#####################################

@deco.get_deco_atleast1d(['average_CNt'])
def get_CN_candidates_type1(average_CNt, *, min_candidate_num=5, candidate_num_fraction=0.5):
    """lower and upper candidates are given
    """
    assert average_CNt.ndim == 1

    num_candidates, max_dimsize, CN_shifts = get_CN_candidate_shifts(
        average_CNt, min_candidate_num, candidate_num_fraction
    )
    upperCN_borders, lowerCN_borders = tools.get_nearest_integer_bounds(average_CNt, elevate=True)

    upper_CNt = upperCN_borders[:, np.newaxis] + CN_shifts
    lower_CNt = lowerCN_borders[:, np.newaxis] - CN_shifts
    lower_CNt[lower_CNt < 0] = np.nan

    return lower_CNt, upper_CNt


@deco.get_deco_squeeze_atleast1d(['average_CNt'])
def get_CN_candidates_type2(average_CNt, *, clonal_CNt=None, min_candidate_num=5, candidate_num_fraction=0.5):
    """Candidates are labeled as clonal or subclonal. Clonal values can be optionally given.
    """
    assert average_CNt.ndim == 1

    num_candidates, max_dimsize, CN_shifts = get_CN_candidate_shifts(
        average_CNt, min_candidate_num, candidate_num_fraction
    )

    if clonal_CNt is None:
        clonal_CNt = np.rint(average_CNt)
    else:
        assert average_CNt.shape == clonal_CNt.shape

    gt_clonal_selector = (average_CNt > clonal_CNt)

    gt_clonal_subCNt = clonal_CNt[gt_clonal_selector][:, np.newaxis] + CN_shifts[gt_clonal_selector]
    lt_clonal_subCNt = clonal_CNt[~gt_clonal_selector][:, np.newaxis] - CN_shifts[~gt_clonal_selector]
    lt_clonal_subCNt[lt_clonal_subCNt < 0] = np.nan

    subclonal_CNt = np.empty((len(average_CNt), max_dimsize), dtype=float)
    subclonal_CNt[gt_clonal_selector, :] = gt_clonal_subCNt
    subclonal_CNt[~gt_clonal_selector, :] = lt_clonal_subCNt

    return clonal_CNt, subclonal_CNt


def get_CN_candidate_shifts(average_CNt, min_candidate_num, candidate_num_fraction):
    num_candidates = np.maximum(
        min_candidate_num, 
        np.rint(average_CNt * candidate_num_fraction),
    )
    max_dimsize = np.nanmax(num_candidates, axis=None).astype(int)

    CN_shifts = np.tile(
        np.arange(max_dimsize)[np.newaxis, :], 
        (len(average_CNt), 1),
    )
    CN_shifts = CN_shifts.astype(float)
    nan_selector = np.stack(
        [
            (
                np.repeat(True, max_dimsize)
                if np.isnan(x) else
                np.repeat([False, True], (x, max_dimsize - x)) 
            )
            for x in num_candidates
        ], 
        axis=0,
    )
    CN_shifts[nan_selector] = np.nan

    return num_candidates, max_dimsize, CN_shifts


#############################
# helpers - B determination #
#############################

def get_B_from_clonal_CN(clonal_CNt, bafs):
    pass
    




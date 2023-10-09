import inspect
import os
import collections
import functools
import warnings

import numpy as np
import pandas as pd

import handygenome.tools as tools
import handygenome.deco as deco
import handygenome.refgenome.refgenome as refgenome
from handygenome.refgenome.refgenome import NoXYError
import handygenome.refgenome.refverfile as refverfile
import handygenome.genomedf.genomedf as libgenomedf
import handygenome.cnv.baf as libbaf
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

@functools.cache
def get_default_CNg_Bg_diploid(refver, is_female):
    """Assumptions
        1) X, Y is sex chromosome
        2) Male has one copy of each sex chromosomes
        3) pseudoautosomal region is available
    Results:
        CN: 0 in female Y, male Y PAR
        B: 0 where CN is 0 or 1
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
        sex_gdf.loc[~is_X, DEFAULT_BG_COLNAME] = float(0)
    else:
        all_sex_gdf = all_chrom_gdf.subset_chroms([X_chrom, Y_chrom])
        nonpar_sex_gdf = all_sex_gdf.subtract(par_gdf)
        nonpar_sex_gdf[DEFAULT_CNG_COLNAME] = float(1)
        nonpar_sex_gdf[DEFAULT_BG_COLNAME] = float(0)

        par_sex_gdf = par_gdf.copy()
        is_X = (par_sex_gdf['Chromosome'] == X_chrom)
        par_sex_gdf.loc[is_X, DEFAULT_CNG_COLNAME] = float(2)
        par_sex_gdf.loc[is_X, DEFAULT_BG_COLNAME] = float(1)
        par_sex_gdf.loc[~is_X, DEFAULT_CNG_COLNAME] = float(0)
        par_sex_gdf.loc[~is_X, DEFAULT_BG_COLNAME] = float(0)

        sex_gdf = GDF.concat([nonpar_sex_gdf, par_sex_gdf])

    # result
    result = GDF.concat([autosomal_gdf, sex_gdf])
    result.sort()
    return result


def add_default_CNg_Bg(gdf, is_female, inplace=True, nproc=1):
    if DEFAULT_CNG_COLNAME not in gdf.columns:
        CNg_Bg_gdf = get_default_CNg_Bg_diploid(gdf.refver, is_female)
        annotated_gdf = gdf.join(
            CNg_Bg_gdf, 
            how='left',
            right_gdf_cols=[DEFAULT_CNG_COLNAME, DEFAULT_BG_COLNAME],
            merge='longest',
            overlapping_length=True,
            omit_N=True,
            suffixes={'longest': ''},
            nproc=nproc,
        )
        if inplace:
            gdf.assign_df(annotated_gdf.df)
        else:
            return annotated_gdf
    else:
        if inplace:
            pass
        else:
            return gdf


def get_default_CNg(seg_gdf, is_female, nproc=1):
    if DEFAULT_CNG_COLNAME not in seg_gdf.columns:
        add_default_CNg_Bg(seg_gdf, is_female, inplace=True, nproc=nproc)
    return seg_gdf[DEFAULT_CNG_COLNAME]


def get_default_Bg(seg_gdf, is_female, nproc=1):
    if DEFAULT_BG_COLNAME not in seg_gdf.columns:
        add_default_CNg_Bg(seg_gdf, is_female, inplace=True, nproc=nproc)
    return seg_gdf[DEFAULT_BG_COLNAME]


def clonal_solution_from_germline(cnv_seg_gdf, is_female, ploidy, nproc=1):
    cnv_seg_gdf_copy = cnv_seg_gdf.copy()

    # find clonal CNt
    onecopy_depth = cnv_seg_gdf_copy.find_onecopy_depth(is_female, ploidy)
    K = get_K_from_onecopy_depth_germline(onecopy_depth)
    clonal_CNt = find_clonal_CNt(
        corrected_depth=cnv_seg_gdf_copy.norm_depth_mean, 
        cellularity=1, 
        K=K,
        CNg=None, 
    )

    # replace CNt with 0 where default germline CN is 0
    try:
        default_CNg = get_default_CNg(cnv_seg_gdf_copy, is_female, nproc=nproc)
    except DefaultCNgUnavailableError:
        pass
    else:
        clonal_CNt[(default_CNg == 0)] = 0

    # find clonal B for each baf index
    try:
        default_Bg = get_default_Bg(cnv_seg_gdf_copy, is_female)
        #B_isnan_selector = np.isnan(
        #    get_default_Bg(cnv_seg_gdf_copy, is_female)
        #)
    except DefaultCNgUnavailableError:
        invalid_B_selector = np.repeat(False, cnv_seg_gdf_copy.nrow)
        #B_isnan_selector = np.repeat(False, cnv_seg_gdf_copy.nrow)
    else:
        invalid_B_selector = (default_Bg == 0)
        #np.repeat(False, cnv_seg_gdf_copy.nrow)
    valid_B_selector = np.logical_not(invalid_B_selector)
    #B_notnan_selector = np.logical_not(B_isnan_selector)

    clonal_Bt = np.empty((cnv_seg_gdf_copy.nrow, ploidy - 1), dtype=float)
    clonal_Bt[invalid_B_selector, :] = 0

    relevant_CNt = clonal_CNt[valid_B_selector]
    for idx, baf_index in enumerate(libbaf.get_baf_indexes(ploidy)):
        relevant_bafs = cnv_seg_gdf_copy.get_corrected_baf(baf_index)[valid_B_selector]
        #assert not np.isnan(relevant_bafs).any()  
            # this may be nan in non-assembled chroms like GL* where baf rawdata is unavailable
        clonal_Bt[valid_B_selector, idx] = find_clonal_Bt(
            baf=relevant_bafs, 
            CNt=relevant_CNt, 
            cellularity=1, 
            CNg=None, 
            Bg=None,
        )

    return clonal_CNt, clonal_Bt


################################
# non-germline sample solution #
################################

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


###########
# helpers #
###########

def get_deco_CNg_Bg(CNg_fillvalue, Bg_fillvalue):
    def decorator(func):
        sig = inspect.signature(func)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()
            if 'CNg' in ba.arguments:
                old_val = ba.arguments['CNg']
                if old_val is not None:
                    new_val = old_val.copy()
                    new_val[np.isnan(new_val)] = CNg_fillvalue
                    ba.arguments['CNg'] = new_val
            if 'Bg' in ba.arguments:
                old_val = ba.arguments['Bg']
                if old_val is not None:
                    new_val = old_val.copy()
                    new_val[np.isnan(new_val)] = Bg_fillvalue
                    ba.arguments['Bg'] = new_val

            return func(*ba.args, **ba.kwargs)

        return wrapper

    return decorator

    

@get_deco_CNg_Bg(CNg_fillvalue=2, Bg_fillvalue=1)
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
    CNg_given = (CNg is not None)
    if (cellularity != 1) and (not CNg_given):
        raise Exception(f'When cellularity is not 1, "CNg" must be given.')

    # main
    if cellularity == 1:
        CNt = (corrected_depth / K)
    else:
        CNt = (
            (corrected_depth / K) 
            - (CNg * (1 - cellularity))
        ) / cellularity

    # correction
    CNt = np.clip(np.rint(CNt), 0, None)
    if CNg_given:
        CNt[CNg == 0] = 0

    return CNt


@get_deco_CNg_Bg(CNg_fillvalue=2, Bg_fillvalue=1)
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
    # sanitycheck
    CNg_given = (CNg is not None)
    Bg_given = (Bg is not None)
    if (
        (cellularity != 1) 
        and (
            (not CNg_given) or (not Bg_given)
        )
    ):
        raise Exception(f'When cellularity is not 1, "CNg" and "Bg" must be given.')

    # main
    c = (1 / cellularity) - 1
    if c == 0:
        Bt = baf * CNt 
    else:
        Bt = (baf * (CNt + CNg * c)) - (Bg * c)

    # correction
    Bt = np.clip(np.rint(Bt), 0, np.floor(CNt / 2))
    if CNg_given or Bg_given:
        if CNg_given:
            zero_selector_CNg = (CNg == 0)
        else:
            zero_selector_CNg = np.repeat(False, len(CNg))

        if Bg_given:
            zero_selector_Bg = (Bg == 0)
        else:
            zero_selector_Bg = np.repeat(False, len(Bg))

        zero_selector = np.logical_or(zero_selector_CNg, zero_selector_Bg)
        Bt[zero_selector] = 0

    return Bt


def get_K_from_onecopy_depth_germline(onecopy_depth):
    """Only make sense with germline wgs sample

    derivation:
        corrected_depths = K * (CNg * (1 - cellularity) + CNt * cellularity)
        corrected_depths = K * (CNt)  (since cellularity == 1)
        K = corrected_depths / CNt
        K = onecopy_depth / 1
    """
    return onecopy_depth


@get_deco_CNg_Bg(CNg_fillvalue=2, Bg_fillvalue=1)
def find_K_and_cellularity(*, depth1, CNt1, CNg1, depth2, CNt2, CNg2):
    """Assume 'depth1' and 'depth2' are corrected so that position-dependent bias has been removed.

    derivation:
        assume depth1, depth2, K are all nonzero

        depth1 = K * (CNg1 * (1 - cellularity) + CNt1 * cellularity)
        depth2 = K * (CNg2 * (1 - cellularity) + CNt2 * cellularity)
        ------------------------------------------------------------
        depth1 = K * (CNg1 * (1 - cellularity) + CNt1 * cellularity)
        1 / K = (1 / depth1) * (CNg1 + cellularity * (CNt1 - CNg1))
        1 / K = (1 / depth2) * (CNg2 + cellularity * (CNt2 - CNg2))
        ------------------------------------------------------------
        (1 / depth1) * (CNg1 + cellularity * (CNt1 - CNg1)) = (1 / depth2) * (CNg2 + cellularity * (CNt2 - CNg2))
        depth2 * (CNg1 + cellularity * (CNt1 - CNg1)) = depth1 * (CNg2 + cellularity * (CNt2 - CNg2))
        depth2 * cellularity * (CNt1 - CNg1) - depth1 * cellularity * (CNt2 - CNg2) = depth1 * CNg2 - depth2 * CNg1
        cellularity * (depth2 * (CNt1 - CNg1) - depth1 * (CNt2 - CNg2)) = depth1 * CNg2 - depth2 * CNg1
        cellularity = (
            (depth1 * CNg2 - depth2 * CNg1)
            / (depth2 * (CNt1 - CNg1) - depth1 * (CNt2 - CNg2))
        )
        ------------------------------------------------------------
        K = 1 / ((1 / depth1) * (CNg1 + cellularity * (CNt1 - CNg1)))
    """
    # cellularity
    cellularity_numer = depth1 * CNg2 - depth2 * CNg1
    cellularity_denom = depth2 * (CNt1 - CNg1) - depth1 * (CNt2 - CNg2)
    if cellularity_denom == 0:
        raise Exception(f'Invalid input value combination: cellularity denominator becomes zero')

    cellularity = cellularity_numer / cellularity_denom

    # K
    K_inverse = (1 / depth1) * (CNg1 + cellularity * (CNt1 - CNg1))
    if K_inverse == 0:
        raise Exception(f'Invalid input value combination: inverse of K becomes zero')
    K = 1 / K_inverse

    return K, cellularity


@get_deco_CNg_Bg(CNg_fillvalue=2, Bg_fillvalue=1)
def get_predicted_depth(CNt, CNg, cellularity, K):
    return K * (CNg * (1 - cellularity) + CNt * cellularity)


@get_deco_CNg_Bg(CNg_fillvalue=2, Bg_fillvalue=1)
def get_predicted_baf(CNt, Bt, CNg, Bg, cellularity):
    """Where CNt == 0 and CNg == 0, result is np.nan"""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return (
            (Bt * cellularity + Bg * (1 - cellularity))
            / (CNt * cellularity + CNg * (1 - cellularity))
        )


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
    




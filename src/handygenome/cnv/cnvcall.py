import inspect
import os
import collections
import functools
import warnings
import itertools

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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


class KCError(Exception):
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

    import handygenome.cnv.baf as libbaf
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


def solution_from_segments(
    depths, 
    bafs, 

    num_ndiv_cand=10,
    num_CNt0_cand=6,
    ndivs_offset_step=0.01,

    depthfit_cutoff=0.15,
    bafdiff_frac=0.1,
    CNt_handicap_factor=1,

    CNg=2,
    Bg=1,
):
    # make grid variables
    CNt0 = np.arange(num_CNt0_cand)
    ndivs = np.arange(num_ndiv_cand) + 1

    ## add perturbations to ndivs
    if ndivs_offset_step is not None:
        assert ndivs_offset_step < depthfit_cutoff
        ndivs_offset = np.arange(0, depthfit_cutoff, ndivs_offset_step)
        ndivs_offset = np.unique(np.concatenate([ndivs_offset, -ndivs_offset]))
        ndivs = (ndivs[:, np.newaxis] + ndivs_offset).ravel()
        ndivs = ndivs[ndivs > 0]

    dd = depths.ptp() / ndivs

    CNt0, dd = np.meshgrid(CNt0, dd, indexing='ij')

    # make K, c
    depth0 = depths.min()
    K, cellularity = find_K_and_cellularity(depth1=depth0, CNt1=CNt0, CNg1=2, depth2=(depth0 + dd), CNt2=(CNt0 + 1), CNg2=2)

    # expand along a new axis which represent different segments
    n_segs = len(depths)
    def expand_helper(arr):
        return np.repeat(np.expand_dims(arr, 0), n_segs, axis=0)
    CNt0 = expand_helper(CNt0)
    dd = expand_helper(dd)
    cellularity = expand_helper(cellularity)

    # expand depths and bafs
    depths = np.expand_dims(depths, (1, 2))
    bafs = np.expand_dims(bafs, (1, 2))

    # make CNt1 for each segment
    depthdiff_quotients = (depths - depth0) / dd
    CNt1 = (CNt0 + np.rint(depthdiff_quotients))

    # determine Bt for each segment
    half_CNt1 = np.floor(CNt1 / 2).astype(int)
    Bt_offsets = np.arange(half_CNt1.max())
    Bt = (np.expand_dims(half_CNt1, -1) - Bt_offsets).astype(float)
    Bt[Bt < 0] = np.nan

    pred_bafs = get_predicted_baf(
        CNt=np.expand_dims(CNt1, -1), 
        Bt=Bt, 
        CNg=CNg, 
        Bg=Bg, 
        cellularity=np.expand_dims(cellularity, -1),
    )
    ## choose indexes with the least difference along the last axis
    baf_diffs = np.abs(pred_bafs - np.expand_dims(bafs, -1))
    # argmins = np.nanargmin(baf_diffs, axis=-1, keepdims=True)
    argmins = nanargmin(baf_diffs, axis=-1, keepdims=True)
    Bt = np.squeeze(np.take_along_axis(Bt, argmins, axis=-1), axis=-1)
    pred_bafs = np.squeeze(np.take_along_axis(pred_bafs, argmins, axis=-1), axis=-1)

    # make onecopy bafdiff
    onecopy_bafdiffs = get_predicted_baf(
        CNt=CNt1, 
        Bt=1, 
        CNg=CNg, 
        Bg=Bg, 
        cellularity=cellularity,
    )

    # get diffs between pred and obs bafs
    baf_errors = np.abs(pred_bafs - bafs)

    # get valid selectors
    depthfit_valid_flags = np.abs(depthdiff_quotients - np.rint(depthdiff_quotients)) < depthfit_cutoff
    baffit_valid_flags = baf_errors < (onecopy_bafdiffs * bafdiff_frac)
    notnan_flags = ~np.isnan(pred_bafs)

    valid_flags = functools.reduce(np.logical_and, [depthfit_valid_flags, baffit_valid_flags, notnan_flags])
    valid_flags_flat = valid_flags.reshape(valid_flags.shape[0], -1).all(axis=0)

    baf_errors_hdcp = baf_errors + CNt1 * CNt_handicap_factor
    baf_errors_flat = baf_errors_hdcp.reshape(baf_errors_hdcp.shape[0], -1).sum(axis=0)
    baf_errors_flat[~valid_flags_flat] = np.nan
    argmin = np.nanargmin(baf_errors_flat)

    unraveled_idx = np.unravel_index(argmin, CNt0.shape[1:])

    return {
        'CNt0': CNt0[0, :][unraveled_idx],
        'onecopy_depth': dd[0, :][unraveled_idx],
        'cellularity': cellularity[0, :][unraveled_idx],
        'Bt': Bt[np.index_exp[:,] + unraveled_idx],
        'CNt1': CNt1[np.index_exp[:,] + unraveled_idx],
        'predicted_baf': pred_bafs[np.index_exp[:,] + unraveled_idx],
    }


def show_solution_from_segments(depths, bafs, solution):
    nseg = len(depths)

    fig, axd = plt.subplot_mosaic(
        [['CN'], ['baf'], ['depth']], 
        gridspec_kw=dict(hspace=0.4),
        figsize=(8, 10),
    )

    xs = np.arange(nseg)
    width = 1

    axd['CN'].set_xlim(-1, nseg)
    axd['CN'].hlines(solution['CNt1'], xs - width/2, xs + width/2, alpha=0.4, color='black')
    axd['CN'].hlines(solution['Bt'], xs - width/2, xs + width/2, alpha=0.4, color='blue')
    axd['CN'].set_title('CN')

    axd['baf'].set_xlim(-1, nseg)
    axd['baf'].set_ylim(-0.02, 0.52)
    axd['baf'].hlines(bafs, xs - width/2, xs + width/2, alpha=0.4)
    axd['baf'].hlines(solution['predicted_baf'], xs - width/2, xs + width/2, alpha=0.4, color='tab:red')
    axd['baf'].set_title('BAF')

    axd['depth'].set_xlim(-1, nseg)
    axd['depth'].set_ylim(0, 2)
    axd['depth'].hlines(depths, xs - width/2, xs + width/2, alpha=0.4)
    axd['depth'].set_title('depth')

    depth0 = depths.min()
    for y in depth0 + solution['onecopy_depth'] * np.arange(20):
        axd['depth'].axhline(y, color='tab:red', alpha=0.4, linestyle='dotted')

    return fig, axd


###########
# helpers #
###########

def get_deco_CNg_Bg(CNg_fillvalue, Bg_fillvalue):
    """To fill in values for non-assembled contigs (like GL...)"""
    def decorator(func):
        sig = inspect.signature(func)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()
            if 'CNg' in ba.arguments:
                old_val = ba.arguments['CNg']
                if (old_val is not None) and (not np.isscalar(old_val)):
                    new_val = np.asarray(old_val).copy()
                    new_val[np.isnan(new_val)] = CNg_fillvalue
                    ba.arguments['CNg'] = new_val
            if 'Bg' in ba.arguments:
                old_val = ba.arguments['Bg']
                if (old_val is not None) and (not np.isscalar(old_val)):
                    new_val = np.asarray(old_val).copy()
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

    # rint
#    diff_from_rint = CNt - np.floor(CNt)
#    midzone_selector = np.logical_and((diff_from_rint > 0.4), (diff_from_rint < 0.6))
#    CNt[midzone_selector] = np.floor(CNt[midzone_selector])
#    CNt[~midzone_selector] = np.rint(CNt[~midzone_selector])

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
def find_Kc(*, depth1, CNt1, CNg1, depth2, CNt2, CNg2):
    """Assume 'depth1' and 'depth2' are corrected so that position-dependent bias has been removed.

    Assert that cellularity must not be zero

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
        raise KCError(f'Invalid input value combination: cellularity denominator becomes zero')

    cellularity = cellularity_numer / cellularity_denom

    # K
    K_inverse = (1 / depth1) * (CNg1 + cellularity * (CNt1 - CNg1))
    if K_inverse == 0:
        raise KCError(f'Invalid input value combination: inverse of K becomes zero')
    K = 1 / K_inverse

    # sanity check
    if K <= 0:
        raise KCError(f'Invalid input value combination: K has a non-positive value')
    if not ((cellularity > 0) and (cellularity <= 1)):
        raise KCError(f'Invalid input value combination: cellularity is not within the valid range: (0, 1]')

    return K, cellularity


def find_Kc_withdelta(*, depth, CNt, CNg, delta_depth):
    """
    - delta_depth: depth difference corresponding to one copy difference
    - CNg values are assumed to be the same between locations used to calculate delta_depth
    - depth, CNt, CNg: values corresponding to one of the locations used to calculate delta_depth
    """
    depth2 = depth + delta_depth
    CNt2 = CNt + 1
    CNg2 = CNg
    return find_Kc(
        depth1=depth, 
        CNt1=CNt, 
        CNg1=CNg, 
        depth2=depth2, 
        CNt2=CNt2, 
        CNg2=CNg2,
    )


@get_deco_CNg_Bg(CNg_fillvalue=2, Bg_fillvalue=1)
def get_predicted_depth(CNt, K, cellularity, CNg=None):
    if cellularity == 1:
        return K * CNt
    else:
        return K * (CNg * (1 - cellularity) + CNt * cellularity)


@get_deco_CNg_Bg(CNg_fillvalue=2, Bg_fillvalue=1)
def get_predicted_baf(CNt, Bt, cellularity, CNg=None, Bg=None):
    """Where CNt == 0 and CNg == 0, result is np.nan"""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        if cellularity == 1:
            return Bt / CNt
        else:
            return (
                (Bt * cellularity + Bg * (1 - cellularity))
                / (CNt * cellularity + CNg * (1 - cellularity))
            )


def get_possible_bafs(*, CNt, cellularity, CNg=2, Bg=1):
    half_CNt = np.floor(CNt / 2)
    Bt_candidates = half_CNt - np.arange(half_CNt + 1)
    predicted_bafs = get_predicted_baf(
        CNt=CNt, 
        Bt=Bt_candidates, 
        CNg=CNg, 
        Bg=Bg, 
        cellularity=cellularity,
    )
    return predicted_bafs


def get_possible_bafs_v2(*, depth0, delta_depth, CNg=2, Bg=1):
    data = dict()
    for CNt in itertools.count(0):
        try:
            K, cellularity = find_Kc_withdelta(depth=depth0, CNt=CNt, CNg=CNg, delta_depth=delta_depth)
        except KCError:
            break
        else:
            data[CNt] = {
                'K': K, 
                'cellularity': cellularity,
                'possible_bafs': get_possible_bafs(CNt=CNt, cellularity=cellularity, CNg=CNg, Bg=Bg),
            }

    return data


def get_possible_bafs_v3(*, cellularity, num_CNt=10, CNg=2, Bg=1):
    data = dict()
    for CNt in range(num_CNt):
        data[CNt] = get_possible_bafs(CNt=CNt, cellularity=cellularity, CNg=CNg, Bg=Bg)

    return data


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
    




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

        cellularity = np.atleast_1d(ba.arguments['cellularity'])
        notna_vals = cellularity[~np.isnan(cellularity)]
        if not np.logical_and(
            notna_vals > 0,
            notna_vals <= 1,
        ).all():
        #if ba.arguments['cellularity'] == 0:
            raise Exception(f'"cellularity" must be within interval (0, 1]')

        return func(*ba.args, **ba.kwargs)

    return wrapper


def deco_return_item(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        def helper(x):
            try:
                return x.item()
            except ValueError:
                return x

        result = func(*args, **kwargs)
        if isinstance(result, np.ndarray):
            return helper(result)
        elif isinstance(result, (tuple, list)):
            return tuple(helper(x) for x in result)
        else:
            raise Exception(f'Return value of decorated function must be either an np.ndarray or tuple or list.')

    return wrapper


def deco_return_squeezed(func):
    sig = inspect.signature(func)
    assert 'squeeze' not in sig.parameters.keys()

    @functools.wraps(func)
    def wrapper(*args, squeeze=True, **kwargs):
        def helper(x):
            result = np.squeeze(x)
            if result.ndim == 0:
                result = result.item()
            return result

        result = func(*args, **kwargs)
        if squeeze:
            if isinstance(result, np.ndarray):
                return helper(result)
            elif isinstance(result, (tuple, list)):
                return tuple(helper(x) for x in result)
            else:
                raise Exception(f'Return value of decorated function must be either an np.ndarray or tuple or list.')
        else:
            return result

    return wrapper


################
# miscellanous #
################

def check_haploid(is_female, chrom):
    return (
        (not is_female) 
        and refgenome.check_xy_chrom(chrom)
    )


def get_average_ploidy(CNt, lengths, axis=None):
    result = np.average(CNt, axis=axis, weights=lengths)
    return result


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
        depth=cnv_seg_gdf_copy.norm_depth_mean, 
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


def solution_from_segments_arghandler(
    depths, 
    bafs, 
    weights,

    all_depths,
    all_bafs,
    all_weights,
):
    def helper(data):
        if data is None:
            return None
        else:
            result = np.atleast_1d(
                np.squeeze(
                    np.asarray(data)
                )
            )
            assert result.ndim == 1
            return result

    if weights is not None:
        weights = np.asarray(weights)
        assert weights.shape == depths.shape

    depths = helper(depths)
    bafs = helper(bafs)
    weights = helper(weights)
    if weights is not None:
        assert weights.shape == depths.shape

    all_depths = helper(all_depths)
    all_bafs = helper(all_bafs)
    all_weights = helper(all_weights)
    if all_weights is not None:
        assert all_weights.shape == all_depths.shape

    return depths, bafs, weights, all_depths, all_bafs, all_weights


def solution_from_segments_common_params(
    depths,
    num_ndiv_cand,
    num_CNt0_cand,
    ndivs_offset_step,
    depthfit_cutoff,
    CNg=2,
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
    ndivs = np.tile(ndivs[np.newaxis, :], (CNt0.shape[0], 1))

    # make K, c
    depth0 = depths.min()
    K, cellularity = find_Kc(
        depth1=depth0, 
        CNt1=CNt0, 
        CNg1=CNg, 
        depth2=(depth0 + dd), 
        CNt2=(CNt0 + 1), 
        CNg2=CNg,
        raise_with_nan=False,
        squeeze=False,
    )

    # expand along a new axis which represent different segments
    #CNt0 = np.expand_dims(CNt0, 0)
    #dd = np.expand_dims(dd, 0)
    #ndivs = np.expand_dims(ndivs, 0)
    #cellularity = np.expand_dims(cellularity, 0)
    #K = np.expand_dims(K, 0)

    return CNt0, dd, ndivs, cellularity, K, depth0


def solution_from_segments_fit_to_data_new(
    depths, 
    bafs, 

    CNg,
    Bg,

    cellularity, 
    K,

    weights=None,
    depthfit_cutoff=None,
    baffit_cutoff=None,
):
    """Replacement of 'baf_targetval' with 0 where nan is inevitable due to condition of genomic position, not due to combination of cellularity and K:
        1) baf is unavailable (e.g. male X) => Bt is nan
        2) CNg == 0 => pred_bafs is nan
        3) CNt == 0 & cellularity == 1 => pred_bafs is nan
    """
    assert depths.ndim == 1
    assert bafs.ndim == 1
    assert CNg.ndim == 1
    assert Bg.ndim == 1

    # depth
    CNt1 = find_clonal_CNt(depths, cellularity, K, CNg=CNg, squeeze=False)
    pred_depths = get_predicted_depth(CNt1, K, cellularity, CNg=CNg, squeeze=False)
    depths_expand = np.expand_dims(
        depths, tuple(np.arange(pred_depths.ndim - 1) + 1),
    )
    depthdiffs = np.abs(depths_expand - pred_depths)
    onecopy_depthdiffs = get_onecopy_depthdiff(K, cellularity, squeeze=False)
    depthdiff_fractions = depthdiffs / onecopy_depthdiffs

    # baf
    Bt = find_clonal_Bt(bafs, CNt1, cellularity, CNg=CNg, Bg=Bg, squeeze=False)
    pred_bafs = get_predicted_baf(CNt1, Bt, cellularity, CNg=CNg, Bg=Bg, squeeze=False)
    bafs_expand = np.expand_dims(
        bafs, tuple(np.arange(pred_bafs.ndim - 1) + 1),
    )
    bafdiffs = np.abs(bafs_expand - pred_bafs)
    onecopy_bafdiffs = get_onecopy_bafdiff(CNt1, cellularity, CNg=CNg, squeeze=False)
    bafdiff_fractions = bafdiffs / onecopy_bafdiffs

    # get valid selectors
    if (depthfit_cutoff is None) and (baffit_cutoff is None):
        valid_flags = None
    else:
        depthfit_valid_flags = (
            np.full(True, depthdiff_fractions.shape)
            if depthfit_cutoff is None else
            depthdiff_fractions < depthfit_cutoff
        )
        baffit_valid_flags = (
            np.full(True, bafdiff_fractions.shape)
            if baffit_cutoff is None else
            bafdiff_fractions < baffit_cutoff
        )
        valid_flags = np.logical_and(depthfit_valid_flags, baffit_valid_flags)
        valid_flags = valid_flags.all(axis=0)

    # depth target function
    depth_targetval = depthdiff_fractions.copy()

    # baf target function
    baf_targetval = bafdiff_fractions.copy()

    # replacement of inevitable nan baf_targetval
    added_dims = tuple(np.arange(CNt1.ndim - 1) + 1)
    flag1 = np.expand_dims(np.isnan(bafs), added_dims)
    flag2 = np.expand_dims((CNg == 0), added_dims)
    flag3 = np.logical_and(
        (CNt1 == 0), 
        np.broadcast_to((cellularity == 1), CNt1.shape),
    )
    total_flag = functools.reduce(
        np.logical_or, np.broadcast_arrays(flag1, flag2, flag3),
    )
    baf_targetval[total_flag] = 0

    # total target function
    targetval = (
        depth_targetval 
        + baf_targetval 
    )

    def process_targetval(targetval, valid_flags, weights):
        if weights is not None:
            weights_expand = np.expand_dims(
                weights, tuple(np.arange(targetval.ndim - 1) + 1),
            )
            targetval = targetval * weights_expand
        targetval_sum = targetval.sum(axis=0)
        if valid_flags is not None:
            targetval_sum[~valid_flags] = np.nan

        #notna_vals = targetval_sum[~np.isnan(targetval_sum)]
        #assert len(notna_vals) > 0
        #mean = np.mean(notna_vals)
        #std = np.std(notna_vals)
        #norm_targetval_sum = (targetval_sum - mean) / std

        return targetval, targetval_sum

    #depth_targetval, depth_targetval_sum, norm_depth_targetval_sum = process_targetval(depth_targetval, valid_flags, weights)
    #baf_targetval, baf_targetval_sum, norm_baf_targetval_sum = process_targetval(baf_targetval, valid_flags, weights)
    #targetval, targetval_sum, norm_targetval_sum = process_targetval(targetval, valid_flags, weights)

    depth_targetval, depth_targetval_sum = process_targetval(depth_targetval, valid_flags, weights)
    baf_targetval, baf_targetval_sum = process_targetval(baf_targetval, valid_flags, weights)
    targetval, targetval_sum = process_targetval(targetval, valid_flags, weights)

    return (
        CNt1, 
        Bt, 
        depth_targetval_sum,
        baf_targetval_sum,
        targetval_sum,
    )


def solution_from_segments_fit_to_data(
    depth0,

    depths, 
    bafs, 
    weights,

    depthfit_cutoff,
    baffit_cutoff,
    CNt_handicap_factor,

    CNg,
    Bg,

    CNt0, 
    dd, 
    cellularity, 

    valid_flags=None,
):
    # expand depths and bafs
    depths = np.expand_dims(depths, (1, 2))
    bafs = np.expand_dims(bafs, (1, 2))

    ### from here ndims of CNt0, dd, cellularity, depths, bafs are all same ###

    # make CNt1 for each segment
    dd_quotients = (depths - depth0) / dd
    dd_quotients[dd_quotients < 0] = 0
    CNt_diffs = np.rint(dd_quotients)
    depthdiff_fractions = np.abs(dd_quotients - CNt_diffs)
    CNt1 = CNt0 + CNt_diffs

    # determine Bt for each segment (axis for Bt candidates is the last one)
    half_CNt1 = (np.floor(CNt1 / 2)).astype(int)
    Bt_offsets = np.arange(half_CNt1.max() + 1)
    Bt_candidates = (np.expand_dims(half_CNt1, -1) - Bt_offsets).astype(float)
    Bt_candidates[Bt_candidates < 0] = np.nan

    # choose Bt candidate indexes with the least difference
    pred_bafs_candidates = get_predicted_baf(
        CNt=np.expand_dims(CNt1, -1), 
        Bt=Bt_candidates, 
        CNg=CNg, 
        Bg=Bg, 
        cellularity=np.expand_dims(cellularity, -1),
        squeeze=False,
    )
    bafdiff_candidates = np.abs(pred_bafs_candidates - np.expand_dims(bafs, -1))
    argmins = tools.nanargmin(bafdiff_candidates, axis=-1, keepdims=True)

    bafdiffs = np.squeeze(np.take_along_axis(bafdiff_candidates, argmins, axis=-1), axis=-1)
    Bt = np.squeeze(np.take_along_axis(Bt_candidates, argmins, axis=-1), axis=-1)
    pred_bafs = np.squeeze(np.take_along_axis(pred_bafs_candidates, argmins, axis=-1), axis=-1)

    # make onecopy bafdiff
    onecopy_bafdiffs = get_onecopy_bafdiff(
        CNt=CNt1, 
        cellularity=cellularity,
        CNg=CNg, 
        squeeze=False,
    )
    bafdiff_fractions = bafdiffs / onecopy_bafdiffs

    # get valid selectors
    if valid_flags is None:
        depthfit_valid_flags = depthdiff_fractions < depthfit_cutoff
        baffit_valid_flags = bafdiff_fractions < baffit_cutoff
        notnan_flags = ~np.isnan(pred_bafs)

        valid_flags = functools.reduce(
            np.logical_and, 
            [depthfit_valid_flags, baffit_valid_flags, notnan_flags],
        )
        valid_flags = valid_flags.all(axis=0)

    # make target function
    depth_targetval = depthdiff_fractions.copy()

    # for regions where baf is unavailable (e.g. male X)
    baf_targetval = bafdiff_fractions.copy()
    baf_isnan = np.isnan(np.squeeze(bafs, axis=(1, 2)))
    baf_targetval[baf_isnan] = 0

    targetval = (
        depth_targetval 
        + baf_targetval 
        #+ (CNt1 * CNt_handicap_factor)
    )

    def process_targetval(targetval, valid_flags, weights):
        if weights is not None:
            targetval = targetval * np.expand_dims(weights, (1, 2))
        targetval_sum = targetval.sum(axis=0)
        targetval_sum[~valid_flags] = np.nan

        notna_vals = targetval_sum[~np.isnan(targetval_sum)]
        #assert len(notna_vals) > 0
        mean = np.mean(notna_vals)
        std = np.std(notna_vals)
        norm_targetval_sum = (targetval_sum - mean) / std

        return targetval, targetval_sum, norm_targetval_sum

    depth_targetval, depth_targetval_sum, norm_depth_targetval_sum = process_targetval(depth_targetval, valid_flags, weights)
    baf_targetval, baf_targetval_sum, norm_baf_targetval_sum = process_targetval(baf_targetval, valid_flags, weights)
    targetval, targetval_sum, norm_targetval_sum = process_targetval(targetval, valid_flags, weights)

    return (
        CNt1, Bt, valid_flags,
        norm_depth_targetval_sum,
        norm_baf_targetval_sum,
        norm_targetval_sum,
        pred_bafs,
    )


def solution_from_segments(
    depths, 
    bafs, 
    weights=None,

    all_depths=None,
    all_bafs=None,
    all_weights=None,

    num_ndiv_cand=20,
    num_CNt0_cand=10,
    ndivs_offset_step=0.01,

    depthfit_cutoff=0.5,
    baffit_cutoff=0.5,
    CNt_handicap_factor=0,

    CNg=2,
    Bg=1,
):
    depths, bafs, weights, all_depths, all_bafs, all_weights = solution_from_segments_arghandler(
        depths, 
        bafs, 
        weights,

        all_depths,
        all_bafs,
        all_weights,
    )

    CNt0, dd, ndivs, cellularity, K, depth0 = solution_from_segments_common_params(
        depths,
        num_ndiv_cand,
        num_CNt0_cand,
        ndivs_offset_step,
        depthfit_cutoff,
    )

    (
        CNt1, Bt, valid_flags,
        depth_targetval_sum,
        baf_targetval_sum,
        targetval_sum,
        pred_bafs,
    ) = solution_from_segments_fit_to_data(
        depth0,

        depths, 
        bafs, 
        weights,

        depthfit_cutoff,
        baffit_cutoff,
        CNt_handicap_factor,

        CNg,
        Bg,

        CNt0, 
        dd, 
        cellularity, 

        valid_flags=None,
    )

    argmin = np.nanargmin(targetval_sum)
    unraveled_idx = np.unravel_index(argmin, shape=targetval_sum.shape)

    if all_depths is not None:
        (
            all_CNt1, all_Bt, _,
            all_depth_targetval_sum,
            all_baf_targetval_sum,
            all_targetval_sum,
            _,
        ) = solution_from_segments_fit_to_data(
            depth0,

            all_depths, 
            all_bafs, 
            all_weights,

            depthfit_cutoff,
            baffit_cutoff,
            CNt_handicap_factor,

            CNg,
            Bg,

            CNt0, 
            dd, 
            cellularity, 

            valid_flags=valid_flags,
        )
    else:
        all_CNt1 = None 
        all_Bt = None 
        all_depth_targetval_sum = None
        all_baf_targetval_sum = None
        all_targetval_sum = None


    #targetval_flat = targetval.reshape(targetval.shape[0], -1).sum(axis=0)
    #baf_errors_hdcp = bafdiffs + CNt1 * CNt_handicap_factor
    #baf_errors_flat = baf_errors_hdcp.reshape(baf_errors_hdcp.shape[0], -1).sum(axis=0)
    #baf_errors_flat[~valid_flags_flat] = np.nan
    #argmin = np.nanargmin(baf_errors_flat)

#    return {
#        'CNt0': CNt0[0, :][unraveled_idx],
#        'onecopy_depth': dd[0, :][unraveled_idx],
#        'cellularity': cellularity[0, :][unraveled_idx],
#        'K': K[0, :][unraveled_idx],
#        'CNt1': CNt1[np.index_exp[:,] + unraveled_idx],
#        'Bt': Bt[np.index_exp[:,] + unraveled_idx],
#        'predicted_baf': pred_bafs[np.index_exp[:,] + unraveled_idx],
#    }

    return {
        'CNt0': CNt0[0, :],
        'onecopy_depth': dd[0, :],
        'num_division': ndivs[0, :],
        'cellularity': cellularity[0, :],
        'K': K[0, :],

        'targetval': targetval_sum,
        'depth_targetval': depth_targetval_sum,
        'baf_targetval': baf_targetval_sum,
        'CNt': CNt1,
        'Bt': Bt,

        'argmin': unraveled_idx,
        'predicted_baf': pred_bafs,

        'all_targetval': all_targetval_sum,
        'all_depth_targetval': all_depth_targetval_sum,
        'all_baf_targetval': all_baf_targetval_sum,
        'all_CNt': all_CNt1,
        'all_Bt': all_Bt,
    }


@deco.get_deco_atleast1d([
    'picked_depths', 
    'depths', 
    'bafs', 
    'lengths',
    'target_depths',
    'target_bafs',
    'target_weights',
])
def solution_from_depth_limits(
    picked_depths,
    depths, 
    bafs, 
    CNg,
    Bg,
    lengths,
    #weights=None,

    target_depths=None,
    target_bafs=None,
    target_weights=None,

    num_ndiv_cand=20,
    num_CNt0_cand=10,
    ndivs_offset_step=0.01,

    depthfit_cutoff=0.3,
    baffit_cutoff=0.3,

    major_CNg=2,
):
    CNt0, dd, ndivs, cellularity, K, _ = solution_from_segments_common_params(
        picked_depths,
        num_ndiv_cand,
        num_CNt0_cand,
        ndivs_offset_step,
        depthfit_cutoff,
        CNg=major_CNg,
    )

    (
        CNt1, 
        Bt, 
        depth_targetval_sum,
        baf_targetval_sum,
        targetval_sum,
    ) = solution_from_segments_fit_to_data_new(
        depths, 
        bafs, 

        CNg,
        Bg,

        cellularity, 
        K,

        weights=lengths,
        depthfit_cutoff=None,
        baffit_cutoff=None,
    )

    argsort = np.argsort(targetval_sum, axis=None)
    unraveled_indexes = [np.unravel_index(x, shape=targetval_sum.shape) for x in argsort]

    #argmin = np.nanargmin(targetval_sum)
    #unraveled_idx = np.unravel_index(argmin, shape=targetval_sum.shape)
    ploidies = get_average_ploidy(CNt1, lengths, axis=0)

    if target_depths is not None:
        (
            target_CNt1, 
            target_Bt, 
            target_depth_targetval_sum,
            target_baf_targetval_sum,
            target_targetval_sum,
        ) = solution_from_segments_fit_to_data_new(
            depths=target_depths, 
            bafs=target_bafs, 

            CNg=CNg,
            Bg=Bg,

            cellularity=cellularity, 
            K=K,

            weights=target_weights,
            depthfit_cutoff=depthfit_cutoff,
            baffit_cutoff=baffit_cutoff,
        )
    else:
        target_CNt1 = None 
        target_Bt = None 
        target_depth_targetval_sum = None
        target_baf_targetval_sum = None
        target_targetval_sum = None

    return {
        'CNt0': CNt0,
        'onecopy_depth': dd,
        'num_division': ndivs,
        'cellularity': cellularity,
        'K': K,

        'CNt': CNt1,
        'Bt': Bt,

        'targetval': targetval_sum,
        'depth_targetval': depth_targetval_sum,
        'baf_targetval': baf_targetval_sum,

        #'argmin': unraveled_idx,
        'sorted_indexes': unraveled_indexes,
        'ploidy': ploidies,

        'target_targetval': target_targetval_sum,
        'target_depth_targetval': target_depth_targetval_sum,
        'target_baf_targetval': target_baf_targetval_sum,
        'target_CNt': target_CNt1,
        'target_Bt': target_Bt,
    }


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


#@get_deco_CNg_Bg(CNg_fillvalue=2, Bg_fillvalue=1)
#@deco_cellularity_sanitycheck
#def find_clonal_CNt_nonvectorized(corrected_depth, cellularity, K, CNg=None):
#    """Assume 'corrected_depth' is corrected so that position-dependent bias has been removed.
#    It can be formulated like:
#        'corrected_depth' = K * (CNg * (1 - cellularity) + CNt * cellularity)
#
#    derivation:
#        corrected_depth = K * (CNg * (1 - cellularity) + CNt * cellularity)
#        corrected_depth / K = CNg * (1 - cellularity) + CNt * cellularity
#        (corrected_depth / K) - (CNg * (1 - cellularity)) = CNt * cellularity
#        CNt = ( (corrected_depth / K) - (CNg * (1 - cellularity)) ) / cellularity
#    """
#    CNg_given = (CNg is not None)
#    if (cellularity != 1) and (not CNg_given):
#        raise Exception(f'When cellularity is not 1, "CNg" must be given.')
#
#    # main
#    if cellularity == 1:
#        CNt = (corrected_depth / K)
#    else:
#        CNt = (
#            (corrected_depth / K) 
#            - (CNg * (1 - cellularity))
#        ) / cellularity
#
#    # rint
##    diff_from_rint = CNt - np.floor(CNt)
##    midzone_selector = np.logical_and((diff_from_rint > 0.4), (diff_from_rint < 0.6))
##    CNt[midzone_selector] = np.floor(CNt[midzone_selector])
##    CNt[~midzone_selector] = np.rint(CNt[~midzone_selector])
#
#    CNt = np.clip(np.rint(CNt), 0, None)
#    if CNg_given:
#        CNt[CNg == 0] = 0
#
#    return CNt


#@deco_return_item
@deco_return_squeezed
@deco.get_deco_atleast1d(['depth', 'cellularity', 'K', 'CNg'], keep_none=True)
@get_deco_CNg_Bg(CNg_fillvalue=2, Bg_fillvalue=1)
@deco_cellularity_sanitycheck
def find_clonal_CNt(depth, cellularity, K, CNg=None):
    """Assume 'depth' is corrected so that position-dependent bias has been removed.
    It can be formulated like:
        'depth' = K * (CNg * (1 - cellularity) + CNt * cellularity)

    derivation:
        depth = K * (CNg * (1 - cellularity) + CNt * cellularity)
        depth / K = CNg * (1 - cellularity) + CNt * cellularity
        (depth / K) - (CNg * (1 - cellularity)) = CNt * cellularity
        CNt = ( (depth / K) - (CNg * (1 - cellularity)) ) / cellularity

    Args:
        depth: ndim must be 1.

    Axis of parameters:
        Axis 0: Different genomic positions. The axis of "depth" and "CNg" must be this.
        Axis 1, 2, ...: Candidates of cellularity and K.

        "cellularity", "K" must be broadcastable to each other.
    """
    # sanitycheck
    assert depth.ndim == 1

    CNg_given = (CNg is not None)
    if CNg_given:
        assert CNg.ndim == 1
    if (cellularity != 1).any() and (not CNg_given):
        raise Exception(f'When any cellularity value is not 1, "CNg" must be given.')

    # non-position shape
    nonpos_shape = np.broadcast_shapes(cellularity.shape, K.shape)
    added_dims = tuple(np.arange(len(nonpos_shape)) + 1)

    # main
    depth_expand = np.expand_dims(depth, added_dims)
    if CNg_given:
        CNg_expand = np.expand_dims(CNg, added_dims)

    if (cellularity == 1).all():
        CNt = depth_expand / K
    else:
        CNt = (
            (depth_expand / K)
            - CNg_expand * (1 - cellularity)
        ) / cellularity

    CNt = np.clip(np.rint(CNt), 0, None)
    if CNg_given:
        CNt[np.broadcast_to(CNg_expand == 0, CNt.shape)] = 0

    return CNt


#@get_deco_CNg_Bg(CNg_fillvalue=2, Bg_fillvalue=1)
#@deco_cellularity_sanitycheck
#def find_clonal_Bt_nonvectorized(baf, CNt, cellularity, CNg=None, Bg=None):
#    """
#    derivation:
#        baf = (
#            (Bt * cellularity + Bg * (1 - cellularity))
#            / (CNt * cellularity + CNg * (1 - cellularity))
#        )
#        Bt * cellularity + Bg * (1 - cellularity) = (
#            baf * (CNt * cellularity + CNg * (1 - cellularity))
#        )
#        Bt * cellularity = (
#            (baf * (CNt * cellularity + CNg * (1 - cellularity)))
#            - (Bg * (1 - cellularity))
#        )
#        Bt = (
#            (baf * (CNt * cellularity + CNg * (1 - cellularity)))
#            - (Bg * (1 - cellularity))
#        ) / cellularity
#        Bt = (
#            (baf * (CNt + CNg * (1 / cellularity - 1)))
#            - (Bg * (1 / cellularity - 1))
#        )
#    """
#    # sanitycheck
#    CNg_given = (CNg is not None)
#    Bg_given = (Bg is not None)
#    if (
#        (cellularity != 1) 
#        and (
#            (not CNg_given) or (not Bg_given)
#        )
#    ):
#        raise Exception(f'When cellularity is not 1, "CNg" and "Bg" must be given.')
#
#    # main
#    c = (1 / cellularity) - 1
#    if c == 0:
#        Bt = baf * CNt 
#    else:
#        Bt = (baf * (CNt + CNg * c)) - (Bg * c)
#
#    # correction
#    Bt = np.clip(np.rint(Bt), 0, np.floor(CNt / 2))
#    if CNg_given or Bg_given:
#        if CNg_given:
#            zero_selector_CNg = (CNg == 0)
#        else:
#            zero_selector_CNg = np.repeat(False, len(CNg))
#
#        if Bg_given:
#            zero_selector_Bg = (Bg == 0)
#        else:
#            zero_selector_Bg = np.repeat(False, len(Bg))
#
#        zero_selector = np.logical_or(zero_selector_CNg, zero_selector_Bg)
#        Bt[zero_selector] = 0
#
#    return Bt


#@deco_return_item
@deco_return_squeezed
@deco.get_deco_atleast1d(['baf', 'CNt', 'cellularity', 'CNg', 'Bg'], keep_none=True)
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

    Args:
        baf: ndim must be 1

    Axis of parameters:
        Axis 0: Different genomic positions. 
            The only axis of "baf" must be this.
            Axis 0 of CNt must correspond to this.
            Axis 0 of CNt and baf must be broadcastable to each other.
        Axis 1, 2, ...: Must represent candidates of cellularity, CNg, and Bg.

        "cellularity", "CNg", and "Bg" must be broadcastable to each other.
    """
    # sanitycheck
    assert baf.ndim == 1
    assert baf.shape[0] == CNt.shape[0]

    CNg_given = (CNg is not None)
    if CNg_given:
        assert CNg.ndim == 1
    Bg_given = (Bg is not None)
    if Bg_given:
        assert Bg.ndim == 1

    assert CNg_given == Bg_given
    ginfo_given = CNg_given

    c_all_full = (cellularity == 1).all()
    if (not c_all_full) and (not ginfo_given):
        raise Exception(f'When any of cellularity is not 1, "CNg" and "Bg" must be given.')

    # non-position shape
    if CNt.ndim == 1:
        nonpos_shape = (1,)
        CNt = np.expand_dims(CNt, 1)
    else:
        nonpos_shape = CNt.shape[1:]

    assert nonpos_shape == np.broadcast_shapes(nonpos_shape, cellularity.shape)
    added_dims = tuple(np.arange(len(nonpos_shape)) + 1)

    # main
    baf_expand = np.expand_dims(baf, added_dims)
    c_prime = (1 / cellularity) - 1
    if c_all_full:
        Bt = baf_expand * CNt 
    else:
        CNg_expand = np.expand_dims(CNg, added_dims)
        Bg_expand = np.expand_dims(Bg, added_dims)
        Bt = (baf_expand * (CNt + CNg_expand * c_prime)) - (Bg_expand * c_prime)

    # correction
    Bt = np.clip(np.rint(Bt), 0, np.floor(CNt / 2))
    if ginfo_given:
        zero_selector_CNg = (CNg_expand == 0)
        zero_selector_Bg = (Bg_expand == 0)
        zero_selector = np.logical_or(
            np.broadcast_to(zero_selector_CNg, Bt.shape), 
            np.broadcast_to(zero_selector_Bg, Bt.shape), 
        )
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


#@deco_return_item
@deco_return_squeezed
@get_deco_CNg_Bg(CNg_fillvalue=2, Bg_fillvalue=1)
@deco.get_deco_atleast1d(['depth1', 'CNt1', 'CNg1', 'depth2', 'CNt2', 'CNg2'])
def find_Kc(*, depth1, CNt1, CNg1, depth2, CNt2, CNg2, raise_with_nan=False):
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
    if (cellularity_denom == 0).any():
        raise KCError(f'Invalid input value combination: cellularity denominator becomes zero')
    cellularity = cellularity_numer / cellularity_denom

    # cellularity sanitycheck
    invalid_selector = np.logical_or((cellularity <= 0), (cellularity > 1))
    if invalid_selector.any() and raise_with_nan:
        raise KCError(f'Invalid input value combination: cellularity is out of range (0, 1]')
    cellularity[invalid_selector] = np.nan

    # K
    K_inverse = (1 / depth1) * (CNg1 + cellularity * (CNt1 - CNg1))
        # since "cellularity" is at least 1d, K_inverse is also at least 1d
    if (K_inverse == 0).any():
        raise KCError(f'Invalid input value combination: inverse of K becomes zero')
    K = 1 / K_inverse

    # K sanitycheck
    invalid_K_selector = (K <= 0)
    if invalid_K_selector.any() and raise_with_nan:
        raise KCError(f'Invalid input value combination: K <= 0')
    K[invalid_K_selector] = np.nan
        
    # result
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


#@deco_return_item
@deco_return_squeezed
@get_deco_CNg_Bg(CNg_fillvalue=2, Bg_fillvalue=1)
@deco.get_deco_atleast1d(['CNt', 'K', 'cellularity', 'CNg'])
def get_predicted_depth(CNt, K, cellularity, CNg=None):
    c_all_full = (cellularity == 1).all()
    CNg_given = (CNg is not None)

    if (not c_all_full) and (not CNg_given):
        raise Exception(f'When any of cellularity is not 1, CNg must be given.')

    # non-position shape
    if CNt.ndim == 1:
        nonpos_shape = (1,)
        CNt = np.expand_dims(CNt, 1)
    else:
        nonpos_shape = CNt.shape[1:]

    assert nonpos_shape == np.broadcast_shapes(nonpos_shape, cellularity.shape)
    assert nonpos_shape == np.broadcast_shapes(nonpos_shape, K.shape)
    added_dims = tuple(np.arange(len(nonpos_shape)) + 1)

    # main
    if c_all_full:
        result = K * CNt
    else:
        CNg_expand = np.expand_dims(CNg, added_dims)
        result = K * (CNg_expand * (1 - cellularity) + CNt * cellularity)

    return result


#@deco_return_item
@deco_return_squeezed
@get_deco_CNg_Bg(CNg_fillvalue=2, Bg_fillvalue=1)
@deco.get_deco_atleast1d(['K', 'cellularity'])
def get_onecopy_depthdiff(K, cellularity):
    result = K * cellularity
    return result


#@deco_return_item
@deco_return_squeezed
@get_deco_CNg_Bg(CNg_fillvalue=2, Bg_fillvalue=1)
@deco.get_deco_atleast1d(['CNt', 'Bt', 'cellularity', 'CNg', 'Bg'])
def get_predicted_baf(CNt, Bt, cellularity, CNg=None, Bg=None):
    """Where CNg == 0 or (cellularity == 1 and CNt == 0), result is np.nan"""
 
    # sanitycheck
    assert CNt.shape == Bt.shape

    c_all_full = (cellularity == 1).all()
    CNg_given = (CNg is not None)
    if CNg_given:
        assert CNg.ndim == 1
    Bg_given = (Bg is not None)
    if Bg_given:
        assert Bg.ndim == 1

    assert CNg_given == Bg_given
    ginfo_given = CNg_given

    if (not c_all_full) and (not ginfo_given):
        raise Exception(f'When any of cellularity is not 1, CNg and Bg must be given.')

    # non-position shape
    if CNt.ndim == 1:
        nonpos_shape = (1,)
        CNt = np.expand_dims(CNt, 1)
        Bt = np.expand_dims(Bt, 1)
    else:
        nonpos_shape = CNt.shape[1:]

    assert nonpos_shape == np.broadcast_shapes(nonpos_shape, cellularity.shape)
    added_dims = tuple(np.arange(len(nonpos_shape)) + 1)

    # main
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)  # in order to ignore zerodivisionwarning
        if c_all_full:
            result = Bt / CNt  # where CNg == 0, CNt ought to be 0
        else:
            CNg_expand = np.expand_dims(CNg, added_dims)
            Bg_expand = np.expand_dims(Bg, added_dims)
            one_minus_c = (1 - cellularity)
            result = (
                (Bt * cellularity + Bg_expand * one_minus_c)
                / (CNt * cellularity + CNg_expand * one_minus_c)
            )

    return result


#@deco_return_item
@deco_return_squeezed
@get_deco_CNg_Bg(CNg_fillvalue=2, Bg_fillvalue=1)
@deco.get_deco_atleast1d(['CNt', 'cellularity', 'CNg'])
def get_onecopy_bafdiff(CNt, cellularity, CNg=None):
    """Where CNg == 0 or (cellularity == 1 and CNt == 0), result is np.nan"""

    CNg_given = (CNg is not None)
    if CNg_given:
        assert CNg.ndim == 1
    c_all_full = (cellularity == 1).all()
    if (not c_all_full) and (not CNg_given):
        raise Exception(f'When any cellularity value is not 1, "CNg" must be given.')

    # non-position shape
    if CNt.ndim == 1:
        nonpos_shape = (1,)
        CNt = np.expand_dims(CNt, 1)
    else:
        nonpos_shape = CNt.shape[1:]
    assert nonpos_shape == np.broadcast_shapes(nonpos_shape, cellularity.shape)
    added_dims = tuple(np.arange(len(nonpos_shape)) + 1)

    # main
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        if c_all_full:
            result = 1 / CNt
        else:
            CNg_expand = np.expand_dims(CNg, added_dims)
            result = (
                cellularity
                / (CNt * cellularity + CNg_expand * (1 - cellularity))
            )

    return result


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
    




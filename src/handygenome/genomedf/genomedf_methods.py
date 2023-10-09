import re
import os
import datetime
import inspect
import itertools
import contextlib
import multiprocessing
import functools

import pandas as pd
import numpy as np
import pyranges as pr
import scipy.stats

import handygenome.deco as deco
import handygenome.refgenome.refgenome as refgenome
import handygenome.tools as tools
import handygenome.logutils as logutils
import handygenome.cnv.rdnacopy as rdnacopy

import handygenome.genomedf.genomedf_utils as genomedf_utils


MERGETYPES_SIMPLEAGG = ['mean', 'std', 'first', 'last', 'max', 'min', 'median', 'skew']
MERGETYPES_WEIGHTED_SIMPLEAGG = ['weighted_mean']
MERGETYPES_BYLENGTH = ['longest', 'shortest']
RESERVED_COLNAMES = [
    'Distance',
    'Overlap',
    'Start_b', 
    'End_b',
]


#def split_left_gdf(gdf, width=1000):
#    partial_left_gdf_list = list()
#    for chrom, subgdf in gdf.group_bychrom(sort=True).items():
#        partial_left_gdf_list.extend(subgdf.equal_nrow_split(width=width))
#
#    return partial_left_gdf_list


def targetfunc(partial_left_gdf, right_gdf, join_kwargs):
    subset_args = (
        partial_left_gdf.chroms[0], 
        partial_left_gdf.start0s[0], 
        partial_left_gdf.end0s[-1], 
    )
    partial_right_gdf = right_gdf.subset(*subset_args)
    result = join_base(partial_left_gdf, partial_right_gdf, **join_kwargs)
    return result


###################


@deco.get_deco_arg_choices({'how': [None, 'left']})
def join_base(
    left_gdf, right_gdf, 
    right_gdf_cols=None,
    how='left', 
    merge=['mean', 'std'],
    ddof=0,
    overlapping_length=False,
    N_colname='N',
    suffixes=None,
    merge_args=tuple(),
    merge_kwargs=dict(),
    winsorize=None,
    omit_N=False,
    find_nearest=False,
):
    # merge arghandler
    merge, mergetype_groups = merge_arghandler(merge)
    # suffixes arghandler
    suffixes = suffixes_arghandler(suffixes, merge)
    # right_gdf-related argument handlers
    right_gdf_cols, cols_tobe_added, winsorize = sanitycheck_arghandler(
        left_gdf, right_gdf, right_gdf_cols, merge, N_colname, suffixes, winsorize, omit_N,
    )
        # cols_tobe_added: later will be used for merged_right column ordering

    # handle non-overlap case
    no_overlap = precheck_nonoverlap(left_gdf, right_gdf)
    if no_overlap:
        return make_nonoverlap_result(left_gdf, how, cols_tobe_added)

    # make parameters
    left_columns = left_gdf.columns
    report_overlap = (
        overlapping_length 
        and (merge is not None)
        and set(merge).intersection(MERGETYPES_BYLENGTH)
    )

    # run pyranges join
    left_gr = left_gdf.gr
    right_gr = right_gdf.gr[right_gdf_cols]
    joined_gr = left_gr.join(
        right_gr, 
        how=how, 
        report_overlap=report_overlap,
    )
    if joined_gr.empty:
        return make_nonoverlap_result(left_gdf, how, cols_tobe_added)

    joined_df = joined_gr.df
    joined_df.reset_index(drop=True, inplace=True)

    # handle unmatched rows
    if find_nearest:
        # this function modifies "joined_df" inplace
        add_nearest(joined_df, right_gdf, left_columns)

    # replace missing values with pd.NA
    unmatched_selector = (joined_df['Start_b'] == -1).to_numpy()
    joined_df.loc[
        unmatched_selector, 
        ~joined_df.columns.isin(left_columns)
    ] = pd.NA  # this becomes np.nan in float columns

    # merge
    if merge is None:
        result_df = joined_df
    else:
        # unmatched rows
        unmatched_df = joined_df.loc[unmatched_selector, left_columns]
        if not omit_N:
            unmatched_df.insert(unmatched_df.shape[1], N_colname, 0)

        # matched left
        matched_left = joined_df.loc[~unmatched_selector, left_columns]
        counts, groupkey = genomedf_utils.get_coord_groupkey(matched_left, left_gdf.chromdict)
            # "counts" is aligned with "groupkey"

        merged_left = matched_left.drop_duplicates()
        merged_left.reset_index(drop=True, inplace=True)
        if not omit_N:
            merged_left.insert(merged_left.shape[1], N_colname, counts)

        # matched right
        matched_right = joined_df.loc[
            ~unmatched_selector, 
            ~joined_df.columns.isin(left_columns)
        ]
        matched_right.reset_index(drop=True, inplace=True)
        if winsorize is not None:
            winsorize_matched_right(
                matched_right, 
                groupkey, 
                winsorize,
            )

        merged_right_list = [
            merge_right(
                matched_right, 
                groupkey, 
                right_gdf_cols, 
                mergetype, 
                suffixes, 
                ddof,
                merge_args, 
                merge_kwargs,
            )
            for mergetype in merge
        ]
        merged_right = pd.concat(merged_right_list, axis=1)

        # concat all
        matched = pd.concat([merged_left, merged_right], axis=1)
        result_df = pd.concat([unmatched_df, matched], axis=0)

        # coerce dtype
        if not omit_N:
            result_df = result_df.astype({N_colname: int})

    # return
    result = left_gdf.spawn(result_df)
    result.sort()
    return result


def precheck_nonoverlap(left_gdf, right_gdf):
    return (
        left_gdf.is_empty
        or right_gdf.is_empty
        or (not set(left_gdf.chroms).intersection(right_gdf.chroms))
        or (left_gdf.end0s.max() <= right_gdf.start0s.min())
        or (left_gdf.start0s.min() >= right_gdf.end0s.max())
    )


def make_nonoverlap_result(left_gdf, how, cols_tobe_added):
    if how is None:
        result = left_gdf.iloc[[], :]
    elif how == 'left':
        result = left_gdf.copy()

    for x in cols_tobe_added:
        result[x] = np.nan

    return result


def add_nearest(joined_df, right_gdf, left_columns):
    unmatched_selector = (joined_df['Start_b'] == -1).to_numpy()
    unmatched_df = joined_df.iloc[np.nonzero(unmatched_selector)[0], :3]
    unmatched_df_annot = pr.PyRanges(unmatched_df).nearest(right_gdf.gr).df
    annot_cols = unmatched_df_annot.columns[3:].to_list()  # includes Start_b, End_b, Distance

    joined_df['Distance'] = np.nan  # joined_df already includes columns "Start_b" and "End_b"
    joined_df.loc[unmatched_selector, annot_cols] = unmatched_df_annot.loc[:, annot_cols].to_numpy()


def group_mergetypes(merge):
    mergetype_groups = {
        'simpleagg': list(),
        'weighted_simpleagg': list(),
        'bylength': list(),
        'custom': list(),
    }
    for mergetype in merge:
        if mergetype in MERGETYPES_SIMPLEAGG:
            mergetype_groups['simpleagg'].append(mergetype)
        elif mergetype in MERGETYPES_WEIGHTED_SIMPLEAGG:
            mergetype_groups['weighted_simpleagg'].append(mergetype)
        elif mergetype in MERGETYPES_BYLENGTH:
            mergetype_groups['bylength'].append(mergetype)
        else:
            mergetype_groups['custom'].append(mergetype)

    return mergetype_groups


def merge_arghandler(merge):
    if merge is None:
        mergetype_groups = None
    else:
        # modification
        merge = np.atleast_1d(merge)
        assert len(merge) > 0
        mergetype_groups = group_mergetypes(merge)

        # sanitycheck
        if not all(map(callable, mergetype_groups['custom'])):
            raise Exception(
                f'Custom merge methods must be non-lambda callables'
            )
        if any(map(tools.check_is_lambda, mergetype_groups['custom'])):
            raise Exception(f'Custom merge methods must not include lambda')

    return merge, mergetype_groups


def get_suffixes_key(mergetype):
    if isinstance(mergetype, str):
        return mergetype
    elif callable(mergetype):
        return mergetype.__name__
    else:
        raise Exception(f'Unavailable mergetype')


def suffixes_arghandler(suffixes, merge):
    if merge is None:
        pass
    else:
        # make default suffixes
        default_suffixes_keys = list()
        default_suffixes_vals = list()

        for mergetype in merge:
            key = get_suffixes_key(mergetype)
            if isinstance(mergetype, str):
                val = f'_{mergetype}'
            elif callable(mergetype):
                val = f'_{mergetype.__name__}'

            default_suffixes_keys.append(key)
            default_suffixes_vals.append(val)

        if not tools.check_unique(default_suffixes_keys):
            raise Exception(f'Overlapping merge method names')
        if not tools.check_unique(default_suffixes_vals):
            raise Exception(f'Overlapping merge method suffixes')

        default_suffixes = dict(
            zip(default_suffixes_keys, default_suffixes_vals)
        )

        # apply custom suffixes
        if suffixes is None:
            suffixes = default_suffixes
        else:
            assert isinstance(suffixes, dict)
            suffixes = default_suffixes | suffixes
            if not tools.check_unique(suffixes.values()):
                raise Exception(f'Overlapping merge method suffixes')

    return suffixes


def sanitycheck_arghandler(
    left_gdf, right_gdf, right_gdf_cols, merge, N_colname, suffixes, winsorize, omit_N,
):
    ########################################
    # sanitycheck BEFORE args modification #
    ########################################
    assert len(left_gdf.columns) == len(set(left_gdf.columns))
    assert len(right_gdf.columns) == len(set(right_gdf.columns))

    assert not left_gdf.check_duplicate_coords()
    assert not right_gdf.check_duplicate_coords()

    #####################
    # args modification #
    #####################

    # right_gdf_cols is a list
    if right_gdf_cols is None:
        right_gdf_cols = list(right_gdf.annot_cols)
    else:
        right_gdf_cols = list(np.atleast_1d(right_gdf_cols))

    if winsorize is not None:
        winsorize = tuple(winsorize)
        assert len(winsorize) == 2

    #######################################
    # sanitycheck AFTER args modification #
    #######################################

    assert len(right_gdf_cols) == len(set(right_gdf_cols))
    assert set(right_gdf_cols).issubset(right_gdf.annot_cols)
    assert not set(left_gdf.annot_cols).intersection(right_gdf_cols), (
        f'There must be no overlapping annotation columns between '
        f'two {left_gdf.__class__.__name__} objects.'
    )

    reserved_colnames = list(RESERVED_COLNAMES)
    if not omit_N:
        reserved_colnames.append(N_colname)

    cols_after_pyranges_join = list(left_gdf.columns) + right_gdf_cols
    assert not set(reserved_colnames).intersection(cols_after_pyranges_join)

    if merge is None:
        cols_tobe_added = right_gdf_cols
    else:
        cols_tobe_added = list(
            (x + y) for (x, y) 
            in itertools.product(right_gdf_cols, suffixes.values())
        )

    assert not set(left_gdf.annot_cols).intersection(cols_tobe_added)
    assert not set(reserved_colnames).intersection(cols_tobe_added)

    return right_gdf_cols, cols_tobe_added, winsorize


def winsorize_matched_right(
    matched_right, 
    groupkey, 
    winsorize,
):
    new_values_dict = dict()
    for colname in matched_right.columns:
        new_values = list()
        ser = matched_right[colname]
        for key, val in ser.groupby(groupkey):
            new_val = scipy.stats.mstats.winsorize(
                val, 
                limits=winsorize,
                nan_policy='omit',
            )
            new_values.append(np.ma.getdata(new_val))

        new_values_dict[colname] = np.concatenate(new_values)

    for colname, new_values in new_values_dict.items():
        matched_right[colname] = new_values


def merge_right(
    matched_right, 
    groupkey, 
    right_gdf_cols, 
    mergetype, 
    suffixes, 
    ddof,
    merge_args, 
    merge_kwargs,
):
    if mergetype in MERGETYPES_SIMPLEAGG:
        merged_right = merge_right_simpleagg(
            matched_right, groupkey, right_gdf_cols, mergetype, ddof, 
        )
    elif mergetype in MERGETYPES_WEIGHTED_SIMPLEAGG:
        merged_right = merge_right_weighted_simpleagg(
            matched_right, groupkey, right_gdf_cols, mergetype, ddof,
        )
    elif mergetype in MERGETYPES_BYLENGTH:
        merged_right = merge_right_bylength(
            matched_right, groupkey, right_gdf_cols, mergetype, ddof,
        )
    else:  # custom merge method
        merged_right = merge_right_custom(
            matched_right, 
            groupkey, 
            right_gdf_cols, 
            mergetype, 
            merge_args, 
            merge_kwargs,
        )

    merged_right.reset_index(drop=True, inplace=True)
    suffix = suffixes[get_suffixes_key(mergetype)]
    merged_right = merged_right.add_suffix(suffix, axis=1)

    return merged_right


def merge_right_custom(
    matched_right, 
    groupkey, 
    right_gdf_cols, 
    mergetype, 
    merge_args, 
    merge_kwargs,
):
    '''When *args or **kwargs are given to "agg", "func" receives
    DataFrame; otherwise "func" receives Series.
    '''
    aggresult_list = list()
    for x in right_gdf_cols:
        ser = matched_right.loc[:, x]
        aggresult = ser.groupby(groupkey, sort=False).agg(
            mergetype, 
            *merge_args, 
            engine=None, 
            engine_kwargs=None,
            **merge_kwargs,
        )
        aggresult.reset_index(drop=True, inplace=True)
        aggresult_list.append(aggresult)

    merged_right = pd.concat(aggresult_list, axis=1)
    return merged_right


def merge_right_simpleagg(
    matched_right, groupkey, right_gdf_cols, mergetype, ddof,
):
    assert mergetype in MERGETYPES_SIMPLEAGG

    groupby_obj = matched_right.groupby(groupkey, sort=False)[right_gdf_cols]

    if mergetype == 'first':
        merged_right = groupby_obj.first()
    elif mergetype == 'last':
        merged_right = groupby_obj.last()
    elif mergetype == 'mean':
        merged_right = groupby_obj.mean()
    elif mergetype == 'std':
        merged_right = groupby_obj.std(ddof=ddof)
    elif mergetype == 'max':
        merged_right = groupby_obj.max()
    elif mergetype == 'min':
        merged_right = groupby_obj.min()
    elif mergetype == 'median':
        merged_right = groupby_obj.median()
    elif mergetype == 'skew':
        merged_right = groupby_obj.skew()
    else:
        raise Exception(f'Unknown mergetype: {mergetype}')

    return merged_right


def merge_right_weighted_simpleagg(matched_right, groupkey, right_gdf_cols, mergetype, ddof):
    assert mergetype in MERGETYPES_WEIGHTED_SIMPLEAGG
    #assert (matched_right.index == pd.RangeIndex(start=0, end=matched_right.shape[0])).all()

    new_matched_right = make_weighted_right(matched_right, right_gdf_cols)
    groupby_obj = new_matched_right.groupby(groupkey, sort=False)[right_gdf_cols]

    if mergetype == 'weighted_mean':
        merged_right = groupby_obj.sum()
    else:
        raise Exception(f'Unknown mergetype: {mergetype}')

    return merged_right


def make_weighted_right(matched_right, right_gdf_cols):
    if 'Overlap' in matched_right.columns:
        lengths = matched_right['Overlap']
    else:
        lengths = matched_right['End_b'] - matched_right['Start_b']

    lengthsum = lengths.sum()

    src_data = {
        key: ((matched_right[key] * lengths) / lengthsum)
        for key in right_gdf_cols
    }
    weighted_matched_right = pd.DataFrame(src_data)

    return weighted_matched_right


def merge_right_bylength(matched_right, groupkey, right_gdf_cols, mergetype, ddof):
    assert mergetype in MERGETYPES_BYLENGTH
    #assert (matched_right.index == pd.RangeIndex(start=0, end=matched_right.shape[0])).all()

    if 'Overlap' in matched_right.columns:
        lengths = matched_right['Overlap']
    else:
        lengths = matched_right['End_b'] - matched_right['Start_b']

    groupby_obj = lengths.reset_index(drop=True).groupby(groupkey, sort=False)
    if mergetype == 'longest':
        indexes = groupby_obj.idxmax()
    elif mergetype == 'shortest':
        indexes = groupby_obj.idxmin()
    else:
        raise Exception(f'Unknown mergetype: {mergetype}')

    right_gdf_cols_idxs = [matched_right.columns.get_loc(x) for x in right_gdf_cols]
    merged_right = matched_right.iloc[indexes.to_numpy(), right_gdf_cols_idxs]

    return merged_right


########
# isec #
########

def isec_base(left_gdf, right_gdf):
    """Assume left_gdf and right_gdf has single identical chromosome"""

    left_gdf.sort()
    right_gdf.sort()

    left_start0s = left_gdf.start0s
    left_end0s = left_gdf.end0s
    right_start0s = right_gdf.start0s
    right_end0s = right_gdf.end0s


# isec boolean array makers

def check_intersect_num_leading_true(arr):
    result = tools.array_index(arr, False)
    if result is None:
        result = arr.shape[0]
    return result


def check_intersect_num_leading_false(arr):
    result = tools.array_index(arr, True)
    if result is None:
        result = arr.shape[0]
    return result


def check_intersect(
    left_start0s, left_end0s, right_start0s, right_end0s, 
    no_inner_overlap=False,
    sanitycheck=True,
):
    if (no_inner_overlap and sanitycheck):
        assert not genomedf_utils.check_interval_overlap(left_start0s, left_end0s)
        assert not genomedf_utils.check_interval_overlap(right_start0s, right_end0s)

    left_start0s = np.asarray(left_start0s)
    left_end0s = np.asarray(left_end0s)
    right_start0s = np.asarray(right_start0s)
    right_end0s = np.asarray(right_end0s)

    nonovlp_selector = np.logical_or(
        (left_end0s <= right_start0s[0]),
        (left_start0s >= right_end0s[-1]),
    )

    result = np.empty((len(left_start0s), len(right_start0s)), dtype=bool)
    result[nonovlp_selector, :] = False

    main_func = (
        check_intersect_relevant_noinneroverlap
        if no_inner_overlap else
        check_intersect_relevant
    )
    result[~nonovlp_selector, :] = main_func(
        left_start0s[~nonovlp_selector],
        left_end0s[~nonovlp_selector],
        right_start0s,
        right_end0s,
    )

    return result
    

def check_intersect_relevant(left_start0s, left_end0s, right_start0s, right_end0s):
    # less than left end
    rstart_lt_lend_list = list()
    for end0 in left_end0s:
        current_selector = (right_start0s < end0)
        rstart_lt_lend_list.append(current_selector)

    rstart_lt_lend = np.stack(rstart_lt_lend_list, axis=0)

    # greater than left start
    for lt_lend_selector, start0 in zip(rstart_lt_lend, left_start0s):
        lt_lend_selector[lt_lend_selector] = (
            right_end0s[lt_lend_selector] > start0
        )

    return rstart_lt_lend


def check_intersect_relevant_noinneroverlap(left_start0s, left_end0s, right_start0s, right_end0s):
    # less than left end
    rstart_lt_lend_list = list()
    numtrue_list = list()
    for end0 in left_end0s:
        if len(numtrue_list) == 0:
            last_numtrue = 0
        else:
            last_numtrue = numtrue_list[-1]

        current_selector = np.concatenate(
            [
                np.repeat(True, last_numtrue),
                (right_start0s[last_numtrue:] < end0),
            ]
        )
        rstart_lt_lend_list.append(current_selector)
        numtrue_list.append(check_intersect_num_leading_true(current_selector))

    rstart_lt_lend = np.stack(rstart_lt_lend_list, axis=0)

    #print(rstart_lt_lend)

    # greater than left start
    #last_num_leadingfalse_list = list()
    #last_leadingfalse = 0
    last_idx = len(left_start0s) - 1
    for idx, (start0, numtrue) in enumerate(zip(left_start0s, numtrue_list)):
        #if len(last_num_leadingfalse_list) == 0:
        #    last_leadingfalse = 0
        #else:
        #    last_leadingfalse = last_num_leadingfalse_list[-1]
        #    last_leadingfalse = check_intersect_num_leading_false(rstart_lt_lend[idx - 1, :])

        subarr = rstart_lt_lend[idx, :]
        subarr[subarr] = right_end0s[subarr] > start0
        if idx != last_idx:
            rstart_lt_lend[idx + 1, :numtrue] = subarr[:numtrue]

        #last_leadingfalse = check_intersect_num_leading_false(subsubarr)

    return rstart_lt_lend



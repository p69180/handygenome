import itertools
import warnings
import operator
import multiprocessing

import numpy as np
import pandas as pd
import pyranges as pr

import handygenome.common as common
import handygenome.workflow as workflow


LOGGER_INFO = workflow.get_debugging_logger(verbose=False)
LOGGER_DEBUG = workflow.get_debugging_logger(verbose=True)


def sort_df_by_coord(df, chromdict):
    coord_cols = df.loc[:, ['Chromosome', 'Start', 'End']]
    sorted_rows = sorted(
        coord_cols.iterrows(),
        key=(lambda x: (chromdict.contigs.index(x[1][0]), x[1][1], x[1][2])),
    )
    sorted_rows_indexes = [x[0] for x in sorted_rows]
    return df.iloc[sorted_rows_indexes, :].reset_index(drop=True)


def isec_union(gr1, gr2):
    isec = gr1.intersect(gr2)
    gr1_diff_isec = gr1.subtract(isec)
    gr2_diff_isec = gr2.subtract(isec)
    result = pr.concat([gr1_diff_isec, isec, gr2_diff_isec])
    result = result[[]]  # remove annotation columns
    result.sort()
    return result


# join #


#def inner_join(left_gr, right_gr):
#    """For handling grs with lots of columns"""
#    new_right_df = right_gr.df.loc[:, ['Chromosome', 'Start', 'End']]
#    new_right_df.insert(3, 'index', range(new_right_df.shape[0]))
#    new_right_gr = pr.PyRanges(new_right_df)
#    joined_gr = left_gr.join(new_right_gr)
#    right_values_df = right_gr.df.iloc[joined_gr.index.to_list(), 3:].reset_index(drop=True)
#    right_values_df.reset_index(drop=True, inplace=True)
#    
#    old_columns = left_gr.columns.to_list()
#    new_columns = [
#        (f'{x}_b' if x in old_columns else x)
#        for x in right_values_df.columns
#    ]
#    right_values_df.columns = new_columns
#
#    result_df = pd.concat([joined_gr.df.iloc[:, :-1], right_values_df], axis=1)
#    return pr.PyRanges(result_df)
#
#
#def add_new_coordinates(joined_gr, how):
#    """Accepts a result of 'join' method, which has columns 'Start_b' and 'End_b'"""
#    assert how in ("union", "intersection", "swap")
#    assert {'Start_b', 'End_b'}.issubset(joined_gr.columns)
#
#    Start = joined_gr.Start
#    End = joined_gr.End
#
#    _ = joined_gr.new_position(how)
#    new_Start = joined_gr.Start
#    new_End = joined_gr.End
#
#    joined_gr.Start = Start
#    joined_gr.End = End
#    joined_gr.new_Start = new_Start
#    joined_gr.new_End = new_End
#
#    return joined_gr
#
#
## this is of similar speed to add_new_coordinates
#def add_new_coordinates_new(joined_gr, how):
#    """Accepts a result of 'join' method, which has columns 'Start_b' and 'End_b'"""
#    assert how in ("union", "intersection")
#    assert {'Start_b', 'End_b'}.issubset(joined_gr.columns)
#
#    df = joined_gr.df
#    starts_subarr = df.loc[:, ['Start', 'Start_b']].to_numpy()
#    ends_subarr = df.loc[:, ['End', 'End_b']].to_numpy()
#
#    if how == 'union':
#        new_starts = starts_subarr.min(axis=1)
#        new_ends = ends_subarr.max(axis=1)
#    elif how == 'intersection':
#        new_starts = starts_subarr.max(axis=1)
#        new_ends = ends_subarr.min(axis=1)
#
#    joined_gr.new_Start = new_starts
#    joined_gr.new_End = new_ends
#
#    return joined_gr
#
#
## join helper
#def extract_unmatched_rows(left_gr, right_gr, joined_gr, find_nearest, added_columns):
#    unmatched_rows_gr = left_gr.subtract(joined_gr)
#    if unmatched_rows_gr.empty:
#        return None
#
#    if find_nearest:
#        unmatched_rows_gr = unmatched_rows_gr.nearest(
#            right_gr, overlap=False, how=None,
#        )
#    else:
#        nrow = unmatched_rows_gr.df.shape[0]
#        for colname in added_columns:
#            setattr(unmatched_rows_gr, colname, np.repeat(np.nan, nrow))
#
#    # remove unused columns
#    cols_to_drop = set(unmatched_rows_gr.columns).intersection(
#        {'Start_b', 'End_b', 'Distance'}
#    )
#    unmatched_rows_gr = unmatched_rows_gr[
#        unmatched_rows_gr.columns.drop(cols_to_drop).to_list()
#    ]
#    return unmatched_rows_gr
#
#
#def join_sanity_check(left_gr, right_gr):
#    assert isinstance(left_gr, pr.PyRanges)
#    assert isinstance(right_gr, pr.PyRanges)
#
#    common_cols = set.intersection(set(left_gr.columns), set(right_gr.columns))
#    assert not common_cols.difference({'Chromosome', 'Start', 'End'})
#
#    assert len(set(right_gr.columns).difference(['Chromosome', 'Start', 'End'])) > 0, (
#        f'"right_gr" Must have columns other than "Chromosome", "Start", "End".'
#    )
#
#
#def join_new(left_gr, right_gr, how='inner', find_nearest=False, merge='mean', as_gr=True):
#    assert merge in {'mean'}
#    assert how in {'inner', 'left'}
#    join_sanity_check(left_gr, right_gr)
#
#    # join
#    with warnings.catch_warnings(): 
#        warnings.simplefilter('ignore', category=FutureWarning)
#        joined_gr = left_gr.join(right_gr)
#
#    joined_df = joined_gr.df
#    added_cols = [x for x in right_gr.columns.to_list() if x not in ('Chromosome', 'Start', 'End')]
#    #if joined_gr.empty:
#        #added_cols = pd.Series(joined_gr.columns).drop(left_gr.columns).drop(['Start_b', 'End_b']).to_list()
#    #else:
#        #added_cols = joined_gr.columns.drop(left_gr.columns).drop(['Start_b', 'End_b']).to_list()
#
#    # unmatched rows
#    if how == 'left':
#        unmatched_rows_gr = extract_unmatched_rows(
#            left_gr, right_gr, joined_gr, find_nearest, added_cols,
#        )
#    elif how == 'inner':
#        unmatched_rows_gr = None
#
#    # merge
#    if joined_gr.empty:
#        matched_rows_gr = joined_gr
#    else:
#        coord_cols = np.concatenate(
#            [
#                joined_df['Chromosome'].cat.codes.values[:, np.newaxis],
#                joined_df.loc[:, ['Start', 'End']],
#            ],
#            axis=1,
#        )
#        values, counts = common.array_grouper(coord_cols, omit_values=True)
#        indexer = list(
#            itertools.chain.from_iterable(
#                itertools.repeat(val, count) for val, count in enumerate(counts)
#            )
#        )
#        grouper = joined_df.groupby(indexer)
#        with warnings.catch_warnings(): 
#            warnings.simplefilter('ignore', category=FutureWarning)
#            aggresult = grouper[added_cols].mean().reset_index(drop=True)
#
#        matched_rows = pd.concat(
#            [
#                joined_df.loc[:, left_gr.columns.to_list()].drop_duplicates(
#                    ['Chromosome', 'Start', 'End'], keep='first',
#                ).reset_index(drop=True), 
#                aggresult,
#            ], 
#            axis=1,
#        )
#        matched_rows_gr = pr.PyRanges(matched_rows)
#
#    # merge matched and nonmatched grs
#    if unmatched_rows_gr is None:
#        result = matched_rows_gr
#    else:
#        if matched_rows_gr.empty:
#            result = unmatched_rows_gr
#        else:
#            result = pr.concat([unmatched_rows_gr, matched_rows_gr]).sort()
#
#    # return
#    if not as_gr:
#        result = result.df
#
#    return result
#
#
#def join(left_gr, right_gr, how='inner', merge=None, find_nearest=False, as_gr=True, verbose=False):
#    """Args:
#        find_nearest: Match a nearest row from "right_gr" for each
#            non-overlapping row of "left_gr". Only relevant when "how" is "left".
#    Returns:
#        A pyranges.PyRanges object 
#    """
#    def merge_helper(joined_gr, added_columns, merge):
#        # modify joined_gr
#        joined_gr = joined_gr.sort()
#        joined_gr = add_new_coordinates(joined_gr, 'intersection')
#            # now 'new_Start' and 'new_End' columns are added
#        joined_gr.isec_length = joined_gr.new_End - joined_gr.new_Start
#
#        # prepare parameters
#        added_value_columns = list(set(added_columns).difference({'Start_b', 'End_b'}))
#        new_columns = ['Chromosome', 'Start', 'End', 'isec_length'] + added_value_columns
#        joined_gr_subset = joined_gr[new_columns]
#
#        # make merged annotation data
#        data = _merge_helper_make_data_vectorize(added_value_columns, joined_gr_subset, merge)
#        #data = _merge_helper_make_data_iterrows(added_value_columns, joined_gr_subset, merge)
#
#        # prepare final result
#        result_df = joined_gr.df.drop_duplicates(
#            subset=['Chromosome', 'Start', 'End'], keep='first', ignore_index=True,
#        )
#        with warnings.catch_warnings(): 
#            warnings.simplefilter('ignore', category=FutureWarning)
#            result_df.loc[:, added_value_columns] = data
#        result = pr.PyRanges(result_df, int64=False)
#        return result
#
#    def _merge_helper_make_data_iterrows(added_value_columns, joined_gr_subset, merge):
#        def _merge_values_weighted_mean(value_rowlist, lengths):
#            return np.average(np.array(value_rowlist), axis=0, weights=lengths)
#
#        def _merge_values_mean(value_rowlist, lengths):
#            return np.average(np.array(value_rowlist), axis=0, weights=None)
#
#        def _merge_values_longest_nonmerge(value_rowlist, lengths):
#            return value_rowlist[np.argmax(lengths)]
#
#        def _merge_values_longest(value_rowlist, lengths):
#            len_val_dict = dict()
#            for vals, l in zip(value_rowlist, lengths):
#                vals = tuple(vals)
#                len_val_dict.setdefault(vals, 0)
#                len_val_dict[vals] += l
#            return max(len_val_dict.items(), key=operator.itemgetter(1))[0]
#
#        if merge == 'mean':
#            values_merger = _merge_values_mean
#        elif merge == 'weighted_mean':
#            values_merger = _merge_values_weighted_mean
#        elif merge == 'longest_nonmerge':
#            values_merger = _merge_values_longest_nonmerge
#        elif merge == 'longest':
#            values_merger = _merge_values_longest
#
#        data = list()
#        for key, subiter in itertools.groupby(
#            (x[1] for x in joined_gr_subset.df.iterrows()),
#            key=(lambda row: tuple(row.iloc[:3])),
#        ):
#            value_rowlist = list()
#            lengths = list()
#            for row in subiter:
#                value_rowlist.append(row.iloc[4:])
#                lengths.append(row.iloc[3])
#
#            if len(value_rowlist) == 1:
#                data.append(value_rowlist[0])
#            else:
#                data.append(values_merger(value_rowlist, lengths))
#
#        return np.array(data)
#
#    # too slow
#    def _merge_helper_make_data_pdgroupby(added_value_columns, joined_gr_subset, merge):
#        data = list()
#        for key, subdf in joined_gr_subset.df.groupby(['Chromosome', 'Start', 'End']):
#            if subdf.shape[0] == 1:
#                data.append(subdf.iloc[0, 4:])
#            else:
#                if merge == 'mean':
#                    data.append(subdf.iloc[:, 4:].mean())
#                elif merge == 'weighted_mean':
#                    data.append(np.average(subdf.iloc[:, 4:], axis=0, weights=subdf.iloc[:, 3]))
#                elif merge == 'longest_nonmerge':
#                    data.append(subdf.iloc[np.argmax(subdf.iloc[:, 3]), 4:])
#                elif merge == 'longest':
#                    len_val_dict = {
#                        tuple(subkey): subsubdf.iloc[:, 3].sum()
#                        for subkey, subsubdf in subdf.groupby(added_value_columns)
#                    }
#                    data.append(
#                        max(len_val_dict.items(), key=operator.itemgetter(1))[0]
#                    )
#
#        return np.array(data)
#
#    # this is the fastest!
#    def _merge_helper_make_data_vectorize(added_value_columns, joined_gr_subset, merge):
#        def _merge_annotvals_weighted_mean(annotvals, lengths):
#            return np.average(np.array(annotvals), axis=0, weights=lengths)
#
#        def _merge_annotvals_mean(annotvals, lengths):
#            #return np.average(np.array(annotvals), axis=0, weights=None)
#            return np.nanmean(np.array(annotvals), axis=0)
#
#        def _merge_annotvals_longest_nonmerge(annotvals, lengths):
#            return annotvals[np.argmax(lengths)]
#
#        def _merge_annotvals_longest(annotvals, lengths):
#            len_val_dict = dict()
#            for vals, l in zip(annotvals, lengths):
#                len_val_dict.setdefault(vals, 0)
#                len_val_dict[vals] += l
#            return max(len_val_dict.items(), key=operator.itemgetter(1))[0]
#
#        coord_list = joined_gr_subset.df.apply(lambda row: tuple(row.iloc[:3]), axis=1)
#        length_list = joined_gr_subset.df.iloc[:, 3]
#        annotval_list = joined_gr_subset.df.apply((lambda row: tuple(row.iloc[4:])), axis=1)
#
#        if merge == 'mean':
#            annotvals_merger = _merge_annotvals_mean
#        elif merge == 'weighted_mean':
#            annotvals_merger = _merge_annotvals_weighted_mean
#        elif merge == 'longest_nonmerge':
#            annotvals_merger = _merge_annotvals_longest_nonmerge
#        elif merge == 'longest':
#            annotvals_merger = _merge_annotvals_longest
#
#        data = list()
#        for key, subiter in itertools.groupby(
#            zip(coord_list, length_list, annotval_list), key=operator.itemgetter(0),
#        ):
#            coords = list()
#            lengths = list()
#            annotvals = list()
#            for x in subiter:
#                coords.append(x[0])
#                lengths.append(x[1])
#                annotvals.append(x[2])
#
#            if len(coords) == 1:
#                data.extend(annotvals)
#            else:
#                data.append(annotvals_merger(annotvals, lengths))
#
#        return np.array(data)
#
#    # main
#    assert merge in ('mean', 'weighted_mean', 'first', 'longest', 'longest_nonmerge', None)
#    assert how in ('left', 'inner')
#    join_sanity_check(left_gr, right_gr)
#
#    # do join
#    #common.funclogger(1)
#    with warnings.catch_warnings(): 
#        warnings.simplefilter('ignore', category=FutureWarning)
#        joined_gr = left_gr.join(right_gr)  
#            # this is inner join
#            # setting "how='left'" does not work and results in inner join (230209)
#    #common.funclogger(2)
#    final_added_columns = [x for x in right_gr.columns.to_list() if x not in ('Chromosome', 'Start', 'End')]
#    added_columns = ['Start_b', 'End_b'] + final_added_columns
#    #if joined_gr.empty:
#        #added_columns = pd.Series(joined_gr.columns).drop(left_gr.columns).to_list()
#    #else:
#        #added_columns = joined_gr.columns.drop(left_gr.columns).to_list()
#        # This includes "Start_b" and "End_b"
#    #common.funclogger(3)
#
#    # handle unmatched rows
#    if how == 'left':
#        unmatched_rows_gr = extract_unmatched_rows(
#            left_gr, right_gr, joined_gr, find_nearest, final_added_columns,
#        )
#    elif how == 'inner':
#        unmatched_rows_gr = None
#    #common.funclogger(4)
#
#    if joined_gr.empty:
#        matched_rows_gr = joined_gr
#    else:
#        # handle matched rows - merge rows with identical (chrom, start, end)
#        if merge is None:
#            matched_rows_gr = joined_gr
#        elif merge == 'first':
#            matched_rows_gr = pr.PyRanges(
#                joined_gr.df.drop_duplicates(
#                    subset=['Chromosome', 'Start', 'End'], keep='first', ignore_index=True,
#                )
#            )
#        elif merge in ('mean', 'weighted_mean', 'longest', 'longest_nonmerge'):
#            matched_rows_gr = merge_helper(joined_gr, added_columns, merge)
#        #common.funclogger(5)
#
#        # handle matched rows - remove unused columns
#        cols_to_drop = set(matched_rows_gr.columns).intersection(
#            {'Start_b', 'End_b', 'new_Start', 'new_End', 'isec_length'}
#        )
#        matched_rows_gr = matched_rows_gr[
#            matched_rows_gr.columns.drop(cols_to_drop).to_list()
#        ]
#
#    # concat matched and unmatched rows
#    if unmatched_rows_gr is None:
#        result = matched_rows_gr
#    else:
#        if matched_rows_gr.empty:
#            result = unmatched_rows_gr
#        else:
#            result = pr.concat([unmatched_rows_gr, matched_rows_gr]).sort()
#
#    # return
#    if as_gr:
#        return result
#    else:
#        return result.df


###############################

def join_preprocess_df(df):
    df = common.arg_into_df(df)

    leading_cols = ['Chromosome', 'Start', 'End']
    if df.columns.to_list()[:3] != leading_cols:
        df = df.loc[
            :, 
            (leading_cols + df.columns.drop(leading_cols).to_list())
        ]
    return df.reset_index(drop=True, inplace=False)


def join_df_sanitycheck(left_df, right_df, index_col):
    common_cols = ["Chromosome", "Start", "End"]

    for df in (left_df, right_df):
        if not set(common_cols).issubset(df.columns):
            raise Exception(f'Input dataframe must include columns {common_cols}')
        if index_col in df.columns:
            raise Exception(f'Input dataframe must not have a column named "{index_col}"')

    left_annot_cols = left_df.columns.drop(common_cols)
    right_annot_cols = right_df.columns.drop(common_cols)

    assert not set(left_annot_cols).intersection(set(right_annot_cols)), f'Input DataFrames has overlapping annotation columns'


def check_intv_invtlist_overlap(intv_start0, intv_end0, intvlist_start0s, intvlist_end0s):
    return np.logical_and(
        intv_end0 > intvlist_start0s,
        intv_start0 < intvlist_end0s,
    )


def compare_coords(left_df, right_df):
    pos_compare_1 = (
        right_df.End.to_numpy() 
        > left_df.Start.to_numpy()[:, np.newaxis]
    )
    pos_compare_2 = (
        right_df.Start.to_numpy() 
        < left_df.End.to_numpy()[:, np.newaxis]
    )
    compare_result = np.logical_and(pos_compare_1, pos_compare_2)

    return compare_result


def compare_coords_new(left_df, right_df):
    """Assumes input dataframes are sorted by [Start, End]"""
    left_start = left_df.Start.to_numpy()
    left_end = left_df.End.to_numpy()
    right_start = right_df.Start.to_numpy()
    right_end = right_df.End.to_numpy()
    left_nrow = left_df.shape[0]
    right_nrow = right_df.shape[0]

    compare_result = np.tile(False, (left_nrow, right_nrow))
    if left_end[-1] <= right_start[0]:
        pass
    elif right_end[-1] <= left_start[0]:
        pass
    else:
        # get matched left_df indexes
        left_idx_start0 = np.where(left_end > right_start[0])[0][0]
        where_result = np.where(left_start >= right_end[-1])[0]
        if len(where_result) == 0:
            left_idx_end0 = left_nrow
        else:
            left_idx_end0 = where_result[0]

        right_idx = 0
        for left_idx in range(left_idx_start0, left_idx_end0):
            compare1 = (left_start[left_idx] < right_end[right_idx:])
            compare2 = (left_end[left_idx] > right_start[right_idx:])
            compare = np.logical_and(compare1, compare2)

            compare_result[left_idx, right_idx:] = compare
            where_result = np.where(compare)[0]
            if len(where_result) > 0:
                right_idx += where_result[0]

    return compare_result


def fetch_nearest(joined_df, right_df, left_rowidxs, annot_columns):
    unmatched_row_selector = ~joined_df.index.isin(left_rowidxs)
    unmatched_row_coords = joined_df.loc[unmatched_row_selector, :]

    distances = calc_distances(
        unmatched_row_coords['Start'].to_numpy(), 
        unmatched_row_coords['End'].to_numpy(), 
        right_df['Start'].to_numpy(), 
        right_df['End'].to_numpy(),
    )
    nearest_right_df_indexes = np.argmin(distances, axis=1)
    joined_df.loc[unmatched_row_selector, annot_columns] = (
        right_df.iloc[nearest_right_df_indexes, :].loc[:, annot_columns].to_numpy()
    )
    return joined_df


def join_singlechrom_dfs_new(
    left_df, right_df, how, keep_right_coords, find_nearest, annot_columns, merge, logger,
):
    logger.debug(f'Sorting dataframes')
    left_df.sort_values(['Start', 'End'], inplace=True)
    right_df.sort_values(['Start', 'End'], inplace=True)

    logger.debug(f'Comparing coordinates')
    compare_result = compare_coords(left_df, right_df)
    left_rowidxs, right_rowidxs = np.where(compare_result)

    logger.debug(f'Subsetting right df')
    if keep_right_coords:
        left_subdf = left_df.iloc[left_rowidxs, :]
        right_subdf = right_df.iloc[right_rowidxs, 3:]
        right_subdf['overlap_length'] = calc_overlap_lengths(
            left_subdf['Start'].to_numpy(), 
            left_subdf['End'].to_numpy(), 
            right_subdf['Start'].to_numpy(), 
            right_subdf['End'].to_numpy(),
        )
    else:
        right_subdf = right_df.iloc[right_rowidxs, 3:]
    right_subdf.index = left_rowidxs
    
    # do merge
    logger.debug(f'Doing merge')
    merged_right_subdf = merge_right_subdf(right_subdf, merge, annot_columns)

    # join with left_df
    logger.debug(f'Joining with left df')
    joined_df = left_df.join(merged_right_subdf, how=how)

    # fetch nearest
    if (how == 'left') and find_nearest:
        logger.debug(f'Doing fetch nearest')
        joined_df = fetch_nearest(joined_df, right_df, left_rowidxs, annot_columns)

    return joined_df


def join_singlechrom_dfs(left_df, right_df, how, keep_right_coords, find_nearest, annot_columns, merge, logger):
    current_chrom = left_df.Chromosome[0]
    logger.debug(f'Beginning chromosome {current_chrom}')

    #logger.debug(f'Comparing coordinates')
    compare_result = compare_coords(left_df, right_df)
    left_rowidxs, right_rowidxs = np.where(compare_result)

    #logger.debug(f'Subsetting right df')
    if keep_right_coords:
        right_subdf = right_df.iloc[right_rowidxs, 3:]
        right_subdf['overlap_length'] = calc_overlap_lengths(
            left_df['Start'].iloc[left_rowidxs].to_numpy(), 
            left_df['End'].iloc[left_rowidxs].to_numpy(), 
            right_df['Start'].iloc[right_rowidxs].to_numpy(), 
            right_df['End'].iloc[right_rowidxs].to_numpy(),
        )
    else:
        right_subdf = right_df.iloc[right_rowidxs, 3:]
    right_subdf.index = left_rowidxs
    
    # do merge
    #logger.debug(f'Doing merge')
    right_subdf_merged = merge_right_subdf(right_subdf, merge, annot_columns)

    # join with left_df
    #logger.debug(f'Joining with left df')
    joined_df = left_df.join(right_subdf_merged, how=how)

    # fetch nearest
    if (how == 'left') and find_nearest:
        #logger.debug(f'Doing fetch nearest')
        joined_df = fetch_nearest(joined_df, right_df, left_rowidxs, annot_columns)

    logger.debug(f'Finished chromosome {current_chrom}')

    return joined_df


def dedup_by_index(df, index_col, keep):
    indexes = df.index.to_numpy()
    diff_indexes = np.nonzero(np.diff(indexes))[0]
    if keep == 'last':
        selector = np.concatenate([diff_indexes, [df.shape[0] - 1]])
    elif keep == 'first':
        selector = np.concatenate([[-1], diff_indexes]) + 1

    result = df.iloc[selector, :]
    result.reset_index(inplace=True, drop=True)
    return result


def merge_right_subdf(right_subdf, merge, annot_columns):
    if merge is None:
        return right_subdf
    elif merge == 'first':
        return dedup_by_index(right_subdf, index_col=None, keep='first')
    elif merge == 'last':
        return dedup_by_index(right_subdf, index_col=None, keep='last')
    else:
        indexes = right_subdf.index.to_numpy()

        if merge in ('mean',):
            right_subdf_merged = right_subdf.groupby(by=indexes, axis=0, sort=False).mean()
        elif merge in ('weighted_mean', 'longest'):
            if merge == 'weighted_mean':
                def subdf_handler(subdf):
                    weights = subdf['overlap_length'].to_numpy()
                    if np.isnan(weights.sum()):
                        #return subdf.drop(columns='overlap_length').iloc[0, :]
                        return subdf.iloc[0, :]
                    else:
                        #return subdf.drop(columns='overlap_length').apply(
                        #    lambda x: common.nanaverage(x.to_numpy(), weights=weights)
                        #)
                        return subdf.apply(
                            lambda x: common.nanaverage(x.to_numpy(), weights=weights)
                        )
            elif merge == 'longest':
                def subdf_handler(subdf):
                    idx = np.argmax(subdf['overlap_length'])
                    return subdf.iloc[idx, :]

            data_gen = (
                subdf_handler(subdf)
                for key, subdf in right_subdf.groupby(by=indexes, axis=0, sort=False)
            )
            right_subdf_merged = pd.DataFrame.from_records(data_gen)

        return right_subdf_merged


def calc_distances(starts_left, ends_left, starts_right, ends_right):
    result = np.empty((len(starts_left), len(starts_right)), dtype=int)

    left_on_left = ends_left[:, np.newaxis] <= starts_right
    right_on_left = starts_left[:, np.newaxis] >= ends_right
    overlaps = np.logical_and(
        np.logical_not(left_on_left), np.logical_not(right_on_left)
    )

    left_indexes, right_indexes = np.where(left_on_left)
    result[left_on_left] = starts_right[right_indexes] - ends_left[left_indexes]

    left_indexes, right_indexes = np.where(right_on_left)
    result[right_on_left] = starts_left[left_indexes] - ends_right[right_indexes]

    result[overlaps] = 0

    return result


def calc_overlap_lengths(starts_left, ends_left, starts_right, ends_right):
    assert len(starts_left) == len(starts_right)
    result = (
        np.min(np.array([ends_left, ends_right]), axis=0)
        - np.max(np.array([starts_left, starts_right]), axis=0)
    )
    result[result < 0] = 0
    return result


def group_df_bychrom(df):
    result = dict()
    for key, subdf in df.groupby(df['Chromosome'].to_numpy()):
        result[key] = subdf.reset_index(drop=True)
    return result


def concat_dfs_bychrom(how, sort, common_chroms, leftonly_chroms, chromdict, joined_bychrom, left_bychrom):
    if how == 'inner':
        if sort:
            sorted_chroms = sorted(
                common_chroms, 
                key=(lambda x: chromdict.contigs.index(x)),
            )
            concat_dfs = (
                joined_bychrom[x] for x in sorted_chroms
            )
        else:
            concat_dfs = iter(joined_bychrom.values())
    elif how == 'left':
        if sort:
            sorted_chroms = sorted(
                itertools.chain(common_chroms, leftonly_chroms),
                key=(lambda x: chromdict.contigs.index(x)),
            )
            concat_dfs = (
                (
                    left_bychrom[x]  
                    if x in leftonly_chroms else
                    joined_bychrom[x]
                )
                for x in sorted_chroms
            )
        else:
            concat_dfs = itertools.chain(
                iter(joined_bychrom.values()),
                (left_bychrom[chrom] for chrom in leftonly_chroms),
            )

    result = pd.concat(concat_dfs, axis=0)
    return result


#def join_main_newest_nonmulti(
#    left_df, right_df, how, merge, find_nearest, index_col, sort, chromdict, logger,
#    nproc,
#):
#    ###
#    annot_columns = right_df.columns[3:].to_list()
#
#    ###
#    logger.debug(f'Grouping dataframes by chromosome')
#    left_bychrom = group_df_bychrom(left_df)
#    right_bychrom = group_df_bychrom(right_df)
#
#    left_chroms = set(left_bychrom.keys())
#    right_chroms = set(right_bychrom.keys())
#    common_chroms = left_chroms.intersection(right_chroms)
#    leftonly_chroms = left_chroms.difference(right_chroms)
#    #rightonly_chroms = right_chroms.difference(left_chroms)
#
#    ###
#    logger.debug(f'Doing join by chromosome')
#    keep_right_coords = (merge in ('weighted_mean', 'longest'))
#    joined_bychrom = dict()
#    for chrom in common_chroms:
#        logger.debug(f'chrom {chrom}')
#        joined_df = join_singlechrom_dfs(
#            left_bychrom[chrom], 
#            right_bychrom[chrom], 
#            how, 
#            keep_right_coords,
#            find_nearest,
#            annot_columns,
#            merge,
#            logger,
#        )
#        joined_bychrom[chrom] = joined_df
#
#    ###
#    logger.debug(f'Concatenating dataframes of all chromosomes')
#    result = concat_dfs_bychrom(how, sort, common_chroms, leftonly_chroms, chromdict, joined_bychrom, left_bychrom)
#    result.reset_index(inplace=True, drop=True)
#    if keep_right_coords:
#        result.drop(columns=['Start_right', 'End_right'], inplace=True)
#
#    return result


def join_main_newest(
    left_df, right_df, how, merge, find_nearest, index_col, sort, chromdict, logger,
    nproc,
):
    ###
    annot_columns = right_df.columns[3:].to_list()

    ###
    logger.debug(f'Grouping dataframes by chromosome')
    left_bychrom = group_df_bychrom(left_df)
    right_bychrom = group_df_bychrom(right_df)

    left_chroms = set(left_bychrom.keys())
    right_chroms = set(right_bychrom.keys())
    common_chroms = left_chroms.intersection(right_chroms)
    leftonly_chroms = left_chroms.difference(right_chroms)
    #rightonly_chroms = right_chroms.difference(left_chroms)

    ###
    logger.debug(f'Doing join by chromosome')
    keep_right_coords = (merge in ('weighted_mean', 'longest'))
    with multiprocessing.Pool(nproc) as pool:
        args = (
            (
                left_bychrom[chrom], 
                right_bychrom[chrom], 
                how, 
                keep_right_coords,
                find_nearest,
                annot_columns,
                merge,
                logger,
            )
            for chrom in common_chroms
        )
        joined_df_list = pool.starmap(join_singlechrom_dfs, args)
    joined_bychrom = dict(zip(common_chroms, joined_df_list))

    ###
    logger.debug(f'Concatenating dataframes of all chromosomes')
    result = concat_dfs_bychrom(
        how, sort, common_chroms, leftonly_chroms, chromdict, joined_bychrom, left_bychrom,
    )
    result.reset_index(inplace=True, drop=True)

    return result


def join_newest(
    left_df, right_df, how='inner', find_nearest=False, merge=None, as_gr=False,
    nproc=1,
    sort=False, chromdict=None, refver=None,
    verbose=False,
):
    index_col = '__index'

    # sanity check
    assert merge in (
        'mean', 'weighted_mean', 'longest', #'longest_nonmerge', 
        'first', 'last',
        None,
    )
    assert how in ('left', 'inner')
    if sort and ((chromdict is None) and (refver is None)):
        raise Exception(f'If "sort" is True, "chromdict" or "refver" must be given')
    join_df_sanitycheck(left_df, right_df, index_col)

    # set logger
    logger = (LOGGER_DEBUG if verbose else LOGGER_INFO)

    # arg handling
    left_df = join_preprocess_df(left_df)
    right_df = join_preprocess_df(right_df)
    if (chromdict is None) and (refver is not None):
        chromdict = common.DEFAULT_CHROMDICTS[refver]

    # join
    result = join_main_newest(left_df, right_df, how, merge, find_nearest, index_col, sort, chromdict, logger, nproc)

    # result
    if as_gr:
        result = pr.PyRanges(result)

    return result


join = join_newest



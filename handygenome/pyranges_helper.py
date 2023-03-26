import itertools
import warnings
import operator

import numpy as np
import pandas as pd
import pyranges as pr

import handygenome.common as common


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


def inner_join(left_gr, right_gr):
    """For handling grs with lots of columns"""
    new_right_df = right_gr.df.loc[:, ['Chromosome', 'Start', 'End']]
    new_right_df.insert(3, 'index', range(new_right_df.shape[0]))
    new_right_gr = pr.PyRanges(new_right_df)
    joined_gr = left_gr.join(new_right_gr)
    right_values_df = right_gr.df.iloc[joined_gr.index.to_list(), 3:].reset_index(drop=True)
    right_values_df.reset_index(drop=True, inplace=True)
    
    old_columns = left_gr.columns.to_list()
    new_columns = [
        (f'{x}_b' if x in old_columns else x)
        for x in right_values_df.columns
    ]
    right_values_df.columns = new_columns

    result_df = pd.concat([joined_gr.df.iloc[:, :-1], right_values_df], axis=1)
    return pr.PyRanges(result_df)


def add_new_coordinates(joined_gr, how):
    """Accepts a result of 'join' method, which has columns 'Start_b' and 'End_b'"""
    assert how in ("union", "intersection", "swap")
    assert {'Start_b', 'End_b'}.issubset(joined_gr.columns)

    Start = joined_gr.Start
    End = joined_gr.End

    _ = joined_gr.new_position(how)
    new_Start = joined_gr.Start
    new_End = joined_gr.End

    joined_gr.Start = Start
    joined_gr.End = End
    joined_gr.new_Start = new_Start
    joined_gr.new_End = new_End

    return joined_gr


# this is of similar speed to add_new_coordinates
def add_new_coordinates_new(joined_gr, how):
    """Accepts a result of 'join' method, which has columns 'Start_b' and 'End_b'"""
    assert how in ("union", "intersection")
    assert {'Start_b', 'End_b'}.issubset(joined_gr.columns)

    df = joined_gr.df
    starts_subarr = df.loc[:, ['Start', 'Start_b']].to_numpy()
    ends_subarr = df.loc[:, ['End', 'End_b']].to_numpy()

    if how == 'union':
        new_starts = starts_subarr.min(axis=1)
        new_ends = ends_subarr.max(axis=1)
    elif how == 'intersection':
        new_starts = starts_subarr.max(axis=1)
        new_ends = ends_subarr.min(axis=1)

    joined_gr.new_Start = new_starts
    joined_gr.new_End = new_ends

    return joined_gr


# join helper
def extract_unmatched_rows(left_gr, right_gr, joined_gr, find_nearest, added_columns):
    unmatched_rows_gr = left_gr.subtract(joined_gr)
    if unmatched_rows_gr.empty:
        return None

    if find_nearest:
        unmatched_rows_gr = unmatched_rows_gr.nearest(
            right_gr, overlap=False, how=None,
        )
    else:
        nrow = unmatched_rows_gr.df.shape[0]
        for colname in added_columns:
            setattr(unmatched_rows_gr, colname, np.repeat(np.nan, nrow))

    # remove unused columns
    cols_to_drop = set(unmatched_rows_gr.columns).intersection(
        {'Start_b', 'End_b', 'Distance'}
    )
    unmatched_rows_gr = unmatched_rows_gr[
        unmatched_rows_gr.columns.drop(cols_to_drop).to_list()
    ]
    return unmatched_rows_gr


def join_sanity_check(left_gr, right_gr):
    assert isinstance(left_gr, pr.PyRanges)
    assert isinstance(right_gr, pr.PyRanges)

    common_cols = set.intersection(set(left_gr.columns), set(right_gr.columns))
    assert not common_cols.difference({'Chromosome', 'Start', 'End'})

    assert len(set(right_gr.columns).difference(['Chromosome', 'Start', 'End'])) > 0, (
        f'"right_gr" Must have columns other than "Chromosome", "Start", "End".'
    )


def join_new(left_gr, right_gr, how='inner', find_nearest=False, merge='mean', as_gr=True):
    assert merge in {'mean'}
    assert how in {'inner', 'left'}
    join_sanity_check(left_gr, right_gr)

    # join
    with warnings.catch_warnings(): 
        warnings.simplefilter('ignore', category=FutureWarning)
        joined_gr = left_gr.join(right_gr)

    joined_df = joined_gr.df
    added_cols = [x for x in right_gr.columns.to_list() if x not in ('Chromosome', 'Start', 'End')]
    #if joined_gr.empty:
        #added_cols = pd.Series(joined_gr.columns).drop(left_gr.columns).drop(['Start_b', 'End_b']).to_list()
    #else:
        #added_cols = joined_gr.columns.drop(left_gr.columns).drop(['Start_b', 'End_b']).to_list()

    # unmatched rows
    if how == 'left':
        unmatched_rows_gr = extract_unmatched_rows(
            left_gr, right_gr, joined_gr, find_nearest, added_cols,
        )
    elif how == 'inner':
        unmatched_rows_gr = None

    # merge
    if joined_gr.empty:
        matched_rows_gr = joined_gr
    else:
        coord_cols = np.concatenate(
            [
                joined_df['Chromosome'].cat.codes.values[:, np.newaxis],
                joined_df.loc[:, ['Start', 'End']],
            ],
            axis=1,
        )
        values, counts = common.array_grouper(coord_cols, omit_values=True)
        indexer = list(
            itertools.chain.from_iterable(
                itertools.repeat(val, count) for val, count in enumerate(counts)
            )
        )
        grouper = joined_df.groupby(indexer)
        with warnings.catch_warnings(): 
            warnings.simplefilter('ignore', category=FutureWarning)
            aggresult = grouper[added_cols].mean().reset_index(drop=True)

        matched_rows = pd.concat(
            [
                joined_df.loc[:, left_gr.columns.to_list()].drop_duplicates(
                    ['Chromosome', 'Start', 'End'], keep='first',
                ).reset_index(drop=True), 
                aggresult,
            ], 
            axis=1,
        )
        matched_rows_gr = pr.PyRanges(matched_rows)

    # merge matched and nonmatched grs
    if unmatched_rows_gr is None:
        result = matched_rows_gr
    else:
        if matched_rows_gr.empty:
            result = unmatched_rows_gr
        else:
            result = pr.concat([unmatched_rows_gr, matched_rows_gr]).sort()

    # return
    if not as_gr:
        result = result.df

    return result


def join(left_gr, right_gr, how='inner', merge=None, find_nearest=False, as_gr=True, verbose=False):
    """Args:
        find_nearest: Match a nearest row from "right_gr" for each
            non-overlapping row of "left_gr". Only relevant when "how" is "left".
    Returns:
        A pyranges.PyRanges object 
    """
    def merge_helper(joined_gr, added_columns, merge):
        # modify joined_gr
        joined_gr = joined_gr.sort()
        joined_gr = add_new_coordinates(joined_gr, 'intersection')
            # now 'new_Start' and 'new_End' columns are added
        joined_gr.isec_length = joined_gr.new_End - joined_gr.new_Start

        # prepare parameters
        added_value_columns = list(set(added_columns).difference({'Start_b', 'End_b'}))
        new_columns = ['Chromosome', 'Start', 'End', 'isec_length'] + added_value_columns
        joined_gr_subset = joined_gr[new_columns]

        # make merged annotation data
        data = _merge_helper_make_data_vectorize(added_value_columns, joined_gr_subset, merge)
        #data = _merge_helper_make_data_iterrows(added_value_columns, joined_gr_subset, merge)

        # prepare final result
        result_df = joined_gr.df.drop_duplicates(
            subset=['Chromosome', 'Start', 'End'], keep='first', ignore_index=True,
        )
        with warnings.catch_warnings(): 
            warnings.simplefilter('ignore', category=FutureWarning)
            result_df.loc[:, added_value_columns] = data
        result = pr.PyRanges(result_df, int64=False)
        return result

    def _merge_helper_make_data_iterrows(added_value_columns, joined_gr_subset, merge):
        def _merge_values_weighted_mean(value_rowlist, lengths):
            return np.average(np.array(value_rowlist), axis=0, weights=lengths)

        def _merge_values_mean(value_rowlist, lengths):
            return np.average(np.array(value_rowlist), axis=0, weights=None)

        def _merge_values_longest_nonmerge(value_rowlist, lengths):
            return value_rowlist[np.argmax(lengths)]

        def _merge_values_longest(value_rowlist, lengths):
            len_val_dict = dict()
            for vals, l in zip(value_rowlist, lengths):
                vals = tuple(vals)
                len_val_dict.setdefault(vals, 0)
                len_val_dict[vals] += l
            return max(len_val_dict.items(), key=operator.itemgetter(1))[0]

        if merge == 'mean':
            values_merger = _merge_values_mean
        elif merge == 'weighted_mean':
            values_merger = _merge_values_weighted_mean
        elif merge == 'longest_nonmerge':
            values_merger = _merge_values_longest_nonmerge
        elif merge == 'longest':
            values_merger = _merge_values_longest

        data = list()
        for key, subiter in itertools.groupby(
            (x[1] for x in joined_gr_subset.df.iterrows()),
            key=(lambda row: tuple(row.iloc[:3])),
        ):
            value_rowlist = list()
            lengths = list()
            for row in subiter:
                value_rowlist.append(row.iloc[4:])
                lengths.append(row.iloc[3])

            if len(value_rowlist) == 1:
                data.append(value_rowlist[0])
            else:
                data.append(values_merger(value_rowlist, lengths))

        return np.array(data)

    # too slow
    def _merge_helper_make_data_pdgroupby(added_value_columns, joined_gr_subset, merge):
        data = list()
        for key, subdf in joined_gr_subset.df.groupby(['Chromosome', 'Start', 'End']):
            if subdf.shape[0] == 1:
                data.append(subdf.iloc[0, 4:])
            else:
                if merge == 'mean':
                    data.append(subdf.iloc[:, 4:].mean())
                elif merge == 'weighted_mean':
                    data.append(np.average(subdf.iloc[:, 4:], axis=0, weights=subdf.iloc[:, 3]))
                elif merge == 'longest_nonmerge':
                    data.append(subdf.iloc[np.argmax(subdf.iloc[:, 3]), 4:])
                elif merge == 'longest':
                    len_val_dict = {
                        tuple(subkey): subsubdf.iloc[:, 3].sum()
                        for subkey, subsubdf in subdf.groupby(added_value_columns)
                    }
                    data.append(
                        max(len_val_dict.items(), key=operator.itemgetter(1))[0]
                    )

        return np.array(data)

    # this is the fastest!
    def _merge_helper_make_data_vectorize(added_value_columns, joined_gr_subset, merge):
        def _merge_annotvals_weighted_mean(annotvals, lengths):
            return np.average(np.array(annotvals), axis=0, weights=lengths)

        def _merge_annotvals_mean(annotvals, lengths):
            #return np.average(np.array(annotvals), axis=0, weights=None)
            return np.nanmean(np.array(annotvals), axis=0)

        def _merge_annotvals_longest_nonmerge(annotvals, lengths):
            return annotvals[np.argmax(lengths)]

        def _merge_annotvals_longest(annotvals, lengths):
            len_val_dict = dict()
            for vals, l in zip(annotvals, lengths):
                len_val_dict.setdefault(vals, 0)
                len_val_dict[vals] += l
            return max(len_val_dict.items(), key=operator.itemgetter(1))[0]

        coord_list = joined_gr_subset.df.apply(lambda row: tuple(row.iloc[:3]), axis=1)
        length_list = joined_gr_subset.df.iloc[:, 3]
        annotval_list = joined_gr_subset.df.apply((lambda row: tuple(row.iloc[4:])), axis=1)

        if merge == 'mean':
            annotvals_merger = _merge_annotvals_mean
        elif merge == 'weighted_mean':
            annotvals_merger = _merge_annotvals_weighted_mean
        elif merge == 'longest_nonmerge':
            annotvals_merger = _merge_annotvals_longest_nonmerge
        elif merge == 'longest':
            annotvals_merger = _merge_annotvals_longest

        data = list()
        for key, subiter in itertools.groupby(
            zip(coord_list, length_list, annotval_list), key=operator.itemgetter(0),
        ):
            coords = list()
            lengths = list()
            annotvals = list()
            for x in subiter:
                coords.append(x[0])
                lengths.append(x[1])
                annotvals.append(x[2])

            if len(coords) == 1:
                data.extend(annotvals)
            else:
                data.append(annotvals_merger(annotvals, lengths))

        return np.array(data)

    # main
    assert merge in ('mean', 'weighted_mean', 'first', 'longest', 'longest_nonmerge', None)
    assert how in ('left', 'inner')
    join_sanity_check(left_gr, right_gr)

    # do join
    #common.funclogger(1)
    with warnings.catch_warnings(): 
        warnings.simplefilter('ignore', category=FutureWarning)
        joined_gr = left_gr.join(right_gr)  
            # this is inner join
            # setting "how='left'" does not work and results in inner join (230209)
    #common.funclogger(2)
    final_added_columns = [x for x in right_gr.columns.to_list() if x not in ('Chromosome', 'Start', 'End')]
    added_columns = ['Start_b', 'End_b'] + final_added_columns
    #if joined_gr.empty:
        #added_columns = pd.Series(joined_gr.columns).drop(left_gr.columns).to_list()
    #else:
        #added_columns = joined_gr.columns.drop(left_gr.columns).to_list()
        # This includes "Start_b" and "End_b"
    #common.funclogger(3)

    # handle unmatched rows
    if how == 'left':
        unmatched_rows_gr = extract_unmatched_rows(
            left_gr, right_gr, joined_gr, find_nearest, final_added_columns,
        )
    elif how == 'inner':
        unmatched_rows_gr = None
    #common.funclogger(4)

    if joined_gr.empty:
        matched_rows_gr = joined_gr
    else:
        # handle matched rows - merge rows with identical (chrom, start, end)
        if merge is None:
            matched_rows_gr = joined_gr
        elif merge == 'first':
            matched_rows_gr = pr.PyRanges(
                joined_gr.df.drop_duplicates(
                    subset=['Chromosome', 'Start', 'End'], keep='first', ignore_index=True,
                )
            )
        elif merge in ('mean', 'weighted_mean', 'longest', 'longest_nonmerge'):
            matched_rows_gr = merge_helper(joined_gr, added_columns, merge)
        #common.funclogger(5)

        # handle matched rows - remove unused columns
        cols_to_drop = set(matched_rows_gr.columns).intersection(
            {'Start_b', 'End_b', 'new_Start', 'new_End', 'isec_length'}
        )
        matched_rows_gr = matched_rows_gr[
            matched_rows_gr.columns.drop(cols_to_drop).to_list()
        ]

    # concat matched and unmatched rows
    if unmatched_rows_gr is None:
        result = matched_rows_gr
    else:
        if matched_rows_gr.empty:
            result = unmatched_rows_gr
        else:
            result = pr.concat([unmatched_rows_gr, matched_rows_gr]).sort()

    # return
    if as_gr:
        return result
    else:
        return result.df


#def join_grs_left_groupby(left_gr, right_gr, merge='mean'):
#    assert merge in ('mean', 'first')
#
#    joined_gr = left_gr.join(right_gr)
#
#    joined_df = joined_gr[list(joined_gr.columns[len(left_gr.columns) + 2:])].df 
#    data = list()
#
#    for key, subdf in joined_df.groupby(['Chromosome', 'Start', 'End']):
#        if merge == 'mean':
#            data.append(list(key) + list(subdf.iloc[:, 3:].mean()))
#        elif merge == 'first':
#            data.append(list(key) + list(subdf.iloc[0:, 3:]))
#
#    averaged_joined_df = pd.DataFrame.from_records(data, columns=list(joined_df.columns))
#    result_df = left_gr.df.join(
#        averaged_joined_df.set_index(['Chromosome', 'Start', 'End'], drop=True), 
#        on=['Chromosome', 'Start', 'End'],
#    )
#
#    return pr.PyRanges(result_df, int64=False)
#
#
#def join_grs_left_itertools(left_gr, right_gr, merge='mean'):
#    assert merge in ('mean', 'first')
#
#    joined_gr = left_gr.join(right_gr)
#    joined_gr = joined_gr.sort()
#
#    added_value_columns = joined_gr.columns.drop(left_gr.columns.to_list() + ['Start_b', 'End_b']).to_list()
#    assert len(added_value_columns) > 0
#    new_columns = ['Chromosome', 'Start', 'End'] + added_value_columns
#    joined_df = joined_gr.df.loc[:, new_columns]
#
#    data = list()
#    for key, subiter in itertools.groupby(
#        (x[1] for x in joined_df.iterrows()),
#        key=(lambda row: tuple(row.iloc[:3])),
#    ):
#        data_items = list()
#        data_items.extend(key)  # chrom, start, end
#
#        if merge == 'mean':
#            rowlist = list()
#            for row in subiter:
#                rowlist.append(row.iloc[3:])
#
#            if len(rowlist) == 1:
#                data_items.extend(rowlist[0])
#            else:
#                data_items.extend(np.array(rowlist).mean(axis=0))
#
#        elif merge == 'first':
#            data_items.extend(next(subiter).iloc[3:])
#
#        data.append(data_items)
#
#    averaged_joined_df = pd.DataFrame.from_records(data, columns=list(joined_df.columns))
#    result_df = left_gr.df.join(
#        averaged_joined_df.set_index(['Chromosome', 'Start', 'End'], drop=True), 
#        on=['Chromosome', 'Start', 'End'],
#    )
#
#    return pr.PyRanges(result_df, int64=False)



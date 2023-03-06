import itertools
import warnings
import operator

import numpy as np
import pandas as pd
import pyranges as pr


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


def inner_join(gr_left, gr_right):
    """For handling grs with lots of columns"""
    new_right_df = gr_right.df.loc[:, ['Chromosome', 'Start', 'End']]
    new_right_df.insert(3, 'index', range(new_right_df.shape[0]))
    new_gr_right = pr.PyRanges(new_right_df)
    joined_gr = gr_left.join(new_gr_right)
    right_values_df = gr_right.df.iloc[joined_gr.index.to_list(), 3:].reset_index(drop=True)
    right_values_df.reset_index(drop=True, inplace=True)
    
    old_columns = gr_left.columns.to_list()
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


def join(gr_left, gr_right, how='inner', merge=None, find_nearest=False, as_gr=True):
    """Args:
        find_nearest: Match a nearest row from "gr_right" for each
            non-overlapping row of "gr_left". Only relevant when "how" is "left".
    Returns:
        A pyranges.PyRanges object 
    """
    def join_left_helper(gr_left, gr_right, joined_gr, find_nearest, added_columns):
        unmatched_rows_gr = gr_left.subtract(joined_gr)
        if unmatched_rows_gr.empty:
            return None

        if find_nearest:
            unmatched_rows_gr = unmatched_rows_gr.nearest(
                gr_right, overlap=False, how=None,
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
    assert isinstance(gr_left, pr.PyRanges)
    assert isinstance(gr_right, pr.PyRanges)
    assert len(set(gr_right.columns).difference(['Chromosome', 'Start', 'End'])) > 0, (
        f'"gr_right" Must have columns other than "Chromosome", "Start", "End".'
    )

    # do join
    with warnings.catch_warnings(): 
        warnings.simplefilter('ignore', category=FutureWarning)
        joined_gr = gr_left.join(gr_right)  
        #joined_gr = inner_join(gr_left, gr_right)
            # this is inner join
            # setting "how='left'" does not work and results in inner join (230209)
    added_columns = joined_gr.columns.drop(gr_left.columns).to_list()
        # This includes "Start_b" and "End_b"

    # handle unmatched rows
    if how == 'left':
        unmatched_rows_gr = join_left_helper(
            gr_left, gr_right, joined_gr, find_nearest, added_columns,
        )
    elif how == 'inner':
        unmatched_rows_gr = None

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
        result = pr.concat([unmatched_rows_gr, matched_rows_gr])
    result = result.sort()

    # return
    if as_gr:
        return result
    else:
        return result.df


#def join_grs_left_groupby(gr_left, gr_right, merge='mean'):
#    assert merge in ('mean', 'first')
#
#    joined_gr = gr_left.join(gr_right)
#
#    joined_df = joined_gr[list(joined_gr.columns[len(gr_left.columns) + 2:])].df 
#    data = list()
#
#    for key, subdf in joined_df.groupby(['Chromosome', 'Start', 'End']):
#        if merge == 'mean':
#            data.append(list(key) + list(subdf.iloc[:, 3:].mean()))
#        elif merge == 'first':
#            data.append(list(key) + list(subdf.iloc[0:, 3:]))
#
#    averaged_joined_df = pd.DataFrame.from_records(data, columns=list(joined_df.columns))
#    result_df = gr_left.df.join(
#        averaged_joined_df.set_index(['Chromosome', 'Start', 'End'], drop=True), 
#        on=['Chromosome', 'Start', 'End'],
#    )
#
#    return pr.PyRanges(result_df, int64=False)
#
#
#def join_grs_left_itertools(gr_left, gr_right, merge='mean'):
#    assert merge in ('mean', 'first')
#
#    joined_gr = gr_left.join(gr_right)
#    joined_gr = joined_gr.sort()
#
#    added_value_columns = joined_gr.columns.drop(gr_left.columns.to_list() + ['Start_b', 'End_b']).to_list()
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
#    result_df = gr_left.df.join(
#        averaged_joined_df.set_index(['Chromosome', 'Start', 'End'], drop=True), 
#        on=['Chromosome', 'Start', 'End'],
#    )
#
#    return pr.PyRanges(result_df, int64=False)



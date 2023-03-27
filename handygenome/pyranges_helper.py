import itertools
import warnings
import operator

import numpy as np
import pandas as pd
import pyranges as pr

import handygenome.common as common
import handygenome.cnv.misc as cnvmisc


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


###############################

def join_preprocess_df(df):
    df = cnvmisc.arg_into_df(df)

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


def compare_coords(left_df, right_df):
    #chrom_compare = (
    #    left_df.Chromosome.to_numpy()[:, np.newaxis] 
    #    == right_df.Chromosome.to_numpy()
    #)
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


def join_singlechrom_dfs(left_df, right_df, how, keep_right_coords, find_nearest, annot_columns, merge):
    print(f'Comparing coordinates')
    compare_result = compare_coords(left_df, right_df)
    left_rowidxs, right_rowidxs = np.where(compare_result)

    print(f'Subsetting right df')
    if keep_right_coords:
        left_subdf = left_df.iloc[left_rowidxs, :]
        right_subdf = right_df.iloc[right_rowidxs, 3:]
        right_subdf['length'] = calc_overlap_lengths(
            left_subdf['Start'].to_numpy(), 
            left_subdf['End'].to_numpy(), 
            right_subdf['Start'].to_numpy(), 
            right_subdf['End'].to_numpy(),
        )
    else:
        right_subdf = right_df.iloc[right_rowidxs, 3:]
    right_subdf.index = left_rowidxs
    
    # do merge
    print(f'Doing merge')
    merged_right_subdf = merge_right_subdf(right_subdf, merge, annot_columns)

    # join with left_df
    print(f'Joining with left df')
    joined_df = left_df.join(merged_right_subdf, how=how)

    # fetch nearest
    if (how == 'left') and find_nearest:
        print(f'Doing fetch nearest')
        joined_df = fetch_nearest(joined_df, right_df, left_rowidxs, annot_columns)

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
            merged_right_subdf = right_subdf.groupby(by=indexes, axis=0, sort=False).mean()
        elif merge in ('weighted_mean', 'longest'):
            if merge == 'weighted_mean':
                def subdf_handler(subdf):
                    weights = subdf['length'].to_numpy()
                    if np.isnan(weights.sum()):
                        return subdf.drop(columns='length').iloc[0, :]
                    else:
                        return subdf.drop(columns='length').apply(lambda x: common.nanaverage(x.to_numpy(), weights=weights))
            elif merge == 'longest':
                def subdf_handler(subdf):
                    idx = np.argmax(subdf['length'])
                    return subdf.drop(columns='length').iloc[idx, :]

            data_gen = (
                subdf_handler(subdf)
                for key, subdf in right_subdf.groupby(by=indexes, axis=0, sort=False)
            )
            merged_right_subdf = pd.DataFrame.from_records(data_gen)

        merged_right_subdf.reset_index(inplace=True, drop=True)
        return merged_right_subdf


def merge_joined_df(
    joined_df, right_df, merge, index_col, how, annot_columns,
):
#    def dedup_by_index_old(df, index_col, keep):
#        df.reset_index(inplace=True, drop=False, names=index_col)
#        df.drop_duplicates(subset=index_col, keep=keep, inplace=True)
#        df.drop(columns=index_col, inplace=True)
#        df.reset_index(inplace=True, drop=True)
#        return df

    if merge is None:
        return joined_df
    elif merge == 'first':
        return dedup_by_index(joined_df, index_col, keep='first')
    elif merge == 'last':
        return dedup_by_index(joined_df, index_col, keep='last')
    else:
        indexes = joined_df.index.to_numpy()
        annot_subdf = joined_df.loc[:, annot_columns]
        nonannot_subdf = joined_df.loc[:, joined_df.columns.difference(annot_columns)]
        nonannot_subdf = dedup_by_index(nonannot_subdf, index_col, keep='first')

        if merge in ('mean',):
            merged_annot_subdf = annot_subdf.groupby(by=indexes, axis=0, sort=False).mean()
        elif merge in ('weighted_mean', 'longest'):
            overlap_legnths = calc_overlap_lengths(
                joined_df['Start'].to_numpy(), 
                joined_df['End'].to_numpy(), 
                joined_df['Start_right'].to_numpy(), 
                joined_df['End_right'].to_numpy(),
            )
            annot_subdf['length'] = overlap_legnths

            if merge == 'weighted_mean':
                def subdf_handler(subdf):
                    weights = subdf['length'].to_numpy()
                    if np.isnan(weights.sum()):
                        return subdf.drop(columns='length').iloc[0, :]
                    else:
                        return subdf.drop(columns='length').apply(lambda x: common.nanaverage(x.to_numpy(), weights=weights))
            elif merge == 'longest':
                def subdf_handler(subdf):
                    idx = np.argmax(subdf['length'])
                    return subdf.drop(columns='length').iloc[idx, :]

            data_gen = (
                subdf_handler(subdf)
                for key, subdf in annot_subdf.groupby(by=indexes, axis=0, sort=False)
            )
            merged_annot_subdf = pd.DataFrame.from_records(data_gen)

        merged_annot_subdf.reset_index(inplace=True, drop=True)
        result = pd.concat([nonannot_subdf, merged_annot_subdf], axis=1)
        return result


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
    for key, subdf in df.groupby('Chromosome'):
        result[key] = subdf.reset_index(drop=True)
    return result


def join_main_newest(left_df, right_df, how, merge, find_nearest, index_col):
    ###
    annot_columns = right_df.columns[3:].to_list()

    ###
    print(f'Grouping dataframes by chromosome')
    left_bychrom = group_df_bychrom(left_df)
    right_bychrom = group_df_bychrom(right_df)

    left_chroms = set(left_bychrom.keys())
    right_chroms = set(right_bychrom.keys())
    common_chroms = left_chroms.intersection(right_chroms)
    leftonly_chroms = left_chroms.difference(right_chroms)
    rightonly_chroms = right_chroms.difference(left_chroms)

    ###
    print(f'Doing join by chromosome')
    keep_right_coords = (merge in ('weighted_mean', 'longest'))
    joined_bychrom = dict()
    for chrom in common_chroms:
        print(f'chrom {chrom}')
        joined_df = join_singlechrom_dfs(
            left_bychrom[chrom], 
            right_bychrom[chrom], 
            how, 
            keep_right_coords,
            find_nearest,
            annot_columns,
            merge,
        )
        #joined_df = merge_joined_df(
        #    joined_df, right_df, merge, index_col, how, annot_columns,
        #)
        joined_bychrom[chrom] = joined_df

    ###
    print(f'Concatenating dataframes of all chromosomes')
    if how == 'inner':
        result = pd.concat(list(joined_bychrom.values()), axis=0)
    elif how == 'left':
        result = pd.concat(
            (
                list(joined_bychrom.values()) 
                + [left_bychrom[chrom] for chrom in leftonly_chroms]
            ), 
            axis=0,
        )

    result.reset_index(inplace=True, drop=True)
    if keep_right_coords:
        result.drop(columns=['Start_right', 'End_right'], inplace=True)

    return result


def join_newest(
    left_df, right_df, how='inner', find_nearest=False, merge=None, as_gr=False,
    sort=False, chromdict=None, 
):
    index_col = '__index'

    # sanity check
    assert merge in (
        'mean', 'weighted_mean', 'longest', #'longest_nonmerge', 
        'first', 'last',
        None,
    )
    assert how in ('left', 'inner')
    join_df_sanitycheck(left_df, right_df, index_col)

    # arg handling
    left_df = join_preprocess_df(left_df)
    right_df = join_preprocess_df(right_df)

    # join
    result = join_main_newest(left_df, right_df, how, merge, find_nearest, index_col)

    return result






















































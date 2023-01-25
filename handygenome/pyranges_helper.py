import itertools
import warnings

import numpy as np
import pandas as pd
import pyranges as pr


def join(gr_left, gr_right, how=None, merge='mean'):
    assert merge in ('mean', 'first', None)

    # do join
    with warnings.catch_warnings(): 
        warnings.simplefilter('ignore', category=FutureWarning)
        joined_gr = gr_left.join(gr_right, how=how)
    joined_gr = joined_gr[joined_gr.columns.drop(['Start_b', 'End_b']).to_list()]
    added_value_columns = joined_gr.columns.drop(gr_left.columns).to_list()
    assert len(added_value_columns) > 0

    if merge is None:
        return joined_gr
    elif merge == 'first':
        joined_df = joined_gr.df.drop_duplicates(subset=['Chromosome', 'Start', 'End'], keep='first', ignore_index=True)
        subdf = joined_df.loc[:, added_value_columns]
        subdf.where(subdf != -1, np.nan, inplace=True)
        joined_df.loc[:, added_value_columns] = subdf
        return pr.PyRanges(joined_df, int64=True)
    elif merge == 'mean':
        # merge duplicated columns
        joined_gr = joined_gr.sort()
        new_columns = ['Chromosome', 'Start', 'End'] + added_value_columns
        joined_subdf = joined_gr.df.loc[:, new_columns]

        data = list()
        for key, subiter in itertools.groupby(
            (x[1] for x in joined_subdf.iterrows()),
            key=(lambda row: tuple(row.iloc[:3])),
        ):
            if merge == 'mean':
                rowlist = list()
                for row in subiter:
                    rowlist.append(row.iloc[3:])

                if len(rowlist) == 1:
                    added_values = rowlist[0]
                else:
                    added_values = np.array(rowlist).mean(axis=0)
            elif merge == 'first':
                added_values = next(subiter).iloc[3:]

            added_values = tuple(
                np.nan if x == -1 else x
                for x in added_values
            )
            data.append(added_values)  # key: chrom, start, end

        result_df = joined_gr.df.drop_duplicates(subset=['Chromosome', 'Start', 'End'], keep='first', ignore_index=True)
        result_df.loc[:, added_value_columns] = np.array(data)
        result = pr.PyRanges(result_df, int64=True)

        return result


def join_grs_left_groupby(gr_left, gr_right, merge='mean'):
    assert merge in ('mean', 'first')

    joined_gr = gr_left.join(gr_right)

    joined_df = joined_gr[list(joined_gr.columns[len(gr_left.columns) + 2:])].df 
    data = list()

    for key, subdf in joined_df.groupby(['Chromosome', 'Start', 'End']):
        if merge == 'mean':
            data.append(list(key) + list(subdf.iloc[:, 3:].mean()))
        elif merge == 'first':
            data.append(list(key) + list(subdf.iloc[0:, 3:]))

    averaged_joined_df = pd.DataFrame.from_records(data, columns=list(joined_df.columns))
    result_df = gr_left.df.join(
        averaged_joined_df.set_index(['Chromosome', 'Start', 'End'], drop=True), 
        on=['Chromosome', 'Start', 'End'],
    )

    return pr.PyRanges(result_df, int64=True)


def join_grs_left_itertools(gr_left, gr_right, merge='mean'):
    assert merge in ('mean', 'first')

    joined_gr = gr_left.join(gr_right)
    joined_gr = joined_gr.sort()

    added_value_columns = joined_gr.columns.drop(gr_left.columns.to_list() + ['Start_b', 'End_b']).to_list()
    assert len(added_value_columns) > 0
    new_columns = ['Chromosome', 'Start', 'End'] + added_value_columns
    joined_df = joined_gr.df.loc[:, new_columns]

    data = list()
    for key, subiter in itertools.groupby(
        (x[1] for x in joined_df.iterrows()),
        key=(lambda row: tuple(row.iloc[:3])),
    ):
        data_items = list()
        data_items.extend(key)  # chrom, start, end

        if merge == 'mean':
            rowlist = list()
            for row in subiter:
                rowlist.append(row.iloc[3:])

            if len(rowlist) == 1:
                data_items.extend(rowlist[0])
            else:
                data_items.extend(np.array(rowlist).mean(axis=0))

        elif merge == 'first':
            data_items.extend(next(subiter).iloc[3:])

        data.append(data_items)

    averaged_joined_df = pd.DataFrame.from_records(data, columns=list(joined_df.columns))
    result_df = gr_left.df.join(
        averaged_joined_df.set_index(['Chromosome', 'Start', 'End'], drop=True), 
        on=['Chromosome', 'Start', 'End'],
    )

    return pr.PyRanges(result_df, int64=True)


def sort_df_by_coord(df, chromdict):
    coord_cols = df.loc[:, ['Chromosome', 'Start', 'End']]
    sorted_rows = sorted(
        coord_cols.iterrows(),
        key=(lambda x: (chromdict.contigs.index(x[1][0]), x[1][1], x[1][2])),
    )
    sorted_rows_indexes = [x[0] for x in sorted_rows]
    return df.iloc[sorted_rows_indexes, :].reset_index(drop=True)


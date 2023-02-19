import itertools
import warnings

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


def add_new_coordinates(joined_gr, how):
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


def join(gr_left, gr_right, how='inner', merge='mean', as_gr=True):
    """Returns:
        A pyranges.PyRanges object 
    """

    assert merge in ('mean', 'weighted_mean', 'first', None)
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
            # this is inner join
            # setting "how='left'" does not work and results in inner join (230209)
    added_value_columns = joined_gr.columns.drop(gr_left.columns).to_list()

    # add non-matching rows of "gr_left"
    if how == 'left':
        unmatched_rows_gr = gr_left.subtract(joined_gr)
        for colname in added_value_columns:
            setattr(
                unmatched_rows_gr, 
                colname, 
                np.repeat(np.nan, unmatched_rows_gr.df.shape[0])
            )
        joined_gr = pr.concat([unmatched_rows_gr, joined_gr]).sort()

    # merge rows with identical (chrom, start, end)
    if merge is None:
        result = joined_gr
    elif merge == 'first':
        joined_df = joined_gr.df.drop_duplicates(subset=['Chromosome', 'Start', 'End'], keep='first', ignore_index=True)
        #subdf = joined_df.loc[:, added_value_columns]
        #subdf.where(subdf != -1, np.nan, inplace=True)
        #joined_df.loc[:, added_value_columns] = subdf
        result = pr.PyRanges(joined_df, int64=False)
    elif merge in ('mean', 'weighted_mean'):
        # merge duplicated columns
        joined_gr = joined_gr.sort()
        new_columns = ['Chromosome', 'Start', 'End'] + added_value_columns
        if merge == 'weighted_mean':
            joined_gr = add_new_coordinates(joined_gr, 'intersection')
                # now 'new_Start' and 'new_End' columns are added
            new_columns += ['new_Start', 'new_End']
        joined_subdf = joined_gr.df.loc[:, new_columns]

        data = list()
        for key, subiter in itertools.groupby(
            (x[1] for x in joined_subdf.iterrows()),
            key=(lambda row: tuple(row.loc[['Chromosome', 'Start', 'End']])),
        ):
            value_rowlist = list()
            lengths = list()
            for row in subiter:
                value_rowlist.append(row.loc[added_value_columns])
                if merge == 'weighted_mean':
                    lengths.append(row.loc['new_End'] - row.loc['new_Start'])

            if len(value_rowlist) == 1:
                added_values = value_rowlist[0]
            else:
                if merge == 'weighted_mean':
                    weights = lengths
                else:
                    weights = None
                added_values = np.average(np.array(value_rowlist), axis=0, weights=weights)

            data.append(added_values)  # key: chrom, start, end

        result_df = joined_gr.df.drop_duplicates(subset=['Chromosome', 'Start', 'End'], keep='first', ignore_index=True)
        with warnings.catch_warnings(): 
            warnings.simplefilter('ignore', category=FutureWarning)
            result_df.loc[:, added_value_columns] = np.array(data)
        result = pr.PyRanges(result_df, int64=False)

    # remove columns 'Start_b', 'End_b', 'new_Start', 'new_End'
    cols_to_drop = set(result.columns).intersection({'Start_b', 'End_b', 'new_Start', 'new_End'})
    result = result[result.columns.drop(cols_to_drop).to_list()]

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



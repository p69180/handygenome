import numpy as np
import pandas as pd
from pandas.util import hash_pandas_object

import handygenome.tools as tools


COMMON_COLUMNS = ['Chromosome', 'Start', 'End']
DEFAULT_DTYPES = {
    'Chromosome': 'string', 
    'Start': int, 
    'End': int,
}


def hash_gr(gr):
    #df_for_hash = pd.concat([getattr(gr, x) for x in gr.columns], axis=1)
    assert not gr.empty

    df_for_hash = gr.df
    df_for_hash.sort_values(COMMON_COLUMNS, inplace=True)
    df_for_hash.reset_index(drop=True, inplace=True)
    return hash_pandas_object(df_for_hash, index=False)


def hash_df(df):
    assert set(COMMON_COLUMNS).issubset(df.columns)
    df_for_hash = df.sort_values(COMMON_COLUMNS, inplace=False)
    df_for_hash.reset_index(drop=True, inplace=True)
    return hash_pandas_object(df_for_hash, index=False)


def check_duplicate_coords(df):
    return df.loc[:, ['Chromosome', 'Start', 'End']].duplicated().any()


def get_coord_groupkey(df, chromdict):
    """Does not sort before grouping, like itertools.groupby"""
    chrom_indexes = df['Chromosome'].apply(chromdict.contigs.index)
    coord_arr = np.stack(
        [chrom_indexes.to_numpy(), df['Start'].to_numpy(), df['End'].to_numpy()], 
        axis=1,
    )
    _, counts, groupkey = tools.array_grouper(coord_arr, omit_values=True)
    return counts, groupkey


def check_interval_overlap(start0s, end0s):
    start0s = np.asarray(start0s)
    end0s = np.asarray(end0s)
    assert (start0s.ndim == 1) and (end0s.ndim == 1)

    idxs = np.argsort(start0s)
    start0s_sort = start0s[idxs]
    end0s_sort = end0s[idxs]

    return not (start0s_sort[1:] >= end0s_sort[:-1]).all()


def compact_coords(start0s, end0s):
    # sanity check
    if check_interval_overlap(start0s, end0s):
        raise Exception(f'Intervals must not overlap.')

    intervals = list(zip(start0s, end0s))
    #encoder = dict()  # noncompact (start, end) -> compact single pos
    #decoder = dict()  # compact single pos -> noncompact (start, end)
    compact_pos0_list = np.argsort(start0s)
    decoder = dict(zip(compact_pos0_list, intervals))
    #for compact_pos0, intv in zip(compact_pos0_list, intervals):
        #encoder[intv] = compact_pos0
    #    decoder[compact_pos0] = intv
        #compact_pos0_list.append(compact_pos0)

    #return encoder, decoder
    return compact_pos0_list, decoder


def uncompact_coords(compact_start0s, compact_end0s, decoder):
    start0s = list()
    end0s = list()
    for compact_s0, compact_e0 in zip(compact_start0s, compact_end0s):
        s0 = decoder[compact_s0][0]
        start0s.append(s0)

        last_pos0 = compact_e0 - 1
        e0 = decoder[last_pos0][1]
        end0s.append(e0)
        
    return start0s, end0s



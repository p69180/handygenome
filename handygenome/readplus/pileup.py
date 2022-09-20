import re
import collections

import numpy as np
import pandas as pd

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
readhandler = importlib.import_module(".".join([top_package_name, "readplus", "readhandler"]))


class Pileup:
    def __init__(self, df=None, uids=None, chrom=None):
        self.df = df
        self.uids = uids
        self.chrom = chrom
        self.split_cigars = None

    def get_max_vaf(self, col_idx):
        col = self.df.iloc[:, col_idx]
        return get_max_vaf_from_col(col)

    def get_seq_specs(self):
        seq_specs = dict()
        for uid, row in self.df.iterrows():
            seq_spec = {
                'seq': ''.join(
                    re.sub('[()]', '', x) for x in row 
                    if x not in ('*', None)
                ),
                'match_length': sum(x is not None for x in row),
                'left_filled': row.iloc[0] is not None,
                'right_filled': row.iloc[-1] is not None,
                'uid': uid,
            }
            seq_specs[seq_spec['uid']] = seq_spec

        return seq_specs


def get_max_vaf_from_col(col):
    counts = collections.Counter(x for x in col if x is not None)
    total = sum(counts.values())
    if total == 0:
        return None
    else:
        return max(counts.values()) / total


def pileup_readfilter(read):
    """True if a bad read"""
    return (
        readhandler.check_bad_read(read) or
        read.is_supplementary or
        read.is_secondary
    )


def get_pileup(bam, chrom, start0, end0, as_array=False, truncate=True, with_reads=False):
    """Returns:
        A numpy.ndarray object if "as_array" argument is True
        A Pileup object if "as_array" argument is False
    """

    read_dict = {
        readhandler.get_uid(x): x 
        for x in readhandler.get_fetch(bam, chrom, start0, end0, filter_fun=pileup_readfilter)
    }
    readlist = tuple(iter(read_dict.values()))
    pileup, pileup_range = get_pileup_from_reads(readlist, as_array=as_array, return_range=True)
    
    if truncate:
        if as_array:
            start0_idx = pileup_range.index(start0)
            end0_idx = pileup_range.index(end0)
            pileup = pileup[:, start0_idx:end0_idx]
        else:
            pileup.df = pileup.df.loc[:, list(range(start0, end0))]

    if with_reads:
        return pileup, read_dict
    else:
        return pileup


def get_pileup_from_reads(readlist, as_array=False, pileup_range=None, return_range=False, del_char='*'):
    def raise_err(read):
        raise Exception(f"Unexpected cigar pattern:\n{read.to_string()}")

    def cigar_sanitycheck(read):
        """Assumes M.* or [IS]M.*"""
        error = False
        first_cigarop = read.cigartuples[0][0]
        if first_cigarop != 0:
            if first_cigarop not in (1, 4):
                error = True
            else:
                if read.cigartuples[1][0] != 0:
                    error = True
        if error:
            raise_err(read)

    def get_pileup_range(readlist):
        start0s = list()
        end0s = list()
        for read in readlist:
            start0s.append(read.reference_start)
            end0s.append(read.reference_end)
        start0 = min(start0s)
        end0 = max(end0s)

        return range(start0, end0)

    def check_overlaps_pileup_range(read, pileup_range):
        if read.cigartuples[-1][0] in readhandler.CIGAROPS_QUERYONLY:
            overlaps_upstream = (read.reference_end >= pileup_range.start)
        else:
            overlaps_upstream = (read.reference_end > pileup_range.start)

        overlaps_downstream = (read.reference_start < pileup_range.stop)

        return overlaps_upstream and overlaps_downstream

    def handle_leading_insclip(read):
        leading_insclip_len = 0
        for cigartup_idx, (cigarop, cigarlen) in enumerate(read.cigartuples):
            if cigarop in (1, 4):
                leading_insclip_len += cigarlen
            else:
                break
                
        if leading_insclip_len == 0:
            leading_insseq = None
        else:
            leading_insseq = read.query_sequence[:leading_insclip_len]
            
        query_idx = leading_insclip_len

        return leading_insseq, cigartup_idx, query_idx

    def handle_last_match(
        read, arr, cigarlen, query_idx, arr_col_idx, arr_row_idx, cigartup_idx
    ):
        sl_query = slice(query_idx, query_idx + cigarlen)
        seq_buffer = tuple(read.query_sequence[sl_query])
        sl_arr_col = slice(arr_col_idx, arr_col_idx + cigarlen)
        arr[arr_row_idx, sl_arr_col] = seq_buffer
        cigartup_idx += 1

        return cigartup_idx

    def add_matches_before_last(
        read, arr, cigarlen, query_idx, arr_col_idx, arr_row_idx
    ):
        sl_query = slice(query_idx, query_idx + cigarlen - 1)
        seq = tuple(read.query_sequence[sl_query])
        query_idx += cigarlen - 1

        sl_arr_col = slice(arr_col_idx, arr_col_idx + cigarlen - 1)
        arr[arr_row_idx, sl_arr_col] = seq
        arr_col_idx += cigarlen - 1

        return query_idx, arr_col_idx

    def handle_del_after_match(
        read, arr, query_idx, arr_col_idx, arr_row_idx, next_cigarlen
    ):
        seq = read.query_sequence[query_idx]
        query_idx += 1
        arr[arr_row_idx, arr_col_idx] = seq
        arr_col_idx += 1

        sl_arr_col = slice(arr_col_idx, arr_col_idx + next_cigarlen)
        arr[arr_row_idx, sl_arr_col] = del_char
        arr_col_idx += next_cigarlen

        return query_idx, arr_col_idx

    def handle_insclip_after_match(
        read, arr, query_idx, arr_col_idx, arr_row_idx, next_cigarlen
    ):
        match_seq = read.query_sequence[query_idx]
        query_idx += 1
        insseq = read.query_sequence[query_idx : (query_idx + next_cigarlen)]
        query_idx += next_cigarlen
        seq = match_seq + f"({insseq})"

        arr[arr_row_idx, arr_col_idx] = seq
        arr_col_idx += 1

        return query_idx, arr_col_idx

    def add_seqs(read, arr, cigartup_idx, query_idx, arr_row_idx, arr_col_idx):
        while True:
            if cigartup_idx == len(read.cigartuples):
                break

            cigarop, cigarlen = read.cigartuples[cigartup_idx]
            if cigarop != 0:
                raise_err(read)

            if cigartup_idx == len(read.cigartuples) - 1:
                cigartup_idx = handle_last_match(
                    read,
                    arr,
                    cigarlen,
                    query_idx,
                    arr_col_idx,
                    arr_row_idx,
                    cigartup_idx,
                )
                continue
            else:
                if cigarlen > 1:
                    query_idx, arr_col_idx = add_matches_before_last(
                        read, arr, cigarlen, query_idx, arr_col_idx, arr_row_idx
                    )

                next_cigarop, next_cigarlen = read.cigartuples[cigartup_idx + 1]
                if next_cigarop == 2:  # deletion
                    query_idx, arr_col_idx = handle_del_after_match(
                        read, arr, query_idx, arr_col_idx, arr_row_idx, next_cigarlen
                    )
                elif next_cigarop in (1, 4):
                    query_idx, arr_col_idx = handle_insclip_after_match(
                        read, arr, query_idx, arr_col_idx, arr_row_idx, next_cigarlen
                    )

                cigartup_idx += 2
                continue

    def add_leading_insseq(read, leading_insseq, arr, arr_row_idx, arr_col_idx_init):
        if leading_insseq is not None:
            old_seq = arr[arr_row_idx, arr_col_idx_init]
            new_seq = f"(l+{leading_insseq})" + old_seq
            arr[arr_row_idx, arr_col_idx_init] = f"({leading_insseq})" + old_seq

    # main
    chrom = readlist[0].reference_name
    if pileup_range is None:
        pileup_range = get_pileup_range(readlist)

    ncol = len(pileup_range)
    nrow = len(readlist)
    #arr = np.empty((nrow, ncol), dtype=object)
    arr = np.full((nrow, ncol), None, dtype=object)

    for arr_row_idx, read in enumerate(
        read for read in readlist 
        if check_overlaps_pileup_range(read, pileup_range)
    ):
        try:
            arr_col_idx_init = read.reference_start - pileup_range.start
            arr_col_idx = arr_col_idx_init
            leading_insseq, cigartup_idx, query_idx = handle_leading_insclip(read)
            add_seqs(read, arr, cigartup_idx, query_idx, arr_row_idx, arr_col_idx)
            add_leading_insseq(read, leading_insseq, arr, arr_row_idx, arr_col_idx_init)
        except Exception as exc:
            try:
                cigar_sanitycheck(read)
            except Exception as exc_cigarpattern:
                raise exc_cigarpattern from exc
            else:
                raise exc

    if as_array:
        if return_range:
            return (arr, pileup_range)
        else:
            return arr
    else:
        uids = [readhandler.get_uid(read) for read in readlist]
        if len(uids) != len(set(uids)):
            reads_as_string = '\n'.join(x.to_string() for x in readlist)
            raise Exception(f'Duplicate read uid. Input reads are as following:\n{reads_as_string}')

        df = pd.DataFrame(arr, columns=list(pileup_range))
        pileup = Pileup(df=df, uids=uids, chrom=chrom)
        if return_range:
            return (pileup, pileup_range)
        else:
            return pileup

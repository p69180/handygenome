import re
import itertools
import collections

import numpy as np
import pandas as pd

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
readhandler = importlib.import_module(
    ".".join([top_package_name, "readplus", "readhandler"])
)
libcigar = importlib.import_module(".".join([top_package_name, "readplus", "cigar"]))


DEL_VALUE = "*"
EMPTY_VALUE = ""
DEFAULT_ACTIVE_THRESHOLD = 0.9
DEFAULT_EXTEND_FETCHEDREADS_BY = 300
DEFAULT_EXTEND_PILEUP_BY = 30


class FetchedReads:
    def __init__(self, bam, chrom, start0, end0, readfilter=None):
        # set basic attrs
        self.chrom = chrom
        self.fetch_start0 = start0
        self.fetch_end0 = end0
        # set readranges, list, and dict
        if readfilter is None:
            readfilter = pileup_readfilter
        self.readranges = list()
        self.uids = list()
        self.list = list()
        for read in readhandler.get_fetch(
            bam, chrom, start0, end0, filter_fun=readfilter
        ):
            self.readranges.append((read.reference_start, read.reference_end))
            self.uids.append(readhandler.get_uid(read))
            self.list.append(read)
        self.dict = dict(zip(self.uids, self.list))
        # check uid uniqueness
        if len(set(self.uids)) != len(self.uids):
            reads_string = "\n".join(x.to_string() for x in self.list)
            uids_string = "\n".join(x for x in self.uids)
            raise Exception(
                f"Duplicate read uids.\n"
                f"Read list:\n{reads_string}.\n"
                f"uid list:\n{uids_string}."
            )
        # set readspan_start0, readspan_end0
        if len(self.readranges) == 0:
            self.readspan_start0 = None
            self.readspan_end0 = None
        else:
            self.readspan_start0 = min(x[0] for x in self.readranges)
            self.readspan_end0 = max(x[1] for x in self.readranges)

    def _fetch_forward(self, fetch_start0, fetch_end0, with_uid=False):
        read_end_ok = False
        for read_spanning_range, uid, read in zip(
            self.readranges, self.uids, self.list
        ):
            if not read_end_ok:
                if read_spanning_range[1] > fetch_start0:
                    read_end_ok = True
                else:
                    continue

            if read_spanning_range[0] < fetch_end0:
                if with_uid:
                    yield read, uid
                else:
                    yield read
            else:
                return

    def _fetch_reverse(self, fetch_start0, fetch_end0, with_uid=False):
        read_start_ok = False
        for read_spanning_range, uid, read in zip(
            reversed(self.readranges), reversed(self.uids), reversed(self.list)
        ):
            if not read_start_ok:
                if read_spanning_range[0] < fetch_end0:
                    read_start_ok = True
                else:
                    continue

            if read_spanning_range[1] > fetch_start0:
                if with_uid:
                    yield read, uid
                else:
                    yield read
            else:
                return

    def fetch(self, fetch_start0=None, fetch_end0=None, with_uid=False):
        if len(self.list) == 0:
            return iter(())
        else:
            if fetch_start0 is None:
                fetch_start0 = self.readspan_start0
            else:
                fetch_start0 = max(fetch_start0, self.readspan_start0)

            if fetch_end0 is None:
                fetch_end0 = self.readspan_end0
            else:
                fetch_end0 = min(fetch_end0, self.readspan_end0)

            if fetch_start0 >= fetch_end0:
                return iter(())
            else:
                dist_from_start = fetch_start0 - self.readspan_start0
                dist_from_end = self.readspan_end0 - fetch_end0

                if dist_from_start <= dist_from_end:
                    return self._fetch_forward(
                        fetch_start0, fetch_end0, with_uid=with_uid
                    )
                else:
                    return reversed(
                        tuple(
                            self._fetch_reverse(
                                fetch_start0, fetch_end0, with_uid=with_uid
                            )
                        )
                    )
                    # return self._fetch_reverse(fetch_start0, fetch_end0)


class Pileup:
    pat_insclip = re.compile("(\(.+\))?([^()]+)(\(.+\))?")
    pat_parenthesis = re.compile("[()]")
    seq_spec_exluded_vals = {DEL_VALUE, EMPTY_VALUE}

    def __init__(
        self, df, chrom, fasta, bam, active_threshold=DEFAULT_ACTIVE_THRESHOLD
    ):
        self.df = df
        self.chrom = chrom
        self.fasta = fasta
        self.bam = bam
        self.active_threshold = active_threshold
        self.fetchedreads_list = list()

        #self.seq_specs = dict()
        # self.uids = uids
        # self.split_cigars = None

    @property
    def start0(self):
        return self.df.columns[0]

    @property
    def end0(self):
        return self.df.columns[-1] + 1

    @property
    def range(self):
        return range(self.start0, self.end0)

    def get_max_vaf_colindex(self, col_idx, with_counts=False):
        col = self.df.iloc[:, col_idx]
        return get_max_vaf_from_pileupcol(col, with_counts=with_counts)

    def get_max_vaf_colpos0(self, col_pos0, with_counts=False):
        col = self.df.loc[:, col_pos0]
        return get_max_vaf_from_pileupcol(col, with_counts=with_counts)

    def subset(self, start0, end0):
        self._coord_arg_sanitycheck(start0, end0)
        # subset dataframe
        new_df = self.df.loc[:, start0:(end0 - 1)]
        # init result
        result = self.__class__(
            df=new_df, chrom=self.chrom, fasta=self.fasta, bam=self.bam,
            active_threshold=self.active_threshold,
        )
        # fetchedreads_list
        result.fetchedreads_list = self.fetchedreads_list
        # active_info
        active_info_idx_start = start0 - self.start0
        active_info_idx_end = end0 - self.start0
        result._active_info = self._active_info[active_info_idx_start:active_info_idx_end]

        return result

    ##################

    def _coord_arg_sanitycheck(self, start0, end0):
        if start0 is not None:
            if start0 < self.start0:
                raise Exception('Input "start0" argument is out of pileup range.')
        if end0 is not None:
            if end0 > self.end0:
                raise Exception('Input "end0" argument is out of pileup range.')

    def iter_active_info(self, start0, end0, reverse=False):
        self._coord_arg_sanitycheck(start0, end0)

        ref_seq = self.fasta.fetch(self.chrom, start0, end0)
        pos0_range = range(start0, end0)
        if reverse:
            ref_seq = ref_seq[::-1]
            pos0_range = range(pos0_range)

        for ref_base, pos0 in zip(ref_seq, pos0_range):
            max_allele, max_vaf = self.get_max_vaf_colpos0(pos0, with_counts=False)
            if max_vaf >= self.active_threshold:
                is_active = not (max_allele == ref_base)
            else:
                is_active = True
            yield is_active

    def _set_active_info(self, start0=None, end0=None):
        self._coord_arg_sanitycheck(start0, end0)

        if not hasattr(self, '_active_info'):
            setattr(self, '_active_info', list())

        if start0 is None:
            start0 = self.start0
        if end0 is None:
            end0 = self.end0
        active_info_idx_start = start0 - self.start0
        active_info_idx_end = end0 - self.start0

        self._active_info[
            active_info_idx_start:active_info_idx_end
        ] = self.iter_active_info(start0, end0, reverse=False)

    def get_active_info(self, start0=None, end0=None):
        self._coord_arg_sanitycheck(start0, end0)

        if start0 is None:
            start0 = self.start0
        if end0 is None:
            end0 = self.end0
        active_info_idx_start = start0 - self.start0
        active_info_idx_end = end0 - self.start0

        return pd.Series(
            self._active_info[active_info_idx_start:active_info_idx_end], 
            index=self.df.columns[active_info_idx_start:active_info_idx_end],
        )

    def get_active_positions(self):
        return tuple(
            itertools.compress(self.df.columns, self._active_info)
        )

    ### extend ###
    def extend_rightward(self, width, extend_fetchedreads_by=DEFAULT_EXTEND_FETCHEDREADS_BY):
        # keep original starts and ends
        border_col_idx_left = self.df.shape[1] - 1
        # extend fetchedreads if needed
        self._extend_helper_prepare_fetchedreads_rightward(
            width, extend_fetchedreads_by
        )
        # make new pileup
        new_pileup = self._extend_helper_make_new_pileup(
            pileup_coord=(self.chrom, self.end0, self.end0 + width),
            fetchedreads=self.fetchedreads_list[-1],
        )
        # join
        self.df = self.df.join(new_pileup.df, how="outer")
        self.df[self.df.isnull()] = EMPTY_VALUE  # turn NaN into EMPTY_VALUE
        self._active_info.extend(new_pileup._active_info)
        # handle insclips facing each other
        self._extend_helper_handle_facing_insclip(border_col_idx_left)

    def extend_leftward(self, width, extend_fetchedreads_by=DEFAULT_EXTEND_FETCHEDREADS_BY):
        # keep original starts and ends
        border_col_idx_left = width - 1
        # extend fetchedreads if needed
        self._extend_helper_prepare_fetchedreads_leftward(width, extend_fetchedreads_by)
        # make new pileup
        new_pileup = self._extend_helper_make_new_pileup(
            pileup_coord=(self.chrom, self.start0 - width, self.start0),
            fetchedreads=self.fetchedreads_list[0],
        )
        # join
        self.df = new_pileup.df.join(self.df, how="outer")
        self.df[self.df.isnull()] = EMPTY_VALUE  # turn NaN into EMPTY_VALUE
        self._active_info = new_pileup._active_info + self._active_info
        # handle insclips facing each other
        self._extend_helper_handle_facing_insclip(border_col_idx_left)

    def _extend_helper_handle_facing_insclip(self, border_col_idx_left):
        border_col_idx_right = border_col_idx_left + 1
        border_col_left = self.df.iloc[:, border_col_idx_left]
        border_col_right = self.df.iloc[:, border_col_idx_right]
        for row_idx, (val_left, val_right) in enumerate(
            zip(border_col_left, border_col_right)
        ):
            if val_right.startswith("(") and val_left.endswith(")"):
                mat_left = self.__class__.pat_insclip.fullmatch(val_left)
                mat_right = self.__class__.pat_insclip.fullmatch(val_right)
                if mat_left.group(3) != mat_right.group(1):
                    raise Exception(
                        f"Insclip seqs of adjacent entries are different.\n{self.df}"
                    )
                self.df.iloc[
                    row_idx, border_col_idx_right
                ] = self.__class__.pat_insclip.sub("\\2\\3", "", val_right)

    def _extend_helper_prepare_fetchedreads_rightward(
        self, width, extend_fetchedreads_by
    ):
        if self.end0 + width > self.fetchedreads_list[-1].fetch_end0:
            new_fetchedreads = FetchedReads(
                bam=self.bam,
                chrom=self.chrom,
                start0=self.end0,
                end0=(self.end0 + extend_fetchedreads_by),
            )
            self.fetchedreads_list.append(new_fetchedreads)

    def _extend_helper_prepare_fetchedreads_leftward(
        self, width, extend_fetchedreads_by
    ):
        if self.start0 - width < self.fetchedreads_list[0].fetch_start0:
            new_fetchedreads = FetchedReads(
                bam=self.bam,
                chrom=self.chrom,
                start0=(self.start0 - extend_fetchedreads_by),
                end0=self.start0,
            )
            self.fetchedreads_list.insert(0, new_fetchedreads)

    def _extend_helper_make_new_pileup(self, pileup_coord, fetchedreads):
        return get_pileup(
            bam=self.bam,
            pileup_coord=pileup_coord,
            fetchedreads=fetchedreads,
            fasta=self.fasta,
            active_threshold=self.active_threshold,
            truncate=True,
            as_array=False,
            return_range=False,
        )
    ##############

    ### seq specs ###
    def get_seq_spec(self, row_idx):
        row = self.df.iloc[row_idx, :]
        return {
            "seq": "".join(
                self.__class__.pat_parenthesis.sub("", x) 
                for x in row if x not in self.__class__.seq_spec_exluded_vals
            ),
            "match_length": sum(x != EMPTY_VALUE for x in row),
            "left_filled": (row.iloc[0] != EMPTY_VALUE),
            "right_filled": (row.iloc[-1] != EMPTY_VALUE),
            "uid": row.name,
        }

    def set_seq_specs(self):
        self.seq_specs = dict()
        for row_idx in range(self.df.shape[0]):
            seq_spec = self.get_seq_spec(row_idx)
            self.seq_specs[seq_spec['uid']] = seq_spec
            
            
def group_seq_specs(seq_specs):
    """Input "seq_spec"s are grouped by whether one is a substring of another."""
    
    def sortkey(seq_spec):
        return -1 * seq_spec['match_length']
    
    def matcher(query, target):
        if query['left_filled']:
            if query['right_filled']:
                if target['left_filled'] and target['right_filled']:
                    return query['seq'] in target['seq']
                else:
                    return False
            else:
                if target['left_filled']:
                    return query['seq'] == target['seq'][:len(query['seq'])]
                else:
                    return False
        else:
            if query['right_filled']:
                if target['right_filled']:
                    return query['seq'] == target['seq'][-len(query['seq']):]
                else:
                    return False
            else:
                return query['seq'] in target['seq']                
        
    # main
    seq_spec_groups = dict()
    for query in sorted(seq_specs.values(), key=sortkey):
        superseq_candidates = list()
        for uid in seq_spec_groups.keys():
            target = seq_specs[uid]
            if matcher(query, target):
                superseq_candidates.append(uid)
                
        if len(superseq_candidates) == 0:                
            seq_spec_groups[query['uid']] = set()
            seq_spec_groups[query['uid']].add(query['uid'])
        else:
            superseq_uid = max(superseq_candidates, key=(lambda uid: seq_specs[uid]['match_length']))
            seq_spec_groups[superseq_uid].add(query['uid'])
            
    return seq_spec_groups


def get_max_vaf_from_pileupcol(col, with_counts=False, empty_value=EMPTY_VALUE):
    counts = collections.Counter(x for x in col if x != empty_value)
    total = sum(counts.values())
    if total == 0:
        raise Exception(f"Pileup column allele count sum is 0.")

    max_allele, max_count = max(counts.items(), key=(lambda x: x[1]))
    max_vaf = max_count / total

    if with_counts:
        return max_allele, max_vaf, counts
    else:
        return max_allele, max_vaf


def pileup_readfilter(read):
    return (
        readhandler.check_bad_read(read) or read.is_supplementary or read.is_secondary
    )


def get_pileup(
    bam,
    pileup_coord=None,
    fetchedreads=None,
    fasta=None,
    active_threshold=DEFAULT_ACTIVE_THRESHOLD,
    truncate=True,
    as_array=False,
    return_range=False,
    del_value=DEL_VALUE,
    empty_value=EMPTY_VALUE,
):
    """Args:
    pileup_coord: A tuple composed of: (chrom, start0, end0)
    """

    def sanity_check(pileup_coord, bam, fetchedreads):
        if fetchedreads is None:
            if pileup_coord is None or bam is None:
                raise Exception(
                    f'If "fetchedreads" is not set, "pileup_coord" and "bam" must be set.'
                )
        if pileup_coord is None:
            if fetchedreads is None:
                raise Exception(
                    f'If "pileup_coord" is not set, "fetchedreads" must be set.'
                )
        if (pileup_coord is not None) and (fetchedreads is not None):
            if (
                pileup_coord[1] < fetchedreads.fetch_start0
                or pileup_coord[2] > fetchedreads.fetch_end0
            ):
                raise Exception(
                    f'Position range specified by "pileup_coord" is not included within fetch range of "fetchedreads".'
                )
            if pileup_coord[0] != fetchedreads.chrom:
                raise Exception(
                    f'chrom is different between "pileup_coord" and "fetchedreads".'
                )
        if (not as_array) and (fasta is None):
            raise Exception(f'If "as_array" is False, "fasta" must be set.')

    def set_params(pileup_coord, bam, fetchedreads):
        if fetchedreads is None:
            chrom, start0, end0 = pileup_coord
            fetchedreads = FetchedReads(bam, chrom, start0, end0)
            pileup_range = range(start0, end0)
        elif pileup_coord is None:
            pileup_range = range(fetchedreads.fetch_start0, fetchedreads.fetch_end0)
        else:
            pileup_range = range(pileup_coord[1], pileup_coord[2])

        return pileup_range, fetchedreads

    def make_readlist(fetchedreads, pileup_range):
        readlist = list()
        start_list = list()
        end_list = list()
        uid_list = list()
        for read, uid in fetchedreads.fetch(
            pileup_range.start, pileup_range.stop, with_uid=True
        ):
            readlist.append(read)
            start_list.append(read.reference_start)
            end_list.append(read.reference_end)
            uid_list.append(uid)
        tmp_pileup_range = range(min(start_list), max(end_list))

        return readlist, tmp_pileup_range, uid_list

    def raise_cigarpattern_error(read):
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
            raise_cigarpattern_error(read)

    def add_seqs(read, arr, arr_row_idx, initial_arr_col_idx, del_value):
        def handle_leading_insclip(read):
            leading_insclip_len = 0
            cigartup_idx = 0
            for cigarop, cigarlen in libcigar.iter_leading_queryonly(read.cigartuples):
                leading_insclip_len += cigarlen
                cigartup_idx += 1

            if leading_insclip_len == 0:
                leading_insseq = None
            else:
                leading_insseq = read.query_sequence[:leading_insclip_len]
            query_idx = leading_insclip_len

            return leading_insseq, cigartup_idx, query_idx

        def handle_queryonly_only_case(read, arr, arr_col_idx, arr_row_idx):
            arr[arr_row_idx, arr_col_idx] = f"({read.query_sequence})"

        def handle_targetonly_after_match(
            read,
            arr,
            query_idx,
            arr_col_idx,
            arr_row_idx,
            trailing_nonM_cigarlen_sum,
            del_value,
        ):
            arr[arr_row_idx, arr_col_idx] = read.query_sequence[
                query_idx
            ]  # query_idx indicates the last base of M
            arr[
                arr_row_idx,
                (arr_col_idx + 1) : (arr_col_idx + 1 + trailing_nonM_cigarlen_sum),
            ] = del_value

        def handle_queryonly_after_match(
            read,
            arr,
            query_idx,
            arr_col_idx,
            arr_row_idx,
            trailing_nonM_cigarlen_sum,
        ):
            last_match_base = read.query_sequence[query_idx]
            insseq = read.query_sequence[
                (query_idx + 1) : (query_idx + 1 + trailing_nonM_cigarlen_sum)
            ]
            arr[arr_row_idx, arr_col_idx] = f"{last_match_base}({insseq})"

        def handle_nonlast_match(
            read,
            arr,
            cigarlen,
            query_idx,
            arr_col_idx,
            arr_row_idx,
            cigartup_idx,
        ):
            # adds all match seqs except the last one
            # query_idx and arr_col_idx updated
            if cigarlen > 1:
                sl_query = slice(query_idx, query_idx + cigarlen - 1)
                sl_arr_col = slice(arr_col_idx, arr_col_idx + cigarlen - 1)
                arr[arr_row_idx, sl_arr_col] = tuple(read.query_sequence[sl_query])
                query_idx += cigarlen - 1
                arr_col_idx += cigarlen - 1
            # get all subsequent non-M cigar units
            # cigartup_idx updated
            trailing_nonMs = list()
            while True:
                cigartup_idx += 1
                if cigartup_idx >= len(read.cigartuples):
                    break
                else:
                    next_cigartup = read.cigartuples[cigartup_idx]
                    if next_cigartup[0] != 0:
                        trailing_nonMs.append(next_cigartup)
                    else:
                        break
            # add trailing non-M cigar units
            # query_idx and arr_col_idx updated
            trailing_nonM_cigarops = set(x[0] for x in trailing_nonMs)
            trailing_nonM_cigarlen_sum = sum(x[1] for x in trailing_nonMs)
            if trailing_nonM_cigarops.issubset(libcigar.CIGAROPS_TARGETONLY):  # D, N
                handle_targetonly_after_match(
                    read,
                    arr,
                    query_idx,
                    arr_col_idx,
                    arr_row_idx,
                    trailing_nonM_cigarlen_sum,
                    del_value,
                )
                query_idx += 1
                arr_col_idx += trailing_nonM_cigarlen_sum + 1
            elif trailing_nonM_cigarops.issubset(libcigar.CIGAROPS_QUERYONLY):  # I, S
                handle_queryonly_after_match(
                    read,
                    arr,
                    query_idx,
                    arr_col_idx,
                    arr_row_idx,
                    trailing_nonM_cigarlen_sum,
                )
                query_idx += trailing_nonM_cigarlen_sum + 1
                arr_col_idx += 1
            else:
                raise Exception(
                    f"Consecutive Non-M cigar units are composed of both target-only and query-only ones:\n{read.to_string()}"
                )

            return cigartup_idx, query_idx, arr_col_idx

        def handle_last_match(
            read,
            arr,
            cigarlen,
            query_idx,
            arr_col_idx,
            arr_row_idx,
        ):
            # adds all match seqs to array
            sl_query = slice(query_idx, query_idx + cigarlen)
            sl_arr_col = slice(arr_col_idx, arr_col_idx + cigarlen)
            arr[arr_row_idx, sl_arr_col] = tuple(read.query_sequence[sl_query])

        # main
        """1) First, leading query-only cigar units are collected.
        2) Then, cigartuple is iterated on. Initial cigartup_idx is at the 
            first cigar unit which is not a query-only one. This cigar unit
            is assumed to be M. (If not, raises an exception)
        3) What is done in a loop cycle:
            The cigar unit which cigartup_idx indicates (which should be M)
            and all subsequent non-M cigarop are treated at single loop cycle.
            The subsequent non-M cigar units are assumed to be composed of only
            query-only ones or target-only ones. (If not, raises an exception)
        """

        # set mutable arr_col_idx. initial one is kept separately.
        arr_col_idx = initial_arr_col_idx

        # handles leading query-only cigar units. cigartup_idx and query_idx are initialized
        leading_insseq, cigartup_idx, query_idx = handle_leading_insclip(read)

        if cigartup_idx == len(read.cigartuples):
            # cigartup is only composed of query-only operations (no M)
            arr[arr_row_idx, arr_col_idx] = "(" + read.query_sequence + ")"
        else:
            if read.cigartuples[cigartup_idx][0] != 0:
                raise Exception(
                    f"The first cigar operation after stripping leading I and S is not M:\n{read.to_string()}"
                )

            # begins loop
            while True:
                if cigartup_idx > len(read.cigartuples) - 1:
                    raise Exception(
                        f'"cigartup_idx" became greater than "len(read.cigartuples) - 1" while looping.'
                    )

                cigarop, cigarlen = read.cigartuples[cigartup_idx]
                assert (
                    cigarop == 0
                ), f'The cigar unit indicated by "cigarop_idx" at the beginning of a loop cycle is not M.'

                if (
                    cigartup_idx == len(read.cigartuples) - 1
                ):  # current cigartup is the last one
                    handle_last_match(
                        read,
                        arr,
                        cigarlen,
                        query_idx,
                        arr_col_idx,
                        arr_row_idx,
                    )
                    break
                else:
                    cigartup_idx, query_idx, arr_col_idx = handle_nonlast_match(
                        read,
                        arr,
                        cigarlen,
                        query_idx,
                        arr_col_idx,
                        arr_row_idx,
                        cigartup_idx,
                    )
                    if cigartup_idx == len(read.cigartuples):
                        break

            # add leading query-only seqs
            if leading_insseq is not None:
                arr[arr_row_idx, initial_arr_col_idx] = (
                    f"({leading_insseq})" + arr[arr_row_idx, initial_arr_col_idx]
                )

    def truncate_array(arr, tmp_pileup_range, pileup_range):
        sl_start = tmp_pileup_range.index(pileup_range.start)
        if pileup_range.stop == tmp_pileup_range.stop:
            sl_end = None
        else:
            sl_end = tmp_pileup_range.index(pileup_range.stop)
        return arr[:, slice(sl_start, sl_end)]

    def make_pileup_from_array(
        arr,
        fetchedreads,
        uid_list,
        pileup_range,
        tmp_pileup_range,
        truncate,
        active_threshold,
        fasta,
        bam,
    ):
        df_columns = pileup_range if truncate else tmp_pileup_range
        df = pd.DataFrame(arr, columns=list(df_columns), index=uid_list)
        pileup = Pileup(
            df=df,
            chrom=fetchedreads.chrom,
            fasta=fasta,
            bam=bam,
            active_threshold=active_threshold,
        )
        pileup.fetchedreads_list.append(fetchedreads)
        pileup._set_active_info()
        return pileup

    # main
    # setup parameters
    sanity_check(pileup_coord, bam, fetchedreads)
    pileup_range, fetchedreads = set_params(pileup_coord, bam, fetchedreads)
    readlist, tmp_pileup_range, uid_list = make_readlist(fetchedreads, pileup_range)
    # create array
    arr = np.full(
        shape=(len(readlist), len(tmp_pileup_range)),
        fill_value=empty_value,
        dtype=object,
    )
    for arr_row_idx, read in enumerate(readlist):
        initial_arr_col_idx = read.reference_start - tmp_pileup_range.start
        try:
            add_seqs(read, arr, arr_row_idx, initial_arr_col_idx, del_value)
        except Exception as exc:
            try:
                cigar_sanitycheck(read)
            except Exception as exc_cigarpattern:
                raise exc_cigarpattern from exc
            else:
                raise exc
    # truncate array
    if truncate:
        arr = truncate_array(arr, tmp_pileup_range, pileup_range)
    # prepare results
    if as_array:
        if return_range:
            return (arr, pileup_range)
        else:
            return arr
    else:
        pileup = make_pileup_from_array(
            arr,
            fetchedreads,
            uid_list,
            pileup_range,
            tmp_pileup_range,
            truncate,
            active_threshold,
            fasta,
            bam,
        )
        if return_range:
            return (pileup, pileup_range)
        else:
            return pileup


# def get_pileup(bam, chrom, start0, end0, as_array=False, truncate=True, with_reads=False):
#    """Returns:
#        A numpy.ndarray object if "as_array" argument is True
#        A Pileup object if "as_array" argument is False
#    """
#
#    read_dict = {
#        readhandler.get_uid(x): x
#        for x in readhandler.get_fetch(bam, chrom, start0, end0, filter_fun=pileup_readfilter)
#    }
#    readlist = tuple(iter(read_dict.values()))
#    pileup, pileup_range = get_pileup_from_reads(readlist, as_array=as_array, return_range=True)
#
#    if truncate:
#        if as_array:
#            start0_idx = pileup_range.index(start0)
#            end0_idx = pileup_range.index(end0)
#            pileup = pileup[:, start0_idx:end0_idx]
#        else:
#            pileup.df = pileup.df.loc[:, list(range(start0, end0))]
#
#    if with_reads:
#        return pileup, read_dict
#    else:
#        return pileup


# def get_pileup_from_reads(readlist, as_array=False, pileup_range=None, return_range=False, del_char='*'):
#    def raise_err(read):
#        raise Exception(f"Unexpected cigar pattern:\n{read.to_string()}")
#
#    def cigar_sanitycheck(read):
#        """Assumes M.* or [IS]M.*"""
#        error = False
#        first_cigarop = read.cigartuples[0][0]
#        if first_cigarop != 0:
#            if first_cigarop not in (1, 4):
#                error = True
#            else:
#                if read.cigartuples[1][0] != 0:
#                    error = True
#        if error:
#            raise_err(read)
#
#    def get_pileup_range(readlist):
#        start0s = list()
#        end0s = list()
#        for read in readlist:
#            start0s.append(read.reference_start)
#            end0s.append(read.reference_end)
#        start0 = min(start0s)
#        end0 = max(end0s)
#
#        return range(start0, end0)
#
#    def check_overlaps_pileup_range(read, pileup_range):
#        if read.cigartuples[-1][0] in readhandler.CIGAROPS_QUERYONLY:
#            overlaps_upstream = (read.reference_end >= pileup_range.start)
#        else:
#            overlaps_upstream = (read.reference_end > pileup_range.start)
#
#        overlaps_downstream = (read.reference_start < pileup_range.stop)
#
#        return overlaps_upstream and overlaps_downstream
#
#    def handle_leading_insclip(read):
#        leading_insclip_len = 0
#        for cigartup_idx, (cigarop, cigarlen) in enumerate(read.cigartuples):
#            if cigarop in (1, 4):
#                leading_insclip_len += cigarlen
#            else:
#                break
#
#        if leading_insclip_len == 0:
#            leading_insseq = None
#        else:
#            leading_insseq = read.query_sequence[:leading_insclip_len]
#
#        query_idx = leading_insclip_len
#
#        return leading_insseq, cigartup_idx, query_idx
#
#    def handle_last_match(
#        read, arr, cigarlen, query_idx, arr_col_idx, arr_row_idx, cigartup_idx
#    ):
#        sl_query = slice(query_idx, query_idx + cigarlen)
#        seq_buffer = tuple(read.query_sequence[sl_query])
#        sl_arr_col = slice(arr_col_idx, arr_col_idx + cigarlen)
#        arr[arr_row_idx, sl_arr_col] = seq_buffer
#        cigartup_idx += 1
#
#        return cigartup_idx
#
#    def add_matches_before_last(
#        read, arr, cigarlen, query_idx, arr_col_idx, arr_row_idx
#    ):
#        sl_query = slice(query_idx, query_idx + cigarlen - 1)
#        seq = tuple(read.query_sequence[sl_query])
#        query_idx += cigarlen - 1
#
#        sl_arr_col = slice(arr_col_idx, arr_col_idx + cigarlen - 1)
#        arr[arr_row_idx, sl_arr_col] = seq
#        arr_col_idx += cigarlen - 1
#
#        return query_idx, arr_col_idx
#
#    def handle_del_after_match(
#        read, arr, query_idx, arr_col_idx, arr_row_idx, next_cigarlen
#    ):
#        seq = read.query_sequence[query_idx]
#        query_idx += 1
#        arr[arr_row_idx, arr_col_idx] = seq
#        arr_col_idx += 1
#
#        sl_arr_col = slice(arr_col_idx, arr_col_idx + next_cigarlen)
#        arr[arr_row_idx, sl_arr_col] = del_char
#        arr_col_idx += next_cigarlen
#
#        return query_idx, arr_col_idx
#
#    def handle_insclip_after_match(
#        read, arr, query_idx, arr_col_idx, arr_row_idx, next_cigarlen
#    ):
#        match_seq = read.query_sequence[query_idx]
#        query_idx += 1
#        insseq = read.query_sequence[query_idx : (query_idx + next_cigarlen)]
#        query_idx += next_cigarlen
#        seq = match_seq + f"({insseq})"
#
#        arr[arr_row_idx, arr_col_idx] = seq
#        arr_col_idx += 1
#
#        return query_idx, arr_col_idx
#
#    def add_seqs(read, arr, cigartup_idx, query_idx, arr_row_idx, arr_col_idx):
#        while True:
#            if cigartup_idx == len(read.cigartuples):
#                break
#
#            cigarop, cigarlen = read.cigartuples[cigartup_idx]
#            if cigarop != 0:
#                raise_err(read)
#
#            if cigartup_idx == len(read.cigartuples) - 1:
#                cigartup_idx = handle_last_match(
#                    read,
#                    arr,
#                    cigarlen,
#                    query_idx,
#                    arr_col_idx,
#                    arr_row_idx,
#                    cigartup_idx,
#                )
#                continue
#            else:
#                if cigarlen > 1:
#                    query_idx, arr_col_idx = add_matches_before_last(
#                        read, arr, cigarlen, query_idx, arr_col_idx, arr_row_idx
#                    )
#
#                next_cigarop, next_cigarlen = read.cigartuples[cigartup_idx + 1]
#                if next_cigarop == 2:  # deletion
#                    query_idx, arr_col_idx = handle_del_after_match(
#                        read, arr, query_idx, arr_col_idx, arr_row_idx, next_cigarlen
#                    )
#                elif next_cigarop in (1, 4):
#                    query_idx, arr_col_idx = handle_insclip_after_match(
#                        read, arr, query_idx, arr_col_idx, arr_row_idx, next_cigarlen
#                    )
#
#                cigartup_idx += 2
#                continue
#
#    def add_leading_insseq(read, leading_insseq, arr, arr_row_idx, arr_col_idx_init):
#        if leading_insseq is not None:
#            old_seq = arr[arr_row_idx, arr_col_idx_init]
#            new_seq = f"(l+{leading_insseq})" + old_seq
#            arr[arr_row_idx, arr_col_idx_init] = f"({leading_insseq})" + old_seq
#
#    # main
#    chrom = readlist[0].reference_name
#    if pileup_range is None:
#        pileup_range = get_pileup_range(readlist)
#
#    ncol = len(pileup_range)
#    nrow = len(readlist)
#    #arr = np.empty((nrow, ncol), dtype=object)
#    arr = np.full((nrow, ncol), None, dtype=object)
#
#    for arr_row_idx, read in enumerate(
#        read for read in readlist
#        if check_overlaps_pileup_range(read, pileup_range)
#    ):
#        try:
#            arr_col_idx_init = read.reference_start - pileup_range.start
#            arr_col_idx = arr_col_idx_init
#            leading_insseq, cigartup_idx, query_idx = handle_leading_insclip(read)
#            add_seqs(read, arr, cigartup_idx, query_idx, arr_row_idx, arr_col_idx)
#            add_leading_insseq(read, leading_insseq, arr, arr_row_idx, arr_col_idx_init)
#        except Exception as exc:
#            try:
#                cigar_sanitycheck(read)
#            except Exception as exc_cigarpattern:
#                raise exc_cigarpattern from exc
#            else:
#                raise exc
#
#    if as_array:
#        if return_range:
#            return (arr, pileup_range)
#        else:
#            return arr
#    else:
#        uids = [readhandler.get_uid(read) for read in readlist]
#        if len(uids) != len(set(uids)):
#            reads_as_string = '\n'.join(x.to_string() for x in readlist)
#            raise Exception(f'Duplicate read uid. Input reads are as following:\n{reads_as_string}')
#
#        df = pd.DataFrame(arr, columns=list(pileup_range), index=uids)
#        pileup = Pileup(df=df, uids=uids, chrom=chrom)
#        if return_range:
#            return (pileup, pileup_range)
#        else:
#            return pileup

import os
import tempfile
import subprocess
import gzip
import collections
import re

import pysam
import Bio.Align

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))


ALIGNER_FILLED = Bio.Align.PairwiseAligner(
    mode='global',
    match_score=2,
    mismatch_score=-3,
    query_internal_open_gap_score=-7,
    query_internal_extend_gap_score=-2,
    target_internal_open_gap_score=-7,
    target_internal_extend_gap_score=-2,
    query_left_open_gap_score=-7,
    query_left_extend_gap_score=-2,
)

ALIGNER_UNFILLED = Bio.Align.PairwiseAligner(
    mode='global',
    match_score=2,
    mismatch_score=-3,
    query_internal_open_gap_score=-7,
    query_internal_extend_gap_score=-2,
    target_internal_open_gap_score=-7,
    target_internal_extend_gap_score=-2,
)

ALIGNER_BLASTN = Bio.Align.PairwiseAligner(
    match_score=2,
    mismatch_score=-3,
    query_internal_open_gap_score=-7,
    query_internal_extend_gap_score=-2,
    target_internal_open_gap_score=-7,
    target_internal_extend_gap_score=-2,
)


def get_aligner_1():
    aligner = Bio.Align.PairwiseAligner(
        match_score=2,
        mismatch_score=-3,
        query_internal_open_gap_score=-7,
        query_internal_extend_gap_score=-2,
        target_internal_open_gap_score=-7,
        target_internal_extend_gap_score=-2,
    )

    return aligner


def get_aligner_2():
    aligner = Bio.Align.PairwiseAligner(
        match_score=2,
        mismatch_score=-3,
        query_internal_open_gap_score=-7,
        query_internal_extend_gap_score=-1,
        target_internal_open_gap_score=-7,
        target_internal_extend_gap_score=-1,
    )

    return aligner


def do_align(query, target, aligner):
    """Arguments:
    query, target: str object
    """
    return aligner.align(target, query)


def get_idx_tuples(alignment):
    idx_tuples = list()
    for idx in range(1, len(alignment.path)):
        idx_tuples.append( tuple(zip(alignment.path[idx-1], alignment.path[idx])) )

    return idx_tuples


def get_cigartuples(idx_tuples):
    cigartuples_list = list()

    for idx, tup_pair in enumerate(idx_tuples):
        target_len = tup_pair[0][1] - tup_pair[0][0]
        query_len = tup_pair[1][1] - tup_pair[1][0]

        if target_len == 0 and query_len != 0: # insertion
            cigartuples_list.append((1, query_len))
        elif target_len != 0 and query_len == 0: # deletion
            if idx != 0 and idx != len(idx_tuples)-1:
                cigartuples_list.append((2, target_len))
        else: # match
            cigartuples_list.append((0, target_len))

    return tuple(cigartuples_list)

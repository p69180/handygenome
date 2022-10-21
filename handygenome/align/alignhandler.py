import re
import collections
import itertools

import Bio
import numpy as np

import importlib
top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
libvcfspec = importlib.import_module(".".join([top_package_name, "variant", "vcfspec"]))


CIGAROPDICT_ITEMS = [
    ("M", 0),
    ("I", 1),
    ("D", 2),
    ("N", 3),
    ("S", 4),
    ("H", 5),
    ("P", 6),
    ("=", 7),
    ("X", 8),
    ("B", 9),
]
CIGAROPDICT = dict(CIGAROPDICT_ITEMS)
CIGAROPDICT_REV = dict((x[1], x[0]) for x in CIGAROPDICT_ITEMS)
CIGARPAT = re.compile(f'([0-9]+)([{"".join(CIGAROPDICT.keys())}])')
CIGAR_WALK_DICT = { # (target, query)
    0: (True, True),
    1: (False, True),
    2: (True, False),
    3: (True, False),
    4: (False, True),
    5: (False, False),
    6: (False, False),
    7: (True, True),
    8: (True, True),
}
CIGAROPS_BOTH = {0, 7, 8}
CIGAROPS_TARGETONLY = {2, 3}
CIGAROPS_QUERYONLY = {1, 4}


class Cigar:
    # constructors #
    @classmethod
    def from_cigarstring(cls, cigarstring):
        result = cls()
        result.cigartuples = cigarstring_to_cigartuples(cigarstring)
        return result
    
    @classmethod
    def from_cigartuples(cls, cigartuples):
        result = cls()
        result.cigartuples = cigartuples
        return result
    ################

    @property
    def cigarstring(self):
        return cigartuples_to_cigarstring(self.cigartuples)

    def iter_trailing_queryonly(self):
        return iter_trailing_queryonly(self.cigartuples)

    def iter_leading_queryonly(self):
        return iter_leading_queryonly(self.cigartuples)

    def check_DN_outside(self):
        # cigar I may be outside!
        opstrings_rmclip = [x for x in self.opstrings if x not in "SH"]

        return opstrings_rmclip[0] in "DN" or opstrings_rmclip[-1] in "DN"

    def check_clip_inside(self):
        if len(self.opstrings) <= 2:
            return False
        else:
            return not {"S", "H"}.isdisjoint(self.opstrings[1:-1])

    def check_unexpected_pattern(self):
        return (
            self.check_DN_outside()
            or self.check_clip_inside()
            or (not {"B", "P"}.isdisjoint(self.opstrings))
        )


#def get_cigartuples(cigarstring):
#    return [
#        (CIGAROPDICT[cigarop], int(count))
#        for (count, cigarop) in CIGARPAT.findall(cigarstring)
#    ]


def cigarstring_to_cigartuples(cigarstring):
    result = list()
    for (count, opstring) in CIGARPAT.findall(cigarstring):
        opcode = CIGAROPDICT[opstring]
        count = int(count)
        result.append((opcode, count))
    return result


def cigartuples_to_cigarstring(cigartuples):
    buffer = list()
    for opcode, count in cigartuples:
        opstring = CIGAROPDICT_REV[opcode]
        buffer.append(str(count))
        buffer.append(opstring)
    return "".join(buffer)


def cigartuples_to_cigarstring_spread(cigartuples):
    buffer = list()
    for opcode, count in cigartuples:
        opstring = CIGAROPDICT_REV[opcode]
        buffer.extend(itertools.repeat(opstring, count))
    return "".join(buffer)


def iter_trailing_queryonly(cigartuples):
    for opcode, count in reversed(cigartuples):
        if opcode in CIGAROPS_QUERYONLY:
            yield opcode, count
        else:
            break


def iter_leading_queryonly(cigartuples):
    for opcode, count in cigartuples:
        if opcode in CIGAROPS_QUERYONLY:
            yield opcode, count
        else:
            break


def get_query_length(cigartuples):
    result = 0
    for opcode, count in cigartuples:
        walk_target, walk_query = CIGAR_WALK_DICT[opcode]
        if walk_query:
            result += count
    return result


def get_target_length(cigartuples):
    result = 0
    for opcode, count in cigartuples:
        walk_target, walk_query = CIGAR_WALK_DICT[opcode]
        if walk_target:
            result += count
    return result


def split_cigar(cigartuples, reference_start0, split_range0):
    tuples_before = list()
    tuples_within = list()
    tuples_after = list()

    for cigarunit_tup, target_range0 in zip(
        cigartuples, 
        walk_cigar_targetonly(cigartuples, reference_start0)
    ):
        if len(target_range0) == 0:
            # query-only cigar operation
            # start and stop of target_range0 are the same
            # start of target_range0 indicates the position next to the insertion site
            pos0 = target_range0.start
            if pos0 < split_range0.start:
                tuples_before.append(cigarunit_tup)
            elif pos0 > split_range0.stop:
                tuples_after.append(cigarunit_tup)
            else:
                # query-only cigarops immediately before or after "split_range0" are all included in "tuples_within"
                tuples_within.append(cigarunit_tup)
        else:
            # current cigarwalk entirely before split_range0
            if target_range0.stop <= split_range0.start:
                tuples_before.append(cigarunit_tup)
            # current cigarwalk entirely after split_range0
            elif target_range0.start >= split_range0.stop:
                tuples_after.append(cigarunit_tup)
            # current cigarwalk overlaps split_range0
            else:
                # add to tuples_before
                if target_range0.start < split_range0.start:
                    tuple_to_add = (
                        cigarunit_tup[0], 
                        split_range0.start - target_range0.start,
                    )
                    tuples_before.append(tuple_to_add)
                # add to tuples_after
                if target_range0.stop > split_range0.stop:
                    tuple_to_add = (
                        cigarunit_tup[0], 
                        target_range0.stop - split_range0.stop,
                    )
                    tuples_after.append(tuple_to_add)
                # add to tuples_within
                tuple_to_add = (
                    cigarunit_tup[0], 
                    (
                        min(target_range0.stop, split_range0.stop)
                        - max(target_range0.start, split_range0.start)
                    )
                )
                tuples_within.append(tuple_to_add)

    return tuples_before, tuples_within, tuples_after


############################################
# Bio.Align.PairwiseAlignment-related ones #
############################################

def walk_cigar(cigartuples, target_start0):
    query_start0 = 0

    target_end0 = target_start0
    query_end0 = query_start0

    for cigartup in cigartuples:
        opcode, count = cigartup
        walk_target, walk_query = CIGAR_WALK_DICT[opcode]
        if walk_target:
            target_end0 += count
        if walk_query:
            query_end0 += count

        target_range0 = range(target_start0, target_end0)
        query_range0 = range(query_start0, query_end0)

        yield target_range0, query_range0

        target_start0 = target_end0
        query_start0 = query_end0


def walk_cigar_targetonly(cigartuples, target_start0):
    #query_start0 = 0

    target_end0 = target_start0
    #query_end0 = query_start0

    for cigartup in cigartuples:
        opcode, count = cigartup
        walk_target, walk_query = CIGAR_WALK_DICT[opcode]
        if walk_target:
            target_end0 += count
        #if walk_query:
        #    query_end0 += count

        target_range0 = range(target_start0, target_end0)
        #query_range0 = range(query_start0, query_end0)

        yield target_range0

        target_start0 = target_end0
        #query_start0 = query_end0


def alignment_to_cigartuples(
    alignment, 
    match_as_78=False, del_as_skip=False, 
    left_ins_as_clip=True, right_ins_as_clip=True,
    remove_left_del=True, remove_right_del=True,
):
    """Args:
        alignment: Bio.Align.PairwiseAlignment object
        assumes global alignment
    """
    # set params
    targetonly_cigarop = 3 if del_as_skip else 2
    match_handler = _match_handler_78 if match_as_78 else _match_handler_0
    # init cigartuples
    cigartuples = list()
    walks = get_walks(alignment, copy=False)
    #list(path_to_walks(alignment.path))
    for target_range0, query_range0 in walks:
        if len(target_range0) == 0:
            cigartuples.append((1, len(query_range0)))
        elif len(query_range0) == 0:
            cigartuples.append((targetonly_cigarop, len(target_range0)))
        else:
            match_handler(alignment, target_range0, query_range0, cigartuples)
    # postprocess - I at borders into clip
    if left_ins_as_clip:
        if cigartuples[0][0] == 1:
            cigartuples[0] = (4, cigartuples[0][1])
    if right_ins_as_clip:
        if cigartuples[-1][0] == 1:
            cigartuples[-1] = (4, cigartuples[-1][1])
    # postprocess - remove flanking targetonly
        # assumes global alignment
    if remove_left_del:
        if cigartuples[0][0] == targetonly_cigarop:
            del cigartuples[0]
    if remove_right_del:
        if cigartuples[-1][0] == targetonly_cigarop:
            del cigartuples[-1]
    # get target offset
    if remove_left_del:
        for target_range0, query_range0 in walks:
            if len(query_range0) > 0:
                target_offset = target_range0.start 
                break
    else:
        target_offset = 0

    return cigartuples, target_offset


def _match_handler_78(alignment, target_range0, query_range0, cigartuples):
    cigarop_list = list()
    for target_idx, query_idx in zip(target_range0, query_range0):
        if alignment.target[target_idx] == alignment.query[query_idx]:
            cigarop_list.append(7)
        else:
            cigarop_list.append(8)

    cigartuples.extend(
        [(x[0], len(tuple(x[1]))) for x in itertools.groupby(cigarop_list)]
    )


def _match_handler_0(alignment, target_range0, query_range0, cigartuples):
    cigartuples.append((0, len(target_range0)))


def walks_to_path(cigarwalks):
    """'path' means 'path' attribute of Bio.Align.PairwiseAlignment object"""
    result = list()
    iterator = iter(cigarwalks)

    target_range0, query_range0 = next(iterator)
    result.append((target_range0.start, query_range0.start))
    result.append((target_range0.stop, query_range0.stop))

    for target_range0, query_range0 in iterator:
        result.append((target_range0.stop, query_range0.stop))

    return tuple(result)


def path_to_walks(alignpath):
    """Args:
        alignpath: 'path' attribute of Bio.Align.PairwiseAlignment object
    """
    for x0, x1 in common.pairwise(alignpath):
        yield range(x0[0], x1[0]), range(x0[1], x1[1])  # target_range, query_range


def set_walks(alignment):
    alignment.walks = list(path_to_walks(alignment.path))


def get_walks(alignment, copy=False):
    if not hasattr(alignment, 'walks'):
        set_walks(alignment)

    if copy:
        return alignment.walks.copy()
    else:
        return alignment.walks


def amend_outer_insdel_left(alignment):
    walks = get_walks(alignment, copy=False)
    if len(walks) < 3:
        return alignment

    if (
        len(walks[0][0]) == 0 and 
        len(walks[1][1]) == 0 and
        (len(walks[2][0]) > 0 and len(walks[2][1]) > 0)
    ):
        new_walks = walks.copy()
        new_walks[0] = (walks[1][0], range(0, 0))
        new_walks[1] = (range(walks[1][0].stop, walks[1][0].stop), walks[0][1])
        return Bio.Align.PairwiseAlignment(
            alignment.target, 
            alignment.query, 
            walks_to_path(new_walks),
            0,
        )
    else:
        return alignment


def amend_outer_insdel_right(alignment):
    walks = get_walks(alignment, copy=False)
    if len(walks) < 3:
        return alignment

    if (
        len(walks[-1][0]) == 0 and 
        len(walks[-2][1]) == 0 and
        (len(walks[-3][0]) > 0 and len(walks[-3][1]) > 0)
    ):
        new_walks = walks.copy()
        new_walks[-1] = (walks[-2][0], range(0, 0))
        new_walks[-2] = (range(walks[-2][0].start, walks[-2][0].start), walks[-1][1])
        return Bio.Align.PairwiseAlignment(
            alignment.target, 
            alignment.query, 
            walks_to_path(new_walks),
            0,
        )
    else:
        return alignment

def amend_outer_insdel_both(alignment):
    return amend_outer_insdel_right(
        amend_outer_insdel_left(alignment)
    )


#def amend_outer_insdel_right(alignment):
#    rev_alignment = reverse_alignment(alignment)
#    rev_alignment_amended = amend_outer_insdel_left(rev_alignment)
#    return reverse_alignment(rev_alignment_amended)


def show_alignment_with_numbers(alignment):
    aln_length = sum(max(len(x), len(y)) for (x, y) in get_walks(alignment))
    string_tens = ''.join(
        str(x) + ' ' * (9 if x == 0 else 9 - np.floor(np.log10(x)).astype('int'))
        for x in range(aln_length // 10 + 1)
    )
    string_ones = ''.join(str(x)[-1] for x in range(aln_length))
    print(string_tens)
    print(string_ones)
    print(alignment, end='')
    print(string_ones)
    print(string_tens)


def reverse_alignment(alignment):
    """Args:
        alignment: Bio.Align.PairwiseAlignment object
    """
    return Bio.Align.PairwiseAlignment(
        alignment.target[::-1], 
        alignment.query[::-1], 
        [
            (len(alignment.target) - x[0], len(alignment.query) - x[1]) 
            for x in alignment.path[::-1]
        ], 
        0,
    )


#def remove_flanking_query_gaps(alignment):
#    offset = 0
#    # strip leading query gap
#    walks = list(path_to_walks(alignment.path))
#    for idx, (target_walk, query_walk) in enumerate(walks):
#        if len(query_walk) != 0:
#            break
#        else:
#            offset += len(target_walk)
#    walks = walks[idx:]
#    # correct coordinates of target_walk's
#    for idx in range(len(walks)):
#        old_target_walk, old_query_walk = walks[idx]
#        new_target_walk = range(old_target_walk.start - offset, old_target_walk.stop - offset)
#        walks[idx] = (new_target_walk, old_query_walk)
#    # strip trailing query gap
#    for idx, (target_walk, query_walk) in enumerate(walks[::-1]):
#        if len(query_walk) != 0:
#            break
#    if idx > 0:
#        walks = walks[:-idx]
#    # make new alignment
#    new_aln = Bio.Align.PairwiseAlignment(
#        alignment.target[offset:], 
#        alignment.query,
#        walks_to_path(walks),
#        0,
#    )
#        
#    return new_aln, offset


def remove_leading_query_gaps_component_input(target, query, walks, as_alignment=True):
    offset = 0
    for idx, (target_walk, query_walk) in enumerate(walks):
        if len(query_walk) == 0:
            offset += len(target_walk)
        else:
            break
    new_walks = walks[idx:]
    new_target = target[offset:]
    # correct coordinates of target_walk's
    for idx in range(len(new_walks)):
        old_target_walk, old_query_walk = new_walks[idx]
        new_target_walk = range(old_target_walk.start - offset, old_target_walk.stop - offset)
        new_walks[idx] = (new_target_walk, old_query_walk)
    # return
    if as_alignment:
        new_aln = Bio.Align.PairwiseAlignment(new_target, query, walks_to_path(new_walks), 0)
        return new_aln, offset
    else:
        return new_target, new_walks, offset


def remove_leading_query_gaps(alignment, as_alignment=True):
    return remove_leading_query_gaps_component_input(alignment.target, alignment.query, get_walks(alignment, copy=False), as_alignment=as_alignment)


def remove_trailing_query_gaps_component_input(target, query, walks, as_alignment=True):
    offset = 0
    for idx, (target_walk, query_walk) in enumerate(walks[::-1]):
        if len(query_walk) == 0:
            offset += len(target_walk)
        else:
            break

    if idx > 0:
        new_walks = walks[:-idx]
    else:
        new_walks = walks.copy()

    if offset == 0:
        new_target = target
    else:
        new_target = target[:-offset]

    # return
    if as_alignment:
        new_aln = Bio.Align.PairwiseAlignment(new_target, query, walks_to_path(new_walks), 0)
        return new_aln, offset
    else:
        return new_target, new_walks, offset


def remove_trailing_query_gaps(alignment, as_alignment=True):
    return remove_trailing_query_gaps_component_input(alignment.target, alignment.query, get_walks(alignment, copy=False), as_alignment=as_alignment)


def remove_flanking_query_gaps(alignment):
    target_lstrip, walks_lstrip, leading_offset = remove_leading_query_gaps(alignment, as_alignment=False)
    target_lrstrip, walks_lrstrip, trailing_offset = remove_trailing_query_gaps_component_input(
        target_lstrip, alignment.query, walks_lstrip, as_alignment=False,
    )
    new_aln = Bio.Align.PairwiseAlignment(
        target_lrstrip, 
        alignment.query, 
        walks_to_path(walks_lrstrip),
        0,
    )
        
    return new_aln, leading_offset


def alignment_to_vcfspec(alignment, target_start0, chrom, fasta, lstrip_query_gaps=True, rstrip_query_gaps=True):
    def groupkey_walks(x):
        if len(x[0]) == 0 or len(x[1]) == 0:
            return 'indel'
        else:
            return 'match'

    def groupkey_matches(x):
        return x[0] != x[1]

    def handle_single_ins(current_pos0, target_idx, query_idx, query_walk, fasta, chrom, new_aln, vcfspec_list):
        pos = current_pos0
        if target_idx == 0:
            ref = fasta.fetch(chrom, current_pos0 - 1, current_pos0)
        else:
            ref = new_aln.target[target_idx - 1]
        inserted_seq = new_aln.query[query_idx:(query_idx + len(query_walk))]
        alt = ref + inserted_seq
        vcfspec_list.append(libvcfspec.Vcfspec(chrom, pos, ref, (alt,), fasta=fasta))

    def handle_single_del(current_pos0, target_idx, target_walk, fasta, chrom, new_aln, vcfspec_list):
        pos = current_pos0
        if target_idx == 0:
            preceding_base = fasta.fetch(chrom, current_pos0 - 1, current_pos0)
        else:
            preceding_base = new_aln.target[target_idx - 1]
        ref = preceding_base + new_aln.target[target_idx:(target_idx + len(target_walk))]
        alt = ref[0]
        vcfspec_list.append(libvcfspec.Vcfspec(chrom, pos, ref, (alt,), fasta=fasta))

    def handle_consecutive_indels(new_aln, walks_subset, target_idx, query_idx, current_pos0, chrom, fasta, vcfspec_list):
        target_walk_length = sum(len(target_walk) for target_walk, query_walk in walks_subset)
        query_walk_length = sum(len(query_walk) for target_walk, query_walk in walks_subset)
        ref = new_aln.target[target_idx:(target_idx + target_walk_length)]
        alt = new_aln.query[query_idx:(query_idx + query_walk_length)]
        pos = current_pos0 + 1
        vcfspec = libvcfspec.Vcfspec(chrom, pos, ref, (alt,), fasta=fasta)
        #vcfspec = vcfspec.parsimonious()
        vcfspec_list.append(vcfspec)

    def handle_match(new_aln, target_idx, query_idx, target_walk, current_pos0, chrom, fasta, vcfspec_list):
        target_seq = new_aln.target[target_idx:(target_idx + len(target_walk))]
        query_seq = new_aln.query[query_idx:(query_idx + len(target_walk))]
        idx_offset = 0
        for is_mismatch, seq_pairs in itertools.groupby(zip(target_seq, query_seq), key=groupkey_matches):
            seq_pairs = tuple(seq_pairs)
            if is_mismatch:
                pos = current_pos0 + idx_offset + 1
                ref = ''.join(x[0] for x in seq_pairs)
                alt = ''.join(x[1] for x in seq_pairs)
                vcfspec_list.append(libvcfspec.Vcfspec(chrom, pos, ref, (alt,), fasta=fasta))
            idx_offset += len(seq_pairs)

    # set paramters
    if lstrip_query_gaps:
        if rstrip_query_gaps:
            new_aln, leading_offset = remove_flanking_query_gaps(alignment)
        else:
            new_aln, leading_offset = remove_leading_query_gaps(alignment, as_alignment=True)
    else:
        if rstrip_query_gaps:
            new_aln, trailing_offset = remove_trailing_query_gaps(alignment, as_alignment=True)
            leading_offset = 0
        else:
            new_aln = alignment
            leading_offset = 0
    
    target_start0 += leading_offset
    current_pos0 = target_start0
    target_idx = 0
    query_idx = 0
    vcfspec_list = list()
    
    # main
    for walktype, walks_subset in itertools.groupby(get_walks(new_aln, copy=False), key=groupkey_walks):
        walks_subset = tuple(walks_subset)

        if walktype == 'match':
            if len(walks_subset) != 1:
                raise Exception(
                    f'Longer than one consecutive runs of match walks.\n'
                    f'Input alignment:\n'
                    f'{alignment}'
                )
            target_walk, query_walk = walks_subset[0]
            handle_match(new_aln, target_idx, query_idx, target_walk, current_pos0, chrom, fasta, vcfspec_list)
        elif walktype == 'indel':
            if len(walks_subset) == 1:
                target_walk, query_walk = walks_subset[0]
                if len(target_walk) == 0:  # ins
                    handle_single_ins(current_pos0, target_idx, query_idx, query_walk, fasta, chrom, new_aln, vcfspec_list)
                elif len(query_walk) == 0:  # del
                    handle_single_del(current_pos0, target_idx, target_walk, fasta, chrom, new_aln, vcfspec_list)
            else:
                if not (
                    any(len(target_walk) > 0 for target_walk, query_walk in walks_subset) and 
                    any(len(query_walk) > 0 for target_walk, query_walk in walks_subset)
                ):
                    raise Exception(
                        f'Consecutive run of indels are composed only of insertion or deletion.\n'
                        f'Input alignment:\n'
                        f'{alignment}'
                    )
                handle_consecutive_indels(new_aln, walks_subset, target_idx, query_idx, current_pos0, chrom, fasta, vcfspec_list)

        for target_walk, query_walk in walks_subset:
            current_pos0 += len(target_walk)
            target_idx += len(target_walk)
            query_idx += len(query_walk)
        
    return tuple(vcfspec_list)


def alignment_hasher(aln):
    return (aln.target, aln.query, tuple(aln.path))


def remove_identical_alignments(alignments):
    hash_store = set()
    for x in alignments:
        hash_val = alignment_hasher(x)
        if hash_val in hash_store:
            continue
        else:
            hash_store.add(hash_val)
            yield x


############################
# alignment post-selection #
############################

class AlignmentTieError(Exception):
    pass


def tiebreaker_scorer_leftmost_gaps(alignment):
    score = 0
    for target_range0, query_range0 in get_walks(alignment, copy=False):
        if len(target_range0) == 0 or len(query_range0) == 0:
            score -= target_range0.start
    return score


def tiebreaker_scorer_gap_length_sum(alignment):
    """Those with shorter gap lengths sum are favored"""
    score = 0
    walks = get_walks(alignment, copy=False)
    for idx, (target_range0, query_range0) in enumerate(walks):
        if len(target_range0) == 0:
            score -= len(query_range0)
        elif len(query_range0) == 0:
            if idx != 0 and idx != (len(walks) - 1):
                score -= len(target_range0)
    return score


def tiebreaker_scorer_gap_ordering(alignment):
    """Favors longer gaps coming earlier than shorter gaps"""
    score = 0
    walks = get_walks(alignment, copy=False)
    for idx, (target_range0, query_range0) in enumerate(walks):
        if len(target_range0) == 0:
            score -= len(query_range0) * target_range0.start
        elif len(query_range0) == 0:
            if idx != 0 and idx != (len(walks) - 1):
                score -= len(target_range0) * query_range0.start
    return score


def tiebreakers_merged_main(alignments):
    # alignment with the highest score is selected
    selected_alns = common.multi_max(alignments, key=tiebreaker_scorer_gap_length_sum, with_target_val=False)
    if len(selected_alns) != 1:
        selected_alns = common.multi_max(selected_alns, key=tiebreaker_scorer_leftmost_gaps, with_target_val=False)
        if len(selected_alns) != 1:
            selected_alns = common.multi_max(selected_alns, key=tiebreaker_scorer_gap_ordering, with_target_val=False)

    return selected_alns


def alignment_tiebreaker(alignments, raise_with_failure=True):
    selected_alns = tiebreakers_merged_main(alignments)
    if len(selected_alns) != 1 and raise_with_failure:
        alignments_string = list()
        for x in selected_alns:
            alignments_string.append(
                f'target: {x.target}\nquery: {x.query}\npath: {x.path}\n{str(x)}'
            )
        alignments_string = '\n'.join(alignments_string)
        raise AlignmentTieError(
            f'Failed to break alignment tie. Alignments are:\n{alignments_string}'
        )
        
    return selected_alns[0]



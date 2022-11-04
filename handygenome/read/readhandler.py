import sys
import re
import itertools
import collections

import pysam
import numpy as np
import pandas as pd

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))


CIGAROPDICT = {
    "M": 0,
    "I": 1,
    "D": 2,
    "N": 3,
    "S": 4,
    "H": 5,
    "P": 6,
    "=": 7,
    "X": 8,
    "B": 9,
}
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
CIGAROPS_TARGETONLY = {2, 3}
CIGAROPS_QUERYONLY = {1, 4}


class NoMDTagError(Exception):
    pass


###################################

ReadUID = collections.namedtuple(
    'ReadUID',
    ('qname', 'flag', 'chrom', 'pos0'),
)


def get_uid(read):
    return ReadUID(read.query_name, read.flag, read.reference_name, read.reference_start)
    #return (read.query_name, read.flag, read.reference_name, read.reference_start)


def get_read_dict(read_iter):
    return {get_uid(x): x for x in read_iter}


def readfilter_bad_read(read):
    """returns True if the read is to be included"""

    return not (
        (not read.is_paired)
        or read.is_unmapped
        or read.mate_is_unmapped
        or read.is_duplicate
    )


def readfilter_pileup(read):
    return (
        readfilter_bad_read(read) and 
        (not read.is_supplementary) and 
        (not read.is_secondary)
    )


def get_fetch(bam, chrom, start, end, readfilter=None):
    """
    - returns a generator
    - similar to bam.fetch iterator except read filtering
    """
    # set default filtering function
    if readfilter is None:
        readfilter = readfilter_bad_read
    # main
    for read in bam.fetch(chrom, start, end):
        if readfilter(read):
            yield read


# cigar-related classes and functions

#class CigarUnit:
#    def __init__(self, opcode, opstring, count):
#        self.opcode = opcode
#        self.opstring = opstring
#        self.count = count

class Cigar:
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
    CIGAROPS_TARGETONLY = {2, 3}
    CIGAROPS_QUERYONLY = {1, 4}

    # constructors #
    @classmethod
    def from_cigarstring(cls, cigarstring):
        result = cls()
        result.cigartuples = cls.cigarstring_to_cigartuples(cigarstring)
        return result
    
    @classmethod
    def from_cigartuples(cls, cigartuples):
        result = cls()
        result.cigartuples = cigartuples
        return result
    ################

    # converters
    @classmethod
    def cigarstring_to_cigartuples(cls, cigarstring):
        result = list()
        for (count, opstring) in cls.CIGARPAT.findall(cigarstring):
            opcode = cls.CIGAROPDICT[opstring]
            count = int(count)
            result.append((opcode, count))
        return result

    @classmethod
    def cigartuples_to_cigarstring(cls, cigartuples):
        buffer = list()
        for opcode, count in cigartuples:
            opstring = cls.CIGAROPDICT_REV[opcode]
            buffer.append(str(count))
            buffer.append(opstring)
        return "".join(buffer)

    # others
    @property
    def cigarstring(self):
        return self.cigartuples_to_cigarstring(self.cigartuples)

    def iter_trailing_queryonly(self):
        for opcode, count in reversed(self.cigartuples):
            if opcode in CIGAROPS_QUERYONLY:
                yield opcode, count
            else:
                break

    def iter_leading_queryonly(self):
        for opcode, count in self.cigartuples:
            if opcode in CIGAROPS_QUERYONLY:
                yield opcode, count
            else:
                break

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


CigarWalk = collections.namedtuple(
    'CigarWalk', 
    ('cigartuple', 'target_range0', 'query_range0'),
)
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

        target_start0 = target_end0
        query_start0 = query_end0

        yield CigarWalk(cigartup, target_range0, query_range0)


def split_cigar(cigartuples, reference_start0, split_range0):
    def add_buffer(tuple_list, queryonly_buffer):
        while queryonly_buffer:
            tuple_list.append(queryonly_buffer.pop(0).cigartuple)

    def add_current(tuple_list, cigarwalk):
        tuple_list.append(cigarwalk.cigartuple)

    tuples_before = list()
    tuples_within = list()
    tuples_after = list()
    reference_pointer = reference_start0
    queryonly_buffer = list()

    for cigarwalk in walk_cigar(cigartuples, reference_start0):
        reference_pointer += len(cigarwalk.target_range0)

        if len(cigarwalk.target_range0) == 0:
            queryonly_buffer.append(cigarwalk)
        else:
            # current cigarwalk entirely before split_range0
            if cigarwalk.target_range0.stop <= split_range0.start:
                add_buffer(tuples_before, queryonly_buffer)
                add_current(tuples_before, cigarwalk)
            # current cigarwalk entirely after split_range0
            elif cigarwalk.target_range0.start >= split_range0.stop:
                add_buffer(tuples_after, queryonly_buffer)
                add_current(tuples_after, cigarwalk)
            # current cigarwalk overlaps split_range0
            else:
                # handle buffer
                if cigarwalk.target_range0.start < split_range0.start:
                    add_buffer(tuples_before, queryonly_buffer)
                else:
                    add_buffer(tuples_within, queryonly_buffer)
                # add to tuples_before
                if cigarwalk.target_range0.start < split_range0.start:
                    tuple_to_add = (
                        cigarwalk.cigartuple[0], 
                        split_range0.start - cigarwalk.target_range0.start,
                    )
                    tuples_before.append(tuple_to_add)
                # add to tuples_after
                if cigarwalk.target_range0.stop > split_range0.stop:
                    tuple_to_add = (
                        cigarwalk.cigartuple[0], 
                        cigarwalk.target_range0.stop - split_range0.stop,
                    )
                    tuples_after.append(tuple_to_add)
                # add to tuples_within
                tuple_to_add = (
                    cigarwalk.cigartuple[0], 
                    (
                        min(cigarwalk.target_range0.stop, split_range0.stop)
                        - max(cigarwalk.target_range0.start, split_range0.start)
                    )
                )
                tuples_within.append(tuple_to_add)
                    
    # now reference_pointer is equal to read.reference_end
    if queryonly_buffer:
        if reference_pointer < split_range0.start:
            add_buffer(tuples_before, queryonly_buffer)
        else:
            if reference_pointer < split_range0.stop:
                add_buffer(tuples_within, queryonly_buffer)
            else:
                add_buffer(tuples_after, queryonly_buffer)

    return tuples_before, tuples_within, tuples_after


def get_aligned_pairs(pos, cigartuples):
    aligned_pairs = list()
    querypos0 = 0
    refpos0 = pos - 1

    for (cigarop, count) in cigartuples:
        if cigarop in (2, 3):  # D, N
            for _ in range(count):
                aligned_pairs.append((None, refpos0))
                refpos0 += 1
        elif cigarop in (1, 4, 5):  # I, S, H
            for _ in range(count):
                aligned_pairs.append((querypos0, None))
                querypos0 += 1
        elif cigarop in (0, 7, 8):  # M, =, X
            for _ in range(count):
                aligned_pairs.append((querypos0, refpos0))
                querypos0 += 1
                refpos0 += 1

    return aligned_pairs


def get_pairs_dict(read, fasta=None, skip_refseq=False):
    """Returns:
    Dict {
        'querypos0': A tuple of query positions
        'refpos0': A tuple of reference genome positions
        'refseq': A tuple of reference sequences (lowercase with mismatch)
        'cigarop': A tuple of cigar operations, as integers
        }
    """

    if read.has_tag("MD"):
        raw = read.get_aligned_pairs(with_seq=True)
        zipped = zip(*raw)
        pairs_dict = dict()
        pairs_dict["querypos0"] = next(zipped)
        pairs_dict["refpos0"] = next(zipped)
        if not skip_refseq:
            pairs_dict["refseq"] = next(zipped)
    else:
        raw = read.get_aligned_pairs(with_seq=False)
        zipped = zip(*raw)
        pairs_dict = dict()
        pairs_dict["querypos0"] = next(zipped)
        pairs_dict["refpos0"] = next(zipped)
        if not skip_refseq:
            pairs_dict["refseq"] = get_pairs_dict_refseq(read, fasta, raw)

    # cigarop
    pairs_dict["cigarop"] = tuple(
        itertools.chain.from_iterable(
            itertools.repeat(cigarop, count)
            for (cigarop, count) in read.cigartuples
            if cigarop != 5
        )
    )

    if len(pairs_dict["cigarop"]) != len(pairs_dict["refpos0"]):
        raise Exception(
            f"The length of cigar operations differ from that of "
            f"aligned pairs:\n"
            f"{read.to_string()}"
        )

    return pairs_dict


def get_pairs_dict_refseq(read, fasta, raw_aligned_pairs):
    refseq = list()

    non_None_refpos0 = [
        refpos0 for (querypos0, refpos0) in raw_aligned_pairs if refpos0 is not None
    ]
    if len(non_None_refpos0) == 0:
        fetched_seq = ""
    else:
        fetched_seq = fasta.fetch(
            read.reference_name, non_None_refpos0[0], non_None_refpos0[-1] + 1
        )
    it = iter(fetched_seq)
    for querypos0, refpos0 in raw_aligned_pairs:
        if refpos0 is None:
            refseq.append(None)
        else:
            if querypos0 is None:
                read_base = None
            else:
                read_base = read.query_sequence[querypos0].upper()

            ref_base = next(it).upper()
            if read_base == ref_base:
                refseq.append(ref_base)
            else:
                refseq.append(ref_base.lower())

    return refseq


# def get_MM(read, pairs_dict):
#    if pairs_dict is None:
#        MM = None
#    else:
#        cigarM_indices = [
#            True
#            if (tup[0] is not None and tup[1] is not None) else
#            False
#            for tup in zip(pairs_dict['querypos0'], pairs_dict['refpos0'])]
#        pairs_dict_seq_subset = itertools.compress(
#            pairs_dict['refseq'], cigarM_indices)
#        MM = len([x for x in pairs_dict_seq_subset if x.islower()])
#
#    return MM


#####################################################


def check_cigar_clip_inside(cigartuples):
    if len(cigartuples) <= 2:
        return False
    else:
        cigarops_inside = [x[0] for x in cigartuples[1:-1]]
        return not {4, 5}.isdisjoint(cigarops_inside)


def check_cigar_DN_outside(cigarops):
    # cigar I can be outside!
    cigarops_rmclip = [x for x in cigarops if x not in (4, 5)]
    return cigarops_rmclip[0] in (2, 3) or cigarops_rmclip[-1] in (2, 3)


def check_unexpected_cigar_pattern(cigartuples):
    cigarops = [x[0] for x in cigartuples]
    if check_cigar_DN_outside(cigarops):
        return True
    elif check_cigar_clip_inside(cigartuples):
        return True
    elif "B" in cigarops:
        return True
    elif "P" in cigarops:
        return True
    else:
        return False


#############################################################


def get_padded_seqs(read, fasta):
    if check_unexpected_cigar_pattern(read.cigartuples):
        raise Exception(
            f"Unexpected cigar pattern with the input read: " f"{read.to_string()}"
        )

    ref_seq_padded = list()
    read_seq_padded = list()

    ref_seq = list(
        fasta.fetch(read.reference_name, read.reference_start, read.reference_end)
    )
    read_seq = list(read.query_sequence)

    for idx_cigar, tup in enumerate(read.cigartuples):
        op, l = tup  # op: operation ; l: length
        if op in (0, 7, 8):  # cigar M, =, X
            for i in range(l):
                ref_seq_padded.append(ref_seq.pop(0))
                read_seq_padded.append(read_seq.pop(0))
        elif op == 1:  # cigar I
            ref_seq_padded.append(None)
            ins_seq_buffer = list()
            for i in range(l):
                ins_seq_buffer.append(read_seq.pop(0))
            read_seq_padded.append("".join(ins_seq_buffer))
        elif op == 2:  # cigar D
            read_seq_padded.append(None)
            del_seq_buffer = list()
            for i in range(l):
                del_seq_buffer.append(ref_seq.pop(0))
            ref_seq_padded.append("".join(del_seq_buffer))
        elif op == 3:  # cigar N
            del ref_seq[:l]
        elif op == 4:  # cigar S
            del read_seq[:l]
        elif op == 5:  # cigar H
            pass

    return ref_seq_padded, read_seq_padded


def get_MD(read, fasta=None, ref_seq_padded=None, read_seq_padded=None):
    if (ref_seq_padded is None) or (read_seq_padded is None):
        ref_seq_padded, read_seq_padded = get_padded_seqs(read, fasta)

    MD_list = list()
    identical = 0
    for ref_base, read_base in zip(ref_seq_padded, read_seq_padded):
        if ref_base is not None and read_base is not None:
            if read_base == "=":  # means a base identical to reference base
                identical += 1
            else:
                if ref_base == read_base:  # identical match
                    identical += 1
                else:  # different match
                    MD_list.append(str(identical))
                    identical = 0
                    MD_list.append(ref_base)
        elif ref_base is None:  # cigar I
            pass
        elif read_base is None:  # cigar D
            MD_list.append(str(identical))
            identical = 0
            MD_list.append("^" + ref_base)

    MD_list.append(str(identical))

    return "".join(MD_list)


def get_NM(read, fasta=None, ref_seq_padded=None, read_seq_padded=None):
    if (ref_seq_padded is None) or (read_seq_padded is None):
        ref_seq_padded, read_seq_padded = get_padded_seqs(read, fasta)

    NM = 0
    for ref_base, read_base in zip(ref_seq_padded, read_seq_padded):
        if ref_base is None:  # cigar I
            NM += len(read_base)
        elif read_base is None:  # cigar D
            NM += len(ref_base)
        else:
            if read_base == "=":  # means a base identical to reference base
                pass
            else:
                if ref_base == read_base:  # identical match
                    pass
                else:  # different match
                    NM += 1

    return NM


def set_NMMD(read, fasta):
    ref_seq_padded, read_seq_padded = get_padded_seqs(read, fasta)
    NM = get_NM(read, ref_seq_padded=ref_seq_padded, read_seq_padded=read_seq_padded)
    MD = get_MD(read, ref_seq_padded=ref_seq_padded, read_seq_padded=read_seq_padded)
    read.set_tag(tag='NM', value_type='i', value=NM)
    read.set_tag(tag='MD', value_type='Z', value=MD)


def get_NMMD(read, fasta):
    ref_seq_padded, read_seq_padded = get_padded_seqs(read, fasta)
    NM = get_NM(read, ref_seq_padded=ref_seq_padded, read_seq_padded=read_seq_padded)
    MD = get_MD(read, ref_seq_padded=ref_seq_padded, read_seq_padded=read_seq_padded)
    return NM, MD


def get_pairorient_substring(read, mate=False):
    if mate:
        orientation = "r" if read.mate_is_reverse else "f"
        read12 = "2" if read.is_read1 else "1"
    else:
        orientation = "r" if read.is_reverse else "f"
        read12 = "1" if read.is_read1 else "2"

    return orientation + read12


def get_pairorient(read):
    if read.mate_is_unmapped:
        pairorient = None
    else:
        if read.template_length == 0:  # This occurs when mate chroms differ
            pairorient = None
        else:
            if read.template_length > 0:
                substring1 = get_pairorient_substring(read, mate=False)
                substring2 = get_pairorient_substring(read, mate=True)
            elif read.template_length < 0:
                substring1 = get_pairorient_substring(read, mate=True)
                substring2 = get_pairorient_substring(read, mate=False)

            pairorient = substring1 + substring2

    return pairorient


#####


def check_TRA(read):
    return read.reference_name != read.next_reference_name


def check_primary_alignment(read):
    return (not read.is_secondary) and (not read.is_supplementary)


#####


def get_fiveprime_end(read):
    if read.is_reverse:
        fiveprime_end = read.reference_end - 1
    else:
        fiveprime_end = read.reference_start

    return fiveprime_end


def get_threeprime_end(read):
    if read.is_reverse:
        threeprime_end = read.reference_start
    else:
        threeprime_end = read.reference_end - 1

    return threeprime_end


#####


def get_template_range0(read, with_softclip=False):
    """Returns a range object representing the interval between 5' ends of
        the read pairs.
    It is in 0-based, half-open format (bed format).
    """

    if read.template_length == 0:
        """
        TLEN == 0 cases include
            - different chromosomes between read pairs
            - two reads are on different strands and have identical 5' end
                base position
        """
        template_range0 = None
    else:
        if read.template_length > 0:
            start = get_fiveprime_end(read)
            end = start + read.template_length
        elif read.template_length < 0:
            end = get_fiveprime_end(read) + 1  # bed format
            start = end + read.template_length

        template_range0 = range(start, end)

    return template_range0


def get_softclip_ends_range0(read):
    if read.cigartuples[0][0] == 4:  # first cigar is S
        start = read.reference_start - read.cigartuples[0][1]
    else:
        start = read.reference_start

    if read.cigartuples[-1][0] == 4:  # last cigar is S
        end = read.reference_end + read.cigartuples[-1][1]
    else:
        end = read.reference_end

    return range(start, end)


Clipspec = collections.namedtuple(
    "Clipspec", ["start1", "is_forward", "seq", "qual", "qname"]
)


def get_softclip_specs(read):
    clipspec_list = list()
    if read.cigartuples[0][0] == 4:
        start1 = read.reference_start
        cliplen = read.cigartuples[0][1]
        seq = read.query_sequence[:cliplen][::-1]
        qual = list(read.query_qualities)[:cliplen][::-1]
        is_forward = False
        qname = read.query_name
        clipspec_list.append(Clipspec(start1, is_forward, seq, qual, qname))

    if read.cigartuples[-1][0] == 4:
        start1 = read.reference_end + 1
        cliplen = read.cigartuples[-1][1]
        seq = read.query_sequence[-cliplen:]
        qual = list(read.query_qualities)[-cliplen:]
        is_forward = True
        qname = read.query_name
        clipspec_list.append(Clipspec(start1, is_forward, seq, qual, qname))

    return clipspec_list


#####


def check_mateless(read):
    return (not read.is_paired) or read.mate_is_unmapped


def get_primary_mate_candidate(bam, read):
    primary_mate_candidate = list()
    for new_read in bam.fetch(
        read.next_reference_name,
        read.next_reference_start,
        read.next_reference_start + 1,
    ):
        if check_primary_alignment(new_read):
            if (new_read.query_name == read.query_name) and (
                read.compare(new_read) != 0
            ):
                # read.compare method returns 0 if the two reads are identical
                primary_mate_candidate.append(new_read)

    return primary_mate_candidate


def get_primary_mate(read, bam):
    """Returns:
        None if the input read is unpaired or the mate is unmapped.
        Else:
            1) Performs fetch with input bam in a position specified by RNEXT
                and PNEXT of the input read.
            2) From the fetch result, picks reads with the same QNAME as the
                input read.
            3) If there are only one picked read, returns it.
                Otherwise, raises an exception.
    """

    if check_mateless(read):
        mate = None
    else:
        primary_mate_candidate = get_primary_mate_candidate(bam, read)

        if len(primary_mate_candidate) == 0:
            raise Exception(
                f"No primary mate candidate found for this read:\n"
                f"{read.to_string()}"
            )
        elif len(primary_mate_candidate) == 1:
            mate = primary_mate_candidate[0]
        elif len(primary_mate_candidate) > 1:
            mate_candidate_strings = "\n".join(
                [read.to_string() for read in primary_mate_candidate]
            )
            raise Exception(
                f"More than one primary mate candidates found for "
                f"the input read:\n"
                f"{read.to_string()}\n"
                f"Detected mate candidates for this read:\n"
                f"{mate_candidate_strings}"
            )

    return mate




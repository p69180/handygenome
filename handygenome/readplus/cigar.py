import re
import collections

import importlib
top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))


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


CigarWalk = collections.namedtuple(
    'CigarWalk', 
    ('cigartuple', 'target_range0', 'query_range0'),
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

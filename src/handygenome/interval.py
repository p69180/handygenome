import functools

import pyranges as pr

import handygenome.tools as tools
import handygenome.deco as deco


class Interval:
    """
    Attributes:
        chrom
        start0
        end0
        start1
        end1
        range0
        length
    """

    def __init__(
        self, chrom, *, start1=None, end1=None, start0=None, end0=None, is_reverse=False,
    ):
        """Args:
            'chrom' is mandatory.
            ('start1' and 'end1') or ('start0' and 'end0') must 
                be given (for coordinate setting).
            'start0' and 'end0' are 0-based half-open system.
            'start1' and 'end1' are 1-based closed system.
        """

        self.chrom = chrom
        self.is_reverse = is_reverse

        # start0, end0, start1, end1
        if (start1 is not None) and (end1 is not None):
            self.start0 = start1 - 1
            self.end0 = end1
        else:
            self.start0 = start0
            self.end0 = end0

        self.start1 = self.start0 + 1
        self.end1 = self.end0

        # range
        self.range0 = range(self.start0, self.end0)

        # length
        self.length = self.end0 - self.start0

    def __repr__(self):
        return (f'<Interval> {self.chrom}:{self.start1:,}-{self.end1:,} '
                f'(1-based)')

    def to_gr(self, **kwargs):
        data = {
            'Chromosome': [self.chrom],
            'Start': [self.start0],
            'End': [self.end0],
        }
        for key, val in kwargs.items():
            data[key] = [val]
        return pr.from_dict(data)

    @classmethod
    def from_gr(cls, gr):
        assert len(gr) == 1
        return cls(chrom=gr.Chromosome[0], start0=gr.Start[0], end0=gr.End[0])

    def includes(self, other):
        return (self.chrom == other.chrom and
                self.start0 <= other.start0 and
                self.end0 >= other.end0)


class IntervalList(list):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._lengths_cumsum = None

    # constructors #
    @classmethod
    def from_gr(cls, gr):
        result = cls()
        for chrom, start0, end0 in zip(gr.Chromosome, gr.Start, gr.End):
            intv = Interval(chrom=chrom, start0=start0, end0=end0)
            result.append(intv)
        return result

    @classmethod
    def from_bed(cls, bedfile):
        gr = pr.read_bed(bedfile)
        return cls.from_gr(gr)

    @classmethod
    def from_chromdict(cls, chromdict):
        result = cls()
        for contig, length in chromdict.items():
            result.append(Interval(chrom=contig, start0=0, end0=length))

        return result

    @classmethod
    def from_vcfspec(cls, vcfspec):
        return cls.from_gr(cls, vcfspec.to_gr())

    @classmethod
    def from_margin(cls, chromdict, chrom_left, start0_left, chrom_right, end0_right):
        #chromdict = DEFAULT_CHROMDICTS[refver]
        cumpos0_left = chromdict.get_cumpos0(chrom_left, start0_left)
        cumpos0_right = chromdict.get_cumpos0(chrom_right, end0_right)
        if cumpos0_left >= cumpos0_right:
            raise Exception(
                f'"left position" comes later than "right_position"; '
                f'chrom_left={chrom_left}, start0_left={start0_left}, '
                f'chrom_right={chrom_right}, end0_right={end0_right}'
            )

        result = cls()
        if chrom_left == chrom_right:
            result.append(Interval(chrom=chrom_left, start0=start0_left, end0=end0_right))
        else:
            chrom_left_idx = chromdict.contigs.index(chrom_left)
            chrom_right_idx = chromdict.contigs.index(chrom_right)
            result.append(
                Interval(chrom=chrom_left, start0=start0_left, end0=chromdict[chrom_left])
            )
            for chrom in chromdict.contigs[(chrom_left_idx + 1):chrom_right_idx]:
                result.append(
                    Interval(chrom=chrom, start0=0, end0=chromdict[chrom])
                )
            result.append(
                Interval(chrom=chrom_right, start0=0, end0=end0_right)
            )

        return result

    @classmethod
    def get_depth_bins(cls, chromdict, width=100_000):
        result = cls()
        for contig, length in zip(chromdict.contigs, chromdict.lengths):
            intvlist_contig = cls()
            intvlist_contig.append(Interval(chrom=contig, start0=0, end0=length))
            for item in intvlist_contig.split(width=width):
                result.extend(item)
        return result
    ################

    def write_bed(self, outfile_path):
        with tools.openfile(outfile_path, 'w') as outfile:
            for intv in self:
                outfile.write(f'{intv.chrom}\t{intv.start0}\t{intv.end0}\n')

    def to_gr(self):
        chroms = list()
        starts = list()
        ends = list()
        for intv in self:
            chroms.append(intv.chrom)
            starts.append(intv.start0)
            ends.append(intv.end0)
        return pr.from_dict({'Chromosome': chroms, 'Start': starts, 'End': ends})

    # properties
    @functools.cached_property
    def lengths_cumsum(self):
        return list(itertools.accumulate(intv.length for intv in self))
    @property
    def length(self):
        return sum(intv.length for intv in self)

    # inclusion check
    def includes_vr(self, vr):
        vr_intv = Interval(chrom=vr.contig, start1=vr.pos, end1=vr.pos)
        return any(intv.includes(vr_intv) for intv in self)

    # calculations
    def sort_intervals(self, chromdict):
        def sortkey(intv):
            return tools.coord_sortkey(intv.chrom, intv.start1, chromdict)
        self.sort(key=sortkey)

    def isec(self, other):
        self_gr = self.to_gr()
        other_gr = other.to_gr()
        isec_gr = self_gr.set_intersect(other_gr)

        return self.__class__.from_gr(isec_gr)

    def subtract(self, other):
        self_gr = self.to_gr()
        other_gr = other.to_gr()
        subtract_gr = self_gr.subtract(other_gr)

        return self.__class__.from_gr(subtract_gr)

    def union(self, other):
        self_gr = self.to_gr()
        other_gr = other.to_gr()
        union_gr = self_gr.set_union(other_gr)

        return self.__class__.from_gr(union_gr)

    def merge(self):
        return self.__class__.from_gr(self.to_gr().merge())

    @deco.get_deco_num_set(('b', 'l', 'r'), 1)
    def slop(self, chromdict, b=None, l=None, r=None):
        def start_handler(start0, width):
            new_start0 = max(0, start0 - width)
            return new_start0

        def end_handler(end0, width, chrom, chromdict):
            new_end0 = min(chromdict[chrom], end0 + width)
            return new_end0

        result = self.__class__()
        for intv in self:
            if b is not None:
                new_start0 = start_handler(intv.start0, b)
                new_end0 = end_handler(intv.end0, b, intv.chrom, chromdict)
            elif l is not None:
                new_start0 = start_handler(intv.start0, l)
                new_end0 = intv.end0
            elif r is not None:
                new_start0 = new_start0
                new_end0 = end_handler(intv.end0, b, intv.chrom, chromdict)

            new_intv = Interval(chrom=intv.chrom, start0=new_start0, end0=new_end0)
            result.append(new_intv)

        return result

    @deco.get_deco_num_set(('num', 'width'), 1)
    def split(self, *, num=None, width=None):
        #self.sort_intervals(chromdict)

        # get result_lengths_cumsum
        total_length = self.lengths_cumsum[-1]
        if num is not None:
            result_lengths_cumsum = list(
                itertools.accumulate(get_interval_lengths_num(total_length, num))
            )
        elif width is not None:
            result_lengths_cumsum = list(itertools.accumulate(
                get_interval_lengths_width(total_length, width)))

        # prepare result values
        result = list()
        for idx in range(len(result_lengths_cumsum)):
            # get start0 and end0 in the single merged coordinate system
            merged_start0 = (0 
                             if idx == 0 else
                             result_lengths_cumsum[idx - 1])
            merged_end0 = result_lengths_cumsum[idx]
            # modify merged coordinates into interval-specific ones
            (chrom_start, 
             pos0_start, 
             self_idx_start) = self._modify_coord(merged_start0, is_end=False)
            (chrom_end, 
             pos0_end, 
             self_idx_end) = self._modify_coord(merged_end0, is_end=True)
            # create an IntervalList corresponding to a split unit
            intvlist = IntervalList()
            if self_idx_start == self_idx_end:
                intv = Interval(chrom_start, start0=pos0_start, end0=pos0_end)
                intvlist.append(intv)
            else:
                intv_first = Interval(chrom_start, 
                                      start0=pos0_start, 
                                      end0=self[self_idx_start].end0)
                intvlist.append(intv_first)

                for self_idx in range(self_idx_start + 1, self_idx_end):
                    intvlist.append(self[self_idx])

                intv_last = Interval(chrom_end, 
                                     start0=self[self_idx_end].start0, 
                                     end0=pos0_end)
                intvlist.append(intv_last)

            result.append(intvlist)

        return result

    def _modify_coord(self, merged_pos0, is_end=False):
        def get_interval_values(merged_pos0, idx, self):
            chrom = self[idx].chrom

            if idx == 0:
                shift_within_interval = merged_pos0
            else:
                shift_within_interval = (merged_pos0 
                                         - self.lengths_cumsum[idx - 1])
            new_pos0 = self[idx].start0 + shift_within_interval

            return chrom, new_pos0

        def handle_current_interval(merged_pos0, idx, length_previous, 
                                    length_current, is_end):
            if merged_pos0 > length_previous and merged_pos0 <= length_current:
                if merged_pos0 == length_current:
                    if is_end:
                        chrom, new_pos0 = get_interval_values(
                            merged_pos0, idx, self)
                    else:
                        chrom, new_pos0 = get_interval_values(
                            merged_pos0, idx + 1, self)
                elif merged_pos0 < length_current:
                    chrom, new_pos0 = get_interval_values(merged_pos0, idx, 
                                                          self)

                to_break = True
            else:
                chrom = None
                new_pos0 = None
                to_break = False

            return chrom, new_pos0, to_break

        # sanity check
        assert merged_pos0 >= 0, f'"merged_pos0" must be non-negative.'
        assert not ((not is_end) and merged_pos0 >= self.lengths_cumsum[-1]), (
            f'If ""is_end" is False, "merged_pos0" must be less than '
            f'the total IntervalList length.')
        assert not (is_end and merged_pos0 > self.lengths_cumsum[-1]), (
            f'If ""is_end" is True, "merged_pos0" must be less than or '
            f'equal to the total IntervalList length.')

        if merged_pos0 == 0:
            chrom = self[0].chrom
            new_pos0 = self[0].start0
            self_idx = 0
        else:
            for idx in range(len(self.lengths_cumsum)):
                length_current = self.lengths_cumsum[idx]
                length_previous = (0 
                                   if idx == 0 else
                                   self.lengths_cumsum[idx - 1])
                (chrom, new_pos0, to_break) = handle_current_interval(
                    merged_pos0, idx, length_previous, length_current, is_end)
                self_idx = idx

                if to_break:
                    break

        return chrom, new_pos0, self_idx


# interval-related functions

def get_interval_lengths_num(total_length, num):
    interval_num = min(num, total_length)
    q, r = divmod(total_length, interval_num)
    interval_width_1 = q + 1
    interval_num_1 = r
    interval_width_2 = q
    interval_num_2 = interval_num - r

    result = list(
        itertools.chain(
            itertools.repeat(interval_width_1, interval_num_1),
            itertools.repeat(interval_width_2, interval_num_2),
        )
    )

    return result


def get_interval_lengths_width(total_length, width):
    interval_width_raw = min(width, total_length)
    q, r = divmod(total_length, interval_width_raw)
    if r == 0:
        interval_width = interval_width_raw
        interval_num = q
        result = list(itertools.repeat(interval_width, interval_num))
    else:
        q2, r2 = divmod(r, q)
        interval_width_1 = interval_width_raw + q2 + 1
        interval_num_1 = r2
        interval_width_2 = interval_width_raw + q2
        interval_num_2 = q - r2
        result = list(itertools.chain(
            itertools.repeat(interval_width_1, interval_num_1),
            itertools.repeat(interval_width_2, interval_num_2)))

    return result


# range operation

def range_included(rng1, rng2):
    """Returns:
        True if rng1 is included in rng2
    """
    return (rng1.start >= rng2.start) and (rng1.stop <= rng2.stop)


def range_intersection(rng1, rng2):
    start = max(rng1.start, rng2.start)
    stop = min(rng1.stop, rng2.stop)
    if stop > start:
        return range(start, stop)
    else:
        return None


def range_subtraction(rng1, rng2):
    isec = range_intersection(rng1, rng2)
    if isec is None:
        return rng1
    else:
        if isec.start == rng1.start:
            if isec.stop == rng1.stop:
                return None
            else:
                return range(isec.stop, rng1.stop)
        else:
            if isec.stop == rng1.stop:
                return range(rng1.start, isec.start)
            else:
                return (range(rng1.start, isec.start), range(isec.stop, rng1.stop))


def reverse_range(rng):
    if rng.step == 1:
        start = rng.stop - 1
        stop = rng.start - 1
        return range(start, stop, -1)
    elif rng.step == -1:
        start = rng.stop + 1
        stop = rng.start + 1
        return range(start, stop, 1)
    else:
        raise Exception(f'Input range obejct must have "step" value of either 1 or -1.')


def check_overlaps_forward_nonzero(rng1, rng2):
    return (rng1.start < rng2.stop) and (rng1.stop > rng2.start)


def check_overlaps(rng1, rng2):
    """
    Meaning of range
        - Assumes range.step is either 1 or -1
        - "current position" travels, beginning from "range.start" toward "range.stop".
        - "cursor" position
            - When range.step == 1: immediately left to "current position"
            - When range.step == -1: immediately right to "current position"
        - Area swept by "cursor" is the area represented by range object.
        - When len(range) == 0:
            - When range.step == 1: 0-length interval between (range.start - 1) and (range.start)
            - When range.step == -1: 0-length interval between (range.start) and (range.start + 1)

    Meaning of overlap with 0-length range
        - consider as overlap when the interspace is contained between ALIGNED read bases (not insertion)

    Args:
        rng1, rng2: python range object.
    """
    # convert into positive-step range
    if rng1.step == -1:
        rng1 = reverse_range(rng1)
    if rng2.step == -1:
        rng2 = reverse_range(rng2)

    # main
    if len(rng1) == 0 and len(rng2) == 0:
        return rng1.start == rng2.start
    else:
        return check_overlaps_forward_nonzero(rng1, rng2)


def check_spans(query_range0, target_range0):
    """Checks if query spans target"""

    if len(query_range0) == 0:
        return False
    else:
        if len(target_range0) == 0:
            return check_overlaps(query_range0, target_range0)
        else:
            return (
                min(query_range0) <= min(target_range0) and
                max(query_range0) >= max(target_range0)
            )



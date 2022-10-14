import re
import itertools
import functools

import Bio

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
alignhandler = importlib.import_module('.'.join([top_package_name, 'align', 'alignhandler']))
realign = importlib.import_module('.'.join([top_package_name, 'align', 'realign']))


DEFAULT_FETCH_EXTEND_LENGTH = 10
SV_ALTS = ('DEL', 'INS', 'DUP', 'INV', 'CNV', 'BND', 'TRA')
CPGMET_ALT = 'CPGMET'


class Vcfspec:
    # constructors #
    def __init__(self, chrom=None, pos=None, ref=None, alts=None, 
                 somaticindex=1, germlineindexes=(0, 0)):
        if alts is not None:
            if not isinstance(alts, (tuple, list)):
                raise Exception(f'"alts" argument must be a tuple or a list.')

        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        if alts is not None:
            self.alts = tuple(alts)
        self.somaticindex = somaticindex
        self.germlineindexes = sorted(germlineindexes)

    @classmethod
    def from_vr(cls, vr):
        return cls(chrom=vr.contig, pos=vr.pos, ref=vr.ref, alts=vr.alts)
    ################

    def __repr__(self):
        if len(self.alts) == 1:
            altstring = str(self.alts[0])
        else:
            altstring = str(list(self.alts))
        return f'<Vcfspec ({self.chrom}:{self.pos} {self.ref}>{altstring})>'

    def __hash__(self):
        return hash(self.get_tuple())

    def __eq__(self, other):
        return all(
            getattr(self, key) == getattr(other, key)
            for key in ('chrom', 'pos', 'ref', 'alts')
        )

    @property
    def pos0(self):
        return self.pos - 1

    @property
    def end0(self):
        return self.pos0 + len(self.ref)

    @property
    def alleles(self):
        return (self.ref,) + self.alts

    @property
    def germline(self):
        alleles = self.alleles
        return tuple(alleles[x] for x in self.germlineindexes)

    @property
    def somatic(self):
        return self.alleles[self.somaticindex]

    def get_id(self):
        return '_'.join([self.chrom, 
                         str(self.pos), 
                         self.ref, 
                         '|'.join(self.alts)])

    def get_mutation_type(self, alt_index=0):
        return get_mttype(self.ref, self.alts[alt_index])

    def get_mttype_firstalt(self):
        return self.get_mutation_type(0)

    def get_tuple(self):
        return (self.chrom, self.pos, self.ref, self.alts)

    ### ranges
    @property
    def REF_range0(self):
        return range(self.pos0, self.end0)

    @functools.cache
    def get_preflank_range0(self, idx=0, flanklen=1):
        assert flanklen >= 1, f'"flanklen" argument must be at least 1.'

        if self.alts[idx][0] == self.ref[0]:
            flanklen = flanklen - 1
        return range(self.pos0 - flanklen, self.pos0)

    @functools.cache
    def get_postflank_range0(self, flanklen=1):
        assert flanklen >= 1, f'"flanklen" argument must be at least 1.'

        return range(self.pos0 + len(self.ref),
                     self.pos0 + len(self.ref) + flanklen)

    # equivalents
    def normalize(self, fasta, fetch_extend_length=DEFAULT_FETCH_EXTEND_LENGTH):
        return leftmost(self, fasta, fetch_extend_length=fetch_extend_length)

    # misc
    def apply_to_vr(self, vr):
        vr.contig = self.chrom
        vr.pos = self.pos
        vr.ref = self.ref
        vr.alts = self.alts

    def iter_monoalts(self):
        for alt in self.alts:
            new_vcfspec = self.__class__(self.chrom, self.pos, self.ref, (alt,))
            yield new_vcfspec

    def get_monoalt(self, alt_index=0):
        return self.__class__(
            self.chrom, self.pos, self.ref, (self.alts[alt_index],)
        )

    def check_without_N(self):
        return (
            without_N(self.ref) and
            all(without_N(x) for x in self.alts)
        )

    def to_hgvsg(self, alt_index=0):
        chrom = self.chrom
        pos = self.pos
        ref = self.ref
        alt = self.alts[alt_index]
        mttype = self.get_mutation_type(alt_index)

        if mttype == 'snv':
            result = f'{chrom}:g.{pos}{ref}>{alt}'
        elif mttype == 'mnv':
            pos2 = pos + (len(ref) - 1)
            result = f'{chrom}:g.{pos}_{pos2}delins{alt}'
        elif mttype == 'ins':
            inserted_seq = alt[1:]
            result = f'{chrom}:g.{pos}_{pos+1}ins{inserted_seq}'
        elif mttype == 'del':
            pos1 = pos + 1
            pos2 = pos + (len(ref) - 1)
            if pos1 == pos2:
                result = f'{chrom}:g.{pos1}del'
            else:
                result = f'{chrom}:g.{pos1}_{pos2}del'
        elif mttype == 'delins':
            pos1 = pos
            pos2 = pos + (len(ref) - 1)
            if pos1 == pos2:
                result = f'{chrom}:g.{pos1}delins{alt}'
            else:
                result = f'{chrom}:g.{pos1}_{pos2}delins{alt}'

        return result

    def to_gr(self):
        ref_range0 = self.REF_range0
        return pr.from_dict(
            {
                'Chromosome': [self.chrom], 
                'Start': [ref_range0.start], 
                'End': [ref_range0.stop],
            }
        )


def check_vcfspec_monoalt(vcfspec):
    if len(vcfspec.alts) != 1:
        raise Exception('The input vcfspec must be with single ALT.')

check_vcfspec_monoallele = check_vcfspec_monoalt


def get_mttype(ref, alt):
    if common.RE_PATS['nucleobases'].fullmatch(alt) is None:
        if any(
            (re.fullmatch(f'<{x}(:.+)?>', alt) is not None) for x in SV_ALTS
        ):
            mttype = 'sv'
        elif (
            (common.RE_PATS['alt_bndstring_1'].fullmatch(alt) is not None) or 
            (common.RE_PATS['alt_bndstring_2'].fullmatch(alt) is not None)
        ):
            mttype = 'sv'
        elif alt == f'<{CPGMET_ALT}>':
            mttype = 'cpgmet'
        else:
            raise Exception(f'Unexpected symbolic ALT allele: {alt}')
    else:
        if len(ref) == len(alt):
            if len(ref) == 1:
                mttype = 'snv'
            else:
                mttype = 'mnv'
        else:
            if len(ref) == 1:
                if ref[0] == alt[0]:
                    mttype = 'ins'
                else:
                    mttype = 'delins'
            elif len(alt) == 1:
                if ref[0] == alt[0]:
                    mttype = 'del'
                else:
                    mttype = 'delins'
            else:
                mttype = 'delins'

    return mttype


def merge_vcfspecs(vcfspec_list, fasta, merging_distance_le=3):
    assert len(set(x.chrom for x in vcfspec_list)) == 1, f'Input vcfspecs belong to more than one chromosomes.'
    assert all(len(x.alts) == 1 for x in vcfspec_list), f'All input vcfspecs must have one ALT.'

    new_vcfspec_list = list()

    def unit_func(vcfspec1, vcfspec2):
        num_intervening_bases = vcfspec2.REF_range0.start - vcfspec1.REF_range0.stop
        if num_intervening_bases <= merging_distance_le:
            new_pos = vcfspec1.pos
            if num_intervening_bases == 0:
                new_ref = vcfspec1.ref + vcfspec2.ref
                new_alt = vcfspec1.alts[0] + vcfspec2.alts[0]
            else:
                intervening_ref = fasta.fetch(vcfspec1.chrom, vcfspec1.REF_range0.stop, vcfspec2.REF_range0.start)
                new_ref = vcfspec1.ref + intervening_ref + vcfspec2.ref
                new_alt = vcfspec1.alts[0] + intervening_ref + vcfspec2.alts[0]
            return Vcfspec(vcfspec1.chrom, new_pos, new_ref, (new_alt,))
        else:
            new_vcfspec_list.append(vcfspec1)
            return vcfspec2

    final = functools.reduce(unit_func, sorted(vcfspec_list, key=(lambda x: x.pos)))
    new_vcfspec_list.append(final)

    return new_vcfspec_list


def split_vcfspec(vcfspec, fasta, splitting_distance_ge=3, show_alignment=False):
    check_vcfspec_monoalt(vcfspec)

    alignment = alignhandler.alignment_tiebreaker(
        realign.ALIGNER_MAIN.align(vcfspec.ref, vcfspec.alts[0]),
        raise_with_failure=True,
    )
    vcfspec_list = alignhandler.alignment_to_vcfspec(
        alignment, target_start0=vcfspec.pos0, chrom=vcfspec.chrom, fasta=fasta, strip_query_gaps=False)
    merged_vcfspec_list = merge_vcfspecs(vcfspec_list, fasta, merging_distance_le=(splitting_distance_ge - 1))
    return merged_vcfspec_list


##############
# equivalent #
##############


def leftmost(vcfspec, fasta, fetch_extend_length=DEFAULT_FETCH_EXTEND_LENGTH):
    return indel_equivalents(vcfspec, fasta, fetch_extend_length)[0]


def rightmost(vcfspec, fasta, fetch_extend_length=DEFAULT_FETCH_EXTEND_LENGTH):
    return indel_equivalents(vcfspec, fasta, fetch_extend_length)[-1]


def check_equivalent(vcfspec1, vcfspec2, fasta):
    return leftmost(vcfspec1, fasta) == leftmost(vcfspec2, fasta)


def indel_equivalents(vcfspec, fasta, parsimoniously=True, fetch_extend_length=DEFAULT_FETCH_EXTEND_LENGTH):
    if parsimoniously:
        vcfspec = make_parsimonious(vcfspec)
    return shifting_equivalents(vcfspec=vcfspec, fasta=fasta, fetch_extend_length=fetch_extend_length)


def make_parsimonious(vcfspec):
    check_vcfspec_monoalt(vcfspec)

    pos = vcfspec.pos
    ref = vcfspec.ref
    alt = vcfspec.alts[0]
    #ref_alt_samelen = (len(ref) == len(alt))

    # left
    min_length = min(len(ref), len(alt))
    action_idx = None
    for idx in range(0, min_length):
        # all identical leading bases are searched for
        if ref[idx] == alt[idx]:
            action_idx = idx
            continue
        else:
            break

    if action_idx is not None:
        ref = ref[action_idx:]
        alt = alt[action_idx:]
        pos += action_idx
        # identical leading bases, except the rightmost one, have been removed

    # right
    min_length = min(len(ref), len(alt))
    action_idx = None
    for idx in range(-1, -min_length, -1):
        # all identical trailing bases, except for the leftmost position, are searched for
        if ref[idx] == alt[idx]:
            action_idx = idx
            continue
        else:
            break

    if action_idx is not None:
        ref = ref[:action_idx]
        alt = alt[:action_idx]
        # all searched identical trailing bases have been removed

    # final modificaiton
    if (
        (len(ref) == len(alt)) or  # e.g. REF = AGCC ; ALT = ATCG
        (len(ref) > 1 and len(alt) > 1)  # e.g. REF = ATTTGA ; ALT = AG
    ):
        ref = ref[1:]
        alt = alt[1:]
        pos += 1

    # result
    return Vcfspec(vcfspec.chrom, pos, ref, (alt,))


def shifting_equivalents(vcfspec, fasta, fetch_extend_length):
    """Args:
        vcfspec: Must have single ALT.

    Returns:
        [vcfspec, vcfspec, ...]
            - each vcfspec is equivalent to the input vcfspec
            - each vcfspec is the most parsimonious form
            - sorted by pos (increasing order)

        Return value includes only input 'vcfspec' itself if one or more
            of the below is true:
            1) ref or alt sequence includes 'n' or 'N'
            2) input vcfspec is not a simple insertion or deletion (e.g. SNV, MNV, complex indel(delins))
    """

    ##### progressive fetch #####

    def progressive_fetch_forward(
        chrom,
        pos,
        fasta,
        fetch_extend_length,
        diff_seq,
        is_ins,
    ):
        extend_count = 0
        out_of_range = False
        chrom_length = fasta.lengths[fasta.references.index(chrom)]
        diff_seq_length = len(diff_seq)
        fetch_range_setter = (
            fetch_range_setter_ins_forward if is_ins else fetch_range_setter_del_forward
        )
        while True:
            if out_of_range:
                break

            extend_count += 1
            fetch_start, fetch_end = fetch_range_setter(
                pos, extend_count, fetch_extend_length, diff_seq_length
            )
            if fetch_end > chrom_length:
                out_of_range = True
                if fetch_start >= chrom_length:
                    continue

            result = fasta.fetch(chrom, fetch_start, fetch_end)
            yield result

    def progressive_fetch_backward(
        chrom,
        pos,
        fasta,
        fetch_extend_length,
    ):
        extend_count = 0
        out_of_range = False
        while True:
            if out_of_range:
                break

            extend_count += 1
            fetch_start, fetch_end = fetch_range_setter_backward(
                pos, extend_count, fetch_extend_length
            )
            if fetch_start < 0:
                out_of_range = True
                fetch_start = 0
                if fetch_end <= 0:
                    continue

            result = reversed(list(fasta.fetch(chrom, fetch_start, fetch_end)))
            yield result

    def fetch_range_setter_ins_forward(
        pos,
        extend_count,
        fetch_extend_length,
        diff_seq_length,
    ):
        fetch_start = pos + ((extend_count - 1) * fetch_extend_length)
        fetch_end = pos + (extend_count * fetch_extend_length)
        return fetch_start, fetch_end

    def fetch_range_setter_del_forward(
        pos,
        extend_count,
        fetch_extend_length,
        diff_seq_length,
    ):
        fetch_start = pos + diff_seq_length + ((extend_count - 1) * fetch_extend_length)
        fetch_end = pos + diff_seq_length + (extend_count * fetch_extend_length)
        return fetch_start, fetch_end

    def fetch_range_setter_backward(
        pos,
        extend_count,
        fetch_extend_length,
    ):
        fetch_start = pos - (extend_count * fetch_extend_length)
        fetch_end = pos - ((extend_count - 1) * fetch_extend_length)
        return fetch_start, fetch_end

    ##### navigate #####

    def navigate(
        chrom, pos, fasta, diff_seq, new_pos_list, fetch_extend_length, is_ins, forward
    ):
        """Args:
        is_ins: True for insertion, False for deletion
        forward: True for forward naviation (to 3' side on the
            positive strand), False for backward navigation
        """
        fetch_reservoir = list()
        fetched_seq = list()
        diff_seq_compared = diff_seq if forward else diff_seq[::-1]

        if forward:
            fetcher = progressive_fetch_forward(
                chrom,
                pos,
                fasta,
                fetch_extend_length,
                diff_seq,
                is_ins,
            )
            addnewpos = (
                navigate_helper_addnewpos_ins_forward
                if is_ins
                else navigate_helper_addnewpos_del_forward
            )
        else:
            fetcher = progressive_fetch_backward(
                chrom,
                pos,
                fasta,
                fetch_extend_length,
            )
            addnewpos = (
                navigate_helper_addnewpos_ins_backward
                if is_ins
                else navigate_helper_addnewpos_del_backward
            )

        for idx in itertools.count(0):
            navigate_helper_extend_fetch(fetch_reservoir, idx, diff_seq, fetcher)

            try:
                seq = fetch_reservoir[idx]
            except IndexError:
                break
            else:
                fetched_seq.append(seq)

            idx_remainder = idx % len(diff_seq)

            if fetched_seq[-1] != diff_seq_compared[idx_remainder]:
                break
            else:
                addnewpos(
                    idx,
                    idx_remainder,
                    fetch_reservoir,
                    pos,
                    diff_seq,
                    diff_seq_compared,
                    new_pos_list,
                    fetched_seq,
                )

    def navigate_helper_extend_fetch(fetch_reservoir, idx, diff_seq, fetcher):
        # "idx + 1 + len(diff_seq)" is the largest index of "fetch_reservoir", used in "navigate_helper_addnewpos_del_backward"
        limit_length = idx + 1 + len(diff_seq)
        while len(fetch_reservoir) <= limit_length:
            try:
                fetcher_output = next(fetcher)
            except StopIteration:
                break
            else:
                fetch_reservoir.extend(fetcher_output)

    def navigate_helper_addnewpos_ins_forward(
        idx,
        idx_remainder,
        fetch_reservoir,
        pos,
        diff_seq,
        diff_seq_compared,
        new_pos_list,
        fetched_seq,
    ):
        try:
            new_ref = fetch_reservoir[idx]
        except IndexError:
            pass
        else:
            new_pos = pos + (idx + 1)
            inserted_seq = (
                diff_seq_compared[idx_remainder + 1 :]
                + diff_seq_compared[: idx_remainder + 1]
            )
            new_alt = new_ref + "".join(inserted_seq)

            new_pos_list.append((new_pos, new_ref, new_alt))

    def navigate_helper_addnewpos_ins_backward(
        idx,
        idx_remainder,
        fetch_reservoir,
        pos,
        diff_seq,
        diff_seq_compared,
        new_pos_list,
        fetched_seq,
    ):
        try:
            new_ref = fetch_reservoir[idx + 1]
        except IndexError:
            pass
        else:
            new_pos = pos - (idx + 1)
            inserted_seq = list(
                reversed(
                    diff_seq_compared[idx_remainder + 1 :]
                    + diff_seq_compared[: idx_remainder + 1]
                )
            )
            new_alt = new_ref + "".join(inserted_seq)

            new_pos_list.append((new_pos, new_ref, new_alt))

    def navigate_helper_addnewpos_del_forward(
        idx,
        idx_remainder,
        fetch_reservoir,
        pos,
        diff_seq,
        diff_seq_compared,
        new_pos_list,
        fetched_seq,
    ):
        # overlapping deletion
        if idx < len(diff_seq):
            new_pos = pos + (idx + 1)
            deleted_seq = diff_seq_compared[idx + 1 :] + fetched_seq
            new_ref = diff_seq[idx] + "".join(deleted_seq)
            new_alt = new_ref[0]

            new_pos_list.append((new_pos, new_ref, new_alt))

        # non-overlapping deletion
        try:
            beyond_seq = get_beyond_seq(fetch_reservoir, idx, diff_seq)
        except IndexError:
            pass
        else:
            if nonoverlapping_deletion_check(
                diff_seq_compared, idx_remainder, beyond_seq
            ):
                new_pos = pos + len(diff_seq) + (idx + 1)
                deleted_seq = beyond_seq
                new_ref = fetch_reservoir[idx] + "".join(deleted_seq)
                # this line is not enclosed in try clause because the index used here is smaller than those used in "get_beyond_seq" function.
                new_alt = new_ref[0]

                new_pos_list.append((new_pos, new_ref, new_alt))

    def navigate_helper_addnewpos_del_backward(
        idx,
        idx_remainder,
        fetch_reservoir,
        pos,
        diff_seq,
        diff_seq_compared,
        new_pos_list,
        fetched_seq,
    ):
        # overlapping deletion
        if idx < len(diff_seq):
            try:
                seq = fetch_reservoir[idx + 1]
            except IndexError:
                pass
            else:
                new_pos = pos - (idx + 1)
                deleted_seq = list(reversed(diff_seq_compared[idx + 1 :] + fetched_seq))
                new_ref = seq + "".join(deleted_seq)
                new_alt = new_ref[0]

                new_pos_list.append((new_pos, new_ref, new_alt))

        # non-overlapping deletion
        try:
            beyond_seq = get_beyond_seq(fetch_reservoir, idx, diff_seq)
        except IndexError:
            pass
        else:
            if nonoverlapping_deletion_check(
                diff_seq_compared, idx_remainder, beyond_seq
            ):
                try:
                    seq = fetch_reservoir[idx + len(diff_seq) + 1]
                except IndexError:
                    pass
                else:
                    new_pos = pos - (idx + 1) - len(diff_seq)
                    deleted_seq = list(reversed(beyond_seq))
                    new_ref = seq + "".join(deleted_seq)
                    new_alt = new_ref[0]

                    new_pos_list.append((new_pos, new_ref, new_alt))

    def get_beyond_seq(fetch_reservoir, idx, diff_seq):
        if len(fetch_reservoir) < idx + 1 + len(diff_seq):
            raise IndexError
        else:
            return fetch_reservoir[(idx + 1) : (idx + 1 + len(diff_seq))]

    def nonoverlapping_deletion_check(diff_seq_compared, idx_remainder, beyond_seq):
        return (
            diff_seq_compared[idx_remainder + 1 :] == beyond_seq[: -(idx_remainder + 1)]
        ) and (
            diff_seq_compared[: idx_remainder + 1] == beyond_seq[-(idx_remainder + 1) :]
        )

    ##################################################

    def get_result(chrom, pos, ref, alt, fasta, diff_seq, fetch_extend_length, is_ins):
        # get new_pos_list
        new_pos_list = list()
        navigate(
            chrom,
            pos,
            fasta,
            diff_seq,
            new_pos_list,
            fetch_extend_length,
            is_ins,
            forward=True,
        )
        navigate(
            chrom,
            pos,
            fasta,
            diff_seq,
            new_pos_list,
            fetch_extend_length,
            is_ins,
            forward=False,
        )
        # make result
        result = list()
        result.append(Vcfspec(chrom, pos, ref, [alt]))
        for new_pos, new_ref, new_alt in new_pos_list:
            result.append(Vcfspec(chrom, new_pos, new_ref, [new_alt]))
        result.sort(key=(lambda x: x.pos))

        return result

    # sanity check
    def sanity_check_indel(vcfspec):
        if vcfspec.ref[0] != vcfspec.alts[0][0]:
            raise Exception(f'Inserted/deleted seq must be on the right side, not left side.')

    # main
    check_vcfspec_monoalt(vcfspec)

    chrom, pos, ref, alts = vcfspec.get_tuple()
    alt = alts[0]

    if "n" not in ref and "N" not in ref and "n" not in alt and "N" not in alt:
        # if set('nN').isdisjoint(ref) and set('nN').isdisjoint(alt):
        mttype = vcfspec.get_mttype_firstalt()
        if mttype in ("ins", "del"):
            sanity_check_indel(vcfspec)
            is_ins = mttype == "ins"
            diff_seq = list(alt[1:]) if is_ins else list(ref[1:])
            result = get_result(
                chrom, pos, ref, alt, fasta, diff_seq, fetch_extend_length, is_ins
            )
        else:
            result = [vcfspec]
    else:
        result = [vcfspec]

    return result


"""
Below functions are kept for demonstration of the basic concept.

def _check_commutative_for_proof(refseq, insseq):
    def cond1(refseq, insseq):
        n = 0
        while True:
            n += 1
            end = n*len(insseq)
            start = (n-1)*len(insseq)
            subseq = refseq[start:end]
            if subseq != insseq:
                return False

            if end in range(len(refseq) - len(insseq) + 1, len(refseq) + 1):
                break

        return end

    def cond2(refseq, insseq):
        return refseq[ -len(insseq) : ] == insseq

    def cond3(refseq, insseq, end):
        #assert end < len(refseq)
        idx = len(refseq) - end
        return refseq[ -idx : ] == insseq[ : idx ]

    def main(refseq, insseq):
        end = cond1(refseq, insseq)
        if end == len(refseq):
            return True
        else:
            if cond2(refseq, insseq):
                return cond3(refseq, insseq, end)
            else:
                return False

    return main(refseq, insseq)
"""

##############

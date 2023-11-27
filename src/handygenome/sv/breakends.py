
import re
import collections
import itertools
import logging
import functools
import operator

import Bio.Seq

import handygenome.refgenome.refgenome as refgenome
import handygenome.tools as tools
import handygenome.deco as deco
import handygenome.workflow as workflow
import handygenome.sv.structvars as structvars
import handygenome.variant.vcfspec as libvcfspec
from handygenome.variant.vcfspec import Vcfspec, NonStandardSVAlt


LOGGER = workflow.get_logger(
    name=__name__,
    formatter=logging.Formatter(
        fmt='[%(asctime)s  %(levelname)s] %(name)s - %(message)s',
        datefmt=workflow.DEFAULT_DATEFMT),
    level='info', stderr=True, filenames=None, append=False)


class BreakendSpec(
    collections.namedtuple('BreakendSpec', ('chrom', 'pos0', 'is5prime'))
):
    @property
    def pos(self):
        return self.pos0 + 1

    @property
    def is3prime(self):
        return not self.is5prime

    def spawn(self, **kwargs):
        kwargs = (
            {key: getattr(self, key) for key in ('chrom', 'pos0', 'is5prime')}
            | kwargs
        )
        return self.__class__(**kwargs)


class Breakends:
    """
    Attributes:
        chrom_bnd1
        pos_bnd1
        is5prime_bnd1
        chrom_bnd2
        pos_bnd2
        is5prime_bnd2
        inserted_seq: 
            - A list composed of the bases of the inserted sequence 
                (e.g. ['A', 'A', 'G', 'C'])
            - As seen in the viewpoint where bnd1-side sequence is on 
                the plus(Crick) strand

        fasta: pysam.FastaFile instance
        chromdict: ChromDict instance
        svtype
        score
        is3prime_bnd1
        is3prime_bnd2

        #_equivs: A list of Breakends objects with maximum reference coverage.
        _pos_range0_bnd1
        _pos_range0_bnd2
        #_homlen = None
        #_homseq = None
    """

    #BASIC_ATTRS = ('chrom_bnd1', 'pos_bnd1', 'is5prime_bnd1', 
    #               'chrom_bnd2', 'pos_bnd2', 'is5prime_bnd2', 
    #               'inserted_seq', 'svtype')

    @deco.get_deco_num_set_differently(
        (
            'inserted_seq', 
            'inserted_seq_view1',
            'inserted_seq_view2',
        ), 
        1, 
        'le',
    )
    @refgenome.deco_standardize
    def __init__(
        self, 
        bndspec1, bndspec2, refver,
        inserted_seq=None, 
        inserted_seq_view1=None,
        inserted_seq_view2=None,
        #chromdict=None, 
        #svtype=None, 
        score=None,
    ):
        """Args:
            bndspec1, bndspec2: tuple (chrom, pos0, is5prime)
            inserted_seq, view12: 
                - Only one of them must be used
                - If used, str or tuple
        """
        arg_bndspec1 = bndspec1
        arg_bndspec2 = bndspec2
        # sanitycheck
        order = tools.compare_coords(
            arg_bndspec1[0], arg_bndspec1[1], 
            arg_bndspec2[0], arg_bndspec2[1], 
            refgenome.get_chromdict(refver),
        )
        assert order != 0

        # basic attributes
        self.refver = refver
        if order == -1:
            arg1_first = True
        elif order == 1:
            arg1_first = False
        else:
            raise Exception()

        if arg1_first:
            self.bndspec1 = BreakendSpec(*arg_bndspec1)
            self.bndspec2 = BreakendSpec(*arg_bndspec2)
        else:
            self.bndspec1 = BreakendSpec(*arg_bndspec2)
            self.bndspec2 = BreakendSpec(*arg_bndspec1)

        # inserted seq
        self.init_inserted_seq(
            inserted_seq, 
            inserted_seq_view1,
            inserted_seq_view2,
            arg1_first,
        )

        # score
        if score is None:
            self.score = 0
        else:
            self.score = score

        # private attributes accesible with getters
        #self._pos_range0_bnd1 = None
        #self._pos_range0_bnd2 = None

    def init_inserted_seq(
        self, 
        inserted_seq, 
        inserted_seq_view1,
        inserted_seq_view2,
        arg1_first,
    ):
        assert all(
            isinstance(x, (str, tuple, list, type(None))) for x in 
            (
                inserted_seq, 
                inserted_seq_view1,
                inserted_seq_view2,
            )
        )

        if all(
            (x is None) for x in (
                inserted_seq, 
                inserted_seq_view1,
                inserted_seq_view2,
            )
        ):
            self.inserted_seq = tuple()
        elif inserted_seq is not None:
            self.inserted_seq = tuple(inserted_seq)
        elif inserted_seq_view1 is not None:
            if arg1_first:
                self.inserted_seq = tuple(inserted_seq_view1)
            else:
                if isinstance(inserted_seq_view1, tuple):
                    inserted_seq_view1 = ''.join(inserted_seq_view1)
                self.inserted_seq = tuple(
                    self._convert_seq_between_bnds(inserted_seq_view1)
                )
        elif inserted_seq_view2 is not None:
            if (not arg1_first):
                self.inserted_seq = tuple(inserted_seq_view2)
            else:
                if isinstance(inserted_seq_view2, tuple):
                    inserted_seq_view2 = ''.join(inserted_seq_view2)
                self.inserted_seq = tuple(
                    self._convert_seq_between_bnds(inserted_seq_view2)
                )

    def __repr__(self):
        def func1(self):
            result = [f'<{self.__class__.__name__}>']
            result.extend(
                [f'{key}: {getattr(self, key)}' 
                 for key in ('chrom_bnd1', 'pos_bnd1', 'is5prime_bnd1', 
                             'inserted_seq',
                             'chrom_bnd2', 'pos_bnd2', 'is5prime_bnd2')])
            return '\n'.join(result)

        def func2(self):
            brkt1 = '[' if self.is5prime_bnd1 else ']'
            brkt2 = '[' if self.is5prime_bnd2 else ']'
            bnd1 = f'{brkt1}{self.chrom_bnd1}:{self.pos_bnd1}{brkt1}'
            bnd2 = f'{brkt2}{self.chrom_bnd2}:{self.pos_bnd2}{brkt2}'
            insseq = ''.join(self.inserted_seq)

            return f'<{self.__class__.__name__} bnd1={bnd1} | bnd2={bnd2} | insseq={repr(insseq)}>'

        def func3(self):
            endtype1 = '5\'' if self.is5prime_bnd1 else '3\''
            endtype2 = '5\'' if self.is5prime_bnd2 else '3\''
            if len(self.inserted_seq) == 0:
                insseq = 'no inserted seq'
            else:
                insseq = ''.join(self.inserted_seq)

            return (f'<{self.__class__.__name__} '
                    f'{self.chrom_bnd1}:{self.pos_bnd1} {endtype1} '
                    f'[{insseq}] '
                    f'{self.chrom_bnd2}:{self.pos_bnd2} {endtype2}>')

        return func2(self)

    def __hash__(self):
        return hash(
            (self.bndspec1, self.bndspec2, self.inserted_seq, self.refver)
        )
        #return hash((self.chrom_bnd1, self.pos_bnd1, self.is5prime_bnd1,
        #             self.chrom_bnd2, self.pos_bnd2, self.is5prime_bnd2,
        #             tuple(self.inserted_seq)))

#    def __hash__(self):
#        bnd1adv_form = get_bnd1adv_form(self)
#        return hash((
#                    bnd1adv_form.chrom_bnd1,
#                    bnd1adv_form.pos_bnd1,
#                    bnd1adv_form.is5prime_bnd1,
#                    bnd1adv_form.chrom_bnd2,
#                    bnd1adv_form.pos_bnd2,
#                    bnd1adv_form.is5prime_bnd2,
#                    bnd1adv_form.inserted_seq,
#                    ))

    def __eq__(self, other):
        #assert isinstance(other, Breakends)
        return hash(self) == hash(other)

    ###########################################

    @classmethod
    def from_vcfspec(cls, vcfspec):
        try:
            (
                t,
                chrom_mate,
                pos_mate,
                current_is5prime,
                mate_is5prime,
            ) = libvcfspec.parse_sv_altstring(vcfspec.alts[0])
        except NonStandardSVAlt as e:
            raise Exception(f'Cannot interpret input vcfspec ALT string') from e

        chromdict = refgenome.get_chromdict(vcfspec.refver)
        is_bnd1 = get_is_bnd1(vcfspec, chrom_mate, pos_mate, chromdict)

        return get_bnds_from_vcfspec_parseinfo(
            chrom=vcfspec.chrom,
            pos=vcfspec.pos,
            ref=vcfspec.ref,
            t=t,
            chrom_mate=chrom_mate,
            pos_mate=pos_mate,
            current_is5prime=current_is5prime,
            mate_is5prime=mate_is5prime,
            is_bnd1=is_bnd1,
            refver=vcfspec.refver,
        )

    @classmethod
    def from_vr(cls, vr):
        assert len(vr.alts) == 1, (
            f'Multiallelic variant record is not allowed:\n{vr}'
        )
        refver = refgenome.infer_refver_vcfheader(vr.header)
        fasta = refgenome.get_fasta(refver)
        chromdict = refgenome.get_chromdict(refver)

        vr_svinfo = get_vr_svinfo_standard_vr(vr, fasta, chromdict)
        bnds = get_bnds_from_vr_svinfo(vr, vr_svinfo, fasta, chromdict)

        return bnds

    ##############
    # properties #
    ##############

    @property
    def chromdict(self):
        return refgenome.get_chromdict(self.refver)

    @property
    def fasta(self):
        return refgenome.get_fasta(self.refver)

    # bnd1 immutable

    @property
    def chrom_bnd1(self):
        return self.bndspec1.chrom

    @property
    def is5prime_bnd1(self):
        return self.bndspec1.is5prime

    @property
    def is3prime_bnd1(self):
        return self.bndspec1.is3prime

    @property
    def pos0_bnd1(self):
        return self.bndspec1.pos0

    @property
    def pos_bnd1(self):
        return self.bndspec1.pos

    # bnd2 immutable

    @property
    def chrom_bnd2(self):
        return self.bndspec2.chrom

    @property
    def is5prime_bnd2(self):
        return self.bndspec2.is5prime

    @property
    def is3prime_bnd2(self):
        return self.bndspec2.is3prime

    @property
    def pos0_bnd2(self):
        return self.bndspec2.pos0

    @property
    def pos_bnd2(self):
        return self.bndspec2.pos

    ###

    @property
    def svtype(self):
        if self.chrom_bnd1 != self.chrom_bnd2:
            self.svtype = 'TRA'
        else:
            if (not self.is5prime_bnd1) and self.is5prime_bnd2:
                self.svtype = 'DEL'
            elif self.is5prime_bnd1 and (not self.is5prime_bnd2):
                self.svtype = 'DUP'
            elif self.is5prime_bnd1 == self.is5prime_bnd2:
                self.svtype = 'INV'

    @functools.cached_property
    def microhomology_info(self):
        bnds_equivs = self.get_equivs()

        homlen = len(bnds_equivs) - 1
        if homlen == 0:
            homseq = ''
        else:
            bnds = bnds_equivs[0]
            if bnds.is5prime_bnd1:
                pos_list = [
                    x.pos_bnd1 for x in 
                    sorted(bnds_equivs, key = lambda x: x.pos_bnd1)[:-1]
                ]
            else:
                pos_list = [
                    x.pos_bnd1 for x in 
                    sorted(bnds_equivs, key = lambda x: x.pos_bnd1)[1:]
                ]
            homseq = bnds.fasta.fetch(bnds.chrom_bnd1, pos_list[0]-1, pos_list[-1])

        return homlen, homseq

    @property
    def homlen(self):
        return self.microhomology_info()[0]

    @property
    def homseq(self):
        return self.microhomology_info()[1]

    ###############
    # equivalents #
    ###############

    @functools.cache
    def get_equivalents(self):
        """Returns: A list of Breakends objects, all equivalent to the 
            input object, sorted such that the first item is the most advanced 
            form with respect to bnd1, and the last item is the most advanced 
            form with respect to bnd2.
        """
        self_copy = self.spawn()

        bnds_list_bnd1adv = list()
        bnds_list_bnd1adv.append(self_copy)
        while True:
            #new = bnds_list_bnd1adv[-1].spawn()
            new_bnds, success = bnds_list_bnd1adv[-1].advance_bnd1()
            if success:
                bnds_list_bnd1adv.append(new_bnds)
                continue
            else:
                break

        bnds_list_bnd2adv = list()
        bnds_list_bnd2adv.append(self_copy)
        while True:
            #new = bnds_list_bnd2adv[-1].spawn()
            new_bnds, success = bnds_list_bnd2adv[-1].advance_bnd2()
            if success:
                bnds_list_bnd2adv.append(new_bnds)
                continue
            else:
                break

        bnds_list_bnd1adv.reverse()
        bnds_list = bnds_list_bnd1adv[:-1] + [self_copy] + bnds_list_bnd2adv[1:]
        max_score = max(x.score for x in bnds_list)

        bnds_equivs = [x for x in bnds_list if x.score == max_score]

        return bnds_equivs

    get_equivs = get_equivalents

    def get_sorted_equivalents(self, mode):
        assert mode in ('bnd1adv_first', 'bnd1adv_last', 'bnd2adv_first', 
                        'bnd2adv_last'), 'Invalid "mode" argument'

        bnds_equivs = self.get_equivalents()

        if mode == 'bnd1adv_first':
            reverse = (not bnds_equivs[0].is5prime_bnd1)
            #key = lambda bnds: bnds.pos_bnd1
            key = operator.attrgetter('pos_bnd1')
        elif mode == 'bnd1adv_last':
            reverse = bnds_equivs[0].is5prime_bnd1
            #key = lambda bnds: bnds.pos_bnd1
            key = operator.attrgetter('pos_bnd1')
        elif mode == 'bnd2adv_first':
            reverse = (not bnds_equivs[0].is5prime_bnd2)
            #key = lambda bnds: bnds.pos_bnd2
            key = operator.attrgetter('pos_bnd2')
        elif mode == 'bnd2adv_last':
            reverse = bnds_equivs[0].is5prime_bnd2
            #key = lambda bnds: bnds.pos_bnd2
            key = operator.attrgetter('pos_bnd2')

        return sorted(bnds_equivs, key=key, reverse=reverse)

    get_sorted_equivs = get_sorted_equivalents

    #################
    # advanced form #
    #################

    def get_bnd1adv_form(self):
        bnds_equivs = self.get_equivalents()
        if bnds_equivs[0].is5prime_bnd1:
            return min(bnds_equivs, key=operator.attrgetter('pos_bnd1'))
        else:
            return max(bnds_equivs, key=operator.attrgetter('pos_bnd1'))

    def get_bnd2adv_form(self):
        bnds_equivs = self.get_equivalents()
        if bnds_equivs[0].is5prime_bnd2:
            return min(bnds_equivs, key=(lambda bnds: bnds.pos_bnd2))
        else:
            return max(bnds_equivs, key=(lambda bnds: bnds.pos_bnd2))

    ###########################################

    def get_params(self, is_bnd1):
        if is_bnd1:
            return (self.chrom_bnd1, 
                    self.get_pos_range0_bnd1(),
                    self.is5prime_bnd1)
        else:
            return (self.chrom_bnd2,
                    self.get_pos_range0_bnd2(),
                    self.is5prime_bnd2)

    def get_vcfspec_bnd1(self):
        chrom = self.chrom_bnd1
        pos = self.pos_bnd1
        ref = self.fasta.fetch(self.chrom_bnd1, self.pos_bnd1 - 1, self.pos_bnd1)

        if self.is5prime_bnd1:
            t = ''.join(self.inserted_seq) + ref
        else:
            t = ref + ''.join(self.inserted_seq)
        bracket = '[' if self.is5prime_bnd2 else ']'
        alt_matestring = f'{bracket}{self.chrom_bnd2}:{self.pos_bnd2}{bracket}'
        if self.is5prime_bnd1:
            alt = alt_matestring + t
        else:
            alt = t + alt_matestring

        return Vcfspec(chrom, pos, ref, [alt], refver=self.refver)

    def get_vcfspec_bnd2(self):
        chrom = self.chrom_bnd2
        pos = self.pos_bnd2
        ref = self.fasta.fetch(self.chrom_bnd2, self.pos_bnd2 - 1, 
                               self.pos_bnd2)

        insseq_converted = self._convert_seq_between_bnds(
            ''.join(self.inserted_seq))
        if self.is5prime_bnd2:
            t = insseq_converted + ref
        else:
            t = ref + insseq_converted

        bracket = '[' if self.is5prime_bnd1 else ']'
        alt_matestring = f'{bracket}{self.chrom_bnd1}:{self.pos_bnd1}{bracket}'
        if self.is5prime_bnd2:
            alt = alt_matestring + t
        else:
            alt = t + alt_matestring

        return Vcfspec(chrom, pos, ref, [alt], refver=self.refver)

    def get_simplesv(self):
        if self.chrom_bnd1 != self.chrom_bnd2:
            raise Exception(f'Translocation cannot be converted to a '
                            f'SimpleStructuralVariant object.')
        else:
            chrom = self.chrom_bnd1
            if self.svtype == 'DEL':
                start1 = self.pos_bnd1 + 1
                end1 = self.pos_bnd2 - 1
                simplesv = structvars.Deletion(chrom=chrom, start1=start1, 
                                               end1=end1, fasta=self.fasta)
            elif self.svtype == 'INV':
                if self.is5prime_bnd1:
                    start1 = self.pos_bnd1
                    end1 = self.pos_bnd2 - 1
                else:
                    start1 = self.pos_bnd1 + 1
                    end1 = self.pos_bnd2
                simplesv = structvars.Inversion(chrom=chrom, start1=start1, 
                                                end1=end1, fasta=self.fasta)
            elif self.svtype == 'DUP':
                start1 = self.pos_bnd1
                end1 = self.pos_bnd2
                simplesv = structvars.TandemDuplication(
                    chrom=chrom, start1=start1, end1=end1, fasta=self.fasta)

            return simplesv

    def get_hgvsg(self):
        return self.get_simplesv().get_hgvsg()

    def get_id(self):
        return '_'.join(
            [
                self.chrom_bnd1,
                str(self.pos_bnd1),
                self.chrom_bnd2,
                str(self.pos_bnd2),
            ]
        )

    def get_id_bnd1(self):
        return self.get_id() + '_1'

    def get_id_bnd2(self):
        return self.get_id() + '_2'

    ###########################################

    @functools.cache
    def get_border_seq(self, is_bnd1):
        (chrom, pos_range0, is5prime) = self.get_params(is_bnd1=is_bnd1)
        fetch_start = min(pos_range0)
        fetch_end = max(pos_range0) + 1
        return self.fasta.fetch(chrom, fetch_start, fetch_end)

    @functools.cache
    def get_seq_beyond_bnd_ref(self, is_bnd1, length, with_border_seq):
        (chrom, pos_range0, is5prime) = self.get_params(is_bnd1=is_bnd1)
        if is5prime:
            fetch_end0 = min(pos_range0)
            fetch_start0 = fetch_end0 - length
        else:
            fetch_start0 = pos_range0.stop
            fetch_end0 = fetch_start0 + length

        seq_beyond = self.fasta.fetch(chrom, fetch_start0, fetch_end0)

        if with_border_seq:
            border_seq = self.get_border_seq(is_bnd1)
            if is5prime:
                return seq_beyond + border_seq
            else:
                return border_seq + seq_beyond
        else:
            return seq_beyond

    @functools.cache
    def get_seq_beyond_bnd_alt(self, is_bnd1, length, with_border_seq):
        # get parameters
        is5prime = self.is5prime_bnd1 if is_bnd1 else self.is5prime_bnd2
        (chrom_mate, pos_range0_mate, is5prime_mate) = self.get_params(
            is_bnd1=(not is_bnd1))
        fetch_length = length - len(self.inserted_seq)

        # fetch_start, fetch_end
        if is5prime_mate:
            fetch_start = pos_range0_mate.start
            fetch_end = fetch_start + fetch_length
        else:
            fetch_end = pos_range0_mate.start + 1
            fetch_start = fetch_end - fetch_length

        # mateside_seq
        if is5prime == is5prime_mate:
            mateside_seq = Bio.Seq.reverse_complement(
                self.fasta.fetch(chrom_mate, fetch_start, fetch_end))
        else:
            mateside_seq = self.fasta.fetch(chrom_mate, fetch_start,
                                            fetch_end)

        # inserted_seq
        if is_bnd1:
            inserted_seq = ''.join(self.inserted_seq)
        else:
            if is5prime == is5prime_mate:
                inserted_seq = Bio.Seq.reverse_complement(''.join(self.inserted_seq))
            else:
                inserted_seq = ''.join(self.inserted_seq)

        # result
        if is5prime:
            seq_beyond = mateside_seq + inserted_seq
        else:
            seq_beyond = inserted_seq + mateside_seq

        # add border_seq
        if with_border_seq:
            border_seq = self.get_border_seq(is_bnd1)
            if is5prime:
                return seq_beyond + border_seq
            else:
                return border_seq + seq_beyond
        else:
            return seq_beyond

        return result

    ###########################################

#    def set_equivs(self):
#        self._equivs = get_bnds_equivalents(self)
#
#    def get_equivs(self):
#        if self._equivs is None:
#            self.set_equivs()
#        return self._equivs
#
#    def get_sorted_equivs(self, mode):
#        return sort_equivs(self.get_equivs(), mode)

    ###############

    @functools.cache
    #def set_pos_range0s(self):
    def get_pos_range0s(self):
        poslist0_bnd1 = list()
        poslist0_bnd2 = list()
        for bnds in self.get_equivs():
            poslist0_bnd1.append(bnds.pos0_bnd1)
            poslist0_bnd2.append(bnds.pos0_bnd2)

        if self.is5prime_bnd1:
            #self._pos_range0_bnd1 = range(max(poslist0_bnd1),
            pos_range0_bnd1 = range(
                max(poslist0_bnd1), min(poslist0_bnd1) - 1, -1,
            )
        else:
            #self._pos_range0_bnd1 = range(min(poslist0_bnd1),
            pos_range0_bnd1 = range(
                min(poslist0_bnd1), max(poslist0_bnd1) + 1, 1,
            )

        if self.is5prime_bnd2:
            #self._pos_range0_bnd2 = range(max(poslist0_bnd2),
            pos_range0_bnd2 = range(
                max(poslist0_bnd2),
                min(poslist0_bnd2) - 1,
                -1,
            )
        else:
            #self._pos_range0_bnd2 = range(min(poslist0_bnd2),
            pos_range0_bnd2 = range(
                min(poslist0_bnd2),
                max(poslist0_bnd2) + 1,
                1,
            )

        return pos_range0_bnd1, pos_range0_bnd2

    def get_pos_range0_bnd1(self):
        """Returns a directional range"""
        #if self._pos_range0_bnd1 is None:
        #    self.set_pos_range0s()
        #return self._pos_range0_bnd1
        return self.get_pos_range0s()[0]

    def get_pos_range0_bnd2(self):
        """Returns a directional range"""
        #if self._pos_range0_bnd2 is None:
        #    self.set_pos_range0s()
        #return self._pos_range0_bnd2
        return self.get_pos_range0s()[1]

    @functools.cache
    def get_flank_range0(self, mode, is_bnd1, flanklen):
        """Args:
            mode: Must be one of "par", "bnd_prox", or "bnd_dist".
                par: partner side
                bnd_prox: breakend side, proximal (adjacent to the most 
                    advanced border)
                bnd_dist: breakend side, distal (adjacent to the most 
                    retracted border)

                             BND side | PARTNER side
                      |               |
        --------------|---|---|---|---|
                      |   |   |   |   | <====>
        --------------|---|---|---|---|  "par" flank 
                      |  borderzone   |
               <=====>         <=====>
              bnd_dist        bnd_prox
                 flank           flank
        """
        assert mode in ('par', 'bnd_prox', 'bnd_dist')

        # modify flanklen for breakend-side flanks
        #if mode in ('bnd_prox', 'bnd_dist'):
        #    flanklen = flanklen - 1

        # get parameters
        pos_range0 = (
            self.get_pos_range0_bnd1()
            if is_bnd1 else
            self.get_pos_range0_bnd2()
        )

        # main
        step = pos_range0.step

        if mode == 'par':
            start = pos_range0.stop
            stop = start + (step * flanklen)
        else:
            if mode == 'bnd_prox':
                stop = pos_range0.stop - step
            elif mode == 'bnd_dist':
                stop = pos_range0.start

            start = stop + (-step * flanklen)

        return range(start, stop, step)

#    @functools.cache
#    def get_parflank_range0_bnd1(self, flanklen=1):
#        return self._get_flankrange_helper(mode='par', is_bnd1=True,
#                                           flanklen=flanklen)
#
#    @functools.cache
#    def get_parflank_range0_bnd2(self, flanklen=1):
#        return self._get_flankrange_helper(mode='par', is_bnd1=False,
#                                           flanklen=flanklen)
#
#    @functools.cache
#    def get_bndproxflank_range0_bnd1(self, flanklen=1):
#        """Returns a directional range"""
#
#        return self._get_flankrange_helper(mode='bnd_prox', is_bnd1=True,
#                                           flanklen=(flanklen-1))
#
#    @functools.cache
#    def get_bndsideflank_range0_bnd2(self, flanklen=1):
#        """Returns a directional range"""
#
#        return self._get_flankrange_helper(mode=False, is_bnd1=False,
#                                           flanklen=(flanklen-1))

    ###################################################################

    #def get_bnd1adv_form(self):
    #    return pick_bnd1adv_form(self.get_equivs())

    #def get_bnd2adv_form(self):
    #    return pick_bnd2adv_form(self.get_equivs())
    
    ###################################################################

    def spawn(self, **kwargs):
        kwargs = (
            dict(
                bndspec1=self.bndspec1,
                bndspec2=self.bndspec2,
                refver=self.refver,
                inserted_seq=self.inserted_seq,
                score=self.score,
            ) | kwargs
        )
        return Breakends(**kwargs)

#        return Breakends(
#            chrom_bnd1=self.chrom_bnd1, 
#            pos_bnd1=self.pos_bnd1, 
#            is5prime_bnd1=self.is5prime_bnd1,
#            chrom_bnd2=self.chrom_bnd2, 
#            pos_bnd2=self.pos_bnd2, 
#            is5prime_bnd2=self.is5prime_bnd2,
#            fasta=self.fasta, chromdict=self.chromdict,
#            inserted_seq=self.inserted_seq.copy(), svtype=self.svtype,
#            score=self.score)

    def identical(self, other):
        return hash(self) == hash(other)

    def sameseq1(self, other):
        return self.get_bnd1adv_form() == other.get_bnd1adv_form()

    def sameseq2(self, other):
        if any(
            (getattr(self, key) != getattr(other, key))
            for key in ('chrom_bnd1', 'is5prime_bnd1', 'chrom_bnd2', 'is5prime_bnd2')
        ):
            return False

        if self.is5prime_bnd1:
            goal_bnd1 = max(self.pos_bnd1, other.pos_bnd1)
        else:
            goal_bnd1 = min(self.pos_bnd1, other.pos_bnd1)

        if self.is5prime_bnd2:
            goal_bnd2 = max(self.pos_bnd2, other.pos_bnd2)
        else:
            goal_bnd2 = min(self.pos_bnd2, other.pos_bnd2)

        def helper(bnds):
            result = bnds

            while True:
                if result.pos_bnd1 == goal_bnd1:
                    break
                else:
                    result = result.retract_bnd1()
                    #bnds.retract_bnd1()
                    continue

            while True:
                if result.pos_bnd2 == goal_bnd2:
                    break
                else:
                    result = result.retract_bnd2()
                    #bnds.retract_bnd2()
                    continue

            return result

        modified_self = helper(self)
        modified_other = helper(other)

        return modified_self.equal(modified_other)

    ###################################################################

    def advance_bnd1(self):
        def get_newpos(self):
            newpos_bnd1 = self._get_advanced_pos_bnd1()
            if len(self.inserted_seq) > 0:
                newpos_bnd2 = self.pos_bnd2
            else:
                newpos_bnd2 = self._get_retracted_pos_bnd2()

            return newpos_bnd1, newpos_bnd2

        def base_match_check(self, newpos_bnd1):
            newbase = self.fasta.fetch(self.chrom_bnd1, newpos_bnd1 - 1, newpos_bnd1)
            if len(self.inserted_seq) > 0:
                oldbase = (self.inserted_seq[-1] 
                           if self.is5prime_bnd1 else 
                           self.inserted_seq[0])
            else:
                oldbase = self._convert_seq_between_bnds(
                    self.fasta.fetch(self.chrom_bnd2, self.pos_bnd2 - 1, self.pos_bnd2)
                )

            return (newbase == oldbase)

        def action(self, newpos_bnd1, newpos_bnd2):
            new_bndspec1 = self.bndspec1.spawn(pos0=(newpos_bnd1 - 1))
            new_bndspec2 = self.bndspec2.spawn(pos0=(newpos_bnd2 - 1))
            new_score = self.score
            if len(self.inserted_seq) > 0:
                new_score += 1
                if self.is5prime_bnd1:
                    new_inserted_seq = self.inserted_seq[:-1]
                else:
                    new_inserted_seq = self.inserted_seq[1:]
            else:
                new_inserted_seq = self.inserted_seq

            #self.pos_bnd1 = newpos_bnd1
            #self.pos_bnd2 = newpos_bnd2
            #if len(self.inserted_seq) > 0:
            #    self.score += 1
            #    if self.is5prime_bnd1:
            #        del self.inserted_seq[-1]
            #    else:
            #        del self.inserted_seq[0]

            return self.spawn(
                bndspec1=new_bndspec1,
                bndspec2=new_bndspec2,
                inserted_seq=new_inserted_seq,
                score=new_score,
            )

        # main
        newpos_bnd1, newpos_bnd2 = get_newpos(self)
        if self._newpos_range_check(newpos_bnd1, newpos_bnd2):
            if base_match_check(self, newpos_bnd1):
                result = action(self, newpos_bnd1, newpos_bnd2)
                success = True
            else:
                result = self.spawn()
                success = False
        else:
            result = self.spawn()
            success = False

        return result, success

    def advance_bnd2(self):
        def get_newpos(self):
            newpos_bnd2 = self._get_advanced_pos_bnd2()
            if len(self.inserted_seq) > 0:
                newpos_bnd1 = self.pos_bnd1
            else:
                newpos_bnd1 = self._get_retracted_pos_bnd1()

            return newpos_bnd1, newpos_bnd2

        def base_match_check(self, newpos_bnd2):
            newbase = self._convert_seq_between_bnds(
                self.fasta.fetch(self.chrom_bnd2, newpos_bnd2 - 1, newpos_bnd2)
            )

            if len(self.inserted_seq) > 0:
                oldbase = (self.inserted_seq[0] 
                           if self.is5prime_bnd1 else 
                           self.inserted_seq[-1])
            else:
                oldbase = self.fasta.fetch(self.chrom_bnd1, self.pos_bnd1 - 1, self.pos_bnd1)

            return (newbase == oldbase)

        def action(self, newpos_bnd1, newpos_bnd2):
            new_score = self.score
            new_bndspec1 = self.bndspec1.spawn(pos0=(newpos_bnd1 - 1))
            new_bndspec2 = self.bndspec2.spawn(pos0=(newpos_bnd2 - 1))
            if len(self.inserted_seq) > 0:
                new_score += 1
                if self.is5prime_bnd1:
                    new_inserted_seq = self.inserted_seq[1:]
                else:
                    new_inserted_seq = self.inserted_seq[:-1]
            else:
                new_inserted_seq = self.inserted_seq

#            self.pos_bnd1 = newpos_bnd1
#            self.pos_bnd2 = newpos_bnd2
#            if len(self.inserted_seq) > 0:
#                self.score += 1
#                if self.is5prime_bnd1:
#                    del self.inserted_seq[0]
#                else:
#                    del self.inserted_seq[-1]

            return self.spawn(
                bndspec1=new_bndspec1,
                bndspec2=new_bndspec2,
                inserted_seq=new_inserted_seq,
                score=new_score,
            )

        # main
        newpos_bnd1, newpos_bnd2 = get_newpos(self)
        if self._newpos_range_check(newpos_bnd1, newpos_bnd2):
            if base_match_check(self, newpos_bnd2):
                result = action(self, newpos_bnd1, newpos_bnd2)
                success = True
            else:
                result = self.spawn()
                success = False
        else:
            result = self.spawn()
            success = False

        return result, success

    def retract_bnd1(self):
        added_base = self.fasta.fetch(
            self.chrom_bnd1, self.pos_bnd1 - 1, self.pos_bnd1,
        )
        if self.is5prime_bnd1:
            new_inserted_seq = self.inserted_seq + (added_base,)
            #self.inserted_seq.append(added_base)
        else:
            new_inserted_seq = (added_base,) + self.inserted_seq 
            #self.inserted_seq.insert(0, added_base)

        new_bndspec1 = self.bndspec1.spawn(
            pos0=(self._get_retracted_pos_bnd1() - 1),
        )
        new_score = self.score - 1
        #self.pos_bnd1 = self._get_retracted_pos_bnd1()
        #self.score -= 1

        return self.spawn(
            inserted_seq=new_inserted_seq,
            bndspec1=new_bndspec1,
            score=new_score,
        )

    def retract_bnd2(self):
        added_base = self._convert_seq_between_bnds(
            self.fasta.fetch(self.chrom_bnd2, self.pos_bnd2 - 1, self.pos_bnd2)
        )
        if self.is5prime_bnd1:
            new_inserted_seq = (added_base,) + self.inserted_seq
            #self.inserted_seq.insert(0, added_base)
        else:
            new_inserted_seq = self.inserted_seq + (added_base,)
            #self.inserted_seq.append(added_base)

        new_bndspec2 = self.bndspec2.spawn(
            pos0=(self._get_retracted_pos_bnd2() - 1),
        )
        new_score = self.score - 1
        #self.pos_bnd2 = self._get_retracted_pos_bnd2()
        #self.score -= 1

        return self.spawn(
            inserted_seq=new_inserted_seq,
            bndspec2=new_bndspec2,
            score=new_score,
        )

    #######################################################

    def _set_svtype(self):
        if self.chrom_bnd1 != self.chrom_bnd2:
            self.svtype = 'TRA'
        else:
            if (not self.is5prime_bnd1) and self.is5prime_bnd2:
                self.svtype = 'DEL'
            elif self.is5prime_bnd1 and (not self.is5prime_bnd2):
                self.svtype = 'DUP'
            elif self.is5prime_bnd1 == self.is5prime_bnd2:
                self.svtype = 'INV'

    def _newpos_range_check(self, newpos_bnd1, newpos_bnd2):
        return (newpos_bnd1 >= 1 and 
                newpos_bnd1 <= self.chromdict[self.chrom_bnd1] and
                newpos_bnd2 >= 1 and 
                newpos_bnd2 <= self.chromdict[self.chrom_bnd2])

    def _get_advanced_pos_bnd1(self):
        if self.is5prime_bnd1:
            return self.pos_bnd1 - 1
        else:
            return self.pos_bnd1 + 1

    def _get_retracted_pos_bnd1(self):
        if self.is5prime_bnd1:
            return self.pos_bnd1 + 1
        else:
            return self.pos_bnd1 - 1

    def _get_advanced_pos_bnd2(self):
        if self.is5prime_bnd2:
            return self.pos_bnd2 - 1
        else:
            return self.pos_bnd2 + 1

    def _get_retracted_pos_bnd2(self):
        if self.is5prime_bnd2:
            return self.pos_bnd2 + 1
        else:
            return self.pos_bnd2 - 1

    def _convert_seq_between_bnds(self, seq):
        return convert_seq_between_bnds(seq, 
                                        self.is5prime_bnd1, 
                                        self.is5prime_bnd2)


#######################################################


class NonStandardSVRecord(Exception):
    pass


def convert_seq_between_bnds(seq, is5prime_bnd1, is5prime_bnd2):
    if is5prime_bnd1 == is5prime_bnd2:
        return Bio.Seq.reverse_complement(seq)
    else:
        return seq


#######################################################


#def get_bnds_equivalents(bnds):
#    """Returns: A list of Breakends objects, all equivalent to the 
#        input object, sorted such that the first item is the most advanced 
#        form with respect to bnd1, and the last item is the most advanced 
#        form with respect to bnd2.
#    """
#
#    input_copy = bnds.spawn()
#
#    bnds_list_bnd1adv = list()
#    bnds_list_bnd1adv.append(input_copy)
#    while True:
#        new = bnds_list_bnd1adv[-1].spawn()
#        success = new.advance_bnd1()
#        if success:
#            bnds_list_bnd1adv.append(new)
#            continue
#        else:
#            break
#
#    bnds_list_bnd2adv = list()
#    bnds_list_bnd2adv.append(input_copy)
#    while True:
#        new = bnds_list_bnd2adv[-1].spawn()
#        success = new.advance_bnd2()
#        if success:
#            bnds_list_bnd2adv.append(new)
#            continue
#        else:
#            break
#
#    bnds_list_bnd1adv.reverse()
#    bnds_list = bnds_list_bnd1adv[:-1] + [input_copy] + bnds_list_bnd2adv[1:]
#    max_score = max(x.score for x in bnds_list)
#
#    bnds_equivs = [x for x in bnds_list if x.score == max_score]
#
#    return bnds_equivs


#def sort_equivs(bnds_equivs, mode):
#    assert mode in ('bnd1adv_first', 'bnd1adv_last', 'bnd2adv_first', 
#                    'bnd2adv_last'), 'Invalid "mode" argument'
#
#    if mode == 'bnd1adv_first':
#        reverse = (not bnds_equivs[0].is5prime_bnd1)
#        key = lambda bnds: bnds.pos_bnd1
#    elif mode == 'bnd1adv_last':
#        reverse = bnds_equivs[0].is5prime_bnd1
#        key = lambda bnds: bnds.pos_bnd1
#    elif mode == 'bnd2adv_first':
#        reverse = (not bnds_equivs[0].is5prime_bnd2)
#        key = lambda bnds: bnds.pos_bnd2
#    elif mode == 'bnd2adv_last':
#        reverse = bnds_equivs[0].is5prime_bnd2
#        key = lambda bnds: bnds.pos_bnd2
#
#    return sorted(bnds_equivs, key=key, reverse=reverse)


#def pick_bnd1adv_form(bnds_equivs):
#    if bnds_equivs[0].is5prime_bnd1:
#        return min(bnds_equivs, key=(lambda bnds: bnds.pos_bnd1))
#    else:
#        return max(bnds_equivs, key=(lambda bnds: bnds.pos_bnd1))
#
#
#def pick_bnd2adv_form(bnds_equivs):
#    if bnds_equivs[0].is5prime_bnd2:
#        return min(bnds_equivs, key=(lambda bnds: bnds.pos_bnd2))
#    else:
#        return max(bnds_equivs, key=(lambda bnds: bnds.pos_bnd2))


#######################################################

#def get_bnds_from_vr(vr, fasta, chromdict):
#    """
#    Raises:
#        If input vr does not conform to a known SV variant record format 
#            (including a valid non-SV variant record)
#    """
#
#    assert len(vr.alts) == 1, (
#        f'Multiallelic variant record is not allowed:\n{vr}')
#
#    vr_svinfo = get_vr_svinfo_standard_vr(vr, fasta, chromdict)
#    bnds = get_bnds_from_vr_svinfo(vr, vr_svinfo, fasta, chromdict)
#
#    return bnds


########################################################


def get_bnds_from_vcfspec_parseinfo(
    chrom,
    pos,
    ref,
    t,
    chrom_mate,
    pos_mate,
    current_is5prime,
    mate_is5prime,
    is_bnd1,
    refver,
):
    if is_bnd1:
        chrom_bnd1 = chrom
        chrom_bnd2 = chrom_mate
        is5prime_bnd1 = current_is5prime
        is5prime_bnd2 = mate_is5prime

        if is5prime_bnd1:
            if t[-1] == ref:
                pos_bnd1 = pos
                inserted_seq = list( t[:-1] )
            else:
                warn()
                pos_bnd1 = pos + 1
                inserted_seq = list( t )
        else:
            if t[0] == ref:
                pos_bnd1 = pos
                inserted_seq = list( t[1:] )
            else:
                warn()
                pos_bnd1 = pos - 1
                inserted_seq = list( t )
            
        pos_bnd2 = pos_mate

    else:
        chrom_bnd1 = chrom_mate
        chrom_bnd2 = chrom
        is5prime_bnd1 = mate_is5prime
        is5prime_bnd2 = current_is5prime

        pos_bnd1 = pos_mate

        if is5prime_bnd2:
            if t[-1] == ref:
                pos_bnd2 = pos
                inserted_seq = list(
                    convert_seq_between_bnds(t[:-1], 
                                             is5prime_bnd1, is5prime_bnd2))
            else:
                warn()
                pos_bnd2 = pos + 1
                inserted_seq = list(
                    convert_seq_between_bnds(t, 
                                             is5prime_bnd1, is5prime_bnd2))
        else:
            if t[0] == ref:
                pos_bnd2 = pos
                inserted_seq = list(
                    convert_seq_between_bnds(t[1:], 
                                             is5prime_bnd1, is5prime_bnd2))
            else:
                warn()
                pos_bnd2 = pos - 1
                inserted_seq = list(
                    convert_seq_between_bnds(t, 
                                             is5prime_bnd1, is5prime_bnd2))

    bnds = Breakends(
        bndspec1=(chrom_bnd1, pos_bnd1 - 1, is5prime_bnd1),
        bndspec2=(chrom_bnd2, pos_bnd2 - 1, is5prime_bnd2),
        refver=refver,
        inserted_seq=inserted_seq,
    )

    return bnds


def get_bnds_from_vr_svinfo(vr, vr_svinfo, fasta, chromdict):
    def warn():
        LOGGER.warning(f'"t" portion of SV ALT string is not an extension '
                       f'of REF string for this variant record:\n{vr}')

    if vr_svinfo['is_bnd1']:
        chrom_bnd1 = vr.contig
        chrom_bnd2 = vr_svinfo['chrom_mate']
        is5prime_bnd1 = vr_svinfo['current_is5prime']
        is5prime_bnd2 = vr_svinfo['mate_is5prime']

        if is5prime_bnd1:
            if vr_svinfo['t'][-1] == vr_svinfo['ref']:
                pos_bnd1 = vr.pos
                inserted_seq = list( vr_svinfo['t'][:-1] )
            else:
                warn()
                pos_bnd1 = vr.pos + 1
                inserted_seq = list( vr_svinfo['t'] )
        else:
            if vr_svinfo['t'][0] == vr_svinfo['ref']:
                pos_bnd1 = vr.pos
                inserted_seq = list( vr_svinfo['t'][1:] )
            else:
                warn()
                pos_bnd1 = vr.pos - 1
                inserted_seq = list( vr_svinfo['t'] )
            
        pos_bnd2 = vr_svinfo['pos_mate']

    else:
        chrom_bnd1 = vr_svinfo['chrom_mate']
        chrom_bnd2 = vr.contig
        is5prime_bnd1 = vr_svinfo['mate_is5prime']
        is5prime_bnd2 = vr_svinfo['current_is5prime']

        pos_bnd1 = vr_svinfo['pos_mate']

        if is5prime_bnd2:
            if vr_svinfo['t'][-1] == vr_svinfo['ref']:
                pos_bnd2 = vr.pos
                inserted_seq = list(
                    convert_seq_between_bnds(vr_svinfo['t'][:-1], 
                                             is5prime_bnd1, is5prime_bnd2))
            else:
                warn()
                pos_bnd2 = vr.pos + 1
                inserted_seq = list(
                    convert_seq_between_bnds(vr_svinfo['t'], 
                                             is5prime_bnd1, is5prime_bnd2))
        else:
            if vr_svinfo['t'][0] == vr_svinfo['ref']:
                pos_bnd2 = vr.pos
                inserted_seq = list(
                    convert_seq_between_bnds(vr_svinfo['t'][1:], 
                                             is5prime_bnd1, is5prime_bnd2))
            else:
                warn()
                pos_bnd2 = vr.pos - 1
                inserted_seq = list(
                    convert_seq_between_bnds(vr_svinfo['t'], 
                                             is5prime_bnd1, is5prime_bnd2))

    bnds = Breakends(
        bndspec1=(chrom_bnd1, pos_bnd1 - 1, is5prime_bnd1),
        bndspec2=(chrom_bnd2, pos_bnd2 - 1, is5prime_bnd2),
        refver=refgenome.infer_refver_vcfheader(vr.header),
        inserted_seq=inserted_seq,
    )
#    bnds = Breakends(chrom_bnd1=chrom_bnd1, pos_bnd1=pos_bnd1,
#                     is5prime_bnd1=is5prime_bnd1,
#                     chrom_bnd2=chrom_bnd2, pos_bnd2=pos_bnd2,
#                     is5prime_bnd2=is5prime_bnd2,
#                     inserted_seq=inserted_seq, fasta=fasta,
#                     chromdict=chromdict)

    return bnds


def get_vr_svinfo_standard_vr(vr, fasta, chromdict):
    vr_svinfo = dict()
    vr_svinfo['ref'] = vr.ref

    try:
        (vr_svinfo['t'],
         vr_svinfo['chrom_mate'],
         vr_svinfo['pos_mate'],
         vr_svinfo['current_is5prime'],
         vr_svinfo['mate_is5prime']) = libvcfspec.parse_sv_altstring(vr.alts[0])
    except NonStandardSVAlt as e:
        e_msg = f'{str(e)}\nInput variant record:\n{vr}'
        raise NonStandardSVRecord(e_msg)

    vr_svinfo['is_bnd1'] = get_is_bnd1_withvr(vr, vr_svinfo, chromdict)

    return vr_svinfo




def get_is_bnd1_withvr(vr, vr_svinfo, chromdict):
    order = tools.compare_coords(
        vr.contig, vr.pos, vr_svinfo['chrom_mate'], vr_svinfo['pos_mate'], chromdict,
    )
    if order < 0:
        is_bnd1 = True
    elif order > 0:
        is_bnd1 = False
    elif order == 0:
        is_bnd1 = True
#        raise Exception(
#            f'Current chrom/pos and mate chrom/pos are identical '
#            f'for vcf record :\n{vr}')

    return is_bnd1


def get_is_bnd1(vcfspec, chrom_mate, pos_mate, chromdict):
    return libvcfspec.get_is_bnd1_base(vcfspec.chrom, vcfspec.pos, chrom_mate, pos_mate, chromdict)



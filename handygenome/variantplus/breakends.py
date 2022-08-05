import itertools
import logging
import functools

import Bio.Seq

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
structvars = importlib.import_module('.'.join([top_package_name, 'svlib', 'structvars']))


LOGGER = workflow.get_logger(
    name=__name__,
    formatter=logging.Formatter(
        fmt='[%(asctime)s  %(levelname)s] %(name)s - %(message)s',
        datefmt=workflow.DEFAULT_DATEFMT),
    level='info', stderr=True, filenames=None, append=False)


class Breakends:
    """
    Attributes:
        chrom_bnd1
        pos_bnd1
        endis5_bnd1
        chrom_bnd2
        pos_bnd2
        endis5_bnd2
        inserted_seq: 
            - A list composed of the bases of the inserted sequence 
                (e.g. ['A', 'A', 'G', 'C'])
            - As seen in the viewpoint where bnd1-side sequence is on 
                the plus(Crick) strand

        fasta: pysam.FastaFile instance
        chromdict: julib.common.ChromDict instance
        svtype
        score
        endis3_bnd1
        endis3_bnd2

        _equivs: A list of Breakends objects with maximum reference coverage.
        _pos_range0_bnd1
        _pos_range0_bnd2
        _homlen = None
        _homseq = None
    """

    BASIC_ATTRS = ('chrom_bnd1', 'pos_bnd1', 'endis5_bnd1', 
                   'chrom_bnd2', 'pos_bnd2', 'endis5_bnd2', 
                   'inserted_seq', 'svtype')

    def __init__(self, chrom_bnd1, pos_bnd1, endis5_bnd1, 
                 chrom_bnd2, pos_bnd2, endis5_bnd2, fasta, 
                 inserted_seq=list(), chromdict=None, 
                 svtype=None, score=None):
        # basic attributes
        self.chrom_bnd1 = chrom_bnd1
        self.pos_bnd1 = pos_bnd1
        self.endis5_bnd1 = endis5_bnd1
        self.chrom_bnd2 = chrom_bnd2
        self.pos_bnd2 = pos_bnd2
        self.endis5_bnd2 = endis5_bnd2
        self.inserted_seq = inserted_seq
        self.fasta = fasta

        if chromdict is None:
            self.chromdict = common.ChromDict(fasta=fasta)
        else:
            self.chromdict = chromdict

        # endis3
        self.endis3_bnd1 = (not self.endis5_bnd1)
        self.endis3_bnd2 = (not self.endis5_bnd2)

        # svtype
        if svtype is None:
            self._set_svtype()
        else:
            self.svtype = svtype

        # score
        if score is None:
            self.score = 0
        else:
            self.score = score

        # private attributes accesible with getters
        self._equivs = None
        self._homlen = None
        self._homseq = None
        self._pos_range0_bnd1 = None
        self._pos_range0_bnd2 = None

    def __repr__(self):
        def func1(self):
            result = ['<Breakends>']
            result.extend(
                [f'{key}: {getattr(self, key)}' 
                 for key in ('chrom_bnd1', 'pos_bnd1', 'endis5_bnd1', 
                             'inserted_seq',
                             'chrom_bnd2', 'pos_bnd2', 'endis5_bnd2')])
            return '\n'.join(result)

        def func2(self):
            brkt1 = '[' if self.endis5_bnd1 else ']'
            brkt2 = '[' if self.endis5_bnd2 else ']'
            str1 = f'{brkt1}{self.chrom_bnd1}:{self.pos_bnd1}{brkt1}'
            str2 = f'{brkt2}{self.chrom_bnd2}:{self.pos_bnd2}{brkt2}'
            insseq = ''.join(self.inserted_seq)

            return f'<Breakends {str1} {insseq} {str2} (1-based coords)>'

        def func3(self):
            endtype1 = '5\'' if self.endis5_bnd1 else '3\''
            endtype2 = '5\'' if self.endis5_bnd2 else '3\''
            if len(self.inserted_seq) == 0:
                insseq = 'no inserted seq'
            else:
                insseq = ''.join(self.inserted_seq)

            return (f'<Breakends '
                    f'{self.chrom_bnd1}:{self.pos_bnd1} {endtype1} '
                    f'[{insseq}] '
                    f'{self.chrom_bnd2}:{self.pos_bnd2} {endtype2}>')

        return func2(self)

    def __hash__(self):
        return hash((self.chrom_bnd1, self.pos_bnd1, self.endis5_bnd1,
                     self.chrom_bnd2, self.pos_bnd2, self.endis5_bnd2,
                     tuple(self.inserted_seq)))

#    def __hash__(self):
#        bnd1adv_form = get_bnd1adv_form(self)
#        return hash((
#                    bnd1adv_form.chrom_bnd1,
#                    bnd1adv_form.pos_bnd1,
#                    bnd1adv_form.endis5_bnd1,
#                    bnd1adv_form.chrom_bnd2,
#                    bnd1adv_form.pos_bnd2,
#                    bnd1adv_form.endis5_bnd2,
#                    bnd1adv_form.inserted_seq,
#                    ))

    def __eq__(self, other):
        #assert isinstance(other, Breakends)
        return hash(self) == hash(other)

    ###########################################

    def get_params(self, is_bnd1):
        if is_bnd1:
            return (self.chrom_bnd1, 
                    self.get_pos_range0_bnd1(),
                    self.endis5_bnd1)
        else:
            return (self.chrom_bnd2,
                    self.get_pos_range0_bnd2(),
                    self.endis5_bnd2)

    def get_vcfspec_bnd1(self):
        chrom = self.chrom_bnd1
        pos = self.pos_bnd1
        ref = self.fasta.fetch(self.chrom_bnd1, self.pos_bnd1 - 1, 
                               self.pos_bnd1)

        if self.endis5_bnd1:
            t = ''.join(self.inserted_seq) + ref
        else:
            t = ref + ''.join(self.inserted_seq)
        bracket = '[' if self.endis5_bnd2 else ']'
        alt_matestring = f'{bracket}{self.chrom_bnd2}:{self.pos_bnd2}{bracket}'
        if self.endis5_bnd1:
            alt = alt_matestring + t
        else:
            alt = t + alt_matestring

        return common.Vcfspec(chrom, pos, ref, [alt])

    def get_vcfspec_bnd2(self):
        chrom = self.chrom_bnd2
        pos = self.pos_bnd2
        ref = self.fasta.fetch(self.chrom_bnd2, self.pos_bnd2 - 1, 
                               self.pos_bnd2)

        insseq_converted = self._convert_seq_between_bnds(
            ''.join(self.inserted_seq))
        if self.endis5_bnd2:
            t = insseq_converted + ref
        else:
            t = ref + insseq_converted

        bracket = '[' if self.endis5_bnd1 else ']'
        alt_matestring = f'{bracket}{self.chrom_bnd1}:{self.pos_bnd1}{bracket}'
        if self.endis5_bnd2:
            alt = alt_matestring + t
        else:
            alt = t + alt_matestring

        return common.Vcfspec(chrom, pos, ref, [alt])

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
                if self.endis5_bnd1:
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
        return '_'.join([self.chrom_bnd1,
                         str(self.pos_bnd1),
                         self.chrom_bnd2,
                         str(self.pos_bnd2)])

    def get_id_bnd1(self):
        return self.get_id() + '_1'

    def get_id_bnd2(self):
        return self.get_id() + '_2'

    ###########################################

    @functools.cache
    def get_border_seq(self, is_bnd1):
        (chrom, pos_range0, endis5) = self.get_params(is_bnd1=is_bnd1)
        fetch_start = min(pos_range0)
        fetch_end = max(pos_range0) + 1
        return self.fasta.fetch(chrom, fetch_start, fetch_end)

    @functools.cache
    def get_seq_beyond_bnd_ref(self, is_bnd1, length, with_border_seq):
        (chrom, pos_range0, endis5) = self.get_params(is_bnd1=is_bnd1)
        if endis5:
            fetch_end0 = min(pos_range0)
            fetch_start0 = fetch_end0 - length
        else:
            fetch_start0 = pos_range0.stop
            fetch_end0 = fetch_start0 + length

        seq_beyond = self.fasta.fetch(chrom, fetch_start0, fetch_end0)

        if with_border_seq:
            border_seq = self.get_border_seq(is_bnd1)
            if endis5:
                return seq_beyond + border_seq
            else:
                return border_seq + seq_beyond
        else:
            return seq_beyond

    @functools.cache
    def get_seq_beyond_bnd_alt(self, is_bnd1, length, with_border_seq):
        # get parameters
        endis5 = self.endis5_bnd1 if is_bnd1 else self.endis5_bnd2
        (chrom_mate, pos_range0_mate, endis5_mate) = self.get_params(
            is_bnd1=(not is_bnd1))
        fetch_length = length - len(self.inserted_seq)

        # fetch_start, fetch_end
        if endis5_mate:
            fetch_start = pos_range0_mate.start
            fetch_end = fetch_start + fetch_length
        else:
            fetch_end = pos_range0_mate.start + 1
            fetch_start = fetch_end - fetch_length

        # mateside_seq
        if endis5 == endis5_mate:
            mateside_seq = Bio.Seq.reverse_complement(
                self.fasta.fetch(chrom_mate, fetch_start, fetch_end))
        else:
            mateside_seq = self.fasta.fetch(chrom_mate, fetch_start,
                                            fetch_end)

        # inserted_seq
        if is_bnd1:
            inserted_seq = ''.join(self.inserted_seq)
        else:
            if endis5 == endis5_mate:
                inserted_seq = Bio.Seq.reverse_complement(
                    ''.join(self.inserted_seq))
            else:
                inserted_seq = ''.join(self.inserted_seq)

        # result
        if endis5:
            seq_beyond = mateside_seq + inserted_seq
        else:
            seq_beyond = inserted_seq + mateside_seq

        # add border_seq
        if with_border_seq:
            border_seq = self.get_border_seq(is_bnd1)
            if endis5:
                return seq_beyond + border_seq
            else:
                return border_seq + seq_beyond
        else:
            return seq_beyond

        return result

    ###########################################

    def set_equivs(self):
        self._equivs = get_bnds_equivalents(self)

    def get_equivs(self):
        if self._equivs is None:
            self.set_equivs()
        return self._equivs

    def get_sorted_equivs(self, mode):
        return sort_equivs(self.get_equivs(), mode)

    ###############

    def set_pos_range0s(self):
        poslist0_bnd1 = list()
        poslist0_bnd2 = list()
        for bnds in self.get_equivs():
            poslist0_bnd1.append(bnds.pos_bnd1 - 1)
            poslist0_bnd2.append(bnds.pos_bnd2 - 1)

        if self.endis5_bnd1:
            self._pos_range0_bnd1 = range(max(poslist0_bnd1),
                                          min(poslist0_bnd1) - 1,
                                          -1)
        else:
            self._pos_range0_bnd1 = range(min(poslist0_bnd1),
                                          max(poslist0_bnd1) + 1,
                                          1)

        if self.endis5_bnd2:
            self._pos_range0_bnd2 = range(max(poslist0_bnd2),
                                          min(poslist0_bnd2) - 1,
                                          -1)
        else:
            self._pos_range0_bnd2 = range(min(poslist0_bnd2),
                                          max(poslist0_bnd2) + 1,
                                          1)

    def get_pos_range0_bnd1(self):
        """Returns a directional range"""

        if self._pos_range0_bnd1 is None:
            self.set_pos_range0s()
        return self._pos_range0_bnd1

    def get_pos_range0_bnd2(self):
        """Returns a directional range"""

        if self._pos_range0_bnd2 is None:
            self.set_pos_range0s()
        return self._pos_range0_bnd2

    @functools.cache
    def get_flank_range0(self, mode, is_bnd1, flanklen):
        """Args:
            mode: Must be one of "par", "bnd_prox", or "bnd_dist".
                par: partner side
                bnd_prox: breakend side, proximal (adjacent to the most 
                    advanced border)
                bnd_dist: breakend side, distal (adjacent to the most 
                    retracted border)
        """
        assert mode in ('par', 'bnd_prox', 'bnd_dist')

        # modify flanklen for breakend-side flanks
        if mode in ('bnd_prox', 'bnd_dist'):
            flanklen = flanklen - 1

        # get parameters
        pos_range0 = (self.get_pos_range0_bnd1()
                      if is_bnd1 else
                      self.get_pos_range0_bnd2())

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

    def get_bnd1adv_form(self):
        return pick_bnd1adv_form(self.get_equivs())

    def get_bnd2adv_form(self):
        return pick_bnd2adv_form(self.get_equivs())

    def set_microhomologies(self):
        (self._homlen, 
         self._homseq) = get_microhomology_spec(self.get_equivs())

    def get_homlen(self):
        if self._homlen is None:
            self.set_microhomologies()
        return self._homlen

    def get_homseq(self):
        if self._homseq is None:
            self.set_microhomologies()
        return self._homseq
    
    ###################################################################

    def copy(self):
        return Breakends(
            chrom_bnd1=self.chrom_bnd1, 
            pos_bnd1=self.pos_bnd1, 
            endis5_bnd1=self.endis5_bnd1,
            chrom_bnd2=self.chrom_bnd2, 
            pos_bnd2=self.pos_bnd2, 
            endis5_bnd2=self.endis5_bnd2,
            fasta=self.fasta, chromdict=self.chromdict,
            inserted_seq=self.inserted_seq.copy(), svtype=self.svtype,
            score=self.score)

    def identical(self, other):
        return hash(self) == hash(other)

    def sameseq1(self, other):
        return self.get_bnd1adv_form() == other.get_bnd1adv_form()

    def sameseq2(self, other):
        for key in ('chrom_bnd1', 'endis5_bnd1', 'chrom_bnd2', 'endis5_bnd2'):
            if getattr(self, key) != getattr(other, key):
                return False

        target = self.copy()
        query = other.copy()

        if target.endis5_bnd1:
            goal_bnd1 = max(target.pos_bnd1, query.pos_bnd1)
        else:
            goal_bnd1 = min(target.pos_bnd1, query.pos_bnd1)

        if target.endis5_bnd2:
            goal_bnd2 = max(target.pos_bnd2, query.pos_bnd2)
        else:
            goal_bnd2 = min(target.pos_bnd2, query.pos_bnd2)

        for bnds in (target, query):
            while True:
                if bnds.pos_bnd1 == goal_bnd1:
                    break
                else:
                    bnds.retract_bnd1()
                    continue
            while True:
                if bnds.pos_bnd2 == goal_bnd2:
                    break
                else:
                    bnds.retract_bnd2()
                    continue

        return target.equal(query)

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
            newbase = self.fasta.fetch(self.chrom_bnd1, newpos_bnd1 - 1, 
                                       newpos_bnd1)
            if len(self.inserted_seq) > 0:
                oldbase = (self.inserted_seq[-1] 
                           if self.endis5_bnd1 else 
                           self.inserted_seq[0])
            else:
                oldbase = self._convert_seq_between_bnds(
                    self.fasta.fetch(self.chrom_bnd2, self.pos_bnd2 - 1, 
                                     self.pos_bnd2))

            return (newbase == oldbase)

        def action(self, newpos_bnd1, newpos_bnd2):
            self.pos_bnd1 = newpos_bnd1
            self.pos_bnd2 = newpos_bnd2
            if len(self.inserted_seq) > 0:
                self.score += 1
                if self.endis5_bnd1:
                    del self.inserted_seq[-1]
                else:
                    del self.inserted_seq[0]

        # main
        newpos_bnd1, newpos_bnd2 = get_newpos(self)
        if self._newpos_range_check(newpos_bnd1, newpos_bnd2):
            if base_match_check(self, newpos_bnd1):
                action(self, newpos_bnd1, newpos_bnd2)
                success = True
            else:
                success = False
        else:
            success = False

        return success

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
                self.fasta.fetch(self.chrom_bnd2, newpos_bnd2 - 1, 
                                 newpos_bnd2))

            if len(self.inserted_seq) > 0:
                oldbase = (self.inserted_seq[0] 
                           if self.endis5_bnd1 else 
                           self.inserted_seq[-1])
            else:
                oldbase = self.fasta.fetch(self.chrom_bnd1, self.pos_bnd1 - 1, 
                                           self.pos_bnd1)

            return (newbase == oldbase)

        def action(self, newpos_bnd1, newpos_bnd2):
            self.pos_bnd1 = newpos_bnd1
            self.pos_bnd2 = newpos_bnd2
            if len(self.inserted_seq) > 0:
                self.score += 1
                if self.endis5_bnd1:
                    del self.inserted_seq[0]
                else:
                    del self.inserted_seq[-1]

        # main
        newpos_bnd1, newpos_bnd2 = get_newpos(self)
        if self._newpos_range_check(newpos_bnd1, newpos_bnd2):
            if base_match_check(self, newpos_bnd2):
                action(self, newpos_bnd1, newpos_bnd2)
                success = True
            else:
                success = False
        else:
            success = False

        return success

    def retract_bnd1(self):
        added_base = self.fasta.fetch(self.chrom_bnd1, self.pos_bnd1 - 1, 
                                      self.pos_bnd1)
        if self.endis5_bnd1:
            self.inserted_seq.append(added_base)
        else:
            self.inserted_seq.insert(0, added_base)

        self.pos_bnd1 = self._get_retracted_pos_bnd1()
        self.score -= 1

    def retract_bnd2(self):
        added_base = self._convert_seq_between_bnds(
            self.fasta.fetch(self.chrom_bnd2, self.pos_bnd2 - 1, 
                             self.pos_bnd2))

        if self.endis5_bnd1:
            self.inserted_seq.insert(0, added_base)
        else:
            self.inserted_seq.append(added_base)

        self.pos_bnd2 = self._get_retracted_pos_bnd2()
        self.score -= 1

    #######################################################

    def _set_svtype(self):
        if self.chrom_bnd1 != self.chrom_bnd2:
            self.svtype = 'TRA'
        else:
            if (not self.endis5_bnd1) and self.endis5_bnd2:
                self.svtype = 'DEL'
            elif self.endis5_bnd1 and (not self.endis5_bnd2):
                self.svtype = 'DUP'
            elif self.endis5_bnd1 == self.endis5_bnd2:
                self.svtype = 'INV'

    def _newpos_range_check(self, newpos_bnd1, newpos_bnd2):
        return (newpos_bnd1 >= 1 and 
                newpos_bnd1 <= self.chromdict[self.chrom_bnd1] and
                newpos_bnd2 >= 1 and 
                newpos_bnd2 <= self.chromdict[self.chrom_bnd2])

    def _get_advanced_pos_bnd1(self):
        if self.endis5_bnd1:
            return self.pos_bnd1 - 1
        else:
            return self.pos_bnd1 + 1

    def _get_retracted_pos_bnd1(self):
        if self.endis5_bnd1:
            return self.pos_bnd1 + 1
        else:
            return self.pos_bnd1 - 1

    def _get_advanced_pos_bnd2(self):
        if self.endis5_bnd2:
            return self.pos_bnd2 - 1
        else:
            return self.pos_bnd2 + 1

    def _get_retracted_pos_bnd2(self):
        if self.endis5_bnd2:
            return self.pos_bnd2 + 1
        else:
            return self.pos_bnd2 - 1

    def _convert_seq_between_bnds(self, seq):
        return convert_seq_between_bnds(seq, 
                                        self.endis5_bnd1, 
                                        self.endis5_bnd2)


#######################################################


class NonStandardSVAlt(Exception):
    pass


class NonStandardSVRecord(Exception):
    pass


def convert_seq_between_bnds(seq, endis5_bnd1, endis5_bnd2):
    if endis5_bnd1 == endis5_bnd2:
        return Bio.Seq.reverse_complement(seq)
    else:
        return seq


#######################################################


def get_bnds_equivalents(bnds):
    """Returns: A list of Breakends objects, all equivalent to the 
        input object, sorted such that the first item is the most advanced 
        form with respect to bnd1, and the last item is the most advanced 
        form with respect to bnd2.
    """

    input_copy = bnds.copy()

    bnds_list_bnd1adv = list()
    bnds_list_bnd1adv.append(input_copy)
    while True:
        new = bnds_list_bnd1adv[-1].copy()
        success = new.advance_bnd1()
        if success:
            bnds_list_bnd1adv.append(new)
            continue
        else:
            break

    bnds_list_bnd2adv = list()
    bnds_list_bnd2adv.append(input_copy)
    while True:
        new = bnds_list_bnd2adv[-1].copy()
        success = new.advance_bnd2()
        if success:
            bnds_list_bnd2adv.append(new)
            continue
        else:
            break

    bnds_list_bnd1adv.reverse()
    bnds_list = bnds_list_bnd1adv[:-1] + [input_copy] + bnds_list_bnd2adv[1:]
    max_score = max(x.score for x in bnds_list)

    bnds_equivs = [x for x in bnds_list if x.score == max_score]

    return bnds_equivs


def sort_equivs(bnds_equivs, mode):
    assert mode in ('bnd1adv_first', 'bnd1adv_last', 'bnd2adv_first', 
                    'bnd2adv_last'), 'Invalid "mode" argument'

    if mode == 'bnd1adv_first':
        reverse = (not bnds_equivs[0].endis5_bnd1)
        key = lambda bnds: bnds.pos_bnd1
    elif mode == 'bnd1adv_last':
        reverse = bnds_equivs[0].endis5_bnd1
        key = lambda bnds: bnds.pos_bnd1
    elif mode == 'bnd2adv_first':
        reverse = (not bnds_equivs[0].endis5_bnd2)
        key = lambda bnds: bnds.pos_bnd2
    elif mode == 'bnd2adv_last':
        reverse = bnds_equivs[0].endis5_bnd2
        key = lambda bnds: bnds.pos_bnd2

    return sorted(bnds_equivs, key=key, reverse=reverse)


def pick_bnd1adv_form(bnds_equivs):
    if bnds_equivs[0].endis5_bnd1:
        return min(bnds_equivs, key=(lambda bnds: bnds.pos_bnd1))
    else:
        return max(bnds_equivs, key=(lambda bnds: bnds.pos_bnd1))


def pick_bnd2adv_form(bnds_equivs):
    if bnds_equivs[0].endis5_bnd2:
        return min(bnds_equivs, key=(lambda bnds: bnds.pos_bnd2))
    else:
        return max(bnds_equivs, key=(lambda bnds: bnds.pos_bnd2))


def get_microhomology_spec(bnds_equivs):
    homlen = len(bnds_equivs) - 1
    if homlen == 0:
        homseq = ''
    else:
        bnds = bnds_equivs[0]
        if bnds.endis5_bnd1:
            pos_list = [
                x.pos_bnd1 for x in 
                sorted(bnds_equivs, key = lambda x: x.pos_bnd1)[:-1]]
        else:
            pos_list = [
                x.pos_bnd1 for x in 
                sorted(bnds_equivs, key = lambda x: x.pos_bnd1)[1:]]
        homseq = bnds.fasta.fetch(bnds.chrom_bnd1, pos_list[0]-1, 
                                  pos_list[-1])

    return homlen, homseq


#######################################################

def get_bnds_from_vr(vr, fasta, chromdict):
    """
    Raises:
        If input vr does not conform to a known SV variant record format 
            (including a valid non-SV variant record)
    """

    assert len(vr.alts) == 1, (
        f'Multiallelic variant record is not allowed:\n{vr}')

    vr_svinfo = get_vr_svinfo_standard_vr(vr, fasta, chromdict)
    bnds = get_bnds_from_vr_svinfo(vr, vr_svinfo, fasta, chromdict)

    return bnds


########################################################

def get_bnds_from_vr_svinfo(vr, vr_svinfo, fasta, chromdict):
    def warn():
        LOGGER.warning(f'"t" portion of SV ALT string is not an extension '
                       f'of REF string for this variant record:\n{vr}')

    if vr_svinfo['is_bnd1']:
        chrom_bnd1 = vr.contig
        chrom_bnd2 = vr_svinfo['chrom_mate']
        endis5_bnd1 = vr_svinfo['current_endis5']
        endis5_bnd2 = vr_svinfo['mate_endis5']

        if endis5_bnd1:
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
        endis5_bnd1 = vr_svinfo['mate_endis5']
        endis5_bnd2 = vr_svinfo['current_endis5']

        pos_bnd1 = vr_svinfo['pos_mate']

        if endis5_bnd2:
            if vr_svinfo['t'][-1] == vr_svinfo['ref']:
                pos_bnd2 = vr.pos
                inserted_seq = list(
                    convert_seq_between_bnds(vr_svinfo['t'][:-1], 
                                             endis5_bnd1, endis5_bnd2))
            else:
                warn()
                pos_bnd2 = vr.pos + 1
                inserted_seq = list(
                    convert_seq_between_bnds(vr_svinfo['t'], 
                                             endis5_bnd1, endis5_bnd2))
        else:
            if vr_svinfo['t'][0] == vr_svinfo['ref']:
                pos_bnd2 = vr.pos
                inserted_seq = list(
                    convert_seq_between_bnds(vr_svinfo['t'][1:], 
                                             endis5_bnd1, endis5_bnd2))
            else:
                warn()
                pos_bnd2 = vr.pos - 1
                inserted_seq = list(
                    convert_seq_between_bnds(vr_svinfo['t'], 
                                             endis5_bnd1, endis5_bnd2))

    bnds = Breakends(chrom_bnd1=chrom_bnd1, pos_bnd1=pos_bnd1,
                     endis5_bnd1=endis5_bnd1,
                     chrom_bnd2=chrom_bnd2, pos_bnd2=pos_bnd2,
                     endis5_bnd2=endis5_bnd2,
                     inserted_seq=inserted_seq, fasta=fasta,
                     chromdict=chromdict)

    return bnds


def get_vr_svinfo_standard_vr(vr, fasta, chromdict):
    vr_svinfo = dict()
    vr_svinfo['ref'] = vr.ref

    try:
        (vr_svinfo['t'],
         vr_svinfo['chrom_mate'],
         vr_svinfo['pos_mate'],
         vr_svinfo['current_endis5'],
         vr_svinfo['mate_endis5']) = parse_sv_altstring(vr.alts[0])
    except NonStandardSVAlt as e:
        e_msg = f'{str(e)}\nInput variant record:\n{vr}'
        raise NonStandardSVRecord(e_msg)

    vr_svinfo['is_bnd1'] = get_is_bnd1(vr, vr_svinfo, chromdict)

    return vr_svinfo


def parse_sv_altstring(sv_altstring):
    mats = [common.RE_PATS['alt_bndstring_1'].match(sv_altstring), 
            common.RE_PATS['alt_bndstring_2'].match(sv_altstring)]
    mats_isNotNone = [(x is not None) for x in mats]

    nTrue = mats_isNotNone.count(True)
    if nTrue == 0: # not a bnd string
        raise NonStandardSVAlt(f'ALT string "{sv_altstring}" does not match '
                               f'the standard SV string pattern.') 

    elif nTrue == 1:
        mat = next(itertools.compress(mats, mats_isNotNone))
        t = mat.group('t') # t : according to VCF spec documentation
        chrom_mate = mat.group('matechrom')
        pos_mate = int(mat.group('matepos'))

        if sv_altstring.startswith('[') or sv_altstring.startswith(']'):
            endtype_current_is5 = True
        else:
            endtype_current_is5 = False
        
        if mat.group('bracket1') == '[':
            endtype_mate_is5 = True
        else:
            endtype_mate_is5 = False

        sv_altstring_parsed = (t, chrom_mate, pos_mate, endtype_current_is5, 
                               endtype_mate_is5)

    elif nTrue == 2: # not a valid bnd string
        raise NonStandardSVAlt(f'ALT string "{sv_altstring}" matches both '
                               f'pat1 and pat2.')

    return sv_altstring_parsed


def get_is_bnd1(vr, vr_svinfo, chromdict):
    order = common.compare_coords(vr.contig, vr.pos, vr_svinfo['chrom_mate'], 
                                  vr_svinfo['pos_mate'], chromdict)
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



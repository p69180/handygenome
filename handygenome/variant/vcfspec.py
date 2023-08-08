import re
import itertools
import functools
import collections

import Bio.Align
import pyranges as pr

import handygenome.refgenome as refgenome
import handygenome.annotation.annotitem as annotitem
import handygenome.variant.repeat as librepeat
from handygenome.annotation.annotitem import AnnotItemInfoALTlist
from handygenome.annotation.annotitem import AnnotItemInfoSingle

#import handygenome.align.alignhandler as alignhandler
#import handygenome.align.realign as realign


DEFAULT_FETCH_EXTEND_LENGTH = 10
SV_ALTS = ('DEL', 'INS', 'DUP', 'INV', 'CNV', 'BND', 'TRA')
CPGMET_ALT = 'CPGMET'

PAT_BND1 = re.compile('^(?P<t>[^\[\]]+)(?P<bracket1>\[|\])(?P<matechrom>[^:]+):(?P<matepos>[0-9]+)(?P<bracket2>\[|\])$')
PAT_BND2 = re.compile('^(?P<bracket1>\[|\])(?P<matechrom>[^:]+):(?P<matepos>[0-9]+)(?P<bracket2>\[|\])(?P<t>[^\[\]]+)$')
PAT_NUCLEOBASES = re.compile('[ACGTNacgtn]+')

ALIGNER_DIFF = Bio.Align.PairwiseAligner(
    mode='global',
    match_score=3,
    mismatch_score=-3,

    query_open_gap_score=-3,
    query_extend_gap_score=0,

    target_open_gap_score=-3,
    target_extend_gap_score=0,
)


'''
My convention:
    Multiple ALTs in a single Vcfspec object is considered to represent
    separate haplotypes. Therefore, a single read (or rp or rpp) is expected 
    to support no more than one of ALTs of a single Vcfspec.
'''


class Vcfspec:
    # I/O
    def __init__(
        self, chrom=None, pos=None, ref=None, alts=None, *,
        #somaticindex=1, germlineindexes=(0, 0),
        refver=None, fasta=None, components=None,
    ):
        # sanity check
        if alts is not None:
            if not isinstance(alts, (tuple, list)):
                raise Exception(f'"alts" argument must be a tuple or a list.')
        # chrom, pos, ref, alt
        self.chrom = chrom
        self.pos = int(pos)  # 221215: There are occasions where this value is np.int64, incompatible with VcfspecComponents.write
        self.ref = ref
        if alts is None:
            self.alts = None
        else:
            self.alts = tuple(alts)
        # others
        #self.somaticindex = somaticindex
        #self.germlineindexes = sorted(germlineindexes)
        self.refver = refver
        if fasta is None:
            if refver is None:
                raise Exception(f'When "fasta" is not set, "refver" must be set.')
            self.fasta = refgenome.get_default_fasta(self.refver)
        else:
            self.fasta = fasta
        # attributes not able to set in __init__
        self.is_leftmost = None
        self.is_rightmost = None
        self._is_parsimonious = None
        self._equivalents = None
        self._repeat_spec = None

        # set components
        if components is None:
            if self.alts is None:
                self.components = None
            else:
                self.components = VcfspecComponents.init_empty(len(self.alts))
        else:
            assert isinstance(components, VcfspecComponents)
            self.components = components

    @classmethod
    def from_vr(cls, vr, refver=None):
        if refver is None:
            refver = refgenome.infer_refver_vcfheader(vr.header)
        components = VcfspecComponents.from_vr(vr, refver)

        result = cls(
            chrom=vr.contig, pos=vr.pos, ref=vr.ref, alts=vr.alts,
            refver=refver,
            components=components,
        )
        return result

    def to_tuple(self):
        return (self.chrom, self.pos, self.ref, self.alts)
        #return annotitem._encode(
        #    (self.chrom, self.pos, self.ref, self.alts),
        #    handle_percent=True,
        #)

    @classmethod
    def from_tuple(cls, tup, refver):
        return cls(*tup, refver=refver)

    ################

    def __repr__(self):
#        if len(self.alts) == 1:
#            altstring = str(self.alts[0])
#        else:
#            altstring = str(list(self.alts))
#        return f'<Vcfspec ({self.chrom}:{self.pos} {self.ref}>{altstring})>'
        return f'Vcfspec(chrom={repr(self.chrom)}, pos={self.pos:,}, ref={repr(self.ref)}, alts={repr(self.alts)})'

    def __hash__(self):
        return hash((self.chrom, self.pos, self.ref, self.alts))

    def __eq__(self, other):
        return all(
            getattr(self, key) == getattr(other, key)
            for key in ('chrom', 'pos', 'ref', 'alts')
        )

    def show_components(self):
        for idx in range(len(self.alts)):
            print(self.get_monoalt(idx))
            print('Concat components:')
            for x in self.concat_components[idx]:
                print('\t', x, sep='')
            print('-----------------')

    def spawn(self, **kwargs):
        default_kwargs = {
            'chrom': self.chrom, 
            'pos': self.pos, 
            'ref': self.ref, 
            'alts': self.alts,
            #'somaticindex': self.somaticindex,
            #'germlineindexes': self.germlineindexes,
            'refver': self.refver,
            'fasta': self.fasta,
        }
        default_kwargs.update(kwargs)
        result = self.__class__(**default_kwargs)

#        for key in (
#            'is_leftmost',
#            'is_rightmost',
#            '_is_parsimonious',
#        ):
#            setattr(result, key, getattr(self, key))
#
#        for key in (
#            '_equivalents',
#            '_concat_components',
#            '_merge_components',
#        ):
#            setattr(result, key, getattr(self, key).copy())
#            # shallow copies of the original

        return result

    @property
    def pos0(self):
        return self.pos - 1

    start0 = pos0

    @property
    def end0(self):
        return self.pos0 + len(self.ref)

    @property
    def length(self):
        return len(self.ref)

    @property
    def alleles(self):
        return (self.ref,) + self.alts

#    @property
#    def germline(self):
#        alleles = self.alleles
#        return tuple(alleles[x] for x in self.germlineindexes)

    #@property
    #def somatic(self):
    #    return self.alleles[self.somaticindex]

    def iter_allele_indexes(self):
        return iter(range(len(self.alts) + 1))

    def iter_alt_indexes(self):
        return iter(range(len(self.alts)))

    def get_id(self):
        return '_'.join([self.chrom, str(self.pos), self.ref, '|'.join(self.alts)])

    def to_string(self):
        return f'{self.chrom}:{self.pos} {self.ref}>{"|".join(self.alts)}'

    def get_mutation_type(self, alt_index=0):
        return get_mutation_type(self.ref, self.alts[alt_index])

    #def get_mttype_firstalt(self):
    #    return self.get_mutation_type(0)

    def get_diff_seq(self):
        self.check_monoalt(raise_with_false=True)

        muttype = self.get_mutation_type(alt_index=0)
        if muttype == 'ins':
            return self.alts[0][1:]
        elif muttype == 'del':
            return self.ref[1:]
        else:
            return None

    def get_tuple(self):
        return (self.chrom, self.pos, self.ref, self.alts)

    # ranges
    @property
    def REF_range0(self):
        return range(self.pos0, self.end0)

    @property
    def range0(self):
        return range(self.pos0, self.end0)

    def get_range0(self):
        return range(self.pos0, self.end0)

    #def get_preflank_range0(self, flanklen=1):
    #    assert flanklen >= 1, f'"flanklen" argument must be at least 1.'

    @functools.cache
    def get_flank_range0s_equivalents(self, flanklen=1):
        preflank_candidates = [
            self.get_preflank_range0_monoalt_leftmost(alt_index=idx, flanklen=flanklen)
            for idx in range(len(self.alts))
        ]
        preflank_range0 = range(
            min(x.start for x in preflank_candidates),
            self.pos0
        )

        postflank_candidates = [
            self.get_postflank_range0_rightmost(alt_index=idx, flanklen=flanklen)
            for idx in range(len(self.alts))
        ]
        postflank_range0 = range(
            self.end0,
            max(x.stop for x in postflank_candidates)
        )

        if len(preflank_range0) == 0 or len(postflank_range0) == 0:
            raise Exception(
                f'Length 0 range.'
                f'self: {self}\n'
                f'preflank_candidates: {preflank_candidates}\n'
                f'postflank_candidates: {postflank_candidates}\n'
            )

        return preflank_range0, postflank_range0

    @functools.cache
    def get_flank_range0s(self, flanklen=1):
#        preflank_candidates = [
#            self.get_preflank_range0_monoalt(alt_index=idx, flanklen=flanklen)
#            for idx in range(len(self.alts))
#        ]
#        preflank_range0 = range(
#            min(x.start for x in preflank_candidates),
#            max(x.stop for x in preflank_candidates)
#        )
        preflank_range0 = self.get_preflank_range0_monoalt(flanklen=flanklen)

        postflank_range0 = self.get_postflank_range0(flanklen=flanklen)

        return preflank_range0, postflank_range0

    def get_preflank_range0_monoalt(self, alt_index=0, flanklen=1):
        assert flanklen >= 1, f'"flanklen" argument must be at least 1.'
        #if self.alts[0][0] == self.ref[0]:
        #    end0 = self.pos0 + 1
        #else:
        #    end0 = self.pos0
        return range(self.pos0 - flanklen, self.pos0)

    def get_preflank_range0_monoalt_leftmost(self, alt_index=0, flanklen=1):
        rng_list = list()

        leftmost_form = self.get_monoalt(alt_index).leftmost()
        rng_list.append(
            leftmost_form.get_preflank_range0_monoalt(flanklen=flanklen)
        )
        if not (self.components.is_empty or self.components.is_missing):
            for component_vcfspec in self.components[alt_index]:
                leftmost_form = component_vcfspec.leftmost()
                rng_list.append(
                    leftmost_form.get_preflank_range0_monoalt(flanklen=flanklen)
                )

        leftmost_range_start = min(x.start for x in rng_list)

        return range(
            min(leftmost_range_start, self.pos0 - flanklen), 
            self.pos0,
        )

    def get_postflank_range0(self, flanklen=1):
        assert flanklen >= 1, f'"flanklen" argument must be at least 1.'
        return range(self.end0, self.end0 + flanklen)

    def get_postflank_range0_rightmost(self, alt_index=0, flanklen=1):
        rng_list = list()

        rightmost_form = self.get_monoalt(alt_index).rightmost()
        rng_list.append(
            rightmost_form.get_postflank_range0(flanklen=flanklen)
        )
        if not (self.components.is_empty or self.components.is_missing):
            for component_vcfspec in self.components[alt_index]:
                rightmost_form = component_vcfspec.rightmost()
                rng_list.append(
                    rightmost_form.get_postflank_range0(flanklen=flanklen)
                )

        rightmost_range_stop = max(x.stop for x in rng_list)

        return range(self.end0, max(rightmost_range_stop, self.end0 + flanklen))

    # conversion
    #@functools.cached_property
    #def equivalents(self):
    #    return equivalents(self, fasta=self.fasta, fetch_extend_length=DEFAULT_FETCH_EXTEND_LENGTH)

    def set_equivalents(self, fetch_extend_length=DEFAULT_FETCH_EXTEND_LENGTH):
        self._equivalents = equivalents(self, fetch_extend_length=fetch_extend_length)

    def get_equivalents(self, fetch_extend_length=DEFAULT_FETCH_EXTEND_LENGTH):
        if self._equivalents is None:
            self.set_equivalents(fetch_extend_length=fetch_extend_length)
        return self._equivalents

    def leftmost(self, fetch_extend_length=DEFAULT_FETCH_EXTEND_LENGTH):
        if self.is_leftmost:
            return self
        else:
            return self.get_equivalents(fetch_extend_length=fetch_extend_length)[0]

    normalize = leftmost

    def check_normalized(self):
        return self == self.normalize()

    def rightmost(self, fetch_extend_length=DEFAULT_FETCH_EXTEND_LENGTH):
        if self.is_rightmost:
            return self
        else:
            return self.get_equivalents(fetch_extend_length=fetch_extend_length)[-1]

    def parsimonious(self):
        if self._is_parsimonious:
            return self
        else:
            return make_parsimonious(self)

    def check_parsimonious(self):
        if self._is_parsimonious is None:
            self._is_parsimonious = (self == self.parsimonious())
        return self._is_parsimonious

    def extend_right(self, length):
        added_seq = self.fasta.fetch(self.chrom, self.end0, self.end0 + length)
        result = self.spawn(pos=self.pos, ref=(self.ref + added_seq), alts=tuple(x + added_seq for x in self.alts))
        result._is_parsimonious = False
        return result

    def extend_left(self, length):
        added_seq = self.fasta.fetch(self.chrom, self.pos0 - length, self.pos0)
        result = self.spawn(pos=(self.pos - length), ref=(added_seq + self.ref), alts=tuple(added_seq + x for x in self.alts))
        result._is_parsimonious = False
        return result

    def merge(self, other):
        return merge(self, other)

    # misc
    def get_sortkey(self, chromdict):
        return (
            (chromdict.contigs.index(self.chrom), self.pos, self.ref) 
            + self.alts
        )

    def apply_to_vr(self, vr):
        vr.contig = self.chrom
        vr.pos = self.pos
        vr.ref = self.ref
        vr.alts = self.alts
        if not self.components.is_missing:
            self.components.write(vr)

    def get_alleleclasses(self):
        return (None,) + tuple(range(-1, len(self.alts) + 1))

    def iter_annotation_forms(self):
        for alt_index in self.iter_alt_indexes():
            component_index = alt_index + 1
            if component_index in self.components:
                component_vcfspecs = self.components[component_index]
                if len(component_vcfspecs) == 0:
                    yield self.spawn(alts=(self.alts[alt_index],))
                elif len(component_vcfspecs) == 1:
                    yield component_vcfspecs[0]
                else:
                    yield concat_list(component_vcfspecs, distance_le=None)[0]
            else:
                yield self.spawn(alts=(self.alts[alt_index],))

    def iter_monoalts(self):
        for alt_index in self.iter_alt_indexes():
            yield self.get_monoalt(alt_index=alt_index)

    def get_monoalt(self, alt_index=0):
        if (
            self.components.is_empty 
            or self.components.is_missing 
        ):
            #new_components = self.components
            new_components = VcfspecComponents.init_empty(num_alt=1)
        else:
            allele_index = alt_index + 1
            new_components = VcfspecComponents.init_nonmissing()
            new_components[1] = self.components[allele_index]

        result = self.spawn(
            alts=(self.alts[alt_index],),
            components=new_components,
        )
        return result

    def check_without_N(self):
        return (
            set('nN').isdisjoint(self.ref) and 
            all(set('nN').isdisjoint(x) for x in self.alts)
        )

    def to_hgvsg(self, alt_index=0):
        chrom = self.chrom
        pos = self.pos
        ref = self.ref
        alt = self.alts[alt_index]
        muttype = self.get_mutation_type(alt_index)

        if muttype == 'snv':
            result = f'{chrom}:g.{pos}{ref}>{alt}'
        elif muttype == 'mnv':
            pos2 = pos + (len(ref) - 1)
            result = f'{chrom}:g.{pos}_{pos2}delins{alt}'
        elif muttype == 'ins':
            inserted_seq = alt[1:]
            result = f'{chrom}:g.{pos}_{pos+1}ins{inserted_seq}'
        elif muttype == 'del':
            pos1 = pos + 1
            pos2 = pos + (len(ref) - 1)
            if pos1 == pos2:
                result = f'{chrom}:g.{pos1}del'
            else:
                result = f'{chrom}:g.{pos1}_{pos2}del'
        elif muttype == 'delins':
            pos1 = pos
            pos2 = pos + (len(ref) - 1)
            if pos1 == pos2:
                result = f'{chrom}:g.{pos1}delins{alt}'
            else:
                result = f'{chrom}:g.{pos1}_{pos2}delins{alt}'

        return result

    def to_gr(self):
        return pr.from_dict(
            {
                'Chromosome': [self.chrom], 
                'Start': [self.pos0], 
                'End': [self.end0],
                'Id': [self.get_id()],
            }
        )

    def check_monoalt(self, raise_with_false=True):
        if raise_with_false:
            if len(self.alts) != 1:
                raise Exception('Vcfspec is not with a single ALT: {self}')
        else:
            return len(self.alts) == 1

    @functools.cached_property
    def repeat_specs(self):
        return [
            x.get_repeat_spec_monoalt()
            for x in self.iter_monoalts()
        ]

    def get_repeat_spec_monoalt(self):
        # sanity check
        self.check_monoalt(raise_with_false=True)

        # main
        normalized = self.normalize()
        diff_seq = normalized.get_diff_seq()
        if diff_seq is None:
            return None

        repeat_unit = librepeat.decompose_repeat(diff_seq)[0][0]
        repeat_range0 = librepeat.search_for_repeat(
            query_seq=repeat_unit, 
            chrom=normalized.chrom,
            start0=(normalized.pos0 + 1),
            fasta=normalized.fasta,
        )
        if repeat_range0 is None:
            return None
        else:
            return librepeat.Repeat(
                chrom=normalized.chrom, 
                start0=repeat_range0.start, 
                unit=repeat_unit,
                count=(len(repeat_range0) // len(repeat_unit)),
                fasta=normalized.fasta,
            )

    def check_same_firstbase(self):
        data = list()
        for alt_index in range(len(self.alts)):
            data.append(self.ref[0] == self.alts[alt_index][0])
        data = set(data)
        if len(data) != 1:
            raise Exception(f'Vcfspec is multi-ALT but ')

    def get_readspan_range0_parsimonious(self):
        #self.check_monoalt(raise_with_false=True)
        start = self.pos0
        muttype = self.get_mutation_type(alt_index=0)
        if muttype in ('ins', 'del'):
            stop = self.end0 + 1
        else:
            stop = self.end0

        return range(start, stop)

    def get_readspan_range0_monoalt(self):
        if len(self.components[1]) == 0:
            return self.get_readspan_range0_parsimonious()
        else:
            start = self.pos0
            rightmost_component = sorted(self.components[1], key=(lambda x: x.pos))[-1]
            rightmost_component_rng0 = rightmost_component.get_readspan_range0_parsimonious()
            stop = max(self.end0, rightmost_component_rng0.stop)
            return range(start, stop)

    @functools.cached_property
    def readspan_range0(self):
        starts = list()
        stops = list()
        for monoalt in self.iter_monoalts():
            rng0 = monoalt.get_readspan_range0_monoalt()
            starts.append(rng0.start)
            stops.append(rng0.stop)
        return range(min(starts), max(stops))

    @functools.cached_property
    def compared_seq_prepad(self):
        if self.readspan_range0.start < self.start0:
            return self.fasta.fetch(self.chrom, self.readspan_range0.start, self.start0)
        else:
            return ''

    @functools.cached_property
    def compared_seq_postpad(self):
        if self.readspan_range0.stop > self.end0:
            return self.fasta.fetch(self.chrom, self.end0, self.readspan_range0.stop)
        else:
            return ''

    def get_compared_seq(self, allele_index):
        return self.alleles[allele_index] + self.compared_seq_postpad

    @functools.cached_property
    def compared_seq_dict(self):
        return {
            allele_index: self.get_compared_seq(allele_index)
            for allele_index in self.iter_allele_indexes()
        }


#class VcfspecConcatComponents(collections.UserList):
#    def encode(self):
#        return annotitem._encode(tuple((x.chrom, x.pos, x.ref, x.alts) for x in self))
#        
#    @classmethod
#    def from_raw_value(cls, vr_value, refver):
#        result = cls()
#        decoded = annotitem._decode(vr_value, handle_percent=True, convert_key=True)
#        for tup in decoded:
#            result.append(Vcfspec(*tup, refver=refver))
#        return result
            

class VcfspecComponents(AnnotItemInfoSingle):
    meta = {
        "ID": "vcfspec_components",
        "Number": "1",
        "Type": "String",
        "Description": "Merge/Concat component information of Vcfspec",
    }

    @classmethod
    def from_vr(cls, vr, refver=None):
        if refver is None:
            refver = refgenome.infer_refver_vcfheader(vr.header)
        return cls.from_vr_base(vr, refver)

    def write(self, vr):
        self.write_base(vr)

    @classmethod
    def from_raw_value(cls, vr_value, refver):
        result = cls.init_nonmissing()
        decoded = annotitem._decode(vr_value, handle_percent=True, convert_key=True)
            # decoded: a dict
                # key: allele_index
                # value: list of vcfspec-encoding-tuples
        for allele_index, l in decoded.items():
            concat_components = list()
            for tup in l:
                concat_components.append(Vcfspec.from_tuple(tup, refver))
            result[allele_index] = concat_components

        return result

    @classmethod
    def init_empty(cls, num_alt):
        result = cls.init_nonmissing()
        for allele_index in range(1, num_alt + 1):
            result[allele_index] = list()
        return result

    @property
    def is_empty(self):
        return (
            #(not self.is_missing)
            all(
                len(val) == 0
                for val in self.values()
            )
        )

    def subset(self, allele_index_list):
        result = self.__class__.init_nonmissing()
        for allele_index in allele_index_list:
            result[allele_index] = self[allele_index]
        return result

    def encode(self):
        data = dict()
        for allele_index, concat_components in self.items():
            data[allele_index] = [x.to_tuple() for x in concat_components]
        return annotitem._encode(data, handle_percent=True)


class IdentityVcfspec:
    pass


#def check_vcfspec_monoalt(vcfspec):
#    if len(vcfspec.alts) != 1:
#        raise Exception('The input vcfspec must be with single ALT.')

#check_vcfspec_monoallele = check_vcfspec_monoalt


def get_mutation_type(ref, alt):
    if PAT_NUCLEOBASES.fullmatch(alt) is None:
        if any(
            (re.fullmatch(f'<{x}(:.+)?>', alt) is not None) for x in SV_ALTS
        ):
            mttype = 'sv'
        elif (
            (PAT_BND1.fullmatch(alt) is not None) or 
            (PAT_BND2.fullmatch(alt) is not None)
        ):
            mttype = 'sv'
        elif alt == f'<{CPGMET_ALT}>':
            mttype = 'cpgmet'
        elif alt == '*':
            mttype = 'overlap_del'
        else:
            raise Exception(f'Unexpected symbolic ALT allele: {repr(alt)}')
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

get_mttype = get_mutation_type


def merge(vcfspec1, vcfspec2):
    """Intended for use with vcfspecs on different contigs(haplotypes).
    Args:
        Input vcfspecs may overlap.
    Returns:
        A new vcfspec with multiple ALTs, containing all ALTs of input vcfspecs
    """
    assert vcfspec1.chrom == vcfspec2.chrom, f'Two vcfspecs must be on the same chromosome.'

    vcfspec1_edit = vcfspec1.spawn()
    vcfspec2_edit = vcfspec2.spawn()
    # equalize left margins
    start_diff = vcfspec1_edit.pos - vcfspec2_edit.pos
    if start_diff > 0:
        vcfspec1_edit = vcfspec1_edit.extend_left(start_diff)
    elif start_diff < 0:
        vcfspec2_edit = vcfspec2_edit.extend_left(-start_diff)
    # equalize right margins
    end_diff = vcfspec1_edit.end0 - vcfspec2_edit.end0
    if end_diff > 0:
        vcfspec2_edit = vcfspec2_edit.extend_right(end_diff)
    elif end_diff < 0:
        vcfspec1_edit = vcfspec1_edit.extend_right(-end_diff)
    # get new alt and new components
    result_vcfspec_components_data = dict()
    for original_vcfspec, edited_vcfspec in (
        (vcfspec1, vcfspec1_edit),
        (vcfspec2, vcfspec2_edit),
    ):
        #for idx, (original_alt, edited_alt) in enumerate(zip(original_vcfspec.alts, edited_vcfspec.alts)):
        for alt_index, edited_alt in enumerate(edited_vcfspec.alts):
            if edited_alt in result_vcfspec_components_data.keys():
                continue

            allele_index = alt_index + 1
            component_unit = original_vcfspec.components[allele_index]
            if len(component_unit) == 0:
                result_vcfspec_components_data[edited_alt] = [original_vcfspec.get_monoalt(alt_index=alt_index)]
            else:
                result_vcfspec_components_data[edited_alt] = component_unit.copy()

    new_alts = sorted(
        result_vcfspec_components_data.keys(),
        key=(lambda x: 
            (
                min(result_vcfspec_components_data[x], key=(lambda y: y.pos)).pos,  
                    # pos of the leftmost component unit
                x,  
                    # edited ALT string itself
            )
        )
    )
    new_components = VcfspecComponents.init_empty(len(new_alts))
    for alt_index, edited_alt in enumerate(new_alts):
        allele_index = alt_index + 1
        new_components[allele_index] = result_vcfspec_components_data[edited_alt]

    # make result
    return vcfspec1_edit.spawn(alts=new_alts, components=new_components)


#def concat_new(vcfspec1, vcfspec2):
#    """Intended for use with vcfspecs on the same contig(haplotype).
#    Args:
#        Must be vcfspecs with single ALT. They must not overlap each other.
#    Returns:
#        A new vcfspec with single ALT, created by filling gap between two input vcfspecs
#    """
#
#    def calc_differing_range(vcfspec):
#        """Calculates actually differing range.
#        example: pos0=50, REF="AT", ALT="ATCCG" -> offset is 2; returns range(52, 55)
#        example: pos0=100, REF="A", ALT="ACC" -> offset is 1; returns range(101, 101)
#        0-length range may be returned!! This means insertion between "range.start" and "range.start - 1"
#        """
#        min_length = min(len(vcfspec.ref, vcfspec.alts[0]))
#        overlapping_ref = vcfspec.ref[:min_length]
#        overlapping_alt = vcfspec.alts[0][:min_length]
#        if overlapping_ref == overlapping_alt:
#            offset = min_length
#        else:
#            for offset, (ref_base, alt_base) in enumerate(zip(overlapping_ref, overlapping_alt)):
#                if ref_base != alt_base:
#                    break
#
#        return range(vcfspec.start0 + offset, vcfspec.end0), offset
#
#    # sanity check
#    assert vcfspec1.chrom == vcfspec2.chrom, f'Two vcfspecs must be on the same chromosome.'
#    vcfspec1.check_monoalt(raise_with_false=True)
#    vcfspec2.check_monoalt(raise_with_false=True)
#
#    # calculate actually differing ranges
#    differing_range0_vcfspec1, offset_vcfspec1 = calc_differing_range(vcfspec1)
#    differing_range0_vcfspec2, offset_vcfspec2 = calc_differing_range(vcfspec2)
#
#    # sanity check 2
#    if differing_range0_vcfspec1 == differing_range0_vcfspec2:
#        raise Exception(f'Input vcfspecs are incompatible to concat; differing ranges are identical')
#
#    # set left/right
#    if differing_range0_vcfspec1.start >= differing_range0_vcfspec2.stop:
#        left_vcfspec = vcfspec2
#        left_vcfspec_differing_range0 = differing_range0_vcfspec2
#        left_offset = offset_vcfspec2
#        right_vcfspec = vcfspec1
#        right_vcfspec_differing_range0 = differing_range0_vcfspec1
#        right_offset = offset_vcfspec1
#    elif differing_range0_vcfspec2.start >= differing_range0_vcfspec1.stop:
#        left_vcfspec = vcfspec1
#        left_vcfspec_differing_range0 = differing_range0_vcfspec1
#        left_offset = offset_vcfspec1
#        right_vcfspec = vcfspec2
#        right_vcfspec_differing_range0 = differing_range0_vcfspec2
#        right_offset = offset_vcfspec2
#    else:
#        raise Exception(f'Input vcfspecs are incompatible to concat; differing ranges overlap with each other')


def concat(vcfspec1, vcfspec2):
    """Intended for use with vcfspecs on the same contig(haplotype).
    Args:
        Must be vcfspecs with single ALT. They must not overlap each other.
    Returns:
        A new vcfspec with single ALT, created by filling gap between two input vcfspecs
    """
    # sanity check
    assert vcfspec1.chrom == vcfspec2.chrom, f'Two vcfspecs must be on the same chromosome.'
    vcfspec1.check_monoalt(raise_with_false=True)
    vcfspec2.check_monoalt(raise_with_false=True)

    # set left and right
    if vcfspec1.pos < vcfspec2.pos:
        left = vcfspec1
        right = vcfspec2
    else:
        left = vcfspec2
        right = vcfspec1

    # main
    if right.start0 < left.end0:
        if right.ref[0] == right.alts[0][0]:
            first_pos_same = right
            other = left
        elif left.ref[0] == left.alts[0][0]:
            first_pos_same = left
            other = right
        else:
            raise Exception(f'Two vcfspecs are overlapping and showing unexpected patterns:\nleft: {left}\nright: {right}')

        if other.end0 != first_pos_same.start0 + 1:
            raise Exception(f'Two vcfspecs are overlapping and showing unexpected patterns:\nleft: {left}\nright: {right}')

        new_pos = other.pos
        new_ref = other.ref + first_pos_same.ref[1:]
        new_alt = other.alts[0] + first_pos_same.alts[0][1:]
    else:
        new_pos = left.pos
        if right.start0 == left.end0:
            new_ref = left.ref + right.ref
            new_alt = left.alts[0] + right.alts[0]
        else:
            intervening_ref = left.fasta.fetch(left.chrom, left.end0, right.start0)
            new_ref = left.ref + intervening_ref + right.ref
            new_alt = left.alts[0] + intervening_ref + right.alts[0]

    result = left.spawn(pos=new_pos, ref=new_ref, alts=(new_alt,))
    result.components = VcfspecComponents.init_empty(1)

    components_element = list()
    for input_vcfspec in (vcfspec1, vcfspec2):
        if len(input_vcfspec.components[1]) == 0:
            components_element.append(input_vcfspec.spawn())
        else:
            components_element.extend(input_vcfspec.components[1])
    result.components[1] = sorted(components_element, key=(lambda x: x.pos))
    return result


def concat_list(vcfspec_list, distance_le=None):
    """Args:
        distance_le: Distance less than or equal to this value should
            be concatenated. When REF range of two vcfspecs are adjacent to 
            each other, distance is 0.
            If None, all vcfspecs are concatenated regardless of distance.
    """
    new_vcfspec_list = list()

    if distance_le is None:
        unit_func = concat
    else:
        def unit_func(vcfspec1, vcfspec2):
            num_intervening_bases = vcfspec2.start0 - vcfspec1.end0
            if num_intervening_bases <= distance_le:
                return concat(vcfspec1, vcfspec2)
            else:
                new_vcfspec_list.append(vcfspec1)
                return vcfspec2

    final = functools.reduce(unit_func, sorted(vcfspec_list, key=(lambda x: x.pos)))
    new_vcfspec_list.append(final)

    return tuple(new_vcfspec_list)


def split(vcfspec, distance_ge, show_alignment=False):
    """Args:
        distance_ge: Distance greater than or equal to this value should
            be split. When REF range of two vcfspecs are adjacent to 
            each other, distance is 0.
    """
    #check_vcfspec_monoalt(vcfspec)
    vcfspec.check_monoalt(raise_with_false=True)

    alignment = alignhandler.alignment_tiebreaker(
        realign.ALIGNER_MAIN.align(vcfspec.ref, vcfspec.alts[0]),
        raise_with_failure=True,
    )
    vcfspec_list = alignhandler.alignment_to_vcfspec(
        alignment, target_start0=vcfspec.pos0, chrom=vcfspec.chrom, fasta=vcfspec.fasta, strip_query_gaps=False
    )
    merged_vcfspec_list = concat(vcfspec_list, merging_distance_le=(distance_ge - 1))
    return merged_vcfspec_list


def diff(vcfspec1, vcfspec2):
    vcfspec1.check_monoalt(raise_with_false=True)
    vcfspec2.check_monoalt(raise_with_false=True)

    merged_vcfspec = merge(vcfspec1, vcfspec2)
    merged_vcfspec.extend_left(2)
    merged_vcfspec.extend_right(2)

    alns = ALIGNER_DIFF.align(merged_vcfspec.alts[0], merged_vcfspec.alts[1])


def diff_tiebreaker(aln):
    pass


##############
# equivalent #
##############

def leftmost(vcfspec, fasta, fetch_extend_length=DEFAULT_FETCH_EXTEND_LENGTH):
    #return equivalents(vcfspec, fasta=fasta, fetch_extend_length=fetch_extend_length)[0]
    return equivalents(vcfspec, fetch_extend_length=fetch_extend_length)[0]


def rightmost(vcfspec, fasta, fetch_extend_length=DEFAULT_FETCH_EXTEND_LENGTH):
    #return equivalents(vcfspec, fasta=fasta, fetch_extend_length=fetch_extend_length)[-1]
    return equivalents(vcfspec, fetch_extend_length=fetch_extend_length)[-1]


def check_equivalent(vcfspec1, vcfspec2, fasta):
    return leftmost(vcfspec1, fasta) == leftmost(vcfspec2, fasta)


def equivalents_monoalt(vcfspec, fetch_extend_length=DEFAULT_FETCH_EXTEND_LENGTH):
    vcfspec.check_monoalt(raise_with_false=True)
    vcfspec = make_parsimonious(vcfspec)
    vcfspec = convert_asterisk(vcfspec)
    return shifting_equivalents(vcfspec=vcfspec, fasta=vcfspec.fasta, fetch_extend_length=fetch_extend_length)


def equivalents(vcfspec, fetch_extend_length=DEFAULT_FETCH_EXTEND_LENGTH):
    if vcfspec.check_monoalt(raise_with_false=False):
        return equivalents_monoalt(vcfspec, fetch_extend_length=fetch_extend_length)
    else:
        monoalt_equivalent_list = [
            equivalents_monoalt(monoalt_vcfspec, fetch_extend_length=fetch_extend_length)
            for monoalt_vcfspec in vcfspec.iter_monoalts()
        ]
        result = [
            functools.reduce(merge, vcfspec_list)
            for vcfspec_list in itertools.product(*monoalt_equivalent_list)
        ]
        result.sort(key=(lambda x: x.pos))
        return result


def convert_asterisk(vcfspec):
    if '*' in vcfspec.alts:
        new_alts = tuple(
            '' if x == '*' else x
            for x in vcfspec.alts
        )
        result = vcfspec.spawn(alts=new_alts)
        result = result.extend_left(1)
        result._is_parsimonious = None
        return result
    else:
        return vcfspec


def make_parsimonious(vcfspec):
    #check_vcfspec_monoalt(vcfspec)
    vcfspec.check_monoalt(raise_with_false=True)

    pos = vcfspec.pos
    ref = vcfspec.ref
    alt = vcfspec.alts[0]

    ########
    # left #
    ########
    min_length = min(len(ref), len(alt))
    for left_action_idx in range(0, min_length + 1):
        if left_action_idx == min_length:
            break
        if ref[left_action_idx] != alt[left_action_idx]:
            break
    # Now "left_action_idx" points to the position next to the last 'leading identical portion'
    # and is equal to the length of the 'leading identical portion'

    if left_action_idx == min_length:  # This means the shorter of ref and alt is entirely 'leading identical portion'
        if len(ref) == len(alt):  # Identity
            new_pos = pos
            new_ref = ref[0]
            new_alt = alt[0]
        else:
            new_ref = ref[(left_action_idx - 1):]
            new_alt = alt[(left_action_idx - 1):]
            new_pos = pos + (left_action_idx - 1)
        result = vcfspec.spawn(pos=new_pos, ref=new_ref, alts=(new_alt,))
        result._is_parsimonious = True
        return result

    #########
    # right #
    #########
    min_length_after_lstrip = min_length - left_action_idx
    for right_action_idx in range(-1, -min_length_after_lstrip - 2, -1):
        if right_action_idx == -min_length_after_lstrip - 1:
            break
        if ref[right_action_idx] != alt[right_action_idx]:
            break
    # Now "right_action_idx" points to the position (counting from the right) next to the last 'trailing identical portion'

    if right_action_idx == -min_length_after_lstrip - 1:  
        # This means after stripping leading identical portion, the shorter of ref and alt is entirely 'trailing identical portion'
        # Now len(ref) != len(alt)
        if left_action_idx > 0:
            new_ref = ref[(left_action_idx - 1):(right_action_idx + 1)]
            new_alt = alt[(left_action_idx - 1):(right_action_idx + 1)]
            new_pos = pos + (left_action_idx - 1)
        else:
            added_seq = vcfspec.fasta.fetch(vcfspec.chrom, vcfspec.pos0 - 1, vcfspec.pos0)
            new_ref = added_seq + ref[:(right_action_idx + 1)]
            new_alt = added_seq + alt[:(right_action_idx + 1)]
            new_pos = pos - 1
    else:
        if right_action_idx == -1:
            slice_stop = None
        else:
            slice_stop = right_action_idx + 1
        new_ref = ref[left_action_idx:slice_stop]
        new_alt = alt[left_action_idx:slice_stop]
        new_pos = pos + left_action_idx

    result = vcfspec.spawn(pos=new_pos, ref=new_ref, alts=(new_alt,))
    result._is_parsimonious = True
    return result


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

    def get_result(original_vcfspec, fasta, diff_seq, fetch_extend_length, is_ins):
        # get new_pos_list
        new_pos_list = list()
        navigate(
            original_vcfspec.chrom,
            original_vcfspec.pos,
            fasta,
            diff_seq,
            new_pos_list,
            fetch_extend_length,
            is_ins,
            forward=True,
        )
        navigate(
            original_vcfspec.chrom,
            original_vcfspec.pos,
            fasta,
            diff_seq,
            new_pos_list,
            fetch_extend_length,
            is_ins,
            forward=False,
        )
        # make result
        result = list()
        result.append(original_vcfspec)
        for new_pos, new_ref, new_alt in new_pos_list:
            result.append(original_vcfspec.spawn(pos=new_pos, ref=new_ref, alts=(new_alt,)))
        result.sort(key=(lambda x: x.pos))

        if len(result) == 1:
            result[0].is_leftmost = True
            result[0].is_rightmost = True
        else:
            result[0].is_leftmost = True
            result[0].is_rightmost = False
            result[-1].is_leftmost = False
            result[-1].is_rightmost = True
            if len(result) >= 3:
                for idx in range(1, len(result) - 1):
                    result[idx].is_leftmost = False
                    result[idx].is_rightmost = False

        return result

    # sanity check
    def sanity_check_indel(vcfspec):
        if vcfspec.ref[0] != vcfspec.alts[0][0]:
            raise Exception(f'Inserted/deleted seq must be on the right side, not left side.')

    def result_postprocess(result):
        # set is_leftmost and is_rightmost
        if len(result) == 1:
            result[0].is_leftmost = True
            result[0].is_rightmost = True
        else:
            result[0].is_leftmost = True
            result[0].is_rightmost = False
            result[-1].is_leftmost = False
            result[-1].is_rightmost = True
            if len(result) >= 3:
                for idx in range(1, len(result) - 1):
                    result[idx].is_leftmost = False
                    result[idx].is_rightmost = False
        # set _equivalents
        for x in result:
            x._equivalents = result.copy()

    # main
    #check_vcfspec_monoalt(vcfspec)
    vcfspec.check_monoalt(raise_with_false=True)

    if (
        "n" not in vcfspec.ref and 
        "N" not in vcfspec.ref and 
        "n" not in vcfspec.alts[0] and 
        "N" not in vcfspec.alts[0]
    ):
        # if set('nN').isdisjoint(ref) and set('nN').isdisjoint(alt):
        mttype = vcfspec.get_mutation_type(0)
        if mttype in ("ins", "del"):
            sanity_check_indel(vcfspec)
            is_ins = (mttype == "ins")
            diff_seq = list(vcfspec.alts[0][1:]) if is_ins else list(vcfspec.ref[1:])
            result = get_result(vcfspec, fasta, diff_seq, fetch_extend_length, is_ins)
        else:
            result = [vcfspec]
    else:
        result = [vcfspec]

    result_postprocess(result)

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

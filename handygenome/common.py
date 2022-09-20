"""
Reference genome version aliases

MGSCv37 == mm9
GRCm38 == mm10
GRCm39 == mm39
NCBI36 == hg18
GRCh37 == hg19
GRCh38 == hg38
"""

"""
pysam file mode string

- Valid mode string pattern : ^[rwa]([bzu]?[0-9]?|[0-9]?[bzu]?)$
- Conflicting characters: (b,z,u)

- 'wb' : Compressed BCF (level 6), regardless of file name. 
- 'wb[0-9]' : BCF with level of compression indicated by the number, regardless of file name. 

- 'wz' : Compressed VCF (level 6), regardless of file name.
- 'wz[0-9]' : VCF with level of compression indicated by the number, regardless of file name. 

- 'wu' : Uncompressed VCF, regardless of file name.
- 'wu[0-9]' : Uncompressed VCF, regardless of file name.

- 'w[0-9]' : Uncompressed VCF, regardless of file name. 

- 'w' :
    *.vcf : uncompressed VCF
    *.vcf.gz : compressed VCF (level 6)
    *.bcf : compressed BCF (level 6)
    *.bcf.gz : compressed VCF (level 6)
"""

import sys
import os
import re
import time
import datetime
import tempfile
import collections
import json
import pprint
import shutil
import urllib.request
import urllib.parse
import urllib.error
import io
import contextlib
import inspect
import itertools
import functools
import gzip
import subprocess

import pysam
import pyranges as pr
import numpy as np
import Bio.Seq

import importlib
TOP_PACKAGE_NAME = __name__.split('.')[0]
TOP_PACKAGE = importlib.import_module(TOP_PACKAGE_NAME)
PROJECT_PATH = os.path.dirname(os.path.dirname(TOP_PACKAGE.__file__))
PACKAGE_LOCATION = PROJECT_PATH
DATA_DIR = os.path.join(PROJECT_PATH, 'data')
UTILS_DIR = os.path.join(PROJECT_PATH, 'utils')

DEFAULT_VCFVER = '4.3'

# re patterns
RE_PATS = {
    'int' : re.compile('-?[0-9]+'),
    'float' : re.compile('(-?[0-9]+\.[0-9]+)|(-?[0-9]+(\.[0-9]+)?e-?[0-9]+)'),
    'nucleobases' : re.compile('[ACGTNacgtn]+'),
    'numbered_chromosome' : re.compile('(chr)?[0-9]+'),
    'assembled_chromosome' : re.compile('(chr)?([0-9]+|X|Y)'),

    'alt_bndstring_1' : re.compile('^(?P<t>[^\[\]]+)(?P<bracket1>\[|\])(?P<matechrom>[^:]+):(?P<matepos>[0-9]+)(?P<bracket2>\[|\])$'),
    'alt_bndstring_2' : re.compile('^(?P<bracket1>\[|\])(?P<matechrom>[^:]+):(?P<matepos>[0-9]+)(?P<bracket2>\[|\])(?P<t>[^\[\]]+)$'),
        # these bndstring patterns assume that "t" portion must not be blank
}

# SV symbolic allele strings
SV_ALTS = ('DEL', 'INS', 'DUP', 'INV', 'CNV', 'BND', 'TRA')
CPGMET_ALT = 'CPGMET'

# executable paths
BASH = '/usr/bin/bash'
BWA = os.path.join(UTILS_DIR, 'bwa')
BEDTOOLS = os.path.join(UTILS_DIR, 'bedtools')
GATK = os.path.join(UTILS_DIR, 'gatk')
PERL = '/home/users/pjh/scripts/conda_wrapper/perl'
PYTHON = '/home/users/pjh/tools/miniconda/210821/miniconda3/envs/genome_v5/bin/python'

VEP_V102 = '/home/users/pjh/scripts/conda_wrapper/vep_v102' # for mm10
#VEP_V104 = '/home/users/pjh/scripts/conda_wrapper/vep_v104'
VEP_V105 = '/home/users/pjh/scripts/conda_wrapper/vep_v105'
VEP = VEP_V105
VEP_MM10 = VEP_V102

CONDABIN_PATH = '/home/users/pjh/conda_bin'
SAMTOOLS = os.path.join(CONDABIN_PATH, 'samtools')
BCFTOOLS = os.path.join(CONDABIN_PATH, 'bcftools')
TABIX = os.path.join(UTILS_DIR, 'tabix')
#BEDTOOLS = os.path.join(CONDABIN_PATH, 'bedtools')

# vep cache directory
VEP_CACHE_DIR = '/home/users/pjh/.vep'

# vcf format constants
BCFTOOLS_FORMAT_DICT = { 'v':'w', 'z':'wz', 'u':'wbu', 'b':'wb' }
CYVCF2_FORMAT_DICT = BCFTOOLS_FORMAT_DICT
PYSAM_FORMAT_DICT = { 'v':'wz0', 'z':'wz', 'u':'wb0', 'b':'wb' }
PYSAM_MODE_DICT = PYSAM_FORMAT_DICT
DEFAULT_MODE_BCFTOOLS = 'z'

# http function constants
HTTP_HEADER_POST = {"Content-Type": "application/json", 
                    "Accept": "application/json"}
HTTP_HEADER_GET = {'Content-Type': 'application/json'}


# colors
COLORS = {
    'red':     '\033[38;5;196m',
    'magenta': '\033[38;5;201m',
    'pink':    '\033[38;5;213m',
    'orange':  '\033[38;5;9m',
    'yellow':  '\033[38;5;214m',
    'gold':    '\033[38;5;11m',
    'green':   '\033[38;5;40m',
    'blue':    '\033[38;5;33m',
    'cyan':    '\033[38;5;14m',
    'purple':  '\033[38;5;93m',
    'gray':    '\033[38;5;8m',
    'white':   '\033[38;5;15m',
    'end':     '\033[0m',
    }


class ColorsQQ:
    mine="\033[48;5;6m"
    busy="\033[48;5;244m"
    free="\033[48;5;238m"
    end="\033[0m"
    nor="\033[48;5;160m"
    nor2="\033[48;5;52m"
    danger1="\033[38;5;208m"
    danger2="\033[38;5;196m"
    darkgray="\033[38;5;240m"


def visualize_colors():
    for i in range(256):
        print(f'{i:<3d} \033[38;5;{i}m\\033[38;5;{i}m\033[0m')


def cpformat(obj, **kwargs):
    result = pprint.pformat(obj, **kwargs)
    result = re.sub('(True)', COLORS['green'] + '\\1' + COLORS['end'], result)
    result = re.sub('(False)', COLORS['red'] + '\\1' + COLORS['end'], result)
    result = re.sub('(None)', COLORS['purple'] + '\\1' + COLORS['end'], result)
    return result


def cpprint(obj):
    print(cpformat(obj))


# timer decorator

def deco_timer(func):
    """Print the runtime of the decorated function"""

    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()    # 1
        value = func(*args, **kwargs)
        end_time = time.perf_counter()      # 2
        run_time = end_time - start_time
        print(run_time)

        return value

    return wrapper_timer


###################################################

# argument sanity check decorators

def get_deco_num_set(names, n):
    def decorator(func):
        sig = inspect.signature(func)
        if not set(names).issubset(sig.parameters.keys()):
            raise Exception(
                f'The names of parameters given to '
                f'"get_deco_num_set" function '
                f'is not included in the parameter names of '
                f'the function "{func.__name__}".')

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            n_set = sum((name in ba.arguments) for name in names)
            if n_set != n:
                raise ValueError(
                    f'For the function "{func.__name__}", the '
                    #f'number of parameters set from arguments '
                    f'number of parameters being set, '
                    f'among {tuple(names)}, must be {n}.')

            return func(*args, **kwargs)

        return wrapper

    return decorator


def get_deco_num_set_differently(names, n):
    def decorator(func):
        sig = inspect.signature(func)
        if not set(names).issubset(sig.parameters.keys()):
            raise Exception(
                f'The names of parameters given to '
                f'"get_deco_num_set_differently" function '
                f'is not included in the parameter names of '
                f'the function "{func.__name__}".')

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()

            n_diff = 0
            for name in names:
                if sig.parameters[name].default is None:
                    if sig.parameters[name].default is not ba.arguments[name]:
                        n_diff += 1
                else:
                    if sig.parameters[name].default != ba.arguments[name]:
                        n_diff += 1

            if n_diff != n:
                raise ValueError(
                    f'For the function "{func.__name__}", the '
                    f'number of parameters, among {tuple(names)}, '
                    f'being set as a value different from the default, '
                    f'must be {n}.')

            return func(*args, **kwargs)

        return wrapper

    return decorator


def get_deco_arg_choices(mapping):
    """Args:
        mapping: {'argname': (valid_value1, valid_value2, ...), ...}
    """

    def decorator(func):
        sig = inspect.signature(func)
        if not set(mapping.keys()).issubset(sig.parameters.keys()):
            raise Exception(
                f'The names of parameters given to '
                f'"get_deco_check_arg_choices" function '
                f'is not included in the parameter names of '
                f'the function "{func.__name__}".')

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()
            for key, val in mapping.items():
                if ba.arguments[key] not in val:
                    raise ValueError(
                        f'For the function "{func.__name__}", '
                        f'the parameter "{key}" must be one of these values: '
                        f'{tuple(val)}.')

            return func(*args, **kwargs)

        return wrapper

    return decorator

# other decorators

#def get_deco_noti(title):
#    def decorator(func):
#        def wrapper(*args, **kwargs):
#            ts = get_timestamp()
#            print_timestamp(f'[{ts}] BEGINNING {title}')
#
#            result = func(*args, **kwargs)
#
#            ts = get_timestamp()
#            print_timestamp(f'[{ts}] FINISHED')
#
#            return result
#
#        return wrapper
#
#    return decorator


def get_deco_timestamp(msg, logger):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            logger.info(f'BEGINNING {msg}')
            result = func(*args, **kwargs)
            logger.info(f'FINISHED {msg}')
            return result
        return wrapper

    return decorator

###################################################


class RefverDict(collections.UserDict):
    standards = ('MGSCv37', 'GRCm38', 'GRCm39', 'NCBI36', 'GRCh37', 'GRCh38')
    aliases = {
        'NCBI36': ('hg18', 'ncbi36'),
        'GRCh37': ('hg19', 'grch37'),
        'GRCh37_hs37d5': tuple(),
        'GRCh38': ('hg38', 'grch38'),
        'MGSCv37': ('mm9',),
        'GRCm38': ('mm10', 'grcm38'),
        'GRCm39': ('mm39', 'grcm39'),
        }

    converter = dict()
    for refver, alias_refvers in aliases.items():
        converter[refver] = refver
        for alias_refver in alias_refvers:
            converter[alias_refver] = refver

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not set(self.keys()).issubset(self.__class__.aliases.keys()):
            raise Exception(
                f'RefverDict construction keys must be restricted to: '
                f'{tuple(self.__class__.aliases.keys())}')

        for key in tuple(self.keys()):
            if key in self.__class__.aliases:
                for new_key in self.__class__.aliases[key]:
                    self[new_key] = self[key]


# chr1 lengths
CHR1_LENGTHS = {
    'MGSCv37': 197_195_432,
    'GRCm38': 195_471_971,
    'GRCm39': 195_154_279,
    'NCBI36': 247_249_719,
    'GRCh37': 249_250_621,
    'GRCh38': 248_956_422,
    }
CHR1_LENGTHS_REV = {val: key for key, val in CHR1_LENGTHS.items()}

# default fasta paths
FASTA_URLS = {
    'GRCh37_1000genomes': 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz',
    'GRCh37_ucsc': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz',
    'GRCh38_1000genomes': 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa',
    }

DEFAULT_FASTA_PATHS = RefverDict({
    'NCBI36': '/home/users/pjh/References/reference_genome/NCBI36/ucsc/custom_files/hg18.fa',
    'GRCh37': '/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta',
    'GRCh37_hs37d5': '/home/users/sypark/02_Reference/15_pcawg/genome.fa',
    'GRCh38': '/home/users/data/01_reference/human_g1k_v38/Homo_sapiens_assembly38.fasta',

    'MGSCv37': '/home/users/pjh/References/reference_genome/mm9/ucsc/custom_files/mm9.fa',
    'GRCm38': '/home/users/pjh/References/reference_genome/GRCm38/ucsc/custom_files/mm10.fa',
    'GRCm39': '/home/users/pjh/References/reference_genome/GRCm39/ucsc/custom_files/mm39.fa',
    })

AVAILABLE_REFVERS = tuple(DEFAULT_FASTA_PATHS.keys())
AVAILABLE_REFVERS_PLUSNONE = AVAILABLE_REFVERS + (None,)


###################################################


class ChromDict(collections.OrderedDict):
    @get_deco_num_set_differently(
        ('fasta_path', 'fasta', 'bam_path', 'bam', 
         'vcfheader', 'bamheader', 'custom', 'refver'), 1)
    @get_deco_arg_choices({'refver': AVAILABLE_REFVERS_PLUSNONE})
    def __init__(self, fasta_path=None, fasta=None, bam_path=None, bam=None, 
                 vcfheader=None, bamheader=None, custom=None, refver=None):
        """
        Args:
            fasta: pysam.FastaFile object
            bam: pysam.AlignmentFile object
            vcfheader: pysam.VariantHeader object
            bamheader: pysam.AlignmentHeader object
            custom: {'contigs': ['contig1', 'contig2', ...], 
                     'lengths': [length1, length2, ...] }
        """

        # set self dict
        if vcfheader is not None:
            for contig in vcfheader.contigs.values():
                self[contig.name] = contig.length
        elif custom is not None:
            for chrom, length in zip(custom['contigs'], custom['lengths']):
                self[chrom] = length
        else:
            if fasta_path is not None:
                wrapper = pysam.FastaFile(fasta_path)
            elif fasta is not None:
                wrapper = fasta
            elif bam_path is not None:
                wrapper = pysam.AlignmentFile(bam_path)
            elif bam is not None:
                wrapper = bam
            elif bamheader is not None:
                wrapper = bamheader
            elif refver is not None:
                wrapper = pysam.FastaFile(DEFAULT_FASTA_PATHS[refver])
    
            for chrom, length in zip(wrapper.references, wrapper.lengths):
                self[chrom] = length
    
            if any(x is not None 
                   for x in (fasta_path, bam_path, refver)):
                wrapper.close()

        # set contigs, lengths
        self.contigs = list(self.keys())
        self.lengths = list(self.values())

    def to_gr(self):
        return pr.PyRanges(
            chromosomes=self.contigs,
            starts=([0] * len(self.contigs)),
            ends=self.lengths
        )

    def to_interval_list(self):
        intvlist = IntervalList()
        for contig, length in zip(self.contigs, self.lengths):
            interval = Interval(contig, start0=0, end0=length)
            intvlist.append(interval)

        return intvlist


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
        return (f'<Vcfspec ({self.chrom}:{self.pos} '
                f'{self.ref}>{altstring})>')

    def __hash__(self):
        return hash(self.get_tuple())

    def __eq__(self, other):
        return all(getattr(self, key) == getattr(other, key)
                   for key in ('chrom', 'pos', 'ref', 'alts'))

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

# alias
check_vcfspec_monoallele = check_vcfspec_monoalt


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

    def __init__(self, chrom, start1=None, end1=None, start0=None, end0=None, is_reverse=False):
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

    def to_gr(self):
        return pr.from_dict({'Chromosome': [self.chrom],
                             'Start': [self.start0],
                             'End': [self.end0]})

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
    def get_depth_bins(cls, refver, width=100_000):
        result = cls()
        chromdict = ChromDict(refver=refver)
        for contig, length in zip(chromdict.contigs, chromdict.lengths):
            intvlist_contig = cls()
            intvlist_contig.append(Interval(chrom=contig, start0=0, end0=length))
            for item in intvlist_contig.split(width=width):
                result.extend(item)
        return result
    ################

    def write_bed(self, outfile_path):
        with openfile(outfile_path, 'w') as outfile:
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
        return list(itertools.accumulate(intv.length 
                                         for intv in self))
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
            return coord_sortkey(intv.chrom, intv.start1, chromdict)
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

    @get_deco_num_set(('b', 'l', 'r'), 1)
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

    @get_deco_num_set(('num', 'width'), 1)
    def split(self, num=None, width=None):
        #self.sort_intervals(chromdict)

        # get result_lengths_cumsum
        total_length = self.lengths_cumsum[-1]
        if num is not None:
            result_lengths_cumsum = list(itertools.accumulate(
                get_interval_lengths_num(total_length, num)))
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


###################################################

def without_N(seq):
    return 'n' not in seq and 'N' not in seq


def zrange(n):
    width = len(str(n-1))
    for idx in range(n):
        zidx = str(idx).zfill(width)
        yield zidx


def zenumerate(iterable):
    iterable_tup = tuple(iterable)
    length = len(iterable_tup)
    width = len(str(length-1))
    for idx, item in enumerate(iterable_tup):
        zidx = str(idx).zfill(width)
        yield zidx, item


def round_sig(num, n):
    """round 'num' with 'n' significant digits"""

    assert n > 0
    from math import log10, floor
    return round(num, -floor(log10(abs(num))) + (n-1))


def get_interval_lengths_num(total_length, num):
    interval_num = min(num, total_length)
    q, r = divmod(total_length, interval_num)
    interval_width_1 = q + 1
    interval_num_1 = r
    interval_width_2 = q
    interval_num_2 = interval_num - r

    result = list(itertools.chain(
        itertools.repeat(interval_width_1, interval_num_1),
        itertools.repeat(interval_width_2, interval_num_2)))

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


def str_to_nonstr(val):
    #assert isinstance(val, str)

    if val.lower() in ('none', 'null'):
        return None
    elif val.lower() == 'nan':
        return np.nan
    elif val.lower() == 'true':
        return True
    elif val.lower() == 'false':
        return False
    elif RE_PATS['int'].fullmatch(val) is not None:
        return int(val)
    elif RE_PATS['float'].fullmatch(val) is not None:
        return float(val)
    else:
        return val


def nonstr_to_str(val):
    assert not isinstance(val, str)

    if val is None:
        return None
    else:
        return str(val)


def get_datestring():
    """
    Returns a string like '2021-12-06 11:49:55'
    """

    return str( datetime.datetime.now() ).split('.')[0]


def get_timestamp():
    """
    Returns a string like 'KST 2021-12-06 11:51:36'
    """

    dt = datetime.datetime.now().astimezone()
    timestamp = f'{str(dt.tzinfo)} {str(dt).split(".")[0]}'

    return timestamp


def print_err(*args, stderr = True, files = None, **kwargs):
    """
    Args:
        stderr: (Bool) Whether to write to stderr
        files: A list of file paths to which message is written
    """

    if stderr:
        print(*args, file = sys.stderr, flush = True, **kwargs)

    if files is not None:
        for fname in files:
            with open(fname, 'a') as f:
                print(*args, file = f, flush = True, **kwargs)


def printerr(*args, **kwargs):
    print_err(*args, **kwargs)


def printerrdate(*args, **kwargs):
    datestring = get_datestring()
    print_err(f'[{datestring}]', *args, **kwargs)


def print_timestamp(*args, **kwargs):
    timestamp = get_timestamp()
    print_err(f'[{timestamp}]', *args, **kwargs)


def coord_sortkey(chrom, pos, chromdict): 
    """
    Args:
        pos: 1-based
        chromdict: ChromDict class instance
    """

    return (chromdict.contigs.index(chrom), pos)


def get_read_sortkey(chromdict):
    def sortkey(read):
        return (chromdict.contigs.index(read.reference_name), read.reference_start)

    return sortkey


def get_vr_sortkey(chromdict):
    def sortkey(vr):
        return (chromdict.contigs.index(vr.contig), vr.pos, vr.ref) + vr.alts
        #return coord_sortkey(vr.contig, vr.pos, chromdict)

    return sortkey


def get_vcfspec_sortkey(chromdict):
    def sortkey(vcfspec):
        return (chromdict.contigs.index(vcfspec.chrom), vcfspec.pos, vcfspec.ref) + vcfspec.alts
        #return coord_sortkey(vcfspec.chrom, vcfspec.pos, chromdict)

    return sortkey


# deprecated
def get_vcfspec_order(vcfspec1, vcfspec2, chromdict):
    """
    Returns:
        0: equal
        negative integer: vcfspec1 comes first
        positive integer: vcfspec2 comes first
    """

    return (coord_sortkey(vcfspec1.chrom, vcfspec1.pos, chromdict) 
            - coord_sortkey(vcfspec2.chrom, vcfspec2.pos, chromdict))


def compare_coords(chrom1, pos1, chrom2, pos2, chromdict):
    """
    Returns:
        0: equal; -1: chrom1/pos1 comes first; 1: chrom2/pos2 comes first
    """

    if chrom1 == chrom2 and pos1 == pos2:
        return 0
    else:
        if (coord_sortkey(chrom1, pos1, chromdict) 
            < coord_sortkey(chrom2, pos2, chromdict)):
            return -1
        else:
            return 1

###################################################


def rm_newline(line):
    return re.sub('(\r)?\n$', '', line)


def rm_newline_byte(byteline):
    return re.sub(b'(\r)?\n$', b'', byteline)


def get_linesp(line, sep='\t'):
    return rm_newline(line).split(sep)


def get_linesp_byte(byteline, sep=b'\t'):
    return rm_newline_byte(byteline).split(sep)


###################################################


def get_mttype(ref, alt):
    if RE_PATS['nucleobases'].fullmatch(alt) is None:
        if any(
                (re.fullmatch(f'<{x}(:.+)?>', alt) is not None)
                for x in SV_ALTS):
            mttype = 'sv'
        elif (
                (RE_PATS['alt_bndstring_1'].fullmatch(alt) is not None) or 
                (RE_PATS['alt_bndstring_2'].fullmatch(alt) is not None)):
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


def get_indelseq(ref, alt):
    mttype = get_mttype(ref, alt)
    if mttype == 'ins':
        indelseq = alt[1:]
    elif mttype == 'del':
        indelseq = ref[1:]
    else:
        indelseq = None
    
    return indelseq


###################################################

def listdir(path):
    return sorted(os.path.join(path, x) for x in os.listdir(path))


def get_padded_indices(n):
    """Begins with 0"""

    width = len(str(n-1))
    result = [str(idx).zfill(width) for idx in range(n)]

    return result


###################################################


def printwidth_get_width_list(df):
    '''
    df: [
    [line1_field1, line1_field2, ... ],
    [line2_field1, line2_field2, ... ],
    ...,
    ]
    '''
    width_list = list()
    for i in range(len(df[0])):
        width_list.append(list())

    for line in df:
        for idx, field in enumerate(line):
            width_list[idx].append(len(field))

    for idx, e in enumerate(width_list):
        width_list[idx] = max(e)

    return width_list


def printwidth_print_line(line, width_list, margin, target):
    printresult = ''
    for idx, e in enumerate(line):
        printresult += f'{e:>{width_list[idx] + margin}s}'
    if target == 'out':
        print(printresult, flush = True)
    elif target == 'err':
        print(printresult, flush = True, file = sys.stderr)


def printwidth(df, margin = 2, target = 'out'):
    for line in df:
        for idx, e in enumerate(line):
            line[idx] = str(e)

    width_list = printwidth_get_width_list(df)
    for line in df:
        printwidth_print_line(line, width_list, margin, target)


###################################################

@get_deco_num_set(('chromdict', 'vcfheader', 'bamheader'), 1)
def infer_refver(chromdict=None, vcfheader=None, bamheader=None):
    if chromdict is not None:
        return infer_refver_chromdict(chromdict)
    elif vcfheader is not None:
        return infer_refver_vcfheader(vcfheader)
    elif bamheader is not None:
        return infer_refver_bamheader(bamheader)


def infer_refver_chromdict(chromdict):
    if '1' in chromdict.contigs:
        chr1_length = chromdict.lengths[chromdict.contigs.index('1')]
    elif 'chr1' in chromdict.contigs:
        chr1_length = chromdict.lengths[chromdict.contigs.index('chr1')]
    else:
        raise Exception(
            f'"1" and "chr1" both absent from the chromosome name list.')
    
    if chr1_length in CHR1_LENGTHS_REV:
        refver = CHR1_LENGTHS_REV[chr1_length]
    else:
        raise Exception(f'Cannot infer refver: unknown chr1 length')
        #refver = None # unknown reference genome
    
    return refver


def infer_refver_vcfheader(vcfheader):
    return infer_refver_chromdict(ChromDict(vcfheader=vcfheader))


def infer_refver_bamheader(bamheader):
    return infer_refver_chromdict(ChromDict(bamheader=bamheader))


def infer_refver_vr(vr):
    return infer_refver_chromdict(ChromDict(vcfheader=vr.header))


###################################################


def http_get(url, params=None, headers=dict(), text=False):
    if params is not None:
        url = url + '?' + urllib.parse.urlencode(params)

    req = urllib.request.Request(url, headers=headers, method='GET')
    try:
        with urllib.request.urlopen(req) as response:
            if text:
                result = str(response.read(), 'utf-8')
            else:
                result = json.loads(response.read())
    except urllib.error.HTTPError as e:
        print(str(e.read(), 'utf-8'))
        raise

    return result


def http_post(url, data, params=None, headers=dict(), text=False):
    data = json.dumps(data).encode('ascii')

    if params is not None:
        url = url + '?' + urllib.parse.urlencode(params)

    req = urllib.request.Request(url, data=data, headers=headers, 
                                 method='POST')
    try:
        with urllib.request.urlopen(req) as response:
            if text:
                result = str(response.read(), 'utf-8')
            else:
                result = json.loads(response.read())
    except urllib.error.HTTPError as e:
        print(str(e.read(), 'utf-8'))
        raise

    return result


def download(url, path):
    with urllib.request.urlopen(url) as response:
        with open(path, 'wb') as outfile:
            shutil.copyfileobj(response, outfile)


def unzip(src, dest, rm_src=False):
    with gzip.open(src, 'rb') as infile:
        with open(dest, 'wb') as outfile:
            shutil.copyfileobj(infile, outfile)      
    if rm_src:
        os.remove(src)


@get_deco_arg_choices({'mode': ('r', 'w')})
def openfile(fname, mode='r'):
    if mode == 'r':
        mode = 'rt'
    elif mode == 'w':
        mode = 'wt'

    if fname.endswith('.gz'):
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)

###################################################

def write_mode_arghandler(mode_bcftools, mode_pysam):
    """
    Args:
        bcftools_mode: bcftools style mode string
        mode_pysam: pysam style mode string

    If 'mode_pysam' is set as a non-None value,
    it overrides 'mode_bcftools'.
    """

    if mode_pysam is None:
        return PYSAM_MODE_DICT[mode_bcftools]
    else:
        return mode_pysam


def fileiter(path, sep='\t', remove_leading_hashes=True):
    infile = openfile(path, 'r')
    headerline = next(infile)
    if remove_leading_hashes:
        headerline = re.sub('^#*', '', headerline)
    header = get_linesp(headerline, sep=sep)

    for line in infile:
        linesp = get_linesp(line, sep=sep)
        if len(linesp) != len(header):
            raise Exception(
                f'Field numbers of the header line and the current '
                f'line are different:\n'
                f'header: {header}\n'
                f'line: {linesp}')
        linedict = dict(zip(header, linesp))
        yield linedict

    infile.close()


def iter_lineno_logging(line_iterator, logger, logging_lineno, msgfunc=None):
    if msgfunc is None:
        def msgfunc(NR):
            return f'Processing {NR:,}th line'

    NR = 0
    for line in line_iterator:
        NR += 1
        if NR % logging_lineno == 0:
            logger.info(msgfunc(NR))
        yield line


def gr_iterrows(gr):
    """Args:
        gr: pyranges.PyRanges object
    """
    indexers = (np.eye(gr.df.shape[0]) == 1)
    for row_indexer in indexers:
        yield gr[row_indexer]


###################################################


def to_pyrimidine(base):
    if base in 'AG':
        return Bio.Seq.reverse_complement(base)
    else:
        return base


def get_different_base(base):
    assert base in 'ACTG'
    if base == 'A':
        return 'C'
    else:
        return 'A'


###################################################

# deprecated; this does not work
def get_vcf_noerr(*args, **kwargs):
    with contextlib.redirect_stderr(io.StringIO()) as err, \
            contextlib.redirect_stdout(io.StringIO()) as out:
        vcf = pysam.VariantFile(*args, **kwargs)

    for buf in (err, out):
        msg = buf.getvalue()
        if not msg.startswith('[E::idx_find_and_load] '
                              'Could not retrieve index file for'):
            print(msg, end='', file=sys.stderr)

    return vcf
        

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
import signal
import logging
import importlib
import operator
import ftplib

import psutil
import pysam
import pandas as pd
import pyranges as pr
import numpy as np
import Bio.Seq
import scipy.stats

import handygenome.deco as deco

TOP_PACKAGE_NAME = __name__.split('.')[0]
TOP_PACKAGE = importlib.import_module(TOP_PACKAGE_NAME)
PROJECT_PATH = os.path.dirname(os.path.dirname(TOP_PACKAGE.__file__))
PACKAGE_LOCATION = PROJECT_PATH
DATA_DIR = os.path.join(PROJECT_PATH, 'data')
UTILS_DIR = os.path.join(PROJECT_PATH, 'externals')
R_DIR = os.path.join(PROJECT_PATH, 'R')

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
#SV_ALTS = ('DEL', 'INS', 'DUP', 'INV', 'CNV', 'BND', 'TRA')
#CPGMET_ALT = 'CPGMET'

# executable paths
BASH = '/usr/bin/bash'
BWA = os.path.join(UTILS_DIR, 'bwa')
BEDTOOLS = os.path.join(UTILS_DIR, 'bedtools')
GATK = os.path.join(UTILS_DIR, 'gatk_wrapper.sh')

PERL = '/home/users/pjh/scripts/conda_wrapper/perl'
PYTHON = '/home/users/pjh/tools/miniconda/221104/miniconda3/envs/genome_v7/bin/python'

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


###################################################


class RefverDict_withAllKeys(collections.UserDict):
    standards = ('MGSCv37', 'GRCm38', 'GRCm39', 'NCBI36', 'GRCh37', 'GRCh38', 'banana')
    aliases = {
        'NCBI36': ('hg18', 'ncbi36'),
        'GRCh37': ('hg19', 'grch37'),
        'GRCh37_hs37d5': tuple(),
        'GRCh38': ('hg38', 'grch38'),
        'MGSCv37': ('mm9',),
        'GRCm38': ('mm10', 'grcm38'),
        'GRCm39': ('mm39', 'grcm39'),
        'banana': tuple(),
        }
    known_refvers = tuple(aliases.keys()) + tuple(itertools.chain.from_iterable(aliases.values()))

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


class RefverDict(collections.UserDict):
    standards_byspecies = {
        'Mus_musculus': (
            'MGSCv3', 
            'MGSCv34', 
            'MGSCv35', 
            'MGSCv36', 
            'MGSCv37', 
            'GRCm38', 
            'GRCm39',
        ),
        'Homo_sapiens': (
            'NCBI33', 
            'NCBI34', 
            'NCBI35', 
            'NCBI36', 
            'GRCh37', 
            'GRCh37_hs37d5', 
            'GRCh38', 
            'T2T-CHM13v2',
        ),
        'Musa_acuminata': (
            'ASM31385v1',
            'ASM31385v2',
        ),
    }
    standards = tuple(itertools.chain.from_iterable(standards_byspecies.values()))
#    (
#        # mouse
#        'MGSCv3', 'MGSCv34', 'MGSCv35', 'MGSCv36', 'MGSCv37', 'GRCm38', 'GRCm39', 
#        # human
#        'NCBI33', 'NCBI34', 'NCBI35', 'NCBI36', 'GRCh37', 'GRCh37_hs37d5', 'GRCh38', 'T2T-CHM13v2',
#        # others
#        'banana',
#    )
    aliases = {
        'NCBI36': ('hg18', 'ncbi36'),
        'GRCh37': ('hg19', 'grch37'),
        'GRCh38': ('hg38', 'grch38'),

        'MGSCv37': ('mm9',),
        'GRCm38': ('mm10', 'grcm38'),
        'GRCm39': ('mm39', 'grcm39'),

        'ASM31385v2': ('banana',),
    }
    assert set(aliases.keys()).issubset(standards)
    known_refvers = (
        tuple(standards)
        + tuple(itertools.chain.from_iterable(aliases.values()))
    )

    @classmethod
    def standardize(cls, refver):
        if refver in cls.standards:
            return refver
        else:
            for standard, alias_list in cls.aliases.items():
                if refver in alias_list:
                    return standard
            raise Exception(f'Input reference version string ({refver}) is unknown one.')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not set(self.keys()).issubset(self.__class__.standards):
            raise Exception(
                f'RefverDict construction keys must be restricted to: '
                f'{self.__class__.standards}'
            )

    def __getitem__(self, key):
        key = self.__class__.standardize(key)
        try:
            result = super().__getitem__(key)
        except KeyError:
            raise Exception(f'Input reference version is not available.')

        return result

    def __contains__(self, key):
        key = self.__class__.standardize(key)
        return super().__contains__(key)

    def get_valid_keys(self):
        result = list()
        for key in self.keys():
            result.append(key)
            if key in self.__class__.aliases.keys():
                result.extend(self.__class__.aliases[key])
        return result


# chr1 lengths
CHR1_LENGTHS = RefverDict({
    'MGSCv37': 197_195_432,
    'GRCm38': 195_471_971,
    'GRCm39': 195_154_279,
    'NCBI36': 247_249_719,
    'GRCh37': 249_250_621,
    'GRCh38': 248_956_422,
    'ASM31385v2': 41_765_374,
})
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

    #'banana': '/home/users/yeonjin/practice/banana/reference/Musa_acuminata_pahang_v4_cp.fasta',
    })

DEFAULT_FASTAS = RefverDict({refver: pysam.FastaFile(path) for refver, path in DEFAULT_FASTA_PATHS.items()})


AVAILABLE_REFVERS = tuple(DEFAULT_FASTA_PATHS.keys())
AVAILABLE_REFVERS_PLUSNONE = AVAILABLE_REFVERS + (None,)


####################################################


class ChromDict(collections.OrderedDict):
    @deco.get_deco_num_set_differently(
        ('fasta_path', 'fasta', 'bam_path', 'bam', 
         'vcfheader', 'bamheader', 'custom', 'refver'), 1)
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
                wrapper = DEFAULT_FASTAS[refver]
    
            for chrom, length in zip(wrapper.references, wrapper.lengths):
                self[chrom] = length
    
            if any((x is not None) for x in (fasta_path, bam_path)):
                wrapper.close()

        # set contigs, lengths
        self.contigs = list(self.keys())
        self.lengths = list(self.values())

    def get_cumpos0(self, chrom, pos0):
        chromidx = self.contigs.index(chrom)
        if chromidx == 0:
            return pos0
        else:
            return pos0 + sum(self.lengths[:chromidx])

    @functools.cached_property
    def is_chr_prefixed(self):
        relevant_chroms = [
            x for x in self.contigs 
            if RE_PATS['assembled_chromosome'].fullmatch(x) is not None
        ]
        startswith_chr = [x.startswith('chr') for x in relevant_chroms]
        if all(startswith_chr):
            return True
        elif not any(startswith_chr):
            return False
        else:
            raise Exception(f'Chromosome names are inconsistent of whether prefixed with "chr"')

    @functools.cached_property
    def XY_names(self):
        def helper(name, xy):
            assert xy in ('X', 'Y')
            if name in self.contigs:
                return name
            else:
                raise Exception(f'No {xy} chromosome name detected')
            
        if self.is_chr_prefixed:
            X = helper('chrX', 'X')
            Y = helper('chrY', 'Y')
        else:
            X = helper('X', 'X')
            Y = helper('Y', 'Y')

        return X, Y

    @functools.cached_property
    def assembled_chroms(self):
        return [x for x in self.contigs if RE_PATS['assembled_chromosome'].fullmatch(x) is not None]

    def to_gr(self, assembled_only=True, as_gr=True):
        result = pd.DataFrame({
            'Chromosome': self.contigs,
            'Start': 0,
            'End': self.lengths,
        })
        if assembled_only:
            selector = result['Chromosome'].apply(
                lambda x: RE_PATS['assembled_chromosome'].fullmatch(x) is not None
            )
            result = result.loc[selector, :]

        if as_gr:
            return pr.PyRanges(result)
        else:
            return result

    def to_interval_list(self):
        intvlist = IntervalList()
        for contig, length in zip(self.contigs, self.lengths):
            interval = Interval(contig, start0=0, end0=length)
            intvlist.append(interval)

        return intvlist

    def get_chrom_indexes(self, chroms):
        if np.isscalar(chroms):
            return self.contigs.index(chroms)
        else:
            where = np.where(
                np.array(chroms)[:, np.newaxis] == np.array(self.contigs)
            )
            assert len(set(where[0])) == len(chroms), f'Unknown chromosome names are included.'
            result = [
                tuple(subiter)[0][1] for key, subiter in itertools.groupby(
                    zip(where[0], where[1]), key=operator.itemgetter(0)
                )
            ]
            return result

    def get_chrompos_sortkey(self, chroms, start0s=None, end0s=None):
        assert not ((start0s is None) and (end0s is not None))

        chrom_indexes = self.get_chrom_indexes(chroms)
        if start0s is None:
            sortkey = np.lexsort([chrom_indexes])
        else:
            if end0s is None:
                sortkey = np.lexsort([start0s, chrom_indexes])
            else:
                sortkey = np.lexsort([end0s, start0s, chrom_indexes])

        return sortkey

    def sort_chrompos(self, chroms, start0s=None, end0s=None):
        assert not ((start0s is None) and (end0s is not None))

        # arg handling
        chroms = np.array(chroms)
        if start0s is not None:
            start0s = np.array(start0s)
        if end0s is not None:
            end0s = np.array(end0s)

        # result
        sortkey = self.get_chrompos_sortkey(chroms, start0s, end0s)
        if start0s is None:
            return chroms[sortkey]
        else:
            if end0s is None:
                return chroms[sortkey], start0s[sortkey]
            else:
                return chroms[sortkey], start0s[sortkey], end0s[sortkey]


DEFAULT_CHROMDICTS = RefverDict({
    key: ChromDict(fasta=val) for key, val in DEFAULT_FASTAS.items()
    if key != 'GRCh37_hs37d5'
})


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
    def from_margin(cls, refver, chrom_left, start0_left, chrom_right, end0_right):
        chromdict = DEFAULT_CHROMDICTS[refver]
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

# logger
def make_funclogger(level, name):
    level = getattr(logging, level.upper())
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False
    
    formatter = logging.Formatter(
        fmt='[%(asctime)s %(levelname)s] %(module)s.%(funcName): %(message)s', 
        datefmt='%Z %Y-%m-%d %H:%M:%S',
    )

    sh = logging.StreamHandler()
    sh.setLevel(level)
    sh.setFormatter(formatter)
    logger.addHandler(sh)

    return logger


FUNCLOGGER_DEBUG = make_funclogger(level='debug', name='FUNCLOGGER_DEBUG')
FUNCLOGGER_INFO = make_funclogger(level='info', name='FUNCLOGGER_INFO')



###################################################


def repr_base(obj, keylist, comma_sep_int=False):
    string_list = list()
    for key in keylist:
        val = getattr(obj, key)
        if isinstance(val, int) and comma_sep_int:
            string_list.append(f'{key}={val:,}')
        else:
            string_list.append(f'{key}={repr(val)}')

    return ', '.join(string_list)



###################################################

def without_N(seq):
    return 'n' not in seq and 'N' not in seq


def zrange(n):
    assert n > 0

    width = len(str(n - 1))
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


def zidx_to_idx(zidx):
    if re.match('^0+$', zidx):
        return 0
    else:
        return int(re.sub(r'^0*', '', zidx))


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

    result = list(
        itertools.chain(
            itertools.repeat(interval_width_1, interval_num_1),
            itertools.repeat(interval_width_2, interval_num_2),
        )
    )

    return result


def get_split_nums(total_length, num):
    arrsp = np.array_split(np.zeros(total_length), num)
    return np.fromiter(
        (x.shape[0] for x in arrsp if x.shape[0] != 0), 
        dtype=int,
    )


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
    """Returns a string like 'KST 2021-12-06 11:51:36'"""
    dt = datetime.datetime.now().astimezone()
    return f'{str(dt.tzinfo)} {str(dt).split(".")[0]}'


def print_err(*args, stderr=True, files=None, **kwargs):
    """
    Args:
        stderr: (Bool) Whether to write to stderr
        files: A list of file paths to which message is written
    """

    if stderr:
        print(*args, file=sys.stderr, flush=True, **kwargs)

    if files is not None:
        for fname in files:
            with open(fname, 'a') as f:
                print(*args, file=f, flush=True, **kwargs)


def printerr(*args, **kwargs):
    print_err(*args, **kwargs)


def printerrdate(*args, **kwargs):
    datestring = get_datestring()
    print_err(f'[{datestring}]', *args, **kwargs)


def print_timestamp(*args, **kwargs):
    print_err(f'[{get_timestamp()}]', *args, **kwargs)


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



#def get_indelseq(ref, alt):
#    mttype = get_mttype(ref, alt)
#    if mttype == 'ins':
#        indelseq = alt[1:]
#    elif mttype == 'del':
#        indelseq = ref[1:]
#    else:
#        indelseq = None
#    
#    return indelseq


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

@deco.get_deco_num_set(('chromdict', 'vcfheader', 'bamheader'), 1)
def infer_refver(chromdict=None, vcfheader=None, bamheader=None):
    if chromdict is not None:
        return infer_refver_chromdict(chromdict)
    elif vcfheader is not None:
        return infer_refver_vcfheader(vcfheader)
    elif bamheader is not None:
        return infer_refver_bamheader(bamheader)


def infer_refver_chromdict(chromdict):
    return infer_refver_base(chromdict.contigs, chromdict.lengths)


def infer_refver_vcfheader(vcfheader):
    contigs = list()
    lengths = list()
    for contig in vcfheader.contigs.values():
        contigs.append(contig.name)
        lengths.append(contig.length)
    return infer_refver_base(contigs, lengths)


def infer_refver_vr(vr):
    return infer_refver_vcfheader(vr.header)


def infer_refver_fasta(fasta):
    return infer_refver_pysamwrapper(fasta)


def infer_refver_bamheader(bamheader):
    return infer_refver_pysamwrapper(bamheader)


def infer_refver_pysamwrapper(wrapper):
    return infer_refver_base(wrapper.references, wrapper.lengths)


def infer_refver_base(contigs, lengths):
    chr1_names_candidates = {'1', 'chr1', 'chr01'}
    intersection = chr1_names_candidates.intersection(contigs)
    if len(intersection) == 0:
        raise Exception(f'There is no chromosome name which looks like "chr1"')
    elif len(intersection) > 1:
        raise Exception(f'There is more than one  chromosome names which looks like "chr1"')
    chr1_name = intersection.pop()
    chr1_length = lengths[contigs.index(chr1_name)]
    
    if chr1_length in CHR1_LENGTHS_REV:
        refver = CHR1_LENGTHS_REV[chr1_length]
    else:
        raise Exception(f'Cannot infer refver: unknown chr1 length')
    
    return refver


###################################################

# network functionalities

def retry_or_raise(n_try, retry_count, exc, retry_interval):
    if n_try > retry_count:
        print(exc.read().decode('utf-8'))
        raise Exception(f'Exceeded maximum retry count({retry_count}).') from exc
    else:
        if isinstance(exc, TimeoutError):
            print_timestamp(f'Retrying due to TimeoutError (Exception: {exc}); n_try={n_try}')
        else:
            print_timestamp(f'Retrying; n_try={n_try}')
        time.sleep(retry_interval)


# ftp

def ftp_login(url, retry_count=10, retry_interval=1, timeout=5):
    n_try = 0
    while True:
        n_try += 1
        try:
            with contextlib.redirect_stdout('/dev/null'):
                ftp = ftplib.FTP(url, timeout=timeout)
                ftp.login()
        except TimeoutError as exc:
            retry_or_raise(n_try, retry_count, exc, retry_interval)
            continue
        except OSError as exc:
            if (str(exc) == '[Errno 101] Network is unreachable'):
                retry_or_raise(n_try, retry_count, exc, retry_interval)
                continue
            else:
                raise
        else:
            break

    return ftp


def trim_path(path):
    path = re.sub('/+$', '', path)
    path = re.sub('/{2,}', '/', path)
    return path


def ftp_listdir(ftp, path):
    path = trim_path(path)
    fname_list = list()
    ftp.cwd(path)
    ftp.retrlines('NLST', callback=(lambda x: fname_list.append(path + '/' + x)))
    return fname_list


# http

def http_run_urlopen(url_or_req, retry_count=10, retry_interval=1, urlopen_timeout=5):
    url_string = (
        url_or_req
        if isinstance(url_or_req, str) else
        url_or_req.full_url
    )
    print_timestamp(f'Trying to open url {repr(url_string)}')

    n_try = 0
    while True:
        n_try += 1
        try:
            response = urllib.request.urlopen(url_or_req, timeout=urlopen_timeout)
        except urllib.error.URLError as exc:
            if isinstance(exc, urllib.error.HTTPError):
                if exc.code == 500:  # internal server error
                    retry_or_raise(n_try, retry_count, exc, retry_interval)
                    continue
                else:
                    print(exc.read().decode('utf-8'))
                    raise
            else: 
                if str(exc) == '<urlopen error [Errno 101] Network is unreachable>':
                    retry_or_raise(n_try, retry_count, exc, retry_interval)
                    continue
                else:
                    print(exc.read().decode('utf-8'))
                    raise
        except TimeoutError as exc:
            retry_or_raise(n_try, retry_count, exc, retry_interval)
            continue
        else:
            break

    print_timestamp(f'Succeeded to open url {repr(url_string)}')

    return response


def http_get(url, params=None, headers=None, text=False, retry_count=10, retry_interval=1):
    # set params
    if params is not None:
        url = url + '?' + urllib.parse.urlencode(params)
    if headers is None:
        if text:
            headers = {'Accept': 'text/plain'}
        else:
            headers = {'Accept': 'application/json'}
    # main
    req = urllib.request.Request(url, headers=headers, method='GET')
    return http_send_request(req, text, retry_count, retry_interval)


def http_post(url, data, params=None, headers=None, text=False, retry_count=10, retry_interval=1):
    # set params
    data = json.dumps(data).encode('ascii')
    if params is not None:
        url = url + '?' + urllib.parse.urlencode(params)
    if headers is None:
        if text:
            headers = {'Accept': 'text/plain'}
        else:
            headers = {'Accept': 'application/json'}
    # main
    req = urllib.request.Request(url, data=data, headers=headers, method='POST')
    return http_send_request(req, text, retry_count, retry_interval)


def http_send_request(req, text, retry_count=10, retry_interval=1, urlopen_timeout=5):
    with http_run_urlopen(
        req, 
        retry_count=retry_count, 
        retry_interval=retry_interval, 
        urlopen_timeout=urlopen_timeout,
    ) as response:
        if text:
            result = response.read().decode('utf-8')
        else:
            result = json.loads(response.read())

    return result


def download(url, path, retry_count=10, retry_interval=1, urlopen_timeout=5):
    while True:
        try:
            with http_run_urlopen(
                url, 
                retry_count=retry_count, 
                retry_interval=retry_interval, 
                urlopen_timeout=urlopen_timeout,
            ) as response:
                with open(path, 'wb') as outfile:
                    shutil.copyfileobj(response, outfile)
        except TimeoutError as exc:
            print_timestamp(f'Retrying due to TimeoutError (Exception: {exc})')
            continue
        else:
            break


def download_wget(url, path):
    subprocess.run(
        ['wget', '-O', path, url],
    )


def unzip(src, dest, rm_src=False):
    with gzip.open(src, 'rb') as infile:
        with open(dest, 'wb') as outfile:
            shutil.copyfileobj(infile, outfile)      
    if rm_src:
        os.remove(src)


@deco.get_deco_arg_choices({'mode': ('r', 'w', 'a')})
def openfile(fname, mode='r'):
    mode = mode + 't'

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


def fileiter(path, sep='\t', remove_leading_hashes=True, skip_double_hashes=True):
    infile = openfile(path, 'r')

    while True:
        line = next(infile)
        if skip_double_hashes:
            if line.startswith('##'):
                continue
            else:
                break
        else:
            break
        
    headerline = line
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

    if logging_lineno is None:
        for line in line_iterator:
            yield line
    else:
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

def pairwise(iterable):
    iterable = iter(iterable)
    try:
        x0 = next(iterable)
        x1 = next(iterable)
    except StopIteration:
        return

    while True:
        yield (x0, x1)
        try:
            x2 = next(iterable)
        except StopIteration:
            return
        x0 = x1
        x1 = x2
        continue


def _multi_selector_base(sequence, comparing_func, key=None, with_target_val=False):
    sequence = list(sequence)
    if key is None:
        values = sequence
    else:
        values = [key(x) for x in sequence]
    target_val = comparing_func(values)
    is_target_val = [x == target_val for x in values]

    if with_target_val:
        return list(itertools.compress(sequence, is_target_val)), target_val
    else:
        return list(itertools.compress(sequence, is_target_val))


def multi_min(sequence, key=None, with_target_val=False):
    return _multi_selector_base(sequence, min, key=key, with_target_val=with_target_val)


def multi_max(sequence, key=None, with_target_val=False):
    return _multi_selector_base(sequence, max, key=key, with_target_val=with_target_val)


def range_overlaps(range1, range2):
    return not ((range1.stop <= range2.start) or (range1.start >= range2.stop))


def range_included(range1, range2):
    """Returns:
        True if range1 is included in range2
    """
    return (range1.start >= range2.start) and (range1.stop <= range2.stop)


def range_intersection(range1, range2):
    start = max(range1.start, range2.start)
    stop = min(range1.stop, range2.stop)
    if stop > start:
        return range(start, stop)
    else:
        return None


def range_subtraction(range1, range2):
    isec = range_intersection(range1, range2)
    if isec is None:
        return range1
    else:
        if isec.start == range1.start:
            if isec.stop == range1.stop:
                return None
            else:
                return range(isec.stop, range1.stop)
        else:
            if isec.stop == range1.stop:
                return range(range1.start, isec.start)
            else:
                return (range(range1.start, isec.start), range(isec.stop, range1.stop))


# timeout (https://daeguowl.tistory.com/139)

class TimeoutError(Exception):
    pass


def timeout(seconds, error_message=''):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)
        
        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            #signal.alarm(seconds)
            signal.setitimer(signal.ITIMER_REAL, seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return functools.wraps(func)(wrapper)

    return decorator


def grouper(iterable, n):
    args = [iter(iterable)] * n
    for subiter in itertools.zip_longest(*args, fillvalue=None):
        yield tuple(x for x in subiter if x is not None)


def bernoulli_rv_generator(p, block_size=int(1e5)):
    while True:
        rvs = scipy.stats.bernoulli.rvs(p=p, loc=0, size=block_size)
        for x in rvs:
            yield x


def bernoulli_iterator(iterable, p, block_size=int(1e5)):
    """Args:
        p: probability of selection
    """
    for x, rv in zip(iter(iterable), bernoulli_rv_generator(p, block_size=block_size)):
        if rv:  # True if 1
            yield x


# from Itertools Recipes
def grouper_Itertools_Recipes(iterable, n, *, incomplete='fill', fillvalue=None):
    "Collect data into non-overlapping fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, fillvalue='x') --> ABC DEF Gxx
    # grouper('ABCDEFG', 3, incomplete='strict') --> ABC DEF ValueError
    # grouper('ABCDEFG', 3, incomplete='ignore') --> ABC DEF
    args = [iter(iterable)] * n
    if incomplete == 'fill':
        return itertools.zip_longest(*args, fillvalue=fillvalue)
    if incomplete == 'strict':
        return zip(*args, strict=True)
    if incomplete == 'ignore':
        return zip(*args)
    else:
        raise ValueError('Expected fill, strict, or ignore')


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
        

# memory usage check
def get_rss(mode='total', unit='b'):
    assert mode in ('self', 'children', 'total')
    assert unit in ('b', 'k', 'm', 'g')

    proc = psutil.Process()
    children = proc.children(recursive=True)
    rss_self = proc.memory_info().rss  # in bytes
    rss_children = sum(p.memory_info().rss for p in children)  # in bytes

    if mode == 'self':
        rss = rss_self
    elif mode == 'children':
        rss = rss_children
    elif mode == 'total':
        rss = rss_self + rss_children

    if unit == 'b':
        rss = rss
    elif unit == 'k':
        rss = rss / 1024
    elif unit == 'm':
        rss = rss / 1024**2
    elif unit == 'g':
        rss = rss / 1024**3

    return rss


# range operations
def check_overlaps(start0_1, end0_1, start0_2, end0_2):
    return (start0_1 < end0_2) and (end0_1 > start0_2)


# misc
def get_indexes_of_array(values, ordered_keys):
    result = np.zeros(len(values))
    for idx, key in enumerate(ordered_keys):
        result[np.where(values == key)[0]] = idx

    return result


def nanaverage(values, weights):
    #assert isinstance(values, np.ndarray)
    #assert isinstance(weights, np.ndarray)

    values = np.array(values)
    weights = np.array(weights)

    selector = ~np.isnan(values)
    new_values = values[selector]
    new_weights = weights[selector]
    return np.average(new_values, weights=new_weights)


def funclogger(msg):
    funcname = inspect.stack()[1].function
    print_timestamp(f'{funcname}: {msg}')
    

def array_grouper(arr, omit_values=False):
    assert arr.ndim in (1, 2)

    diff = np.empty(arr.shape[0], dtype=bool)
    diff[0] = True
    if arr.ndim == 1:
        diff[1:] = np.diff(arr)
    elif arr.ndim == 2:
        diff[1:] = np.diff(arr, axis=0).any(axis=1)

    indexes = np.nonzero(diff)[0]

    if omit_values:
        values = None
    else:
        values = arr[indexes]

    counts = np.empty(indexes.shape, dtype=int)
    counts[:-1] = np.diff(indexes)
    counts[-1] = arr.shape[0] - indexes[-1]

    groupkey = np.repeat(np.arange(len(counts)), counts)
        
    return values, counts, groupkey
                

def get_ranks(arr):
    ranks = scipy.stats.rankdata(data, method='max')
    return ranks / len(ranks)


def arg_into_list(arg):
    if not isinstance(arg, (tuple, list)):
        arg = [arg]
    return arg

arg_to_list = arg_into_list
        

def arg_into_df(arg):
    if arg is None:
        return None
    elif isinstance(arg, pd.DataFrame):
        check_having_coord_cols(arg)
        return arg
    elif isinstance(arg, pr.PyRanges):
        return arg.df
    else:
        raise Exception(f'Argument must be either None, pd.DataFrame, or pr.PyRanges')


def arg_into_gr(arg):
    if arg is None:
        return None
    elif isinstance(arg, pd.DataFrame):
        check_having_coord_cols(arg)
        return pr.PyRanges(arg)
    elif isinstance(arg, pr.PyRanges):
        return arg
    else:
        raise Exception(f'Argument must be either None, pd.DataFrame, or pr.PyRanges')


def check_having_coord_cols(arg):
    if not (
        {'Chromosome', 'Start', 'End'}.issubset(arg.columns)
        or arg.index.names == ['Chromosome', 'Start', 'End']
    ):
        raise Exception(f'Input dataframe must have columns "Chromosome", "Start", and "End".')


def get_mode(data, xs=np.arange(0, 0.51, 0.01)):
    kernel = scipy.stats.gaussian_kde(data)
    ys = kernel(xs)
    return xs[np.argmax(ys)]


def prefix_chr(s):
    if s.startswith('chr'):
        return s
    else:
        return 'chr' + s


def shorten_int(numlist, n_after_dot=3):
    mapping = {0: '', 1: 'K', 2: 'M', 3: 'G'}

    log = np.log10(numlist)
    assert (log < 12).all(), f'Numbers greater than or equal to 10^12 are not allowed.'

    #qs = [int(y) for y in (log / 3)]
    qs = np.floor(log / 3)

    suffixes = [mapping[x] for x in qs]
    #new_numlist = numlist / (10 ** (np.array(qs) * 3))
    new_numlist = numlist / (10 ** (qs * 3))
    formatter = f'.{n_after_dot}f'
    result = [
        (
            f'{int(x)}' 
            if y == '' else 
            f'{{x:{formatter}}} {{y}}'.format(x=x, y=y)
        )
        for x, y in zip(new_numlist, suffixes)
    ]
    return result
    

def mean_mad(values):
    values = np.array(values)
    return np.mean(np.abs(values - np.mean(values)))


def median_mad(values):
    values = np.array(values)
    return np.median(np.abs(values - np.median(values)))


def get_diffmean(values, weights=None):
    assert values.ndim == 1

    if weights is None:
        weights = np.ones_like(values)
    else:
        assert weights.shape == values.shape

    indexes = np.triu_indices(values.shape[0], k=1)
    diffs = values[indexes[0]] - values[indexes[1]]
    diff_weights = weights[indexes[0]] + weights[indexes[1]]

    return np.average(np.abs(diffs), weights=diff_weights)
    

#def broadcast_args(args):
#    return np.broadcast_arrays(*[np.atleast_1d(x) for x in args])


def get_groupby_keys(array):
    array = np.asarray(array)
    rolled = np.roll(array, 1)
    diffs = (arr != rolled)
    diffs[0] = True
    return np.cumsum(diffs)


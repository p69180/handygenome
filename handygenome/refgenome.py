import os
import itertools
import functools
import pprint
import re
import collections

import pysam
import pandas as pd
import pyranges as pr

import handygenome
import handygenome.network as network
import handygenome.tools as tools
import handygenome.logutils as logutils
import handygenome.deco as deco
import handygenome.publicdb.ncbi as libncbi
import handygenome.fastahandler as fastahandler


# constants

MT_CHROMS = ('chrM', 'MT')
CHROM_PATSTRING_CHR = '(?P<chr>chr)'
CHROM_PATSTRING_NUMBER = '(?P<number>0*(?P<number_proper>[1-9][0-9]*))'
CHROM_PATSTRING_XY = '(?P<xy>[XY])'
PAT_NUMERIC_CHROM = re.compile(f'{CHROM_PATSTRING_CHR}?{CHROM_PATSTRING_NUMBER}')
PAT_ASSEMBLED_CHROM = re.compile(f'{CHROM_PATSTRING_CHR}?({CHROM_PATSTRING_NUMBER}|{CHROM_PATSTRING_XY})')
CHR1_NAME_PAT = re.compile(r'(chr)?0*1')

FASTA_URLS = {
    'GRCh37_1000genomes': 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz',
    'GRCh37_ucsc': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz',
    'GRCh38_1000genomes': 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa',
}


class UnknownRefverError(Exception):
    pass


class NoChr1NameError(Exception):
    pass


# Initialize REFVERINFO from the config file

class RefverInfo:
    mandatory_keys = (
        'standard',
        'aliases',
        'parent_standard',
        'fasta_path',
        'species',
    )
    optional_keys = (
        #'assembly_report_path',
        'repeatmasker_path',
        'edited_repeat_path',

        'raw_geneset_path',
        'edited_geneset_path',

        'raw_regulatory_path',
        'edited_regulatory_path',

        'dbsnp_path',
    )
    def __init__(self, initial_data):
        self.data = list()
        for refveritem in initial_data:
            self.add_refveritem(refveritem, reset_caches=False)
        self.reset_caches()

    def reset_caches(self):
        self.set_maps()
        self.set_caches()

    def set_maps(self):
        # aliases
        self.standard_aliases_map = dict((x['standard'], x['aliases']) for x in self.data)
        # species-to-standards
        self.species_standards_map = dict()
        for x in self.data:
            self.species_standards_map.setdefault(x['species'], list())
            self.species_standards_map[x['species']].append(x['standard'])
        # default fasta
        self.standard_fastapath_map = dict(
            (x['standard'], x['fasta_path']) for x in self.data
        )
        # repeat
        self.standard_repeatmasker_map = dict(
            (x['standard'], x['repeatmasker_path']) for x in self.data
        )
        self.standard_edited_repeat_map = dict(
            (x['standard'], x['edited_repeat_path']) for x in self.data
        )
        # geneset
        self.standard_raw_geneset_map = dict(
            (x['standard'], x['raw_geneset_path']) for x in self.data
        )
        self.standard_edited_geneset_map = dict(
            (x['standard'], x['edited_geneset_path']) for x in self.data
        )
        # regulatory
        self.standard_raw_regulatory_map = dict(
            (x['standard'], x['raw_regulatory_path']) for x in self.data
        )
        self.standard_edited_regulatory_map = dict(
            (x['standard'], x['edited_regulatory_path']) for x in self.data
        )

    def set_caches(self):
        self.default_fastas = dict()
        self.default_chromdicts = dict()

        self.assemblyspecs = dict()
        self.refverspecs = dict()

        self.repeat_handles = dict()
        self.geneset_handles = dict()
        self.regulatory_handles = dict()

    def refveritem_sanitycheck(self, refveritem):
        #1
        assert isinstance(refveritem, dict), (
            f'Added reference info item must be a dict object.'
        )
        #2
        assert set(self.__class__.mandatory_keys).issubset(refveritem.keys()), (
            f'Added reference info must include these keys: '
            f'{self.__class__.mandatory_keys}'
        )
        #3
        assert refveritem['standard'] not in self.list_standards(), (
            f'"standard" of input refveritem already exists in the database.'
        )
        #4
        if refveritem['parent_standard'] is not None:
            if refveritem['species'] is None:
                raise Exception(f'When "parent_standard" is set, "species" must be also set.')
            if refveritem['parent_standard'] not in libncbi.REFSEQPATHS.get_available_refvers(refveritem['species']):
                raise Exception(f'"parent_standard" and "species" does not match according to RefSeq database.')
            if refveritem['fasta_path'] is not None:
                refverspec_fasta = RefverSpec.from_fasta_path(refveritem['fasta_path'])
                assemblyspec = AssemblySpec.from_refver(refveritem['parent_standard'], refveritem['species'])
                refverspec_asmblspec = RefverSpec.from_assemblyspec(assemblyspec)
                if not refverspec_fasta.check_length_compatibility(refverspec_asmblspec):
                    raise Exception(
                        f'For the newly added refver item, the fasta file and RefSeq assembly report do not '
                        f'coincide with the length of assembled chromosomes.\n'
                        f'refveritem: {refveritem}'
                    )
        #5
        if (refveritem['fasta_path'] is None) and (refveritem['parent_standard'] is None):
            raise Exception(f'When "parent_standard" is not set, "fasta_path" must be set.')


    # chr1 length related methods
#    @staticmethod
#    def find_chr1_name(chrom_names):
#        chr1_candidates = [x for x in chrom_names if CHR1_NAME_PAT.fullmatch(x)]
#        if len(chr1_candidates) == 1:
#            return chr1_candidates[0]
#        else:
#            raise NoChr1NameError(
#                f'Unique chromosome name which looks like chr1 could not be found.'
#            )
#
#    def find_chr1_length(self, fasta):
#        chr1_name = self.find_chr1_name(fasta.references)
#        chr1_len = fasta.lengths[fasta.references.index(chr1_name)]
#        return chr1_len

    # modifiers
    def change_fasta_path(self, refver, fasta_path):
        refveritem = self.access_refveritem(refver)
        refveritem['fasta_path'] = fasta_path
        self.reset_caches()

    def add_refveritem(self, refveritem, reset_caches=True):
        refveritem = refveritem.copy()
        for key in self.__class__.optional_keys:
            refveritem.setdefault(key, None)
        self.refveritem_sanitycheck(refveritem)
        self.data.append(refveritem)
        if reset_caches:
            self.reset_caches()

    # standardize
    def standardize(self, refver):
        for standard, alias_list in self.standard_aliases_map.items():
            if (refver == standard) or (refver in alias_list):
                return standard
        raise UnknownRefverError(f'Unknown reference version: {refver}')

    # refveritem accessor
    def access_refveritem(self, refver):
        standard = self.standardize(refver)
        for x in self.data:
            if x['standard'] == standard:
                return x
        raise Exception(f'Could not find a refveritem corresponding to the refver "{refver}"')

    # convenience methods
    def list_standards(self):
        return [x['standard'] for x in self.data]

    def list_known_refvers(self):
        result = list()
        for x in self.data:
            result.append(x['standard'])
            result.extend(x['aliases'])
        return result

    def list_known_species(self):
        return tuple(set(x['species'] for x in self.data))

    def get_parent_standard(self, refver):
        refveritem = self.access_refveritem(refver)
        return refveritem['parent_standard']

    # map fetchers
    def get_aliases(self, refver):
        standard = self.standardize(refver)
        alias_list = self.standard_aliases_map[standard]
        return [standard] + alias_list

    def find_species(self, refver):
        standard = self.standardize(refver)
        for key, val in self.species_standards_map.items():
            if standard in val:
                return key
        raise Exception(
            f'A standard refver "{standard}" could not be found from a search over "species_standards_map"'
        )

    # datafile path fetcher
    def get_datafile_path(self, refver, map_name, refseq_pathfunc_name):
        standard = self.standardize(refver)
        config_path = getattr(self, map_name)[standard]
        if config_path is None:
            refveritem = self.access_refveritem(standard)
            assert refveritem['parent_standard'] is not None
            species = refveritem['species']
            result = getattr(libncbi, refseq_pathfunc_name)(standard, species)
        else:
            result = config_path
        return result

    def get_fasta_path(self, refver):
        return self.get_datafile_path(
            refver, 'standard_fastapath_map', 'get_edited_fasta_path',
        )

    def get_repeatmasker_path(self, refver):
        return self.get_datafile_path(
            refver, 'standard_repeatmasker_map', 'get_repeatmasker_path',
        )

    def get_repeat_path(self, refver):
        standard = self.standardize(refver)
        config_path = self.standard_edited_repeat_map[standard]
        if config_path is None:
            return None
        else:
            return config_path

    def get_geneset_path(self, refver):
        standard = self.standardize(refver)
        config_path = self.standard_edited_geneset_map[standard]
        if config_path is None:
            return None
        else:
            return config_path

    def get_regulatory_path(self, refver):
        standard = self.standardize(refver)
        config_path = self.standard_edited_regulatory_map[standard]
        if config_path is None:
            return None
        else:
            return config_path

    # cache fetchers
    def get_default_fasta(self, refver):
        standard = self.standardize(refver)
        if standard not in self.default_fastas.keys():
            self.default_fastas[standard] = pysam.FastaFile(self.get_fasta_path(standard))
        return self.default_fastas[standard]

    def get_default_chromdict(self, refver):
        standard = self.standardize(refver)
        if standard not in self.default_chromdicts.keys():
            self.default_chromdicts[standard] = ChromDict.from_fasta(self.get_default_fasta(standard))
        return self.default_chromdicts[standard]

    def get_assemblyspec(self, refver):
        parent_standard = self.access_refveritem(refver)['parent_standard']
        if parent_standard is None:
            return None
            #raise Exception(f'AssemblySpec not available for this refver.')

        if parent_standard not in self.assemblyspecs.keys():
            species = self.find_species(parent_standard)
            self.assemblyspecs[parent_standard] = AssemblySpec.from_refver(parent_standard, species)

        return self.assemblyspecs[parent_standard]

    def get_refverspec(self, refver):
        standard = self.standardize(refver)

        if standard not in self.refverspecs.keys():
            assemblyspec = self.get_assemblyspec(standard)
            if assemblyspec is not None:
                self.refverspecs[standard] = RefverSpec.from_assemblyspec(assemblyspec)
            else:
                fasta_path = self.get_fasta_path(standard)
                if fasta_path is None:
                    self.refverspecs[standard] = None
                else:
                    self.refverspecs[standard] = RefverSpec.from_fasta_path(fasta_path)

        return self.refverspecs[standard]

    def get_geneset_handle(self, refver):
        standard = self.standardize(refver)
        if standard not in self.geneset_handles.keys():
            self.geneset_handles[standard] = pysam.TabixFile(
                self.get_geneset_path(standard),
                parser=pysam.asGTF(),
            )
        return self.geneset_handles[standard]

    def get_regulatory_handle(self, refver):
        standard = self.standardize(refver)
        if standard not in self.regulatory_handles.keys():
            self.regulatory_handles[standard] = pysam.TabixFile(
                self.get_regulatory_path(standard),
                parser=pysam.asGTF(),
            )
        return self.regulatory_handles[standard]

    def get_repeat_handle(self, refver):
        standard = self.standardize(refver)
        if standard not in self.repeat_handles.keys():
            self.repeat_handles[standard] = pysam.TabixFile(
                self.get_repeat_path(standard),
                parser=pysam.asBed(),
            )
        return self.repeat_handles[standard]



# RefverDict

class RefverDict(collections.UserDict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        standard_list = REFVERINFO.list_standards()
        if not set(self.keys()).issubset(standard_list):
            raise Exception(
                f'RefverDict construction keys must be restricted to: '
                f'{standard_list}'
            )

    def __getitem__(self, key):
        key = REFVERINFO.standardize(key)
        try:
            result = super().__getitem__(key)
        except KeyError:
            raise Exception(f'Input reference version is not available.')

        return result

    def __contains__(self, key):
        key = REFVERINFO.standardize(key)
        return super().__contains__(key)

    def get_valid_keys(self):
        result = list()
        for standard in self.keys():
            result.append(standard)
            result.extend(REFVERINFO.standard_aliases_map[standard])
        return result

    ############################
    # candidates be deprecated #
    @classmethod
    def standardize(cls, refver):
        return REFVERINFO.standardize(refver)

    @classmethod
    def get_aliases(cls, refver):
        return REFVERINFO.get_aliases(refver)
    ############################



#DEFAULT_FASTA_PATHS = RefverDict(handygenome.OPTION['default_fasta'])
#DEFAULT_FASTAS = RefverDict(
#    {
#        refver: pysam.FastaFile(path) 
#        for refver, path in DEFAULT_FASTA_PATHS.items()
#    }
#)



# chromdict

class ChromDict(collections.OrderedDict):
    ################
    # initializers #
    ################

    @classmethod
    def _init_from_pysamobj(cls, pysamobj):
        result = cls()
        for chrom, length in zip(pysamobj.references, pysamobj.lengths):
            result[chrom] = length
        return result

    @classmethod
    def from_fasta_path(cls, fasta_path):
        with pysam.FastaFile(fasta_path) as pysamobj:
            return cls._init_from_pysamobj(pysamobj)

    @classmethod
    def from_fasta(cls, fasta):
        return cls._init_from_pysamobj(fasta)

    @classmethod
    def from_bam_path(cls, bam_path):
        with pysam.AlignmentFile(bam_path) as pysamobj:
            return cls._init_from_pysamobj(pysamobj)

    @classmethod
    def from_bam(cls, bam):
        return cls._init_from_pysamobj(bam)

    @classmethod
    def from_bamheader(cls, bamheader):
        return cls._init_from_pysamobj(bamheader)

    @classmethod
    def from_vcfheader(cls, vcfheader):
        result = cls()
        for contig in vcfheader.contigs.values():
            result[contig.name] = contig.length
        return result

    @classmethod
    def from_refver(cls, refver):
        fasta_path = REFVERINFO.get_fasta_path(refver)
        with pysam.FastaFile(fasta_path) as fasta:
            return cls.from_fasta(fasta)

    @classmethod
    def from_custom(cls, custom):
        """Args:
            custom: {
                'contigs': ['contig1', 'contig2', ...], 
                'lengths': [length1, length2, ...],
            }
        """
        result = cls()
        for chrom, length in zip(custom['contigs'], custom['lengths']):
            result[chrom] = length
        return result

    ##############
    # properties #
    ##############

    @property
    def contigs(self):
        return list(self.keys())

    @property
    def lengths(self):
        return list(self.values())

#    @deco.get_deco_num_set_differently(
#        (
#            'fasta_path', 'fasta', 'bam_path', 'bam', 
#            'vcfheader', 'bamheader', 'custom', 'refver',
#        ), 
#        1,
#    )
#    def __init__(self, fasta_path=None, fasta=None, bam_path=None, bam=None, 
#                 vcfheader=None, bamheader=None, custom=None, refver=None):
#        """
#        Args:
#            fasta: pysam.FastaFile object
#            bam: pysam.AlignmentFile object
#            vcfheader: pysam.VariantHeader object
#            bamheader: pysam.AlignmentHeader object
#            custom: {'contigs': ['contig1', 'contig2', ...], 
#                     'lengths': [length1, length2, ...] }
#        """
#
#        # set self dict
#        if vcfheader is not None:
#            for contig in vcfheader.contigs.values():
#                self[contig.name] = contig.length
#        elif custom is not None:
#            for chrom, length in zip(custom['contigs'], custom['lengths']):
#                self[chrom] = length
#        else:
#            if fasta_path is not None:
#                wrapper = pysam.FastaFile(fasta_path)
#            elif fasta is not None:
#                wrapper = fasta
#            elif bam_path is not None:
#                wrapper = pysam.AlignmentFile(bam_path)
#            elif bam is not None:
#                wrapper = bam
#            elif bamheader is not None:
#                wrapper = bamheader
#            elif refver is not None:
#                wrapper = DEFAULT_FASTAS[refver]
#    
#            for chrom, length in zip(wrapper.references, wrapper.lengths):
#                self[chrom] = length
#    
#            if any((x is not None) for x in (fasta_path, bam_path)):
#                wrapper.close()
#
#        # set contigs, lengths
#        self.contigs = list(self.keys())
#        self.lengths = list(self.values())

    ##########
    # others #
    ##########

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
            if PAT_ASSEMBLED_CHROM.fullmatch(x) is not None
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
        return [x for x in self.contigs if PAT_ASSEMBLED_CHROM.fullmatch(x) is not None]

    def to_gr(self, assembled_only=True, as_gr=True):
        result = pd.DataFrame({
            'Chromosome': self.contigs,
            'Start': 0,
            'End': self.lengths,
        })
        if assembled_only:
            selector = result['Chromosome'].apply(
                lambda x: PAT_ASSEMBLED_CHROM.fullmatch(x) is not None
            )
            result = result.loc[selector, :]

        if as_gr:
            return pr.PyRanges(result)
        else:
            return result

#    def to_interval_list(self):
#        intvlist = libinterval.IntervalList()
#        for contig, length in zip(self.contigs, self.lengths):
#            interval = libinterval.Interval(contig, start0=0, end0=length)
#            intvlist.append(interval)
#
#        return intvlist

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
            start0s = np.asarray(start0s)
        if end0s is not None:
            end0s = np.asarray(end0s)

        # result
        sortkey = self.get_chrompos_sortkey(chroms, start0s, end0s)
        if start0s is None:
            return chroms[sortkey]
        else:
            if end0s is None:
                return chroms[sortkey], start0s[sortkey]
            else:
                return chroms[sortkey], start0s[sortkey], end0s[sortkey]


#DEFAULT_CHROMDICTS = RefverDict({
#    key: ChromDict(fasta=val) for key, val in DEFAULT_FASTAS.items()
#    if key != 'GRCh37_hs37d5'
#})


############################
# Assemblyspec and related #
############################

def normalize_chrom(chrom, strip_chr=False):
    """Examples:
        chr1 -> chr1
        chrX -> chrX
        chr001 -> chr1
        1 -> chr1
        X -> chrX
        MT -> chrM
        NC_000001.10 -> NC_000001.10
        chr19_gl000209_random -> chr19_gl000209_random
    """
    mat = PAT_ASSEMBLED_CHROM.fullmatch(chrom)
    if mat is None:
        if chrom in MT_CHROMS:
            if strip_chr:
                return 'MT'
            else:
                return 'chrM'
        else:  # e.g. chr19_gl000209_random, NC_000001.10
            return chrom
    else:  # e.g. chr1, 1, X, chrY, chr001
        if mat.group('xy') is None:
            if strip_chr:
                return mat.group('number_proper')
            else:
                return 'chr' + mat.group('number_proper')
        else:
            assert mat.group('number') is None
            if strip_chr:
                return mat.group('xy')
            else:
                return 'chr' + mat.group('xy')


class AssemblySpec:
    """Attributes:
        data: { 
            'length': [249250621, 243199373, ...],
            'role': ['assembled-molecule', ... 'unlocalized-scaffold', ... ],
            'default': ['1', '2', ... 'HSCHR1_RANDOM_CTG5', ... ],
            'ucsc': ['chr1', 'chr2', ... 'chr1_gl000191_random', ... ],
            'nochr_plus_genbank': ['1', '2', ... 'GL000191.1', ... ],
            'genbank': ['CM000663.1', 'CM000664.1', ... ],
            'refseq': ['NC_000001.11', 'NC_000002.12', ... ],
        }
        dicts: {('ucsc', 'genbank'): dict, ... }
        chromdicts: {'ucsc': ChromDict object, ... }
        versions: All self.data keys except 'length' and 'role'
    """
    assembly_report_keydict = {
        'default': 'Sequence-Name',
        'ucsc': 'UCSC-style-name',
        'genbank': 'GenBank-Accn',
        'refseq': 'RefSeq-Accn',
        'length': 'Sequence-Length',
        'role': 'Sequence-Role',
    }
    def __init__(self, data):
        """
        Args:
            data: {
                'ucsc': ['chr1', 'chr2', ... ], 
                'refseq': ['NC_000001.11', 'NC_000002.12', ... ],
                ...
                }
        """

        self.data = data
        self.dicts = dict()
        self.chromdicts = dict()

        for key1, key2 in itertools.permutations(data.keys(), 2):
            dic = dict()
            chromdict_data = {'contigs': list(), 'lengths': list()}
            for val1, val2 in zip(self.data[key1], self.data[key2]):
                if val1 is not None:
                    dic[val1] = val2
                    if key2 == 'length':
                        chromdict_data['contigs'].append(val1)
                        chromdict_data['lengths'].append(val2)

            self.dicts[(key1, key2)] = dic
            if key2 == 'length':
                self.chromdicts[key1] = ChromDict(custom=chromdict_data)

        self.versions = list(
            set(self.data.keys()).difference(('length', 'role')))

    @classmethod
    def from_refver(cls, refseq_refver, species):
        #species = REFVERINFO.find_species(refver)
        assemblyfile_path = libncbi.get_assemblyfile_path(
            refseq_refver, species, force_download=False,
        )
        with open(assemblyfile_path, 'r') as infile:
            result = cls.from_assembly_report_file(infile)

        return result

    @classmethod
    def from_assembly_report_file(cls, stream):
        parsed_data = cls.parse_assembly_report(stream)
        result = cls(parsed_data)
        return result

    @classmethod
    def parse_assembly_report(cls, stream):
        data = {'length': list(), 'default': list(), 'ucsc': list(),
                'nochr_plus_genbank': list(), 'role': list(),
                'genbank': list(), 'refseq': list()}

        record_start = False
        for line in stream:
            if not record_start:
                if line.startswith('# Sequence-Name'):
                    keys = tools.get_linesp(re.sub('^#\s*', '', line))
                    record_start = True
                    continue
            else:
                linesp = tools.get_linesp(line)
                linedict = dict(
                    zip(
                        keys, 
                        [(None if x == 'na' else x) for x in linesp]
                    )
                )

                default = linedict[cls.assembly_report_keydict['default']]
                ucsc = linedict[cls.assembly_report_keydict['ucsc']]
                genbank = linedict[cls.assembly_report_keydict['genbank']]
                refseq = linedict[cls.assembly_report_keydict['refseq']]
                length = int(linedict[cls.assembly_report_keydict['length']])
                role = linedict[cls.assembly_report_keydict['role']]

                #assert genbank is not None

                data['default'].append(default)
                data['ucsc'].append(ucsc)
                data['genbank'].append(genbank)
                data['refseq'].append(refseq)
                data['length'].append(length)
                data['role'].append(role)

                if ucsc is None:
                    data['nochr_plus_genbank'].append(genbank)
                elif ucsc == 'chrM':
                    data['nochr_plus_genbank'].append('MT')
                else:
                    mat = PAT_ASSEMBLED_CHROM.fullmatch(ucsc)
                    if mat is None:
                        data['nochr_plus_genbank'].append(genbank)
                    else:
                        data['nochr_plus_genbank'].append(mat.group(2))

        return data

    def __repr__(self):
        tmp_keys = list()
        tmp_vals = list()
        for key, val in self.data.items():
            tmp_keys.append(key)
            tmp_vals.append(val)

        show_result = list()
        for tup in zip(*tmp_vals):
            show_result.append(str(dict(zip(tmp_keys, tup))))

        return '\n'.join(show_result)

    def convert(self, in_chrom, out_version):
        assert out_version in self.versions, (
            f'"out_version" argument must be one of {self.versions}')

        absent = True
        for version, names in self.data.items():
            if version in ('length', 'role'):
                continue
            if in_chrom in names:
                absent = False
                if version == out_version:
                    result = in_chrom
                else:
                    result = self.dicts[(version, out_version)][in_chrom]
                break

        if absent:
            raise Exception(f'Input chrom name "{in_chrom}" is not included '
                            f'in the database.')

        return result

    def get_alias_groups(self):
        data_keys = (
            'ucsc', 
            'genbank', 
            'refseq',
            'default', 
            'nochr_plus_genbank', 
        )
        alias_groups = list()
        for idx in range(len(self.data['length'])):
            group = set()
            for key in data_keys:
                chrom = self.data[key][idx]
                if chrom is None:
                    continue

                norm_chrom = normalize_chrom(chrom, strip_chr=False)
                norm_chrom_nochr = normalize_chrom(chrom, strip_chr=True)
                group.add(chrom)
                group.add(norm_chrom)
                group.add(norm_chrom_nochr)

            alias_groups.append(group)
        return alias_groups

    def get_chrom_converter(self, contigs, lengths):
        # check ChromLenghtHash
        hash_assemblyspec = RefverSpec.from_assemblyspec(self)
        hash_inputdata = RefverSpec.from_data(contigs, lengths)
        if not hash_assemblyspec.check_length_compatibility(hash_inputdata):
            raise Exception(f'Input contig names and lenghts does not match assemblyspec.')

        # main
        alias_groups = self.get_alias_groups()
        converter = dict()
        for chrom in contigs:
            for group in alias_groups:
                if chrom in group:
                    for alias in group:
                        assert alias not in converter
                        converter[alias] = chrom
                    continue
        return converter


#ASSEMBLYFILE_URLS = RefverDict(libncbi.collect_assemblyfile_urls())

#ASSEMBLYFILE_PATHS = {
#    key: os.path.join(ASSEMBLYFILE_DIR, os.path.basename(val))
#    for (key, val) in ASSEMBLYFILE_URLS.items()
#}


#def get_assemblyspec_data(refver):
#    assemblyfile_url = ASSEMBLYFILE_URLS[refver]
#    assembly_report_path = ASSEMBLYFILE_PATHS[refver]
#    if not os.path.exists(assembly_report_path):
#        network.download(assemblyfile_url, assembly_report_path)
#
#    return parse_assembly_report(assembly_report_path)
#
#
#SPECS = RefverDict(
#    {
#        refver: AssemblySpec(get_assemblyspec_data(refver))
#        for refver in ASSEMBLYFILE_URLS.keys()
#    }
#)


class RefverSpec:
    def __init__(self, chrom_list, length_list):
        self.lendict = dict((x, y) for (x, y) in zip(chrom_list, length_list))
        self.normchrom_lendict = self.make_normchrom_lendict(chrom_list, length_list)

    @staticmethod
    def make_normchrom_lendict(chrom_list, length_list):
        normchrom_lendict = dict()
        modified_chrom_groups = dict()
        for chrom, length in zip(chrom_list, length_list):
            modified_chrom = normalize_chrom(chrom, strip_chr=False)
            normchrom_lendict[modified_chrom] = length
            modified_chrom_groups.setdefault(modified_chrom, set())
            modified_chrom_groups[modified_chrom].add(
                (chrom, length)
            )

        # sanity check
        for modified_chrom, tuples in modified_chrom_groups.items():
            lengths = set(x[1] for x in tuples)
            if len(lengths) != 1:
                raise Exception(
                    f'Modified chrom names overlap but lengths are discordant:'
                    f'modified chrom: {modified_chrom}, original chrom - length pairs: {tuples}'
                )

        return normchrom_lendict

    @classmethod
    def from_data(cls, chrom_list, length_list):
        #lendict = cls.make_lendict(chrom_list, length_list)
        return cls(chrom_list, length_list)

    @classmethod
    def from_fasta(cls, fasta):
        return cls.from_data(fasta.references, fasta.lengths)

    @classmethod
    def from_fasta_path(cls, fasta_path):
        with pysam.FastaFile(fasta_path) as fasta:
            return cls.from_fasta(fasta)

    @classmethod
    def from_assemblyspec(cls, assemblyspec):
        assert isinstance(assemblyspec, AssemblySpec)

        chrom_list = list(
            itertools.chain(
                assemblyspec.data['default'],
                assemblyspec.data['ucsc'],
                assemblyspec.data['nochr_plus_genbank'],
                assemblyspec.data['genbank'],
                assemblyspec.data['refseq'],
            )
        )
        length_list = assemblyspec.data['length'] * 5
        # remove None chroms
        notnone_indexes = tuple((x is not None) for x in chrom_list)
        chrom_list = list(itertools.compress(chrom_list, notnone_indexes))
        length_list = list(itertools.compress(length_list, notnone_indexes))

        return cls.from_data(chrom_list, length_list)

    @staticmethod
    def compare_lendicts(self_lendict, other_lendict):
        self_chroms = set(self_lendict.keys())
        other_chroms = set(other_lendict.keys())

        self_only_chroms = self_chroms.difference(other_chroms)
        other_only_chroms = other_chroms.difference(self_chroms)
        common_chroms = self_chroms.intersection(other_chroms)
        length_matches = all(
            (self_lendict[x] == other_lendict[x])
            for x in common_chroms
        )
        return length_matches, common_chroms, self_only_chroms, other_only_chroms

    def check_length_compatibility(self, other):
        """Used for parent standard identification"""
        assert isinstance(other, RefverSpec)
        length_matches, common_chroms, self_only_chroms, other_only_chroms = (
            self.compare_lendicts(self.normchrom_lendict, other.normchrom_lendict)
        )
        return length_matches

    def check_issubset(self, other):
        assert isinstance(other, RefverSpec)
        length_matches, common_chroms, self_only_chroms, other_only_chroms = (
            self.compare_lendicts(self.lendict, other.lendict)
        )
        is_subset = ((len(self_only_chroms) == 0) and length_matches)
        return is_subset, other_only_chroms

    def get_hash(self, chroms=None):
        if chroms is None:
            chrom_len_tuples = iter(self.lendict.items())
        else:
            chrom_len_tuples = ((x, self.lendict[x]) for x in chroms)
        return tuple(sorted(chrom_len_tuples, key=(lambda x: x[0])))

    def compare(self, other):
        assert isinstance(other, RefverSpec)
        common_chroms = set(self.lendict.keys()).intersection(other.lendict.keys())
        return self.get_hash(common_chroms) == other.get_hash(common_chroms)
                
                    


####################
# Setup REFVERINFO #
####################

REFVERINFO = RefverInfo(handygenome.OPTION['refverinfo'])


###########################################################################
# convenience functions for utilizing and modifying the REFVERINFO object #
###########################################################################

def standardize_refver(refver):
    REFVERINFO.standardize(refver)

def add_refveritem(refveritem):
    REFVERINFO.add_refveritem(refveritem, reset_caches=True)
    

def change_fasta_path(refver, fasta_path):
    REFVERINFO.change_fasta_path(refver, fasta_path)


def list_known_refvers():
    return REFVERINFO.list_known_refvers()


def get_fasta_path(refver):
    return REFVERINFO.get_fasta_path(refver)


def get_default_fasta(refver):
    return REFVERINFO.get_default_fasta(refver)


def get_default_chromdict(refver):
    return REFVERINFO.get_default_chromdict(refver)


def get_assemblyspec(refver):
    return REFVERINFO.get_assemblyspec(refver)


def get_geneset_handle(refver):
    return REFVERINFO.get_geneset_handle(refver)
    

def get_regulatory_handle(refver):
    return REFVERINFO.get_regulatory_handle(refver)


def get_repeat_handle(refver):
    return REFVERINFO.get_repeat_handle(refver)


###############################
# reference version inference #
###############################

@deco.get_deco_num_set_differently(('chromdict', 'vcfheader', 'bamheader'), 1)
def infer_refver(chromdict=None, vcfheader=None, bamheader=None):
    if chromdict is not None:
        return infer_refver_chromdict(chromdict)
    elif vcfheader is not None:
        return infer_refver_vcfheader(vcfheader)
    elif bamheader is not None:
        return infer_refver_bamheader(bamheader)


def infer_parent_standard(contigs, lengths):
    query_refverspec = RefverSpec.from_data(contigs, lengths)


def infer_refver_base(contigs, lengths):
    query_hash = RefverSpec.from_data(contigs, lengths)
    compare_results = list()
    for standard in REFVERINFO.list_standards():
        target_hash = REFVERINFO.get_refverspec(standard)
        if target_hash is None:
            continue
        is_subset, other_only_chroms = target_hash.check_issubset(query_hash)
        compare_results.append((standard, is_subset, other_only_chroms))

    subset_refvers = [x for x in compare_results if x[1]]
    if len(subset_refvers) == 1:
        return subset_refvers[0][0]
    elif len(subset_refvers) > 1:
        minimal_descrepancy_refvers = tools.multi_min(
            subset_refvers, key=(lambda x: len(x[2]))
        )
        if len(minimal_descrepancy_refvers) == 1:
            return minimal_descrepancy_refvers[0][0]
        else:
            raise Exception(f'Cannot break tie: {subset_refvers}')
    elif len(subset_refvers) == 0:
        return None
            
        

    

    elif len(candidates) > 1:
        raise Exception(f'Matches assembled chromosome lengths with more than one reference genomes({candidates})')
    else:
        return candidates[0]

    if len(candidates) == 0:
        raise Exception(f'Cannot infer refver: unknown chromosome lengths')
    elif len(candidates) > 1:
        raise Exception(f'Matches assembled chromosome lengths with more than one reference genomes({candidates})')
    else:
        return candidates[0]


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


##################################################
# handling of pseudoautosomal region & runs of N #
##################################################

N_REGIONS_DIR = os.path.join(handygenome.DIRS['data'], 'N_region')
if not os.path.exists(N_REGIONS_DIR):
    os.mkdir(N_REGIONS_DIR)

PARS = RefverDict(
    # 0-based half-open system
    # ref: https://www.ncbi.nlm.nih.gov/grc/human
    {  
        'GRCh37': {
            'PAR1_X': ('X', 60_000, 2_699_520),
            'PAR2_X': ('X', 154_931_043, 155_260_560),
            'PAR1_Y': ('Y', 10_000, 2_649_520),
            'PAR2_Y': ('Y', 59_034_049, 59_373_566),
        }, 
        'GRCh38': {
            'PAR1_X': ('X', 10_000, 2_781_479),
            'PAR2_X': ('X', 155_701_382, 156_030_895),
            'PAR1_Y': ('Y', 10_000, 2_781_479),
            'PAR2_Y': ('Y', 56_887_902, 57_217_415),
        },
    }
)

PAR_GRS = dict()
for key, val in PARS.items():
    df = pd.DataFrame.from_records(
        iter(val.values()), 
        columns=('Chromosome', 'Start', 'End'),
    )
    gr = pr.PyRanges(df).sort()
    PAR_GRS[key] = gr
PAR_GRS = RefverDict(PAR_GRS)


def get_par_gr(refver):
    chromdict = ChromDict.from_refver(refver)
    par_gr = PAR_GRS[refver].copy()  # X, Y
    if chromdict.is_chr_prefixed:
        par_gr.Chromosome = 'chr' + par_gr.Chromosome
    return par_gr


def make_N_region_gr(fasta):
    pat = re.compile('N+')
    chroms = list()
    start0s = list()
    end0s = list()
    for chrom in fasta.references:
        seq = fasta.fetch(chrom)
        offset = 0
        while True:
            mat = pat.search(seq)
            if mat is None:
                break
            else:
                span = mat.span()
                chroms.append(chrom)
                start0s.append(offset + span[0])
                end0s.append(offset + span[1])

                offset += span[1]
                seq = seq[span[1]:]

    return pr.PyRanges(chromosomes=chroms, starts=start0s, ends=end0s)


def get_N_regionfile_path(refver):
    refver = RefverDict.standardize(refver)
    return os.path.join(N_REGIONS_DIR, f'{refver}.N_regions.tsv.gz')


def write_N_regionfile(refver):
    fasta = get_default_fasta(refver)
    gr = make_N_region_gr(fasta)
    N_regionfile_path = get_N_regionfile_path(refver)
    gr.df.to_csv(N_regionfile_path, sep='\t', index=False)


def get_N_region_gr(refver):
    N_regionfile_path = get_N_regionfile_path(refver)
    if not os.path.exists(N_regionfile_path):
        print('Creating a N region file. It may take a few minutes.')
        write_N_regionfile(refver)
    return pr.PyRanges(
        pd.read_csv(N_regionfile_path, sep='\t', header=0, dtype={'Chromosome': str})
    )
        


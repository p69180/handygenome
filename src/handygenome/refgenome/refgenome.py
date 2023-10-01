import os
import itertools
import functools
import pprint
import re
import collections
import inspect
import shutil
import hashlib
import json

import pysam
import pandas as pd
import pyranges as pr

import handygenome
import handygenome.network as network
import handygenome.tools as tools
import handygenome.logutils as logutils
import handygenome.deco as deco
import handygenome.publicdb.ncbi as libncbi
import handygenome.publicdb.ncbi_cache as libncbicache
import handygenome.refgenome.refverfile_wogdf as refverfile_wogdf


# constants

MT_CHROMS = ('chrM', 'MT')

CHROM_PATSTRING_CHR = '(?P<chr>chr)'
CHROM_PATSTRING_NUMBER = '(?P<number>0*(?P<number_proper>[1-9][0-9]*))'
CHROM_PATSTRING_XY = '(?P<xy>[XY])'

PAT_NUMERIC_CHROM = re.compile(f'{CHROM_PATSTRING_CHR}?{CHROM_PATSTRING_NUMBER}', flags=re.I)
PAT_ASSEMBLED_CHROM = re.compile(f'{CHROM_PATSTRING_CHR}?({CHROM_PATSTRING_NUMBER}|{CHROM_PATSTRING_XY})', flags=re.I)
PAT_XY_CHROM = re.compile(f'{CHROM_PATSTRING_CHR}?{CHROM_PATSTRING_XY}', flags=re.I)

CHR1_NAME_PAT = re.compile(r'(chr)?0*1')

FASTA_URLS = {
    'GRCh37_1000genomes': 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz',
    'GRCh37_ucsc': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz',
    'GRCh38_1000genomes': 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa',
}

REFVER_DIR = os.path.join(handygenome.USERDATA_DIR, 'reference_files')
os.makedirs(REFVER_DIR, exist_ok=True)


############################
# chromosome name handlers #
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


def check_assembled_chrom(chrom):
    return PAT_ASSEMBLED_CHROM.fullmatch(chrom) is not None


def choose_assembled_chroms(chroms):
    return list(filter(check_assembled_chrom, chroms))


def check_xy_chrom(chrom):
    return PAT_XY_CHROM.fullmatch(chrom) is not None


##############
# exceptions #
##############

class UnknownRefverError(Exception):
    pass


class NoChr1NameError(Exception):
    pass


class AssemblyspecUnavailableError(Exception):
    pass


class NoXYError(Exception):
    pass


###############################
# RefverEntity and RefverInfo #
###############################

class RefverEntity(collections.UserDict):
    mandatory_keys = (
        'standard',
        'aliases',
        'RefSeq_refver',
        'fasta_path',
        'species',
    )
    optional_keys = (
        'repeatmasker_path',
        'edited_repeat_path',

        'raw_geneset_path',
        'edited_geneset_path',

        'raw_regulatory_path',
        'edited_regulatory_path',

        'dbsnp_path',
    )

    @classmethod
    def from_config(cls, config_data, calc_checksum=True):
        result = cls()
        result.calc_checksum = calc_checksum
        result.update(config_data)
        for key in cls.optional_keys:
            result.setdefault(key, None)
        result.sanitycheck()

        return result

    def sanitycheck(self):
        assert set(self.__class__.mandatory_keys).issubset(self.keys()), (
            f'Data must include these keys: {self.__class__.mandatory_keys}'
        )

        if self['RefSeq_refver'] is not None:
            if self['species'] is None:
                raise Exception(f'When "RefSeq_refver" is set, "species" must be also set.')
            if self['RefSeq_refver'] not in libncbi.GENOME_PATHS.get_available_refvers(self['species']):
                raise Exception(f'"RefSeq_refver" and "species" does not match according to RefSeq database.')
            if self['fasta_path'] is not None:
                chromlenhash_fasta = ChromLengthHash.from_fasta_path(self['fasta_path'])
                assemblyspec = AssemblySpec.from_RefSeq_refver(self['RefSeq_refver'], self['species'])
                chromlenhash_asmblspec = ChromLengthHash.from_assemblyspec(assemblyspec)
                if not chromlenhash_fasta.check_length_compatibility(chromlenhash_asmblspec):
                    raise Exception(
                        f'The fasta file and RefSeq assembly report do not '
                        f'coincide with the length of chromosomes.'
                    )
        ###
        if self['RefSeq_refver'] is None:
            if self['fasta_path'] is None:
                raise Exception(f'When "RefSeq_refver" is not set, "fasta_path" must be set.')

    def set_cached_properties(self):
        logutils.log(
            f'Initializing RefverEntity caches (refver={self["standard"]})', 
            level='info', add_locstring=False,
        )
        for key in (
            'fasta', 
            'chromdict', 
            'assemblyspec', 
            'chrom_converter', 
            'chromlenhash', 
            #'geneset_handle',
        ):
            try:
                _ = getattr(self, key)
            except AssemblyspecUnavailableError:
                pass

    ############################
    # data paths as attributes #
    ############################

    @property
    def datadir(self):
        datadir = os.path.join(REFVER_DIR, self['standard'])
        os.makedirs(datadir, exist_ok=True)
        return datadir

    @property
    def processed_fasta_path(self):
        return os.path.join(self.datadir, 'genome.fasta')

    @property
    def unedited_fasta_checksum_path(self):
        return os.path.join(self.datadir, 'raw_fasta_checksum')

    @property
    def unedited_fasta_wascustom_path(self):
        return os.path.join(self.datadir, 'raw_fasta_wascustom.json')

    @property
    def processed_geneset_path(self):
        return os.path.join(self.datadir, 'geneset.gff3.gz')

    ########################
    # important attributes #
    ########################

    @functools.cached_property
    def fasta(self):
        refverfile_wogdf.prepare_processed_fasta(
            standard_refver=self['standard'],
            edited_path=self.processed_fasta_path,
            custom_path=self['fasta_path'],
            RefSeq_refver=self['RefSeq_refver'],
            species=self['species'],
            chrom_converter=(None if self['RefSeq_refver'] is None else self.chrom_converter),
            unedited_file_checksum_path=self.unedited_fasta_checksum_path,
            unedited_file_wascustom_path=self.unedited_fasta_wascustom_path,
            calc_checksum=self.calc_checksum,
        )
        return pysam.FastaFile(self.processed_fasta_path)
    
    @functools.cached_property
    def chromdict(self):
        return ChromDict.from_fasta(self.fasta, refver=self['standard'])

    @functools.cached_property
    def assemblyspec(self):
        if self['RefSeq_refver'] is None:
            raise AssemblyspecUnavailableError(f'AssemblySpec not available for this refver.')
        else:
            return AssemblySpec.from_RefSeq_refver(self['RefSeq_refver'], self['species'])

    @functools.cached_property
    def chrom_converter(self):
        if self['fasta_path'] is None:
            return self.assemblyspec.get_chrom_converter_from_self()
        else:
            with pysam.FastaFile(self['fasta_path']) as fasta:
                contigs = fasta.references
                lengths = fasta.lengths
            return self.assemblyspec.get_chrom_converter(contigs, lengths)

    @functools.cached_property
    def chromlenhash(self):
        return ChromLengthHash.from_fasta(self.fasta)

    @functools.cached_property
    def geneset_handle(self):
        self.prepare_processed_geneset()
        return pysam.TabixFile(self.processed_geneset_path, parser=pysam.asGTF())


class RefverInfo:
    # init and helpers
    def __init__(self, initial_data, calc_checksum=True):
        self.refver_entities = list()
        for config_data in initial_data:
            refver_entity = RefverEntity.from_config(config_data, calc_checksum=calc_checksum)
            self.add_refver_entity(
                refver_entity, 
                reset_cache=False,
                set_refver_entity_cache=False,
            )

        #for refver_entity in self.refver_entities:
        #    refver_entity.set_cached_properties()

        self.reset_cache()

    def add_refver_entity(self, refver_entity, reset_cache=True, set_refver_entity_cache=True):
        # sanity check
        assert refver_entity['standard'] not in self.list_standards(), (
            f'"standard" of input RefverEntity already exists in the database.'
        )
        assert not set(refver_entity['aliases']).intersection(self.list_all_aliases()), (
            f'"aliases" of input RefverEntity has overlap with already existing ones.'
        )
        # main
        self.refver_entities.append(refver_entity)
        #if set_refver_entity_cache:
        #    refver_entity.set_cached_properties()
        if reset_cache:
            self.reset_cache()

    def reset_cache(self):
        # accessor to RefverEntity
        self.refver_entity_accessor = dict()
        for x in self.refver_entities:
            self.refver_entity_accessor[x['standard']] = x
            for alias in x['aliases']:
                self.refver_entity_accessor[alias] = x

        # aliases
        self.standard_aliases_map = dict((x['standard'], x['aliases']) for x in self.refver_entities)

        # species-to-standards
        self.species_standards_map = dict()
        for x in self.refver_entities:
            self.species_standards_map.setdefault(x['species'], list())
            self.species_standards_map[x['species']].append(x['standard'])

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

    def get_refver_entity(self, refver):
        try:
            return self.refver_entity_accessor[refver]
        except KeyError:
            raise UnknownRefverError(f'Unknown reference version: {refver}')

    # modifiers
    def change_fasta_path(self, refver, fasta_path):
        refver_entity = self.get_refver_entity(refver)
        refver_entity['fasta_path'] = fasta_path
        self.reset_cache()

    # convenience methods
    def standardize(self, refver):
        return self.get_refver_entity(refver)['standard']

    def list_standards(self):
        return list(x['standard'] for x in self.refver_entities)

    def list_all_aliases(self):
        return list(
            itertools.chain.from_iterable(
                x['aliases'] for x in self.refver_entities
            )
        )

    def list_known_refvers(self):
        return self.list_standards() + self.list_all_aliases()

    #def list_known_species(self):
    #    return tuple(set(x['species'] for x in self.refver_entities))

    def get_RefSeq_refver(self, refver):
        refver_entity = self.get_refver_entity(refver)
        RefSeq_refver = refver_entity['RefSeq_refver']
        if RefSeq_refver is None:
            raise Exception(f'This reference version does not have a matched standard RefSeq version name.')
        return RefSeq_refver

    def get_RefSeq_standard(self, refver):
        refver_entity = self.get_refver_entity(refver)
        RefSeq_refver = refver_entity['RefSeq_refver']
        species = refver_entity['species']
        if RefSeq_refver is None:
            raise Exception(f'This reference version does not have a matched standard RefSeq version name.')
        return RefSeq_refver, species
        #return {
        #    'RefSeq_refver': refver_entity['RefSeq_refver'], 
        #    'species': refver_entity['species'],
        #}

    def get_aliases(self, refver):
        refver_entity = self.get_refver_entity(refver)
        return [refver_entity['standard']] + refver_entity['aliases']

    def find_species(self, refver):
        standard = self.standardize(refver)
        for key, val in self.species_standards_map.items():
            if standard in val:
                return key
        raise Exception(
            f'A standard refver "{standard}" could not be found from a search over "species_standards_map"'
        )

    # datafile path fetcher
#    def get_datafile_path(self, refver, map_name, refseq_pathfunc_name):
#        standard = self.standardize(refver)
#        config_path = getattr(self, map_name)[standard]
#        if config_path is None:
#            refver_entity = self.get_refver_entity(standard)
#            assert refver_entity['RefSeq_refver'] is not None
#            species = refver_entity['species']
#            result = getattr(libncbi, refseq_pathfunc_name)(standard, species)
#        else:
#            result = config_path
#        return result

    ###########################
    # reference data fetchers #
    ###########################

    def get_fasta_path(self, refver):
        refver_entity = self.get_refver_entity(refver)
        return refver_entity.processed_fasta_path

    def get_unedited_fasta_path(self, refver):
        refver_entity = self.get_refver_entity(refver)
        return refver_entity.get_unedited_fasta_path()

    def get_fasta(self, refver):
        refver_entity = self.get_refver_entity(refver)
        return refver_entity.fasta

    def get_chromdict(self, refver):
        refver_entity = self.get_refver_entity(refver)
        return refver_entity.chromdict

    def get_assemblyspec(self, refver):
        refver_entity = self.get_refver_entity(refver)
        return refver_entity.assemblyspec

    def get_chrom_converter(self, refver):
        refver_entity = self.get_refver_entity(refver)
        return refver_entity.chrom_converter

    def get_repeatmasker_path(self, refver):
    #    return self.get_datafile_path(
    #        refver, 'standard_repeatmasker_map', 'get_repeatmasker_path',
    #    )
        pass

    def get_repeat_path(self, refver):
        refver_entity = self.get_refver_entity(refver)
        return refver_entity['edited_repeat_path']

    def get_geneset_path(self, refver):
        refver_entity = self.get_refver_entity(refver)
        return refver_entity['edited_geneset_path']

    def get_regulatory_path(self, refver):
        refver_entity = self.get_refver_entity(refver)
        return refver_entity['edited_regulatory_path']

    def get_chromlenhash(self, refver):
        refver_entity = self.get_refver_entity(refver)
        return refver_entity.chromlenhash

#        standard = self.standardize(refver)
#
#        if standard not in self.cache['chromlenhash'].keys():
#            assemblyspec = self.get_assemblyspec(standard)
#            if assemblyspec is not None:
#                self.cache['chromlenhash'][standard] = ChromLengthHash.from_assemblyspec(assemblyspec)
#            else:
#                fasta_path = self.get_fasta_path(standard)
#                if fasta_path is None:
#                    self.cache['chromlenhash'][standard] = None
#                else:
#                    self.cache['chromlenhash'][standard] = ChromLengthHash.from_fasta_path(fasta_path)
#
#        return self.cache['chromlenhash'][standard]

    def get_geneset_handle(self, refver):
        standard = self.standardize(refver)
        if standard not in self.cache['geneset_handle'].keys():
            self.cache['geneset_handle'][standard] = pysam.TabixFile(
                self.get_geneset_path(standard),
                parser=pysam.asGTF(),
            )
        return self.cache['geneset_handle'][standard]

    def get_regulatory_handle(self, refver):
        standard = self.standardize(refver)
        if standard not in self.cache['regulatory_handle'].keys():
            self.cache['regulatory_handle'][standard] = pysam.TabixFile(
                self.get_regulatory_path(standard),
                parser=pysam.asGTF(),
            )
        return self.cache['regulatory_handle'][standard]

    def get_repeat_handle(self, refver):
        standard = self.standardize(refver)
        if standard not in self.cache['repeat_handle'].keys():
            self.cache['repeat_handle'][standard] = pysam.TabixFile(
                self.get_repeat_path(standard),
                parser=pysam.asBed(),
            )
        return self.cache['repeat_handle'][standard]

    # chromosome name converter
#    def get_chrom_converter(self, refver):
#        """For refver entity without a custom fasta path, 
#        ucsc-style chromosome names of assemblyspec are used.
#        """
#        # fetch refver_entity
#        refver_entity = self.access_refver_entity(refver)
#        if refver_entity['RefSeq_refver'] is None:
#            raise Exception(f'Reference version must be based on a RefSeq refver to get a chrom converter.')
#
#        # get assemblyspec
#        assemblyspec = self.get_assemblyspec(refver)
#
#        # get contigs, lengths
#        if refver_entity['fasta_path'] is None:
#            return assemblyspec.get_chrom_converter_from_self()
#        else:
#            with pysam.FastaFile(refver_entity['fasta_path']) as fasta:
#                contigs = fasta.references
#                lengths = fasta.lengths
#            return assemblyspec.get_chrom_converter(contigs, lengths)


##############
# RefverDict #
##############

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
            result.extend(REFVERINFO.get_aliases(standard))
            #result.append(standard)
            #result.extend(REFVERINFO.standard_aliases_map[standard])
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


#############
# chromdict #
#############

class ChromDict(collections.OrderedDict):
    ################
    # initializers #
    ################

    def __init__(self, refver, *args, **kwargs):
        self.refver = refver
        super().__init__(*args, **kwargs)

    @classmethod
    def _init_from_pysamobj(cls, pysamobj, refver=None):
        if refver is None:
            refver = infer_refver_pysamwrapper(pysamobj)
        result = cls(refver)
        for chrom, length in zip(pysamobj.references, pysamobj.lengths):
            result[chrom] = length
        return result

    @classmethod
    def from_fasta_path(cls, fasta_path, refver=None):
        with pysam.FastaFile(fasta_path) as pysamobj:
            return cls._init_from_pysamobj(pysamobj, refver=refver)

    @classmethod
    def from_fasta(cls, fasta, refver=None):
        return cls._init_from_pysamobj(fasta, refver=refver)

    @classmethod
    def from_bam_path(cls, bam_path, refver=None):
        with pysam.AlignmentFile(bam_path) as pysamobj:
            return cls._init_from_pysamobj(pysamobj, refver=refver)

    @classmethod
    def from_bam(cls, bam, refver=None):
        return cls._init_from_pysamobj(bam, refver=refver)

    @classmethod
    def from_bamheader(cls, bamheader, refver=None):
        return cls._init_from_pysamobj(bamheader, refver=refver)

    @classmethod
    def from_vcfheader(cls, vcfheader, refver=None):
        if refver is None:
            refver = infer_refver_vcfheader(vcfheader)
        result = cls(refver)
        for contig in vcfheader.contigs.values():
            result[contig.name] = contig.length
        return result

    @classmethod
    def from_refver(cls, refver):
        fasta_path = REFVERINFO.get_fasta_path(refver)
        with pysam.FastaFile(fasta_path) as fasta:
            return cls.from_fasta(fasta, refver=refver)

    @classmethod
    def from_custom(cls, custom, refver=None):
        """Args:
            custom: {
                'contigs': ['contig1', 'contig2', ...], 
                'lengths': [length1, length2, ...],
            }
        """
        if refver is None:
            refver = infer_refver_base(custom['contigs'], custom['lengths'])
        result = cls(refver)
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
        relevant_chroms = choose_assembled_chroms(self.contigs)
        startswith_chr = [x.startswith('chr') for x in relevant_chroms]
        if all(startswith_chr):
            return True
        elif not any(startswith_chr):
            return False
        else:
            raise Exception(f'Chromosome names are inconsistent of whether prefixed with "chr"')

    #@functools.cached_property
    def get_XY_names(self):
        def helper(name, xy):
            assert xy in ('X', 'Y')
            if name in self.contigs:
                return name
            else:
                raise NoXYError(f'No {xy} chromosome name detected')
            
        if self.is_chr_prefixed:
            X = helper('chrX', 'X')
            Y = helper('chrY', 'Y')
        else:
            X = helper('X', 'X')
            Y = helper('Y', 'Y')

        return X, Y

    @property
    def XY_names(self):
        return self.get_XY_names()

    def check_has_XY(self):
        try:
            X, Y = self.get_XY_names()
        except NoXYError:
            return False
        else:
            return True

    @functools.cached_property
    def assembled_chroms(self):
        return choose_assembled_chroms(self.contigs)

#    def to_gdf(self, assembled_only=True):
#        df = pd.DataFrame({
#            'Chromosome': self.contigs,
#            'Start': 0,
#            'End': self.lengths,
#        })
#        if assembled_only:
#            selector = df['Chromosome'].apply(
#                lambda x: PAT_ASSEMBLED_CHROM.fullmatch(x) is not None
#            )
#            df = df.loc[selector, :]
#
#        return GDF.from_frame(df, refver)

    #def to_gr(self, assembled_only=True, as_gr=True):
    #    gdf = self.to_gdf(assembled_only=assembled_only)
    #    if as_gr:
    #        return gdf.gr
    #    else:
    #        return gdf.df

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


############################
# Assemblyspec and related #
############################

class AssemblySpec:
    """Attributes:
        data: { 
            'length': [249250621, 243199373, ...],
            'role': ['assembled-molecule', ... 'unlocalized-scaffold', ... ],
            'default': ['1', '2', ... 'HSCHR1_RANDOM_CTG5', ... ],
            'ucsc': ['chr1', 'chr2', ... 'chr1_gl000191_random', ... ],
            #'nochr_plus_genbank': ['1', '2', ... 'GL000191.1', ... ],
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
        #self.chromdicts = dict()

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
            #if key2 == 'length':
            #    self.chromdicts[key1] = ChromDict(custom=chromdict_data)

        self.versions = list(
            set(self.data.keys()).difference(['length', 'role'])
        )

    @classmethod
    def from_RefSeq_refver(cls, refseq_refver, species):
        assembly_report_path = libncbicache.get_assembly_report_path(
            refseq_refver, species, force_download=False,
        )
        with open(assembly_report_path, 'r') as infile:
            result = cls.from_assembly_report_file(infile)

        return result

    @classmethod
    def from_assembly_report_file(cls, stream):
        parsed_data = cls.parse_assembly_report(stream)
        result = cls(parsed_data)
        return result

    @classmethod
    def parse_assembly_report(cls, stream):
        data = {
            'length': list(), 
            'default': list(), 
            'ucsc': list(),
            #'nochr_plus_genbank': list(), 
            'role': list(),
            'genbank': list(), 
            'refseq': list(),
        }

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

#                if ucsc is None:
#                    data['nochr_plus_genbank'].append(genbank)
#                elif ucsc == 'chrM':
#                    data['nochr_plus_genbank'].append('MT')
#                else:
#                    mat = PAT_ASSEMBLED_CHROM.fullmatch(ucsc)
#                    if mat is None:
#                        data['nochr_plus_genbank'].append(genbank)
#                    else:
#                        data['nochr_plus_genbank'].append(mat.group(2))

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
            #'nochr_plus_genbank', 
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
        hash_assemblyspec = ChromLengthHash.from_assemblyspec(self)
        hash_inputdata = ChromLengthHash(contigs, lengths)
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

    def get_chrom_converter_from_self(self, key='default'):
        """key was chosen to 'default' since some reference 
        version does not have ucsc names.
        Normalize 'default' names with 'strip_chr=False'
        """

        valid_indexes = tuple((x is not None) for x in self.data[key])
        contigs = list(itertools.compress(self.data[key], valid_indexes))
        lengths = list(itertools.compress(self.data['length'], valid_indexes))

        contigs = [normalize_chrom(x, strip_chr=False) for x in contigs]

        return self.get_chrom_converter(contigs, lengths)


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


class ChromLengthHash:
    def __init__(self, chrom_list, length_list):
        self.lendict = dict(zip(chrom_list, length_list))
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

    #@classmethod
    #def from_RefSeq_refver(cls, refseq_refver):
    #    pass

    #@classmethod
    #def from_refver(cls, refver):
    #    pass

    @classmethod
    def from_fasta(cls, fasta):
        return cls(fasta.references, fasta.lengths)

    @classmethod
    def from_fasta_path(cls, fasta_path):
        with pysam.FastaFile(fasta_path) as fasta:
            return cls.from_fasta(fasta)

    @classmethod
    def from_assemblyspec(cls, assemblyspec):
        assert isinstance(assemblyspec, AssemblySpec)

        used_keys = ('default', 'ucsc', 'genbank', 'refseq')
        chrom_list = list(
            itertools.chain.from_iterable(
                assemblyspec.data[key] for key in used_keys
            )
        )
        length_list = assemblyspec.data['length'] * len(used_keys)
        # remove None chroms
        notnone_indexes = tuple((x is not None) for x in chrom_list)
        chrom_list = list(itertools.compress(chrom_list, notnone_indexes))
        length_list = list(itertools.compress(length_list, notnone_indexes))

        return cls(chrom_list, length_list)

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
        assert isinstance(other, self.__class__)
        length_matches, common_chroms, self_only_chroms, other_only_chroms = (
            self.compare_lendicts(self.normchrom_lendict, other.normchrom_lendict)
        )
        return length_matches

    def check_issubset(self, other):
        assert isinstance(other, self.__class__)
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
        assert isinstance(other, self.__class__)
        common_chroms = set(self.lendict.keys()).intersection(other.lendict.keys())
        return self.get_hash(common_chroms) == other.get_hash(common_chroms)
                
                    


####################
# Setup REFVERINFO #
####################

REFVERINFO = RefverInfo(handygenome.OPTION['refverinfo'], calc_checksum=False)


###########################################################################
# convenience functions for utilizing and modifying the REFVERINFO object #
###########################################################################

def standardize_refver(refver):
    return REFVERINFO.standardize(refver)


standardize = standardize_refver


def deco_standardize(func):
    sig = inspect.signature(func)
    if 'refver' not in sig.parameters.keys():
        raise Exception(
            f'parameter names must include "refver".'
        )

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        ba = sig.bind(*args, **kwargs)
        ba.apply_defaults()
        ba.arguments['refver'] = standardize_refver(ba.arguments['refver'])
        return func(**ba.arguments)

    return wrapper


def add_refver_entity(refver_entity):
    REFVERINFO.add_refver_entity(refver_entity, reset_cache=True)
    

def change_fasta_path(refver, fasta_path):
    REFVERINFO.change_fasta_path(refver, fasta_path)


def list_known_refvers():
    return REFVERINFO.list_known_refvers()


def get_fasta_path(refver):
    return REFVERINFO.get_fasta_path(refver)


def get_fasta(refver):
    return REFVERINFO.get_fasta(refver)


def get_chromdict(refver):
    return REFVERINFO.get_chromdict(refver)


def get_assemblyspec(refver):
    return REFVERINFO.get_assemblyspec(refver)


def get_geneset_handle(refver):
    return REFVERINFO.get_geneset_handle(refver)
    

def get_regulatory_handle(refver):
    return REFVERINFO.get_regulatory_handle(refver)


def get_repeat_handle(refver):
    return REFVERINFO.get_repeat_handle(refver)


def get_datadir(refver):
    return REFVERINFO.get_refver_entity(refver).datadir


def get_RefSeq_standard(refver):
    return REFVERINFO.get_RefSeq_standard(refver)


def get_chrom_converter(refver):
    return REFVERINFO.get_chrom_converter(refver)


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


#def infer_RefSeq_refver(contigs, lengths):
#    query_chromlenhash = ChromLengthHash(contigs, lengths)


def infer_refver_chromdict(chromdict):
    return infer_refver_base(chromdict.contigs, chromdict.lengths)


def infer_refver_vcfheader(vcfheader):
    contigs = list()
    lengths = list()
    for contig in vcfheader.contigs.values():
        contigs.append(contig.name)
        lengths.append(contig.length)
    return infer_refver_base(contigs, lengths)


def infer_refver_vcfpath(vcfpath):
    with pysam.VariantFile(vcfpath) as vcf:
        return infer_refver_vcfheader(vcf.header)


def infer_refver_vr(vr):
    return infer_refver_vcfheader(vr.header)


def infer_refver_fasta(fasta):
    return infer_refver_pysamwrapper(fasta)


def infer_refver_fastapath(fastapath):
    with pysam.FastaFile(fastapath) as fasta:
        return infer_refver_fasta(fasta)


def infer_refver_bamheader(bamheader):
    return infer_refver_pysamwrapper(bamheader)


def infer_refver_bampath(bampath):
    with pysam.AlignmentFile(bampath) as bam:
        return infer_refver_bamheader(bam.header)


def infer_refver_base(contigs, lengths):
    query_hash = ChromLengthHash(contigs, lengths)
    compare_results = list()
    for standard in REFVERINFO.list_standards():
        target_hash = REFVERINFO.get_chromlenhash(standard)
        #if target_hash is None:
        #    continue
        is_subset, other_only_chroms = query_hash.check_issubset(target_hash)
        refverinfo = {
            'standard': standard, 
            'is_subset': is_subset, 
            'other_only_chroms': other_only_chroms,
        }
        compare_results.append(refverinfo)

    subset_refverinfos = [x for x in compare_results if x['is_subset']]
    if len(subset_refverinfos) == 0:
        return None
    else:
        if len(subset_refverinfos) == 1:
            chosen_info = subset_refverinfos[0]
        elif len(subset_refverinfos) > 1:
            minimal_descrepancy_ones = tools.multi_min(
                subset_refverinfos, 
                key=(lambda x: len(x['other_only_chroms'])),
            )
            if len(minimal_descrepancy_ones) == 1:
                chosen_info = minimal_descrepancy_ones[0]
            else:
                raise Exception(f'Cannot break tie: {subset_refverinfos}')

        return chosen_info['standard']


def infer_refver_pysamwrapper(wrapper):
    return infer_refver_base(wrapper.references, wrapper.lengths)



##################################################
# handling of pseudoautosomal region & runs of N #
##################################################


#PARS = RefverDict(
#    # 0-based half-open system
#    # ref: https://www.ncbi.nlm.nih.gov/grc/human
#    {  
#        'GRCh37': {
#            'PAR1_X': ('X', 60_000, 2_699_520),
#            'PAR2_X': ('X', 154_931_043, 155_260_560),
#            'PAR1_Y': ('Y', 10_000, 2_649_520),
#            'PAR2_Y': ('Y', 59_034_049, 59_373_566),
#        }, 
#        'GRCh38': {
#            'PAR1_X': ('X', 10_000, 2_781_479),
#            'PAR2_X': ('X', 155_701_382, 156_030_895),
#            'PAR1_Y': ('Y', 10_000, 2_781_479),
#            'PAR2_Y': ('Y', 56_887_902, 57_217_415),
#        },
#    }
#)
#
#PAR_GRS = dict()
#for key, val in PARS.items():
#    df = pd.DataFrame.from_records(
#        iter(val.values()), 
#        columns=('Chromosome', 'Start', 'End'),
#    )
#    gr = pr.PyRanges(df).sort()
#    PAR_GRS[key] = gr
#PAR_GRS = RefverDict(PAR_GRS)
#
#
#def get_par_gr(refver):
#    chromdict = ChromDict.from_refver(refver)
#    par_gr = PAR_GRS[refver].copy()  # X, Y
#    if chromdict.is_chr_prefixed:
#        par_gr.Chromosome = 'chr' + par_gr.Chromosome
#    return par_gr



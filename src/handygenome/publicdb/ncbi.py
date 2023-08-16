import os
import re
import contextlib
import pickle
import tempfile

import pysam

import handygenome
import handygenome.tools as tools
import handygenome.network as network
import handygenome.logutils as logutils
#import handygenome.publicdb.ncbi_postprocess as ncbi_postprocess
#import handygenome.fastahandler as fastahandler


URL_AUTHORITY = 'ftp.ncbi.nlm.nih.gov'
URL_PATH_REFSEQ_GENOME = '/genomes/refseq'
URL_PATH_DBSNP = "/snp/latest_release/VCF"
REFSEQ_ACC_PATSTR = r'(?P<refseq_acc>(?P<refseq_acc_prefix>[^_]+)_(?P<refseq_acc_main>[^_]+)\.(?P<refseq_acc_version>[^_]+))'
REFVER_PATSTR = r'(?P<refver>(?P<refver_main>[^.]+)(\.(?P<refver_sub>[^.]+))?)'
GENOME_DIR_PAT = re.compile(REFSEQ_ACC_PATSTR + '_' + REFVER_PATSTR)

#SPECIES_PARDIR_MAP = {
#    'Homo_sapiens': 'vertebrate_mammalian',
#    'Mus_musculus': 'vertebrate_mammalian',
#    'Musa_acuminata': 'plant',
#}


class UnavailableSpeciesError(Exception):
    pass


#################################
# RefSeq ftp genome directories #
#################################

class RefseqPaths:
    species_toppaths_pickle = os.path.join(
        handygenome.DIRS['data'], 'refseq_genome_toppaths.pickle',
    )
    all_genome_supersets = [
        #'assembly_summary_refseq.txt',
        #'assembly_summary_refseq_historical.txt',
        #'README.txt',

        'vertebrate_mammalian',
        'vertebrate_other',
        'plant',
        'bacteria',
        'viral',
        'fungi',
        'invertebrate',
        'protozoa',
        'archaea',

        #'unknown',
        #'metagenomes',
    ]

    def __init__(self):
        self.initialize(force_update=False)

    def remove_species_toppaths_pickle(self):
        os.remove(self.__class__.species_toppaths_pickle)

    def initialize(self, force_update=False):
        if force_update:
            self.remove_species_toppaths_pickle()
        self.init_species_toppaths()
        self.refver_toppaths = dict()
        self.latest_refver_toppaths = dict()
        self.refver_details = dict()

    def init_species_toppaths(self):
        if os.path.exists(self.__class__.species_toppaths_pickle):
            self.load_species_toppaths_pickle()
        else:
            self.species_toppaths = dict()

    #####

    def load_species_toppaths_pickle(self):
        with open(self.__class__.species_toppaths_pickle, 'rb') as f:
            self.species_toppaths = pickle.load(f)

    def save_species_toppaths_pickle(self):
        with open(self.__class__.species_toppaths_pickle, 'wb') as f:
            pickle.dump(self.species_toppaths, f)

    def search_for_species_toppath(self, query_species):
        for superset, subdic in self.species_toppaths.items():
            for species in subdic.keys():
                if species == query_species:
                    return subdic[species]
        return None

    def get_species_toppath(self, species):
        toppath = self.search_for_species_toppath(species)
        if toppath is None:  # now query species is not in the data
            for superset in self.__class__.all_genome_supersets:
                success = self.set_species_toppaths(
                    species_superset=superset,
                    retry_count=None, 
                    retry_interval=5, 
                    timeout=180,  # set a large value for "bacteria"
                )
                toppath = self.search_for_species_toppath(species)
                if toppath is None:
                    continue
                else:
                    break

        if toppath is None:
            raise UnavailableSpeciesError(f'Species "{species}" could not be found from the list of RefSeq genomes')

        return toppath

    def set_species_toppaths(
        self, species_superset,
        retry_count=None, retry_interval=5, timeout=60,
    ):
        if species_superset in self.species_toppaths.keys():
            success = True
        else:
            superset_path = network.join_ftp_paths(URL_PATH_REFSEQ_GENOME, species_superset)
            logutils.print_timestamp(f'Fetching from "{superset_path}"')
            try:
                nlst_result = network.ftp_nlst(
                    URL_AUTHORITY, superset_path,
                    retry_count=retry_count, 
                    retry_interval=retry_interval, 
                    timeout=timeout,
                )
            except network.MaxRetryError as exc:
                logutils.print_timestamp(
                    f'Skipping "{superset_path}" because maximum retry count was exceeded'
                )
                success = False
            else:
                # do update
                self.species_toppaths[species_superset] = dict()
                for species_path in nlst_result:
                    species = species_path.split('/')[-1]
                    self.species_toppaths[species_superset][species] = species_path
                success = True
                # save updated dict to pickle file
                self.save_species_toppaths_pickle()

        return success

    #####

    def set_refver_toppaths(
        self, species,
        retry_count=None, retry_interval=5, timeout=60,
    ):
        species_toppath = self.get_species_toppath(species)
        nlst_path = network.join_ftp_paths(species_toppath, 'all_assembly_versions')

        fname_list = network.ftp_nlst(
            URL_AUTHORITY, nlst_path,
            retry_count=retry_count, 
            retry_interval=retry_interval, 
            timeout=timeout,
        )

        data = dict()
        for path in fname_list:
            mat = GENOME_DIR_PAT.fullmatch(os.path.basename(path))
            if mat is None:
                continue
            refver = mat.group('refver_main')
            acc_ver = int(mat.group('refseq_acc_version'))
            data.setdefault(refver, dict())
            data[refver][acc_ver] = path

        self.refver_toppaths[species] = data

    def set_latest_refver_toppaths(
        self, species,
        retry_count=None, retry_interval=5, timeout=60,
    ):
        self.set_refver_toppaths(
            species,
            retry_count=retry_count, 
            retry_interval=retry_interval, 
            timeout=timeout,
        )

        data = dict() 
        for refver, subdic in self.refver_toppaths[species].items():
            data[refver] = max(subdic.items(), key=lambda x: x[0])[1]
        self.latest_refver_toppaths[species] = data

    def get_refver_toppath(
        self, refver, species,
        retry_count=None, retry_interval=5, timeout=60,
    ):
        """Returns: "/~/~/all_assembly_versions"
        """
        if (species not in self.latest_refver_toppaths):
            self.set_latest_refver_toppaths(
                species, 
                retry_count=retry_count, 
                retry_interval=retry_interval, 
                timeout=timeout,
            )
        if refver not in self.latest_refver_toppaths[species]:
            raise Exception(f'There is no refver named "{refver}" in the species "{species}"')
        return self.latest_refver_toppaths[species][refver]

    def get_latest_refver_toppaths(
        self, species,
        retry_count=None, retry_interval=5, timeout=60,
    ):
        if species not in self.latest_refver_toppaths:
            self.set_latest_refver_toppaths(
                species, 
                retry_count=retry_count, 
                retry_interval=retry_interval, 
                timeout=timeout,
            )
        return self.latest_refver_toppaths[species]

    #####

    def set_refver_details(
        self, refver, species,
        retry_count=None, retry_interval=5, timeout=60,
    ):
        refver_toppath = self.get_refver_toppath(
            refver, species,
            retry_count=retry_count, 
            retry_interval=retry_interval, 
            timeout=timeout,
        )
        data = dict()
        for path in network.ftp_nlst(
            URL_AUTHORITY, refver_toppath,
            retry_count=retry_count, 
            retry_interval=retry_interval, 
            timeout=timeout,
        ):
            if path.endswith('assembly_report.txt'):
                data['assembly_report'] = path
            elif (
                path.endswith('genomic.fna.gz')
                and (not path.endswith('rna_from_genomic.fna.gz'))
            ):
                data['fasta'] = path
            elif path.endswith('genomic.gff.gz'):
                data['gff'] = path
            elif path.endswith('protein.faa.gz'):
                data['protein_fasta'] = path
            elif path.endswith('rm.out.gz'):
                data['repeatmasker'] = path

        self.refver_details[refver] = data

    def get_refver_details(
        self, refver, species,
        retry_count=None, retry_interval=5, timeout=60,
    ):
        if (refver not in self.refver_details):
            self.set_refver_details(
                refver, species,
                retry_count=retry_count, 
                retry_interval=retry_interval, 
                timeout=timeout,
            )
        return self.refver_details[refver]

    ###############
    # url getters #
    ###############

    def get_url(self, refver, species, key, force_update):
        if force_update:
            self.initialize(force_update=True)
        path = self.get_refver_details(refver, species)[key]
        return f'https://{URL_AUTHORITY}{path}'

    def get_fasta_url(self, refver, species, force_update=False):
        return self.get_url(refver, species, 'fasta', force_update)

    def get_assembly_report_url(self, refver, species, force_update=False):
        return self.get_url(refver, species, 'assembly_report', force_update)

    def get_gff_url(self, refver, species, force_update=False):
        return self.get_url(refver, species, 'gff', force_update)

    def get_repeatmasker_url(self, refver, species, force_update=False):
        return self.get_url(refver, species, 'repeatmasker', force_update)

    #####################################
    # species/refver availability check #
    #####################################

    def check_species_availability(self, species):
        try:
            toppath = self.get_species_toppath(species)
        except UnavailableSpeciesError:
            return False
        else:
            return True

    def get_available_refvers(
        self, species,
        retry_count=None, retry_interval=5, timeout=60,
    ):
        refver_toppaths = self.get_latest_refver_toppaths(
            species,
            retry_count=retry_count, 
            retry_interval=retry_interval, 
            timeout=timeout,
        )
        return tuple(refver_toppaths.keys())


REFSEQPATHS = RefseqPaths()


# convenience functions

#def get_fasta_url(refver, species, force_update=False):
#    return REFSEQPATHS.get_fasta_url(refver, species, force_update)
#
#
#def get_assembly_report_url(refver, species, force_update=False):
#    return REFSEQPATHS.get_assembly_report_url(refver, species, force_update)
#
#
#def get_gff_url(refver, species, force_update=False):
#    return REFSEQPATHS.get_gff_url(refver, species, force_update)


########################################
# RefSeq file downloaders and fetchers #
########################################

def refseqfile_cacher_base(
    topdir, standard_refver, species, path_basename, force_download, get_url_key,
    download_msg=None,
):
    species_subdir = os.path.join(topdir, species)
    os.makedirs(species_subdir, exist_ok=True)
    filepath = os.path.join(species_subdir, path_basename)

    if (not os.path.exists(filepath)) or force_download:
        did_download = True
        if download_msg is None:
            download_msg = (
                f'Downloading {get_url_key} file for reference version "{standard_refver}"'
            )
        logutils.print_timestamp(download_msg)

        url = REFSEQPATHS.get_url(
            refver=standard_refver, 
            species=species, 
            key=get_url_key, 
            force_update=False,
        )
        network.download(url, filepath)
    else:
        did_download = False

    return filepath, did_download


# fasta

FASTAFILE_DIR = os.path.join(handygenome.DIRS['data'], 'refseq_fasta')
os.makedirs(FASTAFILE_DIR, exist_ok=True)

def get_unedited_fasta_path(standard_refver, species, force_download=False):
    filepath, did_download = refseqfile_cacher_base(
        topdir=FASTAFILE_DIR, 
        standard_refver=standard_refver, 
        species=species, 
        path_basename=f'{standard_refver}.unedited.fna.gz', 
        force_download=force_download, 
        get_url_key='fasta',
    )

    return filepath


def get_edited_fasta_path(standard_refver, species, force_download=False):
    edited_fasta_path = os.path.join(FASTAFILE_DIR, f'{standard_refver}.edited.fasta.gz')
    if not os.path.exists(edited_fasta_path):
        unedited_fasta_path = get_unedited_fasta_path(
            standard_refver, species, force_download=force_download,
        )
        from handygenome.publicdb.ncbi_postprocess import postprocess_fasta
        postprocess_fasta(unedited_fasta_path, edited_fasta_path)
    return edited_fasta_path


# assemblyspec

ASSEMBLYFILE_DIR = os.path.join(handygenome.DIRS['data'], 'assembly_reports')
os.makedirs(ASSEMBLYFILE_DIR, exist_ok=True)

def get_assemblyfile_path(standard_refver, species, force_download=False):
    filepath, did_download = refseqfile_cacher_base(
        topdir=ASSEMBLYFILE_DIR, 
        standard_refver=standard_refver, 
        species=species, 
        path_basename=f'{standard_refver}.txt', 
        force_download=force_download, 
        get_url_key='assembly_report',
    )

    return filepath


# repeatmasker

REPEATMASKER_DIR = os.path.join(handygenome.DIRS['data'], 'refseq_repeatmasker')
os.makedirs(REPEATMASKER_DIR, exist_ok=True)

def get_repeatmasker_path(standard_refver, species, force_download=False):
    filepath, did_download = refseqfile_cacher_base(
        topdir=REPEATMASKER_DIR, 
        standard_refver=standard_refver, 
        species=species, 
        path_basename=f'{standard_refver}.rm.out.gz', 
        force_download=force_download, 
        get_url_key='repeatmasker',
    )

    return filepath


# geneset

GENESET_DIR = os.path.join(handygenome.DIRS['data'], 'refseq_geneset')
os.makedirs(GENESET_DIR, exist_ok=True)

def get_geneset_path(standard_refver, species, force_download=False):
    filepath, did_download = refseqfile_cacher_base(
        topdir=GENESET_DIR, 
        standard_refver=standard_refver, 
        species=species, 
        path_basename=f'{standard_refver}.unedited.gff.gz', 
        force_download=force_download, 
        get_url_key='gff',
    )

    return filepath





#########################################


#def genome_path_sorting(
#    species_genome_topdir,
#    retry_count=None, retry_interval=5, timeout=60,
#):
#    fname_list = network.ftp_nlst(
#        URL_AUTHORITY, species_genome_topdir,
#        retry_count=retry_count, 
#        retry_interval=retry_interval, 
#        timeout=timeout,
#    )
#
#    genome_path_dict = dict()
#    for path in fname_list:
#        mat = GENOME_DIR_PAT.fullmatch(os.path.basename(path))
#        if mat is None:
#            continue
#
#        refver = mat.group('refver_main')
#        acc_ver = int(mat.group('refseq_acc_version'))
#        genome_path_dict.setdefault(refver, dict())
#        genome_path_dict[refver][acc_ver] = path
#
#    return genome_path_dict
#
#
#def refver_subver_key(refver_sub):
#    if refver_sub is None:
#        refver_sub_key = -1
#    elif refver_sub.startswith('p'):
#        refver_sub_key = int(refver_sub[1:])
#    else:
#        refver_sub_key = int(refver_sub)
#    return refver_sub_key
#
#
#def pick_latest_genome_path(genome_path_dict):
#    result = dict()
#    for refver, subdic in genome_path_dict.items():
#        result[refver] = max(subdic.items(), key=lambda x: x[0])[1]
#    return result
#
#
#def collect_latest_genome_paths(
#    species_list, 
#    ftp=None, 
#    retry_count=10, 
#    retry_interval=1, 
#    timeout=5,
#):
#    def helper(ftp, species_genome_topdir, result):
#        fname_dict = genome_path_sorting(ftp, species_genome_topdir)
#        latest_genome_paths_part = pick_latest_genome_path(fname_dict)
#        result.update(latest_genome_paths_part)
#
#    # access ftp server and collect file paths
#    result = dict()
#    if ftp is None:
#        ftp = network.ftp_login(
#            URL_AUTHORITY, 
#            retry_count=retry_count, 
#            retry_interval=retry_interval, 
#            timeout=timeout,
#        )
#
#    with contextlib.redirect_stdout('/dev/null'):
#        species_pardir_pairs = [
#            (species, SPECIES_PARDIR_MAP[species])
#            for species in species_list
#        ]
#        for species, pardir in species_pardir_pairs:
#            topdir = URL_PATH_REFSEQ_GENOME + f'/{pardir}/{species}/all_assembly_versions'
#            helper(ftp, topdir, result)
#
#    return result
#
#
#def find_assemblyfile_from_genomepath(ftp, genome_path):
#    fname_list = network.ftp_listdir(ftp, genome_path)
#    results = [x for x in fname_list if x.endswith('assembly_report.txt')]
#    if len(results) != 1:
#        raise Exception(f'Cannot find a unique assembly report file')
#    return results[0]
#
#
#def collect_assemblyfile_urls(retry_count=10, retry_interval=1, timeout=5):
#    ftp = network.ftp_login(
#        URL_AUTHORITY, 
#        retry_count=retry_count, 
#        retry_interval=retry_interval, 
#        timeout=timeout,
#    )
#    latest_genome_paths = collect_latest_genome_paths(ftp=ftp)
#    assemblyfile_paths = {
#        refver: find_assemblyfile_from_genomepath(ftp, genome_path)
#        for refver, genome_path in latest_genome_paths.items()
#    }
#
#    # add prefix
#    assemblyfile_urls = {
#        key: f'https://{URL_AUTHORITY}' + val
#        for key, val in assemblyfile_paths.items()
#    }
#
#    # quit ftp
#    ftp.quit()
#
#    return assemblyfile_urls


def get_dbsnp_urls(retry_count=10, retry_interval=1, timeout=5):
    # access ftp server and collect file paths
    ftp = network.ftp_login(
        URL_AUTHORITY, 
        retry_count=retry_count, 
        retry_interval=retry_interval, 
        timeout=timeout,
    )
    with contextlib.redirect_stdout('/dev/null'):
        fname_list = network.ftp_listdir(ftp, URL_PATH_DBSNP)
        ftp.quit()

    # make result
    pat = re.compile(r'(.+)\.gz(.*)')
    url_dict = dict()
    for fname in fname_list:
        mat = pat.fullmatch(os.path.basename(fname))
        if mat is None:  # CHECKSUM
            continue

        key = mat.group(1)
        url_dict.setdefault(key, dict())

        if mat.group(2) == '':
            subkey = 'vcf'
        elif mat.group(2) == '.md5':
            subkey = 'vcf_md5'
        else:
            continue

        url_dict[key][subkey] = f'https://{URL_AUTHORITY}' + fname

    return url_dict

import os
import gzip
import multiprocessing
import shutil

import handygenome
import handygenome.tools as tools
import handygenome.network as network
import handygenome.publicdb.ncbi as libncbi


CACHEDIR = os.path.join(handygenome.USERDATA_DIR, 'ncbi_files')
os.makedirs(CACHEDIR, exist_ok=True)


def get_refver_dir(RefSeq_refver, species):
    """The top directory of ncbi files for a specific RefSeq reference version"""
    species_dir = os.path.join(CACHEDIR, species)
    os.makedirs(species_dir, exist_ok=True)

    refver_dir = os.path.join(species_dir, RefSeq_refver)
    os.makedirs(refver_dir, exist_ok=True)

    return refver_dir


#########
# fasta #
#########

def fasta_path_string(RefSeq_refver, species):
    return os.path.join(get_refver_dir(RefSeq_refver, species), 'genome_fasta')


def get_fasta_path(RefSeq_refver, species, force_download=False):
    fasta_path = fasta_path_string(RefSeq_refver, species)
    if force_download or (not os.path.exists(fasta_path)):
        fasta_url = libncbi.GENOME_PATHS.get_fasta_url(RefSeq_refver, species)
        network.download(fasta_url, fasta_path)
    return fasta_path


###################
# assembly_report #
###################

def assembly_report_path_string(RefSeq_refver, species):
    return os.path.join(get_refver_dir(RefSeq_refver, species), 'assembly_report')


def get_assembly_report_path(RefSeq_refver, species, force_download=False):
    assembly_report_path = assembly_report_path_string(RefSeq_refver, species)
    if force_download or (not os.path.exists(assembly_report_path)):
        assembly_report_url = libncbi.GENOME_PATHS.get_assembly_report_url(
            RefSeq_refver, species, force_update=force_download,
        )
        network.download(assembly_report_url, assembly_report_path)
    return assembly_report_path


####################
# assembly_regions #
####################

def assembly_regions_path_string(RefSeq_refver, species):
    return os.path.join(get_refver_dir(RefSeq_refver, species), 'assembly_regions')


def get_assembly_regions_path(RefSeq_refver, species, force_download=False):
    assembly_regions_path = assembly_regions_path_string(RefSeq_refver, species)
    if force_download or (not os.path.exists(assembly_regions_path)):
        assembly_regions_url = libncbi.GENOME_PATHS.get_assembly_regions_url(
            RefSeq_refver, species, force_update=force_download,
        )
        network.download(assembly_regions_url, assembly_regions_path)
    return assembly_regions_path


#######
# gff #
#######

def gff_path_string(RefSeq_refver, species):
    return os.path.join(get_refver_dir(RefSeq_refver, species), 'gff')


def get_gff_path(RefSeq_refver, species, force_download=False):
    gff_path = gff_path_string(RefSeq_refver, species)
    if force_download or (not os.path.exists(gff_path)):
        gff_url = libncbi.GENOME_PATHS.get_gff_url(RefSeq_refver, species)
        network.download(gff_url, gff_path)
    return gff_path


#######################
# repeatmasker output #
#######################

def repeatmasker_path_string(RefSeq_refver, species):
    return os.path.join(get_refver_dir(RefSeq_refver, species), 'repeatmasker')


def get_repeatmasker_path(RefSeq_refver, species, force_download=False):
    repeatmasker_path = repeatmasker_path_string(RefSeq_refver, species)
    if force_download or (not os.path.exists(repeatmasker_path)):
        repeatmasker_url = libncbi.GENOME_PATHS.get_repeatmasker_url(RefSeq_refver, species)
        network.download(repeatmasker_url, repeatmasker_path)
    return repeatmasker_path




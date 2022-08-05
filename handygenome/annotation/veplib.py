# Only deals with command-line VEP
import os
import subprocess

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
initvcf = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'initvcf']))


VEP_INFO_FIELD = 'CSQ'
DEFAULT_DISTANCE = 5000
REFVER_TO_VEPARGS = {
    'hg19' : {'species': 'homo_sapiens', 'assembly': 'GRCh37'},
    'hg38' : {'species': 'homo_sapiens', 'assembly': 'GRCh38'},
    'mm10' : {'species': 'mus_musculus', 'assembly': 'GRCm38'},
    'mm39' : {'species': 'mus_musculus', 'assembly': 'GRCm39'},
    }

ALLOWED_SPECIES_ASSEMBLIES = {
    'homo_sapiens': ('GRCh37', 'GRCh38'),
    'mus_musculus': ('GRCm38', 'GRCm39'),
    }


def check_species_assembly_sanity(species, assembly):
    if species not in ALLOWED_SPECIES_ASSEMBLIES.keys():
        raise Exception(f'"species" must be one of '
                        f'{ALLOWED_SPECIES_ASSEMBLIES.keys()}.')
    else:
        if assembly not in ALLOWED_SPECIES_ASSEMBLIES[species]:
            raise Exception(
                f'For the species {species}, assembly name must be one of '
                f'{ALLOWED_SPECIES_ASSEMBLIES[species]}.')


def get_vep_args(infile_path, outfile_path, fasta_path, species, assembly,
                 distance=DEFAULT_DISTANCE, force_overwrite=True):
    check_species_assembly_sanity(species, assembly)

    if species == 'mus_musculus' and assembly == 'GRCm38':
        vep_path = common.VEP_MM10
    else:
        vep_path = common.VEP

    with pysam.FastaFile(fasta_path) as fasta:
        max_chrom_size = max(fasta.lengths)

    vep_args = [
        vep_path,
        '-i', f'{infile_path}',
        '-o', f'{outfile_path}',
        '--fasta', f'{fasta_path}',
        '--species', f'{species}',
        '--assembly', f'{assembly}',
        '--max_sv_size', f'{max_chrom_size}',

        '--dir', f'{common.VEP_CACHE_DIR}',
        '--cache',
        '--offline',

        '--vcf',
        '--vcf_info_field', f'{VEP_INFO_FIELD}',
        '--terms', 'SO',
        '--shift_hgvs', '0',
        '--dont_skip',
        '--no_stats',

        '--distance', f'{distance}',

        '--symbol',
        '--numbers',
        '--hgvsg',
        '--hgvs',
        '--ccds',
        #'--mirna',
        '--biotype',
        '--canonical',
        '--xref_refseq',
        '--mane',

        '--regulatory',

        '--protein',
        '--uniprot',
        '--sift', 'b',
        '--polyphen', 'b',

        #'--check_existing',
        #'--pubmed',
        #'--af',

        #'--max_af',
        #'--af_1kg',
        #'--af_esp',
        #'--af_gnomad',
        #'--var_synonyms',

        #'--nearest', 'transcript',
    ]

    if force_overwrite:
        vep_args.append('--force_overwrite')

    return vep_args


def run_vep_with_vcfspec(vcfspec, refver, distance=DEFAULT_DISTANCE):
    chromdict = common.ChromDict(refver = refver)
    vr = initvcf.create_vr(chromdict=chromdict, vcfspec=vcfspec)

    return run_vep_with_vr(vr, refver, distance)


def run_vep_with_interval(interval, refver, distance=DEFAULT_DISTANCE):
    chromdict = common.ChromDict(refver = refver)
    vcfspec = common.Vcfspec(interval.chrom, interval.start0, 'N', ['<DEL>'])
    while True:
        vr = initvcf.create_vr(chromdict=chromdict, vcfspec=vcfspec, 
                               end=interval.end1)
        if vr.stop == interval.end1:
            break
        else:
            print(f'INFO/END of created VariantRecord object is faulty. '
                  f'Retrying.')
            continue

    return run_vep_with_vr(vr, refver, distance)


def run_vep_with_vr(vr, refver=None, distance=DEFAULT_DISTANCE):
    if refver is None:
        refver = common.infer_refver_vr(vr)

    assert refver in REFVER_TO_VEPARGS.keys()

    species = REFVER_TO_VEPARGS[refver]['species']
    assembly = REFVER_TO_VEPARGS[refver]['assembly']

    fasta_path = common.DEFAULT_FASTA_PATHS[refver]
    infile_path = workflow.get_tmpfile_path(delete=True)
    outfile_path = infile_path + '.vep'

    with pysam.VariantFile(infile_path, mode='wz', 
                           header=vr.header) as out_vcf:
        out_vcf.write(vr)

    vep_args = get_vep_args(infile_path, outfile_path, fasta_path, species, 
                            assembly, distance=distance)

    p = subprocess.run(args=vep_args, capture_output=True, text=True,
                       check=True)

    with pysam.VariantFile(outfile_path) as in_vcf:
        out_vr = next(in_vcf.fetch())

    os.remove(infile_path)
    os.remove(outfile_path)

    return out_vr



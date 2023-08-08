import os
import shutil
import subprocess
import tempfile
import re

import pysam

import handygenome
import handygenome.workflow as workflow
import handygenome.variant.variantplus as variantplus
import handygenome.bameditor as bameditor


GATK = os.path.join(handygenome.DIRS['externals'], 'gatk_wrapper.sh')


def run_haplotypecaller(fasta_path, infile_path, outfile_path, tmpdir=None, 
                        incl_region_path=None, excl_region_path=None,
                        mem_gb=6, threads=4, mbq=10, rm_tmpdir=True):
    """Args:
        threads: value for --native-pair-hmm-threads option. 
            gatk default is 4.
    """
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp(prefix='HaplotypeCaller_tmpdir_', dir=os.getcwd())

    args = [
        GATK, 
        '--java-options', f'-Xmx{mem_gb}G',
        'HaplotypeCaller',
        '--reference', fasta_path,
        '--input', infile_path,
        '--output', outfile_path,
        '--native-pair-hmm-threads', str(threads),
        '--min-base-quality-score', str(mbq),
        '--sample-ploidy', '2', # gatk default
        '--active-probability-threshold', '0.002', # gatk default
        '--standard-min-confidence-threshold-for-calling', '30', # gatk default
        '--dont-use-soft-clipped-bases', 'false',
        '--create-output-variant-index', 'false',
        '--tmp-dir', tmpdir,
        ]

    if incl_region_path is not None:
        args.extend(['--intervals', incl_region_path])
    if excl_region_path is not None:
        args.extend(['--exclude-intervals', excl_region_path])

    p = subprocess.run(args, capture_output=True, text=True, check=False, shell=False) 
    if rm_tmpdir:
        shutil.rmtree(tmpdir)

    return p


def run_Mutect2(infile_path_list, outfile_path, fasta_path, tmpdir=None, realigned_bam_path=None, incl_region_path=None, excl_region_path=None, mem_gb=12, threads=4, mbq=10, rm_tmpdir=True):
    """For 4.3.0.0"""

    if tmpdir is None:
        tmpdir = tempfile.mkdtemp(prefix='Mutect2_tmpdir_', dir=os.getcwd())

    args = [
        GATK, 
        '--java-options', f'-Xmx{mem_gb}G',
        'Mutect2',
        '--reference', fasta_path,
        #'--input', infile_path,
        '--output', outfile_path,

        '--force-active', 'true',

        '--min-base-quality-score', str(mbq),
        '--native-pair-hmm-threads', str(threads),
        '--active-probability-threshold', '0.002', # gatk default
        '--dont-use-soft-clipped-bases', 'false',
        '--create-output-variant-index', 'false',
        '--tmp-dir', tmpdir,

        '--min-assembly-region-size', f'{50}',  # gatk default
        '--assembly-region-padding', f'{100}',
        #'--dont-trim-active-regions', 'true',
    ]

    for infile_path in infile_path_list:
        args.extend(['--input', infile_path])

    if realigned_bam_path is not None:
        args.extend(['--bam-output', realigned_bam_path])

    if incl_region_path is not None:
        args.extend(['--intervals', incl_region_path])
    if excl_region_path is not None:
        args.extend(['--exclude-intervals', excl_region_path])

    p = subprocess.run(args, capture_output=True, text=True, check=False, shell=False) 
    if rm_tmpdir:
        shutil.rmtree(tmpdir)

    return p


def realign_with_Mutect2(bam_path_list, chrom, start0, end0, fasta_path, realigned_bam_path):
    tmpdir = tempfile.mkdtemp(prefix='tmpdir_Mutect2_realignment_', dir=os.getcwd())
    out_vcf_path = os.path.join(tmpdir, 'output.vcf.gz')

    region_path = os.path.join(tmpdir, 'region.bed')
    with open(region_path, 'wt') as infile:
        infile.write(f'{chrom}\t{start0}\t{end0}\n')

    p = run_Mutect2(bam_path_list, out_vcf_path, fasta_path, tmpdir=tmpdir, realigned_bam_path=realigned_bam_path, incl_region_path=region_path, rm_tmpdir=False)
    if p.returncode != 0:
        raise Exception(p.stderr)

    result_vplist = variantplus.load_vcf(out_vcf_path)

    shutil.rmtree(tmpdir)

    return result_vplist


def split_realigned_bam(realigned_bam_path):
    in_bam = pysam.AlignmentFile(realigned_bam_path, mode='rb')
    readgroups = [x['ID'] for x in in_bam.header['RG'] if x['ID'] != 'ArtificialHaplotypeRG']
    out_bam_paths = {
        rg: re.sub('.bam$', f'.{rg}.bam', realigned_bam_path)
        for rg in readgroups
    }
    out_bams = {
        rg: pysam.AlignmentFile(out_bam_paths[rg], mode='wb', header=in_bam.header)
        for rg in readgroups
    }

    for read in in_bam.fetch():
        rg_tag = read.get_tag('RG')
        if rg_tag != 'ArtificialHaplotypeRG':
            out_bams[read.get_tag('RG')].write(read)

    for bam in out_bams.values():
        bam.close()
    for bam_path in out_bam_paths.values():
        bameditor.sort_and_index(bam_path)

    split_bams = {
        rg: pysam.AlignmentFile(out_bam_paths[rg], mode='rb')
        for rg in readgroups
    }
    return split_bams

import os
import itertools

import pysam

import handygenome.interval as libinterval
import handygenome.workflow as workflow
import handygenome.logutils as logutils
import handygenome.vcfeditor.misc as vcfmisc
import handygenome.deco as deco


LOGGER = logutils.get_logger()


def sanity_check(vcf_path, outdir):
    vcf_path = workflow.arghandler_infile(vcf_path)
    outdir = workflow.arghandler_outdir(outdir)

    return vcf_path, outdir


@deco.get_deco_num_set_differently(('n_file', 'n_line'), 1)
def get_output_lineno_list(total_lineno, n_file=None, n_line=None):
    """
    Args:
        total_lineno: Total number of variant records of input vcf.
        n_file: Number of output files. If greater than the total line number, 
            reduced to the total line number.
        n_line: Number of variant records per one output file. If greater 
            than the total line number, reduced to the total line number.
        * One and only one of n_file or n_line argument must be set.

    Returns:
        A list composed of the number of variant records per output file.
    """

    def warn():
        LOGGER.warning(
            f'"n_file" is greater than "total_lineno". It will be '
            f'changed to be the same with "total_lineno".')
    
    assert total_lineno > 0, f'"total_lineno" must be a positive integer.'

    if n_file is not None:
        if n_file > total_lineno:
            warn()
        result = libinterval.get_interval_lengths_num(total_lineno, n_file)
    elif n_line is not None:
        if n_line > total_lineno:
            warn()
        result = libinterval.get_interval_lengths_width(total_lineno, n_line)

    return result


def write_split_vcfs(vcf_path, output_lineno_list, split_filenames, 
                     mode_pysam):
    with pysam.VariantFile(vcf_path, 'r') as in_vcf:
        fetch = in_vcf.fetch()
        for lineno, outfile_path in zip(output_lineno_list, split_filenames):
            with pysam.VariantFile(outfile_path, mode_pysam, 
                                   header=in_vcf.header) as out_vcf:
                for _ in range(lineno):
                    out_vcf.write(next(fetch))


@deco.get_deco_num_set_differently(('n_file', 'n_line'), 1)
def main(vcf_path, outdir, n_file=None, n_line=None, mode_bcftools='z', 
         mode_pysam=None, prefix='', suffix='.vcf.gz'):
    """
    Splits input vcf by either 'n_file' (number of output files) or 
        'n_line' (number of lines each split file contains).
    Split file names looks like: <prefix>000<suffix>

    Args:
        vcf_path: Input vcf file path
        outdir: Directory where split vcf files will be saved. May not be 
            existing. Must not be an non-empty directory.
        One of 'n_file' and 'n_line' must be set.
    """

    vcf_path, outdir = sanity_check(vcf_path, outdir)

    total_lineno = vcfmisc.get_vcf_lineno(vcf_path)
    output_lineno_list = get_output_lineno_list(
            total_lineno, n_file = n_file, n_line = n_line,
            )
    mode_pysam = vcfmisc.write_mode_arghandler(mode_bcftools, mode_pysam)
    split_filenames = workflow.get_split_filenames( 
        len(output_lineno_list), outdir, prefix, suffix)
        # "split_filenames" begins with 0

    write_split_vcfs(vcf_path, output_lineno_list, split_filenames, 
                     mode_pysam)

    return split_filenames

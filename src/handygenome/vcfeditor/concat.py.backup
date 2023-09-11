import re
import os
import argparse

import pysam

import handygenome.vcfeditor.misc as vcfmisc
import handygenome.workflow as workflow


def sanity_check_args(infile_path_list, outfile_path, outfile_must_not_exist):
    infile_path_list = [workflow.arghandler_infile(x) 
                        for x in infile_path_list]
    arghandler = workflow.get_arghandler_outfile(outfile_must_not_exist)
    outfile_path = arghandler(outfile_path)

    return infile_path_list, outfile_path


def sanity_check_headers(infile_path_list):
    """Checks if all headers are equal"""

    samples_list = list()
    infokeys_list = list()
    formatkeys_list = list()

    for vcf_path in infile_path_list:
        with pysam.VariantFile(vcf_path, 'r') as in_vcf:
            samples_list.append(tuple(in_vcf.header.samples))
            infokeys_list.append(tuple(sorted(in_vcf.header.info.keys())))
            formatkeys_list.append(tuple(sorted(in_vcf.header.formats.keys())))

    # check samples
    if len(set(samples_list)) != 1:
        raise Exception(f'Headers of input vcf files are different in '
                        'samples (orders must match)')
    elif len(set(infokeys_list)) != 1:
        raise Exception(f'Headers of input vcf files are different in '
                        'INFO keys')
    elif len(set(formatkeys_list)) != 1:
        raise Exception(f'Headers of input vcf files are different in '
                        'FORMAT keys')


def write_outfile(infile_path_list, outfile_path, mode_pysam):
    with pysam.VariantFile(infile_path_list[0]) as in_vcf:
        outfile_header = in_vcf.header.copy()

    with pysam.VariantFile(outfile_path, mode=mode_pysam, 
                           header=outfile_header) as out_vcf:
        for vcf_path in infile_path_list:
            with pysam.VariantFile(vcf_path, 'r') as in_vcf:
                for vr in in_vcf.fetch():
                    out_vcf.write(vr)


def main(
        infile_path_list,
        outfile_path,
        mode_bcftools='z', 
        mode_pysam=None, 
        outfile_must_not_exist='ask',
        ):
    infile_path_list, outfile_path = sanity_check_args(
        infile_path_list, outfile_path, outfile_must_not_exist)
    sanity_check_headers(infile_path_list)

    mode_pysam = vcfmisc.write_mode_arghandler(mode_bcftools, mode_pysam)
    write_outfile(infile_path_list, outfile_path, mode_pysam)

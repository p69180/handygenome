import re
import os
import argparse

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
toolsetup = importlib.import_module('.'.join([top_package_name, 'workflow', 'toolsetup']))
alleleindexes = importlib.import_module('.'.join([top_package_name, 'variant', 'alleleindexes']))
varianthandler = importlib.import_module('.'.join([top_package_name, 'variant', 'varianthandler']))
indexing = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'indexing']))


MODE_CHOICES = ('unpaired',)


def argument_parser(cmdargs):
    def sanity_check(args):
        pass

    def postprocess(args):
        pass

    parser_dict = workflow.init_parser(
        description=(
            f'Annotate SomaticIndex and GermlineIndexes INFO values.\n'
            f'Modes:\n'
            f'    unpaired: For tumor sample without paired normal sample. '
            f'For all variants, SomaticIndex is set to 1. GermlineIndexes '
            f'are not modified.'
            ))

    # required
    workflow.add_infile_arg(parser_dict['required'], required=True)
    workflow.add_outfile_arg(parser_dict['required'], required=True, 
                             must_not_exist='ask')
    parser_dict['required'].add_argument(
        '--mode', dest='mode', required=True, choices=MODE_CHOICES,
        help=f'Must be one of {MODE_CHOICES}.')

    # flag

    # others
    workflow.add_logging_args(parser_dict)
    workflow.add_index_arg(parser_dict)
    workflow.add_outfmt_arg(parser_dict['optional'], required=False)

    args = parser_dict['main'].parse_args(cmdargs)
    sanity_check(args)
    postprocess(args)

    return args


def job_unpaired(infile_path, outfile_path, mode_pysam):
    somidx = alleleindexes.SomaticIndex(index=1)

    in_vcf = pysam.VariantFile(infile_path, 'r')
    out_vcf_hdr = in_vcf.header.copy()
    somidx.write_meta(out_vcf_hdr)
    out_vcf = pysam.VariantFile(outfile_path, mode=mode_pysam, 
                                header=out_vcf_hdr)

    for vr in in_vcf.fetch():
        new_vr = varianthandler.reheader(vr, out_vcf_hdr)
        somidx.write(new_vr)
        out_vcf.write(new_vr)

    out_vcf.close()
    in_vcf.close()


def main(cmdargs):
    args = argument_parser(cmdargs)
    logger = toolsetup.setup_logger(args)

    # main
    if args.mode == 'unpaired':
        job_unpaired(args.infile_path, args.outfile_path, args.mode_pysam)

    if not args.donot_index:
        indexing.index_vcf(args.outfile_path)



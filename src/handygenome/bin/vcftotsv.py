import sys
import re
import os
import gzip

import pysam
import pandas as pd

import handygenome.logutils as logutils
import handygenome.refgenome.refgenome as refgenome
import handygenome.workflow as workflow
import handygenome.workflow.shell_wrapper as shell_wrapper
import handygenome.utils.workflow_utils as workflow_utils
from handygenome.utils.workflow_utils import MultiArgsList
from handygenome.variant.vcfdataframe import VCFDataFrame
from handygenome.variant.variantplus import VariantPlusList, VariantPlus


##############
# argparsing #
##############

def argparsing_postprocess(args):
    # sanitycheck
    if not re.fullmatch(r'.+\.vcf(\.gz)?$', args.infile_path):
        raise Exception(f'Input file must end with ".vcf" or ".vcf.gz".')

    if args.outfile_path is None:
        args.outfile_path = workflow.arghandler_outfile_ask(
            re.sub(r'\.vcf(\.gz)?$', '.tsv.gz', args.infile_path)
        )
    else:
        if not re.fullmatch(r'.+\.tsv\.gz$', args.outfile_path):
            raise Exception(f'Output file must end with ".tsv.gz".')


def argparsing(cmdargs=None):
    parser_dict = workflow.init_parser(
        description=f'Convert a vcf into tsv format.',
    )

    # required #
    workflow.add_infile_arg(parser_dict['required'], required=True)

    # optional #
    workflow.add_outfile_arg(
        parser_dict['optional'], 
        required=False, 
        must_not_exist='ask',
        help='Output file path. Must end with ".tsv.gz". If not set, input file path with "vcf" substituted with "tsv.gz" will be used.',
    )

    # main #
    args = parser_dict['main'].parse_args(cmdargs)
    argparsing_postprocess(args)

    return args

########
# main #
########

def main_core(cmdargs):
    args = argparsing(cmdargs)

    in_vcf = pysam.VariantFile(args.infile_path)
    all_samples = list(in_vcf.header.samples)
    vplist = VariantPlusList.from_vcf_lazy(args.infile_path)
    vp_iterator = vplist.get_vp_iter_from_vcf()
    vr_iterator = in_vcf.fetch()
    refver = refgenome.infer_refver_vcfheader(in_vcf.header)
    vp0 = VariantPlus.from_vr(next(in_vcf.fetch()))

    # args postprocess
    (
        sampleids, 
        pon_samples, 
        nonpon_samples,
    ) = VCFDataFrame.postprocess_sampleids(
        sampleids=None, 
        all_samples=all_samples, 
        vp0=vp0,
    )
    n_allele = 2

    # main
    mode = 'samplewise_1'
    allele_columns = VCFDataFrame.get_allele_columns(n_allele)
    all_columns = VCFDataFrame.make_columns(
        mode=mode, 
        allele_columns=allele_columns, 
        pon_samples=pon_samples, 
        nonpon_samples=nonpon_samples,
    )
    outfile = gzip.open(args.outfile_path, 'wt')
    row_dict_iterator = VCFDataFrame.iter_row_dict(
        vp_iterator=vp_iterator,
        mode=mode,
        allele_columns=allele_columns,
        sampleids=sampleids,
        pon_samples=pon_samples,
        nonpon_samples=nonpon_samples,
    )
    for idx, row_dict in enumerate(row_dict_iterator):
        row_df = pd.DataFrame.from_records(
            [row_dict], 
            columns=all_columns,
        )
        row_df.to_csv(
            outfile, 
            sep='\t', na_rep='NA',
            header=(idx == 0), index=False, 
        )

    outfile.close()
    in_vcf.close()

    logutils.log(f'All successfully finished.', level='info')


def main():
    shell_wrapper.run_func(
        main_core,
        args=(sys.argv[1:],),
        use_condaenv=True,
        remove_funcrun_tmpfile_dir=True,
        remove_condawrap_script_dir=True,
        raise_with_failure=False,
    )



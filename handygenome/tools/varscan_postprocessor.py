import pysam

import handygenome.common as common
import handygenome.workflow as workflow
import handygenome.workflow.toolsetup as toolsetup
import handygenome.variant.svcaller_parser as svcaller_parser
import handygenome.variant.varianthandler as varianthandler
import handygenome.vcfeditor.headerhandler as headerhandler
import handygenome.vcfeditor.indexing as indexing
import handygenome.variant.vcfspec as libvcfspec
import handygenome.vcfeditor.varscan_editing as varscan_editing


def argument_parser(cmdargs):
    def sanity_check(args):
        pass

    parser_dict = workflow.init_parser(
        description=textwrap.dedent(f"""\
            Modifies VarScan2 output file:
            - REF and ALT strings are made VCF-compatible
            - Multi-ALT lines are not split
            - non-SOMATIC, non-PASS records are removed by default
            - 
            - Original annotations are removed.

            - Split multiallelic lines.
            - Change SV caller records into breakends form.
            - Normalize non-SV records into the leftmost form.
            - Normalize SV records into the bnd1-advanced form."""
        )
    )

    # required
    parser_dict['required'].add_argument(
        '--snps',
        dest='infile_snv',
        required=True,
        help=f'The result file named "*.snps.vcf.gz"',
    )
    parser_dict['required'].add_argument(
        '--indels',
        dest='infile_indel',
        required=True,
        help=f'The result file named "*.indels.vcf.gz"',
    )

    workflow.add_outfile_arg(parser_dict['required'], required=True, must_not_exist='ask')
    workflow.add_refver_arg(parser_dict['required'], required=True, choices='all')

    # optional
    workflow.add_outfmt_arg(parser_dict['optional'], required=False)
    workflow.add_logging_args(parser_dict)

    # flag
    workflow.add_index_arg(parser_dict)
    parser_dict['flag'].add_argument(
        '--include-non-somatic',
        dest='include_non_somatic',
        action='store_true',
        help=f'If set, records in which INFO/SOMATIC flag is false are included. They are removed by default.',
    )
    parser_dict['flag'].add_argument(
        '--include-non-pass',
        dest='include_non_pass',
        action='store_true',
        help=f'If set, records in which FILTER value is not PASS are included. They are removed by default.',
    )

    args = parser_dict['main'].parse_args(cmdargs)
    sanity_check(args)

    return args


def main(cmdargs):
    args = argument_parser(cmdargs)
    logger = toolsetup.setup_logger(args)

    logger.info('BEGINNING')

    varscan_editing.modify_varscan_vcf(
        infile_path_list=[args.infile_snv, args.infile_indel], 
        outfile_path=args.outfile_path, 
        refver=args.refver, 
        mode_pysam=args.mode_pysam, 
        somatic_only=(not args.include_non_somatic), 
        pass_only=(not args.include_non_pass),
        index=(not args.donot_index),
    )

    logger.info('ALL SUCCESSFULLY FINISHED')

import re
import os
import argparse
import textwrap

import pysam

import handygenome.workflow as workflow
import handygenome.workflow.toolsetup as toolsetup
import handygenome.vcfeditor.merge as libmerge
import handygenome.vcfeditor.indexing as indexing


def argument_parser(cmdargs):
    def sanity_check(args):
        e_msg = ('Invalid --isec-indices value. Please refer to the help '
                 'message.')
        if args.isec_indices is not None:
            if len(args.isec_indices) != len(args.infile_path_list):
                raise Exception(e_msg)
            elif not set(isec_indices).issubset({'0', '1'}):
                raise Exception(e_msg)

    def postprocess(args):
        if args.isec_indices is not None:
            args.isec_indices = [int(x) for x in args.isec_indices]

    parser_dict = workflow.init_parser()

    # required
    workflow.add_infilelist_arg(parser_dict['required'], required=True)
    workflow.add_outfile_arg(parser_dict['required'], required=True, 
                             must_not_exist='ask')
    workflow.add_fasta_arg(parser_dict['required'], required=True)
    parser_dict['required'].add_argument(
        '--mode', required=True, choices=('isec', 'union'),
        help=f'Must be "isec" (which means intersection) or "union".')

    # optional
    parser_dict['optional'].add_argument(
        '--isec-indices', dest='isec_indices', required=False,
        help=textwrap.dedent(f"""\
            (Only applied for intersection) A string composed of 0 or 1,
            with the same length as the number of input files. Files 
            marked with 0 are excluded and those with 1 are included. 
            If not set, all samples are included in intersection."""))

    # flag
    parser_dict['flag'].add_argument(
        '--remove-infoformat', dest='remove_infoformat', action='store_true',
        help=f'If set, all INFO and FORMAT data are removed from the output.')

    # others
    workflow.add_logging_args(parser_dict)
    workflow.add_index_arg(parser_dict)
    workflow.add_outfmt_arg(parser_dict['optional'], required=False)

    args = parser_dict['main'].parse_args(cmdargs)
    sanity_check(args)
    postprocess(args)

    return args


def main(cmdargs):
    args = argument_parser(cmdargs)
    logger = toolsetup.setup_logger(args)

    libmerge.main_file(
        infile_path_list=args.infile_path_list,
        outfile_path=args.outfile_path, fasta_path=args.fasta_path,
        remove_infoformat=args.remove_infoformat,
        isec=(args.mode == 'isec'), isec_indices=args.isec_indices,
        mode_pysam=args.mode_pysam, outfile_must_not_exist='no',
        logger=logger)

    if not args.donot_index:
        indexing.index_vcf(args.outfile_path)

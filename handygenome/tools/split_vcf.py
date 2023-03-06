import textwrap

import handygenome.common as common
import handygenome.workflow as workflow
import handygenome.vcfeditor.split as split_module
import handygenome.vcfeditor.indexing as indexing


def argument_parser(cmdargs):
    def sanity_check(args):
        if (args.n_line, args.n_file).count(None) != 1:
            e_msg = f"Exactly one of '--nfile' or '--nline' must be used."
            raise SystemExit(e_msg)  # this does not show traceback

    parser_dict = workflow.init_parser(description=textwrap.dedent("""\
        - Splits input vcf by either 'nfile' (number of output files) or 
          'nline' (number of variant record lines each split file contains).
        - Split file names looks like: 
          <prefix>000<suffix>, <prefix>001<suffix>, ...
        - Exactly one of '--nfile' or '--nline' options must be used.
        """))

    workflow.add_infile_arg(parser_dict['required'], required=True)
    workflow.add_outdir_arg(parser_dict['required'], required=True)
    workflow.add_outfmt_arg(parser_dict['optional'], required=False, 
                            default='z')
    workflow.add_index_arg(parser_dict)

    parser_dict['optional'].add_argument(
        '--nfile', dest='n_file', required=False, type=int,
        help=textwrap.dedent(f"""\
            The number of output files. 
            Must not be used at the same time with --nline option."""))
    parser_dict['optional'].add_argument(
        '--nline', dest='n_line', required=False, type=int,
        help=textwrap.dedent(f"""\
            The number of variant record lines per output file.
            Must not be used at the same time with --nfile option."""))
    parser_dict['optional'].add_argument(
        '--prefix', dest='prefix', required=False, default='',
        help=textwrap.dedent(f"""\
            prefix for the names of split files. Default is an empty
            string."""))
    parser_dict['optional'].add_argument(
        '--suffix', dest='suffix', required=False, default='.vcf.gz',
        help=f'suffix for the names of split files.')

    args = parser_dict['main'].parse_args(cmdargs)
    sanity_check(args)

    return args


def main(cmdargs):
    args = argument_parser(cmdargs)
    split_module.main(vcf_path=args.infile_path, outdir=args.outdir_path,
                      n_file=args.n_file, n_line=args.n_line,
                      mode_pysam=args.mode_pysam, prefix=args.prefix,
                      suffix=args.suffix)

    if not args.donot_index:
        indexing.index_vcf(args.outfile_path)

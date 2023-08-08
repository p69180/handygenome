import handygenome.workflow as workflow
#import handygenome.workflow.toolsetup as toolsetup
import handygenome.dbconv.cosmicdb as cosmicconv


def argument_parser(cmdargs):
    parser_dict = workflow.init_parser(
        description=(
            f'COSMIC database files "CosmicNCV.tsv.gz" and '
            f'"CosmicMutantExport.tsv.gz" are used to create a single vcf '
            f'file which can be used for annotation.'))

    workflow.add_outfile_arg(
        parser_dict['required'], required=True, must_not_exist='ask',
        help='Output vcf file path.')

    parser_dict['required'].add_argument(
        '--CosmicNCV', dest='ncv_path', 
        required=True, 
        help=(f'Previously downloaded unmodified CosmicNCV.tsv.gz file'))

    parser_dict['required'].add_argument(
        '--CosmicMutantExport', 
        dest='mutexp_path', 
        required=True, 
        help=(f'Previously downloaded unmodified CosmicMutantExport.tsv.gz file'))

    parser_dict['required'].add_argument(
        '--refver', dest="refver",
        required=True, choices=('GRCh37', 'GRCh38'),
        help=f'Reference genome version. Must be one of "GRCh37" or "GRCh38".')
        
    parser_dict['required'].add_argument(
        '--cosmic-version', dest="cosmicver",
        required=True,
        help=f'COSMIC database version (e.g. v96).')

    workflow.add_outfmt_arg(parser_dict['optional'], required=False)
    workflow.add_logging_args(parser_dict)

    args = parser_dict['main'].parse_args(cmdargs)
    return args


def main(cmdargs):
    args = argument_parser(cmdargs)
    #logger = toolsetup.setup_logger(args)
    cosmicconv.main(args.ncv_path, args.mutexp_path, args.outfile_path, args.refver, args.cosmicver)


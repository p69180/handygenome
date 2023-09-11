import handygenome.workflow as workflow
import handygenome.workflow.toolsetup as toolsetup
import handygenome.dbconv.dbsnp as dbsnpconv


def argument_parser(cmdargs):
    parser_dict = workflow.init_parser(
        description=(
            f"The original dbSNP VCF file is downloaded, then "
            f"modified to a handygenome-compatible VCF file."
        )
    )

    parser_dict["required"].add_argument(
        "--outfile",
        dest="outfile_path",
        required=True,
        help=f"Output vcf file path",
    )

    parser_dict["optional"].add_argument(
        "--downloaded-file",
        dest="download_path",
        required=False,
        help=f"Previously downloaded dbSNP file path",
    )

    workflow.add_refver_arg(
        parser_dict['required'], required=True, choices=['GRCh37', 'GRCh38'],
    )
    workflow.add_logging_args(parser_dict)

    args = parser_dict["main"].parse_args(cmdargs)

    return args


def main(cmdargs):
    args = argument_parser(cmdargs)
    logger = toolsetup.setup_logger(args)
    is_grch37 = (args.refver == 'GRCh37')

    dbsnpconv.main(
        download_path=args.download_path,
        outfile_path=args.outfile_path,
        is_grch37=is_grch37,
        logger=logger,
    )

    logger.info("ALL SUCCESSFULLY FINISHED")


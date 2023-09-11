import textwrap

import pysam

import handygenome.refgenome.refgenome as refgenome
import handygenome.workflow as workflow
import handygenome.workflow.toolsetup as toolsetup
import handygenome.variant.svcaller_parser as svcaller_parser
import handygenome.variant.varianthandler as varianthandler
import handygenome.vcfeditor.headerhandler as headerhandler
import handygenome.vcfeditor.indexing as indexing
import handygenome.variant.vcfspec as libvcfspec


def argument_parser(cmdargs):
    def sanity_check(args):
        pass

    parser_dict = workflow.init_parser(
        description=textwrap.dedent(f"""\
            - Split multiallelic lines.
            - Change SV caller records into breakends form.
            - Normalize non-SV records into the leftmost form.
            - Normalize SV records into the bnd1-advanced form."""
        )
    )

    workflow.add_infile_arg(parser_dict['required'], required=True)
    workflow.add_outfile_arg(parser_dict['required'], required=True, must_not_exist='ask')
    #workflow.add_fasta_arg(parser_dict['required'], required=True)
    workflow.add_refver_arg(parser_dict['required'], required=True, choices='all')
    workflow.add_outfmt_arg(parser_dict['optional'], required=False)

    workflow.add_logging_args(parser_dict)

    workflow.add_index_arg(parser_dict)

    args = parser_dict['main'].parse_args(cmdargs)
    sanity_check(args)

    return args


def process_nonsv(vr, out_vr_list, nonsv_set, logger, refver):
    """
    Normalize the input vr(the leftmost form).
    Check if the normalized vcfspec is duplicated, in which case the vr is
        discarded.
    """

    vcfspec_original = libvcfspec.Vcfspec.from_vr(vr, refver=refver)
    vcfspec_left = vcfspec_original.leftmost()
    ID = vcfspec_left.get_id()

    if ID in nonsv_set:
        logger.info(f'Input variant record is discarded since vcfspec is duplicated:\n{vr}')
    else:
        nonsv_set.add(ID)
        varianthandler.apply_vcfspec(vr, vcfspec_left)
        vr.id = ID
        out_vr_list.append(vr)


def process_sv(vr, out_vr_list, bnds_set, fasta, chromdict, logger):
    """
    From a SV vr, generate two variant records representing bnd1 and bnd2.
    The breakends are in pos1adv form.
    Check if the bndspec is duplicated, in which case the vr is discarded.
    """

    bnds = svcaller_parser.get_bnds_from_caller_vr(vr, fasta, chromdict) 

    ##
    ##

    if bnds is None:
        raise Exception(f'Error while parsing this SV variant record:\n{vr}')

    if bnds in bnds_set:
        logger.info(f'Input variant record is discarded since breakends '
                    f'spec is duplicated:\n{vr}')
        return
    else:
        bnds_set.add(bnds)
        vcfspec_bnd1 = bnds.get_vcfspec_bnd1()
        vcfspec_bnd2 = bnds.get_vcfspec_bnd2()
        ID_bnd1 = bnds.get_id_bnd1()
        ID_bnd2 = bnds.get_id_bnd2()

        varianthandler.apply_vcfspec(vr, vcfspec_bnd1)
        vr.id = ID_bnd1
        vr.info['MATEID'] = ID_bnd2

        vr_pos2 = vr.header.new_record()
        varianthandler.apply_vcfspec(vr_pos2, vcfspec_bnd2)
        vr_pos2.id = ID_bnd2
        vr_pos2.info['MATEID'] = ID_bnd1

        out_vr_list.append(vr)
        out_vr_list.append(vr_pos2)


def process_monoallelic_vr(vr, out_vr_list, nonsv_set, bnds_set, fasta, chromdict, logger, refver):
    # skip delly <INS> records
    if svcaller_parser.check_vr_format_delly(vr):
        if vr.alts[0] == '<INS>':
            return

    if varianthandler.check_SV(vr):
        process_sv(vr, out_vr_list, bnds_set, fasta, chromdict, logger)
    else:
        process_nonsv(vr, out_vr_list, nonsv_set, logger, refver)


def main(cmdargs):
    args = argument_parser(cmdargs)
    logger = toolsetup.setup_logger(args)

    logger.info('BEGINNING')

    fasta = refgenome.get_fasta(args.refver)
    chromdict = refgenome.get_chromdict(args.refver)
    refver = args.refver

    nonsv_set = set()
    bnds_set = set()
    out_vr_list = list()
    
    with pysam.VariantFile(args.infile_path) as in_vcf:
        #headerhandler.addmeta_MATEID(in_vcf.header)
        #headerhandler.addmeta_END(in_vcf.header)
        with pysam.VariantFile(args.outfile_path, mode=args.mode_pysam, header=in_vcf.header) as out_vcf:
            if vr.contig not in chromdict.contigs:
                continue

            for vr in in_vcf.fetch():
                if len(vr.alts) == 1:
                    process_monoallelic_vr(vr, out_vr_list, nonsv_set, bnds_set, fasta, chromdict, logger, refver)
                else:
                    for alt in vr.alts:
                        new_vr = vr.copy()
                        new_vr.alts = [alt]
                        process_monoallelic_vr(new_vr, out_vr_list, nonsv_set, bnds_set, fasta, chromdict, logger, refver)

            out_vr_list.sort(key=varianthandler.get_vr_sortkey(chromdict))

            for vr in out_vr_list:
                out_vcf.write(vr)

    if not args.donot_index:
        indexing.index_vcf(args.outfile_path)

    logger.info('ALL SUCCESSFULLY FINISHED')


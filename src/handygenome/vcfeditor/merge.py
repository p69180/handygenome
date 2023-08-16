import os
import itertools
import warnings
import textwrap

import pysam

import handygenome.refgenome as refgenome
import handygenome.workflow as workflow
import handygenome.variant.varianthandler as varianthandler
import handygenome.vcfeditor.headerhandler as headerhandler
import handygenome.vcfeditor.initvcf as initvcf
import handygenome.vcfeditor.misc as vcfmisc
import handygenome.variant.vcfspec as libvcfspec


LOGGER = workflow.get_logger()


def sanity_check_args(infile_path_list, outfile_path, isec, isec_indices,
                      outfile_must_not_exist):
    infile_path_list = workflow.arghandler_infilelist(infile_path_list)

    arghandler = workflow.get_arghandler_outfile(outfile_must_not_exist)
    outfile_path = arghandler(outfile_path)

    sanity_check_isec_indices(isec_indices, len(infile_path_list))

    return infile_path_list, outfile_path


def sanity_check_isec_indices(isec_indices, n_vcf):
    e_msg = textwrap.dedent(f"""\
        'isec_indices' must be a list composed of 0 or 1, \
        with the length the same as the number of vcf inputs.""")
    if isec_indices is not None:
        if len(isec_indices) != n_vcf:
            raise Exception(e_msg)
        elif not set(isec_indices).issubset({0, 1}):
            raise Exception(e_msg)


############################################


def isec_indices_to_bool(isec_indices, n_vcf):
    if isec_indices is None:
        isec_indices_bool = list(itertools.repeat(True, n_vcf))
    else:
        isec_indices_bool = [bool(x) for x in isec_indices]

    return isec_indices_bool


def load_vcf_data(vcf_list, remove_infoformat):
    """Returns:
        [
            {vcfspec, vcfspec, ... },
            {vcfspec, vcfspec, ... },
            ...,
        ] (same length as len(infile_path_list)) 
    """

    raw_vcfspecs = list()
    if remove_infoformat:
        vr_dict = None
    else:
        vr_dict = dict()
    for in_vcf in vcf_list:
        subset = set()
        for vr in in_vcf.fetch():
            vcfspec = libvcfspec.Vcfspec.from_vr(vr)
            #vcfspec = varianthandler.get_vcfspec(vr)
            subset.add(vcfspec)
            if not remove_infoformat:
                vr_dict.setdefault(vcfspec, list())
                vr_dict[vcfspec].append(vr)

        raw_vcfspecs.append(subset)

    return raw_vcfspecs, vr_dict


#def load_vcf_data_vcfpath(infile_path_list):
#    vcf_list = [pysam.VariantFile(x) for x in infile_path_list]
#    raw_vcfspecs, vr_dict = load_vcf_data(vcf_list)
#    for vcf in vcf_list:
#        vcf.close()
#
#    return raw_vcfspecs, vr_dict


########################################


def get_merged_header(infile_path_list):
    """Any conflicting INFO or FORMAT keys 
    (with regard to any of Type, Number, or Description) are discarded."""

    header_list = list()
    for vcf_path in infile_path_list:
        with pysam.VariantFile(vcf_path, 'r') as in_vcf:
            header_list.append(in_vcf.header)
    #compatible = headerhandler.check_header_compatibility(header_list)
    merged_header = headerhandler.merge_vcfheaders(header_list)

    return merged_header


def get_output_vcfspecs_isec(raw_vcfspecs, isec_indices_bool, chromdict):
    vcfspecs_included = set.intersection( 
        *itertools.compress(raw_vcfspecs, isec_indices_bool)
    )
    for subset in itertools.compress(
        raw_vcfspecs, [not x for x in isec_indices_bool]
    ):
        vcfspecs_included.difference_update(subset)
    output_vcfspecs = sorted(
        vcfspecs_included, 
        #key=common.get_vcfspec_sortkey(chromdict)
        key=(lambda x: x.get_sortkey(chromdict)),
    )
    return output_vcfspecs


def get_output_vcfspecs_union(raw_vcfspecs, chromdict):
    vcfspec_set = set.union(*raw_vcfspecs)
    output_vcfspecs = sorted(
        vcfspec_set, 
        #key=common.get_vcfspec_sortkey(chromdict),
        key=(lambda x: x.get_sortkey(chromdict)),
    )

    return output_vcfspecs


def merge_vr_dict(vr_dict, output_header, output_vcfspecs, logger):
    output_vr_list = list()
    NR = 0
    for vcfspec in output_vcfspecs:
        NR += 1
        merged_vr = varianthandler.merge(vr_dict[vcfspec], output_header)
        output_vr_list.append(merged_vr)
        del vr_dict[vcfspec]
        if NR % 50_000 == 0:
            logger.info(f'{NR:,} merged variant records created')

    return output_vr_list


def get_output_vrs_wo_infoformat(output_vcfspecs, output_header):
    output_vr_list = list()
    for vcfspec in output_vcfspecs:
        new_vr = output_header.new_record()
        varianthandler.apply_vcfspec(new_vr, vcfspec)
        output_vr_list.append(new_vr)

    return output_vr_list


def write(output_vr_list, output_vcfspecs, outfile_path, mode_pysam, 
          output_header):
    with pysam.VariantFile(outfile_path, mode=mode_pysam, 
                           header=output_header) as out_vcf:
        for vr in output_vr_list:
            out_vcf.write(vr)


#########################################


def main_file(
    infile_path_list, 
    outfile_path, 
    fasta_path, 
    remove_infoformat=False,
    isec=False, 
    isec_indices=None, 
    mode_bcftools='z', 
    mode_pysam=None,
    outfile_must_not_exist='ask',
    logger=None,
):
    if logger is None:
        logger = LOGGER

    logger.info('Beginning')

    # set basic parameters
    infile_path_list, outfile_path = sanity_check_args(
        infile_path_list, 
        outfile_path, 
        isec, 
        isec_indices, 
        outfile_must_not_exist,
    )
    isec_indices_bool = isec_indices_to_bool(isec_indices, len(infile_path_list))
    mode_pysam = vcfmisc.write_mode_arghandler(mode_bcftools, mode_pysam)
    fasta = pysam.FastaFile(fasta_path)
    chromdict = refgenome.ChromDict.from_fasta(fasta)
    if remove_infoformat:
        output_header = initvcf.create_header(chromdict=chromdict)
    else:
        output_header = get_merged_header(infile_path_list)

    # load input files data
    vcf_list = [pysam.VariantFile(x) for x in infile_path_list]

    logger.info('Loading input vcf data')
    raw_vcfspecs, vr_dict = load_vcf_data(vcf_list, remove_infoformat)

    # create output_vcfspecs
    logger.info('Getting union/intersection of vcfspecs')
    if isec:
        output_vcfspecs = get_output_vcfspecs_isec(raw_vcfspecs, isec_indices_bool, chromdict)
    else:
        output_vcfspecs = get_output_vcfspecs_union(raw_vcfspecs, chromdict)
    logger.info(f'{len(output_vcfspecs):,} variant records are to be written')

    # create output_vr_list
    if remove_infoformat:
        output_vr_list = get_output_vrs_wo_infoformat(output_vcfspecs, output_header)
    else:
        logger.info('Merging overlapping variant records. This step may take a while.')
        output_vr_list = merge_vr_dict(vr_dict, output_header, output_vcfspecs, logger)

    # write final result
    logger.info('Writing final result')
    write(output_vr_list, output_vcfspecs, outfile_path, mode_pysam, output_header)

    for vcf in vcf_list:
        vcf.close()
    fasta.close()


#def main_vcfplus(vcfplus_list, isec=False, isec_indices=None, merge=True):
#    #assert all(isinstance(x, vcfplus.VcfPlus) for x in vcfplus_list)
#    def get_vp_list(output_vcfspecs, vr_dict, fasta, chromdict):
#        vp_list = list()
#        for vcfspec in output_vcfspecs:
#            for vr in vr_dict[vcfspec]:
#                vp = variantplus.VariantPlus(vr, fasta, chromdict)
#                vp_list.append(vp)
#
#        return vp_list
#
#    sanity_check_isec_indices(isec_indices, len(vcfplus_list))
#    chromdict = vcfplus_list[0].chromdict
#    fasta = vcfplus_list[0].fasta
#    isec_indices_bool = isec_indices_to_bool(isec_indices, len(vcfplus_list))
#
#    merged_header = merge_pysamhdr(vcfp.vcf.header for vcfp in vcfplus_list)
#    raw_vcfspecs = load_vcf_data(vcf_list)
#    output_vcfspecs = vcfmergelib.get_output_vcfspecs(raw_vcfspecs, isec, isec_indices_bool, chromdict)
#    vr_dict = vcfmergelib.get_vr_dict(vcfplus_list, output_vcfspecs, merged_header, isec_indices_bool, merge)
#    vp_list = get_vp_list(output_vcfspecs, vr_dict, fasta, chromdict)
#
#    vcfp = vcfplus.VcfPlus(fasta = fasta)
#    vcfp.set_header(header = merged_header)
#    vcfp.init_vp(vp_list = vp_list)
#
#    return vcfp

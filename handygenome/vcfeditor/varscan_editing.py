import pysam

import handygenome.common as common
from handygenome.variant.vcfspec import Vcfspec
import handygenome.variant.vcfspec as libvcfspec
import handygenome.vcfeditor.initvcf as initvcf
import handygenome.vcfeditor.indexing as indexing
from handygenome.annotation.callerinfo import CallerInfo
from handygenome.variant.vcfspec import VcfspecComponents


PASS_SET = {'PASS'}


def sanity_check_pre_split(vr):
    # sanity check 1
    if len(vr.alts) != 1:
        raise Exception(f'Comma-separated ALT detected:\n{vr}')
    if '*' in vr.ref or '*' in vr.alts[0]:
        raise Exception(f'"*" included in REF or ALT:\n{vr}')


def sanity_check_post_split(vr, REFsp, ALTsp):
    if len(REFsp) > 1 and len(ALTsp) > 1:
        raise Exception(f'Both REF and ALT include "/":\n{vr}')

    if len(REFsp) > 2 or len(ALTsp) > 2:
        raise Exception(f'REF or ALT have more than two "/"-separated fields:\n{vr}')

    if len(ALTsp) > 1:
        # 1    937795  .   C   CATTT/ATTTT
        if len(vr.ref) > 1:
            raise Exception(f'ALT is multi-field but REF is not a single-character:\n{vr}')

        if all(len(x) == 1 for x in ALTsp):  # snps
            pass
        else:  # indels
            if not (
                ALTsp[0][0] == vr.ref and
                len(ALTsp[0]) > 1
            ):
                raise Exception(f'Unexpected multi-ALT pattern:\n{vr}')
                #raise Exception(f'ALT is multi-field but the first character of ALT is not the same as REF:\n{vr}')

    if len(REFsp) > 1:
        # 1    1029791 .   CT/+T   C
        if len(vr.alts[0]) > 1:
            raise Exception(f'REF is multi-field but ALT is not a single-character:\n{vr}')
        if not (
            REFsp[0][0] == vr.alts[0] and
            len(REFsp[0]) > 1
        ):
            raise Exception(f'Unexpected multi-REF pattern:\n{vr}')
        #if REFsp[0][0] != vr.alts[0]:
        #    raise Exception(f'REF is multi-field but the first character of REF is not the same as ALT:\n{vr}')


def handle_multi_ALT_indel(vr, fasta, REFsp, ALTsp):
    # Now REF is composed of a single letter
    # First element of ALTsp starts with a letter the same as REF. Strip it.
        # examples:
        # 1    937795  .   C   CATTT/ATTTT
        # 1    51864   .   C   CA/-A
    assert len(ALTsp[0]) > 1
    ALTsp[0] = ALTsp[0][1:]
    merge_components = list()

    for ALT_unit in ALTsp:
        if ALT_unit.startswith('-'):
            # It is postulated to be a deletion
            ALT_unit = ALT_unit[1:]  # strip leading "-"
            pos0 = vr.pos - 1
            fetched_seq_nextto_ref = fasta.fetch(vr.contig, pos0 + 1, pos0 + 1 + len(ALT_unit))
            if fetched_seq_nextto_ref != ALT_unit:
                raise Exception(
                    f'ALT string after "-" character is different from fetched reference sequence:\n'
                    f'Fetched reference sequence: {fetched_seq_nextto_ref} (chrom={vr.contig}, start0={pos0 + 1}, end0={pos0 + 1 + len(ALT_unit)})\n'
                    f'vr:\n{vr}'
                )
            ref = vr.ref + fetched_seq_nextto_ref
            alt = vr.ref
        else:
            # Postulated to be an insertion
            ref = vr.ref
            alt = vr.ref + ALT_unit

        merge_components.append(Vcfspec(chrom=vr.contig, pos=vr.pos, ref=ref, alts=(alt,), fasta=fasta))

    return libvcfspec.merge(*merge_components)


def handle_multi_REF(vr, fasta, REFsp, ALTsp):
    # Now REF is composed of a single letter
    # First element of REFsp starts with a letter the same as ALT. Strip it.
        # examples:
        # 1    1092725 .   CCCATCCCCGCCATCCCCG/+CCCG   C
        # 1    1029791 .   CT/+T   C
        # 1    1204224 .   ATT/T   A
    assert len(REFsp[0]) > 1
    REFsp[0] = REFsp[0][1:]
    merge_components = list()

    for REF_unit in REFsp:
        if REF_unit.startswith('+'):
            # It is postulated to be an insertion
            REF_unit = REF_unit[1:]
            ref = vr.alts[0]
            alt = vr.alts[0] + REF_unit
        else:
            # It is postulated to be a deletion
            ref = vr.alts[0] + REF_unit
            alt = vr.alts[0]

        merge_components.append(Vcfspec(chrom=vr.contig, pos=vr.pos, ref=ref, alts=(alt,), fasta=fasta))

    return libvcfspec.merge(*merge_components)


def extract_vcfspec(vr, fasta):
    sanity_check_pre_split(vr)

    REFsp = vr.ref.split('/')
    ALTsp = vr.alts[0].split('/')

    sanity_check_post_split(vr, REFsp, ALTsp)

    # main
    if len(REFsp) == 1 and len(ALTsp) == 1:
        return Vcfspec(chrom=vr.contig, pos=vr.pos, ref=vr.ref, alts=vr.alts, fasta=fasta)
    elif len(REFsp) == 1 and len(ALTsp) > 1:
        if all(len(x) == 1 for x in ALTsp):  # snps
            return Vcfspec(chrom=vr.contig, pos=vr.pos, ref=vr.ref, alts=ALTsp, fasta=fasta)
        else:
            return handle_multi_ALT_indel(vr, fasta, REFsp, ALTsp)
    elif len(REFsp) > 1 and len(ALTsp) == 1:
        return handle_multi_REF(vr, fasta, REFsp, ALTsp)


def write_outfile(vcfspec_list, outfile_path, output_header, mode, caller_info, index):
    with pysam.VariantFile(outfile_path, header=output_header, mode=mode) as out_vcf:
        for vcfspec in vcfspec_list:
            new_vr = output_header.new_record()
            vcfspec.apply_to_vr(new_vr)
            caller_info.write(new_vr)
            out_vcf.write(new_vr)

    if index:
        indexing.index_vcf(outfile_path)


def modify_varscan_vcf(
    infile_path_list, 
    somatic_outfile_path, 
    germline_outfile_path=None,
    refver=None, 
    fasta=None,
    mode_bcftools='z',
    mode_pysam=None, 
    pass_only=True,
    index=True,
):
    """Merge all files in "infile_path_list" to make a single sorted output.
    INFO/FORMAT contents are all removed.
    Args:
        germline_outfile_path: If not set, germline variant records are not written.
    """
    # sanity check
    if (
        (refver is None) and
        (fasta is None)
    ):
        raise Exception(f'At least one of "refver" or "fasta" must be set.')

    # set params
    if fasta is None:
        fasta = common.DEFAULT_FASTAS[refver]
    chromdict = common.ChromDict(fasta=fasta)

    output_header = initvcf.create_header(chromdict=chromdict)
    CallerInfo.add_meta(output_header)
    VcfspecComponents.add_meta(output_header)

    caller_info = CallerInfo.from_caller_names(['VarScan2'])

    mode = common.write_mode_arghandler(mode_bcftools, mode_pysam)
    somatic_vcfspec_list = list()
    germline_vcfspec_list = list()

    # collect vcfspecs
    for infile_path in infile_path_list:
        with pysam.VariantFile(infile_path) as in_vcf:
            for vr in in_vcf.fetch():
                if pass_only:
                    if set(vr.filter) != PASS_SET:
                        continue

                vcfspec = extract_vcfspec(vr, fasta)
                if vcfspec.chrom not in chromdict.contigs:
                    continue

                if vr.info['SOMATIC']:
                    somatic_vcfspec_list.append(vcfspec)
                else:
                    if germline_outfile_path is not None:
                        germline_vcfspec_list.append(vcfspec)

    # sort vcfspecs
    sortkey = common.get_vcfspec_sortkey(chromdict)
    somatic_vcfspec_list.sort(key=sortkey)
    germline_vcfspec_list.sort(key=sortkey)

    # write output
    write_outfile(somatic_vcfspec_list, somatic_outfile_path, output_header, mode, caller_info, index)
    if germline_outfile_path is not None:
        write_outfile(germline_vcfspec_list, germline_outfile_path, output_header, mode, caller_info, index)



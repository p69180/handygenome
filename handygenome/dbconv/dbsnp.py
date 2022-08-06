import os
import re
import argparse
import gzip
import tempfile
import gc

import pysam

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
workflow = importlib.import_module(".".join([top_package_name, "workflow"]))
assemblyspec = importlib.import_module(".".join([top_package_name, "assemblyspec"]))
initvcf = importlib.import_module(".".join([top_package_name, "vcfeditor", "initvcf"]))
libpopfreq = importlib.import_module(".".join([top_package_name, "annotation", "popfreq"]))
equivalents = importlib.import_module(".".join([top_package_name, "variantplus", "equivalents"]))
varianthandler = importlib.import_module(".".join([top_package_name, "variantplus", "varianthandler"]))


LOGGER = workflow.get_logger(name='dbSNP converter')


def get_chrname_converter(refver):
    if refver == 'GRCh37':
        asmblspec = assemblyspec.SPECS["grch37"]
        chrname_converter = asmblspec.dicts[("refseq", "nochr_plus_genbank")]
    elif refver == 'GRCh38':
        asmblspec = assemblyspec.SPECS["grch38"]
        chrname_converter = asmblspec.dicts[("refseq", "ucsc")]

    return chrname_converter


def get_tmpfile_paths(tmpdir, original_vcf_path, chromdict, chrname_converter):
    with pysam.VariantFile(original_vcf_path) as in_vcf:
        chromlist = list()
        for chrom in in_vcf.header.contigs.keys():
            new_chrom = chrname_converter[chrom]
            if (
                (new_chrom is not None)
                and (new_chrom in chromdict.contigs)
            ):
                chromlist.append((chrom, new_chrom))

        chromlist.sort(key=(lambda x: chromdict.contigs.index(x[1])))

    tmpfile_paths = {new_chrom: os.path.join(tmpdir, f'{new_chrom}.vcf.gz')
                     for new_chrom in [x[1] for x in chromlist]}

    return chromlist, tmpfile_paths


#####################


def get_tmp_new_header(chromdict):
    header = initvcf.create_header(chromdict=chromdict)
    libpopfreq.PopfreqInfoList.add_meta_info(header)
    
    return header


def handle_freq(freq):
    freq_split = (",".join(freq)).split("|")
    freq_dict = dict()
    for item in freq_split:
        item_split = item.split(":")
        freq_dict[item_split[0]] = [
            (float(0) if x == "." else float(x)) for x in item_split[1].split(",")
        ]

    return freq_dict


def make_popfreqinfo(original_vr, freq_dict, alt_idx):
    popinfo = libpopfreq.PopfreqInfo()

    popinfo['id'] = f'rs{original_vr.info["RS"]}'

    if "dbSNPBuildID" in original_vr.info:
        popinfo['dbSNPBuild'] = original_vr.info["dbSNPBuildID"]
    else:
        popinfo['dbSNPBuild'] = None

    assert isinstance(original_vr.info["COMMON"], bool)
    popinfo['common'] = original_vr.info["COMMON"]

    result_freqs = dict()
    for pop in freq_dict.keys():
        allele_freq = freq_dict[pop][alt_idx + 1]
        result_freqs[pop] = allele_freq
    popinfo['freqs'] = result_freqs

    popinfolist = libpopfreq.PopfreqInfoList()
    popinfolist.append(popinfo)

    return popinfolist


def get_new_vr_list(original_vr, freq_dict, new_chrom, fasta, tmp_new_header):
    vcfspec = varianthandler.get_vcfspec(original_vr)
    vcfspec.chrom = new_chrom

    new_vr_list = list()
    for alt_idx, sub_vcfspec in enumerate(vcfspec.iter_monoalts()):
        sub_vcfspec = equivalents.leftmost(sub_vcfspec, fasta)
        new_vr = tmp_new_header.new_record()
        varianthandler.apply_vcfspec(new_vr, sub_vcfspec)
        popinfolist = make_popfreqinfo(original_vr, freq_dict, alt_idx)
        popinfolist.write_info(new_vr)
        new_vr_list.append(new_vr)

    return new_vr_list


def process_by_chrom(original_vcf, chrom, new_chrom, tmp_outfile_path, popnames, tmp_new_header, fasta, chromdict, logging_lineno):
    # 1. setup vr dictionary to write to tmp output file
    out_vrs = dict()
    
    LOGGER.info(f'Collecting popnames and modifying original variant records')
    for original_vr in common.iter_lineno_logging(original_vcf.fetch(contig=chrom), LOGGER, logging_lineno):
        if "FREQ" not in original_vr.info:
            continue
        freq_dict = handle_freq(original_vr.info["FREQ"])
        popnames.update(freq_dict.keys())
        new_vr_list = get_new_vr_list(original_vr, freq_dict, new_chrom, fasta, tmp_new_header)
        for new_vr in new_vr_list:
            vcfspec = varianthandler.get_vcfspec(new_vr)
            out_vrs[vcfspec] = new_vr

    # 2. write to tmp output file
    LOGGER.info(f'BEGINNING sorting vcfspecs')
    vcfspecs_sorted = sorted(out_vrs.keys(), key=common.get_vcfspec_sortkey(chromdict))
    LOGGER.info(f'FINISHED sorting vcfspecs')

    out_vcf = pysam.VariantFile(tmp_outfile_path, mode='wz', header=tmp_new_header)
    LOGGER.info(f'Writing temporary outfile')
    for vcfspec in common.iter_lineno_logging(vcfspecs_sorted, LOGGER, logging_lineno):
        out_vcf.write(out_vrs[vcfspec])
    out_vcf.close()

    del out_vrs
    gc.collect()


def preprocess(original_vcf_path, chromlist, tmpfile_paths, fasta, chromdict, logging_lineno):
    popnames = set()
    original_vcf = pysam.VariantFile(original_vcf_path)
    tmp_new_header = get_tmp_new_header(chromdict)

    for chrom, new_chrom in chromlist:
        LOGGER.info(f'BEGINNING preprocessing for chromosome {new_chrom}({chrom})')
        tmp_outfile_path = tmpfile_paths[new_chrom]
        process_by_chrom(original_vcf, chrom, new_chrom, tmp_outfile_path, popnames, tmp_new_header, fasta, chromdict, logging_lineno)

    original_vcf.close()
    return popnames


######################


def get_popfreq_metadata(popnames):
    popfreq_metadata = libpopfreq.PopfreqMetadata()
    popfreq_metadata['popnames'] = popnames

    return popfreq_metadata


def get_final_new_header(chromdict, popfreq_metadata):
    header = initvcf.create_header(chromdict=chromdict)
    libpopfreq.PopfreqInfoList.add_meta_info(header)
    popfreq_metadata.write_header(header)
    
    return header


def write_final_outfile(popnames, chromdict, tmpfile_paths, outfile_path, chromlist, logging_lineno):
    popfreq_metadata = get_popfreq_metadata(popnames)
    final_new_header = get_final_new_header(chromdict, popfreq_metadata)
    out_vcf = pysam.VariantFile(outfile_path, mode='wz', header=final_new_header)

    for chrom, new_chrom in chromlist:
        LOGGER.info(f'Writing final output for chromosome {new_chrom}')
        tmp_outfile_path = tmpfile_paths[new_chrom]
        with pysam.VariantFile(tmp_outfile_path, mode='r') as in_vcf:
            for vr in common.iter_lineno_logging(in_vcf.fetch(), LOGGER, logging_lineno):
                out_vcf.write(vr)

    out_vcf.close()


######################


@common.get_deco_arg_choices({'refver': ('GRCh37', 'GRCh38')})
def main(original_vcf_path, outfile_path, refver):
    chromdict = common.ChromDict(refver=refver)
    chrname_converter = get_chrname_converter(refver)
    tmpdir = tempfile.mkdtemp(prefix='dbsnp_converter', dir=os.path.dirname(outfile_path))
    chromlist, tmpfile_paths = get_tmpfile_paths(tmpdir, original_vcf_path, chromdict, chrname_converter)

    popnames = preprocess(original_vcf_path, chromlist, tmpfile_paths, fasta, chromdict, logging_lineno=1_000_000)

    write_final_outfile(popnames, chromdict, tmpfile_paths, outfile_path, chromlist, logging_lineno=1_000_000)

    LOGGER.info(f'ALL SUCCESSFULLY FINISHED')


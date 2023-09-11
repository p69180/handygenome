import os
import re
import argparse
import gzip
import ftplib
import contextlib

import pysam

import handygenome.network as network
import handygenome.workflow as workflow
import handygenome.workflow.toolsetup as toolsetup
import handygenome.refgenome.refgenome as refgenome
import handygenome.vcfeditor.initvcf as initvcf
import handygenome.vcfeditor.indexing as indexing
import handygenome.publicdb.ncbi as libncbi


DEFAULT_OUTFILE_BASENAME_GRCH37 = "dbSNP_b155_GRCh37.p13.vcf.gz"
DEFAULT_OUTFILE_BASENAME_GRCH38 = "dbSNP_b155_GRCh38.p13.vcf.gz"
DEFAULT_OUTFILE_PATH_GRCH37 = os.path.join(os.getcwd(), DEFAULT_OUTFILE_BASENAME_GRCH37)
DEFAULT_OUTFILE_PATH_GRCH38 = os.path.join(os.getcwd(), DEFAULT_OUTFILE_BASENAME_GRCH38)

#URL_GRCH37 = "https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.25.gz"
#URL_GRCH37_MD5 = "https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.25.gz.md5"
#URL_GRCH38 = "https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz"
#URL_GRCH38_MD5 = "https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz.md5"

URL_GRCH37 = "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz"
URL_GRCH37_MD5 = "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.md5"
URL_GRCH38 = "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"
URL_GRCH38_MD5 = "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.md5"


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
        #default=DEFAULT_OUTFILE_PATH_GRCH37,
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


# Downloading original VCF file

def get_refver(url):
    response = network.http_run_urlopen(url)
    gf = gzip.GzipFile(fileobj=response, mode='rb')
    for bline in gf:
        if bline.startswith(b'##reference'):
            refver = bline.decode().split('=')[1]  # GRCh37.p13
            refver = refver.split('.')[0]
            return refver
        if not bline.startswith(b'##'):
            break

    raise Exception(f'No reference version information found from VCF header')


def get_refver_pysam(url):
    with pysam.VariantFile(url) as vcf:
        for record in vcf.header.records:
            if record.key == 'reference':
                refver = record.value  # GRCh37.p13
                refver = refver.split('.')[0]
                return refver

    raise Exception(f'No reference version information found from VCF header')


def pick_vcf_url(url_dict, is_grch37=True):
    if is_grch37:
        key = min(url_dict.keys(), key=(lambda x: int(x.split('.')[1])))
    else:
        key = max(url_dict.keys(), key=(lambda x: int(x.split('.')[1])))

    url = url_dict[key]['vcf']
    refver = get_refver(url)
    #refver = get_refver_pysam(url)
        # not used due to connection failure risk
    if (
        (is_grch37 and (refver != 'GRCh37'))
        or
        ((not is_grch37) and (refver != 'GRCh38'))
    ):
        raise Exception(f'Unexpected filename-refver pairing; url_dict={url_dict}, url={url}, refver={refver}')

    return url


def download(logger, url, download_path):
    logger.info(f'Beginning download of url "{url}" to "{download_path}"')
    network.download(url, download_path)
    logger.info(f"Download finished")


##########################################

# Modification of original VCF file

def get_generic_metadata(vcf):
    result = dict()
    for rec in vcf.header.records:
        if rec.type == 'GENERIC':
            result[rec.key] = rec.value

    return result


def get_dbsnp_build(vcf):
    for record in vcf.header.records:
        if record.key == 'dbSNP_BUILD_ID':
            build = record.value  # class 'str'
            return build

    raise Exception(f'No build information found from VCF header')


def get_freq_classes(vcf, logger):
    result = set()
    NR = 0
    for vr in vcf.fetch():
        NR += 1
        if NR % 10_000_000 == 0:
            logger.info(
                f"get_freq_classes: Processing {NR:,}th variant " f"record:\n{vr}"
            )

        try:
            freq = vr.info["FREQ"]
        except KeyError:
            continue

        freq = ",".join(freq)
        classes = [x.split(":")[0] for x in freq.split("|")]
        result.update(classes)

    result = sorted(result)

    return result


def get_new_header(chromdict, freq_classes, original_vcf):
    new_header = initvcf.create_header(chromdict=chromdict)

    for rec in original_vcf.header.records:
        if rec.key in ('fileDate', 'dbSNP_BUILD_ID', 'reference'):
            new_header.add_record(rec)

    new_header.add_meta(
        key="INFO",
        items=[
            ("ID", "rs"),
            ("Type", "Integer"),
            ("Number", 1),
            ("Description", "dbSNP rs ID"),
        ],
    )
    new_header.add_meta(
        key="INFO",
        items=[
            ("ID", "dbSNPBuildID"),
            ("Type", "Integer"),
            ("Number", 1),
            ("Description", "First dbSNP Build for RS"),
        ],
    )
    new_header.add_meta(
        key="INFO",
        items=[
            ("ID", "COMMON"),
            ("Type", "Flag"),
            ("Number", 0),
            (
                "Description",
                (
                    "RS is a common SNP. A common SNP is one that has at least "
                    "one 1000Genomes population with a minor allele of "
                    "frequency >= 1% and for which 2 or more founders "
                    "contribute to that minor allele frequency."
                ),
            ),
        ],
    )

    for pop in freq_classes:
        new_header.add_meta(
            key="INFO",
            items=[
                ("ID", f"AF_{pop}"),
                ("Type", "Float"),
                ("Number", 1),
                ("Description", f"Allele frequency in the population '{pop}'"),
            ],
        )

    return new_header


def handle_freq(freq):
    freq_split = ",".join(freq).split("|")
    freq_dict = dict()
    for item in freq_split:
        item_split = item.split(":")
        freq_dict[item_split[0]] = [
            (float(0) if x == "." else float(x)) for x in item_split[1].split(",")
        ]

    return freq_dict


def get_new_vr_list(
    vr, chromdict, chrname_converter, new_header, freq_classes, infile_vcfspecs
):
    new_vr_list = list()
    infile_vcfspecs.setdefault(vr.contig, set())

    new_chrom = chrname_converter[vr.contig]
    if (
        (new_chrom is not None)
        and (new_chrom in chromdict.contigs)
        and "FREQ" in vr.info
    ):
        freq = vr.info["FREQ"]
        freq_dict = handle_freq(freq)
        for idx, alt in enumerate(vr.alts):
            vcfspec = (vr.pos, vr.ref, alt)
            if vcfspec in infile_vcfspecs[vr.contig]:
                continue

            infile_vcfspecs[vr.contig].add(vcfspec)

            new_vr = new_header.new_record()
            new_vr.chrom = new_chrom
            new_vr.pos = vr.pos
            new_vr.ref = vr.ref
            new_vr.alts = (alt,)

            new_vr.info["rs"] = vr.info["RS"]
            if "dbSNPBuildID" in vr.info:
                new_vr.info["dbSNPBuildID"] = vr.info["dbSNPBuildID"]
            if vr.info["COMMON"]:
                new_vr.info["COMMON"] = True

            for pop in freq_classes:
                if pop in freq_dict:
                    AF = freq_dict[pop][idx + 1]
                    new_vr.info[f"AF_{pop}"] = AF

            new_vr_list.append(new_vr)

    return new_vr_list


def write_outfile(
    original_vcf,
    outfile_path,
    chromdict,
    chrname_converter,
    new_header,
    freq_classes,
    logger,
):
    infile_vcfspecs = dict()
    with pysam.VariantFile(outfile_path, mode="wz", header=new_header) as out_vcf:
        for idx, vr in enumerate(original_vcf.fetch()):
            if idx % 1_000_000 == 0:
                logger.info(
                    f"write_outfile: Processing {idx + 1:,}th variant record:\n{vr}"
                )

            new_vr_list = get_new_vr_list(
                vr,
                chromdict,
                chrname_converter,
                new_header,
                freq_classes,
                infile_vcfspecs,
            )
            for new_vr in new_vr_list:
                out_vcf.write(new_vr)


#def write_outfile_old(
#    download_path,
#    outfile_path,
#    chromdict,
#    chrname_converter,
#    new_header,
#    freq_classes,
#    logger,
#):
#    infile_vcfspecs = dict()
#    with pysam.VariantFile(
#        outfile_path, mode="wz", header=new_header
#    ) as out_vcf, pysam.VariantFile(download_path, mode="r") as in_vcf:
#        NR = 0
#        for vr in in_vcf.fetch():
#            NR += 1
#            if NR % 1_000_000 == 0:
#                logger.info(
#                    f"write_outfile: Processing {NR:,}th variant " f"record:\n{vr}"
#                )
#
#            new_vr_list = get_new_vr_list(
#                vr,
#                chromdict,
#                chrname_converter,
#                new_header,
#                freq_classes,
#                infile_vcfspecs,
#            )
#            for new_vr in new_vr_list:
#                out_vcf.write(new_vr)


#def main_old(
#    download_path,
#    download_url,
#    outfile_path,
#    chromdict,
#    chrname_converter,
#    refver,
#    build,
#    logger,
#):
#    # download
#    if download_path is None:
#        download_path = workflow.get_tmpfile_path()
#        download(logger, download_url, download_path)
#        rm_downloaded = True
#    else:
#        rm_downloaded = False
#
#    # collect FREQ population names
#    logger.info("Collecting FREQ population names from the original file")
#    freq_classes = get_freq_classes(download_path, logger)
#
#    # write output file
#    new_header = get_new_header(chromdict, freq_classes, refver=refver, build=build)
#    logger.info("Modifying and writing original variant records")
#    write_outfile(
#        download_path,
#        outfile_path,
#        chromdict,
#        chrname_converter,
#        new_header,
#        freq_classes,
#        logger,
#    )
#
#    # indexing
#    indexing.index_vcf(outfile_path, overwrite=True)
#
#    # remove temp files
#    if rm_downloaded:
#        os.remove(download_path)


def main(
    download_path,
    outfile_path,
    is_grch37,
    logger,
):
    # download
    if download_path is None:
        url_dict = libncbi.get_dbsnp_urls()
        ftp_url = pick_vcf_url(url_dict, is_grch37=is_grch37)

        download_path = workflow.get_tmpfile_path(prefix="dbSNP_vcf_download_tmpfile_", suffix=".vcf.gz", dir=os.getcwd(), is_dir=False, delete=False)
        logger.info(f"Downloading original VCF file from NCBI ftp server to {repr(download_path)}")
        network.download_wget(ftp_url, download_path)
        logger.info("Download finished")

        url = download_path
        rm_tmp = True
    else:
        url = download_path
        rm_tmp = False

    original_vcf = pysam.VariantFile(url)

    # collect FREQ population names
    logger.info("Collecting FREQ population names from the original file")
    freq_classes = get_freq_classes(original_vcf, logger)
    logger.info("Finished collecting FREQ population names")

    # write output file
    if is_grch37:
        assembly_spec = refgenome.SPECS["GRCh37"]
        chromdict = assembly_spec.chromdicts["nochr_plus_genbank"]
        chrname_converter = assembly_spec.dicts[("refseq", "nochr_plus_genbank")]
    else:
        assembly_spec = refgenome.SPECS["GRCh38"]
        chromdict = assembly_spec.chromdicts["ucsc"]
        chrname_converter = assembly_spec.dicts[("refseq", "ucsc")]

    new_header = get_new_header(chromdict, freq_classes, original_vcf)
    logger.info("Modifying original variant records and writing to the output file")
    write_outfile(
        original_vcf,
        outfile_path,
        chromdict,
        chrname_converter,
        new_header,
        freq_classes,
        logger,
    )
    logger.info("Finished writing the output file")

    # indexing
    indexing.index_vcf(outfile_path, overwrite=True)

    # remove downloaded file
    if rm_tmp:
        os.remove(download_path)


#def main_grch37(args, logger):
#    asmblspec = refgenome.SPECS["grch37"]
#    main_old(
#        download_path=args.download_path_grch37,
#        download_url=URL_GRCH37,
#        outfile_path=args.outfile_path_grch37,
#        chromdict=asmblspec.chromdicts["nochr_plus_genbank"],
#        chrname_converter=asmblspec.dicts[("refseq", "nochr_plus_genbank")],
#        refver="GRCh37.p13",
#        build="155",
#        logger=logger,
#    )
#
#
#def main_grch38(args, logger):
#    asmblspec = refgenome.SPECS["grch38"]
#    main_old(
#        download_path=args.download_path_grch38,
#        download_url=URL_GRCH38,
#        outfile_path=args.outfile_path_grch38,
#        chromdict=asmblspec.chromdicts["ucsc"],
#        chrname_converter=asmblspec.dicts[("refseq", "ucsc")],
#        refver="GRCh38.p13",
#        build="155",
#        logger=logger,
#    )


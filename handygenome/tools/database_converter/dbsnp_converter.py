import os
import re
import argparse
import gzip

import pysam

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
workflow = importlib.import_module(".".join([top_package_name, "workflow"]))
toolsetup = importlib.import_module(
    ".".join([top_package_name, "workflow", "toolsetup"])
)
assemblyspec = importlib.import_module(".".join([top_package_name, "assemblyspec"]))
initvcf = importlib.import_module(".".join([top_package_name, "vcfeditor", "initvcf"]))
indexing = importlib.import_module(
    ".".join([top_package_name, "vcfeditor", "indexing"])
)


DEFAULT_OUTFILE_BASENAME_GRCH37 = "dbSNP_b155_GRCh37.p13.vcf.gz"
DEFAULT_OUTFILE_BASENAME_GRCH38 = "dbSNP_b155_GRCh38.p13.vcf.gz"
DEFAULT_OUTFILE_PATH_GRCH37 = os.path.join(os.getcwd(), DEFAULT_OUTFILE_BASENAME_GRCH37)
DEFAULT_OUTFILE_PATH_GRCH38 = os.path.join(os.getcwd(), DEFAULT_OUTFILE_BASENAME_GRCH38)

URL_GRCH37 = "https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.25.gz"
URL_GRCH37_MD5 = "https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.25.gz.md5"
URL_GRCH38 = "https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz"
URL_GRCH38_MD5 = "https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz.md5"


def argument_parser(cmdargs):
    parser_dict = workflow.init_parser(
        description=(
            f"The original dbSNP vcf file is downloaded, then "
            f"modified to a vcf file with a modified format."
        )
    )

    parser_dict["optional"].add_argument(
        "--outfile-grch37",
        dest="outfile_path_grch37",
        default=DEFAULT_OUTFILE_PATH_GRCH37,
        required=False,
        help=f"Output vcf file path for grch37",
    )

    parser_dict["optional"].add_argument(
        "--outfile-grch38",
        dest="outfile_path_grch38",
        default=DEFAULT_OUTFILE_PATH_GRCH38,
        required=False,
        help=f"Output vcf file path for grch38",
    )

    parser_dict["optional"].add_argument(
        "--downloaded-file-grch37",
        dest="download_path_grch37",
        required=False,
        help=(
            f"(Optional) previously downloaded original gene set gff3 "
            f"file path for grch37"
        ),
    )

    parser_dict["optional"].add_argument(
        "--downloaded-file-grch38",
        dest="download_path_grch38",
        required=False,
        help=f"(Optional) previously downloaded original gene set gff3 file path for grch38",
    )

    workflow.add_logging_args(parser_dict)

    args = parser_dict["main"].parse_args(cmdargs)

    return args


def download(logger, url, download_path):
    logger.info(f'Beginning download of url "{url}" to "{download_path}"')
    common.download(url, download_path)
    logger.info(f"Download finished")


def get_freq_classes(download_path, logger):
    result = set()
    with pysam.VariantFile(download_path, "r") as in_vcf:
        NR = 0
        for vr in in_vcf.fetch():
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


def get_new_header(chromdict, freq_classes, refver, build="155"):
    new_header = initvcf.create_header(chromdict=chromdict)

    new_header.add_meta(key="dbSNP_BUILD_ID", value=build)
    new_header.add_meta(key="reference", value=refver)
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
    download_path,
    outfile_path,
    chromdict,
    chrname_converter,
    new_header,
    freq_classes,
    logger,
):
    infile_vcfspecs = dict()
    with pysam.VariantFile(
        outfile_path, mode="wz", header=new_header
    ) as out_vcf, pysam.VariantFile(download_path, mode="r") as in_vcf:
        NR = 0
        for vr in in_vcf.fetch():
            NR += 1
            if NR % 1_000_000 == 0:
                logger.info(
                    f"write_outfile: Processing {NR:,}th variant " f"record:\n{vr}"
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


def submain(
    download_path,
    download_url,
    outfile_path,
    chromdict,
    chrname_converter,
    refver,
    build,
    logger,
):
    # download
    if download_path is None:
        download_path = workflow.get_tmpfile_path()
        download(logger, download_url, download_path)
        rm_downloaded = True
    else:
        rm_downloaded = False

    # collect FREQ population names
    logger.info("Collecting FREQ population names from the original file")
    freq_classes = get_freq_classes(download_path, logger)

    # write output file
    new_header = get_new_header(chromdict, freq_classes, refver=refver, build=build)
    logger.info("Modifying and writing original variant records")
    write_outfile(
        download_path,
        outfile_path,
        chromdict,
        chrname_converter,
        new_header,
        freq_classes,
        logger,
    )

    # indexing
    indexing.index_vcf(outfile_path, overwrite=True)

    # remove temp files
    if rm_downloaded:
        os.remove(download_path)


def submain_grch37(args, logger):
    asmblspec = assemblyspec.SPECS["grch37"]
    submain(
        download_path=args.download_path_grch37,
        download_url=URL_GRCH37,
        outfile_path=args.outfile_path_grch37,
        chromdict=asmblspec.chromdicts["nochr_plus_genbank"],
        chrname_converter=asmblspec.dicts[("refseq", "nochr_plus_genbank")],
        refver="GRCh37.p13",
        build="155",
        logger=logger,
    )


def submain_grch38(args, logger):
    asmblspec = assemblyspec.SPECS["grch38"]
    submain(
        download_path=args.download_path_grch38,
        download_url=URL_GRCH38,
        outfile_path=args.outfile_path_grch38,
        chromdict=asmblspec.chromdicts["ucsc"],
        chrname_converter=asmblspec.dicts[("refseq", "ucsc")],
        refver="GRCh38.p13",
        build="155",
        logger=logger,
    )


def main(cmdargs):
    args = argument_parser(cmdargs)
    logger = toolsetup.setup_logger(args)

    # submain_grch37(args, logger)
    submain_grch38(args, logger)

    logger.info("ALL SUCCESSFULLY FINISHED")

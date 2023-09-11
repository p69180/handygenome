import os
import re
import argparse
import gzip

import pysam

import handygenome.tools as tools
import handygenome.network as network
import handygenome.workflow as workflow
import handygenome.workflow.toolsetup as toolsetup
import handygenome.refgenome.refgenome as refgenome


DEFAULT_OUTFILE_BASENAME = 'ucsc_repeatmasker_out_sorted.bed.gz'
URL_GRCH37 = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.out.gz'
URL_GRCH38 = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz'
DEFAULT_OUTFILE_PATH_GRCH37 = os.path.join(os.getcwd(), DEFAULT_OUTFILE_BASENAME)
DEFAULT_OUTFILE_PATH_GRCH38 = os.path.join(os.getcwd(), DEFAULT_OUTFILE_BASENAME)


def argument_parser(cmdargs):
    parser_dict = workflow.init_parser(
        description=(f'The original RepeatMasker output files are '
                     f'downloaded from ucsc, then filtered, sorted, '
                     f'bgzipped, and tabix-indexed.'))

    parser_dict['optional'].add_argument(
        '--outfile-grch37', dest='outfile_path_grch37', 
        default=DEFAULT_OUTFILE_PATH_GRCH37,
        required=False, help=f'Output repeat bed file path for grch37')
    parser_dict['optional'].add_argument(
        '--outfile-grch38', dest='outfile_path_grch38', 
        default=DEFAULT_OUTFILE_PATH_GRCH38,
        required=False, help=f'Output repeat bed file path for grch38')
    parser_dict['optional'].add_argument(
        '--downloaded-file-grch37', dest='download_path_grch37', 
        required=False, 
        help=(f'(Optional) previously downloaded original ucsc file path '
              f'for grch37'))
    parser_dict['optional'].add_argument(
        '--downloaded-file-grch38', dest='download_path_grch38', 
        required=False, 
        help=(f'(Optional) previously downloaded original ucsc file path '
              f'for grch38'))

    workflow.add_logging_args(parser_dict)

    args = parser_dict['main'].parse_args(cmdargs)

    return args


def download(logger, url, download_path):
    logger.info(f'Beginning download of url "{url}" to "{download_path}"')
    network.download(url, download_path)
    logger.info(f'Download finished')


def parse(download_path, asmblspec, output_chrname_version):
    parse_result = list()
    with gzip.open(download_path, 'rt') as infile:
        for _ in range(3):
            line = next(infile)

        for line in infile:
            linestrip = tools.rm_newline(line)
            linesp = linestrip.split()

            chrom = asmblspec.convert(linesp[4], output_chrname_version)
            start = int(linesp[5]) - 1
            end = int(linesp[6])
            name = linesp[10] + ':' + linesp[9]
            score = int(linesp[0])
            if linesp[8] == '+':
                strand = '+'
            elif linesp[8] == 'C':
                strand = '-'
            else:
                raise Exception(f'Unexpected strand string of repeatmasker '
                                f'output: {linesp[8]}')

            line_parse = [chrom, start, end, name, score, strand]
            parse_result.append(line_parse)

    return parse_result


def sort_parse_result(chromdict, parse_result):
    def sortkey(line_parse):
        return tools.coord_sortkey(line_parse[0], line_parse[1], chromdict)
    parse_result_sorted = sorted(parse_result, key=sortkey)

    return parse_result_sorted


def write_uncompressed_file(outfile_path_uncomp, parse_result_sorted):
    with open(outfile_path_uncomp, 'wt') as outfile:
        outfile.write(
            '\t'.join(['#contig', 'start', 'end', 'class/family:name', 
                       'SW_score', 'strand']) + '\n')
        for line_parse in parse_result_sorted:
            outfile.write('\t'.join(str(x) for x in line_parse) + '\n')


def convert_common(original_url, output_chrname_version, refver, logger, 
                   outfile_path, download_path=None):
    # download
    if download_path is None:
        download_path = workflow.get_tmpfile_path()
        download(logger, original_url, download_path)
        rm_downloaded = True
    else:
        rm_downloaded = False

    # set refgenome dbs
    asmblspec = refgenome.SPECS[refver]
    chromdict = asmblspec.chromdicts[output_chrname_version]

    # load
    logger.info(f'Loading the original file')
    parse_result = parse(download_path, asmblspec, output_chrname_version)

    # sort
    logger.info(f'Sorting the original file lines')
    parse_result_sorted = sort_parse_result(chromdict, parse_result)

    # write uncompressed outfile
    logger.info(f'Writing the sorted lines')
    outfile_path_uncomp = outfile_path + '.uncompressed'
    write_uncompressed_file(outfile_path_uncomp, parse_result_sorted)

    # compress and index
    logger.info(f'Compressing with bgzip and indexing with tabix')
    pysam.tabix_compress(outfile_path_uncomp, outfile_path)
    pysam.tabix_index(outfile_path, preset='bed')

    # remove temp files
    os.remove(outfile_path_uncomp)
    if rm_downloaded:
        os.remove(download_path)


def convert_grch37(logger, outfile_path, download_path=None):
    convert_common(URL_GRCH37, 'nochr_plus_genbank', 'hg19', logger, 
                   outfile_path, download_path)


def convert_grch38(logger, outfile_path, download_path=None):
    convert_common(URL_GRCH38, 'ucsc', 'hg38', logger, outfile_path, 
                   download_path)


def main(cmdargs):
    args = argument_parser(cmdargs)
    logger = toolsetup.setup_logger(args)

    convert_grch37(logger, args.outfile_path_grch37, args.download_path_grch37)
    convert_grch38(logger, args.outfile_path_grch38, args.download_path_grch38)
    
    logger.info('ALL SUCCESSFULLY FINISHED')

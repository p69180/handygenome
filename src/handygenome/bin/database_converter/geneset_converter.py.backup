import os
import re
import argparse
import gzip

import pysam

import handygenome.network as network
import handygenome.tools as tools
import handygenome.workflow as workflow
import handygenome.workflow.toolsetup as toolsetup
import handygenome.refgenome as refgenome


DEFAULT_OUTFILE_BASENAME = 'ensembl_geneset_sorted.gff3.gz'
URL_GRCH37 = 'http://ftp.ensembl.org/pub/grch37/release-105/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz'
URL_GRCH38 = 'http://ftp.ensembl.org/pub/current_gff3/homo_sapiens/Homo_sapiens.GRCh38.105.gff3.gz'
DEFAULT_OUTFILE_PATH_GRCH37 = os.path.join(os.getcwd(), re.sub('gff3.gz$', 'grch37.gff3.gz', DEFAULT_OUTFILE_BASENAME))
DEFAULT_OUTFILE_PATH_GRCH38 = os.path.join(os.getcwd(), re.sub('gff3.gz$', 'grch38.gff3.gz', DEFAULT_OUTFILE_BASENAME))


def argument_parser(cmdargs):
    parser_dict = workflow.init_parser(
        description=(
            f'The original gff3 gene set files are downloaded, filtered, '
            f'sorted, bgzipped, and tabix-indexed. These lines are filtered '
            f'out: containing "###" or "feature" value is "chromosome" '
            f'or "supercontig".'))

    parser_dict['optional'].add_argument(
        '--outfile-grch37', dest='outfile_path_grch37', 
        default=DEFAULT_OUTFILE_PATH_GRCH37,
        required=False, help=f'Output gff3 file path for grch37')

    parser_dict['optional'].add_argument(
        '--outfile-grch38', dest='outfile_path_grch38', 
        default=DEFAULT_OUTFILE_PATH_GRCH38,
        required=False, help=f'Output gff3 file path for grch38')

    parser_dict['optional'].add_argument(
        '--downloaded-file-grch37', dest='download_path_grch37', 
        required=False, 
        help=(f'(Optional) previously downloaded original gene set gff3 '
              f'file path for grch37'))

    parser_dict['optional'].add_argument(
        '--downloaded-file-grch38', dest='download_path_grch38', 
        required=False, 
        help=(f'(Optional) previously downloaded original gene set gff3 '
              f'file path for grch38'))

    workflow.add_logging_args(parser_dict)

    args = parser_dict['main'].parse_args(cmdargs)

    return args


def download(logger, url, download_path):
    logger.info(f'Beginning download of url "{url}" to "{download_path}"')
    network.download(url, download_path)
    logger.info(f'Download finished')


def get_line_lists(download_path, chromdict, asmblspec, output_chrname_version):
    line_list_comment = list()
    linesp_list = list()
    with gzip.open(download_path, 'rt') as infile:
        for line in infile:
            linestrip = tools.rm_newline(line)
            if linestrip == '###':
                continue
            else:
                if linestrip.startswith('#'):
                    line_list_comment.append(linestrip)
                else:
                    linesp = linestrip.split('\t')
                    if linesp[2] in ('chromosome', 'supercontig'):
                        continue
                    if linesp[0] not in chromdict:
                        linesp[0] = asmblspec.convert(
                            linesp[0], output_chrname_version)

                    linesp_list.append(linesp)

    return line_list_comment, linesp_list


def sort_linesp_list(chromdict, linesp_list):
    def sortkey(linesp):
        return tools.coord_sortkey(linesp[0], int(linesp[3]), chromdict)
    linesp_list_sorted = sorted(linesp_list, key = sortkey)

    return linesp_list_sorted


def write_uncompressed_file(outfile_path_uncomp, line_list_comment, 
                            linesp_list_sorted): 
    with open(outfile_path_uncomp, 'wt') as outfile:
        for linestrip in line_list_comment:
            outfile.write(linestrip + '\n')
        for linesp in linesp_list_sorted:
            outfile.write('\t'.join(linesp) + '\n')


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
    line_list_comment, linesp_list = get_line_lists(download_path, chromdict, 
                                                    asmblspec,
                                                    output_chrname_version)

    # sort
    logger.info(f'Sorting the original file lines')
    linesp_list_sorted = sort_linesp_list(chromdict, linesp_list)

    # write uncompressed outfile
    logger.info(f'Writing the sorted lines')
    outfile_path_uncomp = outfile_path + '.uncompressed'
    write_uncompressed_file(outfile_path_uncomp, line_list_comment, linesp_list_sorted)

    # compress and index
    logger.info(f'Compressing with bgzip and indexing with tabix')
    pysam.tabix_compress(outfile_path_uncomp, outfile_path)
    pysam.tabix_index(outfile_path, preset = 'gff')

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

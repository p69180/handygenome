import os
import re
import argparse
import gzip
import ftplib
import shutil

import pysam

import handygenome.network as network
import handygenome.tools as tools
import handygenome.workflow as workflow
import handygenome.workflow.toolsetup as toolsetup
import handygenome.annotation.customfile as customfile
import handygenome.refgenome as refgenome


DEFAULT_OUTFILE_BASENAME = 'ensembl_regulatory_sorted.gff.gz'

MAINFILE_URL_GRCH37 = 'http://ftp.ensembl.org/pub/grch37/release-105/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz'
MAINFILE_URL_GRCH38 = 'http://ftp.ensembl.org/pub/current_regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz'

FTP_HOST = 'ftp.ensembl.org'
FTP_ACTIVITY_DIRECTORY_GRCH37 = '/pub/grch37/release-105/regulation/homo_sapiens/RegulatoryFeatureActivity'
FTP_ACTIVITY_DIRECTORY_GRCH38 = '/pub/current_regulation/homo_sapiens/RegulatoryFeatureActivity'

DEFAULT_OUTFILE_PATH_GRCH37 = os.path.join(
    os.getcwd(), re.sub('gff.gz$', 'grch37.gff.gz', DEFAULT_OUTFILE_BASENAME))
DEFAULT_OUTFILE_PATH_GRCH38 = os.path.join(
    os.getcwd(), re.sub('gff.gz$', 'grch38.gff.gz', DEFAULT_OUTFILE_BASENAME))


def argument_parser(cmdargs):
    def sanity_check(args):
        pass

    parser_dict = workflow.init_parser(
        description=(f'Ensembl gff files for regulatory elements and '
                     f'activity are downloaded and merged into a single '
                     f'gff file.'))

    parser_dict['optional'].add_argument(
        '--outfile-grch37', dest='outfile_path_grch37', 
        default=DEFAULT_OUTFILE_PATH_GRCH37,
        type=workflow.arghandler_outfile_ask,
        required=False, help=f'Output gff file path for grch37')

    parser_dict['optional'].add_argument(
        '--outfile-grch38', dest='outfile_path_grch38', 
        default=DEFAULT_OUTFILE_PATH_GRCH38,
        type=workflow.arghandler_outfile_ask,
        required=False, help=f'Output gff file path for grch38')

    parser_dict['optional'].add_argument(
        '--download-dir-grch37', dest='download_dir_grch37', required=False, 
        help=f'Directory of downloaded files for grch37. Not removed.')

    parser_dict['optional'].add_argument(
        '--download-dir-grch38', dest='download_dir_grch38', required=False, 
        help=f'Directory of downloaded files for grch38. Not removed.')

    workflow.add_logging_args(parser_dict)

    args = parser_dict['main'].parse_args(cmdargs)
    sanity_check(args)

    return args


def download(logger, mainfile_url, ftp_activity_directory, download_dir, 
             mainfile_path, activity_dir):
    # main file
    logger.info(f'Beginning download of main gff to "{mainfile_path}"')
    network.download(mainfile_url, mainfile_path)

    # activity files
    logger.info(f'LOGIN to ensembl ftp host "{FTP_HOST}"')
    ftp = ftplib.FTP(FTP_HOST)
    ftp.login()
    ftp.cwd(ftp_activity_directory)
    tissue_list = ftp.nlst()
    for tissue in tissue_list:
        logger.info(f'Beginning download of activity file for "{tissue}"')
        local_fname = os.path.join(activity_dir, f'{tissue}.gff.gz')
        ftp.cwd(tissue)
        fname = [x for x in ftp.nlst() if x.startswith('homo_sapiens')][0]
        with open(local_fname, 'wb') as outfile:
            ftp.retrbinary(f'RETR {fname}', outfile.write)
        ftp.cwd('..')

    ftp.quit()

    logger.info(f'Download finished')


def load_activity_data(activity_dir, logger):
    activity_data = dict()
    tissue_list = list()
    for fname in tools.listdir(activity_dir):
        tissue = re.sub('\.gff\.gz$', '', os.path.basename(fname))
        tissue_list.append(tissue)
        logger.info(f'Loading {tissue} file')
        with gzip.open(fname, 'rt') as infile:
            for line in infile:
                linesp = tools.get_linesp(line)
                attrs = linesp[8]
                attrs_parsed = customfile.parse_gff3_attrs(attrs)
                ID = attrs_parsed['regulatory_feature_stable_id']
                activity = attrs_parsed['activity']

                activity_data.setdefault(ID, dict())
                activity_data[ID][tissue] = activity

    return activity_data, tissue_list


def load_mainfile(mainfile_path, chromdict, asmblspec, 
                  output_chrname_version):
    mainfile_data = list()
    with gzip.open(mainfile_path, 'rt') as infile:
        for line in infile:
            linesp = tools.get_linesp(line)

            # edit chrom
            if linesp[0] not in chromdict:
                linesp[0] = asmblspec.convert(linesp[0], output_chrname_version)

            # get attributes dictionary
            raw_attrs = customfile.parse_gff3_attrs(linesp[8])
            attrs = {
                'id': raw_attrs['ID'].split(':')[1],
                'bound_end': raw_attrs['bound_end'],
                'bound_start': raw_attrs['bound_start'],
            }

            # remove linesp attributes field
            del linesp[8]

            line_data = {'linesp': linesp, 'attrs': attrs}
            mainfile_data.append(line_data)

    return mainfile_data


def sort_mainfile_data(chromdict, mainfile_data):
    def sortkey(line_data):
        return tools.coord_sortkey(line_data['linesp'][0], int(line_data['linesp'][3]), chromdict)
    mainfile_data_sorted = sorted(mainfile_data, key=sortkey)

    return mainfile_data_sorted


def write_uncompressed_file(outfile_path_uncomp, mainfile_data_sorted, activity_data, tissue_list):
    with open(outfile_path_uncomp, 'wt') as outfile:
        for line_data in mainfile_data_sorted:
            activity = activity_data[line_data['attrs']['id']]
            line_data['attrs']['activity'] = ','.join(
                f'{tissue}:{activity.setdefault(tissue, "NA")}' 
                for tissue in tissue_list)
            attrs_string = ';'.join(
                f'{key}={val}' for key, val in line_data['attrs'].items())
            line_data['linesp'].append(attrs_string)
            outfile.write('\t'.join(line_data['linesp']) + '\n')


def convert_common(logger, mainfile_url, ftp_activity_directory, refver, 
                   output_chrname_version, outfile_path, download_dir):
    # download
    if download_dir is None:
        download_dir = workflow.get_tmpfile_path(delete=False, is_dir=True)
        mainfile_path = os.path.join(download_dir, 'main.gff.gz')
        activity_dir = os.path.join(download_dir, 'activity')
        os.mkdir(activity_dir)
        download(logger, mainfile_url, ftp_activity_directory, download_dir, 
                 mainfile_path, activity_dir)

        rm_downloaded = True
    else:
        mainfile_path = os.path.join(download_dir, 'main.gff.gz')
        activity_dir = os.path.join(download_dir, 'activity')

        rm_downloaded = False

    # set refgenome dbs
    asmblspec = refgenome.SPECS[refver]
    chromdict = asmblspec.chromdicts[output_chrname_version]

    # load activity data
    logger.info(f'Loading activity files')
    activity_data, tissue_list = load_activity_data(activity_dir, logger)

    # load main file
    logger.info(f'Loading the main file')
    mainfile_data = load_mainfile(mainfile_path, chromdict, asmblspec, 
                                  output_chrname_version)

    # sort
    logger.info(f'Sorting loaded main file data')
    mainfile_data_sorted = sort_mainfile_data(chromdict, mainfile_data)

    # write uncompressed outfile
    logger.info(f'Writing the sorted lines')
    outfile_path_uncomp = outfile_path + '.uncompressed'
    write_uncompressed_file(outfile_path_uncomp, mainfile_data_sorted, 
                            activity_data, tissue_list)

    # compress and index
    logger.info(f'Compressing with bgzip and indexing with tabix')
    pysam.tabix_compress(outfile_path_uncomp, outfile_path)
    pysam.tabix_index(outfile_path, preset = 'gff')

    # remove temp files
    os.remove(outfile_path_uncomp)
    if rm_downloaded:
        shutil.rmtree(download_dir)


def convert_grch37(logger, outfile_path, download_dir):
    convert_common(logger, MAINFILE_URL_GRCH37, FTP_ACTIVITY_DIRECTORY_GRCH37, 'hg19', 'nochr_plus_genbank', outfile_path, download_dir)


def convert_grch38(logger, outfile_path, download_dir):
    convert_common(logger, MAINFILE_URL_GRCH38, FTP_ACTIVITY_DIRECTORY_GRCH38, 'hg38', 'ucsc', outfile_path, download_dir)


def main(cmdargs):
    args = argument_parser(cmdargs)
    logger = toolsetup.setup_logger(args)

    #convert_grch37(logger, args.outfile_path_grch37, args.download_dir_grch37)
    convert_grch38(logger, args.outfile_path_grch38, args.download_dir_grch38)
    
    logger.info('ALL SUCCESSFULLY FINISHED')

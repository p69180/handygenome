import os
import re
import contextlib

import handygenome.common as common


URL_AUTHORITY = 'ftp.ncbi.nlm.nih.gov'
URL_PATH_REFSEQ_GENOME = '/genomes/refseq'
URL_PATH_DBSNP = "/snp/latest_release/VCF"
REFSEQ_ACC_PATSTR = r'(?P<refseq_acc>(?P<refseq_acc_prefix>[^_]+)_(?P<refseq_acc_main>[^_]+)\.(?P<refseq_acc_version>[^_]+))'
REFVER_PATSTR = r'(?P<refver>(?P<refver_main>[^.]+)(\.(?P<refver_sub>[^.]+))?)'
GENOME_DIR_PAT = re.compile(REFSEQ_ACC_PATSTR + '_' + REFVER_PATSTR)


def genome_path_sorting(ftp, species_genome_topdir):
    fname_list = common.ftp_listdir(ftp, species_genome_topdir)

    genome_path_dict = dict()
    for path in fname_list:
        mat = GENOME_DIR_PAT.fullmatch(os.path.basename(path))
        if mat is None:
            continue

        refver = mat.group('refver_main')
        acc_ver = int(mat.group('refseq_acc_version'))
        genome_path_dict.setdefault(refver, dict())
        genome_path_dict[refver][acc_ver] = path

    return genome_path_dict


def refver_subver_key(refver_sub):
    if refver_sub is None:
        refver_sub_key = -1
    elif refver_sub.startswith('p'):
        refver_sub_key = int(refver_sub[1:])
    else:
        refver_sub_key = int(refver_sub)
    return refver_sub_key


def pick_latest_genome_path(genome_path_dict):
    result = dict()
    for refver, subdic in genome_path_dict.items():
        result[refver] = max(subdic.items(), key=lambda x: x[0])[1]
    return result


def collect_latest_genome_paths(ftp=None, retry_count=10, retry_interval=1, timeout=5):
    def helper(ftp, species_genome_topdir, result):
        fname_dict = genome_path_sorting(ftp, species_genome_topdir)
        latest_genome_paths_part = pick_latest_genome_path(fname_dict)
        result.update(latest_genome_paths_part)

    # set params
    topdir_human = URL_PATH_REFSEQ_GENOME + '/vertebrate_mammalian/Homo_sapiens/all_assembly_versions'
    topdir_mouse = URL_PATH_REFSEQ_GENOME + '/vertebrate_mammalian/Mus_musculus/all_assembly_versions'
    topdir_banana = URL_PATH_REFSEQ_GENOME + '/plant/Musa_acuminata/all_assembly_versions'

    # access ftp server and collect file paths
    result = dict()
    if ftp is None:
        ftp = common.ftp_login(
            URL_AUTHORITY, 
            retry_count=retry_count, 
            retry_interval=retry_interval, 
            timeout=timeout,
        )

    with contextlib.redirect_stdout('/dev/null'):
        helper(ftp, topdir_human, result)
        helper(ftp, topdir_mouse, result)
        helper(ftp, topdir_banana, result)

    return result


def find_assemblyfile_from_genomepath(ftp, genome_path):
    fname_list = common.ftp_listdir(ftp, genome_path)
    results = [x for x in fname_list if x.endswith('assembly_report.txt')]
    if len(results) != 1:
        raise Exception(f'Cannot find a unique assembly report file')
    return results[0]


def collect_assemblyfile_urls(retry_count=10, retry_interval=1, timeout=5):
    ftp = common.ftp_login(
        URL_AUTHORITY, 
        retry_count=retry_count, 
        retry_interval=retry_interval, 
        timeout=timeout,
    )
    latest_genome_paths = collect_latest_genome_paths(ftp=ftp)
    assemblyfile_paths = {
        refver: find_assemblyfile_from_genomepath(ftp, genome_path)
        for refver, genome_path in latest_genome_paths.items()
    }

    # add prefix
    assemblyfile_urls = {
        key: f'https://{URL_AUTHORITY}' + val
        for key, val in assemblyfile_paths.items()
    }

    # quit ftp
    ftp.quit()

    return assemblyfile_urls


def get_dbsnp_urls(retry_count=10, retry_interval=1, timeout=5):
    # access ftp server and collect file paths
    ftp = common.ftp_login(
        URL_AUTHORITY, 
        retry_count=retry_count, 
        retry_interval=retry_interval, 
        timeout=timeout,
    )
    with contextlib.redirect_stdout('/dev/null'):
        fname_list = common.ftp_listdir(ftp, URL_PATH_DBSNP)
        ftp.quit()

    # make result
    pat = re.compile(r'(.+)\.gz(.*)')
    url_dict = dict()
    for fname in fname_list:
        mat = pat.fullmatch(os.path.basename(fname))
        if mat is None:  # CHECKSUM
            continue

        key = mat.group(1)
        url_dict.setdefault(key, dict())

        if mat.group(2) == '':
            subkey = 'vcf'
        elif mat.group(2) == '.md5':
            subkey = 'vcf_md5'
        else:
            continue

        url_dict[key][subkey] = f'https://{URL_AUTHORITY}' + fname

    return url_dict


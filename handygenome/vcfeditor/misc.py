import os
import gzip

import pysam
import numpy as np

import handygenome.common as common


def get_vcf_format(vcf_path):
    with pysam.VariantFile(vcf_path) as vcf:
        is_vcf = vcf.is_vcf
        #is_bcf = vcf.is_bcf
        comp = vcf.compression

    if comp not in ('NONE', 'GZIP', 'BGZF'):
        raise Exception(f'Unexpected compression type')

    is_bgzf = (comp == 'BGZF')

    return is_vcf, comp, is_bgzf


def check_has_index(vcf_path):
    with pysam.VariantFile(vcf_path) as vcf:
        has_index = (vcf.index is not None)
    return has_index


def get_indexfile_path(vcf_path):
    if os.path.exists(f'{vcf_path}.csi'):
        return f'{vcf_path}.csi'
    elif os.path.exists(f'{vcf_path}.tbi'):
        return f'{vcf_path}.tbi'
    else:
        return None


def make_index(vcf_path):
    index_filename = f'{vcf_path}.csi'
    pysam.tabix_index(vcf_path, preset='vcf', csi=True, index=index_filename)
    return index_filename


def get_vcf_positions(vcf_path, as_iter=False, verbose=False):
    """Returns:
        as_iter==True: A list of tuples (chrom, start0, end0) of each variant record
        as_iter==False: np.ndarray
    """
    is_vcf, comp, is_bgzf = get_vcf_format(vcf_path)
    if is_vcf:
        line_iter = get_line_iter(vcf_path, comp, verbose)
        result_iter = parse_line_iter(line_iter)
    else:  # BCF
        result_iter = get_vcf_positions_pysam(vcf_path)

    if as_iter:
        return result_iter
    else:
        return np.fromiter(result_iter, dtype=object)


def get_vcf_fetchregions(vcf_path, n, refver, verbose=False):
    all_position_info = get_vcf_positions(vcf_path, as_iter=False, verbose=verbose)
    split_position_info = [
        x for x in np.array_split(all_position_info, n) if x.shape[0] != 0
    ]

    result = list()
    for position_info in split_position_info:
        chrom_left = position_info[0][0]
        start0_left = position_info[0][1]
        chrom_right = position_info[-1][0]
        end0_right = position_info[-1][1] + 1
        intvlist = common.IntervalList.from_margin(
            refver, chrom_left, start0_left, chrom_right, end0_right,
        )
        fetchregion_list = [
            (intv.chrom, intv.start0, intv.end0)
            for intv in intvlist
        ]
        result.append(fetchregion_list)

    return result


def get_vcf_lineno(vcf_path, verbose=False):
    is_vcf, comp, is_bgzf = get_vcf_format(vcf_path)
    if is_vcf:
        line_iter = get_line_iter(vcf_path, comp, verbose)
        return count_iterator(line_iter)
    else:  # BCF
        return get_vcf_lineno_pysam(vcf_path)

###

def count_iterator(iterator):
    idx = -1
    for idx, x in enumerate(iterator):
        pass
    return idx + 1


def get_line_iter(vcf_path, comp, verbose):
    if comp in ('BGZF', 'GZIP'):
        raw_iter = gzip.open(vcf_path, 'rb')
    elif comp == 'NONE':
        raw_iter = open(vcf_path, 'rb')

    if verbose:
        for idx, line in enumerate(raw_iter):
            if idx % 100000 == 0:
                print(idx)
            if not line.startswith(b'#'):
                yield line

    else:
        for line in raw_iter:
            if not line.startswith(b'#'):
                yield line


def parse_line_iter(line_iter):
    """Accepts byte mode iterator"""
    for line in line_iter:
        linesp = line.split(b'\t')
        pos0 = int(linesp[1]) - 1
        yield (linesp[0].decode(), pos0, pos0 + len(linesp[3]))


def get_vcf_positions_pysam(vcf_path):
    with pysam.VariantFile(vcf_path) as vcf:
        for vr in vcf.fetch():
            yield (vr.contig, vr.pos - 1, vr.pos - 1 + len(vr.ref))


def get_vcf_lineno_pysam(vcf_path):
    with pysam.VariantFile(vcf_path) as vcf:
        lineno = count_iterator(vcf.fetch())
    return lineno



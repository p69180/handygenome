import os
import gzip

import pysam
import numpy as np

import handygenome.refgenome.refgenome as refgenome
import handygenome.interval as libinterval


# pysam VariantFile mode related ones

"""
pysam file mode string

- Valid mode string pattern : ^[rwa]([bzu]?[0-9]?|[0-9]?[bzu]?)$
- Conflicting characters: (b,z,u)

- 'wb' : Compressed BCF (level 6), regardless of file name. 
- 'wb[0-9]' : BCF with level of compression indicated by the number, regardless of file name. 

- 'wz' : Compressed VCF (level 6), regardless of file name.
- 'wz[0-9]' : VCF with level of compression indicated by the number, regardless of file name. 

- 'wu' : Uncompressed VCF, regardless of file name.
- 'wu[0-9]' : Uncompressed VCF, regardless of file name.

- 'w[0-9]' : Uncompressed VCF, regardless of file name. 

- 'w' :
    *.vcf : uncompressed VCF
    *.vcf.gz : compressed VCF (level 6)
    *.bcf : compressed BCF (level 6)
    *.bcf.gz : compressed VCF (level 6)
"""

CYVCF2_FORMAT_DICT = {'v':'w', 'z':'wz', 'u':'wbu', 'b':'wb'}  
    # keys: bcftools format strings
    # values: cyvcf2 format strings
PYSAM_MODE_DICT = {'v':'wz0', 'z':'wz', 'u':'wb0', 'b':'wb'}
    # keys: bcftools format strings
    # values: pysam format strings
#DEFAULT_MODE_BCFTOOLS = 'z'


def write_mode_arghandler(mode_bcftools, mode_pysam):
    """
    Args:
        mode_bcftools: bcftools style mode string
        mode_pysam: pysam style mode string

    If 'mode_pysam' is set as a non-None value,
    it overrides 'mode_bcftools'.
    """
    if mode_pysam is None:
        return PYSAM_MODE_DICT[mode_bcftools]
    else:
        return mode_pysam

###

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
        #return list(result_iter)
        return np.array(
            list(result_iter), 
            dtype=[('Chromosome', object), ('Start', int), ('End', int)],
        )


def split_vcf_positions(vcf_positions, n):
    """Splits such that adjacent sub-positions do not meet each other with
    identical POS values
    """
    diffresult = np.diff(vcf_positions['Start'])
    nondiff_indexes = np.where(diffresult == 0)[0] + 1
    spittable_indexes = np.setdiff1d(np.arange(vcf_positions.shape[0]), nondiff_indexes)
    spittable_indexes_split = [
        x for x in np.array_split(spittable_indexes, n) 
        if x.shape[0] != 0
    ]
    splitpoints = [x[0] for x in spittable_indexes_split if x[0] != 0]
    split_positions = np.split(vcf_positions, splitpoints)
        # this works even when "splitpoints" is an empty list

    return split_positions


def get_fetchargs_from_vcf_positions(vcf_positions, refver):
    chrom_left = vcf_positions[0]['Chromosome']
    start0_left = vcf_positions[0]['Start']
    chrom_right = vcf_positions[-1]['Chromosome']
    end0_right = vcf_positions[-1]['Start'] + 1  
        # only the 1-base position indicated by POS value is used
    chromdict = refgenome.get_chromdict(refver)

    intvlist = libinterval.IntervalList.from_margin(
        chromdict, chrom_left, start0_left, chrom_right, end0_right,
    )
    fetchargs_list = [
        (intv.chrom, intv.start0, intv.end0) for intv in intvlist
    ]
    return fetchargs_list


def get_vcf_fetchregions(vcf_path, n, refver, verbose=False):
    all_position_info = get_vcf_positions(vcf_path, as_iter=False, verbose=verbose)
    split_position_info = split_vcf_positions(all_position_info, n)
    result = [
        get_fetchargs_from_vcf_positions(position_info, refver)
        for position_info in split_position_info
    ]

    return result


def get_vcf_lineno(vcf_path, verbose=False):
    is_vcf, comp, is_bgzf = get_vcf_format(vcf_path)
    if is_vcf:
        line_iter = get_line_iter(vcf_path, comp, verbose)
        return count_iterator(line_iter)
    else:  # BCF
        return get_vcf_lineno_pysam(vcf_path)


def get_vr_fetcher(vcf, refver, chrom=None, start0=None, end0=None, respect_refregion=False):
    if chrom is None:
        for vr in vcf.fetch():
            yield vr
    else:
        if start0 is None:
            start0 = 0
        if end0 is None:
            end0 = refgenome.get_chromdict(refver)[chrom]
        
        fetcher = vcf.fetch(contig=chrom, start=start0, stop=end0)

        if respect_refregion:
            def check_vr_inclusion(vr, chrom, start0, end0):
                vr_start0 = vr.pos - 1
                vr_end0 = vr_start0 + len(vr.ref)
                return (
                    (vr.contig == chrom)
                    and (vr_end0 > start0)
                    and (vr_start0 < end0)
                )
        else:
            def check_vr_inclusion(vr, chrom, start0, end0):
                vr_start0 = vr.pos - 1
                vr_end0 = vr.pos
                return (
                    (vr.contig == chrom)
                    and (vr_end0 > start0)
                    and (vr_start0 < end0)
                )

        for vr in fetcher:
            if check_vr_inclusion(vr, chrom, start0, end0):
                yield vr

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


# this does not work
def get_vcf_noerr(*args, **kwargs):
    with contextlib.redirect_stderr(io.StringIO()) as err, \
            contextlib.redirect_stdout(io.StringIO()) as out:
        vcf = pysam.VariantFile(*args, **kwargs)

    for buf in (err, out):
        msg = buf.getvalue()
        if not msg.startswith('[E::idx_find_and_load] '
                              'Could not retrieve index file for'):
            print(msg, end='', file=sys.stderr)

    return vcf



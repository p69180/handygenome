import os

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))


def create_header(chromdict):
    return pysam.AlignmentHeader.from_references(
        reference_names=chromdict.contigs, 
        reference_lengths=chromdict.lengths,
    )


def sort_and_index(bam_path):
    tmp_bam_path = workflow.get_tmpfile_path(prefix=bam_path)
    pysam.sort('-O', 'BAM', '-o', tmp_bam_path, bam_path)
    os.rename(tmp_bam_path, bam_path)
    _ = pysam.index(bam_path)


def samtools_view(in_bam_path, out_bam_path, region_bed_path=None, index=True):
    args = ['-h', '-b']
    if region_bed_path is not None:
        args.extend(['-M', '-L', region_bed_path])
    args.append(in_bam_path)
    with open(out_bam_path, mode='wb') as out_bam:
        out_bam.write(pysam.view(*args))
    if index:
        pysam.index(out_bam_path)


def write_readlist(bam_path, readlist, chromdict, sort=True, index=True):
    bam_path = workflow.arghandler_outfile_mustbeabsent(bam_path)
    hdr = create_header(chromdict)
    tmp_bam_path = workflow.get_tmpfile_path(prefix=bam_path)

    with pysam.AlignmentFile(tmp_bam_path, mode='wb', header=hdr) as in_bam:
        for read in readlist:
            in_bam.write(read)

    if sort:
        pysam.sort('-O', 'BAM', '-o', bam_path, tmp_bam_path)
        os.remove(tmp_bam_path)
        if index:
            _ = pysam.index(bam_path)
    else:
        os.rename(tmp_bam_path, bam_path)


def get_readgroup(bam):
    rg = bam.header['RG']
    if len(rg) != 1:
        raise Exception(f'The number of RG is not one.')
    else:
        return rg[0]['ID']


def check_header_compatibility(bamheader_list):
    """Only checks contig order and lengths"""
    contigs = set()
    lengths = set()
    for bamheader in bamheader_list:
        contigs.add(tuple(bamheader.references))
        lengths.add(tuple(bamheader.lengths))

    return (
        len(contigs) == 1 and
        len(lengths) == 1
    )
        
    
def get_average_depth(bam, readlen=None, aligned_length=None):
    if bam.is_cram:
        raise Exception(f'Cannot infer average depth from a cram file.')
    if not bam.has_index():
        raise Exception(f'Bam file must be indexed.')

    if aligned_length is None:
        aligned_length = sum(bam.lengths)
    readcounts = sum(x.mapped for x in bam.get_index_statistics())
    if readlen is None:
        readlen = check_read_length(bam)
        
    return (readcounts * readlen) / aligned_length


def check_read_length(bam):
    """Gets maximum of the first 100 reads"""
    values = list()
    for read in bam.fetch():
        if read.is_unmapped:
            continue

        if 'H' in read.cigarstring:
            continue

        length = read.query_length
        if length == 0:
            continue

        values.append(length)
        if len(values) == 100:
            break

    return max(values)


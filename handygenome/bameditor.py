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


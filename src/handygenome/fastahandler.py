import re
import os
import gzip
import multiprocessing
#import shutil

import pysam
import Bio.SeqIO
import Bio.SeqUtils
import numpy as np

import handygenome.tools as tools
import handygenome.logutils as logutils


def make_index(fasta_path):
    """If input file is bgzipped, *.fai and *.gzi files are created."""
    #fai_path = fasta_path + '.fai'
    #if not os.path.exists(fai_path):
    _ = pysam.faidx(fasta_path)


##################
# fasta renaming #
##################

def _chrom_converter_handler(chrom_converter):
    if isinstance(chrom_converter, dict):
        def converter(x):
            if x in chrom_converter:
                return chrom_converter[x]
            else:
                return None
        return converter
    elif callable(chrom_converter):
        return chrom_converter
    else:
        raise Exception(f'Invalid input type')


def rename_fasta(infile_path, outfile_path, chrom_converter):
    chrom_converter = _chrom_converter_handler(chrom_converter)

    with (
        tools.openfile(infile_path, 'rt') as infile,
        open(outfile_path, 'wt') as outfile,
    ):
        def record_iter():
            for record in Bio.SeqIO.parse(infile, 'fasta'):
                new_id = chrom_converter(record.id)
                if new_id is None:
                    continue

                record.id = new_id
                record.name = None
                record.description = ''
                yield record

        Bio.SeqIO.write(record_iter(), outfile, 'fasta')


def rename_fasta_compress(infile_path, outfile_path, chrom_converter, index=True):
    """Not to be used since fetching from bgzipped fasta is much slower than
    from plain text fasta. (1.68 ms vs 20.5 ms)
    """
    # make target functions
    def pipe_out(pipe_path, outfile_path, index):
        """Read from fifo, bgzip-compress, and write to final output file"""
        pysam.tabix_compress(pipe_path, outfile_path, force=True)
        if index:
            make_index(outfile_path)

    # make fifo
    pipe_path = tools.get_tmpfile_path(prefix=f'{infile_path}_', suffix=f'.fifo', delete=True)
    os.mkfifo(pipe_path)

    # run
    subp1 = multiprocessing.Process(
        target=rename_fasta_nocompress, 
        args=(infile_path, pipe_path, chrom_converter),
    )
    subp2 = multiprocessing.Process(
        target=pipe_out, 
        args=(pipe_path, outfile_path, index),
    )
    subp1.start()
    subp2.start()
    subp1.join()
    subp2.join()

    # remove fifo
    os.remove(pipe_path)

    # raise if erroneous
    if (subp1.exitcode != 0) or (subp2.exitcode != 0):
        raise Exception(f'Finished with an error.')




############
# N region #
############

def make_N_region_coords(fasta):
    pat = re.compile('N+')
    chroms = list()
    start0s = list()
    end0s = list()
    for chrom in fasta.references:
        seq = fasta.fetch(chrom)
        offset = 0
        while True:
            mat = pat.search(seq)
            if mat is None:
                break
            else:
                span = mat.span()
                chroms.append(chrom)
                start0s.append(offset + span[0])
                end0s.append(offset + span[1])

                offset += span[1]
                seq = seq[span[1]:]

    return chroms, start0s, end0s


###############
# GC fraction #
###############

def calculate_gc_fractions_old(chroms, start0s, end0s, fasta, window=None, as_array=True):
    logutils.log('Beginning gc fraction calculation')
    get_fetchargs = _make_get_fetchargs(window, fasta)
    val_gen = (
        float(Bio.SeqUtils.gc_fraction(fasta.fetch(*get_fetchargs(chrom, start0, end0))))
        for chrom, start0, end0 in zip(chroms, start0s, end0s)
    )

    if as_array:
        return np.fromiter(val_gen, float)
    else:
        return list(val_gen)


def calculate_gc_fractions(chroms, start0s, end0s, fasta, window=None, as_array=True):
    logutils.log('Beginning gc fraction calculation')
    chrom_seq_cache = dict()
    get_fetchargs = _make_get_fetchargs(window, fasta)

    def gcfrac_generator():
        for chrom, start0, end0 in zip(chroms, start0s, end0s):
            chrom, start0, end0 = get_fetchargs(chrom, start0, end0)
            if chrom not in chrom_seq_cache:
                chrom_seq_cache[chrom] = fasta.fetch(chrom)
            seq = chrom_seq_cache[chrom][start0:end0]
            yield float(Bio.SeqUtils.gc_fraction(seq))

    if as_array:
        result = np.fromiter(gcfrac_generator(), float)
    else:
        result = list(gcfrac_generator())

    logutils.log('Finished gc fraction calculation')
    return result


def _make_get_fetchargs(window, fasta):
    if window is None:
        def get_fetchargs(chrom, start0, end0):
            return chrom, start0, end0
    else:
        chromlens = dict(zip(fasta.references, fasta.lengths))
        pad_left = int(window / 2)
        if window % 2 == 0:
            pad_right = pad_left
        else:
            pad_right = pad_left + 1

        def get_fetchargs(chrom, start0, end0):
            if end0 - start0 >= window:
                return chrom, start0, end0
            else:
                mid = int((start0 + end0 - 1) / 2)
                new_start0 = max(0, mid - pad_left)
                new_end0 = min(chromlens[chrom], mid + pad_right)
                return chrom, new_start0, new_end0

    return get_fetchargs




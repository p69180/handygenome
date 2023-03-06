import os
import itertools

import Bio.SeqUtils
import pyBigWig
import numpy as np
import pandas as pd
import pyranges as pr

import handygenome.common as common


GCWIG_DIR = os.path.join(common.DATA_DIR, 'gcwig')


# making and saving gc fraction bigwig files

def write_gcwig(outfile_path, fasta, binsize=100, verbose=False):
    bw = pyBigWig.open(outfile_path, 'w')
    # make header
    header = list() 
    for chrom, length in zip(fasta.references, fasta.lengths):
        header.append((chrom, length))
    bw.addHeader(header)
    # write lines
    for chrom, length in zip(fasta.references, fasta.lengths):
        if verbose:
            print(chrom)
        rng = range(0, length, binsize)
        start0s = (x for x in rng)
        end0s = (min(length, x + binsize) for x in rng)
        chroms = itertools.repeat(chrom, len(rng))
        values = calculate_gcvals(chroms, start0s, end0s, fasta, window=None, as_array=False)
        bw.addEntries(chrom, 0, values=values, span=binsize, step=binsize)

    bw.close()


def get_gcwig_path(refver, binsize):
    refver = common.RefverDict.standardize(refver)
    return os.path.join(GCWIG_DIR, refver, f'gc_binsize{100}.bigwig')


def get_gcwig(refver, binsize):
    gcfile_path = get_gcwig_path(refver, binsize)
    if not os.path.exists(gcfile_path):
        print(f'There is no pre-existing gc bigwig file. A new one is being created. It may take a few minutes.')
        fasta = common.DEFAULT_FASTAS[refver]
        write_gcwig(gcfile_path, fasta, binsize=binsize, verbose=False)
    return pyBigWig.open(gcfile_path)


###########################################

# getting array of gc fraction values

def make_fetchargs_func(window, fasta):
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
                mid = int((start0 + end0) / 2)
                new_start0 = max(0, mid - pad_left)
                new_end0 = min(chromlens[chrom], mid + pad_right)
                return chrom, new_start0, new_end0

    return get_fetchargs


def calculate_gcvals(chroms, start0s, end0s, fasta, window=None, as_array=False):
    get_fetchargs = make_fetchargs_func(window, fasta)
    val_gen = (
        float(Bio.SeqUtils.gc_fraction(fasta.fetch(*get_fetchargs(chrom, start0, end0))))
        for chrom, start0, end0 in zip(chroms, start0s, end0s)
    )

    if as_array:
        return np.fromiter(val_gen, float)
    else:
        return list(val_gen)


def load_gcvals(chroms, start0s, end0s, refver, binsize, fasta, window=None, as_array=False):
    get_fetchargs = make_fetchargs_func(window, fasta)
    gcwig = get_gcwig(refver, binsize)
    val_gen = (
        gcwig.stats(
            *get_fetchargs(chrom, start0, end0), 
            type='mean', 
            exact=True,
        )[0]
        for chrom, start0, end0 in zip(chroms, start0s, end0s)
    )

    if as_array:
        return np.fromiter(val_gen, float)
    else:
        return list(val_gen)


# add gc fraction values to a depth dataframe

def add_gc_calculating(df, fasta, window=None):
    """Args:
        df: pandas.DataFrame or pyranges.PyRanges
    Changes in-place
    """
    gcvals = calculate_gcvals(df.Chromosome, df.Start, df.End, fasta, window=window)

    if isinstance(df, pd.DataFrame):
        df['GC'] = gcvals
    elif isinstance(df, pr.PyRanges):
        df.GC = gcvals


def add_gc_loading(df, refver, binsize, window=None):
    """Args:
        df: pandas.DataFrame or pyranges.PyRanges
    Changes in-place
    """
    gcvals = load_gcvals(
        df.Chromosome, df.Start, df.End, 
        refver=refver, binsize=binsize, 
        fasta=common.DEFAULT_FASTAS[refver], 
        window=window, as_array=False,
    )
    if isinstance(df, pd.DataFrame):
        df['GC'] = gcvals
    elif isinstance(df, pr.PyRanges):
        df.GC = gcvals



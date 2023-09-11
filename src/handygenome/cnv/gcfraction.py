"""This version abandons gc data saving with bigwig and saves DataFrame as tsv file"""

import os
import itertools
import functools

import Bio.SeqUtils
import numpy as np
import pandas as pd
import pyranges as pr

import handygenome
import handygenome.refgenome.refgenome as refgenome
import handygenome.logutils as logutils
import handygenome.pyranges_helper as pyranges_helper
import handygenome.deco as deco
import handygenome.refgenome.refgenome as refgenome


GCDATA_DIR = os.path.join(handygenome.USERDATA_DIR, 'gcdata')
os.makedirs(GCDATA_DIR, exist_ok=True)
#if not os.path.exists(GCDATA_DIR):
#    os.mkdir(GCDATA_DIR)


# making and saving gc fraction bigwig files

def write_gcfile(outfile_path, fasta, binsize=100):
    chromdict = refgenome.ChromDict(fasta=fasta)
    bin_gr = chromdict.to_gr().window(binsize)
    gc_df = bin_gr.df
    gc_df['GC'] = calculate_gcvals(
        tuple(bin_gr.Chromosome), 
        tuple(bin_gr.Start), 
        tuple(bin_gr.End), 
        fasta, 
        window=None, 
        as_array=True,
    )
    gc_df.to_csv(outfile_path, sep='\t', header=True, index=False)


def get_gcfile_path(refver, binsize):
    refver = refgenome.RefverDict.standardize(refver)
    gcdata_refver_dir = os.path.join(GCDATA_DIR, refver)
    if not os.path.exists(gcdata_refver_dir):
        os.mkdir(gcdata_refver_dir)
    return os.path.join(gcdata_refver_dir, f'gc_binsize_{100}.tsv.gz')


@functools.cache
def get_gc_df(refver, binsize, coords_as_index=True):
    gcfile_path = get_gcfile_path(refver, binsize)

    if not os.path.exists(gcfile_path):
        logutils.print_timestamp(f'There is no pre-existing gc data file. A new one is being created. It may take a few minutes.')
        fasta = refgenome.get_fasta(refver)
        write_gcfile(gcfile_path, fasta, binsize=binsize)
        logutils.print_timestamp(f'Finished making a gc data file.')

    if coords_as_index:
        index_col = ['Chromosome', 'Start', 'End']
    else:
        index_col = None

    return pd.read_csv(
        gcfile_path, 
        sep='\t', 
        header=0, 
        dtype={'Chromosome': 'category', 'Start': int, 'End': int, 'GC': float},
        index_col=index_col,
    )


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
                mid = int((start0 + end0 - 1) / 2)
                new_start0 = max(0, mid - pad_left)
                new_end0 = min(chromlens[chrom], mid + pad_right)
                return chrom, new_start0, new_end0

    return get_fetchargs


def calculate_gcvals_generator(chroms, start0s, end0s, fasta, get_fetchargs):
    prev_chrom = ''
    for chrom, start0, end0 in zip(chroms, start0s, end0s):
        if chrom != prev_chrom:
            print(chrom)
        prev_chrom = chrom
        yield float(Bio.SeqUtils.gc_fraction(fasta.fetch(*get_fetchargs(chrom, start0, end0))))


@functools.cache
def calculate_gcvals(chroms, start0s, end0s, fasta, window=None, as_array=True):
    logutils.print_timestamp('Beginning gc fraction calculation')
    get_fetchargs = make_fetchargs_func(window, fasta)
    #val_gen = calculate_gcvals_generator(chroms, start0s, end0s, fasta, get_fetchargs)
    val_gen = (
        float(Bio.SeqUtils.gc_fraction(fasta.fetch(*get_fetchargs(chrom, start0, end0))))
        for chrom, start0, end0 in zip(chroms, start0s, end0s)
    )

    if as_array:
        return np.fromiter(val_gen, float)
    else:
        return list(val_gen)


def calculate_gcvals_with_df(df, fasta, window=None, as_array=True):
    return calculate_gcvals(
        tuple(df.Chromosome), 
        tuple(df.Start), 
        tuple(df.End), 
        fasta, 
        window=window, 
        as_array=as_array,
    )


def load_gcvals(chroms, start0s, end0s, refver, binsize, window=None, fasta=None):
    """Returns:
        A pandas Series
    """
    assert not ((window is not None) and (fasta is None))

    get_fetchargs = make_fetchargs_func(window, fasta)
    gc_df = get_gc_df(refver, binsize)
    left_gr = pr.PyRanges(chromosomes=chroms, starts=start0s, ends=end0s)
    joined_gr = pyranges_helper.join(left_gr, pr.PyRanges(gc_df), merge='weighted_mean', how='left')

    return joined_gr.GC


def load_gcvals_with_df(df, refver, binsize, fasta, window=None):
    return load_gcvals(
        df.Chromosome, df.Start, df.End, 
        refver=refver, 
        binsize=binsize, 
        fasta=refgenome.get_fasta(refver), 
        window=window,
    )


# add gc fraction values to a depth dataframe

@deco.get_deco_num_set_differently(('fasta', 'refver'), 1)
def add_gc_calculating(df, *, refver=None, fasta=None, window=None):
    """Args:
        df: pandas.DataFrame or pyranges.PyRanges
        refver, fasta: refver is recommended because with refver, default 
            fasta object with the given refver is loaded and subsequently
            used for "calculate_gcvals" function. Using a fasta object with
            a fixed id can make use of "cache" functionality of "calculate_gcvals".

    Modifies input DataFrame
    """
    if fasta is None:
        fasta = refgenome.get_fasta(refver)

    gcvals = calculate_gcvals(
        tuple(df.Chromosome), 
        tuple(df.Start), 
        tuple(df.End), 
        fasta, 
        window=window, 
        as_array=True,
    )

    if isinstance(df, pd.DataFrame):
        df['GC'] = gcvals
    elif isinstance(df, pr.PyRanges):
        df.GC = gcvals

    return df


def add_gc_loading(df, refver, binsize, window=None):
    """Args:
        df: pandas.DataFrame or pyranges.PyRanges
    Changes in-place
    """
    gcvals = load_gcvals(
        df.Chromosome, df.Start, df.End, 
        refver=refver, 
        binsize=binsize, 
        fasta=refgenome.get_fasta(refver), 
        window=window,
    )
    if isinstance(df, pd.DataFrame):
        df['GC'] = gcvals
    elif isinstance(df, pr.PyRanges):
        df.GC = gcvals



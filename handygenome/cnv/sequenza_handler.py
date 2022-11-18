import os
import re
import tempfile
import subprocess

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
import pyranges as pr

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
read_cnvfile = importlib.import_module('.'.join([top_package_name, 'cnv', 'read_cnvfile']))


SCRIPTS_DIR = os.path.join(common.PROJECT_PATH, 'misc')
CONDA_WRAPPER = os.path.join(SCRIPTS_DIR, 'conda_wrapper.sh')
RUN_EXTRACT = os.path.join(SCRIPTS_DIR, 'sequenza_extract.R')
RUN_FIT = os.path.join(SCRIPTS_DIR, 'sequenza_fit.R')
WRITE_EXTRACT_TSV = os.path.join(SCRIPTS_DIR, 'sequenza_extract_to_tsv.R')
WRITE_FIT_TSV = os.path.join(SCRIPTS_DIR, 'sequenza_fit_to_tsv.R')
MAKE_SEGMENTS = os.path.join(SCRIPTS_DIR, 'sequenza_results_segment.R')

CONDA_WRAPPER_ARGS = [
    CONDA_WRAPPER,  
    '--condadir', '/home/users/pjh/tools/miniconda/220718_for_sequenza/miniconda3/',
    '--envname', 'sequenza',
    'Rscript',
]


def make_segments_as_df(extractfile_path, cellularity, ploidy, is_female):
    fd, tsv_path = tempfile.mkstemp(
        prefix='tmpfile_sequenza_segments_', suffix='.tsv', dir=os.getcwd(),
    )
    os.close(fd)
    make_segments(extractfile_path, cellularity, ploidy, is_female, tsv_path)
    df = read_cnvfile.read_sequenza_segment(tsv_path)
    os.remove(tsv_path)
    return df


def make_segments(extractfile_path, cellularity, ploidy, is_female, outfile_path):
    args = list()
    args.extend(CONDA_WRAPPER_ARGS)
    args.extend(
        [
            MAKE_SEGMENTS,
            '--extfile', extractfile_path,
            '--outfile', outfile_path,
            '--cellularity', str(cellularity),
            '--ploidy', str(ploidy),
            '--gender', ('F' if is_female else 'M'),
        ]
    )
    p = subprocess.run(args, capture_output=True, text=True, check=True)


def load_fit(fitfile_path):
    fd, tsv_path = tempfile.mkstemp(
        prefix='tmpfile_sequenza_fit_result_', suffix='.tsv', dir=os.getcwd(),
    )
    os.close(fd)

    args = list()
    args.extend(CONDA_WRAPPER_ARGS)
    args.extend(
        [
            WRITE_FIT_TSV,
            '--fitfile', fitfile_path,
            '--outfile', tsv_path,
        ]
    )
    p = subprocess.run(args, capture_output=True, text=True, check=True)

    df = pd.read_table(tsv_path, sep='\t', index_col=0)
    df.columns = df.columns.astype('float')
    df.columns.name = 'cellularity'
    df.index.name = 'ploidy'

    os.remove(tsv_path)

    return df


def load_extract(extractfile_path):
    # make temp file
    fd, tsv_path = tempfile.mkstemp(
        prefix='tmpfile_sequenza_extract_result_', suffix='.tsv', dir=os.getcwd(),
    )
    os.close(fd)

    # run R script
    args = list()
    args.extend(CONDA_WRAPPER_ARGS)
    args.extend(
        [
            WRITE_EXTRACT_TSV,
            '--extractfile', extractfile_path,
            '--outfile', tsv_path,
        ]
    )
    p = subprocess.run(args, capture_output=True, text=True, check=True)

    # load dataframe
    df = pd.read_table(tsv_path, sep='\t', header=0, dtype={'chrom': 'string'})
    df.rename(columns={'start': 'Start', 'end': 'End', 'chrom': 'Chromosome'}, inplace=True)

    os.remove(tsv_path)

    return df


def get_1d_peaks(row_iterator):
    result = list()
    for row in row_iterator:
        row_peaks = list()
        peaks_result, _ = scipy.signal.find_peaks(row)
        diff = np.diff(row)
        for within_row_idx in peaks_result:
            row_peaks.append(within_row_idx)
            # leftward
            current_idx = within_row_idx
            while True:
                current_idx -= 1
                if diff[current_idx] != 0:
                    break
                else:
                    row_peaks.append(current_idx)

                if current_idx == 0:
                    break
            # rightward
            current_idx = within_row_idx - 1
            while True:
                current_idx += 1
                if diff[current_idx] != 0:
                    break
                else:
                    row_peaks.append(current_idx + 1)

                if current_idx == len(diff) - 1:
                    break

        row_peaks.sort()
        result.append(row_peaks)
        
    return result


def get_fitresult_peaks(fitresult):
    row_peaks = get_1d_peaks((x[1] for x in fitresult.iterrows()))
    col_peaks = get_1d_peaks((x[1] for x in fitresult.items()))

    row_peaks_tuples = list()
    for row_idx, peaks in enumerate(row_peaks):
        for col_idx in peaks:
            row_peaks_tuples.append((row_idx, col_idx))

    col_peaks_tuples = list()
    for col_idx, peaks in enumerate(col_peaks):
        for row_idx in peaks:
            col_peaks_tuples.append((row_idx, col_idx))

    result = set(row_peaks_tuples).intersection(set(col_peaks_tuples))
    result = [
        {'ploidy': fitresult.index[x], 'cellularity': fitresult.columns[y], 'lpp': fitresult.iloc[x, y]}
        for (x, y) in result
    ]
    result.sort(key=(lambda x: x['lpp']), reverse=True)
    
    return result



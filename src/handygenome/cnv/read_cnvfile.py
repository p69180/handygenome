import gzip
import subprocess

import pandas as pd
import pyranges as pr

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


#def read_sequenza_segment(segment_path, as_gr=False):
#    df = pd.read_table(
#        segment_path, 
#        names=(
#            'Chromosome', 'Start', 'End',
#            'Bf', 'N.BAF', 'sd.BAF',
#            'depth.ratio', 'N.ratio', 'sd.ratio',
#            'CNt', 'A', 'B', 'LPP',
#        ),
#        skiprows=1,
#    )
#    df.loc[df['CNt'] == 0, 'A'] = 0
#    df.loc[df['CNt'] == 0, 'B'] = 0
#
#    if as_gr:
#        return pr.PyRanges(df=df)
#    else:
#        return df


def read_depth_bed(bed_path, as_gr=False):
    """Assumes that the 4th column represents average depth"""
    df = pd.read_table(
        bed_path, 
        names=('Chromosome', 'Start', 'End', 'Depth'),
        skiprows=0,
    )

    if as_gr:
        return pr.PyRanges(df=df)
    else:
        return df


def read_seqz_baf(seqz_path):
    pass


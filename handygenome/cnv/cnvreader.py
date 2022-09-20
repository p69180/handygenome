import pandas as pd
import pyranges as pr

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))

def read_sequenza_segment(segment_path):
    df = pd.read_table(
        segment_path, 
        names=(
            'Chromosome', 'Start', 'End',
            'Bf', 'N.BAF', 'sd.BAF',
            'depth.ratio', 'N.ratio', 'sd.ratio',
            'CNt', 'A', 'B', 'LPP',
        ),
        skiprows=1,
    )
    gr = pr.PyRanges(df=df)
    return gr


def read_depth_bed(bed_path):
    """Assumes that the 4th column represents average depth"""
    df = pd.read_table(
        bed_path, 
        names=('Chromosome', 'Start', 'End', 'Depth'),
        skiprows=0,
    )
    gr = pr.PyRanges(df=df)
    return gr


def read_seqz(seqz_path):
    pass

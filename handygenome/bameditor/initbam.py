import os
import itertools
import warnings

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


def create_header(chromdict):
    return pysam.AlignmentHeader.from_references(
        reference_names=chromdict.contigs, 
        reference_lengths=chromdict.lengths,
    )

import os

import pysam

import handygenome.refgenome as refgenome


def make_index(fasta_path):
    fai_path = fasta_path + '.fai'
    if not os.path.exists(fai_path):
        _ = pysam.faidx(fasta_path)


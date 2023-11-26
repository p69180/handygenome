import os
import tempfile
import subprocess

import Bio.SeqIO
import Bio.SeqRecord
import Bio.Seq
import pysam

import handygenome
import handygenome.refgenome.refgenome as refgenome


def run_bwa_withseqs(seqlist, refver, namelist=None):
    ref_path = refgenome.get_fasta_path(refver)
    if namelist is None:
        seqrec_list = [Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq(seq), id='query')
                       for seq in seqlist]
    else:
        seqrec_list = list()
        for seq, seqname in zip(seqlist, namelist):
            seqrec = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq(seq), id=seqname)
            seqrec_list.append(seqrec)
    
    with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdir:
        query_path = os.path.join(tmpdir, 'input.fasta')
        sam_path = os.path.join(tmpdir, 'output.sam')
        Bio.SeqIO.write(seqrec_list, query_path, 'fasta')
        p = subprocess.run(
            [
                handygenome.PARAMS['bwa'], 'mem', 
                '-Y', 
                '-M', 
                '-t', '1',
                '-o', sam_path,
                ref_path,
                query_path,
            ]
        )
        readlist = list(pysam.AlignmentFile(sam_path).fetch())

    return readlist


def make_index(fasta_path):
    subprocess.run(
        [handygenome.PARAMS['bwa'], 'index', '-a', 'bwtsw', fasta_path]
    )


def check_index(fasta_path):
    suffixes = ['pac', 'bwt', 'ann', 'amb', 'sa']
    return all(os.path.exists(f'{fasta_path}.{suf}') for suf in suffixes)



import os
import tempfile
import subprocess

import Bio.SeqIO
import Bio.SeqRecord
import Bio.Seq
import pysam

import handygenome
import handygenome.network as network
import handygenome.refgenome as refgenome


BLAT_URL = 'https://genome.ucsc.edu/cgi-bin/hgBlat'
BWA = os.path.join(handygenome.DIRS['externals'], 'bwa')


def run_bwa(seqlist, refver, namelist=None):
    ref_path = refgenome.get_default_fasta_path(refver)
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
        p = subprocess.run([
            BWA, 
            'mem', 
            '-Y', 
            '-M', 
            '-t', '1',
            '-o', sam_path,
            ref_path,
            query_path,
        ])
        readlist = list(pysam.AlignmentFile(sam_path).fetch())

    return readlist


def blat(query, refver):
    """Reference: https://genome.ucsc.edu/FAQ/FAQblat.html#blat14
    """
    refver = refgenome.RefverDict.converter[refver]

    result = network.http_get(
        url=BLAT_URL,
        params={
            'userSeq': query,
            'type': 'DNA',
            'db': refver,
            'output': 'json'
        }
    )
    
    if len(result['blat']) == 1:
        return dict(zip(result['fields'], result['blat'][0]))
    if len(result['blat'] == 0):
        return None
    elif len(result['blat'] > 1):
        raise Exception(f'Unexpected blat result:\n{result}')

import os
import tempfile
import subprocess

import Bio.SeqIO
import Bio.SeqRecord
import Bio.Seq
import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


BLAT_URL = 'https://genome.ucsc.edu/cgi-bin/hgBlat'


def run_bwa(seqlist, refver, namelist=None):
    ref_path = common.DEFAULT_FASTA_PATHS[refver]
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
        p = subprocess.run([common.BWA, 'mem', 
                            '-Y', 
                            '-M', 
                            '-t', '1',
                            '-o', sam_path,
                            ref_path,
                            query_path])
        readlist = list(pysam.AlignmentFile(sam_path).fetch())

    return readlist


def blat(query, refver):
    """Reference: https://genome.ucsc.edu/FAQ/FAQblat.html#blat14
    """
    refver = common.RefverDict.converter[refver]

    result = common.http_get(
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

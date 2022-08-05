"""Multiple Sequence Alignment"""

import os
import tempfile
import subprocess
import gzip
import collections
import re

import Bio.SeqIO
import Bio.SeqRecord
import Bio.AlignIO
import numpy as np

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


MUSCLE_PATH = '/home/users/pjh/tools/miniconda/210821/miniconda3/envs/genome_v5/bin/muscle'
PAT_ALIGN = re.compile('^([^-]+)((-+)(([^-])(.*))?)?$')

    
class MultipleSequenceAlignment:
    def __init__(self, seqs, names=None):
        assert len(set(len(x) for x in seqs)) == 1
        
        shape = (len(seqs), len(seqs[0]), 1)
        self.array = np.reshape([list(x) for x in seqs], shape)
        if names is None:
            self.names = None
        else:
            self.names = list(names)

    def __repr__(self):
        sep = ' ' * 4

        lines = list()
        idx_width = len(str(self.array.shape[0] - 1))
        for idx in range(self.array.shape[0]):
            seq = ''.join(self.array[idx].ravel())
            zidx = f'{idx:{idx_width}d}'
            if self.names is None:
                line = sep.join([zidx, seq])
            else:
                name = self.names[idx]
                line = sep.join([zidx, seq, name])
            lines.append(line)

        return '\n'.join(lines)

    def select_col(self, col_idx, rows_incl=None, rows_excl=None):
        if rows_incl is not None:
            row_selector = np.repeat(False, self.array.shape[0])
            row_selector[rows_incl] = True
        elif rows_excl is not None:
            row_selector = np.repeat(True, self.array.shape[0])
            row_selector[rows_excl] = False
        else:
            row_selector = np.repeat(True, self.array.shape[0])

        return self.array[row_selector, col_idx, 0]

    def get_consensus_col(self, col_idx, threshold=0.5, ambiguous='X'):
        col = self.array[:, col_idx, 0]
        col = col[col != '-']
        if col.shape[0] == 0:
            result = ('-')
        else:
            counter = collections.Counter(col)
            max_val = max(counter.values())
            max_bases = [x for x in counter.keys() if counter[x] == max_val]
            if len(max_bases) == 1:
                max_base = max_bases[0]
                total = sum(counter.values())
                prop = counter[max_base] / total
                if prop >= threshold:
                    result = max_base
                else:
                    result = ambiguous
            else:
                result = ambiguous

        return result

    def get_consensus(self, threshold=0.5, ambiguous='X'):
        result = list()
        for col_idx in range(self.array.shape[1]):
            result_col = self.get_consensus_col(col_idx, threshold=threshold, ambiguous=ambiguous)
            result.append(result_col)

        return ''.join(result)

    def get_consensus_biopython(self):
        seqr_list = (Bio.SeqRecord.SeqRecord(''.join(self.array[idx, :, 0]))
                     for idx in range(self.array.shape[0]))
        msa_bio = Bio.Align.MultipleSeqAlignment(seqr_list)
        summ = Bio.Align.AlignInfo.SummaryInfo(msa_bio)

        return str(summ.dumb_consensus())


def modify_msa_gapclose(msa, feasibility_threshold=0.5):
    def check_feasible(msa, row_idx, mat, feasibility_threshold):
        if mat.group(2) is None or mat.group(4) is None:
            return False

        moved_base = mat.group(5)
        col_idx = mat.start(2)
        col_bases = msa.select_col(col_idx, rows_excl=[row_idx])
        base_counter = collections.Counter(col_bases[col_bases != '-'])
        
        numerator = base_counter[moved_base]
        denominator = sum(base_counter.values())
        if denominator == 0:
            return False
        else:
            fraction = numerator / denominator
            return (fraction >= feasibility_threshold)

    def modify_row(msa, row_idx, feasibility_threshold):
        while True:
            mat = PAT_ALIGN.match(''.join(msa.array[row_idx, :, 0]))
            if check_feasible(msa, row_idx, mat, feasibility_threshold):
                new_seq = (mat.group(1) 
                           + mat.group(5) 
                           + mat.group(3) 
                           + mat.group(6))
                msa.array[row_idx, :, 0] = list(new_seq)
                continue
            else:
                break

    # main
    for idx in range(msa.array.shape[0]):
        modify_row(msa, idx, feasibility_threshold)


def fasta_to_msa(fasta_path):
    seqs = list()
    names = list()

    with common.openfile(fasta_path) as infile:
        for seqrec in Bio.SeqIO.parse(infile, 'fasta'):
            seqs.append(str(seqrec.seq))
            names.append(seqrec.id)

    msa = MultipleSequenceAlignment(seqs=seqs, names=names)

    return msa


def write_fasta(fname, seqlist, seqnamelist=None):
    if fname.endswith('.gz'):
        f = gzip.open(fname, 'wt')
    else:
        f = open(fname, 'wt')

    for seq, seqname in zip(seqlist, seqnamelist):
        f.write('>' + seqname + '\n')
        f.write(seq + '\n')

    f.close()


def run_muscle(infile_path, outfile_path):
    args = [MUSCLE_PATH, '-align', infile_path, '-output', outfile_path]
    p = subprocess.run(args, capture_output=True, text=True, check=True)
    return p


def align_with_muscle(seqlist, seqnamelist=None, biopython=False):
    if seqnamelist is None:
        seqnamelist = [f'seq{x}' for x in range(len(seqlist))]

    with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdir:
        infile_path = os.path.join(tmpdir, 'input.fasta')
        outfile_path = os.path.join(tmpdir, 'output.afa')
        write_fasta(infile_path, seqlist, seqnamelist=seqnamelist)
        p = run_muscle(infile_path, outfile_path)

        if biopython:
            msa = Bio.AlignIO.read(outfile_path, 'fasta')
        else:
            msa = fasta_to_msa(outfile_path)

    return msa

                

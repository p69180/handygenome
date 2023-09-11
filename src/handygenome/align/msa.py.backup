"""Multiple Sequence Alignment"""

import os
import tempfile
import subprocess
import gzip
import collections
import re
import itertools

import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
import Bio.AlignIO
import Bio.Align
import Bio.Align.AlignInfo
import numpy as np
import pysam

import handygenome.tools as tools
import handygenome.align.alignhandler as alignhandler


GAP_CHAR = '-'
MUSCLE_PATH = '/home/users/pjh/tools/miniconda/210821/miniconda3/envs/genome_v5/bin/muscle'
PAT_ALIGN = re.compile('^([^-]+)((-+)(([^-])(.*))?)?$')


#################################################################################
# custom-made msa (targeted for softclip segments with identical anchor points) #
#################################################################################

def msa_with_target(target, alignments, target_seqrec=None, query_seqrecs=None):
    # arg handling
    assert (target_seqrec is None) == (query_seqrecs is None)
    seqrec_given = (target_seqrec is not None)

    # make target gaps
    target_gaps = summarize_target_gaps(alignments)

    # make msa
    if seqrec_given:
        input_seqrecs = list()
        input_seqrecs.append(
            add_gaps_to_target(target_seqrec, target_gaps)
        )
        for aln, query_rec in zip(alignments, query_seqrecs):
            input_seqrecs.append(
                add_gaps_to_query(aln, target_gaps, query_seqrec=query_rec)
            )
    else:
        padded_seqs = list()
        padded_seqs.append(
            add_gaps_to_target(target, target_gaps)
        )
        for aln in alignments:
            padded_seqs.append(
                add_gaps_to_query(aln, target_gaps)
            )
        input_seqrecs = [Bio.SeqRecord.SeqRecord(x) for x in padded_seqs]

    msa = Bio.Align.MultipleSeqAlignment(input_seqrecs)

    # make consensus
    consensus = str(Bio.Align.AlignInfo.SummaryInfo(msa).dumb_consensus())

    return msa, consensus


def summarize_target_gaps(aln_list):
    """Make a summary of target gaps in each alignment
    Returns: dict 
        keys: pos0 of the base on the right side of the gap. May include a value equal to "target length"
        values: gap length
    """
    target_gaps = dict()
    for aln in aln_list:
        walks = alignhandler.get_walks(aln)
        for walk in walks:
            if walk.check_queryonly():
                gap_right_border = walk.target.stop
                gap_length = len(walk.query)
                target_gaps.setdefault(gap_right_border, 0)
                target_gaps[gap_right_border] = max(target_gaps[gap_right_border], gap_length)

    return target_gaps


def add_gaps_to_target(target, target_gaps, gap_char=GAP_CHAR):
    if isinstance(target, str):
        return add_gaps_to_target_string(target, target_gaps, gap_char=gap_char)
    elif isinstance(target, Bio.SeqRecord.SeqRecord):
        padded_target = add_gaps_to_target_string(str(target.seq), target_gaps, gap_char=gap_char)
        target.seq = Bio.Seq.Seq(padded_target)
        #new_seqrec = Bio.SeqRecord.SeqRecord(
        #    padded_target, 
        #    annotations=target.annotations,
        #)
        return target
    else:
        raise Exception(f'"target" argument must be either str or Bio.SeqRecord.SeqRecord')


def add_gaps_to_target_string(target, target_gaps, gap_char=GAP_CHAR):
    """Makes a modified target sequence with gap added as '-'"""
    target_list = list(target)
    target_end_idx = len(target)

    for gap_end, gap_len in target_gaps.items():
        if gap_end == target_end_idx:
            target_list.append(gap_char * gap_len)
        else:
            target_list[gap_end] = (gap_char * gap_len) + target_list[gap_end]

    return ''.join(target_list)


def add_gaps_to_query(alignment, target_gaps, query_seqrec=None, overlapping_gap_to_right=True, gap_char=GAP_CHAR):
    if query_seqrec is None:
        return add_gaps_to_query_string(alignment, target_gaps, overlapping_gap_to_right=overlapping_gap_to_right, gap_char=gap_char)
    else:
        padded_query = add_gaps_to_query_string(alignment, target_gaps, overlapping_gap_to_right=overlapping_gap_to_right, gap_char=gap_char)
        query_seqrec.seq = Bio.Seq.Seq(padded_query)
        #new_seqrec = Bio.SeqRecord.SeqRecord(padded_query, annotations=query_seqrec.annotations)
        return query_seqrec


def add_gaps_to_query_string(alignment, target_gaps, overlapping_gap_to_right=True, gap_char=GAP_CHAR):
    """Given an alignment between target and query, returns a modified query sequence with gap added as '-'"""
    target_indices = alignment.indices[0]
    target_indices_list = list(target_indices)

    original_gap_info = list()
    for is_gap, subiter in itertools.groupby(
        enumerate(target_indices), key=(lambda x: x[1] == -1)
    ):
        if is_gap:
            subiter = list(subiter)
            original_gapend_alnidx = subiter[-1][0] + 1
            original_gaplen = len(subiter)
            original_gap_info.append((original_gapend_alnidx, original_gaplen))

    raw_query_list = list(alignment[1])
    raw_query_list.append('')
    for (new_gapend_idx, new_gaplen) in target_gaps.items():
        if new_gapend_idx == len(alignment.target):
            new_gapend_alnidx = len(target_indices)
        else:
            new_gapend_alnidx = target_indices_list.index(new_gapend_idx)

        ovlp_original_gaps = [
            (original_gapend_alnidx, original_gaplen) 
            for (original_gapend_alnidx, original_gaplen) in original_gap_info 
            if (new_gapend_alnidx == original_gapend_alnidx)
        ]
        assert len(ovlp_original_gaps) <= 1

        if len(ovlp_original_gaps) == 0:
            query_insertion_idx = new_gapend_alnidx
            added_gaplen = new_gaplen
        else:
            original_gapend_alnidx, original_gaplen = ovlp_original_gaps[0]

            if overlapping_gap_to_right:
                query_insertion_idx = new_gapend_alnidx
            else:
                query_insertion_idx = new_gapend_alnidx - original_gaplen

            added_gaplen = new_gaplen - original_gaplen
            assert added_gaplen >= 0    

        raw_query_list[query_insertion_idx] = (gap_char * added_gaplen) + raw_query_list[query_insertion_idx]

    return ''.join(raw_query_list)

##########################################################




############################################################
############################################################
############################################################
############################################################
############################################################

# old ones #

    
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

    with tools.openfile(fasta_path) as infile:
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

                

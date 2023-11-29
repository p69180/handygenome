import itertools

import handygenome.tools as tools


class Repeat:
    def __init__(self, chrom, start0, unit, count, refver):
        self.chrom = chrom
        self.start0 = start0
        self.unit = unit
        self.count = count
        self.refver = refver

    def __repr__(self):
        string = tools.repr_base(self, ('chrom', 'start0', 'unit', 'count'), comma_sep_int=True)
        return f'<Repeat object({string})>'

    @property
    def chromdict(self):
        return refgenome.get_chromdict(self.refver)

    @property
    def fasta(self):
        return refgenome.get_fasta(self.refver)

    @property
    def end0(self):
        return self.start0 + (len(self.unit) * self.count)


def decompose_repeat(seq):
    if len(seq) == 1:
        return [(seq, 1)]
    else:
        result = list()
        seqlen = len(seq)
        for sublen in itertools.chain(range(1, int(seqlen / 2) + 1), [seqlen]):
            div, mod = divmod(seqlen, sublen)
            if mod != 0:
                continue
            subseq = seq[:sublen]
            if subseq * div == seq:
                result.append((subseq, div))

        return result


def ref_seq_generator_forward(chrom, pos0, fasta, extend_by):
    chrom_len = fasta.lengths[fasta.references.index(chrom)]
    start0 = pos0
    while True:
        end0 = min(start0 + extend_by, chrom_len)
        seq = fasta.fetch(chrom, start0, end0)
        for x in seq:
            yield x

        if end0 == chrom_len:
            break

        start0 = end0


def ref_seq_generator_backward(chrom, pos0, fasta, extend_by):
    end0 = pos0 + 1
    while True:
        start0 = max(end0 - extend_by, 0)
        seq = fasta.fetch(chrom, start0, end0)
        for x in seq[::-1]:
            yield x

        if start0 == 0:
            break

        end0 = start0


def search_for_repeat(query_seq, chrom, start0, fasta, extend_by=50):
    end0 = start0
    query_len = len(query_seq)

    # forward
    refseq_gen = ref_seq_generator_forward(chrom, start0, fasta, extend_by)
    chunk_refseq_gen = itertools.zip_longest(*([refseq_gen] * query_len), fillvalue=None)
    for ref_seq in chunk_refseq_gen:
        if None in ref_seq:
            # Hitting the start or end of reference sequence
            break
        elif ''.join(ref_seq) == query_seq:
            end0 += query_len
        else:
            break

    # backward
    refseq_gen = ref_seq_generator_backward(chrom, start0 - 1, fasta, extend_by)
    chunk_refseq_gen = itertools.zip_longest(*([refseq_gen] * query_len), fillvalue=None)
    for ref_seq in chunk_refseq_gen:
        if None in ref_seq:
            break
        elif ''.join(ref_seq[::-1]) == query_seq:
            start0 -= query_len
        else:
            break

    # result
    if end0 == start0:
        return None
    else:
        return range(start0, end0)
            
        

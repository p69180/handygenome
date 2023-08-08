import Bio.Seq

def to_pyrimidine(base):
    if base in 'AG':
        return Bio.Seq.reverse_complement(base)
    else:
        return base

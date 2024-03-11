

class CNVEvent:
    def __init__(self, chrom, start0, end0, refver, CN, mean_ploidy, B=None):
        self.chrom = chrom
        self.start0 = start0
        self.end0 = end0
        self.refver = refver
        self.CN = CN
        self.mean_ploidy = mean_ploidy
        self.B = B


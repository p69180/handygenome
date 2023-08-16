import pysam



class CNVBrowser:
    def __init__(self, refver, normal_bam_path=None, tumor_bam_path=None, region_gr=None):
        self.refver = refver
        self.region_gr = region_gr
        self.normal_bam_path = normal_bam_path
        self.tumor_bam_path = tumor_bam_path
        if self.normal_bam_path is None:
            self.normal_bam = pysam.AlignmentFile(self.normal_bam_path)
        if self.tumor_bam_path is None:
            self.tumor_bam = pysam.AlignmentFile(self.tumor_bam_path)

    def __del__(self):
        for bam in (self.normal_bam, self.tumor_bam):
            if bam is not None:
                bam.close()

    def set_gcdata(self):
        pass

import pysam

import handygenome.annotation.data as annotdata


class Gene:
    def __init__(self, name, refver):
        self.name = name
        self.refver = refver

        tabixfile = annotdata.TABIXFILES_GENESET[self.refver]

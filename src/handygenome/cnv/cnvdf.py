
import pandas as pd
import numpy as np
import pyranges as pr

from handygenome.genomedf import GenomeDataFrame



class BAFDataFrame(GenomeDataFrame):
    baf_colname = 'baf'


class CNVSegments(GenomeDataFrame):
    #depth_colname = 'depth'

    ################
    # constructors #
    ################

    def __init__(self, refver, use_median=False):
        super().__init__(refver)
        self.use_median = use_median
        #self.is_perbase = is_perbase

    def init_sanitycheck(self):
        #assert self.__class__.depth_colname in self.annot_cols
        #assert self.__class__.baf_colname in self.annot_cols
        pass

    def spawn(self, frame):
        result = self.__class__(
            self.refver, 
            use_median=self.use_median, 
            #is_perbase=self.is_perbase,
        )
        result.populate_frame(frame)
        return result

    @classmethod
    def from_frame(cls, frame, refver, use_median=False):
        result = cls(refver, use_median=use_median)
        result.populate_frame(frame)
        result.init_sanitycheck()
        return result

    @classmethod
    def read_tsv(cls, filename, refver, use_median=False):
        """Column names must be set appropriately"""
        frame = make_frame_for_read(filename)
        return cls.from_frame(frame, refver, use_median=use_median)





import pandas as pd
import numpy as np
import pyranges as pr

from handygenome.genomedf import GenomeDataFrame


class DepthDataFrame(GenomeDataFrame):
    depth_colname = 'depth'
    normalized_depth_colname = 'normalized_depth'

    def init_sanitycheck(self):
        assert self.__class__.depth_colname in self.annot_cols

    def spawn(self, frame):
        result = self.__class__(
            self.refver, 
            use_median=self.use_median, 
        )
        result.assign_frame(frame)
        return result

    def __init__(self, refver, use_median=False):
        super().__init__(refver)
        self.use_median = use_median

    @classmethod
    def from_frame(cls, frame, refver, use_median=False):
        result = cls(refver, use_median=use_median)
        result.assign_frame(frame)
        result.init_sanitycheck()
        return result

    @classmethod
    def read_tsv(cls, filename, refver, use_median=False):
        """Column names must be set appropriately"""
        frame = cls.make_frame_for_read(filename)
        return cls.from_frame(frame, refver, use_median=use_median)

    @classmethod
    def load_mosdepth(cls, filename, refver, use_median=False):
        assert filename.endswith('.regions.bed.gz')
        frame = cls.make_frame_for_read(filename, annot_cols=[cls.depth_colname])
        return cls.from_frame(frame, refver, use_median=use_median)

    @property
    def depth(self):
        return self[self.__class__.depth_colname].to_numpy()

    #################
    # normalization #
    #################

    @property
    def normalized_depth(self):
        return self[self.__class__.normalized_depth_colname].to_numpy()

    norm_depth = normalized_depth

    def set_normalized_depth(self):
        all_depths = self.depth
        valid_idxs = np.isnan(all_depths)
        valid_depths = all_depths[valid_idxs]
        valid_lenghts = self.lengths[valid_idxs]
        average_depth = np.average(valid_depths, weights=valid_lenghts)

        self[self.__class__.normalized_depth_colname] = all_depths / average_depth




import multiprocessing
import functools

import pandas as pd
import numpy as np
import pyranges as pr
import pysam

import handygenome.peakutils as peakutils
import handygenome.logutils as logutils
from handygenome.genomedf.genomedf import GenomeDataFrame
import handygenome.cnv.cncall as cncall
from handygenome.cnv.segment import SegmentDataFrame
from handygenome.cnv.depth import DepthSegmentDataFrame
from handygenome.cnv.baf import BAFSegmentDataFrame


CLONAL_CN_COLNAME = 'CN'
CLONAL_B_COLNAME = 'B'
SUBCLONAL_CN_COLNAME = 'subCN'
SUBCLONAL_B_COLNAME = 'subB'
CCF_COLNAME = 'ccf'


#class CNVSegmentDataFrame(SegmentDataFrame):
class CNVSegmentDataFrame(DepthSegmentDataFrame, BAFSegmentDataFrame):
    @classmethod
    def from_segments(cls, seg_gdfs, nproc=1):
        assert len(seg_gdfs) > 0

        refvers = set(x.refver for x in seg_gdfs)
        if len(refvers) != 1:
            raise Exception(f'refver differs between segment gdfs: {refvers}')
        refver = refvers.pop()

        merged_seg_gdf = functools.reduce(
            lambda seg1, seg2: (
                seg1.drop_annots('N', inplace=False).isec_union(
                    seg2.drop_annots('N', inplace=False), 
                    nproc=nproc,
                    drop_annots=False,
                    sort=True,
                )
            ),
            seg_gdfs,
        )

        return cls.from_frame(merged_seg_gdf.df, refver=refver)

    # assign

    def assign_clonal_CN(self, data):
        self[CLONAL_CN_COLNAME] = data

    def assign_clonal_B(self, data):
        self[CLONAL_B_COLNAME] = data

    def assign_subclonal_CN(self, data):
        self[SUBCLONAL_CN_COLNAME] = data

    def assign_subclonal_B(self, data):
        self[SUBCLONAL_B_COLNAME] = data

    def assign_ccf(self, data):
        self[CCF_COLNAME] = data

    # fetch

    def get_clonal_CN(self):
        return self[CLONAL_CN_COLNAME]

    def get_clonal_B(self):
        return self[CLONAL_B_COLNAME]

    def get_subclonal_CN(self):
        return self[SUBCLONAL_CN_COLNAME]

    def get_subclonal_B(self):
        return self[SUBCLONAL_B_COLNAME]

    def get_ccf(self):
        return self[CCF_COLNAME]



import multiprocessing

import pandas as pd
import numpy as np
import pyranges as pr
import pysam

import handygenome.peakutils as peakutils
import handygenome.logutils as logutils
from handygenome.genomedf.genomedf import GenomeDataFrame
from handygenome.cnv.segment import SegmentDataFrame 
import handygenome.cnv.cncall as cncall
from handygenome.cnv.cncall import DefaultCNgUnavailableError


class DepthRawDataFrame(GenomeDataFrame):
    depth_colname = 'depth'
    norm_depth_colname = 'normalized_depth'
    MQ_colname = 'MQ'
    DEFAULT_DTYPES = {
        'Chromosome': 'string', 
        'Start': int, 
        'End': int,
        depth_colname: float,
        norm_depth_colname: float,
    }

    def init_sanitycheck(self):
        pass
        #assert self.__class__.depth_colname in self.annot_cols

    def spawn(self, frame):
        return super().spawn(frame, use_median=self.use_median)

    def __init__(self, refver, use_median=False):
        super().__init__(refver)
        self.use_median = use_median

    @classmethod
    def from_frame(cls, frame, refver, dtype=dict(), use_median=False):
        result = cls(refver, use_median=use_median)
        result.assign_frame(frame, dtype)
        result.init_sanitycheck()
        return result

    @classmethod
    def read_tsv(cls, filename, refver, use_median=False):
        """Column names must be set appropriately"""
        frame = cls.make_frame_from_tsv(filename)
        return cls.from_frame(frame, refver, use_median=use_median)

    @classmethod
    def load_mosdepth(cls, filename, refver, use_median=False):
        assert filename.endswith('.regions.bed.gz')
        frame = cls.make_frame_from_tsv(filename, annot_cols=[cls.depth_colname])
        return cls.from_frame(frame, refver, use_median=use_median)

    #################
    # normalization #
    #################

    def set_normalized_depth(self):
        all_depths = self.depth
        valid_idxs = ~np.isnan(all_depths)
        valid_depths = all_depths[valid_idxs]
        valid_lenghts = self.lengths[valid_idxs]
        average_depth = np.average(valid_depths, weights=valid_lenghts)

        self[self.__class__.norm_depth_colname] = all_depths / average_depth

    ################
    # segmentation #
    ################

    def get_segment(self, *args, **kwargs):
        seg = super().get_segment(*args, **kwargs)
        return DepthSegmentDataFrame.from_frame(seg.df, refver=self.refver)

    ##############
    # properties #
    ##############

    @property
    def depth(self):
        return self[self.__class__.depth_colname]

    @property
    def normalized_depth(self):
        return self[self.__class__.norm_depth_colname]

    norm_depth = normalized_depth

    @property
    def mq(self):
        return self[self.__class__.MQ_colname]

    ###############
    # MQ addition #
    ###############

    def get_MQ_fetchregion_gdf(self, readlength=151):
        self.sort()

        new_chroms = list()
        new_start0s = list()
        new_end0s = list()
        halflen = int(readlength / 2)

        for chrom, subgdf in self.group_bychrom(sort=False).items():
            start0s = subgdf.start0s
            end0s = subgdf.end0s
            lengths = subgdf.lengths

            long_length_selector = (lengths > readlength)
            start0s[long_length_selector] = start0s[long_length_selector] + halflen
            end0s[long_length_selector] = end0s[long_length_selector] - halflen

            new_start0s.append(start0s)
            new_end0s.append(end0s)
            new_chroms.append(np.repeat(chrom, len(start0s)))

        result_chroms = np.concatenate(new_chroms)
        result_start0s = np.concatenate(new_start0s)
        result_end0s = np.concatenate(new_end0s)

        return self.spawn_from_data(
            chroms=result_chroms, 
            start0s=result_start0s, 
            end0s=result_end0s,
        )

    def add_MQ(self, bam_path, readlength=151, nproc=1, n_split=1000, verbose=False):
        # get MQ values
        fetchregion_gdf = self.get_MQ_fetchregion_gdf(readlength)  # self is sorted in this step
        split_fetchregions = fetchregion_gdf.equal_nrow_split(n=n_split)
        args = (
            (bam_path, gdf.chroms, gdf.start0s, gdf.end0s, verbose)
            for gdf in split_fetchregions
        )
        with multiprocessing.Pool(nproc) as pool:
            map_result = pool.starmap(make_MQ_array, args)

        # apply results
        self[self.__class__.MQ_colname] = np.concatenate(map_result)

    def check_has_MQ(self):
        return self.__class__.MQ_colname in self.columns


def make_MQ_array(bam_path, chroms, start0s, end0s, verbose=False):
    if verbose:
        logutils.log(f'Beginning region chrom={chroms[0]}, start0={start0s[0]}, end0={end0s[0]}')
    with pysam.AlignmentFile(bam_path) as bam:
        def iterator():
            for c, s, e in zip(chroms, start0s, end0s):
                MQs = tuple(read.mapping_quality for read in bam.fetch(c, s, e))
                if len(MQs) == 0:
                    yield np.nan
                else:
                    yield np.mean(MQs)

        return np.fromiter(iterator(), dtype=float)


class DepthSegmentDataFrame(SegmentDataFrame):
    norm_depth_mean_colname = DepthRawDataFrame.norm_depth_colname + '_mean'
    norm_depth_std_colname = DepthRawDataFrame.norm_depth_colname + '_std'
    MQ_mean_colname = DepthRawDataFrame.MQ_colname + '_mean'

    @property
    def norm_depth_mean(self):
        return self[self.__class__.norm_depth_mean_colname]

    @property
    def norm_depth_std(self):
        return self[self.__class__.norm_depth_std_colname]

    def find_onecopy_depth(self, is_female, ploidy):
        try:
            default_CNg = self.get_default_CNg(is_female)
        except DefaultCNgUnavailableError:
            data = self.norm_depth_mean
            weights = self.lengths
        else:
            valid_selector = np.logical_and(
                (default_CNg != 0),
                ~np.isnan(default_CNg),
            )

            correction_factor = ploidy / default_CNg[valid_selector]
            data = self.norm_depth_mean[valid_selector] * correction_factor

            weights = self.lengths[valid_selector]

        # select copy-neutral depth segments
        peakresult = peakutils.find_density_peaks(
            data,
            weights=weights,
            xs=np.arange(0, 2, 0.01),
            bw_method=0.005,
        )
        modal_depth = peakresult['peak_xs'][np.argmax(peakresult['peak_ys'])]
        onecopy_depth = modal_depth / ploidy

        return onecopy_depth

    def add_rawdata_info(
        self, 
        depth_rawdata_gdf, 
        merge_methods=['mean', 'std'],
        rawdepth=False,
        nproc=1,
    ):
        depth_colname = (
            DepthRawDataFrame.depth_colname
            if rawdepth else
            DepthRawDataFrame.norm_depth_colname
        )
        right_gdf_cols = [depth_colname]

        MQ_colname = DepthRawDataFrame.MQ_colname
        if MQ_colname in depth_rawdata_gdf.columns:
            right_gdf_cols.append(MQ_colname)

        joined_gdf = self.drop_annots().join(
            depth_rawdata_gdf,
            right_gdf_cols=right_gdf_cols,
            how='left',
            merge=merge_methods,
            winsorize=(0.05, 0.05),
            nproc=nproc,
        )
        self.assign_frame(joined_gdf.df)

    @staticmethod
    def get_filter_colname():
        return f'depth_selected'

    def add_filter(self, MQ_cutoff=(50, None)):
        selector = self.filter_helper(self[self.__class__.MQ_mean_colname], MQ_cutoff)
        self[self.get_filter_colname()] = selector

    def get_filter(self):
        return self[self.get_filter_colname()]



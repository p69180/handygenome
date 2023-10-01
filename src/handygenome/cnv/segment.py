import functools
import multiprocessing

import pandas as pd
import numpy as np
import pyranges as pr

#from handygenome.refgenome.refgenome import NoXYError
from handygenome.genomedf.genomedf import GenomeDataFrame
import handygenome.cnv.cncall as cncall



class SegmentDataFrame(GenomeDataFrame):
    @staticmethod
    def filter_helper(values, cutoff):
        if cutoff is None:
            selector = True
        else:
            assert len(cutoff) == 2
            if cutoff[1] is None:
                selector1 = True
            else:
                selector1 = (values < cutoff[1])

            if cutoff[0] is None:
                selector2 = True
            else:
                selector2 = (values > cutoff[0])

            selector = np.logical_and(selector1, selector2)

        return selector

    #################################
    # merging of low ndata segments #
    #################################

    @staticmethod
    def merge_low_ndata_segments_targetfunc(subgdf, ndata_colname, cutoff=30):
        high_ndata_selector = (subgdf[ndata_colname] >= cutoff)
        high_ndata_selector_numeric = high_ndata_selector.astype(int)

        if high_ndata_selector.all():
            new_start0s = subgdf.start0s
            new_end0s = subgdf.end0s
            new_chroms = subgdf.chroms
        else:
            # make grouper
            cumsum = np.cumsum(high_ndata_selector_numeric)
            diff = np.diff(np.insert(high_ndata_selector_numeric, 0, 1))
            low_ndata_zone_starts = np.nonzero(diff == -1)[0]

            N = len(high_ndata_selector)
            offsets = functools.reduce(
                np.add, 
                [np.repeat([0, 1], [x, N - x]) for x in low_ndata_zone_starts]
            )
            grouper = cumsum + offsets

            # do groupby and filtering
            start_groupby = subgdf.df['Start'].groupby(grouper, sort=False)
            end_groupby = subgdf.df['End'].groupby(grouper, sort=False)

            new_start0s = start_groupby.first().to_numpy()
            new_end0s = end_groupby.last().to_numpy()
            new_chroms = np.repeat(subgdf.chroms[0], len(new_end0s))

        return new_chroms, new_start0s, new_end0s

    def merge_low_ndata_segments(self, ndata_colname, cutoff=30, nproc=1):
        assert ndata_colname in self.columns

        bychrom = self.group_bychrom(sort=True)
        args = ((subgdf, ndata_colname, cutoff) for subgdf in bychrom.values())
        with multiprocessing.Pool(nproc) as pool:
            mp_result = pool.starmap(
                self.__class__.merge_low_ndata_segments_targetfunc, 
                args,
            )

        result_chroms, result_start0s, result_end0s = zip(*mp_result)
        result_chroms = np.concatenate(result_chroms)
        result_start0s = np.concatenate(result_start0s)
        result_end0s = np.concatenate(result_end0s)

        result = self.spawn_from_data(
            chroms=result_chroms,
            start0s=result_start0s,
            end0s=result_end0s,
        )
        return result

    @staticmethod
    def incoporate_low_ndata_segments_targetfunc(subgdf, ndata_colname, chrom_length, cutoff=30):
        start0s = subgdf.start0s
        end0s = subgdf.end0s
        ndata_array = subgdf[ndata_colname]
        high_ndata_selector = (ndata_array >= cutoff)
        low_ndata_selector = np.logical_not(high_ndata_selector)

        # 1. remove low ndata segments at borders
        if low_ndata_selector[0]:  # the first segment is with low ndata
            start0s = np.delete(start0s, 0)
            end0s = np.delete(end0s, 0)
            ndata_array = np.delete(ndata_array, 0)
            high_ndata_selector = np.delete(high_ndata_selector, 0)
            low_ndata_selector = np.delete(low_ndata_selector, 0)

            start0s[0] = 0

        if low_ndata_selector[-1]:  # the last segment is with low ndata
            start0s = np.delete(start0s, -1)
            end0s = np.delete(end0s, -1)
            ndata_array = np.delete(ndata_array, -1)
            high_ndata_selector = np.delete(high_ndata_selector, -1)
            low_ndata_selector = np.delete(low_ndata_selector, -1)

            #end0s[-1] = self.chromdict[subgdf.chroms[0]]
            end0s[-1] = chrom_length

        # 2. handle intervening low ndata segment
        if low_ndata_selector.any():
            low_ndata_indexes = np.nonzero(low_ndata_selector)[0]
            low_ndata_midpoints = np.floor(
                (start0s[low_ndata_indexes] + end0s[low_ndata_indexes] - 1) / 2
            ).astype(int)

            before_indexes = low_ndata_indexes - 1
            end0s[before_indexes] = low_ndata_midpoints

            after_indexes = low_ndata_indexes + 1
            start0s[after_indexes] = low_ndata_midpoints

        # 3. result
        new_start0s = start0s[high_ndata_selector]
        new_end0s = end0s[high_ndata_selector]
        new_chroms = np.repeat(subgdf.chroms[0], len(new_start0s))

        return new_chroms, new_start0s, new_end0s

    def incoporate_low_ndata_segments(self, ndata_colname, cutoff=30, nproc=1):
        assert ndata_colname in self.columns

        bychrom = self.group_bychrom(sort=True)
        subgdf_list = list(bychrom.values())
        chrom_length_list = [self.chromdict[x.chroms[0]] for x in subgdf_list]
        args = (
            (subgdf, ndata_colname, chrom_length, cutoff)
            for (subgdf, chrom_length) in zip(subgdf_list, chrom_length_list)
        )
        with multiprocessing.Pool(nproc) as pool:
            mp_result = pool.starmap(
                self.__class__.incoporate_low_ndata_segments_targetfunc, 
                args,
            )

        result_chroms, result_start0s, result_end0s = zip(*mp_result)
        result_chroms = np.concatenate(result_chroms)
        result_start0s = np.concatenate(result_start0s)
        result_end0s = np.concatenate(result_end0s)

        result = self.spawn_from_data(
            chroms=result_chroms,
            start0s=result_start0s,
            end0s=result_end0s,
        )
        return result

    ###############
    # default CNg #
    ###############

    def add_default_CNg_Bg(self, is_female):
        if cncall.DEFAULT_CNG_COLNAME not in self.columns:
            CNg_Bg_gdf = cncall.get_default_CNg_Bg_diploid(self.refver, is_female)
            annot_seg = self.join(
                CNg_Bg_gdf, 
                how='left',
                right_gdf_cols=[cncall.DEFAULT_CNG_COLNAME, cncall.DEFAULT_BG_COLNAME],
                merge='longest',
                overlapping_length=True,
                omit_N=True,
                suffixes={'longest': ''},
            )
            self.assign_df(annot_seg.df)

    def get_default_CNg(self, is_female):
        if cncall.DEFAULT_CNG_COLNAME not in self.columns:
            self.add_default_CNg_Bg(is_female)
        return self[cncall.DEFAULT_CNG_COLNAME]

    def get_default_Bg(self, is_female):
        if cncall.DEFAULT_BG_COLNAME not in self.columns:
            self.add_default_CNg_Bg(is_female)
        return self[cncall.DEFAULT_BG_COLNAME]



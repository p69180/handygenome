import re
import os
import datetime
import inspect
import itertools
import contextlib
import multiprocessing
import functools
import warnings

import pandas as pd
import numpy as np
import pyranges as pr

import handygenome.deco as deco
import handygenome.refgenome.refgenome as refgenome
import handygenome.tools as tools
import handygenome.logutils as logutils
import handygenome.cnv.rdnacopy as rdnacopy

import handygenome.genomedf.genomedf_utils as genomedf_utils
import handygenome.genomedf.genomedf_methods as genomedf_methods

from handygenome.genomedf.genomedf_draw import GenomeDataFrameDrawingBase
from handygenome.genomedf.genomedf_base import GenomeDataFrameBase


class GenomeDataFrame(GenomeDataFrameBase, GenomeDataFrameDrawingBase):

    ################
    # segmentation #
    ################

    def make_compact(self):
        compact_pos0_list, decoder = genomedf_utils.compact_coords(self.start0s, self.end0s)
        compact_gdf = self.copy()
        compact_gdf['Start'] = compact_pos0_list
        compact_gdf['End'] = compact_gdf['Start'] + 1

        return compact_gdf, decoder

    @classmethod
    def get_segment_base(
        cls,
        gdf, 
        annot_colname,
        drop_annots=False,
        N_colname='N',
        smoothing=False, 
        smooth_kwargs=dict(), 
        segment_kwargs=dict(),
    ):
        # set kwargs
        segment_kwargs = (
            dict()
            | segment_kwargs
        )

        # make compact and run segmentation
        compact_gdf, decoder = gdf.make_compact()
        compact_seg_df = rdnacopy.run_dnacopy_segment(
            chromosomes=compact_gdf.chroms,
            positions=compact_gdf.start0s,
            values=compact_gdf[annot_colname],

            N_colname='N',
            value_colname='value',
            smoothing=smoothing, 
            verbose=False, 
            smooth_kwargs=smooth_kwargs, 
            segment_kwargs=segment_kwargs,
        )
        #compact_seg_gdf = gdf.spawn(compact_seg_df)
        compact_seg_gdf = SegmentDataFrame.from_frame(compact_seg_df, refver=gdf.refver)

        # recover noncompact coordinates
        noncompact_start0s, noncompact_end0s = genomedf_utils.uncompact_coords(
            compact_seg_gdf.start0s, 
            compact_seg_gdf.end0s, 
            decoder,
        )

        # make result
        result = SegmentDataFrame.from_data(
            refver=gdf.refver,
            **{
                'chroms': compact_seg_gdf.chroms,
                'start0s': noncompact_start0s,
                'end0s': noncompact_end0s,
                annot_colname: compact_seg_gdf['value'],
                N_colname: compact_seg_gdf['N'],
            }
        )
        if drop_annots:
            result = result.drop_annots()

        return result

    def get_segment(
        self, 
        annot_colname,
        nproc=1,
        split_width=1000,

        drop_annots=False,
        N_colname='N',
        smoothing=False, 
        verbose=False, 
        smooth_kwargs=dict(), 
        segment_kwargs=dict(),
    ):
        assert not self.is_empty

        if verbose:
            logutils.log(f'Beginning segmentation')

        # sort
        self.sort()

        # subset gdf into non-nan region
        notnan_selector = np.logical_not(np.isnan(self[annot_colname]))
        notnan_self = self.loc[notnan_selector, :]

        # make split inputs
        grouped_gdfs = notnan_self.equal_nrow_split_keepchrom(width=split_width, sort=False)
        with multiprocessing.Pool(nproc) as pool:
            args = (
                (
                    gdf, 
                    annot_colname,
                    drop_annots,
                    N_colname,
                    smoothing,
                    smooth_kwargs,
                    segment_kwargs,
                )
                for gdf in grouped_gdfs
            )
            pool_result = pool.starmap(GenomeDataFrame.get_segment_base, args)

        result = SegmentDataFrame.concat(pool_result)
        result.sort()

        if verbose:
            logutils.log(f'Finished segmentation')

        return result


class SegmentDataFrame(GenomeDataFrame):

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




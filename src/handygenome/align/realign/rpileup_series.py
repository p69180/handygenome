import sys
import collections
import itertools
import functools
#import logging
import inspect
import random
import uuid

import pysam
import Bio.Align
import numpy as np
import pandas as pd
#import pyranges as pr

from handygenome.genomedf.genomedf_base import GenomeDataFrameBase
import handygenome.tools as tools
import handygenome.refgenome.refgenome as refgenome
from handygenome.refgenome.refgenome import RefverObjectBase

import handygenome.workflow as workflow
import handygenome.variant.vcfspec as libvcfspec
import handygenome.read.pileup as libpileup
import handygenome.align.alignhandler as alignhandler
import handygenome.bameditor as bameditor
import handygenome.read.readhandler as readhandler
import handygenome.read.readplus as readplus
import handygenome.align.realign.base as realign_base
from handygenome.align.realign.base import LoggingBase
from handygenome.align.realign.rpileup import RealignerPileup


class RealignerPileupSeries(RefverObjectBase, LoggingBase):
    def __init__(
        self, 
        bam, 
        chrom, 
        start0, 
        end0, 
        _refver=None,
        aligner=realign_base.DEFAULT_ALIGNER,
        verbose=False,
        **kwargs,
    ):
        """Args:
            **kwargs: keyword arguments of RealignerPileup.__init__
                except start0, end0, init_df
        """
        # set params
        self.chrom = chrom
        self.bam = bam
        if _refver is None:
            self.refver = refgenome.infer_refver_bamheader(bam.header)
        else:
            self.refver = _refver
        self.aligner = aligner
        self.verbose = verbose
        self.params = realign_base.parse_rpileup_kwargs(**kwargs)

        # setup pileup_list
        self._init_pileup_list(start0, end0)
        # merge individual read_store's
        self.read_store = dict()
        for pileup in self.pileup_list:
            self.read_store.update(pileup.read_store)

    def _init_pileup_list(self, start0, end0):
        # initiate pileup_list
        self.log_debug(f'Initial pileup generation')

        initial_pileup, result_inactive, result_vcfspec = RealignerPileup.init_and_augment(
            bam=self.bam, 
            chrom=self.chrom, start0=start0, end0=end0, 
            aligner=self.aligner,
            verbose=self.verbose, 
            **self.params,
        )
        self.pileup_list = [initial_pileup]

        #self.log_debug(f'Initial pileup: {initial_pileup}')
        #self.log_debug(f'result_inactive: {result_inactive}')
        #self.log_debug(f'result_vcfspec: {result_vcfspec}')
        #self.log_debug(f'Initial pileup generation Finished\n')

        if (result_inactive is None) and (result_vcfspec is None):
            # When there is no active region in the initial pileup
            return

        #left_insufficient = (not result_inactive.left_okay) or (not result_vcfspec.left_okay)
        #right_insufficient = (not result_inactive.right_okay) or (not result_vcfspec.right_okay)
        left_insufficient = (not result_inactive.left_okay)
        right_insufficient = (not result_inactive.right_okay)

        if left_insufficient:
            self.secure_left()
        if right_insufficient:
            self.secure_right()

        #self.logger.debug(f'Finished initial pileup generation\n')

    def __repr__(self):
        string = '\n'.join(
            f'\t{idx} {pileup}'
            for idx, pileup in enumerate(self.pileup_list)
        )
        return (
            f'<'
            f'{self.__class__.__name__}'
            f'(pileup_list: [\n{string}\n])'
            f'>'
        )

    def __len__(self):
        return len(self.pileup_list)

    @property
    def start0(self):
        return self.pileup_list[0].start0

    @property
    def end0(self):
        return self.pileup_list[-1].end0

    @property
    def range0(self):
        return range(self.start0, self.end0)

    @property
    def width(self):
        return self.pileup_list[-1].end0 - self.pileup_list[0].start0

    #def get_gr(self):
    def get_gdf(self):
        start0s = list()
        end0s = list()
        for pileup in self.pileup_list:
            start0s.append(pileup.start0)
            end0s.append(pileup.end0)
        chroms = [self.chrom] * len(self.pileup_list)
        self_indexes = list(range(len(self.pileup_list)))
        #return pr.from_dict(
        #    {'Chromosome': chroms, 'Start': start0s, 'End': end0s, 'Self_index': self_indexes}
        #)
        return GenomeDataFrameBase.from_data(
            refver=self.refver, chroms=self.chrom, start0s=start0s, end0s=end0s, Self_index=self_indexes,
        )

    def spawn_new_pileup(self, start0, end0):
        return RealignerPileup(
            bam=self.bam, 
            chrom=self.chrom, 
            start0=start0, 
            end0=end0, 
            _refver=self.refver, 
            aligner=self.aligner,
            init_df=True,
            verbose=self.verbose, 
            #logger=self.logger,
            **self.params,
        )

    def left_insert_new_pileup(self, width=1):
        new_pileup = self.spawn_new_pileup(
            start0=(self.start0 - width), end0=self.start0,
        )
        self.pileup_list.insert(0, new_pileup)

    def right_insert_new_pileup(self, width=1):
        new_pileup = self.spawn_new_pileup(
            start0=self.end0, end0=(self.end0 + width),
        )
        self.pileup_list.append(new_pileup)

    def extend_left(self, width):
        """Ignores MQ_limit, depth_limit, start0_limit, end0_limit, max_series_width attributes
        """
        assert width > 0

        target_start0 = self.start0 - width

        original_params = self.params.copy()
        self.params['max_series_width'] = np.inf
        self.params['start0_limit'] = target_start0
        self.params['MQ_limit'] = 0
        self.params['depth_limit'] = 0

        while True:
            current_width = min(self.params['max_pileup_width'], self.start0 - target_start0)
            self.left_insert_new_pileup(width=current_width)
            if self.start0 == target_start0:
                break

        self.params = original_params

#    def extend_left_old(self, width):
#        """Ignores MQ_limit, depth_limit, start0_limit, end0_limit,
#        max_series_width attribute"""
#        original_max_series_width = self.max_series_width
#        self.max_series_width = np.inf
#        original_start0_limit = self.start0_limit
#        self.start0_limit = -np.inf
#
#        target_start0 = self.start0 - width
#
#        rpileup_init_kwargs = self.rpileup_init_kwargs.copy()
#        rpileup_init_kwargs['start0_limit'] = target_start0
#        rpileup_init_kwargs['MQ_limit'] = 0
#        rpileup_init_kwargs['depth_limit'] = 0
#
#        while True:
#            new_pileup = RealignerPileup(
#                bam=self.bam, 
#                chrom=self.chrom, start0=(self.start0 - 1), end0=self.start0, 
#                refver=self.refver, fasta=self.fasta, 
#                verbose=self.verbose, logger=self.logger,
#                **rpileup_init_kwargs,
#            )
#            self.pileup_list.insert(0, new_pileup)
#            result_inactive, result_vcfspec, touched_width_limit = self.secure_left()
#            if self.start0 <= target_start0:
#                break
#
#        self.max_series_width = original_max_series_width
#        self.start0_limit = original_start0_limit

    def extend_right(self, width):
        """Ignores MQ_limit, depth_limit, start0_limit, end0_limit, max_series_width attributes
        """
        assert width > 0

        target_end0 = self.end0 + width

        original_params = self.params.copy()
        self.params['max_series_width'] = np.inf
        self.params['end0_limit'] = target_end0
        self.params['MQ_limit'] = 0
        self.params['depth_limit'] = 0

        while True:
            current_width = min(self.params['max_pileup_width'], target_end0 - self.end0)
            self.right_insert_new_pileup(width=current_width)
            if self.end0 == target_end0:
                break

        self.params = original_params

#    def extend_right_old(self, width):
#        """Ignores MQ_limit, depth_limit, start0_limit, end0_limit,
#        max_series_width attribute
#        """
#        original_max_series_width = self.max_series_width
#        self.max_series_width = np.inf
#        original_end0_limit = self.end0_limit
#        self.end0_limit = np.inf
#
#        target_end0 = self.end0 + width
#
#        rpileup_init_kwargs = self.rpileup_init_kwargs.copy()
#        rpileup_init_kwargs['end0_limit'] = target_end0
#        rpileup_init_kwargs['MQ_limit'] = 0
#        rpileup_init_kwargs['depth_limit'] = 0
#
#        while True:
#            new_pileup = RealignerPileup(
#                bam=self.bam, 
#                chrom=self.chrom, start0=self.end0, end0=(self.end0 + 1), 
#                refver=self.refver, fasta=self.fasta, 
#                verbose=self.verbose, logger=self.logger,
#                **rpileup_init_kwargs,
#            )
#            self.pileup_list.append(new_pileup)
#            result_inactive, result_vcfspec, touched_width_limit = self.secure_right()
#            if self.end0 >= target_end0:
#                break
#
#        self.max_series_width = original_max_series_width
#        self.end0_limit = original_end0_limit

    # init helpers #
    def secure_left(self, inactive_only=True):
        self.log_debug(f'Beginning RealignerPileupSeries secure_left')

        left_pileup = self.pileup_list[0]
        left_pileup.params = self.params.copy()
        left_pileup.params['end0_limit'] = left_pileup.end0
        while True:
            result_inactive, result_vcfspec = left_pileup.secure_margins(inactive_only=inactive_only)

            self.log_debug(self.pileup_list)
            self.log_debug(f'self.start0: {self.start0}')
            self.log_debug(f'result_inactive: {result_inactive}')
            self.log_debug(f'result_vcfspec: {result_vcfspec}')

            if inactive_only:
                okay_cond = (
                    result_inactive.left_okay
                    or (self.width > self.params['max_series_width'])
                )
            else:
                okay_cond = (
                    (result_inactive.left_okay and result_vcfspec.left_okay)
                    or (
                        result_vcfspec.left_low_depth or 
                        result_vcfspec.left_low_MQ or 
                        result_vcfspec.touched_left_limit
                    )
                    or (self.width > self.params['max_series_width'])
                )

            if okay_cond:
                break
            else:
                if not inactive_only:
                    assert result_vcfspec.touched_width_limit

                left_pileup.prepare_vcfspecs()
                split_points = left_pileup.get_split_points()
                if len(split_points) == 0:
                    self.left_insert_new_pileup(width=1)
                else:
                    split_pileups = left_pileup.split(split_points[:1])
                    del self.pileup_list[0]
                    self.pileup_list = split_pileups + self.pileup_list

                left_pileup = self.pileup_list[0]
                left_pileup.params['end0_limit'] = left_pileup.end0

        touched_width_limit = self.width >= self.params['max_series_width']

        self.log_debug(f'Finished RealignerPileupSeries secure_left')

        return result_inactive, result_vcfspec, touched_width_limit

    def secure_right(self, inactive_only=True):
        #self.logger.debug(f'Beginning RealignerPileupSeries secure_right')

        right_pileup = self.pileup_list[-1]
        right_pileup.params = self.params.copy()
        right_pileup.params['start0_limit'] = right_pileup.start0
        while True:
            result_inactive, result_vcfspec = right_pileup.secure_margins(inactive_only=inactive_only)

            #self.logger.debug(self.pileup_list)
            #self.logger.debug(f'self.end0: {self.end0}')
            #self.logger.debug(f'result_inactive: {result_inactive}')
            #self.logger.debug(f'result_vcfspec: {result_vcfspec}')

            if inactive_only:
                okay_cond = (
                    result_inactive.right_okay
                    or (self.width > self.params['max_series_width'])
                )
            else:
                okay_cond = (
                    (result_inactive.right_okay and result_vcfspec.right_okay)
                    or (
                        result_vcfspec.right_low_depth or 
                        result_vcfspec.right_low_MQ or
                        result_vcfspec.touched_right_limit
                    ) 
                    or (self.width > self.params['max_series_width'])
                )

            if okay_cond:
                break
            else:
                if not inactive_only:
                    assert result_vcfspec.touched_width_limit

                right_pileup.prepare_vcfspecs()
                split_points = right_pileup.get_split_points()
                if len(split_points) == 0:
                    self.right_insert_new_pileup(width=1)
                else:
                    split_pileups = right_pileup.split(split_points[-1:])
                    del self.pileup_list[-1]
                    self.pileup_list.extend(split_pileups)

                right_pileup = self.pileup_list[-1]
                right_pileup.params['start0_limit'] = right_pileup.start0

        touched_width_limit = self.width >= self.params['max_series_width']

        #self.logger.debug(f'Finished RealignerPileupSeries secure_right')

        return result_inactive, result_vcfspec, touched_width_limit

#    def secure_left_old(self):
#        self.logger.debug(f'Beginning RealignerPileupSeries secure_left')
#
#        left_pileup = self.pileup_list[0]
#        left_pileup.end0_limit = left_pileup.end0
#        while True:
#            result_inactive, result_vcfspec = left_pileup.secure_margins()
#
#            self.logger.debug(self.pileup_list)
#            self.logger.debug(f'self.start0: {self.start0}')
#            self.logger.debug(f'result_inactive: {result_inactive}')
#            self.logger.debug(f'result_vcfspec: {result_vcfspec}')
#
#            if (
#                (result_inactive.left_okay and result_vcfspec.left_okay)
#                or (
#                    result_vcfspec.left_low_depth or 
#                    result_vcfspec.left_low_MQ or 
#                    result_vcfspec.touched_left_limit
#                )
#                or (self.width > self.max_series_width)
#            ):
#                break
#            else:
#                assert result_vcfspec.touched_width_limit
#                split_points = left_pileup.get_split_points()
#                if len(split_points) == 0:
#                    rpileup_init_kwargs = self.rpileup_init_kwargs.copy()
#                    new_pileup = RealignerPileup(
#                        bam=self.bam, 
#                        chrom=self.chrom, start0=(left_pileup.start0 - 1), end0=left_pileup.start0, 
#                        refver=self.refver, fasta=self.fasta, 
#                        verbose=self.verbose, logger=self.logger,
#                        **rpileup_init_kwargs,
#                    )
#                    self.pileup_list.insert(0, new_pileup)
#                else:
#                    split_pileups = left_pileup.split(split_points[:1])
#                    del self.pileup_list[0]
#                    self.pileup_list = split_pileups + self.pileup_list
#
#                left_pileup = self.pileup_list[0]
#                left_pileup.end0_limit = left_pileup.end0
#
#        touched_width_limit = self.width >= self.max_series_width
#
#        self.logger.debug(f'Finished RealignerPileupSeries secure_left')
#
#        return result_inactive, result_vcfspec, touched_width_limit
#
#    def secure_right_old(self):
#        self.logger.debug(f'Beginning RealignerPileupSeries secure_right')
#
#        right_pileup = self.pileup_list[-1]
#        right_pileup.start0_limit = right_pileup.start0
#        while True:
#            result_inactive, result_vcfspec = right_pileup.secure_margins()
#
#            self.logger.debug(self.pileup_list)
#            self.logger.debug(f'self.end0: {self.end0}')
#            self.logger.debug(f'result_inactive: {result_inactive}')
#            self.logger.debug(f'result_vcfspec: {result_vcfspec}')
#
#            if (
#                (result_inactive.right_okay and result_vcfspec.right_okay)
#                or (
#                    result_vcfspec.right_low_depth or 
#                    result_vcfspec.right_low_MQ or
#                    result_vcfspec.touched_right_limit
#                ) or
#                or (self.width > self.max_series_width)
#            ):
#                break
#            else:
#                assert result_vcfspec.touched_width_limit
#                split_points = right_pileup.get_split_points()
#                if len(split_points) == 0:
#                    rpileup_init_kwargs = self.rpileup_init_kwargs.copy()
#                    new_pileup = RealignerPileup(
#                        bam=self.bam, 
#                        chrom=self.chrom, start0=right_pileup.end0, end0=(right_pileup.end0 + 1), 
#                        refver=self.refver, fasta=self.fasta, 
#                        verbose=self.verbose, logger=self.logger,
#                        **rpileup_init_kwargs,
#                    )
#                    self.pileup_list.append(new_pileup)
#                else:
#                    split_pileups = right_pileup.split(split_points[-1:])
#                    del self.pileup_list[-1]
#                    self.pileup_list.extend(split_pileups)
#
#                right_pileup = self.pileup_list[-1]
#                right_pileup.start0_limit = right_pileup.start0
#
#        touched_width_limit = self.width >= self.max_series_width
#
#        self.logger.debug(f'Finished RealignerPileupSeries secure_right')
#
#        return result_inactive, result_vcfspec, touched_width_limit

    def get_splittable_region_best(self):
        # without splitting contig vcfspec
        return GenomeDataFrameBase.concat(
            [
                rpileup.get_vcfspec_margins_gdf(
                    subseq_portion_threshold=(
                        rpileup.params['allele_portion_threshold'] * 2
                    ),
                    inverse=True,
                    split_contig_vcfspec=False,
                )
                for rpileup in self.pileup_list
            ]
        ).merge()

    def get_splittable_region_split_contig(self):
        return GenomeDataFrameBase.concat(
            [
                rpileup.get_vcfspec_margins_gdf(
                    subseq_portion_threshold=(
                        rpileup.params['allele_portion_threshold'] * 2
                    ),
                    inverse=True,
                    split_contig_vcfspec=True,
                )
                for rpileup in self.pileup_list
            ]
        ).merge()

    def get_splittable_region_inactive_runs(self):
        #inactive_runs_gr_list = list()
        inactive_runs_gdf_list = list()
        for rpileup in self.pileup_list:
            #active_info_gr = rpileup.get_active_info_gr()
            #inactive_gr = active_info_gr[~active_info_gr.Active]
            #inactive_runs_gr_list.append(inactive_gr)
            active_info_gdf = rpileup.get_active_info_gdf()
            inactive_gdf = active_info_gdf.loc[~active_info_gdf['Active'], :]
            inactive_runs_gdf_list.append(inactive_gdf)
        result = GenomeDataFrameBase.concat(inactive_runs_gdf_list).merge()
        if result.is_empty:
            raise Exception(f'No inactive runs in this RealignerPileupSeries object.')
        return result

    def prepare_vcfspecs(self):
        for rpileup in self.pileup_list:
            rpileup.prepare_vcfspecs()

    @staticmethod
    def rearrange_merger(left_pileup, right_pileup):
        left_pileup.merge(right_pileup, other_on_left=False)
        return left_pileup

    def rearrange(self, start0_list, prepare_vcfspecs=True):
        """Mutates in-place"""
        #def merger(left_pileup, right_pileup):
        #    left_pileup.merge(right_pileup, other_on_left=False)
        #    return left_pileup
        
        # sanity check
        if not (start0_list[0] >= self.start0 and start0_list[-1] <= self.end0):
            raise Exception(f'Input splitting range is out of PileupSeries object range.')
        
        # make new pileup_list
#        split_ranges_gr = pr.from_dict(
#            {
#                'Chromosome': [self.chrom] * (len(start0_list) - 1),
#                'Start': start0_list[:-1],
#                'End': start0_list[1:],
#                'Ranges_index': list(range(len(start0_list) - 1)),
#            }
#        )
        split_ranges_gdf = GenomeDataFrameBase.from_data(
            refver=self.refver,
            chroms=self.chrom,
            start0s=start0_list[:-1],
            end0s=start0_list[1:],
            Ranges_index=list(range(len(start0_list) - 1)),
        )
        #joined_gr = split_ranges_gr.join(self.get_gr())
        joined_gdf = split_ranges_gdf.join(self.get_gdf(), merge=None, how='left')

        #new_pileup_list = list()
        #for ranges_index, subiter in itertools.groupby(
        #    (x[1] for x in joined_gr.df.iterrows()),
        #    key=(lambda x: x.Ranges_index),
        #):
#            partial_pileups = list()
#            for row in subiter:
#                new_start0 = max(row.Start, row.Start_b)
#                new_end0 = min(row.End, row.End_b)
#                partial_pileups.append(
#                    self.pileup_list[row.Self_index].subset(new_start0, new_end0)
#                )
#            new_pileup_list.append(functools.reduce(merger, partial_pileups))

        new_pileup_list = list()
        for key, subdf in joined_gdf.df.groupby('Ranges_index'):
            partial_pileups = list()
            for self_index, new_start0, new_end0 in zip(
                subdf['Self_index'].astype(int),
                np.maximum(subdf['Start'], subdf['Start_b']),
                np.minimum(subdf['End'], subdf['End_b']),
            ):
                partial_pileups.append(
                    self.pileup_list[self_index].subset(new_start0, new_end0)
                )
            new_pileup_list.append(functools.reduce(self.rearrange_merger, partial_pileups))

        # result
        self.pileup_list = new_pileup_list
        if prepare_vcfspecs:
            self.prepare_vcfspecs()

    @staticmethod
    def subseq_aln_to_cigartuples(alignment, left_is_empty, right_is_empty):
        # get active region cigartuples
        active_region_cigartuples, active_region_offset = alignhandler.alignment_to_cigartuples(
            alignment,
            match_as_78=False, del_as_skip=False, 
            left_ins_as_clip=left_is_empty, right_ins_as_clip=right_is_empty,
            remove_left_del=left_is_empty, remove_right_del=right_is_empty,
        )
        return active_region_cigartuples, active_region_offset

    def reduce_cigartuples_list_decorated(self, cigartuples_list_decorated):
        while True:
            cigartuples_list_decorated, reduced = self.softclip_reduce_cigartuples_list_decorated(cigartuples_list_decorated)
            if not reduced:
                break

        return functools.reduce(
            self.merge_two_cigartuples_plain,
            cigartuples_list_decorated,
        )[1]

    def softclip_reduce_cigartuples_list_decorated(self, cigartuples_list_decorated):
        reduced = False
        for idx, (left_item, right_item) in enumerate(tools.pairwise(cigartuples_list_decorated)):
            ref_range0_left, cigartuples_left, read_seq_left = left_item
            ref_range0_right, cigartuples_right, read_seq_right = right_item

            if (
                (ref_range0_right is not None)
                and (len(cigartuples_left) == 0)
                and (cigartuples_right[0][0] == 4)
                and (alignhandler.get_target_length(cigartuples_right) == len(ref_range0_right))
            ):  # Exceptional case where softclip faces empty side (221221)
                reduced = True
                merged_items = self.merge_two_cigartuples_right_clip_left_empty(
                    left_item, right_item
                )
                break
            elif (
                (ref_range0_left is not None)
                and (len(cigartuples_right) == 0)
                and (cigartuples_left[-1][0] == 4)
                and (alignhandler.get_target_length(cigartuples_left) == len(ref_range0_left))
            ):  # Exceptional case where softclip faces empty side (221221)
                reduced = True
                merged_items = self.merge_two_cigartuples_left_clip_right_empty(
                    left_item, right_item
                )
                break

        if reduced:
            result_cigartuples_list_decorated = (
                cigartuples_list_decorated[:idx]
                + [merged_items]
                + cigartuples_list_decorated[(idx + 2):]
            )
        else:
            result_cigartuples_list_decorated = cigartuples_list_decorated

        return result_cigartuples_list_decorated, reduced

    def merge_two_cigartuples_helper(self, read_seq_left, read_seq_right, ref_range0_left, ref_range0_right):
        new_read_seq = read_seq_left + read_seq_right
        if (ref_range0_left is None) or (ref_range0_right is None):
            new_ref_range0 = None
        else:
            new_ref_range0 = range(ref_range0_left.start, ref_range0_right.stop)

        return new_read_seq, new_ref_range0

    def merge_two_cigartuples_right_clip_left_empty(self, left_item, right_item):
        ref_range0_left, cigartuples_left, read_seq_left = left_item
        ref_range0_right, cigartuples_right, read_seq_right = right_item
        
        ref_seq = self.fasta.fetch(
            self.chrom, 
            ref_range0_right.start - int(cigartuples_right[0][1] + 20),
            ref_range0_right.start,
        )
        alignment = realign_base.main_aligner(
            target=ref_seq,
            query=read_seq_right[:cigartuples_right[0][1]],
            aligner=self.aligner,
            reverse_align=True,
        )
        new_cigartuples, offset = self.subseq_aln_to_cigartuples(
            alignment, left_is_empty=True, right_is_empty=False,
        )

        new_cigartuples = new_cigartuples + cigartuples_right[1:]

        new_read_seq, new_ref_range0 = self.merge_two_cigartuples_helper(
            read_seq_left, read_seq_right, ref_range0_left, ref_range0_right,
        )

        return new_ref_range0, new_cigartuples, new_read_seq

    def merge_two_cigartuples_left_clip_right_empty(self, left_item, right_item):
        ref_range0_left, cigartuples_left, read_seq_left = left_item
        ref_range0_right, cigartuples_right, read_seq_right = right_item
        
        ref_seq = self.fasta.fetch(
            self.chrom, 
            ref_range0_left.stop,
            ref_range0_left.stop + int(cigartuples_left[-1][1] + 20),
        )
        alignment = realign_base.main_aligner(
            target=ref_seq,
            query=read_seq_left[-cigartuples_left[-1][1]:],
            aligner=self.aligner,
            reverse_align=False,
        )

        new_cigartuples, offset = self.subseq_aln_to_cigartuples(
            alignment, left_is_empty=False, right_is_empty=True,
        )

        new_cigartuples = cigartuples_left[:-1] + new_cigartuples

        new_read_seq, new_ref_range0 = self.merge_two_cigartuples_helper(
            read_seq_left, read_seq_right, ref_range0_left, ref_range0_right,
        )

        return new_ref_range0, new_cigartuples, new_read_seq

    def merge_two_cigartuples_plain(self, left_item, right_item):
        ref_range0_left, cigartuples_left, read_seq_left = left_item
        ref_range0_right, cigartuples_right, read_seq_right = right_item
        
        if len(cigartuples_left) == 0 or len(cigartuples_right) == 0:
            new_cigartuples = cigartuples_left + cigartuples_right
            #return cigartuples_left + cigartuples_right
        else:
            if cigartuples_left[-1][0] == cigartuples_right[0][0]:
                new_cigartuples = (
                    cigartuples_left[:-1]
                    + [(cigartuples_left[-1][0], cigartuples_left[-1][1] + cigartuples_right[0][1])]
                    + cigartuples_right[1:]
                )
            else:
                new_cigartuples = cigartuples_left + cigartuples_right

        new_read_seq, new_ref_range0 = self.merge_two_cigartuples_helper(
            read_seq_left, read_seq_right, ref_range0_left, ref_range0_right,
        )

        return new_ref_range0, new_cigartuples, new_read_seq

    def save_subseq_alignments(self):
        for pileup in self.pileup_list:
            pileup.save_subseq_alignments()

    def set_realigned_reads(self):
        self.realigned_reads = dict()
        for read_uid, read in self.read_store.items():
            # select relevant pileups 
            relevant_pileups = list()
            for pileup in self.pileup_list:
                if read_uid in pileup.seq_df.index:
                    relevant_pileups.append(pileup)
            if len(relevant_pileups) == 0:
                # During 'subset' method of rpileup, out-of-range reads remain in 'read_store'
                continue

            # split cigartuples
            before_cigartuples, within_cigartuples, after_cigartuples = alignhandler.split_cigar(
                cigartuples=read.cigartuples, 
                reference_start0=read.reference_start,
                split_range0=range(self.start0, self.end0),
                trailing_queryonly_to_right=False,
            )

            before_empty = not bool(before_cigartuples)
            after_empty = not bool(after_cigartuples)

            # make cigartuples from alignment within each pileup
            offset_list = list()
            cigartuples_list = list()  # before, within1, within2, ..., after
            ref_range0_list = list()  # before, within1, within2, ..., after

            def helper(pileup, left_is_empty, right_is_empty):
                cigartuples, offset = self.subseq_aln_to_cigartuples(
                    alignment=pileup._subseq_alignments[read_uid], 
                    left_is_empty=left_is_empty,
                    right_is_empty=right_is_empty,
                )
                offset_list.append(offset)
                cigartuples_list.append(cigartuples)
                ref_range0_list.append(pileup.range0)

            cigartuples_list.append(before_cigartuples)
            ref_range0_list.append(None)  # ref_range0 is None for before_cigartuples

            if len(relevant_pileups) == 1:
                helper(relevant_pileups[0], before_empty, after_empty)
            else:
                helper(relevant_pileups[0], before_empty, False)
                for pileup in relevant_pileups[1:-1]:
                    helper(pileup, False, False)
                helper(relevant_pileups[-1], False, after_empty)

            cigartuples_list.append(after_cigartuples)
            ref_range0_list.append(None)

            # make partial read sequences
            partial_read_seqs = list()
            end_idx = 0
            for cigartuples in cigartuples_list:
                start_idx = end_idx
                end_idx += alignhandler.get_query_length(cigartuples)
                partial_read_seqs.append(read.query_sequence[start_idx:end_idx])

            cigartuples_list_decorated = list(
                zip(ref_range0_list, cigartuples_list, partial_read_seqs)
            )

            realigned_cigartuples = self.reduce_cigartuples_list_decorated(
                cigartuples_list_decorated
            )

            # get new read start position
            if not before_empty:
                new_reference_start = read.reference_start
            else:
                new_reference_start = relevant_pileups[0].start0 + offset_list[0]

            # make new read object
            realigned_read = read.__copy__()
            realigned_read.reference_start = new_reference_start
            realigned_read.cigartuples = realigned_cigartuples

            #####################################################

            readhandler.set_NMMD(realigned_read, fasta=self.fasta)
            self.realigned_cigar_sanity_check(read, realigned_read, realigned_cigartuples)
            self.realigned_reads[read_uid] = realigned_read

    def write_realigned_reads(self, bam_path, padding=0):
        # make header
        hdr = pysam.AlignmentHeader.from_references(
            reference_names=self.fasta.references,
            reference_lengths=self.fasta.lengths,
        )
        with pysam.AlignmentFile(bam_path, mode='wb', header=hdr) as in_bam:
            # write non-realigned reads
            for read in readhandler.get_fetch(
                self.bam,
                self.chrom,
                self.start0 - padding,
                self.end0 + padding,
            ):
                uid = readhandler.get_uid(read)
                if uid not in self.realigned_reads.keys():
                    in_bam.write(read)
            # write realigned reads
            for read in self.realigned_reads.values():
                in_bam.write(read)
        # sort and index
        bameditor.sort_and_index(bam_path)

    def get_result_vcfspecs(self, as_components=True, MQ_threshold=40): 
        return list(
            itertools.chain.from_iterable(
                rpileup.get_result_vcfspecs(as_components=as_components, MQ_threshold=MQ_threshold)
                for rpileup in self.pileup_list
            )
        )

    @staticmethod
    def realigned_cigar_sanity_check(original_read, realigned_read, realigned_cigartuples):
        if len(original_read.query_sequence) != alignhandler.get_query_length(realigned_cigartuples):
            raise realign_base.RealignmentCigarError(
                f'Read query sequence length and cigar length differs.\n'
                f'Original read:\n{original_read.to_string()}\n'
                f'Realigned read:\n{realigned_read.to_string()}\n'
            )

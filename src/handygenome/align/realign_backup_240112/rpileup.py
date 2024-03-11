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
import handygenome.deco as deco
import handygenome.workflow as workflow
import handygenome.variant.vcfspec as libvcfspec
#import handygenome.read.pileup as libpileup
from handygenome.read.pileup import PileupBase
import handygenome.align.alignhandler as alignhandler
import handygenome.bameditor as bameditor
import handygenome.read.readhandler as readhandler
import handygenome.read.readplus as readplus
import handygenome.align.realign.base as realign_base
from handygenome.align.realign.base import RealignerPileupBase


class RealignerPileup(RealignerPileupBase):
    # initializers
    def __init__(
        self, 
        bam, 
        chrom, 
        start0=None, 
        end0=None, 

        _refver=None, # only for spawn

        aligner=realign_base.DEFAULT_ALIGNER,
        init_df=True,

        verbose=False,
        #logger=None,

        **kwargs,
    ):
        """Args:
            start0 and end0 are only required when init_df == True
        """
        # initiation
        PileupBase.__init__(
            self, 
            bam=bam, 
            chrom=chrom, start0=start0, end0=end0, 
            _refver=_refver,
            init_df=init_df, 
            verbose=verbose, 
            #logger=logger,
        )

        # set other parameters
        self.params = realign_base.parse_rpileup_kwargs(**kwargs)
        self.aligner = aligner
        if init_df:
            self._set_active_info()

        self.hit_left_margin = False
        self.hit_right_margin = False

    @classmethod
    def init_and_augment(
        cls, 
        bam, 
        chrom, 
        start0, 
        end0, 
        aligner=realign_base.DEFAULT_ALIGNER,
        verbose=False, 
        #logger=None, 
        **kwargs,
    ):
        result = cls(
            bam=bam, 
            chrom=chrom, start0=start0, end0=end0, 
            aligner=aligner,

            init_df=True,

            verbose=verbose,
            #logger=logger,

            **kwargs,
        )
        secresult_inactive, secresult_vcfspec = result.augment_margins()

        return result, secresult_inactive, secresult_vcfspec

    def spawn(self):
        result = self.__class__(
            bam=self.bam,
            chrom=self.chrom,
            _refver=self.refver,
            aligner=self.aligner,
            init_df=False,
            verbose=self.verbose,
            #logger=self.logger,
            **self.params,
        )
        result.read_store = self.read_store
        return result
    
    def subset(self, start0, end0, inplace=False):
        new_active_info = self.active_info.loc[start0:(end0 - 1)]
        if inplace:
            PileupBase._subset_base(self, start0, end0, inplace=inplace)
            self.active_info = new_active_info
        else:
            result = PileupBase._subset_base(self, start0, end0, inplace=inplace)
            result.active_info = new_active_info
            return result

    def merge(self, other, other_on_left):
        self._merge_base(other, other_on_left=other_on_left)
        if other_on_left:
            self.active_info = pd.concat([other.active_info, self.active_info])
        else:
            self.active_info = pd.concat([self.active_info, other.active_info])

    def split(self, start0_list):
        split_pileups = self._split_base(start0_list)
        for pileup in split_pileups:
            pileup.params = self.params.copy()
            pileup.prepare_vcfspecs()
        return split_pileups

    ################
    def check_touches_left_limit(self, start0=None):
        if start0 is None:
            start0 = self.start0
        return start0 <= self.params['start0_limit']

    def check_touches_right_limit(self, end0=None):
        if end0 is None:
            end0 = self.end0
        return end0 >= self.params['end0_limit']

    def check_touches_width_limit(self, start0=None, end0=None):
        if start0 is None:
            start0 = self.start0
        if end0 is None:
            end0 = self.end0
        return end0 - start0 >= self.params['max_pileup_width']

    def check_left_MQ_low(self):
        return self.MQ.iloc[0] < self.params['MQ_limit']

    def check_right_MQ_low(self):
        return self.MQ.iloc[-1] < self.params['MQ_limit']

    def check_left_depth_low(self):
        return self.get_depth(self.start0) < self.params['depth_limit']

    def check_right_depth_low(self):
        return self.get_depth(self.end0 - 1) < self.params['depth_limit']

    # extend and its helpers
    def extend_left(self, width):
        #assert self.start0 - width >= self.start0_limit
        self._extend_base(width, left=True)

    def extend_right(self, width):
        #assert self.end0 + width <= self.end0_limit
        self._extend_base(width, left=False)

    def _make_extend_pileup(self, width, left):
        if left:
            start0 = self.start0 - width
            end0 = self.start0
        else:
            start0 = self.end0
            end0 = self.end0 + width

        return self.__class__(
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

    def get_extend_generator(self):
        return realign_base.PileupExtendGenerator(self)

    # active info related ones #

    def get_inactive_length_left(self):
        for idx in range(self.width):
            if self.active_info.iloc[idx]:
                break
        if idx == self.width - 1:
            idx += 1
        return idx

    def get_inactive_length_right(self):
        for idx in range(-1, -1 -self.width, -1):
            if self.active_info.iloc[idx]:
                break
        if idx == -self.width:
            idx -= 1
        return -idx - 1

    # extension and reduction in search for realignment area
    def augment_margins(self):
        #self.logger.debug(f'Beginning "augment_margins"')

        # immediately after seeding
        if self.active_info.any():
            #self.logger.debug(f'\tskipping initial search')
            pass
        else:
            #self.logger.debug(f'\tbeginning initial search')

            self.search_for_active_position()

            #self.logger.debug(f'\tfinished initial search')

            if not self.active_info.any():
                #self.logger.debug(f'Finished "augment_margins" without secure_margins')

                result_inactive = None
                result_vcfspec = None
                return result_inactive, result_vcfspec

        # trim excessive inactive padding
        self.reduce_left()
        self.reduce_right()
        # do after-discovery
        #self.logger.debug(f'\tbeginning after-discovery')
        result_inactive, result_vcfspec = self.secure_margins()
        #self.logger.debug(f'\tfinished after-discovery')

        #self.logger.debug(f'Finished "augment_margins"')
        return result_inactive, result_vcfspec

    def search_for_active_position(self):
        gen = self.get_extend_generator()
        if gen.aborted or self.active_info.any():
            return
        
        while True:
            if gen.aborted:
                break

            try:
                gen.iter_left()
            except StopIteration:
                pass
            else:
                if gen.left_is_active:
                    break

            try:
                gen.iter_right()
            except StopIteration:
                pass
            else:
                if gen.right_is_active:
                    break

        # subset
        self.subset(gen.current_start0, gen.current_end0, inplace=True)

    def secure_margins(self):
        #self.logger.debug(f'Began secure_margins')

        while True:
            #self.logger.debug(f'Began secure_inactive_padding')
            result_inactive = self.secure_inactive_padding()
            #self.logger.debug(f'Finished secure_inactive_padding. Margins: {self.range0}')

            #self.logger.debug(f'Began prepare_vcfspecs')
            self.prepare_vcfspecs()
            #self.logger.debug(f'Finished prepare_vcfspecs')

            #self.logger.debug(f'Began secure_vcfspec_margins')
            result_vcfspec = self.secure_vcfspec_margins()
            #self.logger.debug(f'Finished secure_vcfspec_margins. Margins: {self.range0}')

            if not result_vcfspec.edited:
                break
            #if (not result_inactive.left_okay) and (not result_inactive.right_okay):
            #    break
            #else:

        #self.logger.debug(f'Finished secure_margins')

        return result_inactive, result_vcfspec

    # extension reduction helpers
    def reduce_left(self):
        if self.width > self.params['inactive_padding']:
            inactive_margin_length = self.get_inactive_length_left()
            if inactive_margin_length == self.params['inactive_padding']:
                left_fixed = True
            elif inactive_margin_length > self.params['inactive_padding']:
                new_start0 = self.start0 + (inactive_margin_length - self.params['inactive_padding'])
                self.subset(new_start0, self.end0, inplace=True)
                left_fixed = True
            else:
                left_fixed = False
        else:
            left_fixed = False

        return left_fixed

    def reduce_right(self):
        if len(self.range0) > self.params['inactive_padding']:
            inactive_margin_length = self.get_inactive_length_right()
            if inactive_margin_length == self.params['inactive_padding']:
                right_fixed = True
            elif inactive_margin_length > self.params['inactive_padding']:
                new_end0 = self.end0 - (inactive_margin_length - self.params['inactive_padding'])
                self.subset(self.start0, new_end0, inplace=True)
                right_fixed = True
            else:
                right_fixed = False
        else:
            right_fixed = False

        return right_fixed

    @classmethod
    def _securehelper_check_aborted(cls, gen, left_okay, right_okay):
        return (
            (
                cls._securehelper_check_left_blocked(gen, left_okay) and
                cls._securehelper_check_right_blocked(gen, right_okay)
            ) or
            gen.touched_width_limit
        )

    @classmethod
    def _securehelper_check_left_blocked(cls, gen, left_okay):
        return gen.left_is_blocked or gen.touched_width_limit or left_okay

    @classmethod
    def _securehelper_check_right_blocked(cls, gen, right_okay):
        return gen.right_is_blocked or gen.touched_width_limit or right_okay

    def secure_inactive_padding(self):
        # set params
        inactive_length_left = self.get_inactive_length_left()
        inactive_length_right = self.get_inactive_length_right()
        left_okay = (inactive_length_left >= self.params['inactive_padding'])
        right_okay = (inactive_length_right >= self.params['inactive_padding'])
        
        # initial abort check
        if left_okay and right_okay:
            return realign_base.SecureResult(
                touched_left_limit=self.check_touches_left_limit(),
                touched_right_limit=self.check_touches_right_limit(), 
                touched_width_limit=self.check_touches_width_limit(), 

                left_low_depth=self.check_left_depth_low(),
                right_low_depth=self.check_right_depth_low(),
                left_low_MQ=self.check_left_MQ_low(),
                right_low_MQ=self.check_right_MQ_low(),

                left_okay=left_okay, 
                right_okay=right_okay,
                edited=False,
            )

        # begin looping
        initial_start0 = self.start0
        initial_end0 = self.end0
        gen = self.get_extend_generator()
        while True:
            if self._securehelper_check_aborted(gen, left_okay, right_okay):
                break

            if not self._securehelper_check_left_blocked(gen, left_okay):
                gen.iter_left()
                if gen.left_is_active:
                    inactive_length_left = 0
                elif not gen.left_is_active:
                    inactive_length_left += 1

                left_okay = (inactive_length_left >= self.params['inactive_padding'])

            if not self._securehelper_check_right_blocked(gen, right_okay):
                gen.iter_right()
                if gen.right_is_active:
                    inactive_length_right = 0
                elif not gen.right_is_active:
                    inactive_length_right += 1

                right_okay = (inactive_length_right >= self.params['inactive_padding'])

        if (initial_start0 != gen.current_start0) or (initial_end0 != gen.current_end0):
            edited = True
            self.subset(gen.current_start0, gen.current_end0, inplace=True)
        else:
            edited = False

        return realign_base.SecureResult(
            touched_left_limit=gen.touched_left_limit,
            touched_right_limit=gen.touched_right_limit, 
            touched_width_limit=gen.touched_width_limit, 

            left_low_depth=gen.left_low_depth,
            right_low_depth=gen.right_low_depth,
            left_low_MQ=gen.left_low_MQ,
            right_low_MQ=gen.right_low_MQ,

            left_okay=left_okay, 
            right_okay=right_okay,
            edited=edited,
        )

    def prepare_vcfspecs(self):
        self.set_row_specs()
        self.set_row_spec_groups()
        self.save_superseq_alignments(raise_with_tie=False)
        self.set_contig_vcfspecs()

    def secure_vcfspec_margins(self):
        vcfspec_margins_gdf = self.get_vcfspec_margins_gdf(
            subseq_portion_threshold=self.params['allele_portion_threshold'],
            inverse=False,
            split_contig_vcfspec=True,
        )

        if vcfspec_margins_gdf.is_empty:
            # When there is no contig vcfspec
            return realign_base.SecureResult(
                touched_left_limit=self.check_touches_left_limit(),
                touched_right_limit=self.check_touches_right_limit(), 
                touched_width_limit=self.check_touches_width_limit(), 

                left_low_depth=self.check_left_depth_low(),
                right_low_depth=self.check_right_depth_low(),
                left_low_MQ=self.check_left_MQ_low(),
                right_low_MQ=self.check_right_MQ_low(),

                left_okay=True, 
                right_okay=True,
                edited=False,
            )

        # set parameters
        desired_start0 = min(vcfspec_margins_gdf.start0s)  #min(candidate_start0s)
        #desired_end0 = max(candidate_end0s)
        desired_end0 = max(vcfspec_margins_gdf.end0s)
        left_okay = self.start0 <= desired_start0
        right_okay = self.end0 >= desired_end0

        # When requirements are already met
        if left_okay and right_okay:
            return realign_base.SecureResult(
                touched_left_limit=self.check_touches_left_limit(),
                touched_right_limit=self.check_touches_right_limit(), 
                touched_width_limit=self.check_touches_width_limit(), 

                left_low_depth=self.check_left_depth_low(),
                right_low_depth=self.check_right_depth_low(),
                left_low_MQ=self.check_left_MQ_low(),
                right_low_MQ=self.check_right_MQ_low(),

                left_okay=left_okay, 
                right_okay=right_okay,
                edited=False,
            )

        # begin looping
        initial_start0 = self.start0
        initial_end0 = self.end0
        gen = self.get_extend_generator()
        while True:
            if self._securehelper_check_aborted(gen, left_okay, right_okay):
                break

            if not self._securehelper_check_left_blocked(gen, left_okay):
                gen.iter_left()
                left_okay = gen.current_start0 <= desired_start0

            if not self._securehelper_check_right_blocked(gen, right_okay):
                gen.iter_right()
                right_okay = gen.current_end0 >= desired_end0

        if (initial_start0 != gen.current_start0) or (initial_end0 != gen.current_end0):
            edited = True
            self.subset(gen.current_start0, gen.current_end0, inplace=True)
        else:
            edited = False

        return realign_base.SecureResult(
            touched_left_limit=gen.touched_left_limit,
            touched_right_limit=gen.touched_right_limit, 
            touched_width_limit=gen.touched_width_limit, 

            left_low_depth=gen.left_low_depth,
            right_low_depth=gen.right_low_depth,
            left_low_MQ=gen.left_low_MQ,
            right_low_MQ=gen.right_low_MQ,

            left_okay=left_okay, 
            right_okay=right_okay,
            edited=edited,
        )

    def get_margin_from_vcfspec(self, vcfspec):
        leftmost = vcfspec.leftmost()
        rightmost = vcfspec.rightmost()
        padding = int(len(leftmost.REF_range0) * self.params['vcfspec_range_factor'])
        return (leftmost.start0 - padding, rightmost.end0 + padding)

    # pyranges
    #def get_gr(self):
    #    return pr.PyRanges(chromosomes=[self.chrom], starts=[self.start0], ends=[self.end0])
    def get_gdf(self):
        return GenomeDataFrameBase.from_data(
            chroms=self.chrom, start0s=self.start0, end0s=self.end0, refver=self.refver,
        )

    #def get_active_info_gr(self):
    def get_active_info_gdf(self):
        chroms = list()
        start0s = list()
        end0s = list()
        is_active = list()
        for key, subiter in itertools.groupby(
            self.active_info.items(), key=(lambda x: x[1]),
        ):
            subiter = tuple(subiter)
            chroms.append(self.chrom)
            start0s.append(subiter[0][0])
            end0s.append(subiter[-1][0] + 1)
            is_active.append(key)

        return GenomeDataFrameBase.from_data(chroms=chroms, start0s=start0s, end0s=end0s, Active=is_active)
        #return pr.from_dict(
        #    {'Chromosome': chroms, 'Start': start0s, 'End': end0s, 'Active': is_active}
        #)

    #def get_vcfspec_margins_gr(self, subseq_portion_threshold, inverse=False, split_contig_vcfspec=True):
    def get_vcfspec_margins_gdf(self, subseq_portion_threshold, inverse=False, split_contig_vcfspec=True):
        """Args:
            split_contig_vcfspec: If False, spaces between component vcfspecs of a contig vcfspec are incoporated into result pyranges.
            inverse: If True, subtraction from whole pileup region to vcfspec region is returned.
        """
        start0s = list()
        end0s = list()
        names = list()
        # collect vcfspec margins
        for row_id, contig_vcfspec in self.iter_contig_vcfspecs(subseq_portion_threshold):
            if contig_vcfspec is None:
                continue

            # make vcfspec iterator 
                # when contig vcfspec is composed of only one component vcfspec, 
                # iterator must be set differently
            if len(contig_vcfspec.components[1]) == 0:
                iterator = [contig_vcfspec]
            else:
                iterator = contig_vcfspec.components[1]

            if split_contig_vcfspec:
                for vcfspec in iterator:
                    margins = self.get_margin_from_vcfspec(vcfspec)
                    start0s.append(margins[0])
                    end0s.append(margins[1])
                    names.append(vcfspec.get_id())
            else:
                # for each contig vcfspec, spaces between component vcfspecs are
                # included in the result
                start0_candidates = list()
                end0_candidates = list()
                for vcfspec in iterator:
                    margins = self.get_margin_from_vcfspec(vcfspec)
                    start0_candidates.append(margins[0])
                    end0_candidates.append(margins[1])
                start0s.append(min(start0_candidates))
                end0s.append(max(end0_candidates))
                names.append(contig_vcfspec.get_id())

        # make pyranges object
        if len(start0s) == 0:
            #vcfspec_margins_gr = pr.PyRanges()
            vcfspec_margins_gdf = GenomeDataFrameBase.init_empty(refver=self.refver)
        else:
            #chroms = [self.chrom] * len(start0s)
            #vcfspec_margins_gr = pr.from_dict(
            #    {'Chromosome': chroms, 'Start': start0s, 'End': end0s, 'Name': names}
            #)
            vcfspec_margins_gdf = GenomeDataFrameBase.from_data(
                refver=self.refver, 
                chroms=self.chrom, 
                start0s=start0s, 
                end0s=end0s, 
                Name=names,
            )

        # return; inverse if needed
        if inverse:
            return self.get_gdf().subtract(vcfspec_margins_gdf)
        else:
            return vcfspec_margins_gdf


    # secure_inactive_padding
#    def secure_inactive_padding_rightward(self, inactive_padding, extend_pileup_by=DEFAULT_EXTEND_PILEUP_BY):
#        secure_inactive_padding_rightward(
#            pileup=self,
#            inactive_padding=inactive_padding, 
#            extend_pileup_by=extend_pileup_by, 
#            inplace=True,
#        )
#
#    def secure_inactive_padding_leftward(self, inactive_padding, extend_pileup_by=DEFAULT_EXTEND_PILEUP_BY):
#        secure_inactive_padding_leftward(
#            pileup=self, 
#            inactive_padding=inactive_padding, 
#            extend_pileup_by=extend_pileup_by, 
#            inplace=True,
#        )
#
#    def secure_inactive_padding(self, inactive_padding, extend_pileup_by=DEFAULT_EXTEND_PILEUP_BY):
#        self.secure_inactive_padding_rightward(inactive_padding, extend_pileup_by=extend_pileup_by)
#        self.secure_inactive_padding_leftward(inactive_padding, extend_pileup_by=extend_pileup_by)

    ### row specs ###
    def set_row_specs(self):
        self.row_specs = dict()
        for read_uid, row in self.dfs['queryonly_added_main'].iterrows():
            #row = self.df.loc[read_uid, :]
            self.row_specs[read_uid] = self._make_row_spec(row)

    def row_spec_getter(self, rowid):
        return self.row_specs[rowid]

    def set_row_spec_groups(self):
        """row_spec's are grouped by whether one is a substring of another."""
        self.row_spec_groups = collections.OrderedDict()

        for query_rowid, query in sorted(
            self.row_specs.items(), 
            #key=(lambda x: (len(x[1]['seq']), x[1]['span_length'])), 
            #key=(lambda x: (x[1]['span_length'], -len(x[1]['seq']))), 
            key=(lambda x: (x[1]['span_length'], -x[1]['cliplen'])), 
            reverse=True,
        ):
            superseq_candidates = list()
            for superseq_rowid in self.row_spec_groups.keys():
                target = self.row_specs[superseq_rowid]
                if self.row_spec_matcher(query, target):
                    superseq_candidates.append(superseq_rowid)
                #if len(superseq_candidates) >= 2:
                #    break
                    
            if len(superseq_candidates) == 0:
                self.row_spec_groups[query_rowid] = {
                    'subseq_rowids': list(),
                    'subseq_hits': 0,
                }
                self.row_spec_groups[query_rowid]['subseq_rowids'].append(query_rowid)
                self.row_spec_groups[query_rowid]['subseq_hits'] += 1
            else:
                superseq_rowid = self.superseq_candidates_tiebreak(superseq_candidates, self.row_spec_getter, self.row_spec_groups)
                groupinfo = self.row_spec_groups[superseq_rowid]
                groupinfo['subseq_rowids'].append(query_rowid)
                groupinfo['subseq_hits'] += 1

        self.decorate_groupinfo(self.row_spec_groups, self.read_getter)

    def read_getter(self, key):
        return self.read_store[key]

    # alignment #
    def align_row_spec_to_ref(self, row_spec, ref_seq, ref_seq_reversed, raise_with_tie):
        if row_spec['seq'] == '':
            return Bio.Align.Alignment(
                sequences=[ref_seq, ''],
                coordinates=np.array([(0, len(ref_seq)), (0, 0)]),
            )

        reverse_align = (not row_spec['left_filled']) and row_spec['right_filled']
        target = ref_seq
        query = row_spec['seq']
        return realign_base.main_aligner(
            target, 
            query, 
            self.aligner, 
            reverse_align=reverse_align, 
            raise_with_tie=raise_with_tie, 
            #logger=self.logger, 
            row_spec=row_spec, 
            target_reversed=ref_seq_reversed,
        )

    def align_row_spec_to_ref_old(self, row_spec, ref_seq, ref_seq_reversed, aligner, raise_with_tie):
        # set params
        query = row_spec['seq']
        reverse_align = (not row_spec['left_filled']) and row_spec['right_filled']
        # run aligner
        try:
            if reverse_align:
                alns = aligner.align(ref_seq_reversed, query[::-1])
            else:
                alns = aligner.align(ref_seq, query)
        except Exception as exc:
            raise Exception(f'Failed alignment:\nrow_spec: {row_spec}\nReference seq: {ref_seq}') from exc

        #self.logger.debug(f'Num of alignments: {len(alns)}; row_spec: {row_spec}')

        # treat dirty alignments
        if len(alns) > 10000:
            aln = random.choice(alns)
            if reverse_align:
                aln = alignhandler.reverse_alignment(aln)
            #aln = alignhandler.amend_outer_insdel_both(aln)
            #self.logger.debug(f'Skipped dirty alignment; row_spec: {row_spec}')
        else:
            # recover reversed alignment
            if reverse_align:
                alns = [alignhandler.reverse_alignment(x) for x in alns]

            if len(alns) == 1:
                aln = alns[0]
                aln = alignhandler.amend_outer_insdel_both(aln)
            else:
                #alns = [alignhandler.amend_outer_insdel_both(x) for x in alns]
                try:
                    aln = self.align_row_spec_to_ref_helper(alns, raise_with_tie, row_spec)
                except TimeoutError:
                    #self.logger.debug(f'skipping alignments tiebreaking due to timeout;\nrow_spec: {row_spec}')
                    aln = alns[0]
                    aln = alignhandler.amend_outer_insdel_both(aln)
            
        return aln

    @deco.get_deco_timeout(0.05)
    def align_row_spec_to_ref_helper(self, alns, raise_with_tie, row_spec):
        alns = list(
            alignhandler.remove_identical_alignments(
                alignhandler.amend_outer_insdel_both(x) for x in alns
            )
        )
        try:
            aln = alignhandler.alignment_tiebreaker(alns, raise_with_failure=raise_with_tie)
        except alignhandler.AlignmentTieError as exc:
            msg = f'Failed to break alignment tie. row_spec is:\n{row_spec}'
            raise Exception(msg) from exc

        return aln

    def save_superseq_alignments(self, raise_with_tie=False):
        """'_superseq_alignments' attribute is set."""
        self._superseq_alignments = dict()
        ref_seq = self.get_ref_seq()
        ref_seq_reversed = ref_seq[::-1]
        for superseq_rowid in self.row_spec_groups.keys():
            row_spec = self.row_specs[superseq_rowid]
            self._superseq_alignments[superseq_rowid] = self.align_row_spec_to_ref(row_spec, ref_seq, ref_seq_reversed, raise_with_tie)

    # not used #
    def save_subseq_alignments(self):
        self._subseq_alignments = dict()
        for superseq_rowid, groupinfo in self.row_spec_groups.items():
            sup_to_ref_aln = self._superseq_alignments[superseq_rowid]
            for subseq_rowid in groupinfo['subseq_rowids']:
                subseq_row_spec = self.row_specs[subseq_rowid]
                self._subseq_alignments[subseq_rowid] = self.align_subseq_to_ref(subseq_row_spec, sup_to_ref_aln)

    # vcfspec #
    def set_contig_vcfspecs(self):
        """- Must be run after row_spec_groups is set and save_superseq_alignments have been run
        """
        self.contig_vcfspecs = dict()
        for superseq_rowid, groupinfo in sorted(
            self.row_spec_groups.items(),
            key=(lambda x: x[1]['subseq_hits']),
            reverse=True,
        ):
            superseq_row_spec = self.row_specs[superseq_rowid]
            superseq_alignment = self._superseq_alignments[superseq_rowid]
            superseq_vcfspec = self.superseq_to_vcfspecs(superseq_row_spec, superseq_alignment)
            self.contig_vcfspecs[superseq_rowid] = superseq_vcfspec

#    def get_realigned_reads(self):
#        """Must be run after 'save_superseq_alignments' method
#        Returns:
#            dict (keys ReadUID, values pysam.AlignedSegment)
#        """
#        realigned_reads = dict()
#    #    for superseq_rowid, groupinfo in sorted(
#    #        active_region_pileup.row_spec_groups.items(),
#    #        key=(lambda x: x[1]['subseq_hits']), 
#    #        reverse=True,
#    #    ):
#        # sorting is not necessary
#
#        for superseq_rowid, groupinfo in self.row_spec_groups.items():
#            #realigned_read_superseq = realign_superseq(superseq_rowid, active_region_pileup, active_range)
#            #realigned_reads[superseq_rowid] = realigned_read_superseq
#                # (221019) This is commented out because superseq itself is
#                # included in subseq list.
#            superseq_aln = self._superseq_alignments[superseq_rowid]
#            for subseq_rowid in groupinfo['subseq_rowids']:
#                #if subseq_rowid in realigned_reads.keys():
#                #    continue
#                    # (221019) This is commented out because subseqs with multiple superseq hits are
#                    # incoporated into only one superseq group.
#                realigned_read_subseq = realign_subseq(subseq_rowid, self, self.range0, superseq_aln)
#                realigned_reads[subseq_rowid] = realigned_read_subseq
#                
#        return realigned_reads

    # split point generation
    def get_split_points(self):
        """Result may be empty"""
        # without splitting contig vcfspec
        #splittable_region = self.get_vcfspec_margins_gr(
        splittable_region = self.get_vcfspec_margins_gdf(
            subseq_portion_threshold=(self.params['allele_portion_threshold'] * 2),
            inverse=True,
            split_contig_vcfspec=False,
        )

        if splittable_region.is_empty:
            return list()
#            # split contig vcfspec
#            splittable_region = self.get_vcfspec_margins_gr(
#                subseq_portion_threshold=(self.allele_portion_threshold * 2),
#                inverse=True,
#                split_contig_vcfspec=True,
#            )
#            if splittable_region.empty:
#                logging.info(f'Getting Pileup splitting point from inactive runs because splitting with vcfspec margins failed. Pileup object being split: {self}')
#                # from inactive runs
#                #split_points = [self.get_split_point_from_inactive_runs()]
#                return list()
#            else:
#                logging.info(f'Getting Pileup splitting point, splitting contig vcfspecs. Pileup object being split: {self}')
#                pass
        else:
            pass

        split_points = list()
        for idx, row in splittable_region.df.iterrows():
            point = int((row.Start + row.End) / 2)
            if point != self.start0:
                split_points.append(point)

        return split_points

    def get_split_point_from_inactive_runs(self):
        """Midpoint of the widest non-marginal inactive area"""
        #active_info_gr = self.get_active_info_gr()
        #inactive_gr = active_info_gr[~active_info_gr.Active]
        active_info_gdf = self.get_active_info_gdf()
        inactive_gdf = active_info_gdf.loc[~active_info_gdf['Active'], :]

        if inactive_gdf.is_empty:
            raise realign_base.SparseInactiveRegionError(f'{self}')

        max_width_index = inactive_gdf.lengths.argmax()
        max_width_row = inactive_gdf.iloc[max_width_index, :]
        return int((max_width_row['Start'] + max_width_row['End']) / 2)

    def get_result_vcfspecs(self, as_components=False, merge=True, subseq_portion_threshold=None, MQ_threshold=40): 
        if subseq_portion_threshold is None:
            subseq_portion_threshold = self.params['allele_portion_threshold']

        if as_components:
            result = list()
            for row_id, contig_vcfspec in self.iter_contig_vcfspecs(subseq_portion_threshold=subseq_portion_threshold, MQ_threshold=MQ_threshold):
                if contig_vcfspec is None:
                    continue

                if len(contig_vcfspec.components[1]) == 0:
                    result.append(contig_vcfspec)
                else:
                    result.extend(contig_vcfspec.components[1])
        else:
            result = list(
                x[1] for x in self.iter_contig_vcfspecs(subseq_portion_threshold, MQ_threshold)
                if x[1] is not None
            )

        if merge:
            if len(result) > 0:
                result = functools.reduce(lambda x, y: libvcfspec.merge(x, y), result)

        return result

    # for debugging
    def show_row_spec_alignments(self, skip_zero_hits=True, show_subseqs=False):
        self._show_row_spec_alignments_helper(
            self.row_spec_groups, 
            lambda x: self.row_specs.__getitem__(x),
            lambda x: self._superseq_alignments.__getitem__(x),
            skip_zero_hits, 
            show_subseqs,
        )



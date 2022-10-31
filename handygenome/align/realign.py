import sys
import collections
import itertools
import functools
import logging

import pysam
import Bio.Align
import numpy as np
import pandas as pd
import pyranges as pr

import importlib
top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
libvcfspec = importlib.import_module(".".join([top_package_name, "variant", "vcfspec"]))
libpileup = importlib.import_module(".".join([top_package_name, "read", "pileup"]))
alignhandler = importlib.import_module(".".join([top_package_name, "align", "alignhandler"]))
bameditor = importlib.import_module(".".join([top_package_name, "bameditor"]))
readhandler = importlib.import_module(".".join([top_package_name, "read", "readhandler"]))
#fetchcache = importlib.import_module(".".join([top_package_name, "read", "fetchcache"]))


logging.basicConfig(level=logging.INFO)


DEFAULT_ACTIVE_THRESHOLD = 0.05
DEFAULT_INACTIVE_PADDING = 10
DEFAULT_VCFSPEC_EXTRACTION_THRESHOLD = 0.05
DEFAULT_VCFSPEC_CONCAT_DISTANCE = None
DEFAULT_MAX_PILEUP_WIDTH = 100  # Must be equal to or less than this value
DEFAULT_EXTEND_PILEUP_BY = 20
DEFAULT_AUGMENT_FACTOR_VCFSPEC_RANGE = 1.5
DEFAULT_MQ_LIMIT = 20
DEFAULT_DEPTH_LIMIT = 1

DEFAULT_RPILEUPSERIES_MAX_WIDTH = 5000


#DEFAULT_AUGMENT_FACTOR_REPEAT_AREA = 1


#ALIGNER_FILLED = Bio.Align.PairwiseAligner(
#    mode='local',
#    match_score=2,
#    mismatch_score=-3,
#    query_internal_open_gap_score=-7,
#    query_internal_extend_gap_score=-2,
#    target_internal_open_gap_score=-7,
#    target_internal_extend_gap_score=-2,
#    query_left_open_gap_score=-7,
#    query_left_extend_gap_score=-2,
#)

#ALIGNER_UNFILLED = Bio.Align.PairwiseAligner(
#    mode='local',
#    match_score=2,
#    mismatch_score=-3,
#    query_internal_open_gap_score=-7,
#    query_internal_extend_gap_score=-2,
#    target_internal_open_gap_score=-7,
#    target_internal_extend_gap_score=-2,
#)

ALIGNER_BLASTN = Bio.Align.PairwiseAligner(
    match_score=2,
    mismatch_score=-3,
    query_internal_open_gap_score=-7,
    query_internal_extend_gap_score=-2,
    target_internal_open_gap_score=-7,
    target_internal_extend_gap_score=-2,
)

ALIGNER_EQUAL_MM_GAP = Bio.Align.PairwiseAligner(
    mode='global',
    match_score=3,
    mismatch_score=-3,
    query_internal_open_gap_score=-3,
    query_internal_extend_gap_score=0,
    target_internal_open_gap_score=-3,
    target_internal_extend_gap_score=0,
)

ALIGNER_MAIN = ALIGNER_EQUAL_MM_GAP


class RealignerPileupBase(libpileup.PileupBase):
    row_spec_exluded_vals = (
        libpileup.PileupBase.DEL_VALUE, libpileup.PileupBase.EMPTY_VALUE
    )

    # row specs related ones
    @classmethod
    def _make_row_spec(cls, row):
        result = dict()
        result['seq'] = "".join(
            cls.pat_parenthesis.sub("", x) 
            for x in row 
            if x not in cls.row_spec_exluded_vals
        )
        result['id'] = row.name

        row_is_not_empty = row != cls.EMPTY_VALUE
        non_empty_index = row.index[row_is_not_empty]
        result['span_start0'] = non_empty_index.min()
        result['span_end0'] = non_empty_index.max() + 1
        result['span_length'] = result['span_end0'] - result['span_start0']
        result['left_filled'] = row_is_not_empty.iloc[0]
        result['right_filled'] = row_is_not_empty.iloc[-1]

        return result

    @staticmethod
    def row_spec_matcher(query, target):
        if (
            target['span_start0'] <= query['span_start0'] and
            target['span_end0'] >= query['span_end0']
        ):
            if query['left_filled'] and query['right_filled']:
                return target['seq'] == query['seq']
            else:
                return query['seq'] in target['seq']
        else:
            return False

    @staticmethod
    def handle_matching_subseq(superseq_candidates, row_spec_groups, subseq_rowid):
        if len(superseq_candidates) == 1:
            superseq_rowid = superseq_candidates[0]
        elif len(superseq_candidates) > 1:
            superseq_rowid = max(
                superseq_candidates, 
                key=(lambda x: row_spec_groups[x]['subseq_hits']),
            )

        row_spec_groups[superseq_rowid]['subseq_rowids'].append(subseq_rowid)
        row_spec_groups[superseq_rowid]['subseq_hits'] += 1

    @classmethod
    def group_row_specs(cls, row_specs):
        """row_spec's are grouped by whether one is a substring of another."""
        row_spec_groups = collections.OrderedDict()
        for query_rowid, query in sorted(
            row_specs.items(), 
            key=(lambda x: x[1]['span_length']), 
            reverse=True,
        ):
            superseq_candidates = list()
            for superseq_rowid in row_spec_groups.keys():
                target = row_specs[superseq_rowid]
                if cls.row_spec_matcher(query, target):
                    superseq_candidates.append(superseq_rowid)
                #if len(superseq_candidates) >= 2:
                #    break
                    
            if len(superseq_candidates) == 0:
                row_spec_groups[query_rowid] = {
                    'subseq_rowids': list(),
                    'subseq_hits': 0,
                }
                row_spec_groups[query_rowid]['subseq_rowids'].append(query_rowid)
                row_spec_groups[query_rowid]['subseq_hits'] += 1
            else:
                cls.handle_matching_subseq(superseq_candidates, row_spec_groups, query_rowid)

        # postprocess
        subseq_hits_sum = sum(x['subseq_hits'] for x in row_spec_groups.values())
        for groupinfo in row_spec_groups.values():
            groupinfo['subseq_hits_portion'] = groupinfo['subseq_hits'] / subseq_hits_sum

        return row_spec_groups

    def set_refined_vcfspec(self, row_concat_dist_le=np.inf, ):
        """All vcfspecs in a row are considered to be phased and belonging to a haplotype.
        """
        self.refined_vcfspec = None

    def superseq_to_vcfspecs(self, superseq_row_spec, superseq_alignment):
        lstrip_query_gaps = not superseq_row_spec["left_filled"]
        rstrip_query_gaps = not superseq_row_spec["right_filled"]
        superseq_vcfspecs = alignhandler.alignment_to_vcfspec(
            alignment=superseq_alignment, 
            target_start0=self.start0,
            chrom=self.chrom, 
            fasta=self.fasta,
            lstrip_query_gaps=lstrip_query_gaps,
            rstrip_query_gaps=rstrip_query_gaps,
        )
        if len(superseq_vcfspecs) > 0:
            superseq_vcfspecs = libvcfspec.concat_list(superseq_vcfspecs, distance_le=self.concat_dist_le)
        return superseq_vcfspecs

    # for debugging
    def _show_row_spec_alignments_helper(self, row_spec_groups, row_specs, skip_zero_hits, show_subseqs):
        subseq_hits_sum = sum(x['subseq_hits'] for x in row_spec_groups.values())
        for superseq_id, groupinfo in sorted(
            row_spec_groups.items(),
            key=(lambda x: x[1]['subseq_hits']),
            reverse=True,
        ):
            if skip_zero_hits:
                if groupinfo['subseq_hits'] == 0:
                    continue

            superseq_row_spec = row_specs[superseq_id]
            lstrip_query_gaps = not superseq_row_spec["left_filled"]
            rstrip_query_gaps = not superseq_row_spec["right_filled"]

            superseq_aln = self._superseq_alignments[superseq_id]
#            superseq_vcfspec_list = alignhandler.alignment_to_vcfspec(
#                superseq_aln, 
#                self.start0, 
#                self.chrom, 
#                self.fasta,
#                lstrip_query_gaps=lstrip_query_gaps,
#                rstrip_query_gaps=rstrip_query_gaps,
#            )

            superseq_vcfspec_list = self.superseq_to_vcfspecs(superseq_row_spec, superseq_aln)

            print('subseq hits:', groupinfo['subseq_hits'])
            print('subseq hits portion:', groupinfo['subseq_hits'] / subseq_hits_sum)
            print('row_spec:', superseq_row_spec)

            for vcfspec in superseq_vcfspec_list:
                print(f'concatenated vcfspec:\n\t{vcfspec}')
                print(f'concat components:')
                for component in vcfspec.concat_components[0]:
                    print('\t', component)

            print(superseq_aln)
            if show_subseqs:
                for rowid in groupinfo['subseq_rowids']:
                    print(rowid)
            print()


#class SecureResult(
#    collections.namedtuple(
#        'SecureResultBase',
#        ('touched_left_limit', 'touched_right_limit', 'touched_width_limit', 'left_okay', 'right_okay', 'edited'),
#    )
#):
#    pass


SecureResult = collections.namedtuple(
    'SecureResult',
    (
        'touched_width_limit', 
        'touched_left_limit', 
        'touched_right_limit', 

        'left_low_depth',  
        'right_low_depth',  
        'left_low_MQ',  
        'right_low_MQ',  

        'left_okay', 
        'right_okay', 
        'edited',
    ),
)


#ExtendGenResult = collections.namedtuple(
#    'ExtendGenResult',
#    ('aborted', 'no_iter', 'current_start0', 'current_end0', 'left_is_active', 'right_is_active', 'touched_left_limit', 'touched_right_limit', 'touched_width_limit'),
#)


class PileupExtendGenerator:
    def __init__(self, rpileup):
        self.pileup = rpileup

        self.current_start0 = self.pileup.start0
        self.current_end0 = self.pileup.end0
        self.touched_left_limit = self.pileup.check_touches_left_limit(start0=self.current_start0)
        self.touched_right_limit = self.pileup.check_touches_right_limit(end0=self.current_end0)

        self.interim_update(left=True, right=True)

        self.gen_left = self._base_generator_left()
        self.gen_right = self._base_generator_right()

    def __repr__(self):
        buffer = list()
        for key in (
            'current_start0',
            'current_end0',
            'touched_left_limit',
            'touched_right_limit',
            'touched_width_limit',
            'left_is_blocked',
            'right_is_blocked',
            'aborted',
        ):
            val = getattr(self, key)
            buffer.append(f'{key}={val}')
        string = ', '.join(buffer)
        return f'<PileupExtendGenerator({string})>'

    def stopiter(self):
        raise StopIteration()

    def iter_left(self):
        if self.left_is_blocked or self.touched_width_limit:
            self.stopiter()
        else:
            self.current_start0, self.touched_left_limit = next(self.gen_left)
            self.interim_update(left=True, right=False)
            #return self.make_genresult()

    def iter_right(self):
        if self.right_is_blocked or self.touched_width_limit:
            self.stopiter()
        else:
            self.current_end0, self.touched_right_limit = next(self.gen_right)
            self.interim_update(left=False, right=True)
            #return self.make_genresult()

    def interim_update(self, left, right):
        if left:
            self.left_is_active = self.pileup.active_info.loc[self.current_start0]
            self.left_low_MQ = self.pileup.MQ.loc[self.current_start0] < self.pileup.MQ_limit
            self.left_low_depth = self.pileup.get_depth(self.current_start0) < self.pileup.depth_limit

        if right:
            self.right_is_active = self.pileup.active_info.loc[self.current_end0 - 1]
            self.right_low_MQ = self.pileup.MQ.loc[self.current_end0 - 1] < self.pileup.MQ_limit
            self.right_low_depth = self.pileup.get_depth(self.current_end0 - 1) < self.pileup.depth_limit

    @property
    def touched_width_limit(self):
        return self.pileup.check_touches_width_limit(start0=self.current_start0, end0=self.current_end0)

    @property
    def left_is_blocked(self):
        return (
            self.touched_left_limit or 
            self.left_low_MQ or
            self.left_low_depth
        )

    @property
    def right_is_blocked(self):
        return (
            self.touched_right_limit or 
            self.right_low_MQ or
            self.right_low_depth
        )

    @property
    def aborted(self):
        return (
            self.touched_width_limit or 
            (self.left_is_blocked and self.right_is_blocked)
        )


#    def make_genresult(self):
#        return ExtendGenResult(
#            aborted=self.aborted, 
#            no_iter=False,
#            current_start0=self.current_start0, 
#            current_end0=self.current_end0, 
#            left_is_active=self.left_is_active, 
#            right_is_active=self.right_is_active, 
#            touched_left_limit=self.touched_left_limit, 
#            touched_right_limit=self.touched_right_limit, 
#            touched_width_limit=self.touched_width_limit,
#        )

#    @staticmethod
#    def make_genresult(rpileup, current_start0, current_end0, left_is_active, right_is_active, touched_left_limit, touched_right_limit):
#        touched_width_limit = rpileup.check_touches_width_limit(start0=current_start0, end0=current_end0)
#        aborted = (
#            touched_width_limit or 
#            (touched_left_limit and touched_right_limit)
#        )
#        return ExtendGenResult(
#            aborted=aborted, 
#            no_iter=False,
#            current_start0=current_start0, 
#            current_end0=current_end0, 
#            left_is_active=left_is_active, 
#            right_is_active=right_is_active, 
#            touched_left_limit=touched_left_limit, 
#            touched_right_limit=touched_right_limit, 
#            touched_width_limit=touched_width_limit,
#        )

    def _base_generator_left(self):
        """Stops only when current_start0 touches start0_limit"""
        current_start0 = self.pileup.start0
        if current_start0 <= self.pileup.start0_limit:
            return

        while True:
            self.pileup.extend_left(self.pileup.extend_pileup_by)
            for _ in range(self.pileup.extend_pileup_by):
                current_start0 -= 1
                #left_is_active = self.pileup.active_info.loc[current_start0]
                touched_left_limit = current_start0 <= self.pileup.start0_limit
                yield current_start0, touched_left_limit
                if touched_left_limit:
                    return

    def _base_generator_right(self):
        """Stops only when current_end0 touches end0_limit"""
        current_end0 = self.pileup.end0
        if current_end0 >= self.pileup.end0_limit:
            return

        while True:
            self.pileup.extend_right(self.pileup.extend_pileup_by)
            for _ in range(self.pileup.extend_pileup_by):
                current_end0 += 1
                #right_is_active = self.pileup.active_info.loc[current_end0]
                touched_right_limit = current_end0 >= self.pileup.end0_limit
                yield current_end0, touched_right_limit
                if touched_right_limit:
                    return


class RealignerPileup(RealignerPileupBase):
    # initializers
    def __init__(
        self, 
        fasta, 
        bam, 
        chrom, 
        start0=None,
        end0=None,
        init_df=True,
        active_threshold=DEFAULT_ACTIVE_THRESHOLD, 
        aligner=ALIGNER_MAIN,
        inactive_padding=DEFAULT_INACTIVE_PADDING,
        allele_portion_threshold=DEFAULT_VCFSPEC_EXTRACTION_THRESHOLD, 
        concat_dist_le=DEFAULT_VCFSPEC_CONCAT_DISTANCE,
        max_pileup_width=DEFAULT_MAX_PILEUP_WIDTH,
        extend_pileup_by=DEFAULT_EXTEND_PILEUP_BY,
        vcfspec_range_factor=DEFAULT_AUGMENT_FACTOR_VCFSPEC_RANGE,
        start0_limit=-np.inf,  # This value is allowed
        end0_limit=np.inf,  # This value is allowed
        MQ_limit=DEFAULT_MQ_LIMIT,  # This value is allowed
        depth_limit=DEFAULT_DEPTH_LIMIT,  # This value is allowed
    ):
        libpileup.PileupBase.__init__(self, fasta, bam, chrom, start0=start0, end0=end0, init_df=init_df)

        self.active_threshold = active_threshold
        if init_df:
            self._set_active_info()

        self.aligner = aligner
        self.inactive_padding = inactive_padding
        self.allele_portion_threshold = allele_portion_threshold
        self.concat_dist_le = concat_dist_le
        self.max_pileup_width = max_pileup_width
        self.extend_pileup_by = extend_pileup_by
        self.vcfspec_range_factor = vcfspec_range_factor
        self.start0_limit = start0_limit
        self.end0_limit = end0_limit
        self.MQ_limit = MQ_limit
        self.depth_limit = depth_limit

        self._vcfspec_cache = dict()
        self.hit_left_margin = False
        self.hit_right_margin = False
    
    def subset(self, start0, end0, inplace=False):
        new_active_info = self.active_info.loc[start0:(end0 - 1)]
        if inplace:
            libpileup.PileupBase._subset_base(self, start0, end0, inplace=inplace)
            self.active_info = new_active_info
        else:
            result = libpileup.PileupBase._subset_base(self, start0, end0, inplace=inplace)
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
            pileup.prepare_vcfspecs()
        return split_pileups

    ################
    def check_touches_left_limit(self, start0=None):
        if start0 is None:
            start0 = self.start0
        return start0 <= self.start0_limit

    def check_touches_right_limit(self, end0=None):
        if end0 is None:
            end0 = self.end0
        return end0 >= self.end0_limit

    def check_touches_width_limit(self, start0=None, end0=None):
        if start0 is None:
            start0 = self.start0
        if end0 is None:
            end0 = self.end0
        return end0 - start0 >= self.max_pileup_width

    def check_left_MQ_low(self):
        return self.MQ.iloc[0] < self.MQ_limit

    def check_right_MQ_low(self):
        return self.MQ.iloc[-1] < self.MQ_limit

    def check_left_depth_low(self):
        return self.get_depth(self.start0) < self.depth_limit

    def check_right_depth_low(self):
        return self.get_depth(self.end0 - 1) < self.depth_limit

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
            fasta=self.fasta, 
            bam=self.bam, 
            chrom=self.chrom, 
            start0=start0,
            end0=end0,
            init_df=True,
            active_threshold=self.active_threshold,
        )

    def get_extend_generator(self):
        return PileupExtendGenerator(self)

    # active info related ones #
    def _active_info_generator(self, start0, end0, reverse=False):
        self._coord_arg_sanitycheck(start0, end0)

        ref_seq = self.fasta.fetch(self.chrom, start0, end0)
        pos0_range = range(start0, end0)
        if reverse:
            ref_seq = ref_seq[::-1]
            pos0_range = pos0_range[::-1]

        for ref_base, pos0 in zip(ref_seq, pos0_range):
            counter = self.get_allele_counter(pos0)
            non_ref_portion = 1 - (counter[ref_base] / sum(counter.values()))
            is_active = (non_ref_portion >= self.active_threshold)
            yield is_active 

    def _set_active_info(self):
        self.active_info = pd.Series(
            self._active_info_generator(self.start0, self.end0, reverse=False),
            index=self.df.columns,
        )

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

#    def check_inactive_margin_right(self, length):
#        if len(self.active_info) < length:
#            return False
#        else:
#            return not self.active_info[-length:].any()
#
#    def check_inactive_margin_left(self, length):
#        if len(self.active_info) < length:
#            return False
#        else:
#            return not self.active_info[:length].any()

    # extension and reduction in search for realignment area
    def augment_margins(self):
        #logging.debug(f'Beginning "augment_margins"')
        # immediately after seeding
        if self.active_info.any():
            pass
            #logging.debug(f'\tskipping initial search')
        else:
            #logging.debug(f'\tbeginning initial search')
            self.search_for_active_position()
            #logging.debug(f'\tfinished initial search')
            if not self.active_info.any():
                result_inactive = None
                result_vcfspec = None
                #logging.debug(f'Finished "augment_margins" without secure_margins')
                return result_inactive, result_vcfspec

        # trim excessive inactive padding
        self.reduce_left()
        self.reduce_right()
        # do after-discovery
        #logging.debug(f'\tbeginning after-discovery')
        result_inactive, result_vcfspec = self.secure_margins()
        #logging.debug(f'\tfinished after-discovery')
        #logging.debug(f'Finished "augment_margins"')
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
        while True:
            result_inactive = self.secure_inactive_padding()
            self.prepare_vcfspecs()
            result_vcfspec = self.secure_vcfspec_margins()
            if not result_vcfspec.edited:
                break
            #if (not result_inactive.left_okay) and (not result_inactive.right_okay):
            #    break
            #else:
        return result_inactive, result_vcfspec

    # extension reduction helpers
    def reduce_left(self, inactive_padding=None):
        if self.width > self.inactive_padding:
            inactive_margin_length = self.get_inactive_length_left()
            if inactive_margin_length == self.inactive_padding:
                left_fixed = True
            elif inactive_margin_length > self.inactive_padding:
                new_start0 = self.start0 + (inactive_margin_length - self.inactive_padding)
                self.subset(new_start0, self.end0, inplace=True)
                left_fixed = True
            else:
                left_fixed = False
        else:
            left_fixed = False

        return left_fixed

    def reduce_right(self, inactive_padding=None):
        if len(self.range0) > self.inactive_padding:
            inactive_margin_length = self.get_inactive_length_right()
            if inactive_margin_length == self.inactive_padding:
                right_fixed = True
            elif inactive_margin_length > self.inactive_padding:
                new_end0 = self.end0 - (inactive_margin_length - self.inactive_padding)
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
        return gen.left_is_blocked or left_okay

    @classmethod
    def _securehelper_check_right_blocked(cls, gen, right_okay):
        return gen.right_is_blocked or right_okay

    def secure_inactive_padding(self):
        # set params
        inactive_length_left = self.get_inactive_length_left()
        inactive_length_right = self.get_inactive_length_right()
        left_okay = (inactive_length_left >= self.inactive_padding)
        right_okay = (inactive_length_right >= self.inactive_padding)
        
        # initial abort check
        if left_okay and right_okay:
            return SecureResult(
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

                left_okay = (inactive_length_left >= self.inactive_padding)

            if not self._securehelper_check_right_blocked(gen, right_okay):
                gen.iter_right()
                if gen.right_is_active:
                    inactive_length_right = 0
                elif not gen.right_is_active:
                    inactive_length_right += 1

                right_okay = (inactive_length_right >= self.inactive_padding)

        if (initial_start0 != gen.current_start0) or (initial_end0 != gen.current_end0):
            edited = True
            self.subset(gen.current_start0, gen.current_end0, inplace=True)
        else:
            edited = False

        return SecureResult(
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
        vcfspec_margins_gr = self.get_vcfspec_margins_gr()

        if vcfspec_margins_gr.empty:
            # When there is no contig vcfspec
            return SecureResult(
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
        desired_start0 = min(vcfspec_margins_gr.Start)  #min(candidate_start0s)
        #desired_end0 = max(candidate_end0s)
        desired_end0 = max(vcfspec_margins_gr.End)
        left_okay = self.start0 <= desired_start0
        right_okay = self.end0 >= desired_end0

        # When requirements are already met
        if left_okay and right_okay:
            return SecureResult(
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

        return SecureResult(
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
        padding = int(len(leftmost.REF_range0) * self.vcfspec_range_factor)
        return (leftmost.start0 - padding, rightmost.end0 + padding)

    # pyranges
    def get_gr(self):
        return pr.PyRanges(chromosomes=[self.chrom], starts=[self.start0], ends=[self.end0])

    def get_active_info_gr(self):
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
        return pr.from_dict(
            {'Chromosome': chroms, 'Start': start0s, 'End': end0s, 'Active': is_active}
        )

    def get_vcfspec_margins_gr(self, subseq_hits_portion_threshold=None, inverse=False):
        if subseq_hits_portion_threshold is None:
            subseq_hits_portion_threshold = self.allele_portion_threshold

        start0s = list()
        end0s = list()
        names = list()
        for row_id, contig_vcfspec_list in self.contig_vcfspecs.items():
            subseq_hits_portion = self.row_spec_groups[row_id]['subseq_hits_portion']
            if subseq_hits_portion < subseq_hits_portion_threshold:
                continue
            contig_vcfspec = contig_vcfspec_list[0]
            for component in contig_vcfspec.concat_components[0]:
                margins = self.get_margin_from_vcfspec(component)
                start0s.append(margins[0])
                end0s.append(margins[1])
                names.append(component.get_id())

        if len(start0s) == 0:
            vcfspec_margins_gr = pr.PyRanges()
        else:
            chroms = [self.chrom] * len(start0s)
            vcfspec_margins_gr = pr.from_dict(
                {'Chromosome': chroms, 'Start': start0s, 'End': end0s, 'Name': names}
            )

        if inverse:
            return self.get_gr().subtract(vcfspec_margins_gr)
        else:
            return vcfspec_margins_gr


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
        for read_uid, row in self.df.iterrows():
            row = self.df.loc[read_uid, :]
            self.row_specs[read_uid] = self._make_row_spec(row)

    def set_row_spec_groups(self):
        self.row_spec_groups = self.__class__.group_row_specs(self.row_specs)

    # alignment and vcfspec related
    @staticmethod
    def align_row_spec_to_ref(row_spec, ref_seq, ref_seq_reversed, aligner, raise_with_tie):
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
        # recover reversed alignment
        if reverse_align:
            alns = [alignhandler.reverse_alignment(x) for x in alns]

        if len(alns) == 1:
            aln = alns[0]
            aln = alignhandler.amend_outer_insdel_both(aln)
        else:
            #alns = [alignhandler.amend_outer_insdel_both(x) for x in alns]
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

    @staticmethod
    def align_subseq_to_ref(subseq_row_spec, sup_to_ref_aln):
        # get subseq index
        reverse_align = (not subseq_row_spec['left_filled']) and subseq_row_spec['right_filled']
        if reverse_align:
            sub_to_sup_index = sup_to_ref_aln.query.rfind(subseq_row_spec['seq'])
        else:
            sub_to_sup_index = sup_to_ref_aln.query.find(subseq_row_spec['seq'])
        if sub_to_sup_index == -1:
            raise Exception(
                f'"subseq" is not a subseq of "superseq":\n'
                f'subseq_row_spec: {subseq_row_spec}\n'
                f'superseq: {sup_to_ref_aln.query}'
            )
        # convert into sub-to-ref alignment
        # supaln: sup_to_ref_aln
        # subaln: sub_to_sup_aln
        active_region_walks = list()
        target_walk_subaln = range(sub_to_sup_index, sub_to_sup_index + len(subseq_row_spec['seq']))
        query_walk_subaln = range(0, len(subseq_row_spec['seq']))
        found_initial_hit = False
        for target_walk_supaln, query_walk_supaln in alignhandler.get_walks(sup_to_ref_aln, copy=False):
            if not found_initial_hit:
                if target_walk_subaln.start in query_walk_supaln:
                    found_initial_hit = True
                    # superseq walk
                    superseq_walk_start = target_walk_subaln.start
                    superseq_walk_stop = min(query_walk_supaln.stop, target_walk_subaln.stop)
                    len_current_superseq_walk = superseq_walk_stop - superseq_walk_start
                    # subseq walk
                    subseq_walk = range(0, len_current_superseq_walk)
                    # ref walk
                    if len(target_walk_supaln) == 0:
                        ref_walk = target_walk_supaln
                    else:
                        # During initial hit query_walk_supaln cannot be zero-length
                        ref_walk_start = target_walk_supaln.start + (superseq_walk_start - query_walk_supaln.start)
                        ref_walk = range(ref_walk_start, ref_walk_start + len_current_superseq_walk)

                    active_region_walks.append((ref_walk, subseq_walk))

                    if query_walk_supaln.stop >= target_walk_subaln.stop:
                        break
            else:
                # superseq walk
                superseq_walk_start = query_walk_supaln.start
                superseq_walk_stop = min(query_walk_supaln.stop, target_walk_subaln.stop)
                len_current_superseq_walk = superseq_walk_stop - superseq_walk_start
                # subseq walk
                subseq_walk_start = active_region_walks[-1][1].stop
                subseq_walk = range(subseq_walk_start, subseq_walk_start + len_current_superseq_walk)
                # ref walk
                if len(target_walk_supaln) == 0 or len(query_walk_supaln) == 0:
                    ref_walk = target_walk_supaln
                else:
                    ref_walk_start = target_walk_supaln.start
                    ref_walk = range(ref_walk_start, ref_walk_start + len_current_superseq_walk)

                active_region_walks.append((ref_walk, subseq_walk))

                if query_walk_supaln.stop >= target_walk_subaln.stop:
                    break
        # pad with gap: right side
        last_ref_walk = active_region_walks[-1][0]
        last_subseq_walk = active_region_walks[-1][1]
        if last_ref_walk.stop < len(sup_to_ref_aln.target):
            added_ref_walk = range(last_ref_walk.stop, len(sup_to_ref_aln.target))
            added_subseq_walk = range(last_subseq_walk.stop, last_subseq_walk.stop)
            active_region_walks.append((added_ref_walk, added_subseq_walk))
        # pad with gap: left
        first_ref_walk = active_region_walks[0][0]
        first_subseq_walk = active_region_walks[0][1]
        if first_ref_walk.start > 0:
            added_ref_walk = range(0, first_ref_walk.start)
            added_subseq_walk = range(0, 0)
            active_region_walks.insert(0, (added_ref_walk, added_subseq_walk))
        # make result alignment
        return Bio.Align.PairwiseAlignment(
            sup_to_ref_aln.target,
            subseq_row_spec['seq'],
            alignhandler.walks_to_path(active_region_walks),
            0
        )

    def save_superseq_alignments(self, raise_with_tie=False):
        """'_superseq_alignments' attribute is set."""
        self._superseq_alignments = dict()
        ref_seq = self.get_ref_seq()
        ref_seq_reversed = ref_seq[::-1]
        for superseq_rowid in self.row_spec_groups.keys():
            row_spec = self.row_specs[superseq_rowid]
            self._superseq_alignments[superseq_rowid] = self.align_row_spec_to_ref(row_spec, ref_seq, ref_seq_reversed, self.aligner, raise_with_tie)

    def save_subseq_alignments(self):
        self._subseq_alignments = dict()
        for superseq_rowid, groupinfo in self.row_spec_groups.items():
            sup_to_ref_aln = self._superseq_alignments[superseq_rowid]
            for subseq_rowid in groupinfo['subseq_rowids']:
                subseq_row_spec = self.row_specs[subseq_rowid]
                self._subseq_alignments[subseq_rowid] = self.align_subseq_to_ref(subseq_row_spec, sup_to_ref_aln)

    def set_contig_vcfspecs(self):
        """- Must be run after row_spec_groups is set and save_superseq_alignments have been run
        - "_vcfspec_cache" attribute is also updated
        """
        #self.contig_vcfspecs = list()
        self.contig_vcfspecs = dict()
        #self.all_contig_vcfspecs = dict()
        #self.selected_contig_vcfspecs = dict()

        #subseq_hits_sum = sum(groupinfo['subseq_hits'] for groupinfo in self.row_spec_groups.values())
        for superseq_rowid, groupinfo in sorted(
            self.row_spec_groups.items(),
            key=(lambda x: x[1]['subseq_hits']),
            reverse=True,
        ):
            #allele_portion = groupinfo['subseq_hits_portion']
            #if allele_portion < self.allele_portion_threshold:
            #    continue

            superseq_row_spec = self.row_specs[superseq_rowid]
            superseq_alignment = self._superseq_alignments[superseq_rowid]
            superseq_vcfspecs = self.superseq_to_vcfspecs(superseq_row_spec, superseq_alignment)
            self._vcfspec_cache[superseq_rowid] = superseq_vcfspecs
            if len(superseq_vcfspecs) > 0:
                #self.contig_vcfspecs.append(superseq_vcfspecs)
                self.contig_vcfspecs[superseq_rowid] = superseq_vcfspecs
                #self.all_contig_vcfspecs[superseq_rowid] = superseq_vcfspecs
                #if groupinfo['subseq_hits_portion'] >= self.allele_portion_threshold:
                #    self.contig_vcfspecs[superseq_rowid] = superseq_vcfspecs

    def show_contig_vcfspecs(self):
        if len(self.contig_vcfspecs) == 0:
            print('empty')
        else:
            #subseq_hits_sum = sum(groupinfo['subseq_hits'] for groupinfo in self.row_spec_groups.values())

            #for concat_vcfspec in itertools.chain.from_iterable(self.contig_vcfspecs.values()):
            for row_id, concat_vcfspec_list in self.contig_vcfspecs.items():
                subseq_hits = self.row_spec_groups[row_id]['subseq_hits']
                subseq_hits_portion = self.row_spec_groups[row_id]['subseq_hits_portion']
                concat_vcfspec = concat_vcfspec_list[0]

                print(concat_vcfspec)
                print(f'subseq_hits: {subseq_hits}; portion: {subseq_hits_portion}')
                print('components:')
                for component in concat_vcfspec.concat_components[0]:
                    print(f'\t{component}')
                print()

    def get_realigned_reads(self):
        """Must be run after 'save_superseq_alignments' method
        Returns:
            dict (keys ReadUID, values pysam.AlignedSegment)
        """
        return get_realigned_reads_helper(self, self.row_spec_groups, self._superseq_alignments)

    # for debugging
    def show_row_spec_alignments(self, skip_zero_hits=True, show_subseqs=True):
        self._show_row_spec_alignments_helper(self.row_spec_groups, self.row_specs, skip_zero_hits, show_subseqs)


class RealignerPileupSeries():
    def __init__(self, fasta, bam, chrom, start0, end0, rpileup_init_kwargs=dict(), max_width=DEFAULT_RPILEUPSERIES_MAX_WIDTH):
        # set params
        self.max_width = max_width
        self.fasta = fasta
        self.bam = bam
        self.chrom = chrom
        # initiate pileup_list
        initial_pileup, result_inactive, result_vcfspec = get_realigner_pileup(fasta, bam, chrom, start0, end0, **rpileup_init_kwargs)

        #print('Initial pileup generation')
        #print(initial_pileup)
        #print('result_inactive', result_inactive)
        #print('result_vcfspec', result_vcfspec)

        if initial_pileup is None:
            self.pileup_list = list()
            return

        left_insufficient = (not result_inactive.left_okay) or (not result_vcfspec.left_okay)
        right_insufficient = (not result_inactive.right_okay) or (not result_vcfspec.right_okay)

        if (not left_insufficient) and (not right_insufficient):
            self.pileup_list = [initial_pileup]
        else:
            self.pileup_list = self.split_initial_pileup(initial_pileup, left_insufficient, right_insufficient)

            #print('Initial pileup list')
            #print(self.pileup_list)

            self.secure_leftward()
            self.secure_rightward()
        # merge individual read_store's
        self.read_store = dict()
        for pileup in self.pileup_list:
            self.read_store.update(pileup.read_store)

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

    @property
    def left_pileup(self):
        return self.pileup_list[0]

    @property
    def right_pileup(self):
        return self.pileup_list[-1]

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

    # init helpers #
    def secure_leftward(self):
        left_pileup = self.pileup_list[0]
        left_pileup.end0_limit = left_pileup.end0
        while True:
            result_inactive, result_vcfspec = left_pileup.secure_margins()

            #print('secure_leftward loop')
            #print(self.pileup_list)
            #print('result_inactive', result_inactive)
            #print('result_vcfspec', result_vcfspec)
            #print()

            if (
                (result_inactive.left_okay and result_vcfspec.left_okay) or
                (result_vcfspec.left_low_depth or result_vcfspec.left_low_MQ) or
                (self.width > self.max_width)
            ):
                break
            else:
                split_pileups = self.split_left_pileup(left_pileup)
                del self.pileup_list[0]
                self.pileup_list = split_pileups + self.pileup_list
                left_pileup = self.pileup_list[0]
                left_pileup.end0_limit = left_pileup.end0

    def secure_rightward(self):
        #print('BEGINNING of secure_rightward')

        right_pileup = self.pileup_list[-1]
        #print(right_pileup)
        right_pileup.start0_limit = right_pileup.start0
        while True:
            result_inactive, result_vcfspec = right_pileup.secure_margins()

            #print('secure_rightward loop')
            #print(self.pileup_list)
            #print('result_inactive', result_inactive)
            #print('result_vcfspec', result_vcfspec)
            #print()

            if (
                (result_inactive.right_okay and result_vcfspec.right_okay) or
                (result_vcfspec.right_low_depth or result_vcfspec.right_low_MQ) or
                (self.width > self.max_width)
            ):
                break
            else:
                split_pileups = self.split_right_pileup(right_pileup)
                del self.pileup_list[-1]
                self.pileup_list.extend(split_pileups)
                right_pileup = self.pileup_list[-1]
                right_pileup.start0_limit = right_pileup.start0

    @classmethod
    def split_initial_pileup(cls, rpileup, left_insufficient, right_insufficient):
        split_points = cls.get_split_points(rpileup)
        if len(split_points) > 1:
            if (not left_insufficient) and right_insufficient:
                split_points = [split_points[-1]]
            elif left_insufficient and (not right_insufficient):
                split_points = [split_points[0]]
        return rpileup.split(split_points)

    @classmethod
    def split_left_pileup(cls, rpileup):
        split_points = cls.get_split_points(rpileup)
        if len(split_points) > 1:
            split_points = [split_points[0]]
        return rpileup.split(split_points)

    @classmethod
    def split_right_pileup(cls, rpileup):
        split_points = cls.get_split_points(rpileup)
        if len(split_points) > 1:
            split_points = [split_points[-1]]
        return rpileup.split(split_points)

    @classmethod
    def get_split_points(cls, rpileup):
        splittable_region = cls.get_split_region_vcfspec_margin(rpileup)
        if splittable_region.empty:
            logging.info(
                f'Getting RealignerPileup splitting point from inactive runs because splitting with vcfspec margins failed. Pileup object being split: {rpileup}'
            )
            # from inactive runs
            split_points = [cls.get_split_point_from_inactive_runs(rpileup)]
        else:
            # from region free of vcfspec flanks
            split_points = list()
            for idx, row in splittable_region.df.iterrows():
                split_points.append(int((row.Start + row.End) / 2))

        return split_points

    @staticmethod
    def get_split_region_vcfspec_margin(rpileup):
        vcfspec_margins_gr = rpileup.get_vcfspec_margins_gr(subseq_hits_portion_threshold=(rpileup.allele_portion_threshold * 2))
        splittable_region = rpileup.get_gr().subtract(vcfspec_margins_gr)

        return splittable_region

    @staticmethod
    def get_split_point_from_inactive_runs(rpileup):
        """Midpoint of the widest non-marginal inactive area"""
        active_info_gr = rpileup.get_active_info_gr()
        inactive_gr = active_info_gr[~active_info_gr.Active]
        #inactive_gr = inactive_gr[inactive_gr.Start != rpileup.start0]
        #inactive_gr = inactive_gr[inactive_gr.End != rpileup.end0]

        if inactive_gr.empty:
            raise SparseInactiveRegionError(f'{rpileup}')

        max_width_index = inactive_gr.lengths().argmax()
        max_width_row = inactive_gr.df.iloc[max_width_index, :]
        return int((max_width_row.Start + max_width_row.End) / 2)

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

    @staticmethod
    def merge_two_cigartuples(cigartuples_left, cigartuples_right):
        if len(cigartuples_left) == 0 or len(cigartuples_right) == 0:
            return cigartuples_left + cigartuples_right
        else:
            if cigartuples_left[-1][0] == cigartuples_right[0][0]:
                return (
                    cigartuples_left[:-1]
                    + [(cigartuples_left[-1][0], cigartuples_left[-1][1] + cigartuples_right[0][1])]
                    + cigartuples_right[1:]
                )
            else:
                return cigartuples_left + cigartuples_right

    def set_realigned_reads(self):
        # prepare subseq alignments
        for pileup in self.pileup_list:
            pileup.save_subseq_alignments()

        self.realigned_reads = dict()
        for read_uid, read in self.read_store.items():
            before_cigartuples, after_cigartuples = alignhandler.split_read_cigartuples(read, range(self.start0, self.end0))

            # select relevant pileups 
            relevant_pileups = list()
            for pileup in self.pileup_list:
                if read_uid in pileup.df.index:
                    #relevant_pileups.append((pileup.range0, pileup._subseq_alignments[read_uid]))
                    relevant_pileups.append(pileup)
            if len(relevant_pileups) == 0:
                # During 'subset' method of rpileup, out-of-range reads remain in 'read_store'
                continue

            # make cigartuples from alignment within each pileup
            cigartuples_list = list()
            if len(relevant_pileups) == 1:
                cigartuples_list.append(
                    self.subseq_aln_to_cigartuples(
                        alignment=relevant_pileups[0]._subseq_alignments[read_uid], 
                        left_is_empty=(not bool(before_cigartuples)),
                        right_is_empty=(not bool(after_cigartuples)),
                    )
                )
            else:
                cigartuples_list.append(
                    self.subseq_aln_to_cigartuples(
                        alignment=relevant_pileups[0]._subseq_alignments[read_uid], 
                        left_is_empty=(not bool(before_cigartuples)),
                        right_is_empty=False,
                    )
                )
                for pileup in relevant_pileups[1:-1]:
                    cigartuples_list.append(
                        self.subseq_aln_to_cigartuples(
                            alignment=pileup._subseq_alignments[read_uid], 
                            left_is_empty=False,
                            right_is_empty=False,
                        )
                    )
                cigartuples_list.append(
                    self.subseq_aln_to_cigartuples(
                        alignment=relevant_pileups[-1]._subseq_alignments[read_uid], 
                        left_is_empty=False,
                        right_is_empty=(not bool(after_cigartuples)),
                    )
                )

            # merge cigartuples
            realigned_cigartuples = functools.reduce(
                self.merge_two_cigartuples, 
                [before_cigartuples] + [x[0] for x in cigartuples_list] + [after_cigartuples],
            )
            # get new read start position
            if len(before_cigartuples) > 0:
                new_reference_start = read.reference_start
            else:
                new_reference_start = relevant_pileups[0].start0 + cigartuples_list[0][1]
            # make new read object
            realigned_read = read.__copy__()
            realigned_read.reference_start = new_reference_start
            realigned_read.cigartuples = realigned_cigartuples
            readhandler.set_NMMD(realigned_read, self.fasta)

            realigned_cigar_sanity_check(read, realigned_read, realigned_cigartuples)

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


class SparseInactiveRegionError(Exception):
    pass


class RealignmentCigarError(Exception):
    pass


def get_realigner_pileup(fasta, bam, chrom, start0, end0, **kwargs):
    rpileup = RealignerPileup(fasta, bam, chrom, start0, end0, **kwargs)
    result_inactive, result_vcfspec = rpileup.augment_margins()
    if (result_inactive is None) and (result_vcfspec is None):
        # No active region found
        rpileup = None

    return rpileup, result_inactive, result_vcfspec


def get_active_region_pileup_multisample(
    chrom,
    start0,
    end0,
    bam_dict,
    fasta,
    active_threshold=DEFAULT_ACTIVE_THRESHOLD,
    inactive_padding=DEFAULT_INACTIVE_PADDING,
    extend_pileup_by=DEFAULT_EXTEND_PILEUP_BY,
    augment_factor_vcfspec_range=DEFAULT_AUGMENT_FACTOR_VCFSPEC_RANGE,
    #augment_factor_repeat_area=DEFAULT_AUGMENT_FACTOR_REPEAT_AREA,
    aligner=None,
    allele_portion_threshold=None,
    vcfspec_concat_dist=None,
    raise_with_tie=False,
    max_pileup_width=DEFAULT_MAX_PILEUP_WIDTH,
):
    # initialize pileup_dict
    pileup_dict = dict()
    for sampleid, bam in bam_dict.items():
        logging.debug(f'Initiating pileup for single sample: {sampleid}')
        active_range, active_region_pileup = get_active_region_pileup(
            chrom,
            start0,
            end0,
            bam,
            fasta,
            active_threshold=active_threshold,
            inactive_padding=inactive_padding,
            extend_pileup_by=extend_pileup_by,
            augment_factor_vcfspec_range=augment_factor_vcfspec_range,
            augment_factor_repeat_area=augment_factor_repeat_area,
            aligner=aligner,
        )
        logging.debug(f'Initiating pileup for single sample finished: {sampleid}')
        pileup_dict[sampleid] = active_region_pileup

    # create MultisamplePileup object and postprocess
    active_region_mspileup = libpileup.MultisamplePileup(pileup_dict)
    augment_margins_multisample(
        active_region_mspileup, 
        inactive_padding=inactive_padding, 
        extend_pileup_by=extend_pileup_by, 
        factor_vcfspec_range=augment_factor_vcfspec_range, 
        factor_repeat_area=augment_factor_repeat_area, 
        aligner=aligner,
        allele_portion_threshold=allele_portion_threshold,
        vcfspec_concat_dist=vcfspec_concat_dist,
        raise_with_tie=raise_with_tie,
        max_pileup_width=max_pileup_width,
    )
    active_region_mspileup.set_row_spec_groups()
    active_region_mspileup.save_superseq_alignments()
    active_region_mspileup.set_vcfspecs(allele_portion_threshold=allele_portion_threshold, concat_dist_le=vcfspec_concat_dist)

    # result values
    if any(
        len(pileup.get_active_positions()) > 0
        for pileup in active_region_mspileup.pileup_dict.values()
    ):
        active_range = active_region_mspileup.range0
    else:
        active_range = None

    return active_range, active_region_mspileup



############################
# inactive padding related #
############################

def get_secure_margin_from_vcfspec(contig_vcfspec, factor):
    margin_start0_cadidates = list()
    margin_end0_cadidates = list()
    for component_vcfspec in contig_vcfspec.concat_components[0]:
        leftmost = component_vcfspec.leftmost()
        rightmost = component_vcfspec.rightmost()
        padding = int(len(leftmost.REF_range0) * factor)
        margin_start0_cadidates.append(leftmost.start0 - padding)
        margin_end0_cadidates.append(rightmost.end0 + padding)
    return min(margin_start0_cadidates), max(margin_end0_cadidates)




#def secure_inactive_padding_rightward(pileup, inactive_padding, extend_pileup_by=libpileup.DEFAULT_EXTEND_PILEUP_BY, inplace=False):
#    """Returns a modified pileup object with inactive region on the right margin"""
#    # set initial_rightmost_active_pos0
#    initial_active_positions = pileup.get_active_positions()
#    if len(initial_active_positions) == 0:
#        initial_rightmost_active_pos0 = pileup.start0 - 1
#    else:
#        initial_rightmost_active_pos0 = max(initial_active_positions)
#    # initial extend
#    while True:
#        if pileup.end0 - (initial_rightmost_active_pos0 + 1) < inactive_padding:
#            pileup.extend_right(width=extend_pileup_by)
#        else:
#            break
#    # initial search loop
#    active_info = pileup.get_active_info(start0=(initial_rightmost_active_pos0 + 1))
#    found_result, farthest_inactive_pos0 = _find_inactive_run(active_info, inactive_padding)
#    # subsequent search loop
#    if not found_result:
#        while True:
#            pileup.extend_right(width=extend_pileup_by)
#            active_info = pileup.get_active_info(
#                start0=(pileup.end0 - extend_pileup_by + 1 - inactive_padding)
#            )
#            found_result, farthest_inactive_pos0 = _find_inactive_run(active_info, inactive_padding)
#            if found_result:
#                break
#    # clip pileup
#    if inplace:
#        pileup.subset(pileup.start0, farthest_inactive_pos0 + 1, inplace=True)
#    else:
#        pileup_clipped = pileup.subset(pileup.start0, farthest_inactive_pos0 + 1, inplace=False)
#        return pileup_clipped
#
#
#def secure_inactive_padding_leftward(pileup, inactive_padding, extend_pileup_by=libpileup.DEFAULT_EXTEND_PILEUP_BY, inplace=False):
#    """Returns a modified pileup object with inactive region on the left margin"""
#    # set initial_leftmost_active_pos0
#    initial_active_positions = pileup.get_active_positions()
#    if len(initial_active_positions) == 0:
#        initial_leftmost_active_pos0 = pileup.end0
#    else:
#        initial_leftmost_active_pos0 = min(initial_active_positions)
#    # initial extend
#    while True:
#        if initial_leftmost_active_pos0 - pileup.start0 < inactive_padding:
#            pileup.extend_left(width=extend_pileup_by)
#        else:
#            break
#    # initial search loop
#    active_info = pileup.get_active_info(end0=initial_leftmost_active_pos0)[::-1]
#    found_result, farthest_pos0 = _find_inactive_run(active_info, inactive_padding)
#    # subsequent search loop
#    if not found_result:
#        while True:
#            pileup.extend_left(width=extend_pileup_by)
#            active_info = pileup.get_active_info(
#                end0=(pileup.start0 + extend_pileup_by - 1 + inactive_padding)
#            )[::-1]
#            found_result, farthest_pos0 = _find_inactive_run(active_info, inactive_padding)
#            if found_result:
#                break
#    # clip pileup
#    if inplace:
#        pileup.subset(farthest_pos0, pileup.end0, inplace=True)
#    else:
#        pileup_clipped = pileup.subset(farthest_pos0, pileup.end0, inplace=False)
#        return pileup_clipped
#
#
#def _find_inactive_run(active_info, inactive_padding):
#    found_result = False
#    farthest_inactive_pos0 = None
#    for start_idx in range(0, len(active_info) - inactive_padding + 1):
#        if not any(active_info[start_idx:(start_idx + inactive_padding)]):
#            found_result = True
#            farthest_inactive_pos0 = active_info.index[start_idx + inactive_padding -1]
#            break
#    return found_result, farthest_inactive_pos0


############################
# get_active_region_pileup #
############################

def secure_margin_with_vcfspecs(active_region_pileup, factor_vcfspec_range, factor_repeat_area):
    candidate_start0s = set()
    candidate_end0s = set()
    for concat_vcfspec in itertools.chain.from_iterable(active_region_pileup.vcfspecs):
        for component_vcfspec in concat_vcfspec.concat_components[0]:
            leftmost = component_vcfspec.leftmost()
            rightmost = component_vcfspec.rightmost()
            # padding with (vcfspec REF_range length) * factor
            padding = int(len(component_vcfspec.REF_range0) * factor_vcfspec_range)
            candidate_start0s.add(leftmost.start0 - padding)
            candidate_end0s.add(rightmost.end0 + padding)
            # padding considering repetitive reference region
            #repeat_area_start0 = leftmost.start0
            #repeat_area_end0 = rightmost.end0
            #padding = int((repeat_area_end0 - repeat_area_start0) * factor_repeat_area)
            #candidate_start0s.add(repeat_area_start0 - padding)
            #candidate_end0s.add(repeat_area_end0 + padding)

    if len(candidate_start0s) == 0:
        left_edited = False
    else:
        desired_start0 = min(candidate_start0s)
        start_margin = desired_start0 - active_region_pileup.start0
        if start_margin < 0:
            left_edited = True
            active_region_pileup.extend_left(-start_margin)
        else:
            left_edited = False

    if len(candidate_end0s) == 0:
        right_edited = False
    else:
        desired_end0 = max(candidate_end0s)
        end_margin = active_region_pileup.end0 - desired_end0
        if end_margin < 0:
            right_edited = True
            active_region_pileup.extend_right(-end_margin)
        else:
            right_edited = False

    edited = left_edited or right_edited

    return edited


def augment_margins(active_region_pileup, inactive_padding, extend_pileup_by, factor_vcfspec_range, factor_repeat_area, aligner, allele_portion_threshold, vcfspec_concat_dist, raise_with_tie, max_pileup_width):
    edited = False
    while True:

        if len(active_region_pileup.range0) >= max_pileup_width:
            logging.debug('break')
            break

        if active_region_pileup.check_inactive_margin_right(inactive_padding):
            right_inactive_margin_edited = False
        else:
            right_inactive_margin_edited = True
            edited = True
            active_region_pileup.secure_inactive_padding_rightward(inactive_padding, extend_pileup_by=extend_pileup_by)

        logging.debug(f'After secure_inactive_padding_rightward {active_region_pileup.range0} {len(active_region_pileup.range0)}')

        if active_region_pileup.check_inactive_margin_left(inactive_padding):
            left_inactive_margin_edited = False
        else:
            left_inactive_margin_edited = True
            edited = True
            active_region_pileup.secure_inactive_padding_leftward(inactive_padding, extend_pileup_by=extend_pileup_by)

        logging.debug(f'After secure_inactive_padding_leftward {active_region_pileup.range0} {len(active_region_pileup.range0)}')

        # preparation for secure_margin_with_vcfspecs
        active_region_pileup.set_row_specs()
        active_region_pileup.set_row_spec_groups()
        active_region_pileup.save_superseq_alignments(aligner=aligner, raise_with_tie=raise_with_tie)
        active_region_pileup.set_vcfspecs(allele_portion_threshold=allele_portion_threshold, concat_dist_le=vcfspec_concat_dist)

        vcfspec_margin_edited = secure_margin_with_vcfspecs(active_region_pileup, factor_vcfspec_range, factor_repeat_area)

        logging.debug(f'After secure_margin_with_vcfspecs {active_region_pileup.range0} {len(active_region_pileup.range0)}')

        if vcfspec_margin_edited:
            edited = True

        if not vcfspec_margin_edited:
            logging.debug('break')
            break


    return edited


def augment_margins_multisample(active_region_mspileup, inactive_padding, extend_pileup_by, factor_vcfspec_range, factor_repeat_area, aligner, allele_portion_threshold, vcfspec_concat_dist, raise_with_tie, max_pileup_width):
    while True:
        active_region_mspileup.equalize_margins()
        augmented = dict()
        for sampleid, pileup in active_region_mspileup.pileup_dict.items():
            logging.debug(f'augment_margins_multisample: {sampleid}')
            augmented[sampleid] = augment_margins(pileup, inactive_padding, extend_pileup_by, factor_vcfspec_range, factor_repeat_area, aligner, allele_portion_threshold, vcfspec_concat_dist, raise_with_tie, max_pileup_width)
        if not any(augmented.values()):
            logging.debug(f'augment_margins_multisample: break')
            break




# get aligned reads

#def get_realigned_reads(active_region_pileup):
#    """Must be run after 'save_superseq_alignments' method"""
#    return get_realigned_reads_helper(
#        active_region_pileup, 
#        active_region_pileup.row_spec_groups, 
#        active_region_pileup._alignment_cache,
#    )


def get_realigned_reads_multisample(active_region_mspileup):
    """Must be run after 'save_superseq_alignments' method"""
    realigned_reads = dict()
    for sampleid, pileup in active_region_mspileup.pileup_dict.items():
        realigned_reads_singlesample = get_realigned_reads_helper(
            pileup, 
            active_region_mspileup.row_spec_groups[sampleid], 
            active_region_mspileup._alignment_cache,
        )
        realigned_reads[sampleid] = realigned_reads_singlesample
            
    return realigned_reads


def get_realigned_reads_helper(subseq_pileup, row_spec_groups, superseq_alignment_cache):
    """Must be run after 'save_superseq_alignments' method"""
    active_range = subseq_pileup.range0
    realigned_reads = dict()
#    for superseq_rowid, groupinfo in sorted(
#        active_region_pileup.row_spec_groups.items(),
#        key=(lambda x: x[1]['subseq_hits']), 
#        reverse=True,
#    ):
    # sorting is not necessary

    for superseq_rowid, groupinfo in row_spec_groups.items():
        #realigned_read_superseq = realign_superseq(superseq_rowid, active_region_pileup, active_range)
        #realigned_reads[superseq_rowid] = realigned_read_superseq
            # (221019) This is commented out because superseq itself is
            # included in subseq list.
        superseq_aln = superseq_alignment_cache[superseq_rowid]
        for subseq_rowid in groupinfo['subseq_rowids']:
            #if subseq_rowid in realigned_reads.keys():
            #    continue
                # (221019) This is commented out because subseqs with multiple superseq hits are
                # incoporated into only one superseq group.
            realigned_read_subseq = realign_subseq(subseq_rowid, subseq_pileup, active_range, superseq_aln)
            realigned_reads[subseq_rowid] = realigned_read_subseq
            
    return realigned_reads


# get naive vcfspecs

#def get_vcfspecs_from_pileup(active_region_pileup, allele_portion_threshold=None, concat_dist_le=None):
#    """Must be run after row_spec_groups is set and save_superseq_alignments have been run"""
#    if allele_portion_threshold is None:
#        allele_portion_threshold = DEFAULT_VCFSPEC_EXTRACTION_THRESHOLD
#
#    contig_vcfspecs = set()
#
#    chrom = active_region_pileup.chrom
#    active_region_start0 = active_region_pileup.start0
#    fasta = active_region_pileup.fasta
#    subseq_hits_sum = sum(x['subseq_hits'] for x in active_region_pileup.row_spec_groups.values())
#
#    for superseq_rowid, groupinfo in active_region_pileup.row_spec_groups.items():
#        allele_portion = groupinfo['subseq_hits'] / subseq_hits_sum
#        if allele_portion < allele_portion_threshold:
#            continue
#
#        superseq_row_spec = active_region_pileup.row_specs[superseq_rowid]
#        superseq_alignment = active_region_pileup._alignment_cache[superseq_rowid]
#
#        superseq_vcfspecs = superseq_to_vcfspecs(superseq_row_spec, superseq_alignment, chrom, active_region_start0, fasta, concat_dist_le)
#        active_region_pileup._vcfspec_cache[superseq_rowid] = superseq_vcfspecs
#        if len(superseq_vcfspecs) > 0:
#            contig_vcfspecs.add(superseq_vcfspecs)
#            
#    return contig_vcfspecs


def get_vcfspecs_from_pileup_multisample(active_region_mspileup, allele_portion_threshold=None, concat_dist_le=None):
    """Must be run after 'active_region_mspileup' has executed 
    'divide_row_spec_groups' and 'get_realigned_reads' methods.
        'divide_row_spec_groups' is needed for setting 'row_spec_groups' for
            each single sample Pileup.
        'get_realigned_reads' is needed for initiating '_alignment_cache' attribute.
    """
    if allele_portion_threshold is None:
        allele_portion_threshold = DEFAULT_VCFSPEC_EXTRACTION_THRESHOLD

    contig_vcfspecs = set()
    processed_superseq_ids = set()

    chrom = active_region_mspileup.chrom
    active_region_start0 = active_region_mspileup.start0
    fasta = active_region_mspileup.fasta

    for sampleid, pileup in active_region_mspileup.pileup_dict.items():
        row_spec_groups = active_region_mspileup.row_spec_groups[sampleid]
        subseq_hits_sum = sum(x['subseq_hits'] for x in row_spec_groups.values())
        for superseq_rowid, groupinfo in row_spec_groups.items():
            # first filter
            if superseq_rowid in processed_superseq_ids:
                continue
            # second filter
            allele_portion = groupinfo['subseq_hits'] / subseq_hits_sum
            if allele_portion < allele_portion_threshold:
                continue

            processed_superseq_ids.add(superseq_rowid)

            superseq_row_spec = active_region_mspileup.superseq_row_specs[superseq_rowid]
            superseq_alignment = active_region_mspileup._alignment_cache[superseq_rowid]

            read_uid = superseq_rowid[1]
            if read_uid in pileup._vcfspec_cache:
                superseq_vcfspecs = pileup._vcfspec_cache[read_uid]
            else:
                superseq_vcfspecs = superseq_to_vcfspecs(superseq_row_spec, superseq_alignment, chrom, active_region_start0, fasta, concat_dist_le)
            if len(superseq_vcfspecs) > 0:
                contig_vcfspecs.add(superseq_vcfspecs)
            
    return contig_vcfspecs


#def superseq_to_vcfspecs(superseq_row_spec, superseq_alignment, chrom, active_region_start0, fasta, concat_dist_le):
#    if concat_dist_le is None:
#        concat_dist_le = DEFAULT_VCFSPEC_CONCAT_DISTANCE
#
#    lstrip_query_gaps = not superseq_row_spec["left_filled"]
#    rstrip_query_gaps = not superseq_row_spec["right_filled"]
#    superseq_vcfspecs = alignhandler.alignment_to_vcfspec(
#        alignment=superseq_alignment, 
#        target_start0=active_region_start0, 
#        chrom=chrom, 
#        fasta=fasta,
#        lstrip_query_gaps=lstrip_query_gaps,
#        rstrip_query_gaps=rstrip_query_gaps,
#    )
#    if len(superseq_vcfspecs) > 0:
#        superseq_vcfspecs = libvcfspec.concat_list(superseq_vcfspecs, distance_le=concat_dist_le)
#    return superseq_vcfspecs


# process naive vcfspecs into final product

#def process_contig_vcfspecs(contig_vcfspecs):
#    vcfspecs_dict = {
#        x.get_id(): x for x in 
#        itertools.chain.from_iterable(contig_vcfspecs)
#    }
#
#    grlist = list()
#    for idx, row in enumerate(contig_vcfspecs):
#        gr = pr.concat(x.to_gr() for x in row)
#        gr = gr.assign('row_index', lambda df: pd.Series(itertools.repeat(idx, df.shape[0])))
#        grlist.append(gr)
#
#    result = list()
#    clustered_gr = pr.concat(grlist).cluster()
#    for cluster_idx in set(clustered_gr.Cluster):
#        concatenated_vcfspecs = list()
#        gr_subset = clustered_gr[clustered_gr.Cluster == cluster_idx]
#        for row_index in set(gr_subset.row_index):
#            current_vcfspec_list = [
#                vcfspecs_dict[vcfspec_id]
#                for vcfspec_id in gr_subset[gr_subset.row_index == row_index].Id
#            ]
#            vcfspec_concat = libvcfspec.concat_list(current_vcfspec_list, distance_le=None)[0]
#            concatenated_vcfspecs.append(vcfspec_concat)
#        merged_vcfspec = functools.reduce(libvcfspec.merge, concatenated_vcfspecs)
#        result.append(merged_vcfspec)
#
#    result.sort(key=(lambda x: x.pos))
#
#    return result


# write realigned bam

def write_realigned_bam(bam_path, realigned_reads, active_region_pileup, padding_width=0):
    """Args:
        realigned_reads: dict (keys: ReadUID, values: pysam.AlignedSegment)
    """
    written_reads_range = range(
        active_region_pileup.start0 - padding_width, 
        active_region_pileup.end0 + padding_width, 
    )
    # write reads
    hdr = pysam.AlignmentHeader.from_references(
        reference_names=active_region_pileup.fasta.references,
        reference_lengths=active_region_pileup.fasta.lengths,
    )
    with pysam.AlignmentFile(bam_path, mode='wb', header=hdr) as in_bam:
        # write non-realigned reads
        for read in readhandler.get_fetch(
            active_region_pileup.bam,
            active_region_pileup.chrom,
            written_reads_range.start,
            written_reads_range.stop,
        ):
            uid = readhandler.get_uid(read)
            if uid not in realigned_reads.keys():
                in_bam.write(read)
        # write realigned reads
        for read in realigned_reads.values():
            in_bam.write(read)
    # sort and index
    bameditor.sort_and_index(bam_path)


def write_realigned_bam_multisample(bam_paths, realigned_reads_multisample, active_region_pileup_multisample, padding_width=0):
    for sampleid, path in bam_paths.items():
        write_realigned_bam(
            path,
            realigned_reads_multisample[sampleid],
            active_region_pileup_multisample.pileup_dict[sampleid],
            padding_width=padding_width,
        )


###################################
# helpers for get_realigned_reads #
###################################

#def save_superseq_alignments(pileup, aligner, raise_with_tie):
#    if aligner is None:
#        aligner = ALIGNER_MAIN
#
#    pileup._alignment_cache = dict()
#    ref_seq = pileup.get_ref_seq()
#    ref_seq_reversed = ref_seq[::-1]
#
#    for superseq_rowid in pileup.row_spec_groups.keys():
#        row_spec = pileup.row_specs[superseq_rowid]
#        pileup._alignment_cache[superseq_rowid] = align_row_spec_to_ref(row_spec, ref_seq, ref_seq_reversed, aligner, raise_with_tie)


def save_superseq_alignments_multisample(mspileup):
    mspileup._alignment_cache = dict()
    for superseq_rowid in mspileup.superseq_row_specs.keys():
        sampleid, read_uid = superseq_rowid
        superseq_aln = mspileup.pileup_dict[sampleid]._alignment_cache[read_uid]
        mspileup._alignment_cache[superseq_rowid] = superseq_aln


#def realign_superseq(row_id, active_region_pileup, active_range):
#    active_region_alignment = active_region_pileup._alignment_cache[row_id]
#    realigned_read = align_row_helper(row_id, active_region_pileup, active_range, active_region_alignment)
#
#    return realigned_read


def realign_subseq(row_id, active_region_pileup, active_range, superseq_alignment):
    row_spec = active_region_pileup.row_specs[row_id]
    active_region_alignment = align_subseq_to_ref(row_spec, superseq_alignment)
    try:
        realigned_read = align_row_helper(row_id, active_region_pileup, active_range, active_region_alignment)
    except RealignmentCigarError as exc:
        raise Exception(f'Subseq realignment failed.\nrow_spec: {row_spec}\n') from exc

    return realigned_read


def align_row_helper(row_id, active_region_pileup, active_range, active_region_alignment):
    # set params
    original_read = active_region_pileup.read_store[row_id]
    original_cigar_split = alignhandler.split_cigar(original_read.cigartuples, original_read.reference_start, active_range)
    empty_before_active_region = (len(original_cigar_split[0]) == 0)
    empty_after_active_region = (len(original_cigar_split[2]) == 0)
    # get active region cigartuples
    active_region_cigartuples, active_region_offset = alignhandler.alignment_to_cigartuples(
        active_region_alignment,
        match_as_78=False, del_as_skip=False, 
        left_ins_as_clip=empty_before_active_region, right_ins_as_clip=empty_after_active_region,
        remove_left_del=empty_before_active_region, remove_right_del=empty_after_active_region,
    )
    # merge split cigartuples
    realigned_cigartuples = functools.reduce(
        merge_two_cigartuples, 
        [original_cigar_split[0], active_region_cigartuples, original_cigar_split[2]],
    )
    
    # get new read reference_start
    if empty_before_active_region:
        new_reference_start = active_range.start + active_region_offset
    else:
        new_reference_start = original_read.reference_start
    # make new read 
    realigned_read = original_read.__copy__()
    realigned_read.reference_start = new_reference_start
    realigned_read.cigartuples = realigned_cigartuples
    readhandler.set_NMMD(realigned_read, active_region_pileup.fasta)

    realigned_cigar_sanity_check(original_read, realigned_read, realigned_cigartuples)

    return realigned_read


def align_subseq_to_ref(subseq_row_spec, sup_to_ref_aln):
    # get subseq index
    reverse_align = (not subseq_row_spec['left_filled']) and subseq_row_spec['right_filled']
    if reverse_align:
        sub_to_sup_index = sup_to_ref_aln.query.rfind(subseq_row_spec['seq'])
    else:
        sub_to_sup_index = sup_to_ref_aln.query.find(subseq_row_spec['seq'])
    if sub_to_sup_index == -1:
        raise Exception(
            f'"subseq" is not a subseq of "superseq":\n'
            f'subseq_row_spec: {subseq_row_spec}\n'
            f'superseq: {sup_to_ref_aln.query}'
        )
    # convert into sub-to-ref alignment
    # supaln: sup_to_ref_aln
    # subaln: sub_to_sup_aln
    active_region_walks = list()
    target_walk_subaln = range(sub_to_sup_index, sub_to_sup_index + len(subseq_row_spec['seq']))
    query_walk_subaln = range(0, len(subseq_row_spec['seq']))
    found_initial_hit = False
    for target_walk_supaln, query_walk_supaln in alignhandler.get_walks(sup_to_ref_aln, copy=False):
        if not found_initial_hit:
            if target_walk_subaln.start in query_walk_supaln:
                found_initial_hit = True
                # superseq walk
                superseq_walk_start = target_walk_subaln.start
                superseq_walk_stop = min(query_walk_supaln.stop, target_walk_subaln.stop)
                len_current_superseq_walk = superseq_walk_stop - superseq_walk_start
                # subseq walk
                subseq_walk = range(0, len_current_superseq_walk)
                # ref walk
                if len(target_walk_supaln) == 0:
                    ref_walk = target_walk_supaln
                else:
                    # During initial hit query_walk_supaln cannot be zero-length
                    ref_walk_start = target_walk_supaln.start + (superseq_walk_start - query_walk_supaln.start)
                    ref_walk = range(ref_walk_start, ref_walk_start + len_current_superseq_walk)

                active_region_walks.append((ref_walk, subseq_walk))

                if query_walk_supaln.stop >= target_walk_subaln.stop:
                    break
        else:
            # superseq walk
            superseq_walk_start = query_walk_supaln.start
            superseq_walk_stop = min(query_walk_supaln.stop, target_walk_subaln.stop)
            len_current_superseq_walk = superseq_walk_stop - superseq_walk_start
            # subseq walk
            subseq_walk_start = active_region_walks[-1][1].stop
            subseq_walk = range(subseq_walk_start, subseq_walk_start + len_current_superseq_walk)
            # ref walk
            if len(target_walk_supaln) == 0 or len(query_walk_supaln) == 0:
                ref_walk = target_walk_supaln
            else:
                ref_walk_start = target_walk_supaln.start
                ref_walk = range(ref_walk_start, ref_walk_start + len_current_superseq_walk)

            active_region_walks.append((ref_walk, subseq_walk))

            if query_walk_supaln.stop >= target_walk_subaln.stop:
                break
    # pad with gap: right side
    last_ref_walk = active_region_walks[-1][0]
    last_subseq_walk = active_region_walks[-1][1]
    if last_ref_walk.stop < len(sup_to_ref_aln.target):
        added_ref_walk = range(last_ref_walk.stop, len(sup_to_ref_aln.target))
        added_subseq_walk = range(last_subseq_walk.stop, last_subseq_walk.stop)
        active_region_walks.append((added_ref_walk, added_subseq_walk))
    # pad with gap: left
    first_ref_walk = active_region_walks[0][0]
    first_subseq_walk = active_region_walks[0][1]
    if first_ref_walk.start > 0:
        added_ref_walk = range(0, first_ref_walk.start)
        added_subseq_walk = range(0, 0)
        active_region_walks.insert(0, (added_ref_walk, added_subseq_walk))
    # make result alignment
    return Bio.Align.PairwiseAlignment(
        sup_to_ref_aln.target,
        subseq_row_spec['seq'],
        alignhandler.walks_to_path(active_region_walks),
        0
    )


#######################################
# amenders: marginal insdel amendment #
#######################################

#def amender_left_right(alignment):
#    return alignhandler.amend_outer_insdel_right(
#        alignhandler.amend_outer_insdel_left(alignment)
#    )
#
#
#def amender_left(alignment):
#    return alignhandler.amend_outer_insdel_left(alignment)
#
#
#def amender_right(alignment):
#    return alignhandler.amend_outer_insdel_right(alignment)
#
#
#def amender_none(alignment):
#    return alignment
#
#
#def get_amender(empty_before_active_region, empty_after_active_region):
#    if empty_before_active_region:
#        if empty_after_active_region:
#            return amender_left_right
#        else:
#            return amender_left
#    else:
#        if empty_after_active_region:
#            return amender_right
#        else:
#            return amender_none


#######################################


def align_row_spec_to_ref(row_spec, ref_seq, ref_seq_reversed, aligner, raise_with_tie):
    # set params
    if aligner is None:
        aligner = ALIGNER_MAIN
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
    # recover reversed alignment
    if reverse_align:
        alns = [alignhandler.reverse_alignment(x) for x in alns]

    if len(alns) == 1:
        aln = alns[0]
        aln = alignhandler.amend_outer_insdel_both(aln)
    else:
        #alns = [alignhandler.amend_outer_insdel_both(x) for x in alns]
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


def merge_two_cigartuples(cigartuples_left, cigartuples_right):
    if len(cigartuples_left) == 0 or len(cigartuples_right) == 0:
        return cigartuples_left + cigartuples_right
    else:
        if cigartuples_left[-1][0] == cigartuples_right[0][0]:
            return (
                cigartuples_left[:-1]
                + [(cigartuples_left[-1][0], cigartuples_left[-1][1] + cigartuples_right[0][1])]
                + cigartuples_right[1:]
            )
        else:
            return cigartuples_left + cigartuples_right


def realigned_cigar_sanity_check(original_read, realigned_read, realigned_cigartuples):
    if len(original_read.query_sequence) != alignhandler.get_query_length(realigned_cigartuples):
        raise RealignmentCigarError(
            f'Read query sequence length and cigar length differs.\n'
            f'Original read:\n{original_read.to_string()}\n'
            f'Realigned read:\n{realigned_read.to_string()}\n'
        )

    

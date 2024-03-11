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

import handygenome.deco as deco
import handygenome.workflow as workflow
import handygenome.logutils as logutils
import handygenome.variant.vcfspec as libvcfspec
from handygenome.read.pileup import PileupBase
import handygenome.align.alignhandler as alignhandler


#logging.basicConfig(level=logging.INFO)

#DEFAULT_ACTIVE_THRESHOLD = 0.05
#DEFAULT_INACTIVE_PADDING = 10
#DEFAULT_VCFSPEC_EXTRACTION_THRESHOLD = 0.05
#DEFAULT_VCFSPEC_CONCAT_DISTANCE = None
#DEFAULT_MAX_PILEUP_WIDTH = 100  # Must be equal to or less than this value
#DEFAULT_EXTEND_PILEUP_BY = 20
#DEFAULT_AUGMENT_FACTOR_VCFSPEC_RANGE = 1.5
#DEFAULT_MQ_LIMIT = 20
#DEFAULT_DEPTH_LIMIT = 1

#DEFAULT_RPILEUPSERIES_MAX_WIDTH = 5000

#DEFAULT_AUGMENT_FACTOR_REPEAT_AREA = 1


DEFAULT_ALIGNER = alignhandler.ALIGNER_EQUAL_MM_GAP

DEFAULT_PARAMS = {
    'active_threshold_vaf': 0.05,
    'active_threshold_depth': 8,
    'inactive_padding': 15,  # This value is allowed
    'allele_portion_threshold': 0.05, 
    'concat_dist_le': None,  # Now unused
    'extend_pileup_by': 20,
    'vcfspec_range_factor': 1.5,
    'start0_limit': -np.inf,  # This value is allowed
    'end0_limit': np.inf,  # This value is allowed
    'MQ_limit': 20,  # This value is allowed
    'depth_limit': 1,  # This value is allowed
    'max_pileup_width': 100,  # This value is allowed
    'max_series_width': 5000,  # This value is allowed
}


class LoggingBase:
    def log_debug(self, msg):
        if self.verbose:
            logutils.log(msg, level='debug', add_locstring=True, locstring_mode='last')


class RealignerPileupBase(PileupBase, LoggingBase):
    row_spec_exluded_vals = (PileupBase.DEL_VALUE, PileupBase.EMPTY_VALUE)

#    @staticmethod
#    def get_logger(verbose):
#        formatter = logging.Formatter(
#            fmt='[%(asctime)s.%(msecs)03d RealignerPileup] line %(lineno)d: %(message)s', 
#            datefmt='%Z %Y-%m-%d %H:%M:%S'
#        )
#        return workflow.get_logger(
#            name=str(uuid.uuid4()),
#            level=('debug' if verbose else 'info'),
#            formatter=formatter,
#        )

    def _active_info_generator(self, start0, end0, reverse=False):
        self._coord_arg_sanitycheck(start0, end0)

        ref_seq = self.fasta.fetch(self.chrom, start0, end0)
        pos0_range = range(start0, end0)
        if reverse:
            ref_seq = ref_seq[::-1]
            pos0_range = pos0_range[::-1]

        for ref_base, pos0 in zip(ref_seq, pos0_range):
            counter = self.get_allele_counter(pos0)
            depth = sum(counter.values())
            if depth < self.params['active_threshold_depth']:
                is_active = False
            else:
                non_ref_portion = 1 - (counter[ref_base] / depth)
                is_active = (non_ref_portion >= self.params['active_threshold_vaf'])

            yield is_active

    def _set_active_info(self):
        if self.seq_df.shape[1] == 0:
            active_info = pd.Series([], dtype=bool)
        else:
            if self.seq_df.shape[0] > 0:
                active_info = pd.Series(
                    self._active_info_generator(self.start0, self.end0, reverse=False),
                    index=self.seq_df.columns,
                )
            else:
                active_info = pd.Series(
                    [False] * self.seq_df.shape[1],
                    index=self.seq_df.columns,
                )

        self.active_info = active_info

    # row specs related ones
    def _make_row_spec(self, row):
        result = dict()
        result['seq'] = "".join(
            x for x in row 
            if x not in self.row_spec_exluded_vals
        )
        result['id'] = row.name

        row_is_not_empty = (row != self.EMPTY_VALUE)
        non_empty_index = row.index[row_is_not_empty]
        result['span_start0'] = non_empty_index.min()
        result['span_end0'] = non_empty_index.max() + 1
        result['span_length'] = result['span_end0'] - result['span_start0']
        result['left_filled'] = row_is_not_empty.iloc[0]
        result['right_filled'] = row_is_not_empty.iloc[-1]

        read = self.read_getter(row.name)
        result['cliplen'] = read.get_cigar_stats()[0][4]

        return result

    @staticmethod
    def row_spec_matcher_old(query, target):
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
    def row_spec_matcher(query, target):
        if query['left_filled'] and query['right_filled'] and target['left_filled'] and target['right_filled']:
            return query['seq'] == target['seq']
        else:
            return query['seq'] in target['seq']

    @staticmethod
    def tiebreak_key(superseq_row_spec, groupinfo):
        #return (len(superseq_row_spec['seq']), superseq_row_spec['span_length'], groupinfo['subseq_hits'])
        return (superseq_row_spec['span_length'], -superseq_row_spec['cliplen'], groupinfo['subseq_hits'])

    @classmethod
    def superseq_candidates_tiebreak(cls, superseq_candidates, superseq_row_spec_getter, row_spec_groups):
        if len(superseq_candidates) == 1:
            superseq_rowid = superseq_candidates[0]
        elif len(superseq_candidates) > 1:
            superseq_rowid = max(
                superseq_candidates, 
                key=(lambda x: cls.tiebreak_key(superseq_row_spec_getter(x), row_spec_groups[x])),
            )

        return superseq_rowid

    def decorate_groupinfo(self, row_spec_groups, read_getter):
        subseq_hits_sum = sum(x['subseq_hits'] for x in row_spec_groups.values())
        for groupinfo in row_spec_groups.values():
            if subseq_hits_sum == 0:
                groupinfo['subseq_hits_portion'] = np.nan
            else:
                groupinfo['subseq_hits_portion'] = groupinfo['subseq_hits'] / subseq_hits_sum

            MQs = [
                read_getter(key).mapping_quality
                for key in groupinfo['subseq_rowids']
            ]
            if len(MQs) == 0:
                groupinfo['mean_MQ'] = np.nan
            else:
                groupinfo['mean_MQ'] = np.mean(MQs)

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
        return Bio.Align.Alignment(
            sequences=[sup_to_ref_aln.target, subseq_row_spec['seq']],
            coordinates=alignhandler.walks_to_coordinates(active_region_walks),
        )
#        return Bio.Align.PairwiseAlignment(
#            sup_to_ref_aln.target,
#            subseq_row_spec['seq'],
#            alignhandler.walks_to_path(active_region_walks),
#            0
#        )

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

        # remove those with N in ALT
        superseq_vcfspecs = [x for x in superseq_vcfspecs if 'N' not in x.alts[0]]

        if len(superseq_vcfspecs) == 0:
            # convert identity variant to None
            superseq_vcfspec = None
        else:
            # concat
            superseq_vcfspec = libvcfspec.concat_list(superseq_vcfspecs, distance_le=None)[0]

        return superseq_vcfspec

    def iter_contig_vcfspecs(self, subseq_portion_threshold=0, MQ_threshold=40):
        for row_id, contig_vcfspec in self.contig_vcfspecs.items():
            groupinfo = self.row_spec_groups[row_id]
            if (
                (groupinfo['subseq_hits_portion'] >= subseq_portion_threshold) and
                (groupinfo['mean_MQ'] >= MQ_threshold)
            ):
                yield row_id, contig_vcfspec

#    def get_result_vcfspecs(self, as_components=False, merge=True, subseq_portion_threshold=None, MQ_threshold=40): 
#        if subseq_portion_threshold is None:
#            subseq_portion_threshold = self.allele_portion_threshold
#
#        if as_components:
#            result = list()
#            for row_id, contig_vcfspec in self.iter_contig_vcfspecs(subseq_portion_threshold=subseq_portion_threshold, MQ_threshold=MQ_threshold):
#                if contig_vcfspec is None:
#                    continue
#
#                if len(contig_vcfspec.components[0]) == 0:
#                    result.append(contig_vcfspec)
#                else:
#                    result.extend(contig_vcfspec.components[0])
#        else:
#            result = list(
#                x[1] for x in self.iter_contig_vcfspecs(subseq_portion_threshold, MQ_threshold)
#                if x[1] is not None
#            )
#
#        if merge:
#            if len(result) > 0:
#                result = functools.reduce(lambda x, y: libvcfspec.merge(x, y), result)
#
#        return result

    def show_contig_vcfspecs(self):
        contig_vcfspecs_sorted = sorted(
            self.iter_contig_vcfspecs(subseq_portion_threshold=0, MQ_threshold=0),
            key=(lambda x: self.row_spec_groups[x[0]]['subseq_hits']),
            reverse=True,
        )
        if len(contig_vcfspecs_sorted) == 0:
            print('empty')
        else:
            for row_id, contig_vcfspec in contig_vcfspecs_sorted:
                subseq_hits = self.row_spec_groups[row_id]['subseq_hits']
                subseq_hits_portion = self.row_spec_groups[row_id]['subseq_hits_portion']
                print(contig_vcfspec)
                print(f'subseq_hits: {subseq_hits}; portion: {subseq_hits_portion}')
                if contig_vcfspec is not None:
                    print('components:')
                    for component in contig_vcfspec.components[1]:
                        print(f'\t{component}')
                print('\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n')

    # for debugging
    def _show_row_spec_alignments_helper(self, row_spec_groups, row_spec_getter, aln_getter, skip_zero_hits, show_subseqs):
        subseq_hits_sum = sum(x['subseq_hits'] for x in row_spec_groups.values())
        for superseq_id, groupinfo in sorted(
            row_spec_groups.items(),
            key=(lambda x: x[1]['subseq_hits']),
            reverse=True,
        ):
            if skip_zero_hits:
                if groupinfo['subseq_hits'] == 0:
                    continue

            # set params
            superseq_row_spec = row_spec_getter(superseq_id)
            superseq_aln = aln_getter(superseq_id)
            superseq_vcfspec = self.superseq_to_vcfspecs(superseq_row_spec, superseq_aln)
            # show superseq row spec and alignment
            print('superseq id:', superseq_id)
            print('row_spec:', superseq_row_spec)
            print(f'alignment:\n{superseq_aln}')
            # show vcfspec
            print(f'contig vcfspec:\n\t{superseq_vcfspec}')
            if superseq_vcfspec is not None and len(superseq_vcfspec.components[1]) > 0:
                print(f'concat components:')
                for component in superseq_vcfspec.components[1]:
                    print('\t', component)
            # show other groupinfo attrs
            for key, val in sorted(
                groupinfo.items(),
                key=(lambda x: x[0] in ('subseq_hits_portion', 'subseq_hits')),
                reverse=True,
            ):
                if key == 'subseq_rowids':
                    if show_subseqs:
                        print('subseq row ids:')
                        for rowid in groupinfo['subseq_rowids']:
                            print('\t', rowid)
                else:
                    print(key, val)

            print('\n@@@@@@@@@@@@@@@@@@@@@\n')

    def show_alns(self, **kwargs):
        self.show_row_spec_alignments(**kwargs)


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
            self.left_low_MQ = self.pileup.MQ.loc[self.current_start0] < self.pileup.params['MQ_limit']
            self.left_low_depth = self.pileup.get_depth(self.current_start0) < self.pileup.params['depth_limit']

        if right:
            self.right_is_active = self.pileup.active_info.loc[self.current_end0 - 1]
            self.right_low_MQ = self.pileup.MQ.loc[self.current_end0 - 1] < self.pileup.params['MQ_limit']
            self.right_low_depth = self.pileup.get_depth(self.current_end0 - 1) < self.pileup.params['depth_limit']

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
        if current_start0 <= self.pileup.params['start0_limit']:
            return

        while True:
            self.pileup.extend_left(self.pileup.params['extend_pileup_by'])
            for _ in range(self.pileup.params['extend_pileup_by']):
                current_start0 -= 1
                #left_is_active = self.pileup.active_info.loc[current_start0]
                touched_left_limit = current_start0 <= self.pileup.params['start0_limit']
                yield current_start0, touched_left_limit
                if touched_left_limit:
                    return

    def _base_generator_right(self):
        """Stops only when current_end0 touches end0_limit"""
        current_end0 = self.pileup.end0
        if current_end0 >= self.pileup.params['end0_limit']:
            return

        while True:
            self.pileup.extend_right(self.pileup.params['extend_pileup_by'])
            for _ in range(self.pileup.params['extend_pileup_by']):
                current_end0 += 1
                #right_is_active = self.pileup.active_info.loc[current_end0]
                touched_right_limit = current_end0 >= self.pileup.params['end0_limit']
                yield current_end0, touched_right_limit
                if touched_right_limit:
                    return


class AlignmentFailureError(Exception):
    pass


class SparseInactiveRegionError(Exception):
    pass


class RealignmentCigarError(Exception):
    pass


def parse_rpileup_kwargs(**kwargs):
    if any(x not in DEFAULT_PARAMS.keys() for x in kwargs.keys()):
        raise Exception(f'Allowed kwargs keys are: {list(DEFAULT_PARAMS.keys())}')

    params = kwargs
    for k, v in DEFAULT_PARAMS.items():
        if k not in params:
            params[k] = v

    return params


def main_aligner(
    target, query, aligner, 
    reverse_align=False, 
    raise_with_tie=False, 
    #logger=None, 
    row_spec=None, 
    target_reversed=None,
):
    if reverse_align:
        if target_reversed is None:
            target = target[::-1]
        else:
            target = target_reversed
        query = query[::-1]

    try:
        alns = aligner.align(target, query)
    except Exception as exc:
        #raise Exception(f'Failed alignment:\nrow_spec: {row_spec}\nReference seq: {ref_seq}') from exc
        raise AlignmentFailureError() from exc

    try:
        aln = main_aligner_helper(alns, raise_with_tie, reverse_align)
    except TimeoutError:
        #if logger is not None:
        #    logger.debug(
        #        f'Skipping alignments tiebreaking due to timeout;\nrow_spec: {row_spec}\n'
        #        f'Num of alignments: {len(alns)}; row_spec: {row_spec}'
        #    )
        aln = alns[0]
        aln = alignhandler.amend_outer_insdel_both(aln)
        if reverse_align:
            aln = alignhandler.reverse_alignment(aln)
        
    return aln


@deco.get_deco_timeout(0.05)
def main_aligner_helper(alns, raise_with_tie, reverse_align):
    # amend outer insdel, remove duplicates
    alns = list(
        alignhandler.remove_identical_alignments(
            alignhandler.amend_outer_insdel_both(x) for x in alns
        )
    )
    # select one
    if len(alns) == 1:
        aln = alns[0]
    else:
        # tiebreaking
        aln = alignhandler.alignment_tiebreaker(alns, raise_with_failure=raise_with_tie)
    # reverse - should be done after tiebreaking
    if reverse_align:
        aln = alignhandler.reverse_alignment(aln)

    return aln


#def logging_debug(msg):
#    logutils.log(msg, level='debug', verbose_locstring=False)



import sys
import collections
import itertools
import functools
import logging
import inspect
import random
import uuid

import pysam
import Bio.Align
import numpy as np
import pandas as pd
import pyranges as pr

import importlib
top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
workflow = importlib.import_module(".".join([top_package_name, "workflow"]))
libvcfspec = importlib.import_module(".".join([top_package_name, "variant", "vcfspec"]))
libpileup = importlib.import_module(".".join([top_package_name, "read", "pileup"]))
alignhandler = importlib.import_module(".".join([top_package_name, "align", "alignhandler"]))
bameditor = importlib.import_module(".".join([top_package_name, "bameditor"]))
readhandler = importlib.import_module(".".join([top_package_name, "read", "readhandler"]))
readplus = importlib.import_module(".".join([top_package_name, "read", "readplus"]))


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

DEFAULT_RPILEUP_PARAMS = {
    'active_threshold': 0.05,
    'aligner': ALIGNER_MAIN,
    'inactive_padding': 15,
    'allele_portion_threshold': 0.05, 
    'concat_dist_le': None,  # Now unused
    'max_pileup_width': 100,
    'extend_pileup_by': 20,
    'vcfspec_range_factor': 1.5,
    'start0_limit': -np.inf,  # This value is allowed
    'end0_limit': np.inf,  # This value is allowed
    'MQ_limit': 20,  # This value is allowed
    'depth_limit': 1,  # This value is allowed
}


class AlignmentFailureError(Exception):
    pass


class RealignerPileupBase(libpileup.PileupBase):
    row_spec_exluded_vals = (
        libpileup.PileupBase.DEL_VALUE, libpileup.PileupBase.EMPTY_VALUE
    )

    @staticmethod
    def get_logger(verbose):
        formatter = logging.Formatter(
            fmt='[%(asctime)s.%(msecs)03d RealignerPileup] line %(lineno)d: %(message)s', 
            datefmt='%Z %Y-%m-%d %H:%M:%S'
        )
        return workflow.get_logger(
            name=str(uuid.uuid4()),
            level=('debug' if verbose else 'info'),
            formatter=formatter,
        )

    def _active_info_generator(self, start0, end0, reverse=False):
        self._coord_arg_sanitycheck(start0, end0)

        ref_seq = self.fasta.fetch(self.chrom, start0, end0)
        pos0_range = range(start0, end0)
        if reverse:
            ref_seq = ref_seq[::-1]
            pos0_range = pos0_range[::-1]

        for ref_base, pos0 in zip(ref_seq, pos0_range):
            counter = self.get_allele_counter(pos0)
            if len(counter) == 0:
                is_active = False
            else:
                non_ref_portion = 1 - (counter[ref_base] / sum(counter.values()))
                is_active = (non_ref_portion >= self.active_threshold)

            yield is_active

    def _set_active_info(self):
        if self.df.shape[0] > 0:
            self.active_info = pd.Series(
                self._active_info_generator(self.start0, self.end0, reverse=False),
                index=self.df.columns,
            )
        else:
            self.active_info = pd.Series(
                [False] * self.df.shape[1],
                index=self.df.columns,
            )

    # row specs related ones
    def _make_row_spec(self, row):
        result = dict()
        result['seq'] = "".join(
            self.pat_parenthesis.sub("", x) 
            for x in row 
            if x not in self.row_spec_exluded_vals
        )
        result['id'] = row.name

        row_is_not_empty = row != self.EMPTY_VALUE
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


        #groupinfo = superseq_groupinfo_getter(superseq_rowid)
        #groupinfo['subseq_rowids'].append(subseq_rowid)
        #groupinfo['subseq_hits'] += 1

#    @classmethod
#    def group_row_specs(cls, subseq_row_specs):
#        """row_spec's are grouped by whether one is a substring of another."""
#        initial_row_spec_groups = collections.OrderedDict()
#
#        for query_rowid, query in sorted(
#            subseq_row_specs.items(), 
#            key=(lambda x: x[1]['span_length']), 
#            reverse=True,
#        ):
#            superseq_candidates = list()
#            for superseq_rowid in initial_row_spec_groups.keys():
#                target = subseq_row_specs[superseq_rowid]
#                if cls.row_spec_matcher(query, target):
#                    superseq_candidates.append(superseq_rowid)
#                #if len(superseq_candidates) >= 2:
#                #    break
#                    
#            if len(superseq_candidates) == 0:
#                initial_row_spec_groups[query_rowid] = {
#                    'subseq_rowids': list(),
#                    'subseq_hits': 0,
#                }
#                initial_row_spec_groups[query_rowid]['subseq_rowids'].append(query_rowid)
#                initial_row_spec_groups[query_rowid]['subseq_hits'] += 1
#            else:
#                cls.superseq_candidates_tiebreak(superseq_candidates, initial_row_spec_groups, query_rowid, subseq_row_specs)
#
#        return initial_row_spec_groups

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
        bam, 
        chrom, start0=None, end0=None, 
        refver=None, fasta=None,

        init_df=True,

        verbose=False,
        logger=None,

        **kwargs,

#        active_threshold=DEFAULT_ACTIVE_THRESHOLD, 
#        aligner=ALIGNER_MAIN,
#        inactive_padding=DEFAULT_INACTIVE_PADDING,
#        allele_portion_threshold=DEFAULT_VCFSPEC_EXTRACTION_THRESHOLD, 
#        concat_dist_le=DEFAULT_VCFSPEC_CONCAT_DISTANCE,  # Now unused
#        max_pileup_width=DEFAULT_MAX_PILEUP_WIDTH,
#        extend_pileup_by=DEFAULT_EXTEND_PILEUP_BY,
#        vcfspec_range_factor=DEFAULT_AUGMENT_FACTOR_VCFSPEC_RANGE,
#        start0_limit=-np.inf,  # This value is allowed
#        end0_limit=np.inf,  # This value is allowed
#        MQ_limit=DEFAULT_MQ_LIMIT,  # This value is allowed
#        depth_limit=DEFAULT_DEPTH_LIMIT,  # This value is allowed
    ):
        # initiation
        libpileup.PileupBase.__init__(
            self, 
            bam=bam, 
            chrom=chrom, start0=start0, end0=end0, 
            refver=refver, fasta=fasta,
            init_df=init_df, 
            verbose=verbose, logger=logger,
        )

        # set other parameters
        for k, v in parse_rpileup_kwargs(**kwargs).items():
            setattr(self, k, v)

        if init_df:
            self._set_active_info()

#        self.aligner = aligner
#        self.inactive_padding = inactive_padding
#        self.allele_portion_threshold = allele_portion_threshold
#        self.concat_dist_le = concat_dist_le
#        self.max_pileup_width = max_pileup_width
#        self.extend_pileup_by = extend_pileup_by
#        self.vcfspec_range_factor = vcfspec_range_factor
#        self.start0_limit = start0_limit
#        self.end0_limit = end0_limit
#        self.MQ_limit = MQ_limit
#        self.depth_limit = depth_limit

        self.hit_left_margin = False
        self.hit_right_margin = False

    @classmethod
    def init_and_augment(
        cls, 
        bam, 
        chrom, start0, end0, 
        refver=None, fasta=None,
        verbose=False, 
        logger=None, 
        **kwargs,
    ):
        result = cls(
            bam=bam, 
            chrom=chrom, start0=start0, end0=end0, 
            refver=refver, fasta=fasta,

            init_df=True,

            verbose=verbose,
            logger=logger,

            **kwargs,
        )
        secresult_inactive, secresult_vcfspec = result.augment_margins()

        return result, secresult_inactive, secresult_vcfspec

    def spawn(self):
        attrs_to_copy = (
            ('bam', 'chrom', 'refver', 'fasta', 'verbose', 'logger')
            + tuple(DEFAULT_RPILEUP_PARAMS.keys())
        )
        return self._spawn_base(attrs_to_copy)
    
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
        self.logger.debug(f'Beginning "augment_margins"')

        # immediately after seeding
        if self.active_info.any():
            self.logger.debug(f'\tskipping initial search')
        else:
            self.logger.debug(f'\tbeginning initial search')

            self.search_for_active_position()

            self.logger.debug(f'\tfinished initial search')

            if not self.active_info.any():
                self.logger.debug(f'Finished "augment_margins" without secure_margins')

                result_inactive = None
                result_vcfspec = None
                return result_inactive, result_vcfspec

        # trim excessive inactive padding
        self.reduce_left()
        self.reduce_right()
        # do after-discovery
        self.logger.debug(f'\tbeginning after-discovery')
        result_inactive, result_vcfspec = self.secure_margins()
        self.logger.debug(f'\tfinished after-discovery')

        self.logger.debug(f'Finished "augment_margins"')
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
        self.logger.debug(f'Began secure_margins')

        while True:
            self.logger.debug(f'Began secure_inactive_padding')
            result_inactive = self.secure_inactive_padding()
            self.logger.debug(f'Finished secure_inactive_padding. Margins: {self.range0}')

            self.logger.debug(f'Began prepare_vcfspecs')
            self.prepare_vcfspecs()
            self.logger.debug(f'Finished prepare_vcfspecs')

            self.logger.debug(f'Began secure_vcfspec_margins')
            result_vcfspec = self.secure_vcfspec_margins()
            self.logger.debug(f'Finished secure_vcfspec_margins. Margins: {self.range0}')

            if not result_vcfspec.edited:
                break
            #if (not result_inactive.left_okay) and (not result_inactive.right_okay):
            #    break
            #else:

        self.logger.debug(f'Finished secure_margins')

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
        return gen.left_is_blocked or gen.touched_width_limit or left_okay

    @classmethod
    def _securehelper_check_right_blocked(cls, gen, right_okay):
        return gen.right_is_blocked or gen.touched_width_limit or right_okay

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
        vcfspec_margins_gr = self.get_vcfspec_margins_gr(
            subseq_portion_threshold=self.allele_portion_threshold,
            inverse=False,
            split_contig_vcfspec=True,
        )

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

    def get_vcfspec_margins_gr(self, subseq_portion_threshold, inverse=False, split_contig_vcfspec=True):
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
            vcfspec_margins_gr = pr.PyRanges()
        else:
            chroms = [self.chrom] * len(start0s)
            vcfspec_margins_gr = pr.from_dict(
                {'Chromosome': chroms, 'Start': start0s, 'End': end0s, 'Name': names}
            )
        # return; inverse if needed
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
#            return Bio.Align.PairwiseAlignment(
#                target=ref_seq,
#                query='',
#                path=[(0, 0), (len(ref_seq), 0)],
#                score=0,
#            )

        reverse_align = (not row_spec['left_filled']) and row_spec['right_filled']
        target = ref_seq
        query = row_spec['seq']
        #if reverse_align:
        #    target = ref_seq_reversed
        #    query = row_spec['seq'][::-1]
        #else:
        #    target = ref_seq
        #    query = row_spec['seq']

        return self.main_aligner(target, query, aligner=self.aligner, reverse_align=reverse_align, raise_with_tie=raise_with_tie, logger=self.logger, row_spec=row_spec, target_reversed=ref_seq_reversed)

    @classmethod
    def main_aligner(cls, target, query, aligner=ALIGNER_MAIN, reverse_align=False, raise_with_tie=False, logger=None, row_spec=None, target_reversed=None):
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
            aln = cls.main_aligner_helper(alns, raise_with_tie, reverse_align)
        except common.TimeoutError:
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

    @staticmethod
    @common.timeout(0.05)
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

        self.logger.debug(f'Num of alignments: {len(alns)}; row_spec: {row_spec}')

        # treat dirty alignments
        if len(alns) > 10000:
            aln = random.choice(alns)
            if reverse_align:
                aln = alignhandler.reverse_alignment(aln)
            #aln = alignhandler.amend_outer_insdel_both(aln)
            self.logger.debug(f'Skipped dirty alignment; row_spec: {row_spec}')
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
                except common.TimeoutError:
                    self.logger.debug(f'skipping alignments tiebreaking due to timeout;\nrow_spec: {row_spec}')
                    aln = alns[0]
                    aln = alignhandler.amend_outer_insdel_both(aln)
            
        return aln

    @common.timeout(0.05)
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
        splittable_region = self.get_vcfspec_margins_gr(
            subseq_portion_threshold=(self.allele_portion_threshold * 2),
            inverse=True,
            split_contig_vcfspec=False,
        )

        if splittable_region.empty:
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
        active_info_gr = self.get_active_info_gr()
        inactive_gr = active_info_gr[~active_info_gr.Active]

        if inactive_gr.empty:
            raise SparseInactiveRegionError(f'{self}')

        max_width_index = inactive_gr.lengths().argmax()
        max_width_row = inactive_gr.df.iloc[max_width_index, :]
        return int((max_width_row.Start + max_width_row.End) / 2)

    def get_result_vcfspecs(self, as_components=False, merge=True, subseq_portion_threshold=None, MQ_threshold=40): 
        if subseq_portion_threshold is None:
            subseq_portion_threshold = self.allele_portion_threshold

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


class MultisampleRealignerPileup(RealignerPileupBase):
    def __init__(
        self, 
        pileup_dict, 
        refver=None,
        fasta=None, 
        verbose=False, logger=None,
    ):
        """Args:
            pileup_dict: keys - sampleid; values - Pileup object
        """
        # set params
        self.pileup_dict = pileup_dict
        self.refver, self.fasta = RealignerPileupBase.refver_fasta_arghandler(refver, fasta)
        self.allele_portion_threshold = self.first_pileup.allele_portion_threshold

        # set logger
        self.verbose, self.logger = RealignerPileupBase.logger_arghandler(verbose, logger)

        # setup data structure
        self.set_df()
        self.set_row_spec_groups()
        self.save_subseq_alignments()
        self.allocate_subseq_alignments()
        self.set_row_spec_groups_bysample()

    @property
    def first_pileup(self):
        return next(iter(self.pileup_dict.values()))

    @property
    def chrom(self):
        return self.first_pileup.chrom

    @property
    def start0(self):
        return self.first_pileup.start0

    @property
    def end0(self):
        return self.first_pileup.end0

    @property
    def range0(self):
        return self.first_pileup.range0

    def set_df(self):
        assert len(set(x.range0 for x in self.pileup_dict.values())) == 1, f'Genomic ranges of pileup objects are different.'

        self.df = pd.concat(
            {key: val.df for key, val in self.pileup_dict.items()},
            names=['SampleID', 'ReadUID'],
        )

    def row_spec_getter(self, key):
        return self.pileup_dict[key[0]].row_specs[key[1]]

    def set_row_spec_groups(self):
        # initialize row_spec_groups
        self.row_spec_groups = dict()
        candidates = dict()
        for superseq_sampleid, sub_pileup in self.pileup_dict.items():
            for superseq_readuid, groupinfo in sub_pileup.row_spec_groups.items():
                superseq_key = (superseq_sampleid, superseq_readuid)
                supserseq_row_spec = sub_pileup.row_specs[superseq_readuid]
                candidates.setdefault(supserseq_row_spec['seq'], list())
                candidates[supserseq_row_spec['seq']].append((superseq_key, supserseq_row_spec))
        for key, val in candidates.items():
            selected_superseq_key = max(
                val, 
                key=(lambda x: (x[1]['span_length'], -x[1]['cliplen'])),
            )[0]
            self.row_spec_groups[selected_superseq_key] = {
                'subseq_rowids': list(),
                'subseq_hits': 0,
            }
        # add subseq entries
        for subseq_sampleid, sub_pileup in self.pileup_dict.items():
            for subseq_readuid, subseq_row_spec in sub_pileup.row_specs.items():
                superseq_candidates = list()
                for superseq_key in self.row_spec_groups.keys():
                    superseq_sampleid, superseq_readuid = superseq_key
                    superseq_row_spec = self.pileup_dict[superseq_sampleid].row_specs[superseq_readuid]
                    if self.row_spec_matcher(subseq_row_spec, superseq_row_spec):
                        superseq_candidates.append(superseq_key)
                        
                if len(superseq_candidates) == 0:
                    raise Exception(f'Subseq row_spec does not match with any of superseq row_specs')
                else:
                    superseq_key = self.superseq_candidates_tiebreak(superseq_candidates, self.row_spec_getter, self.row_spec_groups)
                    groupinfo = self.row_spec_groups[superseq_key]
                    groupinfo['subseq_rowids'].append((subseq_sampleid, subseq_readuid))
                    groupinfo['subseq_hits'] += 1

        # discard 0-hit groups
        no_hits = [
            subseq_key for subseq_key, groupinfo in self.row_spec_groups.items()
            if groupinfo['subseq_hits'] == 0
        ]
        for subseq_key in no_hits:
            del self.row_spec_groups[subseq_key]

        self.decorate_groupinfo(self.row_spec_groups, self.read_getter)

    def read_getter(self, key):
        return self.pileup_dict[key[0]].read_store[key[1]]

    def save_subseq_alignments(self):
        self._subseq_alignments = dict()
        for superseq_key, groupinfo in self.row_spec_groups.items():
            sup_to_ref_aln = self.pileup_dict[superseq_key[0]]._superseq_alignments[superseq_key[1]]
            for subseq_key in groupinfo['subseq_rowids']:
                subseq_row_spec = self.pileup_dict[subseq_key[0]].row_specs[subseq_key[1]]
                try:
                    self._subseq_alignments[subseq_key] = self.align_subseq_to_ref(subseq_row_spec, sup_to_ref_aln)
                except:
                    print(superseq_key)
                    print(subseq_key)
                    raise

    def allocate_subseq_alignments(self):
        subseq_alignments_bysample = {sampleid: dict() for sampleid in self.pileup_dict.keys()}
        for subseq_key, aln in self._subseq_alignments.items():
            subseq_alignments_bysample[subseq_key[0]][subseq_key[1]] = aln

        for sampleid, sub_pileup in self.pileup_dict.items():
            sub_pileup._subseq_alignments = subseq_alignments_bysample[sampleid]

    def set_row_spec_groups_bysample(self):
        self.row_spec_groups_bysample = {sampleid: dict() for sampleid in self.pileup_dict.keys()}
        for sampleid, subgroups in self.row_spec_groups_bysample.items():
            for superseq_key, groupinfo in self.row_spec_groups.items():
                subgroups[superseq_key] = {
                    'subseq_rowids': list(),
                    'subseq_hits': 0,
                }
            
        for superseq_key, groupinfo in self.row_spec_groups.items():
            for subseq_key in groupinfo['subseq_rowids']:
                subgroups = self.row_spec_groups_bysample[subseq_key[0]]
                subgroups[superseq_key]['subseq_rowids'].append(subseq_key)
                subgroups[superseq_key]['subseq_hits'] += 1

        for sampleid, subgroups in self.row_spec_groups_bysample.items():
            self.decorate_groupinfo(subgroups, self.read_getter)

    def iter_contig_vcfspecs(self, subseq_portion_threshold, MQ_threshold, threshold_bysample=True, verbose=False):
        if threshold_bysample:
            for superseq_key, groupinfo in sorted(
                self.row_spec_groups.items(),
                key=(lambda x: x[1]['subseq_hits']),
                reverse=True,
            ):
                if any(
                    (subgroups[superseq_key]['subseq_hits_portion'] >= subseq_portion_threshold) and
                    (subgroups[superseq_key]['mean_MQ'] >= MQ_threshold)
                    for subgroups in self.row_spec_groups_bysample.values()
                ):
                    if verbose:
                        print('vcfspec:', self.pileup_dict[superseq_key[0]].contig_vcfspecs[superseq_key[1]])
                        print('vcfspec components:', self.pileup_dict[superseq_key[0]].contig_vcfspecs[superseq_key[1]].components)
                        for sampleid, subgroups in self.row_spec_groups_bysample.items():
                            print('sample ID:', sampleid, 'portion:', subgroups[superseq_key]['subseq_hits_portion'], 'meanMQ:', subgroups[superseq_key]['mean_MQ'])
                        print()

                    yield superseq_key, self.pileup_dict[superseq_key[0]].contig_vcfspecs[superseq_key[1]]
        else:
            for superseq_key, groupinfo in sorted(
                self.row_spec_groups.items(),
                key=(lambda x: x[1]['subseq_hits']),
                reverse=True,
            ):
                if (
                    (groupinfo['subseq_hits_portion'] >= subseq_portion_threshold) and
                    (groupinfo['mean_MQ'] >= MQ_threshold)
                ):
                    if verbose:
                        print('vcfspec:', self.pileup_dict[superseq_key[0]].contig_vcfspecs[superseq_key[1]])
                        print('vcfspec components:', self.pileup_dict[superseq_key[0]].contig_vcfspecs[superseq_key[1]].components)
                        print('portion:', groupinfo['subseq_hits_portion'], 'meanMQ:', groupinfo['mean_MQ'])
                        print()

                    yield superseq_key, self.pileup_dict[superseq_key[0]].contig_vcfspecs[superseq_key[1]]

    # get_result_vcfspec BEGIN #

    def get_result_vcfspec_recalc_hits(self, vaf_cutoff=0.01, verbose=False):
        candidate_monoalts = list()
        for superseq_key in self.row_spec_groups.keys():
            contig_vcfspec = self.pileup_dict[superseq_key[0]].contig_vcfspecs[superseq_key[1]]
            if contig_vcfspec is not None:
                candidate_monoalts.append(contig_vcfspec.normalize())
        if len(candidate_monoalts) == 0:
            return None

        # create a range which covers all candidate mono-ALT vcfspecs
        merged_ref_range0 = range(
            min(x.start0 for x in candidate_monoalts),
            max(x.end0 for x in candidate_monoalts),
        )
        # make rpplist using the range with each sample
        rpplist_dict = dict()
        for sampleid, pileup in self.pileup_dict.items():
            rpplist_dict[sampleid] = readplus.get_rpplist_nonsv(
                bam=pileup.bam, 
                fasta=self.fasta, 
                chromdict=common.DEFAULT_CHROMDICTS[self.refver],
                chrom=self.chrom,
                start0=merged_ref_range0.start,
                end0=merged_ref_range0.stop,
            )
        # calculate vaf for each mono-ALT vcfspec from rpplists
        # choose if vaf is greater than cutoff in any sample
        selected_monoalt_vcfspecs = set()
        for monoalt_vcfspec in candidate_monoalts:
            if verbose:
                print()
                print(monoalt_vcfspec)

            for sampleid, rpplist in rpplist_dict.items():
                rpplist.update_alleleclass(monoalt_vcfspec)
                alleleclass_counter = collections.Counter(
                    rpp.alleleclass[monoalt_vcfspec] for rpp in rpplist
                )
                total_reads = sum(
                    v for (k, v) in alleleclass_counter.items()
                    if k is not None
                )
                var_reads = alleleclass_counter[1]
                if total_reads == 0:
                    vaf = None
                else:
                    vaf = var_reads / total_reads

                if verbose:
                    print(f'sampleid {sampleid}')
                    print(f'vaf {vaf}')
                    
                if (vaf is not None) and (vaf >= vaf_cutoff):
                    selected_monoalt_vcfspecs.add(monoalt_vcfspec)
                    break

        if len(selected_monoalt_vcfspecs) == 0:
            return None

        # sort selected vcfspecs
        selected_monoalt_vcfspecs = sorted(
            selected_monoalt_vcfspecs,
            key=(lambda x: (x.pos, x.alts)),
        )
        # merge and return
        return functools.reduce(libvcfspec.merge, selected_monoalt_vcfspecs)

    def get_result_vcfspec_recalc_hits_dual_filtering(
        self, vaf_cutoff=0.01, subseq_portion_threshold=None, 
        MQ_threshold=40, threshold_bysample=True,
        verbose=False,
    ):
        # set params
        if subseq_portion_threshold is None:
            subseq_portion_threshold = self.allele_portion_threshold

        # 1st filtering
        # make candidate mono-ALT vcfspecs
        candidate_monoalts = list()
        for x in self.iter_contig_vcfspecs(
            subseq_portion_threshold=subseq_portion_threshold, 
            MQ_threshold=MQ_threshold, 
            threshold_bysample=threshold_bysample,
        ):
            if x[1] is not None:
                candidate_monoalts.append(x[1].normalize())

        if len(candidate_monoalts) == 0:
            return None

        # 2nd filtering
        # create rpplists for each component sample
        merged_ref_range0 = range(
            min(x.start0 for x in candidate_monoalts),
            max(x.end0 for x in candidate_monoalts),
        )
        rpplist_dict = dict()
        for sampleid, pileup in self.pileup_dict.items():
            rpplist_dict[sampleid] = readplus.get_rpplist_nonsv(
                bam=pileup.bam, 
                fasta=self.fasta, 
                chromdict=common.DEFAULT_CHROMDICTS[self.refver],
                chrom=self.chrom,
                start0=merged_ref_range0.start,
                end0=merged_ref_range0.stop,
            )
        # calculate vaf for each mono-ALT vcfspec from rpplists
        # choose if vaf is greater than cutoff in any sample
        selected_monoalt_vcfspecs = set()
        for monoalt_vcfspec in candidate_monoalts:
            if verbose:
                print()
                print(monoalt_vcfspec)

            for sampleid, rpplist in rpplist_dict.items():
                rpplist.update_alleleclass(monoalt_vcfspec)
                alleleclass_counter = collections.Counter(
                    rpp.alleleclass[monoalt_vcfspec] for rpp in rpplist
                )
                total_reads = sum(
                    v for (k, v) in alleleclass_counter.items()
                    if k is not None
                )
                var_reads = alleleclass_counter[1]
                if total_reads == 0:
                    vaf = None
                else:
                    vaf = var_reads / total_reads

                if verbose:
                    print(f'sampleid {sampleid}')
                    print(f'vaf {vaf}')
                    
                if (vaf is not None) and (vaf >= vaf_cutoff):
                    selected_monoalt_vcfspecs.add(monoalt_vcfspec)
                    break

        if len(selected_monoalt_vcfspecs) == 0:
            return None

        # sort selected vcfspecs
        selected_monoalt_vcfspecs = sorted(
            selected_monoalt_vcfspecs,
            key=(lambda x: (x.pos, x.alts)),
        )
        # merge and return
        return functools.reduce(libvcfspec.merge, selected_monoalt_vcfspecs)

    def get_result_vcfspec_use_rowspecgroup_hits(
        self, 
        #as_components=False, 
        #merge=True, 
        subseq_portion_threshold=None, 
        MQ_threshold=30, 
        threshold_bysample=True,
        verbose=False,
    ):
        if subseq_portion_threshold is None:
            subseq_portion_threshold = self.allele_portion_threshold

#        if as_components:
#            result = list()
#            for row_id, contig_vcfspec in self.iter_contig_vcfspecs(subseq_portion_threshold=subseq_portion_threshold, MQ_threshold=MQ_threshold, threshold_bysample=threshold_bysample):
#                if contig_vcfspec is None:
#                    continue
#
#                if len(contig_vcfspec.components[0]) == 0:
#                    result.append(contig_vcfspec)
#                else:
#                    result.extend(contig_vcfspec.components[0])
        result = list(
            x[1] for x in self.iter_contig_vcfspecs(
                subseq_portion_threshold=subseq_portion_threshold, 
                MQ_threshold=MQ_threshold, 
                threshold_bysample=threshold_bysample,
                verbose=verbose,
            )
            if x[1] is not None
        )

        if len(result) == 0:
            return None
        else:
            return functools.reduce(libvcfspec.merge, result)

    # SELECT AN APPROPRIATE FUNCTION
    get_result_vcfspec = get_result_vcfspec_use_rowspecgroup_hits

    # get_result_vcfspec END #

    def show_row_spec_alignments(self, skip_zero_hits=True, show_subseqs=False):
        self._show_row_spec_alignments_helper(
            self.row_spec_groups, 
            lambda x: self.pileup_dict[x[0]].row_specs.__getitem__(x[1]),
            lambda x: self.pileup_dict[x[0]]._superseq_alignments.__getitem__(x[1]),
            skip_zero_hits, 
            show_subseqs,
        )


class RealignerPileupSeries:
    def __init__(
        self, 
        bam, 
        chrom, start0, end0, 
        refver=None, 
        fasta=None,
        max_series_width=DEFAULT_RPILEUPSERIES_MAX_WIDTH,
        verbose=False,
        logger=None,
        **kwargs,
    ):
        """Args:
            **kwargs: keyword arguments of RealignerPileup.__init__
                except start0, end0, init_df
        """
        # set params
        self.chrom = chrom
        self.bam = bam
        self.refver, self.fasta = RealignerPileupBase.refver_fasta_arghandler(refver, fasta)
        self.max_series_width = max_series_width

        # set logger
        self.verbose, self.logger = RealignerPileupBase.logger_arghandler(verbose, logger)

        # set RealignerPileup kwargs
        self.rpileup_init_kwargs = parse_rpileup_kwargs(**kwargs)
        for k, v in self.rpileup_init_kwargs.items():
            setattr(self, k, v)

        # setup pileup_list
        self._init_pileup_list(start0, end0)
        # merge individual read_store's
        self.read_store = dict()
        for pileup in self.pileup_list:
            self.read_store.update(pileup.read_store)

    def _init_pileup_list(self, start0, end0):
        # initiate pileup_list
        self.logger.debug(f'Initial pileup generation')

        initial_pileup, result_inactive, result_vcfspec = RealignerPileup.init_and_augment(
            bam=self.bam, 
            chrom=self.chrom, start0=start0, end0=end0, 
            refver=self.refver, fasta=self.fasta,
            verbose=self.verbose, 
            logger=self.logger, 
            **self.rpileup_init_kwargs,
        )
        self.pileup_list = [initial_pileup]

        self.logger.debug(f'Initial pileup: {initial_pileup}')
        self.logger.debug(f'result_inactive: {result_inactive}')
        self.logger.debug(f'result_vcfspec: {result_vcfspec}')
        self.logger.debug(f'Initial pileup generation Finished\n')

        if (result_inactive is None) and (result_vcfspec is None):
            # When there is no active region in the initial pileup
            return

        left_insufficient = (not result_inactive.left_okay) or (not result_vcfspec.left_okay)
        right_insufficient = (not result_inactive.right_okay) or (not result_vcfspec.right_okay)

        if left_insufficient:
            self.secure_left()
        if right_insufficient:
            self.secure_right()

        self.logger.debug(f'Finished initial pileup generation\n')

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

    def get_gr(self):
        start0s = list()
        end0s = list()
        for pileup in self.pileup_list:
            start0s.append(pileup.start0)
            end0s.append(pileup.end0)
        chroms = [self.chrom] * len(self.pileup_list)
        self_indexes = list(range(len(self.pileup_list)))
        return pr.from_dict(
            {'Chromosome': chroms, 'Start': start0s, 'End': end0s, 'Self_index': self_indexes}
        )

    def extend_left(self, width):
        """Ignores MQ_limit, depth_limit, start0_limit, end0_limit,
        max_series_width attribute"""
        original_max_series_width = self.max_series_width
        self.max_series_width = np.inf
        original_start0_limit = self.start0_limit
        self.start0_limit = -np.inf

        target_start0 = self.start0 - width

        rpileup_init_kwargs = self.rpileup_init_kwargs.copy()
        rpileup_init_kwargs['start0_limit'] = target_start0
        rpileup_init_kwargs['MQ_limit'] = 0
        rpileup_init_kwargs['depth_limit'] = 0

        while True:
            new_pileup = RealignerPileup(
                bam=self.bam, 
                chrom=self.chrom, start0=(self.start0 - 1), end0=self.start0, 
                refver=self.refver, fasta=self.fasta, 
                verbose=self.verbose, logger=self.logger,
                **rpileup_init_kwargs,
            )
            self.pileup_list.insert(0, new_pileup)
            result_inactive, result_vcfspec, touched_width_limit = self.secure_left()
            if self.start0 <= target_start0:
                break

        self.max_series_width = original_max_series_width
        self.start0_limit = original_start0_limit

    def extend_right(self, width):
        """Ignores MQ_limit, depth_limit, start0_limit, end0_limit,
        max_series_width attribute"""
        original_max_series_width = self.max_series_width
        self.max_series_width = np.inf
        original_end0_limit = self.end0_limit
        self.end0_limit = np.inf

        target_end0 = self.end0 + width

        rpileup_init_kwargs = self.rpileup_init_kwargs.copy()
        rpileup_init_kwargs['end0_limit'] = target_end0
        rpileup_init_kwargs['MQ_limit'] = 0
        rpileup_init_kwargs['depth_limit'] = 0

        while True:
            new_pileup = RealignerPileup(
                bam=self.bam, 
                chrom=self.chrom, start0=self.end0, end0=(self.end0 + 1), 
                refver=self.refver, fasta=self.fasta, 
                verbose=self.verbose, logger=self.logger,
                **rpileup_init_kwargs,
            )
            self.pileup_list.append(new_pileup)
            result_inactive, result_vcfspec, touched_width_limit = self.secure_right()
            if self.end0 >= target_end0:
                break

        self.max_series_width = original_max_series_width
        self.end0_limit = original_end0_limit

    # init helpers #
    def secure_left(self):
        self.logger.debug(f'Beginning RealignerPileupSeries secure_left')

        left_pileup = self.pileup_list[0]
        left_pileup.end0_limit = left_pileup.end0
        while True:
            result_inactive, result_vcfspec = left_pileup.secure_margins()

            self.logger.debug(self.pileup_list)
            self.logger.debug(f'self.start0: {self.start0}')
            self.logger.debug(f'result_inactive: {result_inactive}')
            self.logger.debug(f'result_vcfspec: {result_vcfspec}')

            if (
                (result_inactive.left_okay and result_vcfspec.left_okay) or
                (
                    result_vcfspec.left_low_depth or 
                    result_vcfspec.left_low_MQ or 
                    result_vcfspec.touched_left_limit
                ) or
                (self.width > self.max_series_width)
            ):
                break
            else:
                assert result_vcfspec.touched_width_limit
                split_points = left_pileup.get_split_points()
                if len(split_points) == 0:
                    rpileup_init_kwargs = self.rpileup_init_kwargs.copy()
                    new_pileup = RealignerPileup(
                        bam=self.bam, 
                        chrom=self.chrom, start0=(left_pileup.start0 - 1), end0=left_pileup.start0, 
                        refver=self.refver, fasta=self.fasta, 
                        verbose=self.verbose, logger=self.logger,
                        **rpileup_init_kwargs,
                    )
                    self.pileup_list.insert(0, new_pileup)
                else:
                    split_pileups = left_pileup.split(split_points[:1])
                    del self.pileup_list[0]
                    self.pileup_list = split_pileups + self.pileup_list

                left_pileup = self.pileup_list[0]
                left_pileup.end0_limit = left_pileup.end0

        touched_width_limit = self.width >= self.max_series_width

        self.logger.debug(f'Finished RealignerPileupSeries secure_left')

        return result_inactive, result_vcfspec, touched_width_limit

    def secure_right(self):
        self.logger.debug(f'Beginning RealignerPileupSeries secure_right')

        right_pileup = self.pileup_list[-1]
        right_pileup.start0_limit = right_pileup.start0
        while True:
            result_inactive, result_vcfspec = right_pileup.secure_margins()

            self.logger.debug(self.pileup_list)
            self.logger.debug(f'self.end0: {self.end0}')
            self.logger.debug(f'result_inactive: {result_inactive}')
            self.logger.debug(f'result_vcfspec: {result_vcfspec}')

            if (
                (result_inactive.right_okay and result_vcfspec.right_okay) or
                (
                    result_vcfspec.right_low_depth or 
                    result_vcfspec.right_low_MQ or
                    result_vcfspec.touched_right_limit
                ) or
                (self.width > self.max_series_width)
            ):
                break
            else:
                assert result_vcfspec.touched_width_limit
                split_points = right_pileup.get_split_points()
                if len(split_points) == 0:
                    rpileup_init_kwargs = self.rpileup_init_kwargs.copy()
                    new_pileup = RealignerPileup(
                        bam=self.bam, 
                        chrom=self.chrom, start0=right_pileup.end0, end0=(right_pileup.end0 + 1), 
                        refver=self.refver, fasta=self.fasta, 
                        verbose=self.verbose, logger=self.logger,
                        **rpileup_init_kwargs,
                    )
                    self.pileup_list.append(new_pileup)
                else:
                    split_pileups = right_pileup.split(split_points[-1:])
                    del self.pileup_list[-1]
                    self.pileup_list.extend(split_pileups)

                right_pileup = self.pileup_list[-1]
                right_pileup.start0_limit = right_pileup.start0

        touched_width_limit = self.width >= self.max_series_width

        self.logger.debug(f'Finished RealignerPileupSeries secure_right')

        return result_inactive, result_vcfspec, touched_width_limit

    def get_splittable_region_best(self):
        # without splitting contig vcfspec
        return pr.concat(
            [
                rpileup.get_vcfspec_margins_gr(
                    subseq_portion_threshold=(rpileup.allele_portion_threshold * 2),
                    inverse=True,
                    split_contig_vcfspec=False,
                )
                for rpileup in self.pileup_list
            ]
        ).merge()

    def get_splittable_region_split_contig(self):
        return pr.concat(
            [
                rpileup.get_vcfspec_margins_gr(
                    subseq_portion_threshold=(rpileup.allele_portion_threshold * 2),
                    inverse=True,
                    split_contig_vcfspec=True,
                )
                for rpileup in self.pileup_list
            ]
        ).merge()

    def get_splittable_region_inactive_runs(self):
        inactive_runs_gr_list = list()
        for rpileup in self.pileup_list:
            active_info_gr = rpileup.get_active_info_gr()
            inactive_gr = active_info_gr[~active_info_gr.Active]
            inactive_runs_gr_list.append(inactive_gr)
        result = pr.concat(inactive_runs_gr_list).merge()
        if result.empty:
            raise Exception(f'No inactive runs in this RealignerPileupSeries object.')
        return result

    def rearrange(self, start0_list, prepare_vcfspecs=True):
        """Mutates in-place"""
        def merger(left_pileup, right_pileup):
            left_pileup.merge(right_pileup, other_on_left=False)
            return left_pileup
        
        # sanity check
        if not (start0_list[0] >= self.start0 and start0_list[-1] <= self.end0):
            raise Exception(f'Input splitting range is out of PileupSeries object range.')
        
        # make new pileup_list
        split_ranges_gr = pr.from_dict(
            {
                'Chromosome': [self.chrom] * (len(start0_list) - 1),
                'Start': start0_list[:-1],
                'End': start0_list[1:],
                'Ranges_index': list(range(len(start0_list) - 1)),
            }
        )
        joined_gr = split_ranges_gr.join(self.get_gr())

        new_pileup_list = list()
        for ranges_index, subiter in itertools.groupby(
            (x[1] for x in joined_gr.df.iterrows()),
            key=(lambda x: x.Ranges_index),
        ):
            partial_pileups = list()
            for row in subiter:
                new_start0 = max(row.Start, row.Start_b)
                new_end0 = min(row.End, row.End_b)
                partial_pileups.append(
                    self.pileup_list[row.Self_index].subset(new_start0, new_end0)
                )
            new_pileup_list.append(functools.reduce(merger, partial_pileups))

        # prepare vcfspecs
        if prepare_vcfspecs:
            for pileup in new_pileup_list:
                pileup.prepare_vcfspecs()

        # result
        self.pileup_list = new_pileup_list

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
        for idx, (left_item, right_item) in enumerate(common.pairwise(cigartuples_list_decorated)):
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
        alignment = RealignerPileup.main_aligner(
            target=ref_seq,
            query=read_seq_right[:cigartuples_right[0][1]],
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
        alignment = RealignerPileup.main_aligner(
            target=ref_seq,
            query=read_seq_left[-cigartuples_left[-1][1]:],
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
                if read_uid in pileup.df.index:
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
            raise RealignmentCigarError(
                f'Read query sequence length and cigar length differs.\n'
                f'Original read:\n{original_read.to_string()}\n'
                f'Realigned read:\n{realigned_read.to_string()}\n'
            )


class MultisampleRealignerPileupSeries:
    def __repr__(self):
        return str(self.pileupseries_dict)

    def iter_series(self):
        return iter(self.pileupseries_dict.values())

    def iter_mspileups(self):
        return iter(self.mspileup_list)

    @functools.cached_property
    def df(self):
        data = list()
        index = list()
        for idx, (sampleid, rpileup_series) in enumerate(self.pileupseries_dict.items()):
            index.append(sampleid)
            data.append(list(rpileup_series.pileup_list))
            if idx == 0:
                breaks = list()
                breaks.append(rpileup_series.pileup_list[0].start0)
                for rpileup in rpileup_series.pileup_list:
                    breaks.append(rpileup.end0)
                columns = pd.IntervalIndex.from_breaks(breaks, closed='left')
        return pd.DataFrame.from_records(data, index=index, columns=columns)

#    @property
#    def rows(self):
#        return pd.DataFrame.from_dict(self.pileupseries_dict, orient='index')
#
#    @property
#    def cols(self):
#        return pd.DataFrame.from_records(
#            [
#                {
#                    (x.start0, x.end0): x 
#                    for x in self.mspileup_list
#                }
#            ]
#        )

    #@property
    #def series_length(self):
        #return len(self.rows.iloc[0, :])
    #    return len(next(iter(self.pileupseries_dict.values())))

    @property
    def ranges(self):
        #return [range(*x) for x in self.cols.columns]
        return [x.range0 for x in next(iter(self.pileupseries_dict.values())).pileup_list]

    @property
    def start0(self):
        #return self.rows.iloc[0, :].start0
        return next(self.pileupseries_dict.values()).start0

    @property
    def end0(self):
        #return self.rows.iloc[0, :].end0
        return next(self.pileupseries_dict.values()).end0

    def write_realigned_reads(self, bam_path_dict, padding=0):
        for sampleid, rpileup_ser in self.pileupseries_dict.items():
            rpileup_ser.write_realigned_reads(bam_path_dict[sampleid], padding=padding)

    def iter_contig_vcfspecs(self, subseq_portion_threshold=0, MQ_threshold=40):
        return itertools.chain.from_iterable(
            msrpileup.iter_contig_vcfspecs(subseq_portion_threshold=subseq_portion_threshold, MQ_threshold=MQ_threshold)
            for msrpileup in self.mspileup_list
        )

    #def get_result_vcfspecs(self, vaf_cutoff=0.01, verbose=False):  # for MultisampleRealignerPileup.get_result_vcfspec_recalc_hits
    def get_result_vcfspecs(self, portion_cutoff=0.1, MQ_cutoff=30, verbose=False):
        return list(
            msrpileup.get_result_vcfspec(
                subseq_portion_threshold=portion_cutoff, 
                MQ_threshold=MQ_cutoff, 
                threshold_bysample=True,
                verbose=verbose,
            )
            for msrpileup in self.mspileup_list
        )

    def get_result_vcfspecs_old(
        self, vaf_cutoff=0.01, subseq_portion_threshold=None, 
        MQ_threshold=40, threshold_bysample=True,
        verbose=False,
    ):
        return list(
            msrpileup.get_result_vcfspec_old(
                vaf_cutoff=vaf_cutoff, 
                subseq_portion_threshold=subseq_portion_threshold, 
                MQ_threshold=MQ_threshold, 
                threshold_bysample=threshold_bysample,
                verbose=verbose,
            )
            for msrpileup in self.mspileup_list
        )

    def get_result_vcfspecs_old(
        self, as_components=False, merge=True, subseq_portion_threshold=None, MQ_threshold=40, threshold_bysample=True,
    ): 
        return list(
            msrpileup.get_result_vcfspecs(
                as_components=as_components, 
                merge=merge, 
                subseq_portion_threshold=subseq_portion_threshold, 
                MQ_threshold=MQ_threshold, 
                threshold_bysample=threshold_bysample,
            )
            for msrpileup in self.mspileup_list
        )

    ###########################
    # initializer and helpers #
    ###########################
    def __init__(
        self, 
        bam_dict, 
        chrom, start0, end0, 
        refver=None,
        fasta=None,
        max_series_width=DEFAULT_RPILEUPSERIES_MAX_WIDTH, 
        verbose=False, logger=None, 
        **kwargs,
    ):
        # set params
        self.chrom = chrom
        self.bam_dict = bam_dict
        self.refver, self.fasta = RealignerPileupBase.refver_fasta_arghandler(refver, fasta)
        self.max_series_width = max_series_width

        # set RealignerPileup kwargs
        self.rpileup_init_kwargs = parse_rpileup_kwargs(**kwargs)
        for k, v in self.rpileup_init_kwargs.items():
            setattr(self, k, v)

        # set logger
        self.verbose, self.logger = RealignerPileupBase.logger_arghandler(verbose, logger)

        # setup data
        self.set_series_dict(start0, end0)  # self.pileupseries_dict
        self.set_multisample_pileups()  # self.mspileup_list

        #for sampleid, series in self.pileupseries_dict.items():
        #    print(sampleid)
        #    for rpup in series.pileup_list:
        #        print(rpup)
        #        print(rpup._subseq_alignments)

        self.set_realigned_reads()

    def set_series_dict(self, seed_start0, seed_end0):
        self.pileupseries_dict = dict()
        self.no_variant = False

        # initialize
        for sampleid, bam in self.bam_dict.items():
            self.logger.debug(f'@@@ Initializing RealignerPileupSeries of sample {sampleid} @@@\n')
            self.pileupseries_dict[sampleid] = RealignerPileupSeries(
                bam=bam, 
                chrom=self.chrom, start0=seed_start0, end0=seed_end0, 
                refver=self.refver,
                fasta=self.fasta, 
                max_series_width=self.max_series_width,
                verbose=self.verbose,
                logger=self.logger,
                **self.rpileup_init_kwargs,
            )
            self.logger.debug(f'@@@ Finished initialization of RealignerPileupSeries of sample {sampleid} @@@\n\n')

        # equalize whole margins
        self.equalize_left()
        self.equalize_right()

        # equalize sub-pileup margins
        self.logger.debug(f'@@@ Beginning inner margin equalization @@@\n')
        self.equalize_inner_margins()
        self.logger.debug(f'@@@ Finished inner margin equalization @@@\n')

    def equalize_left(self):
        self.logger.debug(f'@@@ Beginning equalize_left @@@\n')

        # equalize left
        while True:
            if len(set(x.start0 for x in self.pileupseries_dict.values())) == 1:
                break
            # extend
            target_start0 = min(x.start0 for x in self.pileupseries_dict.values())
            self.logger.debug(f'target_start0: {target_start0}')
            for sampleid, pileup_ser in self.pileupseries_dict.items():
                width = pileup_ser.start0 - target_start0
                if width > 0:
                    self.logger.debug(f'Beginning EXTEND of {sampleid}')
                    self.logger.debug(f'current start0 of {sampleid}: {pileup_ser.start0}')
                    pileup_ser.extend_left(width)
                    self.logger.debug(f'Finished EXTEND of {sampleid}\n')

            # check if hit width limit
            if any(pileup_ser.width >= self.max_series_width for pileup_ser in self.pileupseries_dict.values()):
                break
            # secure
            for sampleid, pileup_ser in self.pileupseries_dict.items():
                self.logger.debug(f'Beginning SECURE of {sampleid}')
                pileup_ser.secure_left()
                self.logger.debug(f'Finished SECURE of {sampleid}\n')

        self.logger.debug(f'@@@ Finished equalize_left @@@\n')

    def equalize_right(self):
        self.logger.debug(f'@@@ Beginning equalize_right @@@\n')

        while True:
            if len(set(x.end0 for x in self.pileupseries_dict.values())) == 1:
                break
            # extend
            target_end0 = max(x.end0 for x in self.pileupseries_dict.values())
            self.logger.debug(f'target_end0: {target_end0}')
            for sampleid, pileup_ser in self.pileupseries_dict.items():
                width = target_end0 - pileup_ser.end0
                if width > 0:
                    self.logger.debug(f'Beginning EXTEND of {sampleid}')
                    self.logger.debug(f'current end0 of {sampleid}: {pileup_ser.end0}')
                    pileup_ser.extend_right(width)
                    self.logger.debug(f'Finished EXTEND of {sampleid}\n')

            # check if hit width limit
            if any(pileup_ser.width >= self.max_series_width for pileup_ser in self.pileupseries_dict.values()):
                break
            # secure
            for pileup_ser in self.pileupseries_dict.values():
                self.logger.debug(f'Beginning SECURE of {sampleid}')
                pileup_ser.secure_right()
                self.logger.debug(f'Finished SECURE of {sampleid}\n')

        self.logger.debug(f'@@@ Finished equalize_right @@@\n')

    def equalize_inner_margins(self):
        # set interim parameters
        first_pileupseries = next(self.iter_series())
        series_start0 = first_pileupseries.start0
        series_end0 = first_pileupseries.end0
        series_width = series_end0 - series_start0
        max_pileup_width = first_pileupseries.pileup_list[0].max_pileup_width

        # When there is no need to further split the series range
        if series_end0 - series_start0 <= max_pileup_width:
            start0_list = [series_start0, series_end0]
            for pileup_ser in self.pileupseries_dict.values():
                if len(pileup_ser.pileup_list) > 1:
                    pileup_ser.rearrange(start0_list, prepare_vcfspecs=True)
            self.rearrange_mode = 'not_done'
            return

        # best case - without splitting contig vcfspec
        start0_list, all_splittable = self.get_rearrangement_points(
            split_region_gr=functools.reduce(
                lambda x, y: x.intersect(y), 
                (pileup_ser.get_splittable_region_best() for pileup_ser in self.pileupseries_dict.values())
            ), 
            trim_margins=True, 
            max_pileup_width=max_pileup_width, 
            series_start0=series_start0, 
            series_end0=series_end0,
        )

        if start0_list is not None:
            for pileup_ser in self.pileupseries_dict.values():
                pileup_ser.rearrange(start0_list, prepare_vcfspecs=True)
            self.rearrange_mode = 'preserve_contig'
            return
        else:
            if all_splittable:  # no variant case
                # split series range into evenly-spaced subsets
                start0_list = np.linspace(
                    series_start0,
                    series_end0,
                    np.ceil(series_width / max_pileup_width).astype('int') + 1,
                    endpoint=True,
                ).astype('int')
                for pileup_ser in self.pileupseries_dict.values():
                    pileup_ser.rearrange(start0_list, prepare_vcfspecs=True)

                self.no_variant = True
                self.rearrange_mode = None
                return

        # splitting contig vcfspec
        start0_list, all_splittable = self.get_rearrangement_points(
            split_region_gr=functools.reduce(
                lambda x, y: x.intersect(y), 
                (pileup_ser.get_splittable_region_split_contig() for pileup_ser in self.pileupseries_dict.values())
            ), 
            trim_margins=True, 
            max_pileup_width=max_pileup_width, 
            series_start0=series_start0, 
            series_end0=series_end0,
        )

        if start0_list is not None:
            for pileup_ser in self.pileupseries_dict.values():
                pileup_ser.rearrange(start0_list, prepare_vcfspecs=True)
            self.rearrange_mode = 'split_contig'
            return
        else:
            if all_splittable:
                raise Exception(f'all-splittable is True in split-contig after False in non-split-contig')

        # using inactive runs
        start0_list, all_splittable = self.get_rearrangement_points(
            split_region_gr=functools.reduce(
                lambda x, y: x.intersect(y), 
                (pileup_ser.get_splittable_region_inactive_runs() for pileup_ser in self.pileupseries_dict.values())
            ), 
            trim_margins=False, 
            max_pileup_width=max_pileup_width, 
            series_start0=series_start0, 
            series_end0=series_end0,
        )

        if start0_list is not None:
            for pileup_ser in self.pileupseries_dict.values():
                pileup_ser.rearrange(start0_list, prepare_vcfspecs=True)
            self.rearrange_mode = 'inactive_runs'
            return
        else:
            if all_splittable:
                raise Exception(f'All PileupSeries region is inactive.')
            else:
                self.rearrange_mode = 'failed'
                return

    @staticmethod
    def get_rearrangement_points(split_region_gr, trim_margins, max_pileup_width, series_start0, series_end0):
        """Returns:
            start0_list: list of 0-based positions. 
                With pairwise iteration, for each iteratation result (x, y),
                x is used as start0 and y as end0 for a RealignerPileup object.
        """
        if split_region_gr.empty:
            # abort this split strategy
            start0_list = None
            all_splittable = False
            return start0_list, all_splittable

        # when entire range is splittable
        if split_region_gr.df.shape[0] == 1:
            row = split_region_gr.df.iloc[0, :]
            if (row.Start == series_start0) and (row.End == series_end0):
                start0_list = None
                all_splittable = True
                return start0_list, all_splittable

        # make into ranges
        start_candidate_ranges = [range(row.Start, row.End + 1) for (idx, row) in split_region_gr.df.iterrows()]

        # trim or add margins
        if start_candidate_ranges[0].start == series_start0:
            if trim_margins:
                start_candidate_ranges[0] = range(start_candidate_ranges[0].stop - 1, start_candidate_ranges[0].stop)
        else:
            start_candidate_ranges.insert(0, range(series_start0, series_start0 + 1))

        if start_candidate_ranges[-1].stop - 1 == series_end0:
            if trim_margins:
                start_candidate_ranges[-1] = range(start_candidate_ranges[-1].start, start_candidate_ranges[-1].start + 1)
        else:
            start_candidate_ranges.append(range(series_end0, series_end0 + 1))

        # now len(start_candidate_ranges) is at least 2 

        # If candidate range is smaller than max_pileup_width after trimming
        if trim_margins:
            new_start0 = start_candidate_ranges[0].start
            new_end0 = start_candidate_ranges[-1].stop - 1
            if new_end0 - new_start0 <= max_pileup_width:
                start0_list = [new_start0, new_end0]
                all_splittable = False
                return start0_list, all_splittable

        # validity check
        if any(
            rng2.start - (rng1.stop - 1) > max_pileup_width
            for rng1, rng2 in common.pairwise(start_candidate_ranges)
        ):
            # abort this split strategy
            start0_list = None
            all_splittable = False
            return start0_list, all_splittable

        # search for actual split points
        iterator = iter(start_candidate_ranges)
        current_rng = next(iterator)
        next_rng = next(iterator)
        start0_list = [current_rng.start]
        stopiter = False

        while True:
            if stopiter:
                break

            last_start0 = start0_list[-1]
            while True:
                if next_rng.start <= last_start0 + max_pileup_width:
                    current_rng = next_rng

                    try:
                        next_rng = next(iterator)
                    except StopIteration:
                        stopiter = True
                        break
                else:
                    if current_rng.stop - 1 <= last_start0:
                        raise Exception(f'Cannot make pileup margins within max_pileup_width')
                    start0_list.append(
                        min(current_rng.stop - 1, last_start0 + max_pileup_width)
                    )
                    break

        start0_list.append(next_rng.start)

        all_splittable = False
        return start0_list, all_splittable

    def set_multisample_pileups(self):
        self.mspileup_list = list()
        num_col = len(next(self.iter_series()).pileup_list)
        for idx in range(num_col):
            pileup_dict = {
                sampleid: self.pileupseries_dict[sampleid].pileup_list[idx]
                for sampleid in self.pileupseries_dict.keys()
            }
            self.mspileup_list.append(
                MultisampleRealignerPileup(
                    pileup_dict=pileup_dict, 
                    refver=self.refver,
                    fasta=self.fasta, 
                    verbose=self.verbose, logger=self.logger,
                )
            )

    def set_realigned_reads(self):
        """Assumes all mspileups have executed 'save_subseq_alignments' and 'allocate_subseq_alignments'"""
        #for msrpileup in self.mspileup_list:
        #    msrpileup.save_subseq_alignments()
        #    msrpileup.allocate_subseq_alignments()
        self.realigned_reads = {sampleid: dict() for sampleid in self.pileupseries_dict.keys()}
        for sampleid, rpileup_ser in self.pileupseries_dict.items():
            rpileup_ser.set_realigned_reads()
            self.realigned_reads[sampleid] = rpileup_ser.realigned_reads

    ###########################
    ###########################


class SparseInactiveRegionError(Exception):
    pass


class RealignmentCigarError(Exception):
    pass


def parse_rpileup_kwargs(**kwargs):
    if any(x not in DEFAULT_RPILEUP_PARAMS.keys() for x in kwargs.keys()):
        raise Exception(f'Allowed kwargs keys are: {list(DEFAULT_RPILEUP_PARAMS.keys())}')

    rpileup_params = kwargs
    for k, v in DEFAULT_RPILEUP_PARAMS.items():
        if k not in rpileup_params:
            rpileup_params[k] = v

    return rpileup_params


#def parse_rpileup_kwargs_old(**kwargs):
#    assert not any(x in kwargs for x in ('init_df',))
#
#    ba = inspect.signature(RealignerPileup.__init__).bind_partial(**kwargs)
#    ba.apply_defaults()
#
#    rpileup_init_kwargs = dict()
#    for k, v in ba.arguments.items():  # these include default kwargs of RealignerPileup.__init__
#        if k in ('start0', 'end0', 'init_df', 'verbose', 'logger'):
#            continue
#        rpileup_init_kwargs[k] = v
#
#    return rpileup_init_kwargs



#def get_realigner_pileup(fasta, bam, chrom, start0, end0, verbose=False, logger=None, **kwargs):
#    rpileup = RealignerPileup(
#        fasta, bam, chrom, start0, end0, 
#        verbose=verbose, 
#        logger=logger, 
#        **kwargs,
#    )
#    result_inactive, result_vcfspec = rpileup.augment_margins()
#    #if (result_inactive is None) and (result_vcfspec is None):
#        # No active region found
#    #    rpileup = None
#
#    return rpileup, result_inactive, result_vcfspec



###################################
# helpers for get_realigned_reads #
###################################



#def realign_subseq(row_id, active_region_pileup, active_range, superseq_alignment):
#    row_spec = active_region_pileup.row_specs[row_id]
#    active_region_alignment = align_subseq_to_ref(row_spec, superseq_alignment)
#    try:
#        realigned_read = align_row_helper(row_id, active_region_pileup, active_range, active_region_alignment)
#    except RealignmentCigarError as exc:
#        raise Exception(f'Subseq realignment failed.\nrow_spec: {row_spec}\n') from exc
#
#    return realigned_read
#
#
#def align_row_helper(row_id, active_region_pileup, active_range, active_region_alignment):
#    # set params
#    original_read = active_region_pileup.read_store[row_id]
#    original_cigar_split = alignhandler.split_cigar(original_read.cigartuples, original_read.reference_start, active_range)
#    empty_before_active_region = (len(original_cigar_split[0]) == 0)
#    empty_after_active_region = (len(original_cigar_split[2]) == 0)
#    # get active region cigartuples
#    active_region_cigartuples, active_region_offset = alignhandler.alignment_to_cigartuples(
#        active_region_alignment,
#        match_as_78=False, del_as_skip=False, 
#        left_ins_as_clip=empty_before_active_region, right_ins_as_clip=empty_after_active_region,
#        remove_left_del=empty_before_active_region, remove_right_del=empty_after_active_region,
#    )
#    # merge split cigartuples
#    realigned_cigartuples = functools.reduce(
#        merge_two_cigartuples, 
#        [original_cigar_split[0], active_region_cigartuples, original_cigar_split[2]],
#    )
#    
#    # get new read reference_start
#    if empty_before_active_region:
#        new_reference_start = active_range.start + active_region_offset
#    else:
#        new_reference_start = original_read.reference_start
#    # make new read 
#    realigned_read = original_read.__copy__()
#    realigned_read.reference_start = new_reference_start
#    realigned_read.cigartuples = realigned_cigartuples
#    readhandler.set_NMMD(realigned_read, active_region_pileup.fasta)
#
#    realigned_cigar_sanity_check(original_read, realigned_read, realigned_cigartuples)
#
#    return realigned_read
#
#
#def align_subseq_to_ref(subseq_row_spec, sup_to_ref_aln):
#    # get subseq index
#    reverse_align = (not subseq_row_spec['left_filled']) and subseq_row_spec['right_filled']
#    if reverse_align:
#        sub_to_sup_index = sup_to_ref_aln.query.rfind(subseq_row_spec['seq'])
#    else:
#        sub_to_sup_index = sup_to_ref_aln.query.find(subseq_row_spec['seq'])
#    if sub_to_sup_index == -1:
#        raise Exception(
#            f'"subseq" is not a subseq of "superseq":\n'
#            f'subseq_row_spec: {subseq_row_spec}\n'
#            f'superseq: {sup_to_ref_aln.query}'
#        )
#    # convert into sub-to-ref alignment
#    # supaln: sup_to_ref_aln
#    # subaln: sub_to_sup_aln
#    active_region_walks = list()
#    target_walk_subaln = range(sub_to_sup_index, sub_to_sup_index + len(subseq_row_spec['seq']))
#    query_walk_subaln = range(0, len(subseq_row_spec['seq']))
#    found_initial_hit = False
#    for target_walk_supaln, query_walk_supaln in alignhandler.get_walks(sup_to_ref_aln, copy=False):
#        if not found_initial_hit:
#            if target_walk_subaln.start in query_walk_supaln:
#                found_initial_hit = True
#                # superseq walk
#                superseq_walk_start = target_walk_subaln.start
#                superseq_walk_stop = min(query_walk_supaln.stop, target_walk_subaln.stop)
#                len_current_superseq_walk = superseq_walk_stop - superseq_walk_start
#                # subseq walk
#                subseq_walk = range(0, len_current_superseq_walk)
#                # ref walk
#                if len(target_walk_supaln) == 0:
#                    ref_walk = target_walk_supaln
#                else:
#                    # During initial hit query_walk_supaln cannot be zero-length
#                    ref_walk_start = target_walk_supaln.start + (superseq_walk_start - query_walk_supaln.start)
#                    ref_walk = range(ref_walk_start, ref_walk_start + len_current_superseq_walk)
#
#                active_region_walks.append((ref_walk, subseq_walk))
#
#                if query_walk_supaln.stop >= target_walk_subaln.stop:
#                    break
#        else:
#            # superseq walk
#            superseq_walk_start = query_walk_supaln.start
#            superseq_walk_stop = min(query_walk_supaln.stop, target_walk_subaln.stop)
#            len_current_superseq_walk = superseq_walk_stop - superseq_walk_start
#            # subseq walk
#            subseq_walk_start = active_region_walks[-1][1].stop
#            subseq_walk = range(subseq_walk_start, subseq_walk_start + len_current_superseq_walk)
#            # ref walk
#            if len(target_walk_supaln) == 0 or len(query_walk_supaln) == 0:
#                ref_walk = target_walk_supaln
#            else:
#                ref_walk_start = target_walk_supaln.start
#                ref_walk = range(ref_walk_start, ref_walk_start + len_current_superseq_walk)
#
#            active_region_walks.append((ref_walk, subseq_walk))
#
#            if query_walk_supaln.stop >= target_walk_subaln.stop:
#                break
#    # pad with gap: right side
#    last_ref_walk = active_region_walks[-1][0]
#    last_subseq_walk = active_region_walks[-1][1]
#    if last_ref_walk.stop < len(sup_to_ref_aln.target):
#        added_ref_walk = range(last_ref_walk.stop, len(sup_to_ref_aln.target))
#        added_subseq_walk = range(last_subseq_walk.stop, last_subseq_walk.stop)
#        active_region_walks.append((added_ref_walk, added_subseq_walk))
#    # pad with gap: left
#    first_ref_walk = active_region_walks[0][0]
#    first_subseq_walk = active_region_walks[0][1]
#    if first_ref_walk.start > 0:
#        added_ref_walk = range(0, first_ref_walk.start)
#        added_subseq_walk = range(0, 0)
#        active_region_walks.insert(0, (added_ref_walk, added_subseq_walk))
#    # make result alignment
#    return Bio.Align.PairwiseAlignment(
#        sup_to_ref_aln.target,
#        subseq_row_spec['seq'],
#        alignhandler.walks_to_path(active_region_walks),
#        0
#    )


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


#def align_row_spec_to_ref(row_spec, ref_seq, ref_seq_reversed, aligner, raise_with_tie):
#    # set params
#    if aligner is None:
#        aligner = ALIGNER_MAIN
#    query = row_spec['seq']
#    reverse_align = (not row_spec['left_filled']) and row_spec['right_filled']
#    # run aligner
#    try:
#        if reverse_align:
#            alns = aligner.align(ref_seq_reversed, query[::-1])
#        else:
#            alns = aligner.align(ref_seq, query)
#    except Exception as exc:
#        raise Exception(f'Failed alignment:\nrow_spec: {row_spec}\nReference seq: {ref_seq}') from exc
#    # recover reversed alignment
#    if reverse_align:
#        alns = [alignhandler.reverse_alignment(x) for x in alns]
#
#    if len(alns) == 1:
#        aln = alns[0]
#        aln = alignhandler.amend_outer_insdel_both(aln)
#    else:
#        #alns = [alignhandler.amend_outer_insdel_both(x) for x in alns]
#        alns = list(
#            alignhandler.remove_identical_alignments(
#                alignhandler.amend_outer_insdel_both(x) for x in alns
#            )
#        )
#        try:
#            aln = alignhandler.alignment_tiebreaker(alns, raise_with_failure=raise_with_tie)
#        except alignhandler.AlignmentTieError as exc:
#            msg = f'Failed to break alignment tie. row_spec is:\n{row_spec}'
#            raise Exception(msg) from exc
#        
#    return aln


#def merge_two_cigartuples(cigartuples_left, cigartuples_right):
#    if len(cigartuples_left) == 0 or len(cigartuples_right) == 0:
#        return cigartuples_left + cigartuples_right
#    else:
#        if cigartuples_left[-1][0] == cigartuples_right[0][0]:
#            return (
#                cigartuples_left[:-1]
#                + [(cigartuples_left[-1][0], cigartuples_left[-1][1] + cigartuples_right[0][1])]
#                + cigartuples_right[1:]
#            )
#        else:
#            return cigartuples_left + cigartuples_right


#def realigned_cigar_sanity_check(original_read, realigned_read, realigned_cigartuples):
#    if len(original_read.query_sequence) != alignhandler.get_query_length(realigned_cigartuples):
#        raise RealignmentCigarError(
#            f'Read query sequence length and cigar length differs.\n'
#            f'Original read:\n{original_read.to_string()}\n'
#            f'Realigned read:\n{realigned_read.to_string()}\n'
#        )

    

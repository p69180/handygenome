import sys
import collections
import itertools
import functools
import warnings

import pysam
import Bio.Align

import importlib
top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
libvcfspec = importlib.import_module(".".join([top_package_name, "variant", "vcfspec"]))
libpileup = importlib.import_module(".".join([top_package_name, "read", "pileup"]))
alignhandler = importlib.import_module(".".join([top_package_name, "align", "alignhandler"]))
bameditor = importlib.import_module(".".join([top_package_name, "bameditor"]))
readhandler = importlib.import_module(".".join([top_package_name, "read", "readhandler"]))


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



###########################
# secure_inactive_padding #
###########################

def secure_inactive_padding_rightward(pileup, inactive_padding, extend_pileup_by=libpileup.DEFAULT_EXTEND_PILEUP_BY, extend_fetchedreads_by=libpileup.DEFAULT_EXTEND_FETCHEDREADS_BY, inplace=False):
    """Returns a modified pileup object with inactive region on the right margin"""
    # set initial_rightmost_active_pos0
    initial_active_positions = pileup.get_active_positions()
    if len(initial_active_positions) == 0:
        initial_rightmost_active_pos0 = pileup.start0 - 1
    else:
        initial_rightmost_active_pos0 = max(initial_active_positions)
    # initial extend
    while True:
        if pileup.end0 - (initial_rightmost_active_pos0 + 1) < inactive_padding:
            pileup.extend_rightward(width=extend_pileup_by, extend_fetchedreads_by=extend_fetchedreads_by)
        else:
            break
    # initial search loop
    active_info = pileup.get_active_info(start0=(initial_rightmost_active_pos0 + 1))
    found_result, farthest_inactive_pos0 = _find_inactive_run(active_info, inactive_padding)
    # subsequent search loop
    if not found_result:
        while True:
            pileup.extend_rightward(width=extend_pileup_by, extend_fetchedreads_by=extend_fetchedreads_by)
            active_info = pileup.get_active_info(
                start0=(pileup.end0 - extend_pileup_by + 1 - inactive_padding)
            )
            found_result, farthest_inactive_pos0 = _find_inactive_run(active_info, inactive_padding)
            if found_result:
                break
    # clip pileup
    if inplace:
        pileup.subset(pileup.start0, farthest_inactive_pos0 + 1, inplace=True)
    else:
        pileup_clipped = pileup.subset(pileup.start0, farthest_inactive_pos0 + 1, inplace=False)
        return pileup_clipped


def secure_inactive_padding_leftward(pileup, inactive_padding, extend_pileup_by=libpileup.DEFAULT_EXTEND_PILEUP_BY, extend_fetchedreads_by=libpileup.DEFAULT_EXTEND_FETCHEDREADS_BY, inplace=False):
    """Returns a modified pileup object with inactive region on the left margin"""
    # set initial_leftmost_active_pos0
    initial_active_positions = pileup.get_active_positions()
    if len(initial_active_positions) == 0:
        initial_leftmost_active_pos0 = pileup.end0
    else:
        initial_leftmost_active_pos0 = min(initial_active_positions)
    # initial extend
    while True:
        if initial_leftmost_active_pos0 - pileup.start0 < inactive_padding:
            pileup.extend_leftward(width=extend_pileup_by, extend_fetchedreads_by=extend_fetchedreads_by)
        else:
            break
    # initial search loop
    active_info = pileup.get_active_info(end0=initial_leftmost_active_pos0)[::-1]
    found_result, farthest_pos0 = _find_inactive_run(active_info, inactive_padding)
    # subsequent search loop
    if not found_result:
        while True:
            pileup.extend_leftward(width=extend_pileup_by, extend_fetchedreads_by=extend_fetchedreads_by)
            active_info = pileup.get_active_info(
                end0=(pileup.start0 + extend_pileup_by - 1 + inactive_padding)
            )[::-1]
            found_result, farthest_pos0 = _find_inactive_run(active_info, inactive_padding)
            if found_result:
                break
    # clip pileup
    if inplace:
        pileup.subset(farthest_pos0, pileup.end0, inplace=True)
    else:
        pileup_clipped = pileup.subset(farthest_pos0, pileup.end0, inplace=False)
        return pileup_clipped


def _find_inactive_run(active_info, inactive_padding):
    found_result = False
    farthest_inactive_pos0 = None
    for start_idx in range(0, len(active_info) - inactive_padding + 1):
        if not any(active_info[start_idx:(start_idx + inactive_padding)]):
            found_result = True
            farthest_inactive_pos0 = active_info.index[start_idx + inactive_padding -1]
            break
    return found_result, farthest_inactive_pos0


############################
# get_active_region_pileup #
############################

def get_active_region_pileup(
    bam,
    chrom,
    start0,
    end0,
    fasta,
    threshold=libpileup.DEFAULT_ACTIVE_THRESHOLD,
    inactive_padding=10,
    extend_fetchedreads_by=libpileup.DEFAULT_EXTEND_FETCHEDREADS_BY,
    extend_pileup_by=libpileup.DEFAULT_EXTEND_PILEUP_BY,
):
    """Returns:
        "active_range" and "pileup" are returned.
        When there is no active position around input range (start0, end0), 
            returned "active_range" is None. Pileup object, created as a 
            byproduct, is returned.
        Otherwise, "active_range" is python range object. Pileup object is
            clipped by "subset" method to fit "active_range", 
            then "set_row_specs" and "set_row_spec_groups" methods are run, 
            then returned.
    """

    def get_initial_fetch_range_old(start0, end0, extend_fetchedreads_by):
        query_range = range(start0, end0)
        diff_range_len = extend_fetchedreads_by - len(query_range)
        if diff_range_len > 0:
            pad = int(diff_range_len / 2)
            if diff_range_len % 2 == 0:
                initial_fetch_range = range(start0 - pad, end0 + pad)
            else:
                initial_fetch_range = range(start0 - pad, end0 + pad + 1)
        else:
            initial_fetch_range = query_range

        return initial_fetch_range

    def get_initial_fetch_range(start0, end0, inactive_padding):
        return range(start0 - inactive_padding * 3, end0 + inactive_padding * 3)

    def inspect_initial_range(
        chrom, start0, end0, fasta, bam, threshold, extend_fetchedreads_by, inactive_padding,
    ):
        initial_fetch_range = get_initial_fetch_range(start0, end0, inactive_padding)
        fetchedreads_initial = libpileup.FetchedReads.from_fetch(
            bam, chrom, initial_fetch_range.start, initial_fetch_range.stop, readfilter=readhandler.readfilter_pileup,
        )
        active_region_pileup = libpileup.get_pileup(
            chrom=chrom, 
            start0=start0, 
            end0=end0, 
            fetchedreads=fetchedreads_initial,
            fasta=fasta,
            active_threshold=threshold,
            truncate=True,
            as_array=False,
            return_range=False,
        )
        return active_region_pileup

    # main
    active_region_pileup = inspect_initial_range(
        chrom, start0, end0, fasta, bam, threshold, extend_fetchedreads_by, inactive_padding,
    )
    active_region_pileup.secure_inactive_padding_rightward(inactive_padding, extend_pileup_by, extend_fetchedreads_by)
    active_region_pileup.secure_inactive_padding_leftward(inactive_padding, extend_pileup_by, extend_fetchedreads_by)
    active_region_pileup.set_row_specs()
    active_region_pileup.set_row_spec_groups()

    active_poslist = active_region_pileup.get_active_positions()
    if len(active_poslist) == 0:
        active_range = None
    else:
        active_range = active_region_pileup.range

    return active_range, active_region_pileup


####################


def get_realigned_reads(active_region_pileup, active_range, ref_seq, ref_seq_reversed, aligner=ALIGNER_MAIN):
    realigned_read_ids = set()
    realigned_reads = list()
    for superseq_id, groupinfo in sorted(
        active_region_pileup.row_spec_groups.items(),
        key=(lambda x: x[1]['subseq_hits']), 
        reverse=True,
    ):
        realigned_read_ids.add(superseq_id)

        realigned_read_super, superseq_aln = realign_superseq(superseq_id, active_region_pileup, active_range, ref_seq, ref_seq_reversed, aligner)
        realigned_reads.append(realigned_read_super)
        active_region_pileup._alignment_store[superseq_id] = superseq_aln

        superseq_row_spec = active_region_pileup.row_specs[superseq_id]
        superseq = superseq_row_spec['seq']
        for subseq_uid in groupinfo['subseq_rowids']:
            if subseq_uid in realigned_read_ids:
                continue
            realigned_read_ids.add(subseq_uid)
            realigned_read_sub, subseq_aln = realign_subseq(subseq_uid, active_region_pileup, active_range, superseq, superseq_aln)
            realigned_reads.append(realigned_read_sub)
            
    return realigned_reads, realigned_read_ids


def get_vcfspecs(active_region_pileup, active_range, ref_seq, subseq_portion_threshold=0.05):
    vcfspec_list = list()
    subseq_hits_sum = sum(x['subseq_hits'] for x in active_region_pileup.row_spec_groups.values())
    for superseq_id, dic in active_region_pileup.row_spec_groups.items():
        if dic['subseq_hits'] / subseq_hits_sum > subseq_portion_threshold:
            superseq_aln = active_region_pileup._alignment_store[superseq_id]
            superseq_vcfspec_list = alignhandler.alignment_to_vcfspec(
                superseq_aln, 
                active_region_pileup.start0, 
                active_region_pileup.chrom, 
                active_region_pileup.fasta,
                strip_query_gaps=True,
            )
            if (superseq_vcfspec_list not in vcfspec_list) and (len(superseq_vcfspec_list) > 0):
                vcfspec_list.append(superseq_vcfspec_list)
            
    return vcfspec_list


def get_vcfspec_list_old(active_region_pileup, active_range, ref_seq, subseq_portion_threshold=0.05):
    # step 1: get selected_superseq_uids
    # selected_superseq_uids = list()
    vcfspec_list = list()
    subseq_hits_sum = sum(x['subseq_hits'] for x in active_region_pileup.row_spec_groups.values())
    for superseq_uid, dic in active_region_pileup.row_spec_groups.items():
        if dic['subseq_hits'] / subseq_hits_sum > subseq_portion_threshold:
            row_spec = active_region_pileup.row_specs[superseq_uid]
            if ref_seq == row_spec['seq']:
                continue
            vcfspec = libvcfspec.Vcfspec(active_region_pileup.chrom, active_range.start + 1, ref_seq, (row_spec['seq'],))
            vcfspec = vcfspec.normalize(active_region_pileup.fasta)
            vcfspec_list.append(vcfspec)
            
    return vcfspec_list


def write_realigned_bam(bam_path, realigned_reads, realigned_read_uids, active_region_pileup, padding_width=0):
    written_reads_range = range(
        active_region_pileup.start0 - padding_width, 
        active_region_pileup.end0 + padding_width, 
    )
    # extend fetchedreads if it does not contain the range of reads to be written
    active_region_pileup.fetchedreads.add_reads(active_region_pileup.chrom, written_reads_range.start, written_reads_range.stop)
    # write reads
    hdr = pysam.AlignmentHeader.from_references(
        reference_names=active_region_pileup.bam.references,
        reference_lengths=active_region_pileup.bam.lengths,
    )
    with pysam.AlignmentFile(bam_path, mode='wb', header=hdr) as in_bam:
        # write non-realigned reads
        for uid, read in active_region_pileup.fetchedreads.fetch(
            active_region_pileup.chrom,
            written_reads_range.start,
            written_reads_range.stop,
            with_uid=True,
        ):
            if uid not in realigned_read_uids:
                in_bam.write(read)
        # write realigned reads
        for read in realigned_reads:
            in_bam.write(read)
    # sort and index
    bameditor.sort_and_index(bam_path)


###################################
# helpers for get_realigned_reads #
###################################

#def realign_superseq_old(uid, active_region_pileup, active_range, ref_seq, ref_seq_reversed):
#    read = active_region_pileup.get_read(uid)
#    row_spec = active_region_pileup.row_specs[uid]
#    active_region_alignment = align_row_spec(row_spec, ref_seq, ref_seq_reversed)
#    realigned_read = active_aln_to_full_read(active_region_alignment, read, active_range)
#    return realigned_read, active_region_alignment


def realign_superseq(uid, active_region_pileup, active_range, ref_seq, ref_seq_reversed, aligner):
    # set params
    original_read = active_region_pileup.get_read(uid)
    original_cigar_split = alignhandler.split_cigar(original_read.cigartuples, original_read.reference_start, active_range)
    row_spec = active_region_pileup.row_specs[uid]
    empty_before_active_region = (len(original_cigar_split[0]) == 0)
    empty_after_active_region = (len(original_cigar_split[2]) == 0)
    # do alignment
    active_region_alignment = align_row_spec(row_spec, ref_seq, ref_seq_reversed, aligner, empty_before_active_region, empty_after_active_region)
    # get active region cigartuples
    active_region_cigartuples, active_region_offset = alignhandler.alignment_to_cigartuples(
        active_region_alignment,
        match_as_78=False, del_as_skip=False, 
        left_ins_as_clip=empty_before_active_region, right_ins_as_clip=empty_after_active_region,
        remove_left_del=empty_before_active_region, remove_right_del=empty_after_active_region,
    )
    # split cigartuples are merged into one
    realigned_cigartuples = merge_cigartuples_list(
        [original_cigar_split[0], active_region_cigartuples, original_cigar_split[2]]
    )
    # get new read reference_start
    if empty_before_active_region:
        new_reference_start = active_range.start + active_region_offset
    else:
        new_reference_start = original_read.reference_start
        
    realigned_read = original_read.__copy__()
    realigned_read.reference_start = new_reference_start
    realigned_read.cigartuples = realigned_cigartuples

    realigned_cigar_sanity_check(original_read, realigned_read, realigned_cigartuples)

    return realigned_read, active_region_alignment


#def realign_subseq_old(uid, active_region_pileup, active_range, superseq, superseq_reversed, superseq_alignment):
#    # set params
#    read = active_region_pileup.get_read(uid)
#    row_spec = active_region_pileup.row_specs[uid]
#        
#    sub_to_sup_aln = align_row_spec(row_spec, superseq, superseq_reversed)
#
#    # sanity check
#    active_region_cigartuples, active_region_offset = alignhandler.alignment_to_cigartuples(sub_to_sup_aln)
#    if set(x[0] for x in active_region_cigartuples) != alignhandler.CIGAROPS_BOTH:
#        raise Exception(
#            f'subseq-to-superseq alignment is not match-only.\n'
#            f'subseq row_spec: {row_spec}\n'
#            f'subseq-to-superseq alignment: {sub_to_sup_aln}'
#        )
#    
#    active_region_alignment = alignhandler.map_alignments(superseq_alignment, sub_to_sup_aln)   
#    # superseq_alignment.map(sub_to_sup_aln)
#    
#    realigned_read = active_aln_to_full_read(active_region_alignment, read, active_range)
#    return realigned_read, active_region_alignment


def realign_subseq(uid, active_region_pileup, active_range, superseq, superseq_alignment):
    # set params
    original_read = active_region_pileup.get_read(uid)
    original_cigar_split = alignhandler.split_cigar(original_read.cigartuples, original_read.reference_start, active_range)
    row_spec = active_region_pileup.row_specs[uid]
    empty_before_active_region = (len(original_cigar_split[0]) == 0)
    empty_after_active_region = (len(original_cigar_split[2]) == 0)

    # do initial alignment
    active_region_alignment = align_subseq_to_ref(row_spec, superseq_alignment)

    # get active region cigartuples
    active_region_cigartuples, active_region_offset = alignhandler.alignment_to_cigartuples(
        active_region_alignment,
        match_as_78=False, del_as_skip=False, 
        left_ins_as_clip=empty_before_active_region, right_ins_as_clip=empty_after_active_region,
        remove_left_del=empty_before_active_region, remove_right_del=empty_after_active_region,
    )
    # merge split cigartuples
    realigned_cigartuples = merge_cigartuples_list(
        [original_cigar_split[0], active_region_cigartuples, original_cigar_split[2]]
    )
    # get new read reference_start
    if empty_before_active_region:
        new_reference_start = active_range.start + active_region_offset
    else:
        new_reference_start = original_read.reference_start
        
    realigned_read = original_read.__copy__()
    realigned_read.reference_start = new_reference_start
    realigned_read.cigartuples = realigned_cigartuples

    realigned_cigar_sanity_check(original_read, realigned_read, realigned_cigartuples)

    return realigned_read, active_region_alignment


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
    #superseq_alignment_walks = list(alignhandler.path_to_walks(sup_to_ref_aln.path))

    found_initial_hit = False
    for target_walk_supaln, query_walk_supaln in alignhandler.path_to_walks(sup_to_ref_aln.path):
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
    path = alignhandler.walks_to_path(active_region_walks)
    return Bio.Align.PairwiseAlignment(
        sup_to_ref_aln.target,
        subseq_row_spec['seq'],
        path,
        0
    )


#######################################
# amenders: marginal insdel amendment #
#######################################

def amender_left_right(alignment):
    return alignhandler.amend_outer_insdel_right(
        alignhandler.amend_outer_insdel_left(alignment)
    )


def amender_left(alignment):
    return alignhandler.amend_outer_insdel_left(alignment)


def amender_right(alignment):
    return alignhandler.amend_outer_insdel_right(alignment)


def amender_none(alignment):
    return alignment


def get_amender(empty_before_active_region, empty_after_active_region):
    if empty_before_active_region:
        if empty_after_active_region:
            return amender_left_right
        else:
            return amender_left
    else:
        if empty_after_active_region:
            return amender_right
        else:
            return amender_none


#######################################


def align_row_spec(row_spec, target, target_reversed, aligner, empty_before_active_region, empty_after_active_region):
    # set params
    query = row_spec['seq']
    reverse_align = (not row_spec['left_filled']) and row_spec['right_filled']
    amender = get_amender(empty_before_active_region, empty_after_active_region)
    # run aligner
    try:
        if reverse_align:
            alns = aligner.align(target_reversed, query[::-1])
        else:
            alns = aligner.align(target, query)
    except Exception as exc:
        raise Exception(f'Failed alignment:\nrow_spec: {row_spec}\ntarget: {target}') from exc
    # recover reversed alignment
    if reverse_align:
        alns = [alignhandler.reverse_alignment(x) for x in alns]

    if len(alns) == 1:
        aln = alns[0]
        aln = amender(aln)
    else:
        alns = [amender(x) for x in alns]
        try:
            aln = alignhandler.alignment_tiebreaker(alns, raise_with_failure=True)
        except alignhandler.AlignmentTieError as exc:
            msg = 'Failed to break alignment tie. row_spec is:\n{row_spec}'
            raise Exception(msg) from exc
        
    return aln


#def active_aln_to_full_read(active_region_alignment, read, active_range):
#    original_cigar_split = alignhandler.split_cigar(read.cigartuples, read.reference_start, active_range)
#
#    active_region_cigartuples, active_region_offset = alignhandler.alignment_to_cigartuples(active_region_alignment)
#    
#    realigned_cigartuples = merge_cigartuples_list(
#        [original_cigar_split[0], active_region_cigartuples, original_cigar_split[2]]
#    )
#    
#    if len(original_cigar_split[0]) > 0:
#        #assert active_region_offset == 0
#        new_reference_start = read.reference_start
#    else:
#        new_reference_start = active_range.start + active_region_offset
#        
#    realigned_read = read.__copy__()
#    realigned_read.reference_start = new_reference_start
#    realigned_read.cigartuples = realigned_cigartuples
#
#    # sanity check
#    # if len(read.query_sequence) != alignhandler.get_query_length(realigned_cigartuples):
#    #     raise Exception(
#    #         f'Read query sequence length and cigar length differs:\n'
#    #         f'Original read: {read.to_string()}\n'
#    #         f'Realigned read: {realigned_read.to_string()}\n'
#    #     )
#
#    return realigned_read


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

    
def merge_cigartuples_list(cigartuples_list):
    return list(functools.reduce(merge_two_cigartuples, cigartuples_list))


def realigned_cigar_sanity_check(original_read, realigned_read, realigned_cigartuples):
    if len(original_read.query_sequence) != alignhandler.get_query_length(realigned_cigartuples):
        raise Exception(
            f'Read query sequence length and cigar length differs.\n'
            f'Original read:\n{original_read.to_string()}\n'
            f'Realigned read:\n{realigned_read.to_string()}\n'
        )

    
###################################

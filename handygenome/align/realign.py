import sys
import collections
import itertools
import functools
import warnings

import pysam
import Bio.Align
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


DEFAULT_EXTEND_PILEUP_BY = 20
DEFAULT_INACTIVE_PADDING = 10
DEFAULT_ACTIVE_THRESHOLD = 0.05
DEFAULT_VCFSPEC_EXTRACTION_THRESHOLD = 0.05
DEFAULT_VCFSPEC_CONCAT_DISTANCE = None


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



############################
# inactive padding related #
############################

def secure_inactive_padding_rightward(pileup, inactive_padding, extend_pileup_by=libpileup.DEFAULT_EXTEND_PILEUP_BY, inplace=False):
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
            pileup.extend_rightward(width=extend_pileup_by)
        else:
            break
    # initial search loop
    active_info = pileup.get_active_info(start0=(initial_rightmost_active_pos0 + 1))
    found_result, farthest_inactive_pos0 = _find_inactive_run(active_info, inactive_padding)
    # subsequent search loop
    if not found_result:
        while True:
            pileup.extend_rightward(width=extend_pileup_by)
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


def secure_inactive_padding_leftward(pileup, inactive_padding, extend_pileup_by=libpileup.DEFAULT_EXTEND_PILEUP_BY, inplace=False):
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
            pileup.extend_leftward(width=extend_pileup_by)
        else:
            break
    # initial search loop
    active_info = pileup.get_active_info(end0=initial_leftmost_active_pos0)[::-1]
    found_result, farthest_pos0 = _find_inactive_run(active_info, inactive_padding)
    # subsequent search loop
    if not found_result:
        while True:
            pileup.extend_leftward(width=extend_pileup_by)
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

def secure_margin_with_vcfspecs(active_region_pileup, factor_vcfspec_range, factor_repeat_area):
    candidate_start0s = set()
    candidate_end0s = set()
    for vcfspec in itertools.chain.from_iterable(active_region_pileup.vcfspecs):
        # padding with (vcfspec REF_range length) * factor
        vcfspec_range = vcfspec.REF_range0
        padding = int(len(vcfspec_range) * factor_vcfspec_range)
        candidate_start0s.add(vcfspec_range.start - padding)
        candidate_end0s.add(vcfspec_range.stop + padding)
        # padding considering repetitive reference region
        equivalents = vcfspec.get_equivalents()
        repeat_area_start0 = equivalents[0].pos0
        repeat_area_end0 = equivalents[-1].end0
        padding = int((repeat_area_end0 - repeat_area_start0) * factor_repeat_area)
        candidate_start0s.add(repeat_area_start0 - padding)
        candidate_end0s.add(repeat_area_end0 + padding)

    if len(candidate_start0s) == 0:
        left_edited = False
    else:
        desired_start0 = min(candidate_start0s)
        start_margin = desired_start0 - active_region_pileup.start0
        if start_margin < 0:
            left_edited = True
            active_region_pileup.extend_leftward(-start_margin)
        else:
            left_edited = False

    if len(candidate_end0s) == 0:
        right_edited = False
    else:
        desired_end0 = max(candidate_end0s)
        end_margin = active_region_pileup.end0 - desired_end0
        if end_margin < 0:
            right_edited = True
            active_region_pileup.extend_rightward(-end_margin)
        else:
            right_edited = False

    edited = left_edited or right_edited

    return edited


def augment_margins(active_region_pileup, inactive_padding, extend_pileup_by, factor_vcfspec_range, factor_repeat_area, aligner, allele_portion_threshold, vcfspec_concat_dist):
    edited = False
    while True:
        if active_region_pileup.check_inactive_margin_right(inactive_padding):
            right_inactive_margin_edited = False
        else:
            right_inactive_margin_edited = True
            edited = True
            active_region_pileup.secure_inactive_padding_rightward(inactive_padding, extend_pileup_by=extend_pileup_by)

        if active_region_pileup.check_inactive_margin_left(inactive_padding):
            left_inactive_margin_edited = False
        else:
            left_inactive_margin_edited = True
            edited = True
            active_region_pileup.secure_inactive_padding_leftward(inactive_padding, extend_pileup_by=extend_pileup_by)

        # preparation for secure_margin_with_vcfspecs
        active_region_pileup.set_row_specs()
        active_region_pileup.set_row_spec_groups()
        active_region_pileup.save_superseq_alignments(aligner=aligner)
        active_region_pileup.set_vcfspecs(allele_portion_threshold=allele_portion_threshold, concat_dist_le=vcfspec_concat_dist)

        vcfspec_margin_edited = secure_margin_with_vcfspecs(active_region_pileup, factor_vcfspec_range, factor_repeat_area)

        if vcfspec_margin_edited:
            edited = True

        if not vcfspec_margin_edited:
            break

    return edited


def augment_margins_multisample(active_region_mspileup, inactive_padding, extend_pileup_by, factor_vcfspec_range, factor_repeat_area, aligner, allele_portion_threshold, vcfspec_concat_dist):
    while True:
        active_region_mspileup.equalize_margins()
        augmented = dict()
        for sampleid, pileup in active_region_mspileup.pileup_dict.items():
            augmented[sampleid] = augment_margins(pileup, inactive_padding, extend_pileup_by, factor_vcfspec_range, factor_repeat_area, aligner, allele_portion_threshold, vcfspec_concat_dist)
        if not any(augmented.values()):
            break


def get_active_region_pileup(
    chrom,
    start0,
    end0,
    bam,
    fasta,
    active_threshold=DEFAULT_ACTIVE_THRESHOLD,
    inactive_padding=DEFAULT_INACTIVE_PADDING,
    extend_pileup_by=libpileup.DEFAULT_EXTEND_PILEUP_BY,
    augment_factor_vcfspec_range=1.5,
    augment_factor_repeat_area=1.5,
    aligner=None,
    allele_portion_threshold=None,
    vcfspec_concat_dist=None,
):
    """Returns:
        "active_range" and "pileup" are returned.
        When there is no active position around input range (start0, end0), 
            returned "active_range" is None. Pileup object, created as a 
            byproduct, is returned.
    """
    # main
    active_region_pileup = libpileup.get_pileup(
        chrom=chrom, 
        start0=start0, 
        end0=end0, 
        bam=bam,
        fasta=fasta,
        active_threshold=active_threshold,
        truncate=True,
        as_array=False,
        return_range=False,
    )
    edited = augment_margins(
        active_region_pileup, 
        inactive_padding, 
        extend_pileup_by, 
        factor_vcfspec_range=augment_factor_vcfspec_range, 
        factor_repeat_area=augment_factor_repeat_area,
        aligner=aligner,
        allele_portion_threshold=allele_portion_threshold,
        vcfspec_concat_dist=vcfspec_concat_dist,
    )
        # now set_row_specs(), set_row_spec_groups(), save_superseq_alignments(), set_vcfspecs() are done

    active_poslist = active_region_pileup.get_active_positions()
    if len(active_poslist) == 0:
        active_range = None
    else:
        active_range = active_region_pileup.range0

    return active_range, active_region_pileup


def get_active_region_pileup_multisample(
    chrom,
    start0,
    end0,
    bam_dict,
    fasta,
    active_threshold=DEFAULT_ACTIVE_THRESHOLD,
    inactive_padding=DEFAULT_INACTIVE_PADDING,
    extend_pileup_by=libpileup.DEFAULT_EXTEND_PILEUP_BY,
    augment_factor_vcfspec_range=1.5,
    augment_factor_repeat_area=1.5,
    aligner=None,
    allele_portion_threshold=None,
    vcfspec_concat_dist=None,
):
    # initialize pileup_dict
    pileup_dict = dict()
    for sampleid, bam in bam_dict.items():
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


# get aligned reads

def get_realigned_reads(active_region_pileup):
    """Must be run after 'save_superseq_alignments' method"""
    return get_realigned_reads_helper(
        active_region_pileup, 
        active_region_pileup.row_spec_groups, 
        active_region_pileup._alignment_cache,
    )


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

def get_vcfspecs_from_pileup(active_region_pileup, allele_portion_threshold=None, concat_dist_le=None):
    """Must be run after row_spec_groups is set and save_superseq_alignments have been run"""
    if allele_portion_threshold is None:
        allele_portion_threshold = DEFAULT_VCFSPEC_EXTRACTION_THRESHOLD

    naive_vcfspecs = set()

    chrom = active_region_pileup.chrom
    active_region_start0 = active_region_pileup.start0
    fasta = active_region_pileup.fasta
    subseq_hits_sum = sum(x['subseq_hits'] for x in active_region_pileup.row_spec_groups.values())

    for superseq_rowid, groupinfo in active_region_pileup.row_spec_groups.items():
        allele_portion = groupinfo['subseq_hits'] / subseq_hits_sum
        if allele_portion < allele_portion_threshold:
            continue

        superseq_row_spec = active_region_pileup.row_specs[superseq_rowid]
        superseq_alignment = active_region_pileup._alignment_cache[superseq_rowid]

        superseq_vcfspecs = get_vcfspecs_helper(superseq_row_spec, superseq_alignment, chrom, active_region_start0, fasta, concat_dist_le)
        active_region_pileup._vcfspec_cache[superseq_rowid] = superseq_vcfspecs
        if len(superseq_vcfspecs) > 0:
            naive_vcfspecs.add(superseq_vcfspecs)
            
    return naive_vcfspecs


def get_vcfspecs_from_pileup_multisample(active_region_mspileup, allele_portion_threshold=None, concat_dist_le=None):
    """Must be run after 'active_region_mspileup' has executed 
    'divide_row_spec_groups' and 'get_realigned_reads' methods.
        'divide_row_spec_groups' is needed for setting 'row_spec_groups' for
            each single sample Pileup.
        'get_realigned_reads' is needed for initiating '_alignment_cache' attribute.
    """
    if allele_portion_threshold is None:
        allele_portion_threshold = DEFAULT_VCFSPEC_EXTRACTION_THRESHOLD

    naive_vcfspecs = set()
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
                superseq_vcfspecs = get_vcfspecs_helper(superseq_row_spec, superseq_alignment, chrom, active_region_start0, fasta, concat_dist_le)
            if len(superseq_vcfspecs) > 0:
                naive_vcfspecs.add(superseq_vcfspecs)
            
    return naive_vcfspecs


def get_vcfspecs_helper(superseq_row_spec, superseq_alignment, chrom, active_region_start0, fasta, concat_dist_le):
    if concat_dist_le is None:
        concat_dist_le = DEFAULT_VCFSPEC_CONCAT_DISTANCE

    lstrip_query_gaps = not superseq_row_spec["left_filled"]
    rstrip_query_gaps = not superseq_row_spec["right_filled"]
    superseq_vcfspecs = alignhandler.alignment_to_vcfspec(
        alignment=superseq_alignment, 
        target_start0=active_region_start0, 
        chrom=chrom, 
        fasta=fasta,
        lstrip_query_gaps=lstrip_query_gaps,
        rstrip_query_gaps=rstrip_query_gaps,
    )
    if len(superseq_vcfspecs) > 0:
        superseq_vcfspecs = libvcfspec.concat_list(superseq_vcfspecs, distance_le=concat_dist_le)
    return superseq_vcfspecs


# process naive vcfspecs into final product

def process_naive_vcfspecs(naive_vcfspecs):
    vcfspecs_dict = {
        x.get_id(): x for x in 
        itertools.chain.from_iterable(naive_vcfspecs)
    }

    grlist = list()
    for idx, row in enumerate(naive_vcfspecs):
        gr = pr.concat(x.to_gr() for x in row)
        gr = gr.assign('row_index', lambda df: pd.Series(itertools.repeat(idx, df.shape[0])))
        grlist.append(gr)

    result = list()
    clustered_gr = pr.concat(grlist).cluster()
    for cluster_idx in set(clustered_gr.Cluster):
        concatenated_vcfspecs = list()
        gr_subset = clustered_gr[clustered_gr.Cluster == cluster_idx]
        for row_index in set(gr_subset.row_index):
            current_vcfspec_list = [
                vcfspecs_dict[vcfspec_id]
                for vcfspec_id in gr_subset[gr_subset.row_index == row_index].Id
            ]
            vcfspec_concat = libvcfspec.concat_list(current_vcfspec_list, distance_le=None)[0]
            concatenated_vcfspecs.append(vcfspec_concat)
        merged_vcfspec = functools.reduce(libvcfspec.merge, concatenated_vcfspecs)
        result.append(merged_vcfspec)

    result.sort(key=(lambda x: x.pos))

    return result


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

def save_superseq_alignments(pileup, aligner=None):
    if aligner is None:
        aligner = ALIGNER_MAIN

    pileup._alignment_cache = dict()
    ref_seq = pileup.get_ref_seq()
    ref_seq_reversed = ref_seq[::-1]
    for superseq_rowid in pileup.row_spec_groups.keys():
        row_spec = pileup.row_specs[superseq_rowid]
        pileup._alignment_cache[superseq_rowid] = align_row_spec_to_ref(row_spec, ref_seq, ref_seq_reversed, aligner)


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
    realigned_read = align_row_helper(row_id, active_region_pileup, active_range, active_region_alignment)

    return realigned_read


def align_row_helper(row_id, active_region_pileup, active_range, active_region_alignment):
    # set params
    original_read = active_region_pileup.read_cache[row_id]
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


def align_row_spec_to_ref(row_spec, ref_seq, ref_seq_reversed, aligner):
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
            aln = alignhandler.alignment_tiebreaker(alns, raise_with_failure=True)
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
        raise Exception(
            f'Read query sequence length and cigar length differs.\n'
            f'Original read:\n{original_read.to_string()}\n'
            f'Realigned read:\n{realigned_read.to_string()}\n'
        )

    

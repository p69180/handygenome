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
fetchcache = importlib.import_module(".".join([top_package_name, "read", "fetchcache"]))


DEFAULT_INACTIVE_PADDING = 10
DEFAULT_ACTIVE_THRESHOLD = 0.05
DEFAULT_VCFSPEC_EXTRACTION_THRESHOLD = 0.05


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


###########################
# pileup range correction #
###########################

def augment_pileup_range(active_region_pileup, factor_vcfspec_range=2, factor_repeat_area=1.5):
    active_region_pileup.set_row_specs()
    active_region_pileup.set_row_spec_groups()
    pileup_vcfspecs = get_vcfspecs_from_pileup(active_region_pileup)

    candidate_start0s = list()
    candidate_end0s = list()
    for vcfspec in itertools.chain.from_iterable(pileup_vcfspecs):
        # padding with (vcfspec REF_range length) * factor
        vcfspec_range = vcfspec.REF_range0
        spanning_length = len(vcfspec_range)
        padding = int(len(vcfspec_range) * factor_vcfspec_range)
        candidate_start0s.append(vcfspec_range.start - padding)
        candidate_end0s.append(vcfspec_range.stop + padding)
        # padding considering repetitive reference region
        equivalents = vcfspec.get_equivalents(active_region_pileup.fasta)


############################
# get_active_region_pileup #
############################

def get_active_region_pileup(
    chrom,
    start0,
    end0,
    bam,
    fasta,
    active_threshold=DEFAULT_ACTIVE_THRESHOLD,
    inactive_padding=DEFAULT_INACTIVE_PADDING,
    fetchedreads_addition_padding=fetchcache.DEFAULT_FETCHEDREADS_ADDITION_PADDING,
    extend_pileup_by=libpileup.DEFAULT_EXTEND_PILEUP_BY,
    set_row_specs=True,
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
        chrom, start0, end0, fasta, bam, active_threshold, inactive_padding, fetchedreads_addition_padding,
    ):
        initial_fetch_range = get_initial_fetch_range(start0, end0, inactive_padding)
        fetchedreads_initial = fetchcache.FetchedReads.from_fetch(
            bam, chrom, initial_fetch_range.start, initial_fetch_range.stop, 
            readfilter=readhandler.readfilter_pileup, 
            addition_padding=fetchedreads_addition_padding,
        )
        active_region_pileup = libpileup.get_pileup(
            chrom=chrom, 
            start0=start0, 
            end0=end0, 
            fetchedreads=fetchedreads_initial,
            fasta=fasta,
            active_threshold=active_threshold,
            truncate=True,
            as_array=False,
            return_range=False,
        )
        return active_region_pileup

#    def augment_inactive_padding(active_region_pileup, inactive_padding, extend_pileup_by, extend_fetchedreads_by, factor=2):
#        active_region_length = len(active_region_pileup.range0) - (2 * inactive_padding)
#        if inactive_padding < factor * active_region_length:
#            new_inactive_padding = factor * active_region_length
#            active_region_pileup.secure_inactive_padding_rightward(new_inactive_padding, extend_pileup_by, extend_fetchedreads_by)
#            active_region_pileup.secure_inactive_padding_leftward(new_inactive_padding, extend_pileup_by, extend_fetchedreads_by)

    def extend_pileup_by_active_area_length(active_region_pileup, inactive_padding, extend_pileup_by, extend_fetchedreads_by, factor=1.3):
        active_area_end0 = active_region_pileup.end0 - inactive_padding
        active_area_start0 = active_region_pileup.start0 + inactive_padding
        active_area_length = active_area_end0 - active_area_start0
        # correct right margin
        desired_pileup_end0 = active_area_end0 + int(active_area_length * factor)
        if active_region_pileup.end0 < desired_pileup_end0:
            active_region_pileup.extend_rightward(width=(desired_pileup_end0 - active_region_pileup.end0), extend_fetchedreads_by=extend_fetchedreads_by)
        # correct left margin
        desired_pileup_start0 = active_area_start0 - int(active_area_length * factor)
        if desired_pileup_start0 < active_region_pileup.start0:
            active_region_pileup.extend_leftward(width=(active_region_pileup.start0 - desired_pileup_start0), extend_fetchedreads_by=extend_fetchedreads_by)
        # secure inactive margin
        active_region_pileup.secure_inactive_padding(inactive_padding, extend_pileup_by, extend_fetchedreads_by)

    def extend_pileup_by_reference_repeat(active_region_pileup, extend_pileup_by, extend_fetchedreads_by, factor=1.3):
        pass

    # main
    active_region_pileup = inspect_initial_range(
        chrom, start0, end0, fasta, bam, active_threshold, inactive_padding, fetchedreads_addition_padding,
    )
    active_region_pileup.secure_inactive_padding(inactive_padding, extend_pileup_by)

    if set_row_specs:
        active_region_pileup.set_row_specs()
        active_region_pileup.set_row_spec_groups()

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
    fetchedreads_addition_padding=fetchcache.DEFAULT_FETCHEDREADS_ADDITION_PADDING,
    extend_pileup_by=libpileup.DEFAULT_EXTEND_PILEUP_BY,
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
            fetchedreads_addition_padding=fetchedreads_addition_padding,
            extend_pileup_by=extend_pileup_by,
            set_row_specs=False,
        )
        pileup_dict[sampleid] = active_region_pileup

    # create MultisamplePileup object and postprocess
    active_region_mspileup = libpileup.MultisamplePileup(pileup_dict)
    active_region_mspileup.secure_inactive_padding(inactive_padding, extend_pileup_by=extend_pileup_by)
    active_region_mspileup.set_df()
    active_region_mspileup.set_row_specs()
    active_region_mspileup.set_row_spec_groups()
    active_region_mspileup.divide_row_spec_groups()

    # result values
    active_poslist = active_region_mspileup.get_active_positions()
    if len(active_poslist) == 0:
        active_range = None
    else:
        active_range = active_region_mspileup.range0

    return active_range, active_region_mspileup


####################


def get_realigned_reads(active_region_pileup, active_range, ref_seq, ref_seq_reversed, aligner):
    realigned_reads = dict()
    for superseq_rowid, groupinfo in sorted(
        active_region_pileup.row_spec_groups.items(),
        key=(lambda x: x[1]['subseq_hits']), 
        reverse=True,
    ):
        realigned_read_super, superseq_aln = realign_superseq(superseq_rowid, active_region_pileup, active_range, ref_seq, ref_seq_reversed, aligner)
        realigned_reads[superseq_rowid] = realigned_read_super
        active_region_pileup._alignment_cache[superseq_rowid] = superseq_aln

        superseq_row_spec = active_region_pileup.row_specs[superseq_rowid]
        superseq = superseq_row_spec['seq']
        for subseq_rowid in groupinfo['subseq_rowids']:
            if subseq_rowid in realigned_reads.keys():
                continue
            realigned_read_sub, subseq_aln = realign_subseq(subseq_rowid, active_region_pileup, active_range, superseq, superseq_aln)
            realigned_reads[subseq_rowid] = realigned_read_sub
            
    return realigned_reads


def get_vcfspecs_from_pileup(active_region_pileup, allele_portion_threshold=DEFAULT_VCFSPEC_EXTRACTION_THRESHOLD):
    """Must be run after row_spec_groups is set"""
    # prepare alignments for superseq rows
    unprepared_superseq_rowids = set(active_region_pileup.row_spec_groups.keys()).difference(
        active_region_pileup._alignment_cache.keys()
    )
    active_region_pileup.save_alignments(unprepared_superseq_rowids)
    # main
    vcfspecs = list()
    subseq_hits_sum = sum(x['subseq_hits'] for x in active_region_pileup.row_spec_groups.values())
    for superseq_rowid, groupinfo in active_region_pileup.row_spec_groups.items():
        allele_portion = groupinfo['subseq_hits'] / subseq_hits_sum
        if allele_portion < allele_portion_threshold:
            continue

        superseq_row_spec = active_region_pileup.row_specs[superseq_rowid]
        lstrip_query_gaps = not superseq_row_spec["left_filled"]
        rstrip_query_gaps = not superseq_row_spec["right_filled"]
        superseq_aln = active_region_pileup._alignment_cache[superseq_rowid]

        superseq_vcfspec_tuple = alignhandler.alignment_to_vcfspec(
            alignment=superseq_aln, 
            target_start0=active_region_pileup.start0, 
            chrom=active_region_pileup.chrom, 
            fasta=active_region_pileup.fasta,
            lstrip_query_gaps=lstrip_query_gaps,
            rstrip_query_gaps=rstrip_query_gaps,
        )
        if len(superseq_vcfspec_tuple) > 0:
            vcfspecs.append(superseq_vcfspec_tuple)
            
    return vcfspecs


def get_vcfspecs_multisample_pileup(active_region_mspileup, allele_portion_threshold):
    """Must be run after 'active_region_mspileup' has executed 
    'divide_row_spec_groups' and 'get_realigned_reads' methods.
        'divide_row_spec_groups' is needed for setting 'row_spec_groups' for
            each single sample Pileup.
        'get_realigned_reads' is needed for initiating '_alignment_cache' attribute.
    """
    vcfspecs_by_sample = dict()
    calculated_vcfspecs_cache = dict()
    #{sampleid: list() for sampleid in active_region_mspileup.pileup_dict.keys()}
    for sampleid, pileup in active_region_mspileup.pileup_dict.items():
        vcfspecs_onesample = set()
        subseq_hits_sum = sum(x['subseq_hits'] for x in pileup.row_spec_groups.values())

        for superseq_rowid, groupinfo in pileup.row_spec_groups.items():
            allele_portion = groupinfo['subseq_hits'] / subseq_hits_sum
            if allele_portion < allele_portion_threshold:
                continue

            if superseq_rowid in calculated_vcfspecs_cache.keys():
                superseq_vcfspec_tuple = calculated_vcfspecs_cache[superseq_rowid]
            else:
                superseq_row_spec = active_region_mspileup.row_specs[superseq_rowid]
                lstrip_query_gaps = not superseq_row_spec["left_filled"]
                rstrip_query_gaps = not superseq_row_spec["right_filled"]
                superseq_aln = active_region_mspileup._alignment_cache[superseq_rowid]
                superseq_vcfspec_tuple = alignhandler.alignment_to_vcfspec(
                    alignment=superseq_aln, 
                    target_start0=active_region_mspileup.start0, 
                    chrom=active_region_mspileup.chrom, 
                    fasta=active_region_mspileup.fasta,
                    lstrip_query_gaps=lstrip_query_gaps,
                    rstrip_query_gaps=rstrip_query_gaps,
                )
                calculated_vcfspecs_cache[superseq_rowid] = superseq_vcfspec_tuple

            if len(superseq_vcfspec_tuple) > 0:
                vcfspecs_onesample.add(superseq_vcfspec_tuple)

        vcfspecs_by_sample[sampleid] = vcfspecs_onesample
    # final result
    vcfspecs_allsamples = set(itertools.chain.from_iterable(vcfspecs_by_sample.values()))
            
    return vcfspecs_allsamples


#def get_vcfspec_list_old(active_region_pileup, active_range, ref_seq, allele_portion_threshold=0.05):
#    # step 1: get selected_superseq_uids
#    # selected_superseq_uids = list()
#    vcfspec_list = list()
#    subseq_hits_sum = sum(x['subseq_hits'] for x in active_region_pileup.row_spec_groups.values())
#    for superseq_uid, dic in active_region_pileup.row_spec_groups.items():
#        if dic['subseq_hits'] / subseq_hits_sum > allele_portion_threshold:
#            row_spec = active_region_pileup.row_specs[superseq_uid]
#            if ref_seq == row_spec['seq']:
#                continue
#            vcfspec = libvcfspec.Vcfspec(active_region_pileup.chrom, active_range.start + 1, ref_seq, (row_spec['seq'],))
#            vcfspec = vcfspec.normalize(active_region_pileup.fasta)
#            vcfspec_list.append(vcfspec)
#            
#    return vcfspec_list


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
        for uid, read in active_region_pileup.fetchedreads.fetch(
            active_region_pileup.chrom,
            written_reads_range.start,
            written_reads_range.stop,
            with_uid=True,
        ):
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

#def realign_superseq_old(uid, active_region_pileup, active_range, ref_seq, ref_seq_reversed):
#    read = active_region_pileup.get_read(uid)
#    row_spec = active_region_pileup.row_specs[uid]
#    active_region_alignment = align_row_spec(row_spec, ref_seq, ref_seq_reversed)
#    realigned_read = active_aln_to_full_read(active_region_alignment, read, active_range)
#    return realigned_read, active_region_alignment


def get_alignment_from_rowid(row_id, active_region_pileup, ref_seq, ref_seq_reversed, aligner):
    row_spec = active_region_pileup.row_specs[row_id]
    return align_row_spec(row_spec, ref_seq, ref_seq_reversed, aligner)


def realign_superseq(row_id, active_region_pileup, active_range, ref_seq, ref_seq_reversed, aligner):
    # set params
    original_read = active_region_pileup.get_read(row_id)
    original_cigar_split = alignhandler.split_cigar(original_read.cigartuples, original_read.reference_start, active_range)
    row_spec = active_region_pileup.row_specs[row_id]
    empty_before_active_region = (len(original_cigar_split[0]) == 0)
    empty_after_active_region = (len(original_cigar_split[2]) == 0)
    # do alignment
    active_region_alignment = align_row_spec(row_spec, ref_seq, ref_seq_reversed, aligner)
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
    readhandler.set_NMMD(realigned_read, active_region_pileup.fasta)

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


def realign_subseq(row_id, active_region_pileup, active_range, superseq, superseq_alignment):
    # set params
    original_read = active_region_pileup.get_read(row_id)
    original_cigar_split = alignhandler.split_cigar(original_read.cigartuples, original_read.reference_start, active_range)
    row_spec = active_region_pileup.row_specs[row_id]
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
    readhandler.set_NMMD(realigned_read, active_region_pileup.fasta)

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


def align_row_spec(row_spec, target, target_reversed, aligner):
    # set params
    query = row_spec['seq']
    reverse_align = (not row_spec['left_filled']) and row_spec['right_filled']
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
        aln = amender_left_right(aln)
    else:
        alns = [amender_left_right(x) for x in alns]
        try:
            aln = alignhandler.alignment_tiebreaker(alns, raise_with_failure=True)
        except alignhandler.AlignmentTieError as exc:
            msg = f'Failed to break alignment tie. row_spec is:\n{row_spec}'
            raise Exception(msg) from exc
        
    return aln


def align_row_spec_old(row_spec, target, target_reversed, aligner, empty_before_active_region, empty_after_active_region):
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
            msg = f'Failed to break alignment tie. row_spec is:\n{row_spec}'
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

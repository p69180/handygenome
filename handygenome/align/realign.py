import collections
import itertools

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
libpileup = importlib.import_module(".".join([top_package_name, "readplus", "pileup"]))


# get_active_range
def active_info_generator(bam, chrom, start0, end0, fasta, threshold, reverse=False):
    pileup = libpileup.get_pileup(
        bam, chrom, start0, end0, as_array=False, truncate=True
    )
    if reverse:
        pileup_iterator = pileup.df.iloc[:, ::-1].iteritems()
    else:
        pileup_iterator = pileup.df.iteritems()

    ref_seq = fasta.fetch(chrom, start0, end0)
    if reverse:
        ref_seq = ref_seq[::-1]

    for (pos0, col), ref_base in zip(pileup_iterator, ref_seq):
        max_allele, max_vaf = libpileup.get_max_vaf_from_pileupcol(
            col, with_counts=False
        )
        if max_vaf >= threshold:
            is_active = not (max_allele == ref_base)
        else:
            is_active = True

        yield (pos0, is_active)


def check_active_col(
    bam, chrom, current_pos0, fasta, threshold, extend_pileup_by, leftward=False
):
    if leftward:
        end0 = current_pos0 + 1
        start0 = end0 - extend_pileup_by
        reverse = True

        def update_start_end(start0, end0, extend_pileup_by):
            new_end0 = start0
            new_start0 = new_end0 - extend_pileup_by
            return new_start0, new_end0

    else:
        start0 = current_pos0
        end0 = start0 + extend_pileup_by
        reverse = False

        def update_start_end(start0, end0, extend_pileup_by):
            new_start0 = end0
            new_end0 = new_start0 + extend_pileup_by
            return new_start0, new_end0

    while True:
        for pos0, is_active in active_info_generator(
            bam,
            chrom,
            start0,
            end0,
            fasta,
            threshold,
            reverse,
        ):
            yield pos0, is_active

        start0, end0 = update_start_end(start0, end0, extend_pileup_by)
        continue


def get_active_range(
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
    """Args:
    inactive_padding: not applied when all positions are inactive between
        "start0" and "end0" argument values.
    """

    def get_initial_fetch_range(start0, end0, extend_fetchedreads_by):
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

    def inspect_initial_range(
        chrom, start0, end0, fasta, bam, threshold, extend_fetchedreads_by
    ):
        initial_fetch_range = get_initial_fetch_range(start0, end0, extend_fetchedreads_by)
        fetchedreads_initial = libpileup.FetchedReads(
            bam, chrom, initial_fetch_range.start, initial_fetch_range.stop
        )
        pileup = libpileup.get_pileup(
            bam=bam,
            pileup_coord=(chrom, start0, end0),
            fetchedreads=fetchedreads_initial,
            fasta=fasta,
            active_threshold=threshold,
            truncate=True,
            as_array=False,
            return_range=False,
        )
        return pileup

    def find_inactive_run(active_info, inactive_padding):
        found_result = False
        farthest_pos0 = None
        for start_idx in range(0, len(active_info) - inactive_padding + 1):
            if not any(active_info[start_idx:(start_idx + inactive_padding)]):
                found_result = True
                farthest_pos0 = active_info.index[start_idx + inactive_padding -1]
                break
        return found_result, farthest_pos0

    def search_rightward(pileup, max_initial_active_pos0, inactive_padding, extend_pileup_by, extend_fetchedreads_by):
        # initial extend
        while True:
            if pileup.end0 - (max_initial_active_pos0 + 1) < inactive_padding:
                pileup.extend_rightward(width=extend_pileup_by, extend_fetchedreads_by=extend_fetchedreads_by)
            else:
                break
        # initial search loop
        active_info = pileup.get_active_info(start0=(max_initial_active_pos0 + 1))
        found_result, farthest_pos0 = find_inactive_run(active_info, inactive_padding)
        # subsequent search loop
        if not found_result:
            while True:
                pileup.extend_rightward(width=extend_pileup_by, extend_fetchedreads_by=extend_fetchedreads_by)
                active_info = pileup.get_active_info(
                    start0=(pileup.end0 - extend_pileup_by + 1 - inactive_padding)
                )
                found_result, farthest_pos0 = find_inactive_run(active_info, inactive_padding)
                if found_result:
                    break
        # prepare result
        result_stop = farthest_pos0 + 1
        return result_stop

    def search_leftward(pileup, min_initial_active_pos0, inactive_padding, extend_pileup_by, extend_fetchedreads_by):
        # initial extend
        while True:
            if (min_initial_active_pos0 - 1) - pileup.start0 + 1 < inactive_padding:
                pileup.extend_leftward(width=extend_pileup_by, extend_fetchedreads_by=extend_fetchedreads_by)
            else:
                break
        # initial search loop
        active_info = pileup.get_active_info(end0=min_initial_active_pos0)[::-1]
        found_result, farthest_pos0 = find_inactive_run(active_info, inactive_padding)
        # subsequent search loop
        if not found_result:
            while True:
                pileup.extend_leftward(width=extend_pileup_by, extend_fetchedreads_by=extend_fetchedreads_by)
                active_info = pileup.get_active_info(
                    end0=(pileup.start0 + extend_pileup_by - 1 + inactive_padding)
                )[::-1]
                found_result, farthest_pos0 = find_inactive_run(active_info, inactive_padding)
                if found_result:
                    break
        # prepare result
        result_start = farthest_pos0
        return result_start

    # main
    pileup = inspect_initial_range(
        chrom, start0, end0, fasta, bam, threshold, extend_fetchedreads_by
    )
    initial_active_poslist = pileup.get_active_positions()
    if len(initial_active_poslist) == 0:
        active_range = None
    else:
        max_initial_active_pos0 = max(initial_active_poslist)
        min_initial_active_pos0 = min(initial_active_poslist)

        result_stop = search_rightward(pileup, max_initial_active_pos0, inactive_padding, extend_pileup_by, extend_fetchedreads_by)
        result_start = search_leftward(pileup, min_initial_active_pos0, inactive_padding, extend_pileup_by, extend_fetchedreads_by)
        active_range = range(result_start, result_stop)

    return active_range, pileup

##################

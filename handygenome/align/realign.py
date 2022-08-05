import collections

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
readhandler = importlib.import_module(
    ".".join([top_package_name, "readplus", "readhandler"])
)


# get_active_range


def get_max_vaf(col):
    counts = collections.Counter(x for x in col if x is not None)
    max_vaf = max(counts.values()) / sum(counts.values())
    return max_vaf


def get_iterator(pileup, start0, end0):
    rng = range(start0, end0)
    return ((rng[idx], pileup[:, idx]) for idx in range(pileup.shape[1]))


def get_active_info(bam, chrom, start0, end0, threshold):
    pileup = readhandler.get_pileup(
        bam, chrom, start0, end0, as_array=True, truncate=True
    )
    result = list()
    for pos0, col in get_iterator(pileup, start0, end0):
        is_active = get_max_vaf(col) <= threshold
        result.append((pos0, is_active))

    return result


def check_active_col_rightward(bam, chrom, pos0_init, threshold, add_pileup_by):
    start0 = pos0_init
    end0 = start0 + add_pileup_by
    while True:
        pileup = readhandler.get_pileup(
            bam, chrom, start0, end0, as_array=True, truncate=True
        )
        for pos0, col in get_iterator(pileup, start0, end0):
            is_active = get_max_vaf(col) <= threshold
            yield pos0, is_active

        start0 = end0
        end0 = start0 + add_pileup_by
        continue


def check_active_col_leftward(bam, chrom, pos0_init, threshold, add_pileup_by):
    end0 = pos0_init + 1
    start0 = end0 - add_pileup_by
    while True:
        pileup = readhandler.get_pileup(
            bam, chrom, start0, end0, as_array=True, truncate=True
        )
        for pos0, col in reversed(tuple(get_iterator(pileup, start0, end0))):
            is_active = get_max_vaf(col) <= threshold
            yield pos0, is_active

        end0 = start0
        start0 = end0 - add_pileup_by
        continue


def get_active_range(
    bam, chrom, start0, end0, threshold=0.9, inactive_padding=10, add_pileup_by=10
):
    naive_active_info = get_active_info(bam, chrom, start0, end0, threshold)
    naive_active_positions = [
        pos0 for (pos0, is_active) in naive_active_info if is_active
    ]
    if len(naive_active_positions) == 0:
        return None

    naive_active_positions_min = min(naive_active_positions)
    naive_active_positions_max = max(naive_active_positions)

    active_info_rightward = list()
    active_info_leftward = list()

    # rightward
    current_pos0 = naive_active_positions_max + 1
    for pos0, is_active in check_active_col_rightward(
        bam,
        chrom,
        current_pos0,
        threshold=threshold,
        add_pileup_by=add_pileup_by,
    ):
        active_info_rightward.append((pos0, is_active))
        if len(active_info_rightward) >= inactive_padding:
            if not any(x[1] for x in active_info_rightward[-inactive_padding:]):
                break

    # leftward
    current_pos0 = naive_active_positions_min - 1
    for pos0, is_active in check_active_col_leftward(
        bam,
        chrom,
        current_pos0,
        threshold=threshold,
        add_pileup_by=add_pileup_by,
    ):
        active_info_leftward.append((pos0, is_active))
        if len(active_info_leftward) >= inactive_padding:
            if not any(x[1] for x in active_info_leftward[-inactive_padding:]):
                break

    active_range = range(active_info_leftward[-1][0], active_info_rightward[-1][0] + 1)

    return active_range


##################

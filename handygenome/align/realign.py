import collections

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
libpileup = importlib.import_module(".".join([top_package_name, "readplus", "pileup"]))


# get_active_range
def active_info_generator(bam, chrom, start0, end0, threshold, reverse=False):
    pileup = libpileup.get_pileup(
        bam, chrom, start0, end0, as_array=False, truncate=True
    )
    if reverse:
        iterator = pileup.df.iloc[:, ::-1].iteritems()
    else:
        iterator = pileup.df.iteritems()

    for pos0, col in iterator:
        max_vaf = libpileup.get_max_vaf_from_col(col)
        if max_vaf is None:
            is_active = False
        else:
            is_active = max_vaf <= threshold

        yield (pos0, is_active)


def get_active_info(bam, chrom, start0, end0, threshold):
    return list(active_info_generator(bam, chrom, start0, end0, threshold))


def check_active_col_rightward(bam, chrom, pos0_init, threshold, add_pileup_by):
    start0 = pos0_init
    end0 = start0 + add_pileup_by
    while True:
        for pos0, is_active in active_info_generator(
            bam, chrom, start0, end0, threshold, reverse=False,
        ):
            yield pos0, is_active

        start0 = end0
        end0 = start0 + add_pileup_by
        continue


def check_active_col_leftward(bam, chrom, pos0_init, threshold, add_pileup_by):
    end0 = pos0_init + 1
    start0 = end0 - add_pileup_by
    while True:
        for pos0, is_active in active_info_generator(
            bam, chrom, start0, end0, threshold, reverse=True,
        ):
            yield pos0, is_active

        end0 = start0
        start0 = end0 - add_pileup_by
        continue


def get_active_range(
    bam, chrom, start0, end0, threshold=0.9, inactive_padding=10, add_pileup_by=10
):
    naive_active_info = list(active_info_generator(bam, chrom, start0, end0, threshold, reverse=False))
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

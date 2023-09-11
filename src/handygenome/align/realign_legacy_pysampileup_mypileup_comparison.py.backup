import collections

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
readhandler = importlib.import_module('.'.join([top_package_name, 'readplus', 'readhandler']))


def call_pileup(bam, chrom, start0, end0):
    pup = bam.pileup(chrom, start0, end0, 
                     truncate=True, max_depth=10**6, stepper='all', 
                     ignore_overlaps=False, ignore_orphans=False, 
                     min_base_quality=0, min_mapping_quality=0, 
                     adjust_capq_threshold=0, compute_baq=False, 
                     redo_baq=False)    

    return pup

    
def get_max_vaf_pysampileup(pupcol):
    counts = collections.Counter(pupcol)
    max_vaf = max(counts.values()) / sum(counts.values())
    return max_vaf


def get_is_active_pupcol(pupcol, threshold):
    max_vaf = get_max_vaf_pysampileup(
        x.upper() for x in
        pupcol.get_query_sequences(mark_matches=False, mark_ends=False, add_indels=True)
    )
    is_active = max_vaf <= threshold
    return pupcol.reference_pos, is_active


def get_active_info_from_pileup(pup, threshold):
    return [get_is_active_pupcol(pupcol, threshold) for pupcol in pup]


def check_active_col_pysampileup_rightward(bam, chrom, pos0_init, threshold, add_pileup_by):
    start0 = pos0_init
    end0 = start0 + add_pileup_by
    while True:
        pup = call_pileup(bam, chrom, start0, end0)
        for pupcol in pup:
            pos0, is_active = get_is_active_pupcol(pupcol, threshold)
            yield pos0, is_active

        start0 = end0
        end0 = start0 + add_pileup_by
        continue

        
def check_active_col_pysampileup_leftward(bam, chrom, pos0_init, threshold, add_pileup_by):
    end0 = pos0_init + 1
    start0 = end0 - add_pileup_by
    while True:
        pup = call_pileup(bam, chrom, start0, end0)
        active_info_buffer = get_active_info_from_pileup(pup, threshold)
        for pos0, is_active in reversed(active_info_buffer):
            yield pos0, is_active

        end0 = start0
        start0 = end0 - add_pileup_by
        continue

            
def get_active_range_pysampileup(bam, chrom, start0, end0, threshold=0.9, inactive_padding=5, add_pileup_by=10):
#     readlist = list(bam.fetch(chrom, start0, end0))
#     pileup_df = get_pileup(readlist, as_df=True)
    pup = call_pileup(bam, chrom, start0, end0)
    naive_active_info = get_active_info_from_pileup(pup, threshold)        
    naive_active_positions = [pos0 for (pos0, is_active) in naive_active_info if is_active]    
    print(naive_active_info)
    if len(naive_active_positions) == 0:
        return
    
    naive_active_positions_min = min(naive_active_positions)
    naive_active_positions_max = max(naive_active_positions)
    
#     idx_min = naive_active_info.index((naive_active_positions_min, True))
#     idx_max = naive_active_info.index((naive_active_positions_max, True))
#     naive_active_info_filtered = naive_active_info[idx_min:(idx_min + 1)]
    
    active_info_rightward = list()
    active_info_leftward = list()
    
    # rightward
    current_pos0 = naive_active_positions_max + 1
    for pos0, is_active in check_active_col_pysampileup_rightward(
        bam, chrom, current_pos0, threshold=threshold, add_pileup_by=add_pileup_by
    ):
        active_info_rightward.append((pos0, is_active))
        if len(active_info_rightward) >= inactive_padding:
            if not any(x[1] for x in active_info_rightward[-inactive_padding:]):
                break

    # leftward
    current_pos0 = naive_active_positions_min - 1
    for pos0, is_active in check_active_col_pysampileup_leftward(
        bam, chrom, current_pos0, threshold=threshold, add_pileup_by=add_pileup_by
    ):
        active_info_leftward.append((pos0, is_active))
        if len(active_info_leftward) >= inactive_padding:
            if not any(x[1] for x in active_info_leftward[-inactive_padding:]):
                break
                
#     active_info = active_info_leftward[::-1] + naive_active_info_filtered + active_info_rightward
#     positions = [pos0 for (pos0, is_active) in active_info]    
    active_range = range(active_info_leftward[-1][0], active_info_rightward[-1][0] + 1)
    
    return active_range


###############################

def get_max_vaf_mypileup(col):
    counts = collections.Counter(x for x in col if x is not None)
    max_vaf = max(counts.values()) / sum(counts.values())    
    return max_vaf

def get_iterator(pileup, as_array, start0, end0):
    if as_array:
        rng = range(start0, end0)
        return ((rng[idx], pileup[:, idx]) for idx in range(pileup.shape[1]))
    else:
        return pileup.items()

    
def get_active_info_mypileup(bam, chrom, start0, end0, threshold, as_array):
    pileup = readhandler.get_pileup(bam, chrom, start0, end0, as_array=as_array, truncate=True)
    result = list()

    for pos0, col in get_iterator(pileup, as_array, start0, end0):
        is_active = (get_max_vaf_mypileup(col) <= threshold)
        result.append((pos0, is_active))
    return result


def check_active_col_mypileup_rightward(bam, chrom, pos0_init, threshold, add_pileup_by, as_array):
    start0 = pos0_init
    end0 = start0 + add_pileup_by
    while True:
        pileup = readhandler.get_pileup(bam, chrom, start0, end0, as_array=as_array, truncate=True)
        for pos0, col in get_iterator(pileup, as_array, start0, end0):
            is_active = (get_max_vaf_mypileup(col) <= threshold)
            yield pos0, is_active

        start0 = end0
        end0 = start0 + add_pileup_by
        continue
        
        
def check_active_col_mypileup_leftward(bam, chrom, pos0_init, threshold, add_pileup_by, as_array):
    end0 = pos0_init + 1
    start0 = end0 - add_pileup_by
    while True:
        pileup = readhandler.get_pileup(bam, chrom, start0, end0, as_array=as_array, truncate=True)
        for pos0, col in reversed(tuple(get_iterator(pileup, as_array, start0, end0))):
            is_active = (get_max_vaf_mypileup(col) <= threshold)
            yield pos0, is_active

        end0 = start0
        start0 = end0 - add_pileup_by
        continue    
    
###
    

def get_active_range_mypileup(bam, chrom, start0, end0, threshold=0.9, inactive_padding=5, add_pileup_by=10, as_array=False):
    naive_active_info = get_active_info_mypileup(bam, chrom, start0, end0, threshold, as_array)
    naive_active_positions = [pos0 for (pos0, is_active) in naive_active_info if is_active]    
    print(naive_active_info)
    if len(naive_active_positions) == 0:
        return
    
    naive_active_positions_min = min(naive_active_positions)
    naive_active_positions_max = max(naive_active_positions)
        
    active_info_rightward = list()
    active_info_leftward = list()
    
    # rightward
    current_pos0 = naive_active_positions_max + 1
#     n = 0
    for pos0, is_active in check_active_col_mypileup_rightward(
        bam, chrom, current_pos0, threshold=threshold, add_pileup_by=add_pileup_by, as_array=as_array,
    ):
#         n += 1
#         if n == 500:
#             print('overflow')
#             break
        active_info_rightward.append((pos0, is_active))
        if len(active_info_rightward) >= inactive_padding:
            if not any(x[1] for x in active_info_rightward[-inactive_padding:]):
                break

    # leftward
    current_pos0 = naive_active_positions_min - 1
#     n = 0
    for pos0, is_active in check_active_col_mypileup_leftward(
        bam, chrom, current_pos0, threshold=threshold, add_pileup_by=add_pileup_by, as_array=as_array,
    ):
#         n += 1
#         if n == 500:
#             print('overflow')
#             break
        active_info_leftward.append((pos0, is_active))
        if len(active_info_leftward) >= inactive_padding:
            if not any(x[1] for x in active_info_leftward[-inactive_padding:]):
                break
                
#     active_info = active_info_leftward[::-1] + naive_active_info_filtered + active_info_rightward
#     positions = [pos0 for (pos0, is_active) in active_info]    
    active_range = range(active_info_leftward[-1][0], active_info_rightward[-1][0] + 1)
    
    return active_range
import collections
import operator

import Bio.Align

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
readplus = importlib.import_module('.'.join([top_package_name, 'read', 'readplus']))


DEFAULT_FLANKLEN = 1

ALIGNER_COMPARE = Bio.Align.PairwiseAligner(
    mode='global',
    match_score=0,
    mismatch_score=-1,
    query_gap_score=-1,
    target_gap_score=-1,
)


"""Meaning of alleleclass:
    None: not informative. When a read does not span flanking areas.
        All other alleleclasses imply a read spans flanking areas.
    -1: other than REF and ALTS
    0: REF
    1: 1st ALT
    2: 2nd ALT
    ...
"""


# merging results of two readplus objects

def merge_monoalt_supports(supports_rp1, supports_rp2):
    """Args:
        results of "get_normalized_monoalt_supports_readplus" function
    """
    merged_monoalt_supports = dict()
    for alt_index in supports_rp1.keys():
        merged_monoalt_supports[alt_index] = merge_alleleclasses(
            supports_rp1[alt_index], supports_rp2[alt_index],
        )

    return merged_monoalt_supports


def merge_alleleclasses(alleleclass_rp1, alleleclass_rp2):
    if (alleleclass_rp1 is not None) and (alleleclass_rp2 is not None):
        if alleleclass_rp1 == alleleclass_rp2:
            return alleleclass_rp1
        else:
            return None
    elif (alleleclass_rp1 is None) and (alleleclass_rp2 is None):
        return None
    else:
        if alleleclass_rp1 is not None:
            return alleleclass_rp1
        elif alleleclass_rp2 is not None:
            return alleleclass_rp2


#############################################

# asis

def get_alleleclass_asis_readplus_new(
    vcfspec, rp, 
    coverage_cutoff=1, 
    partial_coverage_del_len_cutoff=5,
    compare_mode='identity', 
    similarity_params={'factor': 0.1},
):
    """Each ReadPlus can support no more than one allele"""

    #assert compare_mode in ('identity', 'similarity')

    # cigar N full span check
    if rp.check_cigarN_includes_range(vcfspec.range0):
        return None

    # treat as None if not spanning relevant repeat region
    read_range0 = rp.get_range0()
    spans_repeat_list = list()
    for repeat_spec in vcfspec.repeat_specs:
        if repeat_spec is None:
            spans_repeat = True
        else:
            spans_repeat = (
                (repeat_spec.start0 in read_range0)
                and (repeat_spec.end0 in read_range0)  
            )  # The position 1 base right to the repeat end must be spanned by the read
        spans_repeat_list.append(spans_repeat)

    if not all(spans_repeat_list):
        return None

    # vcfspec coverage
    spans_left, spans_right, relevant_read_seq = get_vcfspec_coverage_info(rp, vcfspec)

    # treat as None by vcfspec coverage
    if (not spans_left) and (not spans_right):  # no span
        return None
    if (not spans_left) or (not spans_right):  # partial span
        vcfspec_coverage, modified_vcfspec_range0 = get_vcfspec_coverage_fraction(vcfspec, rp)
        if len(modified_vcfspec_range0) >= partial_coverage_del_len_cutoff:
            # len(modified_vcfspec_range0) equals to the deletion length
            if not (vcfspec_coverage >= coverage_cutoff):
                return None
        else:
            return None
        
    # get seqs to use as target (read sequence is query)
    target_seqs = get_target_seqs(vcfspec, spans_left, spans_right, relevant_read_seq)
    # compare target and query seqs
    alleleclass = compare_targets_and_query(
        target_seqs, relevant_read_seq, mode=compare_mode, similarity_params=similarity_params,
    )

    return alleleclass


def get_alleleclass_asis_readplus(vcfspec, rp, flanklen=DEFAULT_FLANKLEN):
    """Each ReadPlus can support no more than one allele"""

    assert flanklen >= 1, f'"flanklen" argument must be at least 1.'

    # cigar N full span check
    if rp.check_cigarN_includes_range(vcfspec.range0):
        return None

    spans, matches = rp.check_spans_and_matches_vcfspec_flanks(vcfspec, flanklen=flanklen)
    if spans:
        if matches:
            allele_seq = rp.get_seq_from_range0(
                vcfspec.REF_range0,
                flanking_queryonly_default_mode=False,
                include_leading_queryonly=True, 
                include_trailing_queryonly=True,
            )

            if allele_seq in vcfspec.alleles:
                alleleclass = vcfspec.alleles.index(allele_seq)
            else:
                alleleclass = -1
        else:
            alleleclass = -1
    else:
        if 'S' in rp.read.cigarstring:
            if rp.check_softclip_spans_vcfspec_flanks(vcfspec, flanklen=flanklen):
                alleleclass = -1
            else:
                alleleclass = None
        else:
            alleleclass = None

    return alleleclass


################################################################

# helpers

def alleleclass_sortkey(x):
    if x is None:
        return -2
    else:
        return x


def get_vcfspec_coverage_info(rp, vcfspec):
    read_range0 = rp.range0
    spans_left = (min(vcfspec.readspan_range0) in read_range0)
    spans_right = (max(vcfspec.readspan_range0) in read_range0)
    relevant_read_seq = rp.get_seq_from_pairs_indexes(
        rp.get_pairs_indexes(
            vcfspec.readspan_range0, flanking_queryonly_default_mode=True,
        )
    )

    return spans_left, spans_right, relevant_read_seq


def get_vcfspec_coverage_fraction(vcfspec, rp):
    mutation_types = set(
        vcfspec.get_mutation_type(alt_index=alt_index)
        for alt_index in vcfspec.iter_alt_indexes()
    )
    if mutation_types == {'del'}:
        modified_vcfspec_range0 = range(
            vcfspec.readspan_range0.start + 1,
            vcfspec.readspan_range0.stop - 1,
        )
    else:
        modified_vcfspec_range0 = vcfspec.readspan_range0

    read_span_start = max(rp.range0.start, modified_vcfspec_range0.start)
    read_span_stop = min(rp.range0.stop, modified_vcfspec_range0.stop)
    numerator = read_span_stop - read_span_start
    denominator = len(modified_vcfspec_range0)
    vcfspec_coverage = numerator / denominator

    return vcfspec_coverage, modified_vcfspec_range0


def get_target_seqs(vcfspec, spans_left, spans_right, relevant_read_seq):
    target_seqs = vcfspec.compared_seq_dict.copy()
    if spans_left and (not spans_right):
        target_seqs = {key: val[:len(relevant_read_seq)] for key, val in target_seqs.items()}
    elif (not spans_left) and spans_right:
        target_seqs = {key: val[-len(relevant_read_seq):] for key, val in target_seqs.items()}

    return target_seqs


def compare_targets_and_query(target_seqs, relevant_read_seq, mode, similarity_params):
    """Args:
        target_seqs: Result of Vcfspec.compared_seq_dict method. 
            keys are allele_index.
    Result:
        alleleclass
    """
    assert mode in ('identity', 'similarity')

    if mode == 'identity':
        hit_allele_indexes = list()
        for allele_index, target in target_seqs.items():
            if target == relevant_read_seq:
                hit_allele_indexes.append(allele_index)

        if len(hit_allele_indexes) == 1:
            return hit_allele_indexes[0]
        else:
            return None

    elif mode == 'similarity':
        alignment_info = list()
        for allele_index, target in target_seqs.items():
            alns = ALIGNER_COMPARE.align(target, relevant_read_seq)
            if similarity_checker(alns, **similarity_params):
                alignment_info.append((allele_index, alns.score))
        if len(alignment_info) == 0:
            return None
        else:
            alignment_info.sort(key=operator.itemgetter(1))
            return alignment_info[-1][0]


def similarity_checker(alns, factor=0.1):
    """Returns:
        True means passed
    """
    query_len = len(alns.sequences[1])
    return -alns.score < (query_len * factor)
            

################################################################

# monoalt support

class VariantSupportMonoALT(
    collections.namedtuple(
        'VariantSupport', 
        ('supports', 'relevant')
    )
):
    pass


def check_monoalt_support(vcfspec, rp, flanklen=DEFAULT_FLANKLEN):
    alleleclass = get_alleleclass_asis_readplus(vcfspec, rp, flanklen=flanklen)
    if alleleclass is None:
        return None
    else:
        return (alleleclass == 1)


def get_normalized_monoalt_supports_readplus(vcfspec, rp, flanklen=DEFAULT_FLANKLEN):
    monoalt_supports = dict()
    monoalt_list = [x.normalize() for x in vcfspec.merge_components]
    #for alt_index, vcfspec_monoalt in enumerate(monoalt_list):
    for alt_index, vcfspec_monoalt in enumerate(monoalt_list):
        monoalt_supports[alt_index] = check_monoalt_support(vcfspec_monoalt, rp, flanklen=flanklen)
        #monoalt_supports.append(check_monoalt_support(vcfspec_monoalt, rp, flanklen=flanklen))

    return monoalt_supports



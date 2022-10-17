import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
readplus = importlib.import_module('.'.join([top_package_name, 'read', 'readplus']))


DEFAULT_FLANKLEN = 1


def make_alleleinfoitem_readpluspair(vcfspec, aiitem_rp1, aiitem_rp2):
    alleleinfoitem = dict()

    if aiitem_rp1['spans'] and aiitem_rp2['spans']:
        if aiitem_rp1['alleleclass'] == aiitem_rp2['alleleclass']:
            alleleinfoitem['alleleclass'] = aiitem_rp1['alleleclass']
        else:
            alleleinfoitem['alleleclass'] = None
    else:
        if (not aiitem_rp1['spans']) and (not aiitem_rp2['spans']):
            alleleinfoitem['alleleclass'] = None
        else:
            if aiitem_rp1['spans']:
                alleleinfoitem['alleleclass'] = aiitem_rp1['alleleclass']
            elif aiitem_rp2['spans']:
                alleleinfoitem['alleleclass'] = aiitem_rp2['alleleclass']

    return alleleinfoitem


def make_alleleinfoitem_readplus(vcfspec, rp, flanklen=DEFAULT_FLANKLEN):
    assert flanklen >= 1, f'"flanklen" argument must be at least 1.'

    preflank_range0, postflank_range0 = get_flank_ranges(vcfspec, flanklen, rp.fasta)
    spans = (rp.check_spans(preflank_range0) and rp.check_spans(postflank_range0))
    if spans:
        alleleclass = get_alleleclass(vcfspec, rp, preflank_range0, postflank_range0)
    else:
        alleleclass = None

    alleleinfoitem = dict()
    alleleinfoitem['spans'] = spans
    alleleinfoitem['alleleclass'] = alleleclass

    return alleleinfoitem


def get_flank_ranges(vcfspec, flanklen, fasta):
    equivalents = vcfspec.get_equivalents(fasta)
    preflank_range0 = get_preflank_range0(equivalents[0], flanklen)
    postflank_range0 = equivalents[-1].get_postflank_range0(flanklen=flanklen)

    return preflank_range0, postflank_range0


def get_preflank_range0(vcfspec, flanklen):
#    preflank_range0_candidates = set(
#        vcfspec.get_preflank_range0(idx=idx, flanklen=flanklen)
#        for idx in range(len(vcfspec.alts))
#    )
#    if len(preflank_range0_candidates) != 1:
#        raise Exception(
#            f'There are more than one possible pre-flanking ranges '
#            f'for the input vcfspec: {vcfspec}')
#    else:
#        preflank_range0 = preflank_range0_candidates.pop()

    return min(
        (
            vcfspec.get_preflank_range0(idx=idx, flanklen=flanklen)
            for idx in range(len(vcfspec.alts))
        ),
        key=(lambda x: x.start),
    )



def get_alleleclass(vcfspec, rp, preflank_range0, postflank_range0):
    """Returns:
        None: not informative
        -1: other than REF and ALTS
        0: REF
        1: 1st ALT
        2: 2nd ALT
        ...
    """

    if rp.check_matches(preflank_range0) and rp.check_matches(postflank_range0):
        allele_seq = rp.get_seq_from_range0(vcfspec.REF_range0)
        alleles = vcfspec.alleles
        if allele_seq in alleles:
            alleleclass = alleles.index(allele_seq)
        else:
            alleleclass = -1
    else:
        alleleclass = -1

    return alleleclass

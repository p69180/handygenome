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

    preflank_range0, postflank_range0 = vcfspec.get_flank_range0s_equivalents(flanklen=flanklen)
    spans = (rp.check_spans(preflank_range0) and rp.check_spans(postflank_range0))
    if spans:
        alleleclass = get_alleleclass(vcfspec, rp, preflank_range0, postflank_range0)
    else:
        alleleclass = None

    alleleinfoitem = dict()
    alleleinfoitem['spans'] = spans
    alleleinfoitem['alleleclass'] = alleleclass

    return alleleinfoitem


def get_alleleclass(vcfspec, rp, preflank_range0, postflank_range0):
    """Returns:
        None: not informative
        -1: other than REF and ALTS
        0: REF
        1: 1st ALT
        2: 2nd ALT
        ...
    """

    if (
        rp.check_matches(
            preflank_range0,
            flanking_queryonly_default_mode=False,
            include_leading_queryonly=True, 
            include_trailing_queryonly=False,
        ) and 
        rp.check_matches(
            postflank_range0,
            flanking_queryonly_default_mode=False,
            include_leading_queryonly=False, 
            include_trailing_queryonly=True,
        )
    ):
        allele_seq = rp.get_seq_from_range0(
            vcfspec.REF_range0,
            flanking_queryonly_default_mode=False,
            include_leading_queryonly=True, 
            include_trailing_queryonly=True,
        )

        alleles = vcfspec.alleles
        if allele_seq in alleles:
            alleleclass = alleles.index(allele_seq)
        else:
            alleleclass = -1
    else:
        alleleclass = -1

    return alleleclass
"""        border-          
            line
    NONCLIP   |  CLIP
              |
------------- | ------\
| | | | | |*| | | | | |
------------- | ------/
   BND SIDE   |  PARTNER SIDE (par side)
dist.    prox.|

* : endpoint

=========================================

              |           |
--------------|---|---|---|--\
              |   |   |   |  |
--------------|---|---|---|--/
              | borderzone| 

=========================================

- BND side (distal vs. proximal)
- PARTNER side
- endpoint: the aligned base most proximal to borderline
- borderzone: region composed with endpoints of all equivalent forms.
    represented by varibable "pos_range0"
- borderline

- advance vs retract
    - advance: borderline moves such that length of aligned bases increases
    - retract: opposite
"""

DEFAULT_FLANKLEN_PARSIDE = 1
DEFAULT_FLANKLEN_BNDSIDE = 1
BND_DISTANCE_THRESHOLD = 1000


##################
# basic concepts #
##################

def on_bndside_notonborder(read, pos_range0, is5prime):
    """Entire length of read lies on BND SIDE (based on most RETRACTED position!!)
        Does not span borderzone
    """
    if is5prime:
        return read.reference_start > pos_range0.start
    else:
        return read.reference_end <= pos_range0.start


def on_bndside(read, pos_range0, is5prime):
    """Entire length of read lies on BND SIDE (based on most ADVANCED position!!)
    """
    if is5prime:
        return read.reference_start > pos_range0.stop
    else:
        return read.reference_end <= pos_range0.stop


# parside: partner side (named in conformity with the VCF spec)
def on_parside(read, pos_range0, is5prime):
    """Entire length of read lies on PARTNER SIDE"""
    if is5prime:
        return read.reference_end <= pos_range0.stop + 1
    else:
        return read.reference_start >= pos_range0.stop


def clipped_toward_border(read, is5prime):
    if is5prime:
        return read.cigartuples[0][0] in (4, 5)
    else:
        return read.cigartuples[-1][0] in (4, 5)


def faces_border_frombnd(read, is5prime):
    if is5prime:
        return read.is_reverse
    else:
        return read.is_forward
   
    
def faces_border_frompar(read, is5prime):
    if is5prime:
        return read.is_forward
    else:
        return read.is_reverse


################
# for readplus #
################

def handle_bndside_notonborder(rp, is5prime, alleleclass_item):
    if faces_border_frombnd(rp.read, is5prime):
        if clipped_toward_border(rp.read, is5prime):
            alleleclass_item['noninformative'] = True
        else:
            alleleclass_item['facing_frombnd'] = True
    else:
        alleleclass_item['noninformative'] = True


def handle_bndside_onborder_clipped(
    bnds, rp, pos_range0, is5prime, alleleclass_item, is_bnd1,
    flanklen_parside, 
    bndside_distal_flank,
):
    """Only for reads which span 'bndside_distal_flank'"""
    required_queryseq_length = len(pos_range0) + flanklen_parside
    queryseq_borderzone_andbeyond = rp.get_seq_from_coord(
        start0=pos_range0.start, 
        length=required_queryseq_length,
        forward=(not is5prime),
    )

    if len(queryseq_borderzone_andbeyond) < required_queryseq_length:
        alleleclass_item['noninformative'] = True
    else:  # now flanking length conditions are fulfilled
        if rp.check_matches(bndside_distal_flank):
            seq_ref = bnds.get_seq_beyond_bnd_ref(
                is_bnd1=is_bnd1, 
                length=flanklen_parside, 
                with_border_seq=True,
            )
            seq_alt = bnds.get_seq_beyond_bnd_alt(
                is_bnd1=is_bnd1, 
                length=flanklen_parside,
                with_border_seq=True,
            )
            if queryseq_borderzone_andbeyond == seq_ref:
                alleleclass_item['ref_support'] = True
            elif queryseq_borderzone_andbeyond == seq_alt:
                alleleclass_item['alt_support'] = True
            else:
                alleleclass_item['other_support'] = True
        else:
            alleleclass_item['other_support'] = True


def handle_crosser(
    bnds, rp, pos_range0, is5prime, flanklen_parside,
    alleleclass_item, bndside_distal_flank,
    is_bnd1,
):
    queryseq_beyond_bnd = rp.get_seq_from_coord(
        start0=pos_range0.stop, 
        length=flanklen_parside,
        forward=(not is5prime),
    )
    if len(queryseq_beyond_bnd) < flanklen_parside:
        alleleclass_item['noninformative'] = True
    else:  # now flanking length conditions are fulfilled
        if (
            rp.check_matches(bndside_distal_flank) 
            and rp.check_matches(pos_range0)
        ):
            seq_beyond_bnd_ref = bnds.get_seq_beyond_bnd_ref(
                is_bnd1=is_bnd1, 
                length=flanklen_parside, 
                with_border_seq=False,
            )
            seq_beyond_bnd_alt = bnds.get_seq_beyond_bnd_alt(
                is_bnd1=is_bnd1, 
                length=flanklen_parside,
                with_border_seq=False,
            )

            if queryseq_beyond_bnd == seq_beyond_bnd_ref:
                alleleclass_item['ref_support'] = True
            elif queryseq_beyond_bnd == seq_beyond_bnd_alt:
                alleleclass_item['alt_support'] = True
            else:
                alleleclass_item['other_support'] = True
        else:
            alleleclass_item['other_support'] = True


def make_alleleclass_item_readplus(
    bnds, is_bnd1, rp,
    flanklen_parside=DEFAULT_FLANKLEN_PARSIDE, 
    flanklen_bndside=DEFAULT_FLANKLEN_BNDSIDE,
):
    # get parameters
    (chrom, pos_range0, is5prime) = bnds.get_params(is_bnd1)
    bndside_distal_flank = bnds.get_flank_range0(
        mode='bnd_dist', 
        is_bnd1=is_bnd1, 
        flanklen=flanklen_bndside,
    )

    # initialize alleleclass_item
    alleleclass_item = {
        'noninformative': False,
        'ref_support': False,  
        'alt_support': False,
        'other_support': False,
        'facing_frombnd': False,
        'facing_frompar': False,
    }
        # ref_support, alt_support: 
            # crossing border
            # matching sequence beyond borderline
            # no need to match sequence before borderline
        # other_support, alt_support: crossing border and sequence mismatch
        # noninformative:
            # too far away from borderline (BND_DISTANCE_THRESHOLD)
            # not facing toward borderline
            # peculiar (too short?) ones (on BND side & overlapping borderzone & not spanning bndside flank)

    # main
    if chrom == rp.read.reference_name:
        distance = rp.get_distance(pos_range0)
        if distance > BND_DISTANCE_THRESHOLD:
            alleleclass_item['noninformative'] = True  # too far away
        else:
            if on_bndside(rp.read, pos_range0, is5prime):
                if on_bndside_notonborder(rp.read, pos_range0, is5prime):
                    handle_bndside_notonborder(rp, is5prime, alleleclass_item)
                else:  # on bndside and on the border
                    if rp.check_spans(bndside_distal_flank):  # flank region most far from borderline
                        if clipped_toward_border(rp.read, is5prime):
                            handle_bndside_onborder_clipped(  # POTENTIAL FOR ALT SUPPORT
                                bnds, rp, pos_range0, is5prime, alleleclass_item, 
                                is_bnd1, flanklen_parside, 
                                bndside_distal_flank,
                            )
                        else:
                            if faces_border_frombnd(rp.read, is5prime):
                                alleleclass_item['facing_frombnd'] = True
                            else:
                                alleleclass_item['noninformative'] = True  # not facing borderline
                    else:
                        alleleclass_item['noninformative'] = True  # PECULIAR

            elif on_parside(rp.read, pos_range0, is5prime):
                if faces_border_frompar(rp.read, is5prime):
                    alleleclass_item['facing_frompar'] = True
                else:
                    alleleclass_item['noninformative'] = True  # not facing borderline

            else:  # crosses the most advanced border
                if rp.check_spans(bndside_distal_flank):
                    handle_crosser(
                        bnds, rp, pos_range0, is5prime, flanklen_parside,
                        alleleclass_item, bndside_distal_flank,
                        is_bnd1,
                    )
                else:
                    alleleclass_item['noninformative'] = True  # PECULIAR
    else:
        alleleclass_item['noninformative'] = True

    return alleleclass_item


####################
# for readpluspair #
####################

def init_rpp_alleleclass_item():
    return {
        'ref_support_direct': False,
        'ref_support_indirect': False,
        'alt_support_direct': False,
        'alt_support_indirect': False,
        'other_support': False,
        'noninformative': False,
    }
        # alt_support_direct: borderline crossing and seqeunce matching
        # alt_support_indirect: not crossing, but facing each other, on different BND sides


def get_alleleclasses():
    return tuple(init_rpp_alleleclass_item().keys())


def check_rpp_other_support(aiitem_rp1, aiitem_rp2, aiitem_rpp):
    if any(
        (
            aiitem_rp1['other_support'], 
            aiitem_rp2['other_support'], 
        )
    ):
        aiitem_rpp['other_support'] = True


def check_rpp_direct_alt_support(
    aiitem_bnd1_rp1, aiitem_bnd1_rp2, alleleclass_item_bnd1,
    aiitem_bnd2_rp1, aiitem_bnd2_rp2, alleleclass_item_bnd2,
):
    if (
        (
            aiitem_bnd1_rp1['alt_support']
            and (
                aiitem_bnd1_rp2['alt_support']
                or aiitem_bnd2_rp2['alt_support']
                or aiitem_bnd2_rp2['facing_frombnd']
            )
        )
        or (
            aiitem_bnd1_rp2['alt_support']
            and (
                aiitem_bnd1_rp1['alt_support']
                or aiitem_bnd2_rp1['alt_support']
                or aiitem_bnd2_rp1['facing_frombnd']
            )
        )
    ):
        alleleclass_item_bnd1['alt_support_direct'] = True

    if (
        (
            aiitem_bnd2_rp1['alt_support']
            and (
                aiitem_bnd2_rp2['alt_support']
                or aiitem_bnd1_rp2['alt_support']
                or aiitem_bnd1_rp2['facing_frombnd']
            )
        )
        or (
            aiitem_bnd2_rp2['alt_support']
            and (
                aiitem_bnd2_rp1['alt_support']
                or aiitem_bnd1_rp1['alt_support']
                or aiitem_bnd1_rp1['facing_frombnd']
            )
        )
    ):
        alleleclass_item_bnd2['alt_support_direct'] = True


def check_rpp_ref_support(aiitem_rp1, aiitem_rp2, aiitem_rpp):
    if not (
        aiitem_rpp['other_support'] 
        or aiitem_rpp['alt_support_direct'] 
        or aiitem_rpp['alt_support_indirect']
    ):  # regarded not as ref-supporting when other- or alt- supporting
        if (
            (
                aiitem_rp1['facing_frombnd']
                and aiitem_rp2['facing_frompar']
            )
            or (
                aiitem_rp2['facing_frombnd'] 
                and aiitem_rp1['facing_frompar']
            )
        ):
            # the mate reads facing each other, one read on bndside of a 
                # breakend and the other read on parside of the same breakend
            aiitem_rpp['ref_support_indirect'] = True
        else:
            if (
                (
                    aiitem_rp1['ref_support']
                    and (
                        aiitem_rp2['ref_support']
                        or aiitem_rp2['facing_frompar']
                        or aiitem_rp2['facing_frombnd']
                    )
                )
                or (
                    aiitem_rp2['ref_support']
                    and (
                        aiitem_rp1['ref_support']
                        or aiitem_rp1['facing_frompar']
                        or aiitem_rp1['facing_frombnd']
                    )
                )
            ):
            # direct ref support with read crossing
                aiitem_rpp['ref_support_direct'] = True


def check_rpp_noninformative(aiitem_rpp):
    if all(
        (not val) 
        for (key, val) in aiitem_rpp.items()
        if (key != 'noninformative')
    ):
        aiitem_rpp['noninformative'] = True


def make_alleleclass_item_readpluspair_old(
    aiitem_bnd1_rp1, aiitem_bnd1_rp2,
    aiitem_bnd2_rp1, aiitem_bnd2_rp2, 
    bnds,
):
    """Assumes that each rp supports at most 1 bnd (one of bnd1 or bnd2).
        This assumption is monitored in "update_alleleclass_sv" method
        of ReadPlus class.
    """

    alleleclass_item = init_rpp_alleleclass_item()

    # other
    if any(
        (
            aiitem_bnd1_rp1['other_support'], 
            aiitem_bnd2_rp1['other_support'],
            aiitem_bnd1_rp2['other_support'], 
            aiitem_bnd2_rp2['other_support'],
        )
    ):
        alleleclass_item['other_support'] = True

    # alt
    if any(
        (
            aiitem_bnd1_rp1['alt_support'], 
            aiitem_bnd2_rp1['alt_support'],
            aiitem_bnd1_rp2['alt_support'],
            aiitem_bnd2_rp2['alt_support'],
        )
    ):
        alleleclass_item['alt_support_direct'] = True
    elif (
        (
            aiitem_bnd1_rp1['facing_frombnd'] 
            and aiitem_bnd2_rp2['facing_frombnd']
        ) 
        or (
            aiitem_bnd1_rp2['facing_frombnd']
            and aiitem_bnd2_rp1['facing_frombnd']
        )
    ):
        alleleclass_item['alt_support_indirect'] = True

    # ref
    if not (
        alleleclass_item['other_support'] 
        or alleleclass_item['alt_support_direct'] 
        or alleleclass_item['alt_support_indirect']
    ):  # regarded not as ref-supporting when other- or alt- supporting
        if any(
            (
                aiitem_bnd1_rp1['ref_support'], 
                aiitem_bnd2_rp1['ref_support'],
                aiitem_bnd1_rp2['ref_support'],
                aiitem_bnd2_rp2['ref_support'],
            )
        ):  # direct ref support with read crossing
            alleleclass_item['ref_support_direct'] = True
        elif any(
            (
                (
                    aiitem_bnd1_rp1['facing_frombnd']
                    and aiitem_bnd1_rp2['facing_frompar']
                ),
                (
                    aiitem_bnd1_rp2['facing_frombnd'] 
                    and aiitem_bnd1_rp1['facing_frompar']
                ),
                (
                    aiitem_bnd2_rp1['facing_frombnd'] 
                    and aiitem_bnd2_rp2['facing_frompar']
                ),
                (
                    aiitem_bnd2_rp2['facing_frombnd'] 
                    and aiitem_bnd2_rp1['facing_frompar']
                ),
            )
        ):
            # the mate reads facing each other, one read on bndside of a 
                # breakend and the other read on parside of the same breakend
            alleleclass_item['ref_support_indirect'] = True

    # otherwise noninformative
    if not any(
        v for (k, v) in alleleclass_item.items() 
        if k != 'noninformative'
    ):
        alleleclass_item['noninformative'] = True

    return alleleclass_item


def make_alleleclass_item_readpluspair(
    aiitem_bnd1_rp1, aiitem_bnd1_rp2,
    aiitem_bnd2_rp1, aiitem_bnd2_rp2, 
):
    """Assumes that each rp supports at most 1 bnd (one of bnd1 or bnd2).
        This assumption is monitored in "update_alleleclass_sv" method
        of ReadPlus class.
    """

    alleleclass_item_bnd1 = init_rpp_alleleclass_item()
    alleleclass_item_bnd2 = init_rpp_alleleclass_item()

    #########
    # other #
    #########

    check_rpp_other_support(aiitem_bnd1_rp1, aiitem_bnd1_rp2, alleleclass_item_bnd1)
    check_rpp_other_support(aiitem_bnd2_rp1, aiitem_bnd2_rp2, alleleclass_item_bnd2)

    #######
    # ALT #
    #######

    # If an rpp is indirect-ALT-support, it cannot be direct ALT support
    indirect_support = (
        (
            aiitem_bnd1_rp1['facing_frombnd'] 
            and aiitem_bnd2_rp2['facing_frombnd']
        ) 
        or (
            aiitem_bnd1_rp2['facing_frombnd']
            and aiitem_bnd2_rp1['facing_frombnd']
        )
    )
    if indirect_support:
        alleleclass_item_bnd1['alt_support_indirect'] = True
        alleleclass_item_bnd2['alt_support_indirect'] = True
    else:
        check_rpp_direct_alt_support(
            aiitem_bnd1_rp1, aiitem_bnd1_rp2, alleleclass_item_bnd1,
            aiitem_bnd2_rp1, aiitem_bnd2_rp2, alleleclass_item_bnd2,
        )

    #######
    # ref #
    #######

    check_rpp_ref_support(aiitem_bnd1_rp1, aiitem_bnd1_rp2, alleleclass_item_bnd1)
    check_rpp_ref_support(aiitem_bnd2_rp1, aiitem_bnd2_rp2, alleleclass_item_bnd2)

    ##################
    # noninformative #
    ##################

    check_rpp_noninformative(alleleclass_item_bnd1)
    check_rpp_noninformative(alleleclass_item_bnd2)

    return alleleclass_item_bnd1, alleleclass_item_bnd2




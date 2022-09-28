import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
breakends = importlib.import_module('.'.join([top_package_name, 'sv', 'breakends']))
alleleinfosetup = importlib.import_module('.'.join([top_package_name, 'readplus', 'alleleinfosetup']))


DEFAULT_FLANKLEN_PARSIDE = 1
DEFAULT_FLANKLEN_BNDSIDE = 1
BND_DISTANCE_THRESHOLD = 1000


def make_alleleinfoitem_readpluspair(aiitem_bnd1_rp1, aiitem_bnd1_rp2,
                                     aiitem_bnd2_rp1, aiitem_bnd2_rp2, bnds):
    """Assumes that each rp supports at most 1 bnd (one of bnd1 or bnd2).
        This assumption is monitored in "update_alleleinfo_sv" method
        of ReadPlus class.
    """

    alleleinfoitem = {
        'ref_support_direct': False,
        'ref_support_indirect': False,
        'alt_support_direct': False,
        'alt_support_indirect': False,
        'other_support': False,
        'noninformative': False,
        }

    # other
    if any((aiitem_bnd1_rp1['other_support'], 
            aiitem_bnd2_rp1['other_support'],
            aiitem_bnd1_rp2['other_support'], 
            aiitem_bnd2_rp2['other_support'])):
        alleleinfoitem['other_support'] = True

    # alt
    if any((aiitem_bnd1_rp1['alt_support'], 
            aiitem_bnd2_rp1['alt_support'],
            aiitem_bnd1_rp2['alt_support'],
            aiitem_bnd2_rp2['alt_support'])):
        # direct alt support with read crossing
        alleleinfoitem['alt_support_direct'] = True
    elif (
            (
                aiitem_bnd1_rp1['facing_frombnd'] and 
                aiitem_bnd2_rp2['facing_frombnd']) or
            (
                aiitem_bnd1_rp2['facing_frombnd'] and 
                aiitem_bnd2_rp1['facing_frombnd'])):
        # the mate reads facing each other, on bndside of each breakend
        alleleinfoitem['alt_support_indirect'] = True

    # ref
    if not (alleleinfoitem['other_support'] or 
            alleleinfoitem['alt_support_direct'] or
            alleleinfoitem['alt_support_indirect']):
        # regarded not as ref-supporting when other- or alt- supporting
        if any((aiitem_bnd1_rp1['ref_support'], 
                aiitem_bnd2_rp1['ref_support'],
                aiitem_bnd1_rp2['ref_support'],
                aiitem_bnd2_rp2['ref_support'])):
            # direct ref support with read crossing
            alleleinfoitem['ref_support_direct'] = True
        elif (
                (
                    aiitem_bnd1_rp1['facing_frombnd'] and 
                    aiitem_bnd1_rp2['facing_frompar']) or
                (
                    aiitem_bnd1_rp2['facing_frombnd'] and 
                    aiitem_bnd1_rp1['facing_frompar']) or
                (
                    aiitem_bnd2_rp1['facing_frombnd'] and 
                    aiitem_bnd2_rp2['facing_frompar']) or
                (
                    aiitem_bnd2_rp2['facing_frombnd'] and 
                    aiitem_bnd2_rp1['facing_frompar'])):
            # the mate reads facing each other, one read on bndside of a 
                # breakend and the other read on parside of the same breakend
            alleleinfoitem['ref_support_indirect'] = True

    # noninformative
    if not any(v for (k, v) in alleleinfoitem.items() 
               if k != 'noninformative'):
        alleleinfoitem['noninformative'] = True

    return alleleinfoitem


def make_alleleinfoitem_readplus(bnds, is_bnd1, rp,
                                 flanklen_parside=DEFAULT_FLANKLEN_PARSIDE, 
                                 flanklen_bndside=DEFAULT_FLANKLEN_BNDSIDE):
    def handle_bndside_notonborder(rp, endis5, alleleinfoitem):
        if faces_border_frombnd(rp.read, endis5):
            if clipped_toward_border(rp.read, endis5):
                alleleinfoitem['noninformative'] = True
            else:
                alleleinfoitem['facing_frombnd'] = True
        else:
            alleleinfoitem['noninformative'] = True

    def handle_bndside_onborder_clipped(
            rp, pos_range0, endis5, alleleinfoitem, is_bnd1,
            flanklen_parside, bndside_flankrange0):
        queryseq_length = len(pos_range0) + flanklen_parside
        queryseq_beyond_bndsideflank = rp.get_seq_from_coord(
            start0=pos_range0.start, length=queryseq_length,
            forward=(not endis5))

        if len(queryseq_beyond_bndsideflank) < queryseq_length:
            alleleinfoitem['noninformative'] = True
        else:  # now flanking length conditions are fulfilled
            if rp.check_matches(bndside_flankrange0):
                seq_ref = bnds.get_seq_beyond_bnd_ref(is_bnd1=is_bnd1, 
                                                      length=flanklen_parside, 
                                                      with_border_seq=True)
                seq_alt = bnds.get_seq_beyond_bnd_alt(is_bnd1=is_bnd1, 
                                                      length=flanklen_parside,
                                                      with_border_seq=True)
                if queryseq_beyond_bndsideflank == seq_ref:
                    alleleinfoitem['ref_support'] = True
                elif queryseq_beyond_bndsideflank == seq_alt:
                    alleleinfoitem['alt_support'] = True
                else:
                    alleleinfoitem['other_support'] = True
            else:
                alleleinfoitem['other_support'] = True

    def handle_traverser(rp, pos_range0, endis5, flanklen_parside,
                         alleleinfoitem, bndside_flankrange0,
                         is_bnd1):
        queryseq_beyond_bnd = rp.get_seq_from_coord(
            start0=pos_range0.stop, length=flanklen_parside,
            forward=(not endis5))
        if len(queryseq_beyond_bnd) < flanklen_parside:
            alleleinfoitem['noninformative'] = True
        else:  # now flanking length conditions are fulfilled
            if (
                    rp.check_matches(bndside_flankrange0) and
                    rp.check_matches(pos_range0)):
                seq_beyond_bnd_ref = bnds.get_seq_beyond_bnd_ref(
                    is_bnd1=is_bnd1, length=flanklen_parside, 
                    with_border_seq=False)
                seq_beyond_bnd_alt = bnds.get_seq_beyond_bnd_alt(
                    is_bnd1=is_bnd1, length=flanklen_parside,
                    with_border_seq=False)

                if queryseq_beyond_bnd == seq_beyond_bnd_ref:
                    alleleinfoitem['ref_support'] = True
                elif queryseq_beyond_bnd == seq_beyond_bnd_alt:
                    alleleinfoitem['alt_support'] = True
                else:
                    alleleinfoitem['other_support'] = True
            else:
                alleleinfoitem['other_support'] = True

    # get parameters
    (chrom, pos_range0, endis5) = bnds.get_params(is_bnd1)
    bndside_flankrange0 = bnds.get_flank_range0(mode='bnd_dist', 
                                                is_bnd1=is_bnd1, 
                                                flanklen=flanklen_bndside)

    # initialize alleleinfoitem
    alleleinfoitem = {
        'noninformative': False,
        'ref_support': False,
        'alt_support': False,
        'other_support': False,
        'facing_frombnd': False,
        'facing_frompar': False,
        }

    # main
    if chrom == rp.read.reference_name:
        distance = rp.get_distance(pos_range0)
        if distance > BND_DISTANCE_THRESHOLD:
            alleleinfoitem['noninformative'] = True
        else:
            if on_bndside(rp.read, pos_range0, endis5):
                if on_bndside_notonborder(rp.read, pos_range0, endis5):
                    handle_bndside_notonborder(rp, endis5, alleleinfoitem)
                else:  # on bndside and on the border
                    if rp.check_spans(bndside_flankrange0):
                        if clipped_toward_border(rp.read, endis5):
                            handle_bndside_onborder_clipped(
                                rp, pos_range0, endis5, alleleinfoitem, 
                                is_bnd1, flanklen_parside, 
                                bndside_flankrange0)
                        else:
                            if faces_border_frombnd(rp.read, endis5):
                                alleleinfoitem['facing_frombnd'] = True
                            else:
                                alleleinfoitem['noninformative'] = True
                    else:
                        alleleinfoitem['noninformative'] = True

            elif on_parside(rp.read, pos_range0, endis5):
                if faces_border_frompar(rp.read, endis5):
                    alleleinfoitem['facing_frompar'] = True
                else:
                    alleleinfoitem['noninformative'] = True

            else:  # traverses the most advanced border
                if rp.check_spans(bndside_flankrange0):
                    handle_traverser(rp, pos_range0, endis5, flanklen_parside,
                                     alleleinfoitem, bndside_flankrange0,
                                     is_bnd1)
                else:
                    alleleinfoitem['noninformative'] = True
    else:
        alleleinfoitem['noninformative'] = True

    return alleleinfoitem


def on_bndside_notonborder(read, pos_range0, endis5):
    if endis5:
        return read.reference_start > pos_range0.start
    else:
        return read.reference_end <= pos_range0.start


def on_bndside(read, pos_range0, endis5):
    if endis5:
        return read.reference_start > pos_range0.stop
    else:
        return read.reference_end <= pos_range0.stop


# parside: partner side (named in conformity with the VCF spec)
def on_parside(read, pos_range0, endis5):
    if endis5:
        return read.reference_end <= pos_range0.stop + 1
    else:
        return read.reference_start >= pos_range0.stop


def clipped_toward_border(read, endis5):
    if endis5:
        return read.cigartuples[0][0] in (4, 5)
    else:
        return read.cigartuples[-1][0] in (4, 5)


def faces_border_frombnd(read, endis5):
    if endis5:
        return read.is_reverse
    else:
        return read.is_forward
   
    
def faces_border_frompar(read, endis5):
    if endis5:
        return read.is_forward
    else:
        return read.is_reverse










# helper functions

#def spans_bndsideflank(rp, bndsideflank_range0):
#    return (rp.ref_range0.start <= min(bndsideflank_range0) and
#            rp.ref_range0.stop > max(bndsideflank_range0))
#
#
#def spans_retracted_pos(rp, bndsideflank_range0):
#    return bndsideflank_range0.stop in rp.ref_range0
#
#
#def get_querypos0_beyond_bndsideflank(rp, bndsideflank_range0, endis5):
#    idx_retractedpos = rp.pairs_dict['refpos0'].index(
#        bndsideflank_range0.stop)
#
#    if endis5:
#        sl = slice(0, idx_retractedpos + 1)
#    else:
#        sl = slice(0, idx_retractedpos + 1)
#        querypos0_beyond_deep = rp.pairs_dict['querypos0'][idx_retractedpos:]
#
#    querypos0_beyond_deep = rp.pairs_dict['querypos0'][sl]
#    querypos0_beyond_deep = [x for x in querypos0_beyond_deep if x is not None]
#
#    return querypos0_beyond_deep
#
#
##def get_length_beyond_bndsideflank(pos_range0
#
#
#def spans_beyond_bnd(read, parflank_range0, endis5):
#    if endis5:
#        return read.reference_start <= parflank_range0[-1]
#    else:
#        return read.reference_end > parflank_range0[0]
#
#
#def spans_beyond_parflank(read, parflank_range0, endis5):
#    if endis5:
#        return read.reference_start <= parflank_range0[0]
#    else:
#        return read.reference_end > parflank_range0[-1]
#
#
#def get_length_otherside(read, pos_range0, endis5):
#    # assumes that the read does not span beyond the breakend
#    if endis5:
#        if read.cigartuples[0][0] == 4:
#            remaining_length = read.cigartuples[0][1]
#        else:
#            remaining_length = 0
#        length_thisside = pos_range0.stop - read.reference_start - 1
#    else:
#        if read.cigartuples[-1][0] == 4:
#            remaining_length = read.cigartuples[-1][1]
#        else:
#            remaining_length = 0
#        length_thisside = pos_range0.stop - read.reference_end
#
#    length_otherside = remaining_length - length_thisside
#
#    return length_otherside
#
#
#def clippedat_border(read, endis5):
#    if endis5:
#        return read.cigartuples[0][0] == 4
#    else:
#        return read.cigartuples[-1][0] == 4


#def get_supporting_bnds(rp):
#
#    def get_suppl_info(rp, SAitem):
#        """suppl_is5prime: 
#            Indicates whether the aligned portion in the supplementary 
#            alignment belongs to the 5-prime side (toward the first base) 
#            of the read.
#        suppl_isleft:
#            Indicates whether the aligned portion in the supplementary 
#            alignment belongs to the leftmost side, with regard to the
#            plus strand of the reference sequence.
#        """
#
#        if SAitem['is_forward']:
#            if (
#                    SAitem['cigartuples'][0][0] == 4 and 
#                    SAitem['cigartuples'][-1][0] == 0):
#                suppl_is5prime = False
#                suppl_isleft = False
#            elif (
#                    SAitem['cigartuples'][0][0] == 0 and 
#                    SAitem['cigartuples'][-1][0] == 4):
#                suppl_is5prime = True
#                suppl_isleft = True
#        else:
#            if (
#                    SAitem['cigartuples'][0][0] == 4 and 
#                    SAitem['cigartuples'][-1][0] == 0):
#                suppl_is5prime = True
#                suppl_isleft = False
#            elif (
#                    SAitem['cigartuples'][0][0] == 0 and 
#                    SAitem['cigartuples'][-1][0] == 4):
#                suppl_is5prime = False
#                suppl_isleft = True
#
#        return suppl_is5prime, suppl_isleft
#
#    def get_primary_clip_info(rp, suppl_is5prime):
#        """The softclip length in the primary alignment, 
#            corresponding to the supplementary alignment.
#        """
#
#        if suppl_is5prime:
#            clip_cigartuple = (rp.read.cigartuples[0]
#                               if rp.read.is_forward else
#                               rp.read.cigartuples[-1])
#        else:
#            clip_cigartuple = (rp.read.cigartuples[-1]
#                               if rp.read.is_forward else
#                               rp.read.cigartuples[0])
#
#        if clip_cigartuple[0] != 4:
#            raise Exception(
#                f'The cigarop of the primary read on the SA side is not '
#                f'a softclip.')
#
#        primary_cliplen = clip_cigartuple[1]
#
#        return primary_cliplen
#
#    # main
#    if rp.SAlist is None:
#        return None
#    
#    for SAitem in rp.SAlist:
#        suppl_is5prime, suppl_isleft = get_suppl_info(rp, SAitem)
#        primary_cliplen = get_primary_cliplen(rp, suppl_is5prime)
                
                
        
        
        

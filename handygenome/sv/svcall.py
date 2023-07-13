
import pyranges as pr
import numpy as np
import Bio.Seq
import Bio.SeqRecord

import handygenome.common as common
import handygenome.read.readplus as libreadplus
import handygenome.align.msa as libmsa
import handygenome.align.alignhandler as alignhandler


def call_breakends(bam, chrom, start0, end0, refver, SA_region_factor=1):
    # make rpplist
    rpplist = libreadplus.ReadPlusPairList.from_bam(bam, chrom, start0, end0)

    # set params
    chromdict = common.DEFAULT_CHROMDICTS[refver]
    rp_dict = dict()
    clipspec_list = list()
    for rp, clipspec in rpplist.iter_clipspecs():
        rp_dict[rp.uid] = rp
        clipspec_list.append(clipspec)

    # get consensus sequences
    consensus_result = get_consensus_from_clipspecs(clipspec_list)

    # parse SA realignment regions
    SA_region_info = dict()
    for clipkey, sublist in consensus_result.items():
        # sanity check
        all_consensus_list = [x['consensus'] for x in sublist]
        assert len(all_consensus_list) == len(set(all_consensus_list))

        for consensus_item in sublist:
            rp_list = list()
            for seqrec in consensus_item['msa']:
                uid = seqrec.annotations['clipspec'].readuid
                rp = rp_dict[uid]
                rp_list.append(rp)

            extended_SA_regions, SAregion_is_reversed = get_extended_SA_regions(
                rp_list, chromdict, SA_region_factor,
            )

            key = clipkey + (consensus_item['consensus'],)
            SA_region_info[key] = {
                'confident': extended_SA_regions['confident'],
                'inconfident': extended_SA_regions['inconfident'],
                'is_reversed': SAregion_is_reversed,
            }

    return SA_region_info, consensus_result


def get_extended_SA_regions(rp_list, chromdict, SA_region_factor):
    # collect SA region intervals
    confident_SA_regions = list()
    all_SA_regions = list()
    SAregion_is_reversed_data = set()

    for rp in rp_list:
        n_clip = sum((x[0] == 4) for x in rp.read.cigartuples)
        is_confident = ((n_clip == 1) and (len(rp.SAinfo) == 1))
            # This SAitem confidently matches the softclip

        for SAitem in rp.SAinfo:
            aligned_length = alignhandler.get_target_length(SAitem['cigartuples'])
            _start0 = SAitem['pos'] - 1
            _end0 = _start0 + aligned_length
            SA_region = (SAitem['chrom'], _start0, _end0)

            all_SA_regions.append(SA_region)
            if is_confident:
                confident_SA_regions.append(SA_region)
                is_reversed = (SAitem['is_forward'] != rp.read.is_forward)
                SAregion_is_reversed_data.add(is_reversed)

    # determine if SA region is reversed or not
    if len(SAregion_is_reversed_data) == 1:
        SAregion_is_reversed = SAregion_is_reversed_data.pop()
    else:  # 0 or >=2
        SAregion_is_reversed = None

    # merge SA regions, considering confident/nonconfident regions
    merged_SA_regions = {'confident': list(), 'inconfident': list()}
    if len(all_SA_regions) == 0:
        pass
    else:
        # make into clustered grs
        chroms, start0s, end0s = zip(*all_SA_regions)
        all_gr = pr.PyRanges(chromosomes=chroms, starts=start0s, ends=end0s).merge()

        if confident_SA_regions:
            chroms, start0s, end0s = zip(*confident_SA_regions)
            confident_gr = pr.PyRanges(chromosomes=chroms, starts=start0s, ends=end0s).merge()
            for selector in np.eye(confident_gr.df.shape[0]).astype(bool):
                part_confident_gr = confident_gr[selector]
                joined = part_confident_gr.join(all_gr, report_overlap=True)
                if not joined.empty:
                    max_row = max(
                        (x[1] for x in joined.df.iterrows()),
                        key=(lambda x: x['Overlap']),
                    )
                    SA_region = (max_row['Chromosome'], max_row['Start_b'], max_row['End_b'])
                    merged_SA_regions['confident'].append(SA_region)
        else:
            for idx, row in all_gr.df.iterrows():
                SA_region = (row['Chromosome'], row['Start'], row['End'])
                merged_SA_regions['inconfident'].append(SA_region)

    # extend merged SA regions
    def extend_region(region, SA_region_factor, chromdict):
        _chrom, _start0, _end0 = region

        length = _end0 - _start0
        pad = int(length * SA_region_factor)

        new_start0 = max(0, _start0 - pad)
        new_end0 = min(chromdict[_chrom], _end0 + pad)

        return (_chrom, new_start0, new_end0)

    extened_SA_regions = dict()
    extened_SA_regions['confident'] = [
        extend_region(x, SA_region_factor, chromdict)
        for x in merged_SA_regions['confident']
    ]
    extened_SA_regions['inconfident'] = [
        extend_region(x, SA_region_factor, chromdict)
        for x in merged_SA_regions['inconfident']
    ]

    # return
    return extened_SA_regions, SAregion_is_reversed



#############################################
# getting consensus from softclip fragments #
#############################################

def get_consensus_from_clipspecs(clipspec_list, aligner=alignhandler.ALIGNER_EQUAL_MM_GAP):
    # make SeqRecord objects from clipspecs
    seqrec_list = list()
    for clipspec in clipspec_list:
        seqrec_seq = (clipspec.seq if clipspec.is_forward else clipspec.seq[::-1])
        seqrec_id = f'(qname={clipspec.readuid.qname}, pos0={clipspec.pos0:,}, is_forward={clipspec.is_forward})'
        seqrec = Bio.SeqRecord.SeqRecord(
            seq=Bio.Seq.Seq(seqrec_seq), 
            id=seqrec_id, 
            annotations={'clipspec': clipspec},
        )
        seqrec_list.append(seqrec)

    # group seqrecs by clipspec identity
    clipid_groups = dict()
    for seqrec in seqrec_list:
        clipspec = seqrec.annotations['clipspec']
        clipid = (clipspec.readuid.chrom, clipspec.pos0, clipspec.is_forward)
        clipid_groups.setdefault(clipid, list())
        clipid_groups[clipid].append(seqrec)

    # make consensus seqs for each clipspec group
    result = dict()
    for clipid, sublist in clipid_groups.items():
        result[clipid] = list()
        seqrec_groups = group_seqrecs(sublist, aligner, score_limit=0.7)
        for target, subseq_info in seqrec_groups.items():
            target_seqrec = subseq_info['superseq']['seqrec']
            query_seqrecs = [x['seqrec'] for x in subseq_info['subseqs']]
            alignments = [x['alignment'] for x in subseq_info['subseqs']]

            # get consensus
            msa, consensus = libmsa.msa_with_target(
                target, alignments, 
                target_seqrec=target_seqrec, 
                query_seqrecs=query_seqrecs
            )
            true_consensus = (consensus if clipid[2] else consensus[::-1])

            # result
            result[clipid].append(
                {
                    'target': target, 
                    'msa': msa, 
                    'consensus': consensus, 
                    'true_consensus': true_consensus,
                }
            )

    return result


def group_seqrecs(seqrec_list, aligner, score_limit=0.7):
    """Groups SeqRecord objects according to sequence similarity"""
    def init_group(superseq_seqrec):
        return {
            'superseq': {
                'seqrec': superseq_seqrec, 'alignment': None, 'adj_score': None,
            },
            'subseqs': list(),
        }

    groups = dict()
    for seqrec in sorted(seqrec_list, key=len, reverse=True):
        seq = str(seqrec.seq)
        if len(groups) == 0:
            groups[seq] = init_group(seqrec)
            #groups[seqrec.seq].append({'seqrec': seqrec, 'alignment': None, 'adj_score': None})
            continue

        # find superseq candidates
        superseq_matching_info = list()
        for superseq in groups.keys():
            aln = alignhandler.alignment_tiebreaker(aligner.align(superseq, seqrec.seq))
            adj_score = score_alignment(aln)
            superseq_matching_info.append((superseq, aln, adj_score))

        # pick the most similar superseq
        chosen_info = max(superseq_matching_info, key=(lambda x: x[2]))
        superseq, aln, adj_score = chosen_info

        if adj_score >= score_limit:
            groups[superseq]['subseqs'].append(
                {'seqrec': seqrec, 'alignment': aln, 'adj_score': adj_score}
            )
        else:
            assert seq not in groups
            groups[seq] = init_group(seqrec)

    return groups


def score_alignment(alignment, match_score=1, mismatch_score=0, gap_score=0):
    # get walks
    walks = alignhandler.get_walks(alignment, copy=True)
    if not walks[0].check_both():
        del walks[0]
    if not walks[-1].check_both():
        del walks[-1]

    # calc score
    #target_seq = np.array(tuple(alignment.target))
    #query_seq = np.array(tuple(alignment.query))
    target_seq = np.fromiter(iter(alignment.target), dtype='<U1')
    query_seq = np.fromiter(iter(alignment.query), dtype='<U1')
    score = 0
    for walk in walks:
        if walk.check_both():
            target_subseq = target_seq[walk.get_target_slice()]
            query_subseq = query_seq[walk.get_query_slice()]

            n_match = (target_subseq == query_subseq).sum()
            n_mismatch = len(walk) - n_match

            score += n_match * match_score
            score += n_mismatch * mismatch_score
        else:
            score += len(walk) * gap_score

    # make adjusted score
    aligned_length = sum(len(x) for x in walks)
    adj_score = score / aligned_length

    return adj_score



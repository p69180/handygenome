
import pyranges as pr
import numpy as np
import Bio.Seq
import Bio.SeqRecord

import handygenome.refgenome.refgenome as refgenome
import handygenome.read.readplus as libreadplus
import handygenome.align.msa as libmsa
import handygenome.align.alignhandler as alignhandler


def call_breakends(bam, chrom, start0, end0, refver, SA_region_factor=1):
    # make rpplist
    rpplist = libreadplus.ReadPlusPairList.from_bam(bam, chrom, start0, end0)

    # set params
    chromdict = refgenome.get_chromdict(refver)
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

            confident_regions, unconfident_regions = get_SA_region_info(
                rp_list, chromdict, SA_region_factor,
            )

            key = clipkey + (consensus_item['true_consensus'],)
            SA_region_info[key] = {
                'confident': confident_regions,
                'unconfident': unconfident_regions,
                #'is_reversed': SAregion_is_reversed,
            }

    return SA_region_info, consensus_result


def get_SA_region_info(rp_list, chromdict, SA_region_factor):
    """Args:
        rp_list: list of ReadPlus objects which harbor a softclip of a given consensus
    """
    def make_gr_from_data(data):
        data_edit = tuple(zip(*data))
        gr = pr.from_dict(
            dict(Chromosome=data_edit[0], Start=data_edit[1], End=data_edit[2], is_reversed=data_edit[3])
        )

        # cluster
        gr = gr.cluster(slack=0)

        # merge by cluster; handle "is_reversed" values within a cluster
        merged_grs = list()
        for key, subdf in gr.df.groupby('Cluster'):
            flags = set(subdf['is_reversed'])
            if len(flags) == 1:
                is_reversed = flags.pop()
            else:
                is_reversed = None

            subgr = pr.PyRanges(subdf).merge()
            subgr.is_reversed = is_reversed
            merged_grs.append(subgr)

        # concat
        result = pr.concat(merged_grs)
        return result

    def extend_region(region, SA_region_factor, chromdict):
        """Modifies 'region' in-place"""
        length = region['end0'] - region['start0']
        pad = int(length * SA_region_factor)
        region['start0'] = max(0, region['start0'] - pad)
        region['end0'] = min(chromdict[region['chrom']], region['end0'] + pad)

    # collect SA region intervals
    confident_SA_regions = list()
    all_SA_regions = list()

    for rp in rp_list:
        n_clip = sum((x[0] == 4) for x in rp.read.cigartuples)
        is_confident = (n_clip == 1)  # This SAitem confidently matches the softclip

        for SAitem in rp.SAinfo:
            aligned_length = alignhandler.get_target_length(SAitem['cigartuples'])
            _start0 = SAitem['pos'] - 1
            _end0 = _start0 + aligned_length
            is_reversed = (SAitem['is_forward'] != rp.read.is_forward)
            SA_region = (SAitem['chrom'], _start0, _end0, is_reversed)

            all_SA_regions.append(SA_region)
            if is_confident:
                confident_SA_regions.append(SA_region)

    # merge SA regions, considering confident/nonconfident regions
    confident_regions = list()
    unconfident_regions = list()
    if all_SA_regions:
        # make into clustered grs
        all_gr = make_gr_from_data(all_SA_regions)
        if confident_SA_regions:
            confident_gr = make_gr_from_data(confident_SA_regions)
            for selector in np.eye(confident_gr.df.shape[0]).astype(bool):
                part_confident_gr = confident_gr[selector]
                joined = part_confident_gr.join(all_gr, report_overlap=True)
                if not joined.empty:
                    max_row = max(
                        (x[1] for x in joined.df.iterrows()),
                        key=(lambda x: x['Overlap']),
                    )
                    SA_region = {
                        'chrom': max_row['Chromosome'], 
                        'start0': max_row['Start_b'], 
                        'end0': max_row['End_b'], 
                        'is_reversed': max_row['is_reversed'],
                    }
                        # follows "is_reversed" value of the confident gr
                    confident_regions.append(SA_region)
        else:
            for idx, row in all_gr.df.iterrows():
                SA_region = {
                    'chrom': row['Chromosome'], 
                    'start0': row['Start'], 
                    'end0': row['End'], 
                    'is_reversed': row['is_reversed'],
                }
                unconfident_regions.append(SA_region)

    # extend merged SA regions
    for region in confident_regions:
        extend_region(region, SA_region_factor, chromdict)
    for region in unconfident_regions:
        extend_region(region, SA_region_factor, chromdict)

    return confident_regions, unconfident_regions



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



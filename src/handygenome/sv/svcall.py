import collections
import pprint

#import pyranges as pr
import numpy as np
import Bio.Seq
import Bio.SeqRecord

import handygenome.refgenome.refgenome as refgenome
from handygenome.read.readplus import ReadPlusPairList
import handygenome.align.msa as libmsa
import handygenome.align.alignhandler as alignhandler
from handygenome.genomedf.genomedf import GenomeDataFrame
from handygenome.sv.breakends import Breakends


#class ClipConsensusKey(
#    collections.namedtuple(
#        'ConsensusKey',
#        ('chrom', 'pos0', 'is5prime', 'seq')
#    )
#):
#    pass


def call_breakends(bam, chrom, start0, end0, SA_extend=30, mate_match_length=5, aligner=alignhandler.ALIGNER_EQUAL_MM_GAP):

    # make rpplist
    rpplist = ReadPlusPairList.from_bam(bam, chrom, start0, end0)

    # set params
    refver = refgenome.infer_refver_bamheader(bam.header)
    fasta = refgenome.get_fasta(refver)
    #chromdict = refgenome.get_chromdict(refgenome.infer_refver_bamheader(bam.header))

    clipspec_list = list()
    for rp, clipspec in rpplist.iter_clipspecs():
        clipspec_list.append(clipspec)

    # get consensus sequences
    consensus_result = get_consensus_from_clipspecs(
        clipspec_list, SA_extend=SA_extend,
    )

    # find alignment with each consensus
    alignments = list()
    for consensus_item in consensus_result:
        if consensus_item['SA_regions'] is None:
            continue
        aln_results = align_consensus_toSA(consensus_item, fasta, match_length=mate_match_length, aligner=aligner)
        valid_aln_results = [x for x in aln_results if x is not None]
        assert len(valid_aln_results) == 1
        consensus_item.update(valid_aln_results[0])

    # make breakends object
    bnds_list = list()
    for x in consensus_result:
        bnds = Breakends(
            bndspec1=(x['chrom'], x['pos0'], x['is5prime']),
            bndspec2=(x['mate_chrom'], x['mate_bnd_pos0'], x['mate_is5prime']),
            refver=refver,
            inserted_seq_view1=x['insseq_clipsideview'],
        )
        bnds_list.append(bnds)

    return consensus_result, bnds_list


########################################
# find alignment of consensus sequence #
########################################

def align_consensus_toSA(consensus_item, fasta, match_length=5, aligner=alignhandler.ALIGNER_EQUAL_MM_GAP):
    SA_gdf = consensus_item['SA_regions']
    result = list()
    for idx, row in SA_gdf.df.iterrows():
        target_chrom = row['Chromosome']
        target_start0 = row['Start']
        target_end0 = row['End']
        is_reversed = row['is_reversed']
        target = fasta.fetch(target_chrom, target_start0, target_end0)

        query = ( 
            Bio.Seq.reverse_complement(consensus_item['consensus'])
            if is_reversed else
            consensus_item['consensus']
        )
        tmp_flip_aln = (
            (not is_reversed) if consensus_item['is5prime'] else is_reversed
        )
        if tmp_flip_aln:
            aln = alignhandler.tiebreaker(
                aligner.align(target[::-1], query[::-1])
            )[:, ::-1]
        else:
            aln = alignhandler.tiebreaker(aligner.align(target, query))

        ###

        cigartuples, target_offset = alignhandler.alignment_to_cigartuples(
            aln, match_as_78=False, remove_left_del=False, remove_right_del=False,
        )
        long_matches = [(x[0] == 0 and x[1] >= match_length) for x in cigartuples]

        if not any(long_matches):
            subresult = None
        else:
            walks = list(
                alignhandler.walk_cigar(cigartuples, target_start0=target_start0)
            )
            flip_aln = (
                (not is_reversed)
                if consensus_item['is5prime'] else
                is_reversed
            )

            if flip_aln:  # counting from right
                lmidx = (len(long_matches) - 1) - long_matches[::-1].index(True)
                    # lmidx: long match index
                mate_bnd_pos0 = walks[lmidx].target.stop - 1
                querylen_before_longmatch = sum(len(x.query) for x in walks[lmidx + 1:])
                if querylen_before_longmatch == 0:
                    insseq = ''
                else:
                    insseq = query[-querylen_before_longmatch:]

            else:  # counting from left
                lmidx = long_matches.index(True)
                mate_bnd_pos0 = walks[lmidx].target.start
                querylen_before_longmatch = sum(len(x.query) for x in walks[:lmidx])
                insseq = query[:querylen_before_longmatch]
          
            # results
            if is_reversed:
                insseq_clipsideview = Bio.Seq.reverse_complement(insseq)
            else:
                insseq_clipsideview = insseq

            mate_is5prime = (
                is_reversed
                if consensus_item['is5prime'] else
                (not is_reversed)
            )
            subresult = {
                'mate_chrom': target_chrom,
                'mate_bnd_pos0': mate_bnd_pos0,
                'mate_is5prime': mate_is5prime,
                'insseq_clipsideview': insseq_clipsideview,
            }

        result.append(subresult)
    return result


##################
# SA tag parsing #
##################

def get_SA_region_info(rp_list, extend=30):
    """Args:
        rp_list: list of ReadPlus objects which harbor a softclip of a given consensus
    """
    refver = rp_list[0].get_refver()
    # step 1: make gdf from all SA regions
    chroms = list()
    start0s = list()
    end0s = list()
    is_reversed_list = list()

    for rp in rp_list:
        #if len(rp.SAinfo) > 0:
            #print(rp.read.query_name)
        for SAitem in rp.SAinfo:
            aligned_length = alignhandler.get_target_length(SAitem['cigartuples'])
            _start0 = SAitem['pos'] - 1
            _end0 = _start0 + aligned_length

            chroms.append(SAitem['chrom'])
            start0s.append(_start0)
            end0s.append(_end0)
            is_reversed_list.append(
                (SAitem['is_forward'] != rp.read.is_forward)
            )

    SA_gdf = GenomeDataFrame.from_data(
        refver=refver,
        chroms=chroms,
        start0s=start0s,
        end0s=end0s,
        is_reversed=is_reversed_list,
    )
    if SA_gdf.is_empty:
        #raise Exception(f'SA tag could not be found')
        return None
    else:
        SA_gdf.sort()

        # step 2: merge it and make is_reversed
        SA_gdf = SA_gdf.cluster(slack=(2 * extend))
        groupby = SA_gdf.df.groupby(SA_gdf['Cluster'])

        merged_df = groupby[['Chromosome', 'Start']].first()
        new_chroms = merged_df['Chromosome']
        new_start0s = merged_df['Start']
        new_end0s = groupby['End'].last()

        new_is_reversed_list = list()
        for key, subdf in groupby:
            values = subdf['is_reversed']
            if len(set(values)) == 1:
                new_is_reversed_list.append(values[0])
            else:
                new_is_reversed_list.append(None)
                
        merged_SA_gdf = GenomeDataFrame.from_data(
            refver=refver,
            chroms=new_chroms,
            start0s=new_start0s,
            end0s=new_end0s,
            is_reversed=new_is_reversed_list,
        )
        
        # sanitycheck
        if merged_SA_gdf.nrow != 1:
            raise Exception(f'More than 1 SA regions: {merged_SA_gdf}')

        # step 3: extend
        new_start0s = np.maximum(
            merged_SA_gdf.start0s - extend,
            0,
        )
        new_end0s = np.minimum(
            merged_SA_gdf.end0s + extend,
            merged_SA_gdf.get_chrom_lengths(),
        )
        merged_SA_gdf['Start'] = new_start0s
        merged_SA_gdf['End'] = new_end0s

        return merged_SA_gdf


#def get_SA_region_info_old(rp_list, chromdict, SA_region_factor):
#    """Args:
#        rp_list: list of ReadPlus objects which harbor a softclip of a given consensus
#    """
#    def make_gr_from_data(data):
#        data_edit = tuple(zip(*data))
#        gr = pr.from_dict(
#            dict(Chromosome=data_edit[0], Start=data_edit[1], End=data_edit[2], is_reversed=data_edit[3])
#        )
#
#        # cluster
#        gr = gr.cluster(slack=0)
#
#        # merge by cluster; handle "is_reversed" values within a cluster
#        merged_grs = list()
#        for key, subdf in gr.df.groupby('Cluster'):
#            flags = set(subdf['is_reversed'])
#            if len(flags) == 1:
#                is_reversed = flags.pop()
#            else:
#                is_reversed = None
#
#            subgr = pr.PyRanges(subdf).merge()
#            subgr.is_reversed = is_reversed
#            merged_grs.append(subgr)
#
#        # concat
#        result = pr.concat(merged_grs)
#        return result
#
#    def extend_region(region, SA_region_factor, chromdict):
#        """Modifies 'region' in-place"""
#        length = region['end0'] - region['start0']
#        pad = int(length * SA_region_factor)
#        region['start0'] = max(0, region['start0'] - pad)
#        region['end0'] = min(chromdict[region['chrom']], region['end0'] + pad)
#
#    # collect SA region intervals
#    confident_SA_regions = list()
#    all_SA_regions = list()
#
#    for rp in rp_list:
#        is_confident = (
#            sum((x[0] == 4) for x in rp.read.cigartuples) == 1
#        )  # When a read has only one softclip, its SA tag item is considered to confidently match the only softclip
#
#        for SAitem in rp.SAinfo:
#            aligned_length = alignhandler.get_target_length(SAitem['cigartuples'])
#            _start0 = SAitem['pos'] - 1
#            _end0 = _start0 + aligned_length
#            is_reversed = (SAitem['is_forward'] != rp.read.is_forward)
#            SA_region = (SAitem['chrom'], _start0, _end0, is_reversed)
#
#            all_SA_regions.append(SA_region)
#            if is_confident:
#                confident_SA_regions.append(SA_region)
#
#    # merge SA regions, considering confident/nonconfident regions
#    confident_regions = list()
#    unconfident_regions = list()
#    if all_SA_regions:
#        # make into clustered grs
#        all_gr = make_gr_from_data(all_SA_regions)
#        if confident_SA_regions:
#            confident_gr = make_gr_from_data(confident_SA_regions)
#            for selector in np.eye(confident_gr.df.shape[0]).astype(bool):
#                part_confident_gr = confident_gr[selector]
#                joined = part_confident_gr.join(all_gr, report_overlap=True)
#                if not joined.empty:
#                    max_row = max(
#                        (x[1] for x in joined.df.iterrows()),
#                        key=(lambda x: x['Overlap']),
#                    )
#                    SA_region = {
#                        'chrom': max_row['Chromosome'], 
#                        'start0': max_row['Start_b'], 
#                        'end0': max_row['End_b'], 
#                        'is_reversed': max_row['is_reversed'],
#                    }
#                        # follows "is_reversed" value of the confident gr
#                    confident_regions.append(SA_region)
#        else:
#            for idx, row in all_gr.df.iterrows():
#                SA_region = {
#                    'chrom': row['Chromosome'], 
#                    'start0': row['Start'], 
#                    'end0': row['End'], 
#                    'is_reversed': row['is_reversed'],
#                }
#                unconfident_regions.append(SA_region)
#
#    # extend merged SA regions
#    for region in confident_regions:
#        extend_region(region, SA_region_factor, chromdict)
#    for region in unconfident_regions:
#        extend_region(region, SA_region_factor, chromdict)
#
#    return confident_regions, unconfident_regions



#############################################
# getting consensus from softclip fragments #
#############################################

def get_consensus_from_clipspecs(
    clipspec_list, 
    cutoff_gt=1, 
    aligner=alignhandler.ALIGNER_EQUAL_MM_GAP,
    SA_extend=30,
):
    # make SeqRecord objects from clipspecs
    seqrec_list = list()
    for clipspec in clipspec_list:
        seqrec = Bio.SeqRecord.SeqRecord(
            seq=Bio.Seq.Seq(  # beginning of clip is on the left
                (clipspec.seq[::-1] if clipspec.is5prime else clipspec.seq)
            ), 
            id=clipspec.get_uid(), 
            annotations={'clipspec': clipspec},
        )
        seqrec_list.append(seqrec)

    # group seqrecs by clipspec identity (start position and direction of clip)
    clipid_groups = dict()
    for seqrec in seqrec_list:
        clipspec = seqrec.annotations['clipspec']
        clipid = (clipspec.chrom, clipspec.pos0, clipspec.is5prime)
        clipid_groups.setdefault(clipid, list())
        clipid_groups[clipid].append(seqrec)

    # make consensus seqs for each clipspec group
    result = list()
    for clipid, sublist in clipid_groups.items():  # this "sublist" contains clips with same positions and directions
        # filter seqrec groups with cutoff #1
        if len(sublist) <= cutoff_gt:
            continue

        seqrec_groups = group_seqrecs(sublist, aligner, score_limit=0.7)  # group similar sequences together ; multiple separate breakend events may share identical position and direction

        for target, subseq_info in seqrec_groups.items():  # seqrecs contained in "subseq_info" are of similar sequences
            target_seqrec = subseq_info['superseq']['seqrec']
            query_seqrecs = [x['seqrec'] for x in subseq_info['subseqs']]
            alignments = [x['alignment'] for x in subseq_info['subseqs']]

            # filter seqrec groups with cutoff #2
            all_seqrecs = [target_seqrec] + query_seqrecs
            all_qnames = [
                x.annotations['clipspec'].rp.read.query_name
                for x in all_seqrecs
            ]
            if len(set(all_qnames)) <= cutoff_gt:
                continue

            # get consensus
            msa, msa_consensus = libmsa.msa_with_target(
                target, alignments, 
                target_seqrec=target_seqrec, 
                query_seqrecs=query_seqrecs
            )
            consensus = (msa_consensus[::-1] if clipid[2] else msa_consensus)

            # make SA region GDF
            rp_list = [
                seqrec.annotations['clipspec'].rp
                for seqrec in msa
            ]
            SA_region_gdf = get_SA_region_info(rp_list, extend=SA_extend)

            # result
            result.append(
                {
                    'chrom': clipid[0],
                    'pos0': clipid[1],
                    'is5prime': clipid[2],
                    'msa_target': target, 
                    'msa': msa, 
                    'rp_list': rp_list,
                    'consensus': consensus, 
                    'msa_consensus': msa_consensus, 
                    'SA_regions': SA_region_gdf,
                }
            )

    return result


def group_seqrecs(seqrec_list, aligner, score_limit=0.7):
    """Groups SeqRecord objects according to sequence similarity"""
    def init_group(superseq_seqrec):
        return {
            'superseq': {
                'seqrec': superseq_seqrec, 
                'alignment': None, 
                'adj_score': None,
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



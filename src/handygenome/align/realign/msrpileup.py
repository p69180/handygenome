import sys
import collections
import itertools
import functools
import logging
import inspect
import random
import uuid

import pysam
import Bio.Align
import numpy as np
import pandas as pd
import pyranges as pr

import handygenome.refgenome.refgenome as refgenome
import handygenome.workflow as workflow
import handygenome.variant.vcfspec as libvcfspec
import handygenome.read.pileup as libpileup
import handygenome.align.alignhandler as alignhandler
import handygenome.bameditor as bameditor
import handygenome.read.readhandler as readhandler
import handygenome.read.readplus as readplus
import handygenome.align.realign.base as realign_base


class MultisampleRealignerPileup(realign_base.RealignerPileupBase):
    def __init__(
        self, 
        pileup_dict, 
        refver=None,
        fasta=None, 
        aligner=realign_base.DEFAULT_ALIGNER,
        verbose=False, 
        logger=None,
        **kwargs,
    ):
        """Args:
            pileup_dict: keys - sampleid; values - Pileup object
        """
        # set params
        self.pileup_dict = pileup_dict
        self.refver, self.fasta = realign_base.RealignerPileupBase.refver_fasta_arghandler(refver, fasta)
        self.aligner = aligner
        self.verbose, self.logger = realign_base.RealignerPileupBase.logger_arghandler(verbose, logger)
        self.params = realign_base.parse_rpileup_kwargs(**kwargs)

        # setup data structure
        self.set_df()
        self.set_row_spec_groups()
        self.save_subseq_alignments()
        self.allocate_subseq_alignments()
        self.set_row_spec_groups_bysample()

    @property
    def first_pileup(self):
        return next(iter(self.pileup_dict.values()))

    @property
    def chrom(self):
        return self.first_pileup.chrom

    @property
    def start0(self):
        return self.first_pileup.start0

    @property
    def end0(self):
        return self.first_pileup.end0

    @property
    def range0(self):
        return self.first_pileup.range0

    def set_df(self):
        assert len(set(x.range0 for x in self.pileup_dict.values())) == 1, f'Genomic ranges of pileup objects are different.'

        self.df = pd.concat(
            {key: val.seq_df for key, val in self.pileup_dict.items()},
            names=['SampleID', 'ReadUID'],
        )

    def row_spec_getter(self, key):
        return self.pileup_dict[key[0]].row_specs[key[1]]

    def set_row_spec_groups(self):
        # initialize row_spec_groups
        self.row_spec_groups = dict()
        candidates = dict()
        for superseq_sampleid, sub_pileup in self.pileup_dict.items():
            for superseq_readuid, groupinfo in sub_pileup.row_spec_groups.items():
                superseq_key = (superseq_sampleid, superseq_readuid)
                supserseq_row_spec = sub_pileup.row_specs[superseq_readuid]
                candidates.setdefault(supserseq_row_spec['seq'], list())
                candidates[supserseq_row_spec['seq']].append((superseq_key, supserseq_row_spec))

        for key, val in candidates.items():
            selected_superseq_key = max(
                val, 
                key=(lambda x: (x[1]['span_length'], -x[1]['cliplen'])),
            )[0]
            self.row_spec_groups[selected_superseq_key] = {
                'subseq_rowids': list(),
                'subseq_hits': 0,
            }

        # add subseq entries
        for subseq_sampleid, sub_pileup in self.pileup_dict.items():
            for subseq_readuid, subseq_row_spec in sub_pileup.row_specs.items():
                superseq_candidates = list()
                for superseq_key in self.row_spec_groups.keys():
                    superseq_sampleid, superseq_readuid = superseq_key
                    superseq_row_spec = self.pileup_dict[superseq_sampleid].row_specs[superseq_readuid]
                    if self.row_spec_matcher(subseq_row_spec, superseq_row_spec):
                        superseq_candidates.append(superseq_key)
                        
                if len(superseq_candidates) == 0:
                    raise Exception(f'Subseq row_spec does not match with any of superseq row_specs')
                else:
                    superseq_key = self.superseq_candidates_tiebreak(superseq_candidates, self.row_spec_getter, self.row_spec_groups)
                    groupinfo = self.row_spec_groups[superseq_key]
                    groupinfo['subseq_rowids'].append((subseq_sampleid, subseq_readuid))
                    groupinfo['subseq_hits'] += 1

        # discard 0-hit groups
        no_hits = [
            subseq_key for subseq_key, groupinfo in self.row_spec_groups.items()
            if groupinfo['subseq_hits'] == 0
        ]
        for subseq_key in no_hits:
            del self.row_spec_groups[subseq_key]

        self.decorate_groupinfo(self.row_spec_groups, self.read_getter)

    def read_getter(self, key):
        return self.pileup_dict[key[0]].read_store[key[1]]

    def save_subseq_alignments(self):
        self._subseq_alignments = dict()
        for superseq_key, groupinfo in self.row_spec_groups.items():
            sup_to_ref_aln = self.pileup_dict[superseq_key[0]]._superseq_alignments[superseq_key[1]]
            for subseq_key in groupinfo['subseq_rowids']:
                subseq_row_spec = self.pileup_dict[subseq_key[0]].row_specs[subseq_key[1]]
                try:
                    self._subseq_alignments[subseq_key] = self.align_subseq_to_ref(subseq_row_spec, sup_to_ref_aln)
                except:
                    print(superseq_key)
                    print(subseq_key)
                    raise

    def allocate_subseq_alignments(self):
        subseq_alignments_bysample = {sampleid: dict() for sampleid in self.pileup_dict.keys()}
        for subseq_key, aln in self._subseq_alignments.items():
            subseq_alignments_bysample[subseq_key[0]][subseq_key[1]] = aln

        for sampleid, sub_pileup in self.pileup_dict.items():
            sub_pileup._subseq_alignments = subseq_alignments_bysample[sampleid]

    def set_row_spec_groups_bysample(self):
        self.row_spec_groups_bysample = {sampleid: dict() for sampleid in self.pileup_dict.keys()}
        for sampleid, subgroups in self.row_spec_groups_bysample.items():
            for superseq_key, groupinfo in self.row_spec_groups.items():
                subgroups[superseq_key] = {
                    'subseq_rowids': list(),
                    'subseq_hits': 0,
                }
            
        for superseq_key, groupinfo in self.row_spec_groups.items():
            for subseq_key in groupinfo['subseq_rowids']:
                subgroups = self.row_spec_groups_bysample[subseq_key[0]]
                subgroups[superseq_key]['subseq_rowids'].append(subseq_key)
                subgroups[superseq_key]['subseq_hits'] += 1

        for sampleid, subgroups in self.row_spec_groups_bysample.items():
            self.decorate_groupinfo(subgroups, self.read_getter)

    def iter_contig_vcfspecs(self, subseq_portion_threshold, MQ_threshold, threshold_bysample=True, verbose=False):
        if threshold_bysample:
            for superseq_key, groupinfo in sorted(
                self.row_spec_groups.items(),
                key=(lambda x: x[1]['subseq_hits']),
                reverse=True,
            ):
                if any(
                    (subgroups[superseq_key]['subseq_hits_portion'] >= subseq_portion_threshold) and
                    (subgroups[superseq_key]['mean_MQ'] >= MQ_threshold)
                    for subgroups in self.row_spec_groups_bysample.values()
                ):
                    if verbose:
                        print('vcfspec:', self.pileup_dict[superseq_key[0]].contig_vcfspecs[superseq_key[1]])
                        print('vcfspec components:', self.pileup_dict[superseq_key[0]].contig_vcfspecs[superseq_key[1]].components)
                        for sampleid, subgroups in self.row_spec_groups_bysample.items():
                            print('sample ID:', sampleid, 'portion:', subgroups[superseq_key]['subseq_hits_portion'], 'meanMQ:', subgroups[superseq_key]['mean_MQ'])
                        print()

                    yield superseq_key, self.pileup_dict[superseq_key[0]].contig_vcfspecs[superseq_key[1]]
        else:
            for superseq_key, groupinfo in sorted(
                self.row_spec_groups.items(),
                key=(lambda x: x[1]['subseq_hits']),
                reverse=True,
            ):
                if (
                    (groupinfo['subseq_hits_portion'] >= subseq_portion_threshold) and
                    (groupinfo['mean_MQ'] >= MQ_threshold)
                ):
                    if verbose:
                        print('vcfspec:', self.pileup_dict[superseq_key[0]].contig_vcfspecs[superseq_key[1]])
                        print('vcfspec components:', self.pileup_dict[superseq_key[0]].contig_vcfspecs[superseq_key[1]].components)
                        print('portion:', groupinfo['subseq_hits_portion'], 'meanMQ:', groupinfo['mean_MQ'])
                        print()

                    yield superseq_key, self.pileup_dict[superseq_key[0]].contig_vcfspecs[superseq_key[1]]

    # get_result_vcfspec BEGIN #

    def get_result_vcfspec_recalc_hits(self, vaf_cutoff=0.01, verbose=False):
        candidate_monoalts = list()
        for superseq_key in self.row_spec_groups.keys():
            contig_vcfspec = self.pileup_dict[superseq_key[0]].contig_vcfspecs[superseq_key[1]]
            if contig_vcfspec is not None:
                candidate_monoalts.append(contig_vcfspec.normalize())
        if len(candidate_monoalts) == 0:
            return None

        # create a range which covers all candidate mono-ALT vcfspecs
        merged_ref_range0 = range(
            min(x.start0 for x in candidate_monoalts),
            max(x.end0 for x in candidate_monoalts),
        )
        # make rpplist using the range with each sample
        rpplist_dict = dict()
        for sampleid, pileup in self.pileup_dict.items():
            rpplist_dict[sampleid] = readplus.get_rpplist_nonsv(
                bam=pileup.bam, 
                fasta=self.fasta, 
                chromdict=refgenome.get_chromdict(self.refver),
                chrom=self.chrom,
                start0=merged_ref_range0.start,
                end0=merged_ref_range0.stop,
            )
        # calculate vaf for each mono-ALT vcfspec from rpplists
        # choose if vaf is greater than cutoff in any sample
        selected_monoalt_vcfspecs = set()
        for monoalt_vcfspec in candidate_monoalts:
            if verbose:
                print()
                print(monoalt_vcfspec)

            for sampleid, rpplist in rpplist_dict.items():
                rpplist.update_alleleclass(monoalt_vcfspec)
                alleleclass_counter = collections.Counter(
                    rpp.alleleclass[monoalt_vcfspec] for rpp in rpplist
                )
                total_reads = sum(
                    v for (k, v) in alleleclass_counter.items()
                    if k is not None
                )
                var_reads = alleleclass_counter[1]
                if total_reads == 0:
                    vaf = None
                else:
                    vaf = var_reads / total_reads

                if verbose:
                    print(f'sampleid {sampleid}')
                    print(f'vaf {vaf}')
                    
                if (vaf is not None) and (vaf >= vaf_cutoff):
                    selected_monoalt_vcfspecs.add(monoalt_vcfspec)
                    break

        if len(selected_monoalt_vcfspecs) == 0:
            return None

        # sort selected vcfspecs
        selected_monoalt_vcfspecs = sorted(
            selected_monoalt_vcfspecs,
            key=(lambda x: (x.pos, x.alts)),
        )
        # merge and return
        return functools.reduce(libvcfspec.merge, selected_monoalt_vcfspecs)

    def get_result_vcfspec_recalc_hits_dual_filtering(
        self, vaf_cutoff=0.01, subseq_portion_threshold=None, 
        MQ_threshold=40, threshold_bysample=True,
        verbose=False,
    ):
        # set params
        if subseq_portion_threshold is None:
            subseq_portion_threshold = self.params['allele_portion_threshold']

        # 1st filtering
        # make candidate mono-ALT vcfspecs
        candidate_monoalts = list()
        for x in self.iter_contig_vcfspecs(
            subseq_portion_threshold=subseq_portion_threshold, 
            MQ_threshold=MQ_threshold, 
            threshold_bysample=threshold_bysample,
        ):
            if x[1] is not None:
                candidate_monoalts.append(x[1].normalize())

        if len(candidate_monoalts) == 0:
            return None

        # 2nd filtering
        # create rpplists for each component sample
        merged_ref_range0 = range(
            min(x.start0 for x in candidate_monoalts),
            max(x.end0 for x in candidate_monoalts),
        )
        rpplist_dict = dict()
        for sampleid, pileup in self.pileup_dict.items():
            rpplist_dict[sampleid] = readplus.get_rpplist_nonsv(
                bam=pileup.bam, 
                fasta=self.fasta, 
                chromdict=refgenome.get_chromdict(self.refver),
                chrom=self.chrom,
                start0=merged_ref_range0.start,
                end0=merged_ref_range0.stop,
            )
        # calculate vaf for each mono-ALT vcfspec from rpplists
        # choose if vaf is greater than cutoff in any sample
        selected_monoalt_vcfspecs = set()
        for monoalt_vcfspec in candidate_monoalts:
            if verbose:
                print()
                print(monoalt_vcfspec)

            for sampleid, rpplist in rpplist_dict.items():
                rpplist.update_alleleclass(monoalt_vcfspec)
                alleleclass_counter = collections.Counter(
                    rpp.alleleclass[monoalt_vcfspec] for rpp in rpplist
                )
                total_reads = sum(
                    v for (k, v) in alleleclass_counter.items()
                    if k is not None
                )
                var_reads = alleleclass_counter[1]
                if total_reads == 0:
                    vaf = None
                else:
                    vaf = var_reads / total_reads

                if verbose:
                    print(f'sampleid {sampleid}')
                    print(f'vaf {vaf}')
                    
                if (vaf is not None) and (vaf >= vaf_cutoff):
                    selected_monoalt_vcfspecs.add(monoalt_vcfspec)
                    break

        if len(selected_monoalt_vcfspecs) == 0:
            return None

        # sort selected vcfspecs
        selected_monoalt_vcfspecs = sorted(
            selected_monoalt_vcfspecs,
            key=(lambda x: (x.pos, x.alts)),
        )
        # merge and return
        return functools.reduce(libvcfspec.merge, selected_monoalt_vcfspecs)

    def get_result_vcfspec_use_rowspecgroup_hits(
        self, 
        #as_components=False, 
        #merge=True, 
        subseq_portion_threshold=None, 
        MQ_threshold=30, 
        threshold_bysample=True,
        verbose=False,
    ):
        if subseq_portion_threshold is None:
            subseq_portion_threshold = self.params['allele_portion_threshold']

#        if as_components:
#            result = list()
#            for row_id, contig_vcfspec in self.iter_contig_vcfspecs(subseq_portion_threshold=subseq_portion_threshold, MQ_threshold=MQ_threshold, threshold_bysample=threshold_bysample):
#                if contig_vcfspec is None:
#                    continue
#
#                if len(contig_vcfspec.components[0]) == 0:
#                    result.append(contig_vcfspec)
#                else:
#                    result.extend(contig_vcfspec.components[0])
        result = list(
            x[1] for x in self.iter_contig_vcfspecs(
                subseq_portion_threshold=subseq_portion_threshold, 
                MQ_threshold=MQ_threshold, 
                threshold_bysample=threshold_bysample,
                verbose=verbose,
            )
            if x[1] is not None
        )

        if len(result) == 0:
            return None
        else:
            return functools.reduce(libvcfspec.merge, result)

    # SELECT AN APPROPRIATE FUNCTION
    get_result_vcfspec = get_result_vcfspec_use_rowspecgroup_hits

    # get_result_vcfspec END #

    def show_row_spec_alignments(self, skip_zero_hits=True, show_subseqs=False):
        self._show_row_spec_alignments_helper(
            self.row_spec_groups, 
            lambda x: self.pileup_dict[x[0]].row_specs.__getitem__(x[1]),
            lambda x: self.pileup_dict[x[0]]._superseq_alignments.__getitem__(x[1]),
            skip_zero_hits, 
            show_subseqs,
        )

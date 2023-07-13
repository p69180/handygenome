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

import handygenome.common as common
import handygenome.workflow as workflow
import handygenome.variant.vcfspec as libvcfspec
import handygenome.read.pileup as libpileup
import handygenome.align.alignhandler as alignhandler
import handygenome.bameditor as bameditor
import handygenome.read.readhandler as readhandler
import handygenome.read.readplus as readplus
import handygenome.align.realign.base as realign_base
import handygenome.align.realign.msrpileup as realign_msrpileup
import handygenome.align.realign.rpileup_series as realign_rpileup_series


class MultisampleRealignerPileupSeries:
    def __repr__(self):
        return str(self.pileupseries_dict)

    def iter_series(self):
        return iter(self.pileupseries_dict.values())

    def iter_mspileups(self):
        return iter(self.mspileup_list)

    @functools.cached_property
    def df(self):
        data = list()
        index = list()
        for idx, (sampleid, rpileup_series) in enumerate(self.pileupseries_dict.items()):
            index.append(sampleid)
            data.append(list(rpileup_series.pileup_list))
            if idx == 0:
                breaks = list()
                breaks.append(rpileup_series.pileup_list[0].start0)
                for rpileup in rpileup_series.pileup_list:
                    breaks.append(rpileup.end0)
                columns = pd.IntervalIndex.from_breaks(breaks, closed='left')
        return pd.DataFrame.from_records(data, index=index, columns=columns)

#    @property
#    def rows(self):
#        return pd.DataFrame.from_dict(self.pileupseries_dict, orient='index')
#
#    @property
#    def cols(self):
#        return pd.DataFrame.from_records(
#            [
#                {
#                    (x.start0, x.end0): x 
#                    for x in self.mspileup_list
#                }
#            ]
#        )

    #@property
    #def series_length(self):
        #return len(self.rows.iloc[0, :])
    #    return len(next(iter(self.pileupseries_dict.values())))

    @property
    def ranges(self):
        #return [range(*x) for x in self.cols.columns]
        return [x.range0 for x in next(iter(self.pileupseries_dict.values())).pileup_list]

    @property
    def start0(self):
        #return self.rows.iloc[0, :].start0
        return next(self.pileupseries_dict.values()).start0

    @property
    def end0(self):
        #return self.rows.iloc[0, :].end0
        return next(self.pileupseries_dict.values()).end0

    def write_realigned_reads(self, bam_path_dict, padding=0):
        for sampleid, rpileup_ser in self.pileupseries_dict.items():
            rpileup_ser.write_realigned_reads(bam_path_dict[sampleid], padding=padding)

    def iter_contig_vcfspecs(self, subseq_portion_threshold=0, MQ_threshold=40):
        return itertools.chain.from_iterable(
            msrpileup.iter_contig_vcfspecs(subseq_portion_threshold=subseq_portion_threshold, MQ_threshold=MQ_threshold)
            for msrpileup in self.mspileup_list
        )

    def get_result_vcfspecs(self, portion_cutoff=0.1, MQ_cutoff=30, verbose=False):
        return list(
            msrpileup.get_result_vcfspec(
                subseq_portion_threshold=portion_cutoff, 
                MQ_threshold=MQ_cutoff, 
                threshold_bysample=True,
                verbose=verbose,
            )
            for msrpileup in self.mspileup_list
        )

    def get_result_vcfspecs_old(
        self, vaf_cutoff=0.01, subseq_portion_threshold=None, 
        MQ_threshold=40, threshold_bysample=True,
        verbose=False,
    ):
        return list(
            msrpileup.get_result_vcfspec_old(
                vaf_cutoff=vaf_cutoff, 
                subseq_portion_threshold=subseq_portion_threshold, 
                MQ_threshold=MQ_threshold, 
                threshold_bysample=threshold_bysample,
                verbose=verbose,
            )
            for msrpileup in self.mspileup_list
        )

    def get_result_vcfspecs_old(
        self, as_components=False, merge=True, subseq_portion_threshold=None, MQ_threshold=40, threshold_bysample=True,
    ): 
        return list(
            msrpileup.get_result_vcfspecs(
                as_components=as_components, 
                merge=merge, 
                subseq_portion_threshold=subseq_portion_threshold, 
                MQ_threshold=MQ_threshold, 
                threshold_bysample=threshold_bysample,
            )
            for msrpileup in self.mspileup_list
        )

    ###########################
    # initializer and helpers #
    ###########################
    def __init__(
        self, 
        bam_dict, 
        chrom, 
        start0, 
        end0, 
        refver=None,
        fasta=None,
        aligner=realign_base.DEFAULT_ALIGNER,
        verbose=False, 
        logger=None, 
        init_blank=False,
        **kwargs,
    ):
        # set params
        self.chrom = chrom
        self.bam_dict = bam_dict
        self.refver, self.fasta = realign_base.RealignerPileupBase.refver_fasta_arghandler(refver, fasta)
        self.aligner = aligner
        self.verbose, self.logger = realign_base.RealignerPileupBase.logger_arghandler(verbose, logger)
        self.params = realign_base.parse_rpileup_kwargs(**kwargs)

        # setup data
        if not init_blank:
            self.set_series_dict(start0, end0)  # self.pileupseries_dict
            self.set_multisample_pileups()  # self.mspileup_list
            self.set_realigned_reads()

    def set_series_dict(self, seed_start0, seed_end0):
        self.pileupseries_dict = dict()
        self.no_variant = False

        # initialize
        for sampleid, bam in self.bam_dict.items():
            self.logger.debug(f'@@@ Initializing RealignerPileupSeries of sample {sampleid} @@@\n')
            self.pileupseries_dict[sampleid] = realign_rpileup_series.RealignerPileupSeries(
                bam=bam, 
                chrom=self.chrom, start0=seed_start0, end0=seed_end0, 
                refver=self.refver,
                fasta=self.fasta, 
                aligner=self.aligner,
                verbose=self.verbose,
                logger=self.logger,
                **self.params,
            )
            self.logger.debug(f'@@@ Finished initialization of RealignerPileupSeries of sample {sampleid} @@@\n\n')

        # equalize whole margins
        self.equalize_left()
        self.equalize_right()

        # equalize sub-pileup margins
        self.logger.debug(f'@@@ Beginning inner margin equalization @@@\n')
        self.equalize_inner_margins()
        self.logger.debug(f'@@@ Finished inner margin equalization @@@\n')

    def equalize_left(self):
        self.logger.debug(f'@@@ Beginning equalize_left @@@\n')

        # equalize left
        while True:
            if len(set(x.start0 for x in self.pileupseries_dict.values())) == 1:
                break
            # extend
            target_start0 = min(x.start0 for x in self.pileupseries_dict.values())
            self.logger.debug(f'target_start0: {target_start0}')
            for sampleid, pileup_ser in self.pileupseries_dict.items():
                width = pileup_ser.start0 - target_start0
                if width > 0:
                    self.logger.debug(f'Beginning EXTEND of {sampleid}')
                    self.logger.debug(f'current start0 of {sampleid}: {pileup_ser.start0}')
                    pileup_ser.extend_left(width)
                    self.logger.debug(f'Finished EXTEND of {sampleid}\n')

            #for sampleid, pileup_ser in self.pileupseries_dict.items():
                #print(sampleid)
                #print(pileup_ser.pileup_list[0].df)
            #    print(pileup_ser.pileup_list[0].active_info)
            #    print()
            #print('----------------')

            # check if hit width limit
            if any(
                pileup_ser.width >= self.params['max_series_width']
                for pileup_ser in self.pileupseries_dict.values()
            ):
                break
            # secure
            for sampleid, pileup_ser in self.pileupseries_dict.items():
                self.logger.debug(f'Beginning SECURE of {sampleid}')
                pileup_ser.secure_left()
                self.logger.debug(f'Finished SECURE of {sampleid}\n')

        self.logger.debug(f'@@@ Finished equalize_left @@@\n')

    def equalize_right(self):
        self.logger.debug(f'@@@ Beginning equalize_right @@@\n')

        while True:
            if len(set(x.end0 for x in self.pileupseries_dict.values())) == 1:
                break
            # extend
            target_end0 = max(x.end0 for x in self.pileupseries_dict.values())
            self.logger.debug(f'target_end0: {target_end0}')
            for sampleid, pileup_ser in self.pileupseries_dict.items():
                width = target_end0 - pileup_ser.end0
                if width > 0:
                    self.logger.debug(f'Beginning EXTEND of {sampleid}')
                    self.logger.debug(f'current end0 of {sampleid}: {pileup_ser.end0}')
                    pileup_ser.extend_right(width)
                    self.logger.debug(f'Finished EXTEND of {sampleid}\n')

            for sampleid, pileup_ser in self.pileupseries_dict.items():
                print(sampleid)
                #print(pileup_ser.pileup_list[-1].df)
                print(pileup_ser.pileup_list[-1].active_info)
                print()
            print('----------------')

            # check if hit width limit
            if any(
                pileup_ser.width >= self.params['max_series_width'] 
                for pileup_ser in self.pileupseries_dict.values()
            ):
                break
            # secure
            for pileup_ser in self.pileupseries_dict.values():
                self.logger.debug(f'Beginning SECURE of {sampleid}')
                pileup_ser.secure_right()
                self.logger.debug(f'Finished SECURE of {sampleid}\n')

        self.logger.debug(f'@@@ Finished equalize_right @@@\n')

    def equalize_inner_margins(self):
        # set interim parameters
        first_pileupseries = next(self.iter_series())
        series_start0 = first_pileupseries.start0
        series_end0 = first_pileupseries.end0
        series_width = series_end0 - series_start0
        max_pileup_width = first_pileupseries.pileup_list[0].params['max_pileup_width']

        # When there is no need to further split the series range
        if series_end0 - series_start0 <= max_pileup_width:
            start0_list = [series_start0, series_end0]
            for pileup_ser in self.pileupseries_dict.values():
                if len(pileup_ser.pileup_list) > 1:
                    pileup_ser.rearrange(start0_list, prepare_vcfspecs=True)
            self.rearrange_mode = 'not_done'
            return

        # best case - without splitting contig vcfspec
        start0_list, all_splittable = self.get_rearrangement_points(
            split_region_gr=functools.reduce(
                lambda x, y: x.intersect(y), 
                (pileup_ser.get_splittable_region_best() for pileup_ser in self.pileupseries_dict.values())
            ), 
            trim_margins=True, 
            max_pileup_width=max_pileup_width, 
            series_start0=series_start0, 
            series_end0=series_end0,
        )

        if start0_list is not None:
            for pileup_ser in self.pileupseries_dict.values():
                pileup_ser.rearrange(start0_list, prepare_vcfspecs=True)
            self.rearrange_mode = 'preserve_contig'
            return
        else:
            if all_splittable:  # no variant case
                # split series range into evenly-spaced subsets
                start0_list = np.linspace(
                    series_start0,
                    series_end0,
                    np.ceil(series_width / max_pileup_width).astype('int') + 1,
                    endpoint=True,
                ).astype('int')
                for pileup_ser in self.pileupseries_dict.values():
                    pileup_ser.rearrange(start0_list, prepare_vcfspecs=True)

                self.no_variant = True
                self.rearrange_mode = None
                return

        # splitting contig vcfspec
        start0_list, all_splittable = self.get_rearrangement_points(
            split_region_gr=functools.reduce(
                lambda x, y: x.intersect(y), 
                (pileup_ser.get_splittable_region_split_contig() for pileup_ser in self.pileupseries_dict.values())
            ), 
            trim_margins=True, 
            max_pileup_width=max_pileup_width, 
            series_start0=series_start0, 
            series_end0=series_end0,
        )

        if start0_list is not None:
            for pileup_ser in self.pileupseries_dict.values():
                pileup_ser.rearrange(start0_list, prepare_vcfspecs=True)
            self.rearrange_mode = 'split_contig'
            return
        else:
            if all_splittable:
                raise Exception(f'all-splittable is True in split-contig after False in non-split-contig')

        # using inactive runs
        start0_list, all_splittable = self.get_rearrangement_points(
            split_region_gr=functools.reduce(
                lambda x, y: x.intersect(y), 
                (pileup_ser.get_splittable_region_inactive_runs() for pileup_ser in self.pileupseries_dict.values())
            ), 
            trim_margins=False, 
            max_pileup_width=max_pileup_width, 
            series_start0=series_start0, 
            series_end0=series_end0,
        )

        if start0_list is not None:
            for pileup_ser in self.pileupseries_dict.values():
                pileup_ser.rearrange(start0_list, prepare_vcfspecs=True)
            self.rearrange_mode = 'inactive_runs'
            return
        else:
            if all_splittable:
                raise Exception(f'All PileupSeries region is inactive.')
            else:
                self.rearrange_mode = 'failed'
                return

    @staticmethod
    def get_rearrangement_points(split_region_gr, trim_margins, max_pileup_width, series_start0, series_end0):
        """Returns:
            start0_list: list of 0-based positions. 
                With pairwise iteration, for each iteratation result (x, y),
                x is used as start0 and y as end0 for a RealignerPileup object.
        """
        if split_region_gr.empty:
            # abort this split strategy
            start0_list = None
            all_splittable = False
            return start0_list, all_splittable

        # when entire range is splittable
        if split_region_gr.df.shape[0] == 1:
            row = split_region_gr.df.iloc[0, :]
            if (row.Start == series_start0) and (row.End == series_end0):
                start0_list = None
                all_splittable = True
                return start0_list, all_splittable

        # make into ranges
        start_candidate_ranges = [range(row.Start, row.End + 1) for (idx, row) in split_region_gr.df.iterrows()]

        # trim or add margins
        if start_candidate_ranges[0].start == series_start0:
            if trim_margins:
                start_candidate_ranges[0] = range(start_candidate_ranges[0].stop - 1, start_candidate_ranges[0].stop)
        else:
            start_candidate_ranges.insert(0, range(series_start0, series_start0 + 1))

        if start_candidate_ranges[-1].stop - 1 == series_end0:
            if trim_margins:
                start_candidate_ranges[-1] = range(start_candidate_ranges[-1].start, start_candidate_ranges[-1].start + 1)
        else:
            start_candidate_ranges.append(range(series_end0, series_end0 + 1))

        # now len(start_candidate_ranges) is at least 2 

        # If candidate range is smaller than max_pileup_width after trimming
        if trim_margins:
            new_start0 = start_candidate_ranges[0].start
            new_end0 = start_candidate_ranges[-1].stop - 1
            if new_end0 - new_start0 <= max_pileup_width:
                start0_list = [new_start0, new_end0]
                all_splittable = False
                return start0_list, all_splittable

        # validity check
        if any(
            rng2.start - (rng1.stop - 1) > max_pileup_width
            for rng1, rng2 in common.pairwise(start_candidate_ranges)
        ):
            # abort this split strategy
            start0_list = None
            all_splittable = False
            return start0_list, all_splittable

        # search for actual split points
        iterator = iter(start_candidate_ranges)
        current_rng = next(iterator)
        next_rng = next(iterator)
        start0_list = [current_rng.start]
        stopiter = False

        while True:
            if stopiter:
                break

            last_start0 = start0_list[-1]
            while True:
                if next_rng.start <= last_start0 + max_pileup_width:
                    current_rng = next_rng

                    try:
                        next_rng = next(iterator)
                    except StopIteration:
                        stopiter = True
                        break
                else:
                    if current_rng.stop - 1 <= last_start0:
                        raise Exception(f'Cannot make pileup margins within max_pileup_width')
                    start0_list.append(
                        min(current_rng.stop - 1, last_start0 + max_pileup_width)
                    )
                    break

        start0_list.append(next_rng.start)

        all_splittable = False
        return start0_list, all_splittable

    def set_multisample_pileups(self):
        self.mspileup_list = list()
        num_col = len(next(self.iter_series()).pileup_list)
        for idx in range(num_col):
            pileup_dict = {
                sampleid: self.pileupseries_dict[sampleid].pileup_list[idx]
                for sampleid in self.pileupseries_dict.keys()
            }
            self.mspileup_list.append(
                realign_msrpileup.MultisampleRealignerPileup(
                    pileup_dict=pileup_dict, 
                    refver=self.refver,
                    fasta=self.fasta, 
                    aligner=self.aligner,
                    verbose=self.verbose, 
                    logger=self.logger,
                    **self.params,
                )
            )

    def set_realigned_reads(self):
        """Assumes all mspileups have executed 'save_subseq_alignments' and 'allocate_subseq_alignments'"""
        #for msrpileup in self.mspileup_list:
        #    msrpileup.save_subseq_alignments()
        #    msrpileup.allocate_subseq_alignments()
        self.realigned_reads = {sampleid: dict() for sampleid in self.pileupseries_dict.keys()}
        for sampleid, rpileup_ser in self.pileupseries_dict.items():
            rpileup_ser.set_realigned_reads()
            self.realigned_reads[sampleid] = rpileup_ser.realigned_reads

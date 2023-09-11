import re
import itertools
import collections
import functools
import inspect
import logging
import uuid

import numpy as np
import pandas as pd

import handygenome.tools as tools
import handygenome.refgenome.refgenome as refgenome
import handygenome.workflow as workflow
import handygenome.read.readhandler as readhandler
import handygenome.align.alignhandler as alignhandler
import handygenome.read.readplus as readplus


class PileupBase:
    pat_insclip = re.compile("(\(.+\))?([^()]+)(\(.+\))?")
    pat_parenthesis = re.compile("[()]")
    DEL_VALUE = "*"
    EMPTY_VALUE = "_"
    DF_FILL_VALUES = {
        'main': EMPTY_VALUE,
        'bq': -1,
        'queryonly': '',
        'queryonly_bq': np.nan,
        'onleft': False,
        'unedited_main': EMPTY_VALUE,
        'queryonly_added_main': EMPTY_VALUE,
    }
        
    ##########################
    # methods for direct use #
    ##########################
    @property
    def start0(self):
        return self.seq_df.columns[0]

    @property
    def end0(self):
        return self.seq_df.columns[-1] + 1

    @property
    def range0(self):
        return range(self.start0, self.end0)

    @property
    def width(self):
        return self.seq_df.shape[1]

    def get_ref_seq(self):
        return self.fasta.fetch(self.chrom, self.start0, self.end0)

    def get_depth(self, pos0):
        col = self.seq_df.loc[:, pos0]
        return (col != self.__class__.EMPTY_VALUE).sum()
        #return sum(x != self.__class__.EMPTY_VALUE for x in col)

    def get_queryonly_added_df(self):
        main_asnp = self.dfs['main'].to_numpy().astype(str)
        queryonly_asnp = self.dfs['queryonly'].to_numpy().astype(str)

        queryonly_exists = (queryonly_asnp != self.__class__.DF_FILL_VALUES['queryonly'])
        main_replaced = np.where(
            (main_asnp == self.__class__.DF_FILL_VALUES['main']),
            '',
            main_asnp,
        )
        result = np.where(
            queryonly_exists,
            np.where(
                self.dfs['onleft'], 
                np.char.add(queryonly_asnp, main_asnp),
                np.char.add(main_asnp, queryonly_asnp),
            ),
            main_asnp,
        )
        result = pd.DataFrame(
            result, 
            index=self.dfs['main'].index, 
            columns=self.dfs['main'].columns, 
        )
        return result

    def get_allele_counter(self, pos0):
        col = self.dfs['queryonly_added_main'].loc[:, pos0]
        return collections.Counter(col[col != self.__class__.EMPTY_VALUE])
        #return collections.Counter(x for x in col if x != self.__class__.EMPTY_VALUE)

    def get_allele_portions(self, pos0):
        counter = self.get_allele_counter(pos0)
        counter_sum = sum(counter.values())
        portions = dict()
        for key, val in counter.items():
            portions[key] = val / counter_sum
        return portions

    def get_read(self, read_uid):
        return self.read_store[read_uid]

    def get_lowbq_replaced_df(self, main_df, bq_df, cutoff=15):
        ref_seq = self.get_ref_seq()

        bq_exists = (bq_df != self.__class__.DF_FILL_VALUES['bq'])
        bq_is_low = np.logical_and(
            bq_exists,
            bq_df < cutoff,
        )
        bq_is_low = bq_is_low.to_numpy()

        result = main_df.copy()
        for row_idx, col_idx in zip(*bq_is_low.nonzero()):
            result.iloc[row_idx, col_idx] = ref_seq[col_idx]

        return result

    @property
    def seq_df(self):
        try:
            return self.dfs['queryonly_added_main']
        except KeyError:
            return self.dfs['main']
            

    ############
    # backends #
    ############
    def __init__(
        self, 
        bam, chrom, 
        start0=None, end0=None, 
        refver=None, fasta=None,

        init_df=True, 

        lowbq_cutoff=15,

        verbose=False, 
        logger=None,
    ):
        """Args:
            start0 and end0 are only required when init_df == True
        """
        # set params
        self.bam = bam
        self.chrom = chrom
        self.refver, self.fasta = self.refver_fasta_arghandler(refver, fasta)

        # set logger
        self.verbose, self.logger = self.logger_arghandler(verbose, logger)

        # init df and subsequent ones
        if init_df:
            (
                df, bq_df, queryonly_df, queryonly_bq_df, onleft_df, read_store,
            )= make_pileup_components(
                chrom,
                start0,
                end0,
                bam,
                truncate=True,
                as_array=False,
                verbose=self.verbose,
                logger=self.logger,
            )
            self.dfs = {
                'main': df,
                'bq': bq_df,
                'queryonly': queryonly_df,
                'queryonly_bq': queryonly_bq_df,
                'onleft': onleft_df,
            }
            self.read_store = read_store

            # replace lowbq with REF
            lowbq_replaced_df = self.get_lowbq_replaced_df(
                main_df=self.dfs['main'],
                bq_df=self.dfs['bq'],
                cutoff=lowbq_cutoff,
            )
            self.dfs['unedited_main'] = self.dfs['main']
            self.dfs['main'] = lowbq_replaced_df

            # make queryonly-added main
            self.dfs['queryonly_added_main'] = self.get_queryonly_added_df()

            self._set_MQ()
        else:
            self.dfs = None
            self.read_store = None
            self.MQ = None

    def __repr__(self):
        return (
            f'<'
            f'{self.__class__.__name__}'
            f'(chrom={self.chrom}, start0={self.start0:,}, end0={self.end0:,}, shape={self.seq_df.shape})'
            f'>'
        )

    def _coord_arg_sanitycheck(self, start0, end0):
        #if start0 is not None:
        if start0 < self.start0:
            raise Exception('Input "start0" argument is out of pileup range.')
        #if end0 is not None:
        if end0 > self.end0:
            raise Exception('Input "end0" argument is out of pileup range.')

    def _coord_arg_sanitycheck_pos0(self, pos0):
        if pos0 not in self.seq_df.columns:
            raise Exception('Input "pos0" argument is out of pileup range.')

    def _set_MQ(self):
        if self.seq_df.shape[1] == 0:
            self.MQ = pd.Series([], dtype=float)
        else:
            if self.seq_df.shape[0] > 0:
                MQ_data = {pos0: list() for pos0 in self.seq_df.columns}
                #for read_uid, row in self.df.iterrows():
                for read_uid in self.seq_df.index:
                    read = self.read_store[read_uid]
                    start0 = max(self.start0, read.reference_start)
                    end0 = min(self.end0, read.reference_end)
                    for pos0 in range(start0, end0):
                        MQ_data[pos0].append(read.mapping_quality)
                self.MQ = pd.Series(
                    [
                        np.nan if (len(MQ_data[pos0]) == 0) else np.mean(MQ_data[pos0])
                        for pos0 in self.seq_df.columns
                    ],
                    index=self.seq_df.columns,
                )
            else:
                self.MQ = pd.Series(
                    [np.inf] * self.seq_df.shape[1],
                    index=self.seq_df.columns,
                )

    def _spawn_base(self, attrs_to_copy, init_df=False):
        """Make a partial copy of self with essential attributes"""
        #for key in inspect.signature(self.__init__).parameters.keys():
        #    if key not in ('init_df', 'start0', 'end0'):
        #        kwargs[key] = getattr(self, key)

        kwargs = dict()
        for key in attrs_to_copy: 
            kwargs[key] = getattr(self, key)
        kwargs['init_df'] = init_df

        result = self.__class__(**kwargs)
        result.read_store = self.read_store
        return result

    def _split_base(self, start0_list):
        assert all((x in range(self.start0 + 1, self.end0)) for x in start0_list)
        result = list()
        margins = [self.start0] + start0_list + [self.end0]
        for new_start0, new_end0 in tools.pairwise(margins):
            result.append(self.subset(new_start0, new_end0, inplace=False))
        return result

    def _subset_base(self, start0, end0, inplace=False):
        #subset_df, subset_bq_df = self._get_subset_df(start0, end0)
        subset_dfs = self._get_subset_df(start0, end0)
        new_MQ = self.MQ.loc[start0:(end0 - 1)]
        if inplace:
            self.dfs = subset_dfs
            self.MQ = new_MQ
            return None
        else:
            result = self.spawn()
            result.dfs = subset_dfs
            result.MQ = new_MQ
            return result

    def _get_subset_df(self, start0, end0):
        self._coord_arg_sanitycheck(start0, end0)
        # get row selector
        row_within_range = list()
        for row_id in self.seq_df.index:
            read = self.read_store[row_id]
            if read.reference_end <= start0 or read.reference_start >= end0:
                row_within_range.append(False)
            else:
                row_within_range.append(True)

        # subset dataframe
        new_dfs = {
            key: val.loc[row_within_range, start0:(end0 - 1)]
            for key, val in self.dfs.items()
        }
        #new_df = self.df.loc[row_within_range, start0:(end0 - 1)]
        #new_bq_df = self.bq_df.loc[row_within_range, start0:(end0 - 1)]

        return new_dfs

    def _get_subset_df_old(self, start0, end0):
        self._coord_arg_sanitycheck(start0, end0)
        # subset dataframe
        new_df = self.df.loc[:, start0:(end0 - 1)]
        new_bq_df = self.bq_df.loc[:, start0:(end0 - 1)]
        # remove out-of-range rows
        row_within_range = list()
        #for row_id, row in new_df.iterrows():
        for row_id in new_df.index:
            read = self.read_store[row_id]
            if read.reference_end <= start0 or read.reference_start >= end0:
                row_within_range.append(False)
            else:
                row_within_range.append(True)
        new_df = new_df.loc[row_within_range, :]
        new_bq_df = new_bq_df.loc[row_within_range, :]

        return new_df, new_bq_df

    def _merge_base(self, other, other_on_left):
        if other_on_left:
            left = other
            right = self
        else:
            left = self
            right = other

        assert left.end0 == right.start0
        border_col_idx_left = left.seq_df.shape[1] - 1
        new_dfs = {
            key: left.dfs[key].join(right.dfs[key], how="outer")
            for key in left.dfs.keys()
        }
        #for key, df in zip(keys, vals):
        for key in list(new_dfs.keys()):
            df = new_dfs[key]
            try:
                df.where(df.notna(), other=self.__class__.DF_FILL_VALUES[key], inplace=True)
            except TypeError as exc:
                if str(exc) == 'Cannot do inplace boolean setting on mixed-types with a non np.nan value':
                    new_dfs[key] = df.where(
                        df.notna(), other=self.__class__.DF_FILL_VALUES[key], inplace=False,
                    )
                else:
                    raise
                    
        self.dfs = new_dfs
        #self.df = left.df.join(right.df, how="outer")
        #self.bq_df = left.bq_df.join(right.bq_df, how="outer")
        self.MQ = pd.concat([left.MQ, right.MQ])

#        try:
#            self.df[self.df.isnull()] = self.__class__.EMPTY_VALUE  # turn NaN into EMPTY_VALUE
#        except:
#            print(self.df.shape)
#            print(self.df)
#            print(type(self.df.iloc[0, 0]))
#            raise
        #self._handle_facing_insclips(border_col_idx_left)
        self.read_store.update(other.read_store)

    def _handle_facing_insclips(self, border_col_idx_left):
        border_col_idx_right = border_col_idx_left + 1
        border_col_left = self.df.iloc[:, border_col_idx_left]
        border_col_right = self.df.iloc[:, border_col_idx_right]
        for row_idx, (val_left, val_right) in enumerate(
            zip(border_col_left, border_col_right)
        ):
            if val_right.startswith("(") and val_left.endswith(")"):
                mat_left = self.__class__.pat_insclip.fullmatch(val_left)
                mat_right = self.__class__.pat_insclip.fullmatch(val_right)
                if mat_left.group(3) != mat_right.group(1):
                    raise Exception(
                        f"Insclip seqs of adjacent entries are different.\n{self.df}"
                    )
                self.df.iloc[
                    row_idx, border_col_idx_right
                ] = self.__class__.pat_insclip.sub("\\2\\3", val_right)

    def _extend_base(self, width, left):
        new_pileup = self._make_extend_pileup(width, left)
        self.merge(new_pileup, other_on_left=left)

    @staticmethod
    def refver_fasta_arghandler(refver, fasta):
        if (refver is None) and (fasta is None):
            raise Exception(f'At least one of "refver" or "fasta" must be set.')

        if refver is None:
            refver_result = refgenome.infer_refver_fasta(fasta)
        else:
            refver_result = refver

        if fasta is None:
            fasta_result = refgenome.get_fasta(refver_result)
        else:
            fasta_result = fasta

        return refver_result, fasta_result

    @classmethod
    def logger_arghandler(cls, verbose, logger):
        verbose_result = verbose
        if logger is None:
            logger_result = cls.get_logger(verbose_result)
        else:
            logger_result = logger

        return verbose_result, logger_result

    @staticmethod
    def get_logger(verbose):
        formatter = logging.Formatter(
            fmt='[%(asctime)s.%(msecs)03d] Pileup: %(message)s', 
            datefmt='%Z %Y-%m-%d %H:%M:%S'
        )
        return workflow.get_logger(
            name=str(uuid.uuid4()),
            level=('debug' if verbose else 'info'),
            formatter=formatter,
        )


class Pileup(PileupBase):
    # initializers
    def __init__(
        self, 
        bam, chrom, 
        start0=None, end0=None, 
        refver=None, fasta=None,
        init_df=True, 
        verbose=False, logger=None,
    ):
        PileupBase.__init__(
            self, 
            bam=bam, chrom=chrom, 
            start0=start0, end0=end0, 
            refver=refver, fasta=fasta,
            init_df=init_df, 
            verbose=verbose, logger=logger,
        )

    def subset(self, start0, end0, inplace=False):
        return self._subset_base(self, start0, end0, inplace=inplace)

    def merge(self, other, other_on_left):
        self._merge_base(other, other_on_left=other_on_left)

    def split(self, start0_list):
        return self._split_base(start0_list)

    def spawn(self):
        return self._spawn_base(
            ('bam', 'chrom', 'refver', 'fasta', 'verbose', 'logger'),
            init_df=False,
        )

    # extend
    def extend_left(self, width):
        self._extend_base(width, left=True)

    def extend_right(self, width):
        self._extend_base(width, left=False)

    def _make_extend_pileup(self, width, left):
        if left:
            start0 = self.start0 - width
            end0 = self.start0
        else:
            start0 = self.end0
            end0 = self.end0 + width
        return self.__class__(
            fasta=self.fasta, 
            bam=self.bam, 
            chrom=self.chrom, 
            start0=start0,
            end0=end0,
            init_df=True,
        )


class MultisamplePileup(PileupBase):
    def __init__(self, pileup_dict):
        """Args:
            pileup_dict: keys - sampleid; values - Pileup object
        """
        self.pileup_dict = pileup_dict

        self.first_pileup = next(iter(self.pileup_dict.values()))
        self.chrom = self.first_pileup.chrom
        self.fasta = self.first_pileup.fasta

    def _get_active_threshold(self):
        depth_sum = sum(x.df.shape[0] for x in self.pileup_dict.values())
        corrected_active_thresholds = (
            x.active_threshold * (x.df.shape[0] / depth_sum)
            for x in self.pileup_dict.values()
        )
        return min(corrected_active_thresholds)

    def _get_active_threshold_simple(self):
        return sum(x.active_threshold for x in self.pileup_dict.values()) / (len(self.pileup_dict) ** 2)

    @property
    def start0(self):
        return self.first_pileup.start0

    @property
    def end0(self):
        return self.first_pileup.end0

    def set_df(self):
        assert len(set(x.range0 for x in self.pileup_dict.values())) == 1, f'Genomic ranges of pileup objects are different.'
        self.df = pd.concat(
            {key: val.df for key, val in self.pileup_dict.items()},
            names=['SampleID', 'ReadUID'],
        )
        #self._set_active_info()

    def get_read(self, row_id):
        sampleid, read_uid = row_id
        return self.pileup_dict[sampleid].read_store[read_uid]

    ### extend ###
    def extend_rightward(self, width):
        for pileup in self.pileup_dict.values():
            pileup.extend_rightward(width)

    def extend_leftward(self, width):
        for pileup in self.pileup_dict.values():
            pileup.extend_leftward(width)

    # equalize
    def equalize_margins_leftward(self):
        max_range_start0 = min(x.start0 for x in self.pileup_dict.values())
        for pileup in self.pileup_dict.values():
            if pileup.start0 > max_range_start0:
                pileup.extend_leftward(pileup.start0 - max_range_start0)

    def equalize_margins_rightward(self):
        max_range_end0 = max(x.end0 for x in self.pileup_dict.values())
        for pileup in self.pileup_dict.values():
            if pileup.end0 < max_range_end0:
                pileup.extend_rightward(max_range_end0 - pileup.end0)

    def equalize_margins(self):
        self.equalize_margins_leftward()
        self.equalize_margins_rightward()

    ### row specs ###
#    def set_row_specs(self, set_for_each=False):
#        if set_for_each:
#            for pileup in self.pileup_dict.values():
#                pileup.set_row_specs()
#
#        self.row_specs = dict()
#        for row_id, row in self.df.iterrows():
#            sampleid, read_uid = row_id
#            self.row_specs[row_id] = self.pileup_dict[sampleid].row_specs[read_uid]

    def set_row_spec_groups(self):
        """Assumes row_spec_groups is set for each single sample Pileup"""
        # make superseq_row_specs
        candidate_superseq_row_specs = dict()
        for sampleid, pileup in self.pileup_dict.items():
            for read_uid in pileup.row_spec_groups.keys():
                candidate_superseq_row_specs[(sampleid, read_uid)] = pileup.row_specs[read_uid]
        superseq_row_specs = dict()
        for superseq_rowid in self.__class__.group_row_specs(candidate_superseq_row_specs).keys():
            superseq_row_specs[superseq_rowid] = candidate_superseq_row_specs[superseq_rowid]
        # make row_spec_groups for each sample
        row_spec_groups_by_sample = dict()
        for sampleid, pileup in self.pileup_dict.items():
            row_spec_groups = collections.OrderedDict()
            for superseq_rowid in superseq_row_specs.keys():
                row_spec_groups[superseq_rowid] = {
                    'subseq_rowids': set(),
                    'subseq_hits': 0,
                }
            for query_read_uid, query_row_spec in pileup.row_specs.items():
                superseq_candidates = list()
                for superseq_rowid, superseq_row_spec in superseq_row_specs.items():
                    if self.__class__.row_spec_matcher(query_row_spec, superseq_row_spec):
                        superseq_candidates.append(superseq_rowid)
                    #if len(superseq_candidates) >= 2:
                    #    break

                if len(superseq_candidates) == 0:
                    raise Exception(f'row_spec in sub-pileup does not fit to any of superseqs:\nsampleid: {sampleid}\nrow_spec of subseq: {query_row_spec}')
                else:
                    self.__class__.handle_matching_subseq(superseq_candidates, row_spec_groups, query_read_uid)

            row_spec_groups_by_sample[sampleid] = row_spec_groups
        # final    
        self.superseq_row_specs = superseq_row_specs
        self.row_spec_groups = row_spec_groups_by_sample

    def save_superseq_alignments(self):
        librealign.save_superseq_alignments_multisample(self)

    # get realigned reads #
    def get_realigned_reads(self):
        """Must be run by those created from 'get_active_region_pileup_multisample' function

        Returns:
            dict (keys sampleid, values dict (keys ReadUID, values pysam.AlignedSegment))
        """
        return librealign.get_realigned_reads_multisample(self)

    def set_vcfspecs(self, allele_portion_threshold=None, concat_dist_le=None):
        """Must be run by those created from 'get_active_region_pileup_multisample' function"""
        self.vcfspecs = librealign.get_vcfspecs_from_pileup_multisample(self, allele_portion_threshold=allele_portion_threshold, concat_dist_le=concat_dist_le)

    # for debugging
    def show_row_spec_alignments(self, skip_zero_hits=True, show_subseqs=True):
        for sampleid, pileup in self.pileup_dict.items():
            print('@@@', sampleid, '@@@')
            print()
            self._show_row_spec_alignments_helper(self.row_spec_groups[sampleid], self.superseq_row_specs, skip_zero_hits, show_subseqs)


def make_pileup_components(
    chrom,
    start0,
    end0,
    bam,
    truncate=True,
    as_array=False,
    del_value=PileupBase.DEL_VALUE,
    empty_value=PileupBase.EMPTY_VALUE,
    verbose=False,
    logger=None,
    trailing_queryonly_to_right=False,
    preserve_trailing_queryonly=True,
):
    """Args:
        trailing_queryonly_to_right: 
            - If True, insertions or softclips on the right border of the read
                are assigned to "read.reference_end"
            - If False, assigned to "read.reference_end - 1"
        preserve_trailing_queryonly:
            - Relevant only when "truncate" is True and "trailing_queryonly_to_right" is True.
            - If True, if trailing queryonly is placed on "end0" position,
                returned array or dataframe gets to include the position "end0"
                to accommodate the trailing queryonly sequence. In this case,
                returned "ref_span_range" becomes range(start0, end0 + 1).
            - If False, if trailing queryonly is placed on "end0" position,
                they are not included in the returned array or dataframe,
                and "ref_span_range" becomes range(start0, end0).
    """
    # set logger
    if logger is None:
        logger = PileupBase.get_logger(verbose)
    else:
        logger = logger

    def make_read_store(bam, chrom, start0, end0):
        target_range0 = range(start0, end0)
        readlist = list()
        uid_list = list()
        start_list = list()
        end_list = list()
        for read in readhandler.get_fetch(bam, chrom, start0, end0, readfilter=readhandler.readfilter_pileup):
            if readhandler.check_cigarN_includes_range(read, start0, end0):
                continue

#            target_seq = readplus.ReadPlus(read, minimal=True).get_seq_from_range0(
#                target_range0,
#                flanking_queryonly_default_mode=True,
#            )
#            if target_seq == '':
#                print(read)
#                continue

            readlist.append(read)
            uid_list.append(readhandler.get_uid(read))
            start_list.append(read.reference_start)
            end_list.append(read.reference_end)

        if len(start_list) == 0:
            tmp_pileup_range = None
        else:
            tmp_pileup_range = range(
                min(min(start_list), start0), 
                max(max(end_list), end0) + 1,
            )

        read_store = collections.OrderedDict(zip(uid_list, readlist))

        return read_store, tmp_pileup_range

    def raise_cigarpattern_error(read):
        raise Exception(f"Unexpected cigar pattern:\n{read.to_string()}")

    def cigar_sanitycheck(read):
        """Assumes M.* or [IS]M.*"""
        error = False
        first_cigarop = read.cigartuples[0][0]
        if first_cigarop != 0:
            if first_cigarop not in (1, 4):
                error = True
            else:
                if read.cigartuples[1][0] != 0:
                    error = True
        if error:
            raise_cigarpattern_error(read)

    def truncate_array(arr, tmp_pileup_range, ref_span_range):
        sl_start = tmp_pileup_range.index(ref_span_range.start)
        if ref_span_range.stop == tmp_pileup_range.stop:
            sl_end = None
        else:
            sl_end = tmp_pileup_range.index(ref_span_range.stop)
        return arr[:, slice(sl_start, sl_end)]

    def make_array(
        start0, end0, tmp_pileup_range, read_store, empty_value, del_value, trailing_queryonly_to_right,
    ):
        if tmp_pileup_range is None:  # zero depth
            ref_span_range = range(start0, end0)
            arr = np.empty((0, len(ref_span_range)))
            bq_arr = arr.copy()
            queryonly_arr = arr.copy()
            queryonly_bq_arr = arr.copy()
            onleft_arr = arr.copy()
        else:
            # initialize array
            arr = np.full(
                shape=(len(read_store), len(tmp_pileup_range)),
                fill_value=PileupBase.DF_FILL_VALUES['main'],
                dtype=str,
            )
            bq_arr = np.full(
                shape=(len(read_store), len(tmp_pileup_range)),
                fill_value=PileupBase.DF_FILL_VALUES['bq'],
                dtype=int,
            )
            queryonly_arr = np.full(
                shape=(len(read_store), len(tmp_pileup_range)),
                fill_value=PileupBase.DF_FILL_VALUES['queryonly'],
                dtype=object,
            )
            queryonly_bq_arr = np.full(
                shape=(len(read_store), len(tmp_pileup_range)),
                fill_value=PileupBase.DF_FILL_VALUES['queryonly_bq'],
                dtype=object,
            )
            onleft_arr = np.full(
                shape=(len(read_store), len(tmp_pileup_range)),
                fill_value=PileupBase.DF_FILL_VALUES['onleft'],
                dtype=bool,
            )

            for arr_row_idx, read in enumerate(read_store.values()):
                arr_col_idx = read.reference_start - tmp_pileup_range.start
                try:
                    _write_read_to_array_row(
                        read, arr, arr_row_idx, arr_col_idx, del_value, 
                        bq_arr, queryonly_arr, queryonly_bq_arr, onleft_arr, 
                        trailing_queryonly_to_right,
                    )
                except Exception as exc:
                    try:
                        cigar_sanitycheck(read)
                    except Exception as exc_cigarpattern:
                        raise exc_cigarpattern from exc
                    else:
                        raise exc

            last_col_is_empty = (set(arr[:, -1]) == {empty_value})
            if last_col_is_empty:
                arr = arr[:, :-1]
                bq_arr = bq_arr[:, :-1]
                queryonly_arr = queryonly_arr[:, :-1]
                queryonly_bq_arr = queryonly_bq_arr[:, :-1]
                onleft_arr = onleft_arr[:, :-1]
                tmp_pileup_range = range(tmp_pileup_range.start, tmp_pileup_range.stop - 1)

            # truncate array
            if truncate:
                ref_span_range = range(start0, end0)
                arr = truncate_array(arr, tmp_pileup_range, ref_span_range)
                bq_arr = truncate_array(bq_arr, tmp_pileup_range, ref_span_range)
                queryonly_arr = truncate_array(queryonly_arr, tmp_pileup_range, ref_span_range)
                queryonly_bq_arr = truncate_array(queryonly_bq_arr, tmp_pileup_range, ref_span_range)
                onleft_arr = truncate_array(onleft_arr, tmp_pileup_range, ref_span_range)
            else:
                ref_span_range = tmp_pileup_range

        return arr, bq_arr, queryonly_arr, queryonly_bq_arr, onleft_arr, ref_span_range

    # main
    read_store, tmp_pileup_range = make_read_store(bam, chrom, start0, end0)
        # tmp_pileup_range contains end0 of the rightmost fetched read

    arr, bq_arr, queryonly_arr, queryonly_bq_arr, onleft_arr, ref_span_range = make_array(
        start0, end0, tmp_pileup_range, read_store, empty_value, del_value, trailing_queryonly_to_right,
    )

    # final result
    if as_array:
        return arr, bq_arr, queryonly_arr, queryonly_bq_arr, onleft_arr, ref_span_range, read_store
    else:
        df = pd.DataFrame(
            arr, 
            columns=list(ref_span_range), 
            index=pd.Index(read_store.keys(), tupleize_cols=False),
        )
        bq_df = pd.DataFrame(
            bq_arr, 
            columns=list(ref_span_range), 
            index=pd.Index(read_store.keys(), tupleize_cols=False),
        )
        queryonly_df = pd.DataFrame(
            queryonly_arr, 
            columns=list(ref_span_range), 
            index=pd.Index(read_store.keys(), tupleize_cols=False),
        )
        queryonly_bq_df = pd.DataFrame(
            queryonly_bq_arr, 
            columns=list(ref_span_range), 
            index=pd.Index(read_store.keys(), tupleize_cols=False),
        )
        onleft_df = pd.DataFrame(
            onleft_arr, 
            columns=list(ref_span_range), 
            index=pd.Index(read_store.keys(), tupleize_cols=False),
        )
        return df, bq_df, queryonly_df, queryonly_bq_df, onleft_df, read_store


def _write_read_to_array_row(
    read, arr, arr_row_idx, arr_col_idx, del_value, 
    bq_arr, queryonly_arr, queryonly_bq_arr, onleft_arr, 
    trailing_queryonly_to_right=True,
    bq_array_includes_queryonly=False,
):
    """Helper function for 'make_pileup_components'"""
    query_idx = 0
    queryonly_seq = None
    for idx, (key, subiter) in enumerate(
        itertools.groupby(
            read.cigartuples, 
            key=(lambda x: alignhandler.CIGAR_WALK_DICT[x[0]]),
        )
    ):
        consume_target, consume_query = key
        cigarlen_sum = sum(x[1] for x in subiter)

        if consume_target and consume_query:  # match
            if queryonly_seq is None:
                sl_query = slice(query_idx, query_idx + cigarlen_sum)
                sl_arr_col = slice(arr_col_idx, arr_col_idx + cigarlen_sum)
                arr[arr_row_idx, sl_arr_col] = tuple(read.query_sequence[sl_query])

                bq_arr[arr_row_idx, sl_arr_col] = tuple(read.query_qualities[sl_query])  # BQ array
            else:
                # add query-only information
                arr[arr_row_idx, arr_col_idx] = read.query_sequence[query_idx]
                bq_arr[arr_row_idx, arr_col_idx] = read.query_qualities[query_idx]
                queryonly_arr[arr_row_idx, arr_col_idx] = queryonly_seq
                queryonly_bq_arr[arr_row_idx, arr_col_idx] = queryonly_bq
                onleft_arr[arr_row_idx, arr_col_idx] = True

                # add trailing matched bases
                sl_query = slice(query_idx + 1, query_idx + cigarlen_sum)
                sl_arr_col = slice(arr_col_idx + 1, arr_col_idx + cigarlen_sum)
                arr[arr_row_idx, sl_arr_col] = tuple(read.query_sequence[sl_query])
                bq_arr[arr_row_idx, sl_arr_col] = tuple(read.query_qualities[sl_query])  # BQ array

                queryonly_seq = None
                queryonly_bq = None

            arr_col_idx += cigarlen_sum
            query_idx += cigarlen_sum

        elif (not consume_target) and consume_query:  # ins, softclip
            queryonly_seq = read.query_sequence[query_idx:(query_idx + cigarlen_sum)]
            queryonly_bq = tuple(read.query_qualities[query_idx:(query_idx + cigarlen_sum)])  # BQ array
            query_idx += cigarlen_sum

        elif consume_target and (not consume_query):  # del, skip
            if queryonly_seq is not None:
                raise Exception(f'Cigar D or N comes right after I or S. Current read: {read.to_string()}')

            sl_arr_col = slice(arr_col_idx, arr_col_idx + cigarlen_sum)
            arr[arr_row_idx, sl_arr_col] = del_value

            arr_col_idx += cigarlen_sum

    # When the last cigar operation is query-only
    if queryonly_seq is not None:
        if idx == 0:  # query-only-only case
            queryonly_arr[arr_row_idx, arr_col_idx] = queryonly_seq
            queryonly_bq_arr[arr_row_idx, arr_col_idx] = queryonly_bq  # BQ array
        else:  # queryonly_seq is appended to the base on the left
            queryonly_arr[arr_row_idx, arr_col_idx - 1] = queryonly_seq
            onleft_arr[arr_row_idx, arr_col_idx - 1] = False
            queryonly_bq_arr[arr_row_idx, arr_col_idx - 1] = queryonly_bq


def _write_read_to_array_row_old(
    read, arr, arr_row_idx, arr_col_idx, del_value, 
    bq_arr, queryonly_arr, onleft_arr, 
    trailing_queryonly_to_right=True,
    bq_array_includes_queryonly=False,
):
    """Helper function for 'make_pileup_components'"""
    query_idx = 0
    queryonly_seq = None
    for idx, (key, subiter) in enumerate(
        itertools.groupby(
            read.cigartuples, 
            key=(lambda x: alignhandler.CIGAR_WALK_DICT[x[0]]),
        )
    ):
        consume_target, consume_query = key
        cigarlen_sum = sum(x[1] for x in subiter)

        if consume_target and consume_query:  # match
            if queryonly_seq is None:
                sl_query = slice(query_idx, query_idx + cigarlen_sum)
                sl_arr_col = slice(arr_col_idx, arr_col_idx + cigarlen_sum)
                arr[arr_row_idx, sl_arr_col] = tuple(read.query_sequence[sl_query])

                bq_arr[arr_row_idx, sl_arr_col] = tuple(read.query_qualities[sl_query])  # BQ array
            else:
                # add the first match-base, prefixed with query-only sequence buffer
                arr[arr_row_idx, arr_col_idx] = f'({queryonly_seq})' + read.query_sequence[query_idx]
                if bq_array_includes_queryonly:  # BQ array
                    bq_arr[arr_row_idx, arr_col_idx] = (
                        queryonly_bq, 
                        read.query_qualities[query_idx],
                    )  
                else:
                    bq_arr[arr_row_idx, arr_col_idx] = read.query_qualities[query_idx]

                # add trailing matched bases
                sl_query = slice(query_idx + 1, query_idx + cigarlen_sum)
                sl_arr_col = slice(arr_col_idx + 1, arr_col_idx + cigarlen_sum)
                arr[arr_row_idx, sl_arr_col] = tuple(read.query_sequence[sl_query])
                bq_arr[arr_row_idx, sl_arr_col] = tuple(read.query_qualities[sl_query])  # BQ array

                queryonly_seq = None
                queryonly_bq = None

            arr_col_idx += cigarlen_sum
            query_idx += cigarlen_sum

        elif (not consume_target) and consume_query:  # ins, softclip
            queryonly_seq = read.query_sequence[query_idx:(query_idx + cigarlen_sum)]
            queryonly_bq = tuple(read.query_qualities[query_idx:(query_idx + cigarlen_sum)])  # BQ array
            query_idx += cigarlen_sum

        elif consume_target and (not consume_query):  # del, skip
            if queryonly_seq is not None:
                raise Exception(f'Cigar D or N comes right after I or S. Current read: {read.to_string()}')

            sl_arr_col = slice(arr_col_idx, arr_col_idx + cigarlen_sum)
            arr[arr_row_idx, sl_arr_col] = del_value
            bq_arr[arr_row_idx, sl_arr_col] = del_value  # BQ array

            arr_col_idx += cigarlen_sum

    # When the last cigar operation is query-only
    if queryonly_seq is not None:
        if idx == 0:  # query-only-only case
            arr[arr_row_idx, arr_col_idx] = f'({queryonly_seq})'
            if bq_array_includes_queryonly:  # BQ array
                bq_arr[arr_row_idx, arr_col_idx] = queryonly_bq  # BQ array
        else:  # queryonly_seq is appended to the base on the left
            if trailing_queryonly_to_right:
                arr[arr_row_idx, arr_col_idx] = f'({queryonly_seq})'
                if bq_array_includes_queryonly:  # BQ array
                    bq_arr[arr_row_idx, arr_col_idx] = queryonly_bq  # BQ array
            else:
                arr[arr_row_idx, arr_col_idx - 1] = arr[arr_row_idx, arr_col_idx - 1] + f'({queryonly_seq})'
                if bq_array_includes_queryonly:  # BQ array
                    bq_arr[arr_row_idx, arr_col_idx - 1] = (
                        bq_arr[arr_row_idx, arr_col_idx - 1],
                        queryonly_bq,
                    )  # BQ array


def _write_read_to_array_row_deprecated(read, arr, arr_row_idx, initial_arr_col_idx, del_value):
    """NOT USED. THIS VERSION ASSIGNS INSCLIP TO THE POSITION ON THE LEFT."""
    def handle_leading_insclip(read):
        leading_insclip_len = 0
        cigartup_idx = 0
        for cigarop, cigarlen in alignhandler.iter_leading_queryonly(read.cigartuples):
            leading_insclip_len += cigarlen
            cigartup_idx += 1

        if leading_insclip_len == 0:
            leading_insseq = None
        else:
            leading_insseq = read.query_sequence[:leading_insclip_len]
        query_idx = leading_insclip_len

        return leading_insseq, cigartup_idx, query_idx

    def handle_queryonly_only_case(read, arr, arr_col_idx, arr_row_idx):
        arr[arr_row_idx, arr_col_idx] = "(" + read.query_sequence + ")"            

    def handle_targetonly_after_match(
        read,
        arr,
        query_idx,
        arr_col_idx,
        arr_row_idx,
        trailing_nonM_cigarlen_sum,
        del_value,
    ):
        arr[arr_row_idx, arr_col_idx] = read.query_sequence[
            query_idx
        ]  # query_idx indicates the last base of M
        arr[
            arr_row_idx,
            (arr_col_idx + 1) : (arr_col_idx + 1 + trailing_nonM_cigarlen_sum),
        ] = del_value

    def handle_queryonly_after_match(
        read,
        arr,
        query_idx,
        arr_col_idx,
        arr_row_idx,
        trailing_nonM_cigarlen_sum,
    ):
        last_match_base = read.query_sequence[query_idx]
        insseq = read.query_sequence[
            (query_idx + 1) : (query_idx + 1 + trailing_nonM_cigarlen_sum)
        ]
        arr[arr_row_idx, arr_col_idx] = f"{last_match_base}({insseq})"

    def handle_nonlast_match(
        read,
        arr,
        cigarlen,
        query_idx,
        arr_col_idx,
        arr_row_idx,
        cigartup_idx,
    ):
        # adds all match seqs except the last one
        # query_idx and arr_col_idx updated
        if cigarlen > 1:
            sl_query = slice(query_idx, query_idx + cigarlen - 1)
            sl_arr_col = slice(arr_col_idx, arr_col_idx + cigarlen - 1)
            arr[arr_row_idx, sl_arr_col] = tuple(read.query_sequence[sl_query])
            query_idx += cigarlen - 1
            arr_col_idx += cigarlen - 1
        # get all subsequent non-M cigar units
        # cigartup_idx updated
        trailing_nonMs = list()
        while True:
            cigartup_idx += 1
            if cigartup_idx >= len(read.cigartuples):
                break
            else:
                next_cigartup = read.cigartuples[cigartup_idx]
                if next_cigartup[0] != 0:
                    trailing_nonMs.append(next_cigartup)
                else:
                    break
        # add trailing non-M cigar units
        # query_idx and arr_col_idx updated
        trailing_nonM_cigarops = set(x[0] for x in trailing_nonMs)
        trailing_nonM_cigarlen_sum = sum(x[1] for x in trailing_nonMs)
        if trailing_nonM_cigarops.issubset(alignhandler.CIGAROPS_TARGETONLY):  # D, N
            handle_targetonly_after_match(
                read,
                arr,
                query_idx,
                arr_col_idx,
                arr_row_idx,
                trailing_nonM_cigarlen_sum,
                del_value,
            )
            query_idx += 1
            arr_col_idx += trailing_nonM_cigarlen_sum + 1
        elif trailing_nonM_cigarops.issubset(alignhandler.CIGAROPS_QUERYONLY):  # I, S
            handle_queryonly_after_match(
                read,
                arr,
                query_idx,
                arr_col_idx,
                arr_row_idx,
                trailing_nonM_cigarlen_sum,
            )
            query_idx += trailing_nonM_cigarlen_sum + 1
            arr_col_idx += 1
        else:
            raise Exception(f"Consecutive Non-M cigar units are composed of both target-only and query-only ones:\n{read.to_string()}")

        return cigartup_idx, query_idx, arr_col_idx

    def handle_last_match(
        read,
        arr,
        cigarlen,
        query_idx,
        arr_col_idx,
        arr_row_idx,
    ):
        # adds all match seqs to array
        sl_query = slice(query_idx, query_idx + cigarlen)
        sl_arr_col = slice(arr_col_idx, arr_col_idx + cigarlen)
        arr[arr_row_idx, sl_arr_col] = tuple(read.query_sequence[sl_query])

    # main
    """1) First, leading query-only cigar units are collected.
    2) Then, cigartuple is iterated on. Initial cigartup_idx is at the 
        first cigar unit which is not a query-only one. This cigar unit
        is assumed to be M. (If not, raises an exception)
    3) What is done in a loop cycle:
        The cigar unit which cigartup_idx indicates (which should be M)
        and all subsequent non-M cigarop are treated at single loop cycle.
        The subsequent non-M cigar units are assumed to be composed of only
        query-only ones or target-only ones. (If not, raises an exception)
    """

    # set mutable arr_col_idx. initial one is kept separately.
    arr_col_idx = initial_arr_col_idx

    # handles leading query-only cigar units. cigartup_idx and query_idx are initialized
    leading_insseq, cigartup_idx, query_idx = handle_leading_insclip(read)

    if cigartup_idx == len(read.cigartuples):
        # cigartup is only composed of query-only operations (I or S)
        handle_queryonly_only_case(read, arr, arr_col_idx, arr_row_idx)
    else:
        if read.cigartuples[cigartup_idx][0] != 0:
            raise Exception(
                f"The first cigar operation after stripping leading I and S is not M:\n{read.to_string()}"
            )

        # begins loop
        while True:
            if cigartup_idx > len(read.cigartuples) - 1:
                raise Exception(f'"cigartup_idx" became greater than "len(read.cigartuples) - 1" while looping.')

            cigarop, cigarlen = read.cigartuples[cigartup_idx]
            if cigarop != 0:
                raise Exception(f'The cigar unit indicated by "cigarop_idx" at the beginning of a loop cycle is not M.')

            if cigartup_idx == len(read.cigartuples) - 1:  
                # current cigartup is the last one
                handle_last_match(
                    read,
                    arr,
                    cigarlen,
                    query_idx,
                    arr_col_idx,
                    arr_row_idx,
                )
                break
            else:
                cigartup_idx, query_idx, arr_col_idx = handle_nonlast_match(
                    read,
                    arr,
                    cigarlen,
                    query_idx,
                    arr_col_idx,
                    arr_row_idx,
                    cigartup_idx,
                )
                if cigartup_idx == len(read.cigartuples):
                    break

        # add leading query-only seq
        if leading_insseq is not None:
            arr[arr_row_idx, initial_arr_col_idx] = (
                f"({leading_insseq})" + arr[arr_row_idx, initial_arr_col_idx]
            )


#def get_pileup(
#    chrom,
#    start0,
#    end0,
#    bam,
#    fasta,
#    truncate=True,
#    del_value=DEL_VALUE,
#    empty_value=EMPTY_VALUE,
#):
#    df, read_store = make_pileup_components(
#        chrom,
#        start0,
#        end0,
#        bam,
#        truncate=truncate,
#        as_array=False,
#        del_value=del_value,
#        empty_value=empty_value,
#    )
#    pileup = Pileup(df, chrom, fasta, bam, read_store)
#
#    return pileup





#def get_pileup_multisample(
#    chrom,
#    start0,
#    end0,
#    bam_dict,
#    fasta,
#    active_threshold_onesample=None,
#    del_value=PileupBase.DEL_VALUE,
#    empty_value=PileupBase.EMPTY_VALUE,
#):
#    # parameter handling
#    if active_threshold_onesample is None:
#        active_threshold_onesample = librealign.DEFAULT_ACTIVE_THRESHOLD
#    # initialize pileup_dict
#    pileup_dict = dict()
#    for sampleid, bam in bam_dict.items():
#        pileup_dict[sampleid] = get_pileup(
#            chrom,
#            start0,
#            end0,
#            bam=bam,
#            fasta=fasta,
#            active_threshold=active_threshold_onesample,
#            truncate=True,
#            as_array=False,
#            return_range=False,
#            del_value=del_value,
#            empty_value=empty_value,
#        )
#    # create MultisamplePileup object and postprocess
#    mspileup = MultisamplePileup(pileup_dict)
#    mspileup.set_df()
#
#    return mspileup




import re
import itertools
import collections
import functools

import numpy as np
import pandas as pd
#import pyranges as pr

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
readhandler = importlib.import_module(
    ".".join([top_package_name, "read", "readhandler"])
)
alignhandler = importlib.import_module(".".join([top_package_name, "align", "alignhandler"]))
librealign = importlib.import_module(".".join([top_package_name, "align", "realign"]))
#fetchcache = importlib.import_module(".".join([top_package_name, "read", "fetchcache"]))


DEL_VALUE = "*"
EMPTY_VALUE = ""
#DEFAULT_EXTEND_FETCHEDREADS_BY = 300
DEFAULT_EXTEND_PILEUP_BY = 20


class PileupBase:
    def _coord_arg_sanitycheck(self, start0, end0):
        #if start0 is not None:
        if start0 < self.start0:
            raise Exception('Input "start0" argument is out of pileup range.')
        #if end0 is not None:
        if end0 > self.end0:
            raise Exception('Input "end0" argument is out of pileup range.')

    def _coord_arg_sanitycheck_pos0(self, pos0):
        if pos0 not in self.df.columns:
            raise Exception('Input "pos0" argument is out of pileup range.')

    @property
    def range0(self):
        return range(self.start0, self.end0)

    def get_ref_seq(self):
        return self.fasta.fetch(self.chrom, self.start0, self.end0)

#    @staticmethod
#    def get_max_vaf_from_pileupcol(col, with_counts=False, empty_value=EMPTY_VALUE):
#        counts = collections.Counter(x for x in col if x != empty_value)
#        total = sum(counts.values())
#        if total == 0:
#            raise Exception(f"Pileup column allele count sum is 0.")
#
#        max_allele, max_count = max(counts.items(), key=(lambda x: x[1]))
#        max_vaf = max_count / total
#
#        if with_counts:
#            return max_allele, max_vaf, counts
#        else:
#            return max_allele, max_vaf
#
#    def get_max_vaf_colindex(self, col_idx, with_counts=False):
#        col = self.df.iloc[:, col_idx]
#        return self.get_max_vaf_from_pileupcol(col, with_counts=with_counts)
#
#    def get_max_vaf_colpos0(self, col_pos0, with_counts=False):
#        col = self.df.loc[:, col_pos0]
#        return self.get_max_vaf_from_pileupcol(col, with_counts=with_counts)

    def get_allele_counter(self, pos0):
        #self._coord_arg_sanitycheck_pos0(pos0)
        col = self.df.loc[:, pos0]
        return collections.Counter(x for x in col if x != EMPTY_VALUE)

    def get_allele_portions(self, pos0):
        counter = self.get_allele_counter(pos0)
        counter_sum = sum(counter.values())
        portions = dict()
        for key, val in counter.items():
            portions[key] = val / counter_sum
        return portions

    # active info related ones #
    def iter_active_info(self, start0, end0, reverse=False):
        self._coord_arg_sanitycheck(start0, end0)

        ref_seq = self.fasta.fetch(self.chrom, start0, end0)
        pos0_range = range(start0, end0)
        if reverse:
            ref_seq = ref_seq[::-1]
            pos0_range = pos0_range[::-1]

        for ref_base, pos0 in zip(ref_seq, pos0_range):
            counter = self.get_allele_counter(pos0)
            non_ref_portion = 1 - (counter[ref_base] / sum(counter.values()))
            is_active = (non_ref_portion >= self.active_threshold)
            yield is_active 

    def _set_active_info(self, start0=None, end0=None):
        if start0 is None:
            start0 = self.start0
        if end0 is None:
            end0 = self.end0

        self._coord_arg_sanitycheck(start0, end0)
        if not hasattr(self, '_active_info'):
            setattr(self, '_active_info', list())

        active_info_idx_start = start0 - self.start0
        active_info_idx_end = end0 - self.start0

        self._active_info[
            active_info_idx_start:active_info_idx_end
        ] = self.iter_active_info(start0, end0, reverse=False)

    def get_active_info(self, start0=None, end0=None):
        if start0 is None:
            start0 = self.start0
        if end0 is None:
            end0 = self.end0

        self._coord_arg_sanitycheck(start0, end0)

        active_info_idx_start = start0 - self.start0
        active_info_idx_end = end0 - self.start0

        return pd.Series(
            self._active_info[active_info_idx_start:active_info_idx_end], 
            index=self.df.columns[active_info_idx_start:active_info_idx_end],
        )

    def get_active_positions(self):
        return tuple(
            itertools.compress(self.df.columns, self._active_info)
        )

    def check_inactive_margin_right(self, length):
        if len(self._active_info) < length:
            return False
        else:
            return not any(self._active_info[-length:])

    def check_inactive_margin_left(self, length):
        if len(self._active_info) < length:
            return False
        else:
            return not any(self._active_info[:length])

    # row specs related ones
    @staticmethod
    def row_spec_matcher(query, target):
        if query['left_filled']:
            if query['right_filled']:
                if target['left_filled'] and target['right_filled']:
                    return query['seq'] in target['seq']
                else:
                    return False
            else:
                if target['left_filled']:
                    return query['seq'] == target['seq'][:len(query['seq'])]
                else:
                    return False
        else:
            if query['right_filled']:
                if target['right_filled']:
                    return query['seq'] == target['seq'][-len(query['seq']):]
                else:
                    return False
            else:
                return query['seq'] in target['seq']

    @staticmethod
    def handle_matching_subseq(superseq_candidates, row_spec_groups, subseq_rowid):
        if len(superseq_candidates) == 1:
            superseq_rowid = superseq_candidates[0]
        elif len(superseq_candidates) > 1:
            superseq_rowid = max(
                superseq_candidates, 
                key=(lambda x: row_spec_groups[x]['subseq_hits']),
            )

        row_spec_groups[superseq_rowid]['subseq_rowids'].add(subseq_rowid)
        row_spec_groups[superseq_rowid]['subseq_hits'] += 1

    @classmethod
    def group_row_specs(cls, row_specs):
        """row_spec's are grouped by whether one is a substring of another."""
        row_spec_groups = collections.OrderedDict()
        for query_rowid, query in sorted(
            row_specs.items(), 
            key=(lambda x: x[1]['match_length']), 
            reverse=True,
        ):
            superseq_candidates = list()
            for superseq_rowid in row_spec_groups.keys():
                target = row_specs[superseq_rowid]
                if cls.row_spec_matcher(query, target):
                    superseq_candidates.append(superseq_rowid)
                #if len(superseq_candidates) >= 2:
                #    break
                    
            if len(superseq_candidates) == 0:
                row_spec_groups[query_rowid] = {
                    'subseq_rowids': set(),
                    'subseq_hits': 0,
                }
                row_spec_groups[query_rowid]['subseq_rowids'].add(query_rowid)
                row_spec_groups[query_rowid]['subseq_hits'] += 1
            else:
                cls.handle_matching_subseq(superseq_candidates, row_spec_groups, query_rowid)

        return row_spec_groups

    def set_refined_vcfspec(self, row_concat_dist_le=np.inf, ):
        """All vcfspecs in a row are considered to be phased and belonging to a haplotype.
        """
        self.refined_vcfspec = None


class Pileup(PileupBase):
    pat_insclip = re.compile("(\(.+\))?([^()]+)(\(.+\))?")
    pat_parenthesis = re.compile("[()]")
    row_spec_exluded_vals = {DEL_VALUE, EMPTY_VALUE}

    def __init__(self, df, chrom, fasta, bam, read_cache=None, active_threshold=None):
        self.df = df
        self.chrom = chrom
        self.fasta = fasta
        self.bam = bam
        self.read_cache = read_cache

        if active_threshold is None:
            self.active_threshold = librealign.DEFAULT_ACTIVE_THRESHOLD
        else:
            self.active_threshold = active_threshold

        #self._ref_seq_cache = dict()
        #self._alignment_cache = dict()
        self._vcfspec_cache = dict()

        self._set_active_info()

    @property
    def start0(self):
        return self.df.columns[0]

    @property
    def end0(self):
        return self.df.columns[-1] + 1

    def subset(self, start0, end0, inplace=False):
        self._coord_arg_sanitycheck(start0, end0)
        # subset dataframe
        new_df = self.df.loc[:, start0:(end0 - 1)]
        # remove out-of-range rows
        row_within_range = list()
        for row_id, row in new_df.iterrows():
            read = self.read_cache[row_id]
            if read.reference_end <= start0 or read.reference_start >= end0:
                row_within_range.append(False)
            else:
                row_within_range.append(True)
        new_df = new_df.loc[row_within_range, :]
        # active_info
        active_info_idx_start = start0 - self.start0
        active_info_idx_end = end0 - self.start0
        new_active_info = self._active_info[active_info_idx_start:active_info_idx_end]

        # result
        if inplace:
            self.df = new_df
            self._active_info = new_active_info
        else:
            result = self.__class__(
                df=new_df, 
                chrom=self.chrom, 
                fasta=self.fasta,
                bam=self.bam,
                active_threshold=self.active_threshold,
            )
            result._active_info = new_active_info
            #result._ref_seq_cache = self._ref_seq_cache
            return result

    def get_read(self, read_uid):
        return self.read_cache[read_uid]

    ### extend ###  
    def extend_rightward(self, width):
        # keep original starts and ends
        border_col_idx_left = self.df.shape[1] - 1
        # make new pileup
        new_pileup = self._extend_helper_make_new_pileup(
            chrom=self.chrom,
            start0=self.end0,
            end0=(self.end0 + width),
        )
        # join
        self.df = self.df.join(new_pileup.df, how="outer")
        self.df[self.df.isnull()] = EMPTY_VALUE  # turn NaN into EMPTY_VALUE
        self._active_info.extend(new_pileup._active_info)
        self.read_cache.update(new_pileup.read_cache)
        # handle insclips facing each other
        self._extend_helper_handle_facing_insclip(border_col_idx_left)

    def extend_leftward(self, width):
        # keep original starts and ends
        border_col_idx_left = width - 1
        # make new pileup
        new_pileup = self._extend_helper_make_new_pileup(
            chrom=self.chrom,
            start0=(self.start0 - width),
            end0=self.start0,
        )
        # join
        self.df = new_pileup.df.join(self.df, how="outer")
        self.df[self.df.isnull()] = EMPTY_VALUE  # turn NaN into EMPTY_VALUE
        self._active_info = new_pileup._active_info + self._active_info
        self.read_cache.update(new_pileup.read_cache)
        # handle insclips facing each other
        self._extend_helper_handle_facing_insclip(border_col_idx_left)

    def _extend_helper_handle_facing_insclip(self, border_col_idx_left):
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
                ] = self.__class__.pat_insclip.sub("\\2\\3", "", val_right)

    def _extend_helper_make_new_pileup(self, chrom, start0, end0):
        return get_pileup(
            chrom, 
            start0, 
            end0,
            bam=self.bam,
            fasta=self.fasta,
            active_threshold=self.active_threshold,
            truncate=True,
            as_array=False,
            return_range=False,
        )

    ################################
    # realignment-related features #
    ################################

    # secure_inactive_padding
    def secure_inactive_padding_rightward(self, inactive_padding, extend_pileup_by=DEFAULT_EXTEND_PILEUP_BY):
        librealign.secure_inactive_padding_rightward(
            pileup=self,
            inactive_padding=inactive_padding, 
            extend_pileup_by=extend_pileup_by, 
            inplace=True,
        )

    def secure_inactive_padding_leftward(self, inactive_padding, extend_pileup_by=DEFAULT_EXTEND_PILEUP_BY):
        librealign.secure_inactive_padding_leftward(
            pileup=self, 
            inactive_padding=inactive_padding, 
            extend_pileup_by=extend_pileup_by, 
            inplace=True,
        )

    def secure_inactive_padding(self, inactive_padding, extend_pileup_by=DEFAULT_EXTEND_PILEUP_BY):
        self.secure_inactive_padding_rightward(inactive_padding, extend_pileup_by=extend_pileup_by)
        self.secure_inactive_padding_leftward(inactive_padding, extend_pileup_by=extend_pileup_by)

    ### row specs ###
    def get_row_spec(self, read_uid):
        row = self.df.loc[read_uid, :]
        return {
            "seq": "".join(
                self.__class__.pat_parenthesis.sub("", x) 
                for x in row if x not in self.__class__.row_spec_exluded_vals
            ),
            "match_length": sum(x != EMPTY_VALUE for x in row),
            "left_filled": (row.iloc[0] != EMPTY_VALUE),
            "right_filled": (row.iloc[-1] != EMPTY_VALUE),
            "id": row.name,
        }

    def set_row_specs(self):
        self.row_specs = dict()
        for read_uid, row in self.df.iterrows():
            self.row_specs[read_uid] = self.get_row_spec(read_uid)

    def set_row_spec_groups(self):
        self.row_spec_groups = self.__class__.group_row_specs(self.row_specs)

    # others
    def save_superseq_alignments(self, aligner=None):
        librealign.save_superseq_alignments(self, aligner=aligner)

    def set_vcfspecs(self, allele_portion_threshold=None, concat_dist_le=None):
        self.vcfspecs = librealign.get_vcfspecs_from_pileup(self, allele_portion_threshold=allele_portion_threshold, concat_dist_le=concat_dist_le)

    def get_realigned_reads(self):
        """Returns:
            dict (keys ReadUID, values pysam.AlignedSegment)
        """
        return librealign.get_realigned_reads(self)

    # for debugging
    def show_row_spec_alignments(self):
        subseq_hits_sum = sum(x['subseq_hits'] for x in self.row_spec_groups.values())
        for superseq_id, groupinfo in sorted(
            self.row_spec_groups.items(),
            key=(lambda x: x[1]['subseq_hits']),
            reverse=True,
        ):
            superseq_row_spec = self.row_specs[superseq_id]
            lstrip_query_gaps = not superseq_row_spec["left_filled"]
            rstrip_query_gaps = not superseq_row_spec["right_filled"]

            superseq_aln = self._alignment_cache[superseq_id]
            superseq_vcfspec_list = alignhandler.alignment_to_vcfspec(
                superseq_aln, 
                self.start0, 
                self.chrom, 
                self.fasta,
                lstrip_query_gaps=lstrip_query_gaps,
                rstrip_query_gaps=rstrip_query_gaps,
            )

            print('subseq hits:', groupinfo['subseq_hits'])
            print('subseq hits portion:', groupinfo['subseq_hits'] / subseq_hits_sum)
            print(superseq_row_spec)
            print(superseq_vcfspec_list)
            print(superseq_aln)
            print()


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
        self._set_active_info()

    def get_read(self, row_id):
        sampleid, read_uid = row_id
        return self.pileup_dict[sampleid].read_cache[read_uid]

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
    def show_row_spec_alignments(self, show_subseqs=True):
        for sampleid, pileup in self.pileup_dict.items():
            print('@@@', sampleid, '@@@')
            print()
            row_spec_groups = self.row_spec_groups[sampleid]
            subseq_hits_sum = sum(x['subseq_hits'] for x in row_spec_groups.values())
            for superseq_id, groupinfo in sorted(
                row_spec_groups.items(),
                key=(lambda x: x[1]['subseq_hits']),
                reverse=True,
            ):
                if groupinfo['subseq_hits'] == 0:
                    continue

                superseq_row_spec = self.superseq_row_specs[superseq_id]
                lstrip_query_gaps = not superseq_row_spec["left_filled"]
                rstrip_query_gaps = not superseq_row_spec["right_filled"]

                superseq_aln = self._alignment_cache[superseq_id]
                superseq_vcfspec_list = alignhandler.alignment_to_vcfspec(
                    superseq_aln, 
                    self.start0, 
                    self.chrom, 
                    self.fasta,
                    lstrip_query_gaps=lstrip_query_gaps,
                    rstrip_query_gaps=rstrip_query_gaps,
                )

                print('subseq hits:', groupinfo['subseq_hits'])
                print('subseq hits portion:', groupinfo['subseq_hits'] / subseq_hits_sum)
                print(superseq_row_spec)
                print(superseq_vcfspec_list)
                print(superseq_aln)
                if show_subseqs:
                    for readuid in groupinfo['subseq_rowids']:
                        print(readuid)
                print()


def get_pileup(
    chrom,
    start0,
    end0,
    bam,
    fasta=None,
    active_threshold=None,
    truncate=True,
    as_array=False,
    return_range=False,
    del_value=DEL_VALUE,
    empty_value=EMPTY_VALUE,
):
    def sanity_check():
        if (not as_array) and (fasta is None):
            raise Exception(f'If "as_array" is False, "fasta" must be set.')

    def make_readlist(bam, chrom, start0, end0):
        readlist = list()
        start_list = list()
        end_list = list()
        uid_list = list()
        for read in readhandler.get_fetch(bam, chrom, start0, end0, readfilter=readhandler.readfilter_pileup):
            uid = readhandler.get_uid(read)
            readlist.append(read)
            start_list.append(read.reference_start)
            end_list.append(read.reference_end)
            uid_list.append(uid)
        tmp_pileup_range = range(min(start_list), max(end_list))

        return readlist, tmp_pileup_range, uid_list

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

    def write_read_to_array_row(read, arr, arr_row_idx, initial_arr_col_idx, del_value):
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
            arr[arr_row_idx, arr_col_idx] = f"({read.query_sequence})"

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
                raise Exception(
                    f"Consecutive Non-M cigar units are composed of both target-only and query-only ones:\n{read.to_string()}"
                )

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
            # cigartup is only composed of query-only operations (no M)
            arr[arr_row_idx, arr_col_idx] = "(" + read.query_sequence + ")"            
        else:
            if read.cigartuples[cigartup_idx][0] != 0:
                raise Exception(
                    f"The first cigar operation after stripping leading I and S is not M:\n{read.to_string()}"
                )

            # begins loop
            while True:
                if cigartup_idx > len(read.cigartuples) - 1:
                    raise Exception(
                        f'"cigartup_idx" became greater than "len(read.cigartuples) - 1" while looping.'
                    )

                cigarop, cigarlen = read.cigartuples[cigartup_idx]
                assert (
                    cigarop == 0
                ), f'The cigar unit indicated by "cigarop_idx" at the beginning of a loop cycle is not M.'

                if (
                    cigartup_idx == len(read.cigartuples) - 1
                ):  # current cigartup is the last one
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

    def truncate_array(arr, tmp_pileup_range, pileup_range):
        sl_start = tmp_pileup_range.index(pileup_range.start)
        if pileup_range.stop == tmp_pileup_range.stop:
            sl_end = None
        else:
            sl_end = tmp_pileup_range.index(pileup_range.stop)
        return arr[:, slice(sl_start, sl_end)]

    def make_pileup_from_array(
        chrom,
        arr,
        bam,
        read_cache,
        uid_list,
        pileup_range,
        tmp_pileup_range,
        truncate,
        active_threshold,
        fasta,
    ):
        df_columns = pileup_range if truncate else tmp_pileup_range
        df = pd.DataFrame(arr, columns=list(df_columns), index=uid_list)
        pileup = Pileup(
            df=df,
            chrom=chrom,
            fasta=fasta,
            bam=bam,
            read_cache=read_cache,
            active_threshold=active_threshold,
        )
        return pileup

    # main
    # parameter handling
    sanity_check()
    if active_threshold is None:
        active_threshold = librealign.DEFAULT_ACTIVE_THRESHOLD
    # prepare pileup_range
    pileup_range = range(start0, end0)
    readlist, tmp_pileup_range, uid_list = make_readlist(bam, chrom, start0, end0)
    # create array
    arr = np.full(
        shape=(len(readlist), len(tmp_pileup_range)),
        fill_value=empty_value,
        dtype=object,
    )
    for arr_row_idx, read in enumerate(readlist):
        initial_arr_col_idx = read.reference_start - tmp_pileup_range.start
        try:
            write_read_to_array_row(read, arr, arr_row_idx, initial_arr_col_idx, del_value)
        except Exception as exc:
            try:
                cigar_sanitycheck(read)
            except Exception as exc_cigarpattern:
                raise exc_cigarpattern from exc
            else:
                raise exc
    # truncate array
    if truncate:
        arr = truncate_array(arr, tmp_pileup_range, pileup_range)
    # prepare results
    if as_array:
        if return_range:
            return (arr, pileup_range)
        else:
            return arr
    else:
        read_cache = dict(zip(uid_list, readlist))
        pileup = make_pileup_from_array(
            chrom,
            arr,
            bam,
            read_cache,
            uid_list,
            pileup_range,
            tmp_pileup_range,
            truncate,
            active_threshold,
            fasta,
        )
        if return_range:
            return (pileup, pileup_range)
        else:
            return pileup


def get_pileup_multisample(
    chrom,
    start0,
    end0,
    bam_dict,
    fasta,
    active_threshold_onesample=None,
    del_value=DEL_VALUE,
    empty_value=EMPTY_VALUE,
):
    # parameter handling
    if active_threshold_onesample is None:
        active_threshold_onesample = librealign.DEFAULT_ACTIVE_THRESHOLD
    # initialize pileup_dict
    pileup_dict = dict()
    for sampleid, bam in bam_dict.items():
        pileup_dict[sampleid] = get_pileup(
            chrom,
            start0,
            end0,
            bam=bam,
            fasta=fasta,
            active_threshold=active_threshold_onesample,
            truncate=True,
            as_array=False,
            return_range=False,
            del_value=del_value,
            empty_value=empty_value,
        )
    # create MultisamplePileup object and postprocess
    mspileup = MultisamplePileup(pileup_dict)
    mspileup.set_df()

    return mspileup




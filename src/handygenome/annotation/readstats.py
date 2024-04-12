import os
import warnings
import re
import logging
import copy
import pprint
import collections
import functools

import pysam
import numpy as np
import scipy.stats

import handygenome.refgenome.refgenome as refgenome
import handygenome.logutils as logutils
import handygenome.workflow as workflow
from handygenome.annotation.annotitem import AnnotItemFormatSingle, AnnotItemFormatSampledict
import handygenome.variant.infoformat as infoformat
from handygenome.variant.vcfspec import Vcfspec
import handygenome.read.readplus as readplus
from handygenome.read.readplus import ReadPlusPairList
import handygenome.read.alleleclass as liballeleclass
import handygenome.read.alleleclass_sv as liballeleclass_sv
from handygenome.sv.breakends import Breakends


ZERO_ONE_UNIFORM = scipy.stats.uniform(loc=0, scale=1)
DEFAULT_DEPTH_LIMITS = (0, -1)
DEFAULT_MQ_LIMITS = (0, -1)


class AlleleclassError(Exception):
    pass


def get_position_info(bam, chrom, pos0):
    mqlist = list()
    idx = -1
    for idx, read in enumerate(bam.fetch(chrom, pos0, pos0 + 1)):
        mqlist.append(read.mapping_quality)
    depth = idx + 1
    return depth, mqlist


def get_position_info_pileup(bam, chrom, pos0):
    pup = bam.pileup(
        chrom, pos0, pos0 + 1, 
        truncate=True, 
        stepper='nofilter', 
        ignore_overlaps=False, 
        flag_filter=0, 
        ignore_orphans=False, 
        min_base_quality=0, 
        min_mapping_quality=0,
    )
    pupcol = next(pup)
    mqlist = pupcol.get_mapping_qualities()
    depth = pupcol.get_num_aligned()
    return depth, mqlist


class ReadStats(AnnotItemFormatSingle):
    def __init__(self, is_missing=False, is_pon=False):
        super().__init__(is_missing=is_missing)
        self['is_pon'] = is_pon

    @staticmethod
    def handle_limits_arg(limits):
        result = list(limits)
        if result[1] == -1:
            result[1] = np.inf
        return result

    @staticmethod
    def check_position_validity(depth, mqlist, depth_limits, mq_limits):
        mean_mq = np.mean(mqlist) if mqlist else np.nan
        depth_okay = (depth >= depth_limits[0]) and (depth <= depth_limits[1])
        mq_okay = (
            np.isnan(mean_mq)
            or ((mean_mq >= mq_limits[0]) and (mean_mq <= mq_limits[1]))
        )
        return depth_okay and mq_okay

    @property
    def fasta(self):
        return self.vcfspec.fasta

    @property
    def chromdict(self):
        return self.vcfspec.chromdict

    ################
    # initializers #
    ################
        
    @classmethod
    def from_bam(
        cls, 
        vcfspec, 
        bam, 
        rpplist_kwargs=dict(),
        alleleclass_kwargs=dict(),
        countonly=False,
        depth_limits=DEFAULT_DEPTH_LIMITS,
        mq_limits=DEFAULT_MQ_LIMITS,
        include_mNM_items=False,
        is_pon=False,
    ):
        depth, mqlist = get_position_info(bam, vcfspec.chrom, vcfspec.pos0)
        depth_limits = cls.handle_limits_arg(depth_limits)
        mq_limits = cls.handle_limits_arg(mq_limits)

        if cls.check_position_validity(depth, mqlist, depth_limits, mq_limits):
            if not vcfspec.check_is_sv():
                readstats_data = get_readstats_data(
                    vcfspec, 
                    bam, 
                    rpplist_kwargs=rpplist_kwargs,
                    alleleclass_kwargs=alleleclass_kwargs,
                    countonly=countonly,
                )
                result = cls.from_readstats_data(
                    readstats_data, 
                    vcfspec, 
                    countonly=countonly,
                    include_mNM_items=include_mNM_items,
                    is_pon=is_pon,
                )
                del readstats_data
            else:  # SV
                readstats_data_bnd1, readstats_data_bnd2 = get_readstats_data(
                    vcfspec, 
                    bam, 
                    rpplist_kwargs=rpplist_kwargs,
                    alleleclass_kwargs=alleleclass_kwargs,
                    countonly=countonly,
                )
                result_bnd1 = cls.from_readstats_data(
                    readstats_data_bnd1, 
                    vcfspec, 
                    countonly=True,
                    include_mNM_items=include_mNM_items,
                    is_pon=is_pon,
                )
                result_bnd2 = cls.from_readstats_data(
                    readstats_data_bnd2, 
                    vcfspec, 
                    countonly=True,
                    include_mNM_items=include_mNM_items,
                    is_pon=is_pon,
                )
                del readstats_data_bnd1
                del readstats_data_bnd2

                result = cls(is_missing=False)
                result.vcfspec = vcfspec
                result['bnd1'] = result_bnd1
                result['bnd2'] = result_bnd2
        else:
            if not vcfspec.check_is_sv():
                result = cls.init_invalid(vcfspec, countonly=countonly, is_pon=is_pon)
            else:
                result = cls(is_missing=False)
                result['bnd1'] = cls.init_invalid(vcfspec, countonly=countonly, is_pon=is_pon)
                result['bnd2'] = cls.init_invalid(vcfspec, countonly=countonly, is_pon=is_pon)

        return result

    @classmethod
    def from_readstats_data(
        cls, 
        readstats_data,  
        vcfspec, 
        countonly=False,
        include_mNM_items=False,
        is_pon=False,
    ):
        # main
        result = cls(is_missing=False, is_pon=is_pon)
        result.is_invalid = False
        result.vcfspec = vcfspec

        result['rppcounts'] = readstats_data['count'].copy()
        if not countonly:
            for pf in ('clipped', 'f1r2', 'f2r1', 'IDcontext'):
                result[f'{pf}_rppcounts'] = readstats_data[f'{pf}_count'].copy()

        if not countonly:
            result['mean_BQs'] = cls.allele_means(readstats_data['BQ'])
            result['median_BQs'] = cls.allele_medians(readstats_data['BQ'])

            result['mean_MQs'] = cls.allele_means(readstats_data['MQ'])
            result['median_MQs'] = cls.allele_medians(readstats_data['MQ'])

            result['mean_NMs'] = cls.allele_means(readstats_data['NM'])
            result['median_NMs'] = cls.allele_medians(readstats_data['NM'])

            result['mean_cliplens'] = cls.allele_means(readstats_data['cliplen'])
            result['median_cliplens'] = cls.allele_medians(readstats_data['cliplen'])

            (
                result['mNM'], 
                result['recurrent_mNM'], 
                mNMitems_norecur, 
                mNMitems_recur,
            ) = cls.handle_mNM(
                readstats_data['mNM'], 
                readstats_data['count'], 
                readstats_data['rpp_range0'], 
                recurrence_cutoff_fraction=0.9,
                recurrence_cutoff_denominator=3,
            )
            if include_mNM_items:
                result['mNM_items'] = mNMitems_norecur
                result['recurrent_mNM_items'] = mNMitems_recur

            result['pairorient_pvalues'] = cls.handle_orient(readstats_data['pairorient'], mode='pairorient')
            result['readorient_pvalues'] = cls.handle_orient(readstats_data['readorient'], mode='pairorient')

            result['varpos_uniform_pvalues'] = cls.varpos_kstest(
                readstats_data['pos0_left_fraction']
            )
            result['mean_varpos_fractions'] = cls.allele_means(readstats_data['pos0_left_fraction'])
            result['median_varpos_fractions'] = cls.allele_medians(readstats_data['pos0_left_fraction'])
            result['std_varpos_fractions'] = cls.allele_stds(readstats_data['pos0_left_fraction'])

        return result

    @classmethod
    def from_custom_countonly(cls, vcfspec, counts, is_pon=False):
        result = cls(is_missing=False, is_pon=is_pon)
        result.is_invalid = False
        result.vcfspec = vcfspec

        result['rppcounts'] = dict()
        result['rppcounts']['softclip_overlap'] = 0
        for alleleclass in vcfspec.get_alleleclasses():  # this includes -1
            try:
                n_read = counts[alleleclass]
            except KeyError:
                n_read = 0
            result['rppcounts'][alleleclass] = n_read

        return result

    @classmethod
    def init_invalid(cls, vcfspec, countonly=False, verbose=True, is_pon=False):
        if verbose:
            logutils.log(f'Initiating ReadStats object as invalid mode', level='warning')

        # main
        result = cls(is_missing=False, is_pon=is_pon)
        result.is_invalid = True
        result.vcfspec = vcfspec

        result['rppcounts'] = {
            alleleclass: np.nan 
            for alleleclass in vcfspec.get_alleleclasses()
        }
        result['rppcounts']['softclip_overlap'] = np.nan

        if not countonly:
            for key in (
                'mean_BQs',
                'median_BQs',
                'mean_MQs',
                'median_MQs',
                'mean_cliplens',
                'median_cliplens',
                'mNM',
                'recurrent_mNM',
                'pairorient_pvalues',
                'readorient_pvalues',
                'varpos_uniform_pvalues',
                'mean_varpos_fractions',
                'median_varpos_fractions',
            ):
                result[key] = {
                    alleleclass: np.nan 
                    for alleleclass in vcfspec.get_alleleclasses()
                }

        return result

    ####################################
    # helpers of "from_readstats_data" #
    ####################################

    @staticmethod
    def allele_stats(data, func, *args, **kwargs):
        return dict(
            (
                alleleclass, 
                np.nan if len(values) == 0 else func(values, *args, **kwargs),
            )
            for (alleleclass, values) in data.items()
        )

    @classmethod
    def allele_means(cls, data):
        return cls.allele_stats(data, np.mean)

    @classmethod
    def allele_medians(cls, data):
        return cls.allele_stats(data, np.median)

    @classmethod
    def allele_stds(cls, data, ddof=0):
        return cls.allele_stats(data, np.std, ddof=ddof)

    @staticmethod
    def varpos_kstest(data):
        summary = dict()
        for alleleclass, values in data.items():
            if len(values) == 0:
                pval = np.nan
            else:
                pval = scipy.stats.kstest(values, ZERO_ONE_UNIFORM.cdf).pvalue
            summary[alleleclass] = pval

        return summary

    @staticmethod
    def handle_orient(data, mode):
        summary = dict()
        for alleleclass, values in data.items():
            if len(values) == 0:
                summary[alleleclass] = np.nan
            else:
                counter = collections.Counter(values)
                if mode == 'pairorient':
                    x = (counter['f1r2'], counter['f2r1'])
                elif mode == 'readorient':
                    x = (counter['f'], counter['r'])
                summary[alleleclass] = scipy.stats.binom_test(x, p=0.5, alternative='two-sided')

        return summary

    @staticmethod
    def handle_mNM(
        data, 
        rppcounts_dict, 
        rpp_range0_dict, 
        recurrence_cutoff_fraction,
        recurrence_cutoff_denominator,
    ):
        """Difference from ordinary NM: recurrent and non-recurrent mismatches are separated.
        """
        summary_without_recurrence = dict()
        summary_only_recurrence = dict()
        mNMitems_norecur = dict()
        mNMitems_recur = dict()
        for alleleclass, values in data.items():    
            # set values
            mNM_count_sum_with_recurrence = 0
            mNM_count_sum_without_recurrence = 0

            mNMitems_norecur[alleleclass] = list()
            mNMitems_recur[alleleclass] = list()

            if len(values) > 0:
                counter = collections.Counter(values)
                for key, val in counter.items():
                    refpos0 = key[0]
                    relevant_rppcount = sum(
                        any((refpos0 in x) for x in range0_list)
                        for range0_list in rpp_range0_dict[alleleclass]
                    )
                    if relevant_rppcount == 0:
                        raise Exception(
                            f'The number of relevant rpp is 0. mNM element: {key}'
                        )

                    if (
                        ((val / relevant_rppcount) >= recurrence_cutoff_fraction)
                        and (relevant_rppcount >= recurrence_cutoff_denominator)
                    ):  # This means a recurrent MM
                        mNM_count_sum_with_recurrence += 1
                        mNMitems_recur[alleleclass].append(key)
                    else:
                        mNM_count_sum_without_recurrence += 1
                        mNMitems_norecur[alleleclass].append(key)

            # summarize
            if rppcounts_dict[alleleclass] == 0:
                summary_without_recurrence[alleleclass] = np.nan
            else:
                summary_without_recurrence[alleleclass] = (
                    mNM_count_sum_without_recurrence / rppcounts_dict[alleleclass]
                )
            summary_only_recurrence[alleleclass] = mNM_count_sum_with_recurrence

        return (
            summary_without_recurrence, 
            summary_only_recurrence, 
            mNMitems_norecur, 
            mNMitems_recur,
        )

    #def write(self, vr, sampleid):
    #    return self.write_base(vr, sampleid)

    #############
    # utilities #
    #############

    def get_allele_indexes_average(self, key, allele_indexes=None):
        assert key != 'rppcounts'
        if allele_indexes is None:
            allele_indexes = tuple(self[key].keys())

        values = np.fromiter((self[key][x] for x in allele_indexes), dtype=float)
        weights = np.fromiter((self['rppcounts'][x] for x in allele_indexes), dtype=float)

        nan_indexes = np.logical_or(np.isnan(values), np.isnan(weights))
        values = values[~nan_indexes]
        weights = weights[~nan_indexes]

        if sum(weights) == 0:
            return np.nan
        else:
            return np.average(values, weights=weights)

    get_allele_indexes_mean = get_allele_indexes_average
    #get_alleleindexes_mean = get_allele_indexes_average

    def get_total_rppcount(self, exclude_other=False):
        if exclude_other:
            return sum(
                val for (key, val) in self['rppcounts'].items()
                if isinstance(key, int) and key != -1
            )
        else:
            return sum(
                val for (key, val) in self['rppcounts'].items()
                if isinstance(key, int)
            )

    def get_vaf(self, alleleclass=1, exclude_other=False):
        if (alleleclass == -1) and exclude_other:
            raise Exception(f'When "alleleclass" is -1, "exclude_other" must be False.')

        total_rppcount = self.get_total_rppcount(exclude_other=exclude_other)
        if total_rppcount == 0:
            return np.nan
        else:
            try:
                return self['rppcounts'][alleleclass] / total_rppcount
            except KeyError:
                raise AlleleclassError(f'Invalid alleleclass')

    def get_vafs(self, n_allele=None, exclude_other=False):
        if n_allele is None:
            n_allele = sum(
                (isinstance(x, int) and (x >= 0)) 
                for x in self['rppcounts'].keys()
            )
                
        total_rppcount = self.get_total_rppcount(exclude_other=exclude_other)
        if total_rppcount == 0:
            return np.repeat(np.nan, n_allele)
        else:
            data = list()
            for alleleclass in range(n_allele):
                try:
                    vaf = self['rppcounts'][alleleclass] / total_rppcount
                except KeyError:
                    vaf = np.nan
                data.append(vaf)
            return np.array(data)

    def get_sorted_vafs(self, exclude_other=True, reverse=True):
        vafs = list()
        for alleleclass in self['rppcounts'].keys():
            if alleleclass is None:
                continue
            if exclude_other:
                if alleleclass == -1:
                    continue
            vafs.append(
                self.get_vaf(alleleclass=alleleclass, exclude_other=exclude_other)
            )

        return sorted(vafs, reverse=reverse)

    #def get_all_allele_rppcounts(self):
    #    result = dict()
        #total_rppcount = self.get_total_rppcount()
    #    for alleleclass, rppcount in self['rppcounts'].items():
    #        result[alleleclass] = rppcount / total_rppcount
    #    return result


class ReadStatsSampledict(AnnotItemFormatSampledict):
    meta = {
        'ID': 'readstats', 'Number': 1, 'Type': 'String',
        'Description': 'Read-based statistics for the variant.',
    }
    unit_class = ReadStats

    @classmethod
    def init_missing(cls, vr):
        result = super().init_missing(vr)
        result.vcfspec = Vcfspec.from_vr(vr)
        return result

    @staticmethod
    def handle_limits_arg(limits, bam_dict):
        if isinstance(limits, (tuple, list)):
            new_limits = {key: list(limits) for key in bam_dict.keys()}
        elif isinstance(limits, dict):
            #if set(limits.keys()) != set(bam_dict.keys()):
            if not set(bam_dict.keys()).issubset(set(limits.keys())):
                raise Exception(f'Keys of "limits" argument are not a superset of the keys of "bam_dict".')
            new_limits = limits
        else:
            raise Exception(f'"limits" argument must be either a tuple, list, or dict')

        return new_limits

    @classmethod
    def from_bam_dict(
        cls, 
        bam_dict, 
        vcfspec, 
        rpplist_kwargs=dict(),
        alleleclass_kwargs=dict(),
        countonly=False,
        depth_limits=DEFAULT_DEPTH_LIMITS,
        mq_limits=DEFAULT_MQ_LIMITS,
        init_invalid=False,
        include_mNM_items=False,
        verbose=False,
        pon_samples=list(),
    ):
        depth_limits = cls.handle_limits_arg(depth_limits, bam_dict)
        mq_limits = cls.handle_limits_arg(mq_limits, bam_dict)

        result = cls()
        result.vcfspec = vcfspec
        if init_invalid:
            for sampleid, bam in bam_dict.items():
                if not vcfspec.check_is_sv():
                    result[sampleid] = ReadStats.init_invalid(
                        vcfspec, countonly=countonly, is_pon=(sampleid in pon_samples),
                    )
                else:
                    result[sampleid] = ReadStats(is_missing=False)
                    result[sampleid]['bnd1'] = ReadStats.init_invalid(
                        vcfspec, countonly=countonly, is_pon=(sampleid in pon_samples),
                    )
                    result[sampleid]['bnd2'] = ReadStats.init_invalid(
                        vcfspec, countonly=countonly, is_pon=(sampleid in pon_samples),
                    )
        else:
            for sampleid, bam in bam_dict.items():
                if verbose:
                    logutils.log(f'Creating ReadStats for {repr(sampleid)}')
                result[sampleid] = ReadStats.from_bam(
                    vcfspec, 
                    bam, 
                    rpplist_kwargs=rpplist_kwargs,
                    alleleclass_kwargs=alleleclass_kwargs,
                    countonly=countonly,
                    depth_limits=depth_limits[sampleid],
                    mq_limits=mq_limits[sampleid],
                    include_mNM_items=include_mNM_items,
                    is_pon=(sampleid in pon_samples),
                )
        return result

    @classmethod
    def from_vr(cls, vr, sampleid_list=None):
        return cls.from_vr_base(vr, sampleid_list=sampleid_list)

    def get_first_readstats(self):
        return next(iter(self.values()))

    @property
    def refver(self):
        return self.vcfspec.refver

    @property
    def fasta(self):
        return refgenome.get_fasta(self.refver)

    @property
    def chromdict(self):
        return refgenome.get_chromdict(self.refver)

    def update_bams(
        self, 
        bam_dict,
        rpplist_kwargs=dict(),
        alleleclass_kwargs=dict(),
        countonly=False,
        depth_limits=DEFAULT_DEPTH_LIMITS,
        mq_limits=DEFAULT_MQ_LIMITS,
        init_invalid=False,
        include_mNM_items=False,
        pon_samples=list(),
    ):
        if len(self) == 0:
            raise Exception(f'Length of self must be greater than 0')

        depth_limits = self.handle_limits_arg(depth_limits, bam_dict)
        mq_limits = self.handle_limits_arg(mq_limits, bam_dict)

        if init_invalid:
            for sampleid, bam in bam_dict.items():
                self[sampleid] = ReadStats.init_invalid(
                    self.vcfspec, 
                    countonly=countonly, 
                    is_pon=(sampleid in pon_samples),
                )
        else:
            for sampleid, bam in bam_dict.items():
                self[sampleid] = ReadStats.from_bam(
                    self.vcfspec, 
                    bam, 
                    rpplist_kwargs=rpplist_kwargs,
                    alleleclass_kwargs=alleleclass_kwargs,
                    countonly=countonly,
                    depth_limits=depth_limits[sampleid],
                    mq_limits=mq_limits[sampleid],
                    include_mNM_items=include_mNM_items,
                    is_pon=(sampleid in pon_samples),
                )

    def write(self, vr, donot_write_missing=True):
        return self.write_base(vr, donot_write_missing=donot_write_missing)

    def show_rppcounts(self):
        for sampleid, readstats in self.items():
            print(sampleid)
            pprint.pprint(readstats['rppcounts'])

    def show_vafs(self):
        for sampleid, readstats in self.items():
            print(sampleid)
            vaf_dict = {
                alleleclass: readstats.get_vaf(alleleclass)
                for alleleclass in sorted(
                    readstats['rppcounts'].keys(), 
                    key=(
                        lambda x: (
                            -2 if x is None else 
                            (np.inf if x == 'softclip_overlap' else x)
                        )
                    )
                )
            }
            pprint.pprint(vaf_dict)


##################
# readstats data #
##################

def rpplist_to_readstats_data_countonly_sv(
    rpplist, bnds, 
    flanklen_parside=liballeleclass_sv.DEFAULT_FLANKLEN_PARSIDE,
    flanklen_bndside=liballeleclass_sv.DEFAULT_FLANKLEN_BNDSIDE,
):
    data_bnd1 = {
        'count': {x: 0 for x in (None, -1, 0, 1)},
    }
    data_bnd2 = {
        'count': {x: 0 for x in (None, -1, 0, 1)},
    }

    # add data
    for rpp in rpplist:
        if bnds not in rpp.alleleclass.keys():
            rpp.update_alleleclass_sv(
                bnds,
                flanklen_parside=flanklen_parside,
                flanklen_bndside=flanklen_bndside,
            )

    for rpp in rpplist:
        aiitem = rpp.alleleclass[bnds]
        data_bnd1['count'][None] += aiitem['noninformative']
        data_bnd2['count'][None] += aiitem['noninformative']

        data_bnd1['count'][-1] += aiitem['other_support']
        data_bnd2['count'][-1] += aiitem['other_support']

        data_bnd1['count'][0] += (
            aiitem['ref_support_direct_bnd1']
            + aiitem['ref_support_indirect_bnd1']
        )
        data_bnd2['count'][0] += (
            aiitem['ref_support_direct_bnd2']
            + aiitem['ref_support_indirect_bnd2']
        )

        data_bnd1['count'][1] += (
            aiitem['alt_support_direct']
            + aiitem['alt_support_indirect']
        )
        data_bnd2['count'][1] += (
            aiitem['alt_support_direct']
            + aiitem['alt_support_indirect']
        )

    return data_bnd1, data_bnd2


def add_var_pos0s_rp(rp, data, alleleclass_rpp, vcfspec):
    var_querypos0s = rp.get_querypos0_of_range0_allmodes(vcfspec.REF_range0)

    if var_querypos0s['left'] is not None:
        data['pos0_left'][alleleclass_rpp].append(var_querypos0s['left'])
        data['pos0_right'][alleleclass_rpp].append(var_querypos0s['right'])
        data['pos0_5prime'][alleleclass_rpp].append(var_querypos0s['5prime'])
        data['pos0_3prime'][alleleclass_rpp].append(var_querypos0s['3prime'])

        data['pos0_left_fraction'][alleleclass_rpp].append(var_querypos0s['left_fraction'])
        data['pos0_right_fraction'][alleleclass_rpp].append(var_querypos0s['right_fraction'])
        data['pos0_5prime_fraction'][alleleclass_rpp].append(var_querypos0s['5prime_fraction'])
        data['pos0_3prime_fraction'][alleleclass_rpp].append(var_querypos0s['3prime_fraction'])


def rpplist_to_readstats_data(
    rpplist, 
    vcfspec, 
    flanklen=liballeleclass.DEFAULT_FLANKLEN,
    countonly=False,
):
    # initialize
    alleleclass_keys = vcfspec.get_alleleclasses()
    data = dict()

    # fields initialized as integer
    data['count'] = {x: 0 for x in alleleclass_keys}
    if not countonly:
        data['count']['softclip_overlap'] = 0
        data['clipped_count'] = {x: 0 for x in alleleclass_keys}
        data['f1r2_count'] = {x: 0 for x in alleleclass_keys}
        data['f2r1_count'] = {x: 0 for x in alleleclass_keys}
        data['IDcontext_count'] = {x: 0 for x in alleleclass_keys}

    # fields initialized as list
    if not countonly:
        for key in (
            'rpp_range0',
            'MQ', 'BQ', 'NM', 'mNM', 
            #'clipspec', 
            'cliplen', 'pairorient', 'readorient',
            'pos0_left', 'pos0_right', 'pos0_5prime', 'pos0_3prime',
            'pos0_left_fraction', 'pos0_right_fraction', 
            'pos0_5prime_fraction', 'pos0_3prime_fraction',
        ):
            data[key] = {x: list() for x in alleleclass_keys}

    # add data
    for rpp in rpplist:
        if vcfspec not in rpp.alleleclass.keys():
            rpp.update_alleleclass(vcfspec=vcfspec, flanklen=flanklen)
        alleleclass_rpp = rpp.alleleclass[vcfspec]

        # count
        data['count'][alleleclass_rpp] += 1
            
        if not countonly:
            # softclip overlap
            if alleleclass_rpp is None:
                if rpp.check_softclip_overlaps_vcfspec(vcfspec):
                    data['count']['softclip_overlap'] += 1

            # count only softclip-harboring rpps
            if any(
                (rp.read.get_cigar_stats()[0][4] > 0)
                #(4 in [x[0] for x in rp.read.cigartuples])
                for rp in rpp.get_rp_iter()
            ):
                data['clipped_count'][alleleclass_rpp] += 1

            # count only f1r2 configuration
            if rpp.pairorient == 'f1r2':
                data['f1r2_count'][alleleclass_rpp] += 1
            # count only f2r1 configuration
            if rpp.pairorient == 'f2r1':
                data['f2r1_count'][alleleclass_rpp] += 1

            # count only rpps with indel around variant position
            hits = False
            for rp in rpp.get_rp_iter():
                mNM_data = rp.get_mNM_data()
                if any(
                    (
                        (x[1] in (1, 2))
                        and (abs(x[0] - vcfspec.start0) <= 10)
                    )
                    for x in mNM_data
                ):
                    hits = True
                    break
            if hits:
                data['IDcontext_count'][alleleclass_rpp] += 1

            # rpp_range0
            rpp_range0s = list()
            if rpp.rp1.read.reference_name == vcfspec.chrom:
                rpp_range0s.append(rpp.rp1.range0)
            if rpp.rp2 is not None:
                if rpp.rp2.read.reference_name == vcfspec.chrom:
                    rpp_range0s.append(rpp.rp2.range0)
            data['rpp_range0'][alleleclass_rpp].append(rpp_range0s)

            # MQ
            data['MQ'][alleleclass_rpp].append(rpp.get_MQ())

            # BQ
            data['BQ'][alleleclass_rpp].extend(rpp.rp1.get_BQlist(vcfspec))
            if rpp.rp2 is not None:
                data['BQ'][alleleclass_rpp].extend(rpp.rp2.get_BQlist(vcfspec))

            # cliplen
            data['cliplen'][alleleclass_rpp].append(rpp.get_cliplen())

            # pairorient
            data['pairorient'][alleleclass_rpp].append(rpp.pairorient)

            # readorient
            if rpp.rp1.alleleclass[vcfspec] == alleleclass_rpp:
                data['readorient'][alleleclass_rpp].append(
                    'f' if rpp.rp1.read.is_forward else 'r'
                )
            if rpp.rp2 is not None:
                if rpp.rp2.alleleclass[vcfspec] == alleleclass_rpp:
                    data['readorient'][alleleclass_rpp].append(
                        'f' if rpp.rp2.read.is_forward else 'r'
                    )

            # NM, mNM
            data['NM'][alleleclass_rpp].append(rpp.get_NM())
            #mNM_data, clipspec_data = rpp.get_mNM_clipspec_data(vcfspec)
            mNM_data = rpp.get_mNM_data()
            data['mNM'][alleleclass_rpp].extend(mNM_data)
            #data['clipspec'][alleleclass_rpp].extend(clipspec_data)

            # var_pos0s
            add_var_pos0s_rp(rpp.rp1, data, alleleclass_rpp, vcfspec)
            if rpp.rp2 is not None:
                add_var_pos0s_rp(rpp.rp2, data, alleleclass_rpp, vcfspec)

    return data


def get_readstats_data(
    vcfspec, bam, 
    rpplist_kwargs=dict(),
    alleleclass_kwargs=dict(),
    countonly=False,
):
    """Only for non-sv cases"""

    #rpplist_kwargs.update({'view': False})
    rpplist_kwargs['view'] = False
    if vcfspec.check_is_sv():
        bnds = Breakends.from_vcfspec(vcfspec)
        rpplist = ReadPlusPairList.from_bam_sv(
            bam=bam, 
            bnds=bnds,
            **rpplist_kwargs,
        )
        rpplist.update_alleleclass_sv(
            bnds=bnds,
            **alleleclass_kwargs,
        )
        (
            readstats_data_bnd1, 
            readstats_data_bnd2,
        ) = rpplist_to_readstats_data_countonly_sv(rpplist, bnds)
        del rpplist

        return readstats_data_bnd1, readstats_data_bnd2
    else:
        rpplist = ReadPlusPairList.from_bam(
            bam=bam, 
            chrom=vcfspec.chrom, 
            start0=vcfspec.pos0, 
            end0=vcfspec.end0, 
            **rpplist_kwargs,
        )
        rpplist.update_alleleclass(
            vcfspec=vcfspec, 
            **alleleclass_kwargs,
        )
        readstats_data = rpplist_to_readstats_data(rpplist, vcfspec, countonly=countonly)

        del rpplist

        return readstats_data




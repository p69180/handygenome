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

import handygenome.logutils as logutils
import handygenome.workflow as workflow
import handygenome.annotation.annotitem as annotitem
import handygenome.variant.infoformat as infoformat
import handygenome.read.readplus as readplus
import handygenome.read.alleleinfo as liballeleinfo


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


class ReadStats(annotitem.AnnotItemFormatSingle):
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
        
    @classmethod
    def from_bam(
        cls, 
        vcfspec, 
        bam, 
        fasta=None, 
        chromdict=None,
        rpplist_kwargs=dict(),
        alleleinfo_kwargs=dict(),
        countonly=False,
        depth_limits=DEFAULT_DEPTH_LIMITS,
        mq_limits=DEFAULT_MQ_LIMITS,
        include_mNM_items=False,
    ):
        depth, mqlist = get_position_info(bam, vcfspec.chrom, vcfspec.pos0)
        depth_limits = cls.handle_limits_arg(depth_limits)
        mq_limits = cls.handle_limits_arg(mq_limits)

        if cls.check_position_validity(depth, mqlist, depth_limits, mq_limits):
            readstats_data = get_readstats_data(
                vcfspec, bam, fasta, chromdict,
                rpplist_kwargs=rpplist_kwargs,
                alleleinfo_kwargs=alleleinfo_kwargs,
                countonly=countonly,
            )
            result = cls.from_readstats_data(
                readstats_data, 
                vcfspec, 
                fasta, 
                chromdict,
                countonly=countonly,
                include_mNM_items=include_mNM_items,
            )
            del readstats_data
        else:
            result = cls.init_invalid(vcfspec, fasta, chromdict, countonly=countonly)

        return result

    @classmethod
    def from_readstats_data(
        cls, 
        readstats_data,  
        vcfspec, 
        fasta=None, 
        chromdict=None,
        countonly=False,
        include_mNM_items=False,
    ):
        def allele_means(data):
            return dict(
                (
                    alleleclass, 
                    np.nan if len(values) == 0 else np.mean(values),
                )
                for (alleleclass, values) in data.items()
            )

        def allele_medians(data):
            return dict(
                (
                    alleleclass, 
                    np.nan if len(values) == 0 else np.median(values),
                ) 
                for (alleleclass, values) in data.items()
            )

        def varpos_kstest(data):
            summary = dict()
            for alleleclass, values in data.items():
                if len(values) == 0:
                    pval = np.nan
                else:
                    pval = scipy.stats.kstest(values, ZERO_ONE_UNIFORM.cdf).pvalue
                summary[alleleclass] = pval

            return summary

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

        def handle_mNM(
            data, 
            rppcounts_dict, 
            rpp_range0_dict, 
            recurrence_cutoff_fraction,
            recurrence_cutoff_denominator,
        ):
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
                            (val / relevant_rppcount) >= recurrence_cutoff_fraction and
                            relevant_rppcount >= recurrence_cutoff_denominator
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

            return summary_without_recurrence, summary_only_recurrence, mNMitems_norecur, mNMitems_recur

        # main
        result = cls(is_missing=False)
        result.is_invalid = False
        result.vcfspec = vcfspec
        result.fasta = (
            refgenome.get_fasta(vcfspec.refver)
            if fasta is None else
            fasta
        )
        result.chromdict = (
            refgenome.get_chromdict(vcfspec.refver)
            if chromdict is None else
            chromdict
        )

        result['rppcounts'] = readstats_data['count'].copy()

        if not countonly:
            result['mean_BQs'] = allele_means(readstats_data['BQ'])
            result['median_BQs'] = allele_medians(readstats_data['BQ'])

            result['mean_MQs'] = allele_means(readstats_data['MQ'])
            result['median_MQs'] = allele_medians(readstats_data['MQ'])

            result['mean_cliplens'] = allele_means(readstats_data['cliplen'])
            result['median_cliplens'] = allele_medians(readstats_data['cliplen'])

            result['mNM'], result['recurrent_mNM'], mNMitems_norecur, mNMitems_recur = handle_mNM(
                readstats_data['mNM'], 
                readstats_data['count'], 
                readstats_data['rpp_range0'], 
                recurrence_cutoff_fraction=0.9,
                recurrence_cutoff_denominator=3,
            )
            if include_mNM_items:
                result['mNM_items'] = mNMitems_norecur
                result['recurrent_mNM_items'] = mNMitems_recur

            result['pairorient_pvalues'] = handle_orient(readstats_data['pairorient'], mode='pairorient')
            result['readorient_pvalues'] = handle_orient(readstats_data['readorient'], mode='pairorient')

            result['varpos_uniform_pvalues'] = varpos_kstest(readstats_data['pos0_left_fraction'])
            result['mean_varpos_fractions'] = allele_means(readstats_data['pos0_left_fraction'])
            result['median_varpos_fractions'] = allele_medians(readstats_data['pos0_left_fraction'])

        return result

    @classmethod
    def init_invalid(cls, vcfspec, fasta=None, chromdict=None, countonly=False, verbose=True):
        if verbose:
            logutils.print_timestamp(f'Initiating ReadStats object as invalid mode')
        # main
        result = cls(is_missing=False)
        result.is_invalid = True
        result.vcfspec = vcfspec
        result.fasta = (
            refgenome.get_fasta(vcfspec.refver)
            if fasta is None else
            fasta
        )
        result.chromdict = (
            refgenome.get_chromdict(vcfspec.refver)
            if chromdict is None else
            chromdict
        )

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

    #def write(self, vr, sampleid):
    #    return self.write_base(vr, sampleid)

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


class ReadStatsSampledict(annotitem.AnnotItemFormatSampledict):
    meta = {
        'ID': 'readstats', 'Number': 1, 'Type': 'String',
        'Description': 'Read-based statistics for the variant.',
    }
    unit_class = ReadStats

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
        fasta=None, 
        chromdict=None,
        rpplist_kwargs=dict(),
        alleleinfo_kwargs=dict(),
        countonly=False,
        depth_limits=DEFAULT_DEPTH_LIMITS,
        mq_limits=DEFAULT_MQ_LIMITS,
        init_invalid=False,
        include_mNM_items=False,
    ):
        depth_limits = cls.handle_limits_arg(depth_limits, bam_dict)
        mq_limits = cls.handle_limits_arg(mq_limits, bam_dict)

        result = cls()
        if init_invalid:
            for sampleid, bam in bam_dict.items():
                result[sampleid] = ReadStats.init_invalid(vcfspec, fasta, chromdict, countonly=countonly)
        else:
            for sampleid, bam in bam_dict.items():
                result[sampleid] = ReadStats.from_bam(
                    vcfspec, bam, fasta, chromdict,
                    rpplist_kwargs=rpplist_kwargs,
                    alleleinfo_kwargs=alleleinfo_kwargs,
                    countonly=countonly,
                    depth_limits=depth_limits[sampleid],
                    mq_limits=mq_limits[sampleid],
                    include_mNM_items=include_mNM_items,
                )
        return result

    @classmethod
    def from_vr(cls, vr, sampleid_list=None):
        return cls.from_vr_base(vr, sampleid_list=sampleid_list)

    def get_first_readstats(self):
        return next(iter(self.values()))

    @property
    def vcfspec(self):
        return self.get_first_readstats().vcfspec

    @property
    def fasta(self):
        return self.get_first_readstats().fasta

    @property
    def chromdict(self):
        return self.get_first_readstats().chromdict

    def update_bams(
        self, 
        bam_dict,
        rpplist_kwargs=dict(),
        alleleinfo_kwargs=dict(),
        countonly=False,
        depth_limits=DEFAULT_DEPTH_LIMITS,
        mq_limits=DEFAULT_MQ_LIMITS,
        init_invalid=False,
        include_mNM_items=False,
    ):
        if len(self) == 0:
            raise Exception(f'Length of self must be greater than 0')

        depth_limits = self.handle_limits_arg(depth_limits, bam_dict)
        mq_limits = self.handle_limits_arg(mq_limits, bam_dict)

        if init_invalid:
            for sampleid, bam in bam_dict.items():
                self[sampleid] = ReadStats.init_invalid(
                    self.vcfspec, self.fasta, self.chromdict, countonly=countonly,
                )
        else:
            for sampleid, bam in bam_dict.items():
                self[sampleid] = ReadStats.from_bam(
                    self.vcfspec, bam, self.fasta, self.chromdict,
                    rpplist_kwargs=rpplist_kwargs,
                    alleleinfo_kwargs=alleleinfo_kwargs,
                    countonly=countonly,
                    depth_limits=depth_limits[sampleid],
                    mq_limits=mq_limits[sampleid],
                    include_mNM_items=include_mNM_items,
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


def rpplist_to_readstats_data(
    rpplist, vcfspec, flanklen=liballeleinfo.DEFAULT_FLANKLEN,
):
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

    # initialize
    alleleclass_keys = vcfspec.get_alleleclasses()
    data = dict()

    # fields initialized as integer
    data['count'] = {x: 0 for x in alleleclass_keys}
    data['count']['softclip_overlap'] = 0

    # fields initialized as list
    for key in (
        'rpp_range0',
        'MQ', 'BQ', 'mNM', 'clipspec', 'cliplen', 'pairorient', 'readorient',
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
        if alleleclass_rpp is None:
            if rpp.check_softclip_overlaps_vcfspec(vcfspec):
                data['count']['softclip_overlap'] += 1

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

        # mNM
        mNM_data, clipspec_data = rpp.get_mNM_clipspec_data(vcfspec)
        data['mNM'][alleleclass_rpp].extend(mNM_data)
        data['clipspec'][alleleclass_rpp].extend(clipspec_data)

        # var_pos0s
        add_var_pos0s_rp(rpp.rp1, data, alleleclass_rpp, vcfspec)
        if rpp.rp2 is not None:
            add_var_pos0s_rp(rpp.rp2, data, alleleclass_rpp, vcfspec)

    # turn into numpy arrays
#    for key in data.keys():
#        if key != 'count':
#            data[key] = {
#                subkey: np.array(subval)
#                for subkey, subval in data[key].items()
#            }

    return data


def rpplist_to_readstats_data_countonly(
    rpplist, vcfspec, flanklen=liballeleinfo.DEFAULT_FLANKLEN,
):
    # initialize
    alleleclass_keys = vcfspec.get_alleleclasses()
    data = dict()
    data['count'] = {x: 0 for x in alleleclass_keys}
    data['count']['softclip_overlap'] = 0
    # add data
    for rpp in rpplist:
        if vcfspec not in rpp.alleleclass.keys():
            rpp.update_alleleclass(vcfspec=vcfspec, flanklen=flanklen)
        alleleclass_rpp = rpp.alleleclass[vcfspec]
        data['count'][alleleclass_rpp] += 1

    return data


def get_readstats_data(
    vcfspec, bam, 
    fasta=None, 
    chromdict=None, 
    rpplist_kwargs=dict(),
    alleleinfo_kwargs=dict(),
    countonly=False,
):
    """Only for non-sv cases"""

    #rpplist_kwargs.update({'view': False})
    rpplist_kwargs['view'] = False
    rpplist = readplus.get_rpplist_nonsv(
        bam=bam, 
        fasta=fasta, 
        chromdict=chromdict, 
        chrom=vcfspec.chrom, 
        start0=vcfspec.pos0, 
        end0=vcfspec.end0, 
        **rpplist_kwargs,
    )
    rpplist.update_alleleclass(
        vcfspec=vcfspec, 
        **alleleinfo_kwargs,
    )
    if countonly:
        readstats_data = rpplist_to_readstats_data_countonly(rpplist, vcfspec)
    else:
        readstats_data = rpplist_to_readstats_data(rpplist, vcfspec)

    del rpplist

    return readstats_data



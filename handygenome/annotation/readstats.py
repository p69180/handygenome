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

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
annotitem = importlib.import_module('.'.join([top_package_name, 'annotation', 'annotitem']))
infoformat = importlib.import_module('.'.join([top_package_name, 'variant', 'infoformat']))
readplus = importlib.import_module('.'.join([top_package_name, 'read', 'readplus']))
liballeleinfo = importlib.import_module('.'.join([top_package_name, 'read', 'alleleinfo']))


ZERO_ONE_UNIFORM = scipy.stats.uniform(loc=0, scale=1)


class ReadStats(annotitem.AnnotItemFormatSingle):
    @classmethod
    def from_bam(
        cls, 
        vcfspec, bam, fasta, chromdict,
        rpplist_kwargs=dict(),
        alleleinfo_kwargs=dict(),
        countonly=False,
    ):
        readstats_data = get_readstats_data(
            vcfspec, bam, fasta, chromdict,
            rpplist_kwargs=rpplist_kwargs,
            alleleinfo_kwargs=alleleinfo_kwargs,
            countonly=countonly,
        )
        return cls.from_readstats_data(readstats_data, countonly=countonly)

    @classmethod
    def from_readstats_data(cls, readstats_data, countonly=False):
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
            summary = dict()
            for alleleclass, values in data.items():    
                if len(values) == 0:
                    summary[alleleclass] = np.nan
                else:
                    mNM_count_sum = 0
                    counter = collections.Counter(values)
                    for key, val in counter.items():
                        refpos0 = key[0]
                        relevant_rppcount = sum(
                            any((refpos0 in x) for x in range0_list)
                            for range0_list in rpp_range0_dict[alleleclass]
                        )
                        if relevant_rppcount == 0:
                            raise Exception(f'The number of relevant rpp is 0. mNM element: {key}')
                        #assert relevant_rppcount > 0
                        #print(alleleclass, key, (val / relevant_rppcount))
                        if (
                            (val / relevant_rppcount) >= recurrence_cutoff_fraction and
                            relevant_rppcount >= recurrence_cutoff_denominator
                        ):
                            continue
                        else:
                            mNM_count_sum += 1
                    #print('mNM_count_sum', mNM_count_sum)
                            
                    summary[alleleclass] = mNM_count_sum / rppcounts_dict[alleleclass]

            return summary

        # main
        result = cls(is_missing=False)

        result['rppcounts'] = readstats_data['count'].copy()

        if not countonly:
            result['mean_BQs'] = allele_means(readstats_data['BQ'])
            result['median_BQs'] = allele_medians(readstats_data['BQ'])

            result['mean_MQs'] = allele_means(readstats_data['MQ'])
            result['median_MQs'] = allele_medians(readstats_data['MQ'])

            result['mean_cliplens'] = allele_means(readstats_data['cliplen'])
            result['median_cliplens'] = allele_medians(readstats_data['cliplen'])

            result['mNM'] = handle_mNM(
                readstats_data['mNM'], 
                readstats_data['count'], 
                readstats_data['rpp_range0'], 
                recurrence_cutoff_fraction=0.9,
                recurrence_cutoff_denominator=3,
            )

            result['pairorient_pvalues'] = handle_orient(readstats_data['pairorient'], mode='pairorient')
            result['readorient_pvalues'] = handle_orient(readstats_data['readorient'], mode='pairorient')

            result['varpos_uniform_pvalues'] = varpos_kstest(readstats_data['pos0_left_fraction'])
            result['mean_varpos_fractions'] = allele_means(readstats_data['pos0_left_fraction'])
            result['median_varpos_fractions'] = allele_medians(readstats_data['pos0_left_fraction'])

        return result

    #def write(self, vr, sampleid):
    #    return self.write_base(vr, sampleid)

    def get_allele_indexes_mean(self, key, allele_indexes):
        numerator = sum(
            self[key][x] * self['rppcounts'][x]
            for x in allele_indexes
        )
        denominator = sum(self['rppcounts'][x] for x in allele_indexes)

        if denominator == 0:
            return np.nan
        else:
            return numerator / denominator

    get_alleleindexes_mean = get_allele_indexes_mean

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
            return self['rppcounts'][alleleclass] / total_rppcount

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

    @classmethod
    def from_bam_dict(
        cls, bam_dict, vcfspec, fasta, chromdict,
        rpplist_kwargs=dict(),
        alleleinfo_kwargs=dict(),
        countonly=False,
    ):
        result = cls()
        for sampleid, bam in bam_dict.items():
            result[sampleid] = ReadStats.from_bam(
                vcfspec, bam, fasta, chromdict,
                rpplist_kwargs=rpplist_kwargs,
                alleleinfo_kwargs=alleleinfo_kwargs,
                countonly=countonly,
            )
        return result

    @classmethod
    def from_vr(cls, vr, sampleid_list=None):
        return cls.from_vr_base(vr, sampleid_list=sampleid_list)

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
    alleleclass_keys = (None,) + tuple(range(-1, len(vcfspec.alts) + 1))
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
    alleleclass_keys = (None,) + tuple(range(-1, len(vcfspec.alts) + 1))
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
    vcfspec, bam, fasta, chromdict, 
    rpplist_kwargs=dict(),
    alleleinfo_kwargs=dict(),
    countonly=False,
):
    """Only for non-sv cases"""

    rpplist_kwargs.update({'view': False})
    rpplist = readplus.get_rpplist_nonsv(
        bam=bam, 
        fasta=fasta, 
        chromdict=chromdict, 
        chrom=vcfspec.chrom, 
        start0=vcfspec.pos0, 
        end0=vcfspec.end0, 
        **rpplist_kwargs,
    )
    #rpplist.update_alleleinfo(
    #    vcfspec=vcfspec, 
    #    **alleleinfo_kwargs,
    #)
    rpplist.update_alleleclass(
        vcfspec=vcfspec, 
        **alleleinfo_kwargs,
    )
    if countonly:
        readstats_data = rpplist_to_readstats_data_countonly(rpplist, vcfspec)
    else:
        readstats_data = rpplist_to_readstats_data(rpplist, vcfspec)

    return readstats_data


#def summarize_readstats_data(readstats_data):
#    def ref_mean(data):
#        return np.mean(data[0])
#
#    def alt_means(data, alt_keys):
#        return [np.mean(data[x]) for x in alt_keys]
#
#    def allele_means(data):
#        return dict((alleleclass, np.mean(array)) 
#                    for (alleleclass, array) in data.items())
#
#    def allele_medians(data):
#        return dict((alleleclass, np.median(array)) 
#                    for (alleleclass, array) in data.items())
#
#    def varpos_kstest(data):
#        summary = dict()
#        for alleleclass, array in data.items():
#            if len(array) == 0:
#                pval = np.nan
#            else:
#                pval = scipy.stats.kstest(array, 'uniform').pvalue
#            summary[alleleclass] = pval
#
#        return summary
#
#    def setup_old(readstats, readstats_data):
#        alt_keys = sorted(x for x in readstats_data['count'].keys()
#                          if isinstance(x, int) and (x > 0))
#
#        readstats['ref_rppcount'] = readstats_data['count'][0]
#        readstats['alt_rppcounts'] = [readstats_data['count'][x] 
#                                      for x in alt_keys]
#        readstats['other_rppcount'] = readstats_data['count'][-1]
#        readstats['none_rppcount'] = readstats_data['count'][None]
#        readstats['clipoverlap_rppcount'] = readstats_data['count']['softclip_overlap']
#        
#        readstats['ref_mean_BQ'] = ref_mean(readstats_data['BQ'])
#        readstats['alt_mean_BQs'] = alt_means(readstats_data['BQ'], alt_keys)
#
#        readstats['ref_mean_MQ'] = ref_mean(readstats_data['MQ'])
#        readstats['alt_mean_MQs'] = alt_means(readstats_data['MQ'], alt_keys)
#
#        readstats['ref_mean_cliplen'] = ref_mean(readstats_data['cliplen'])
#        readstats['alt_mean_cliplens'] = alt_means(readstats_data['cliplen'], alt_keys)
#
#    def setup_new(readstats, readstats_data):
#        readstats['rppcounts'] = readstats_data['count'].copy()
#
#        readstats['mean_BQs'] = allele_means(readstats_data['BQ'])
#        readstats['median_BQs'] = allele_medians(readstats_data['BQ'])
#
#        readstats['mean_MQs'] = allele_means(readstats_data['MQ'])
#        readstats['median_MQs'] = allele_medians(readstats_data['MQ'])
#
#        readstats['mean_cliplens'] = allele_means(readstats_data['cliplen'])
#        readstats['median_cliplens'] = allele_medians(readstats_data['cliplen'])
#
#        readstats['varpos_uniform_pvalues'] = varpos_kstest(readstats_data['pos0_left_fraction'])
#        readstats['mean_varpos_fractions'] = allele_means(readstats_data['pos0_left_fraction'])
#        readstats['median_varpos_fractions'] = allele_medians(readstats_data['pos0_left_fraction'])
#
#    # main
#    readstats = ReadStats()
#    readstats.set_is_not_missing()
#    setup_new(readstats, readstats_data)
#
#    return readstats


#def add_meta(vcfheader):
#    vcfheader.add_meta(
#        key='FORMAT',
#        items=[('ID', READSTATS_FORMAT_KEY),
#               ('Type', 'String'),
#               ('Number', 1),
#               ('Description', 
#                'Read information statistics for the sample.')])
    


import re
import functools
import collections

import pysam
import numpy as np

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))

import handygenome.workflow as workflow
from handygenome.annotation.annotitem import AnnotItemInfoSingle, AnnotItemFormatSingle, AnnotItemFormatSampledict
import handygenome.variant.ponbams as libponbams


PON_LOGGER_INFO = workflow.get_debugging_logger(title='PONfilter', verbose=False)
PON_LOGGER_DEBUG = workflow.get_debugging_logger(title='PONfilter', verbose=True)
SUFFICIENT_PON_SAMPLE_NUM = 10


# preset filter combinations
def get_preset_filter_germline(
    with_pon=True, 
    pon_samples=None, 
    pon_cohorts=('BGI', 'PCAWG', 'SNULUNG'),
    refver='hg19',
    cutoff_diffbq=-5,
    cutoff_absbq=20,
    cutoff_diffmq=-10,
    cutoff_absmq=40,
    cutoff_clipovlp=2.0,
    cutoff_cliplen=20,
    cutoff_altcount=2,
    cutoff_otherratio=1.5,
    cutoff_totalcount=10,
    cutoff_unifpval=0.01,
):
    # sanity check
    if with_pon:
        if sum([
            pon_samples is not None,
            pon_cohorts is not None,
        ]) != 1:
            raise Exception(f'If "with_pon" is True, exactly 1 of "pon_samples" and "pon_cohorts" must be set.')

    result = FilterList(
        [
#            PopfreqFilter(
#                popnames=("GnomAD", "1000Genomes", "KOREAN", "Korea1K"), 
#                cutoff=0.01,
#            ),
            DiffMeanBQFilter(cutoff=cutoff_diffbq),
            AbsMeanBQFilter(cutoff=cutoff_absbq),
            DiffMeanMQFilter(cutoff=cutoff_diffmq),
            AbsMeanMQFilter(cutoff=cutoff_absmq),
            ClipoverlapFilter(cutoff=cutoff_clipovlp),
            CliplenFilter(cutoff=cutoff_cliplen),
            ReadcountFilter(cutoff=cutoff_altcount),
            OthercountRatioFilter(cutoff=cutoff_otherratio),
            TotaldepthGTFilter(cutoff=cutoff_totalcount),
            VarposUniformFilter(cutoff=cutoff_unifpval),
        ]
    )
    if with_pon:
        result.append(
            get_ponfilter(
                pon_samples=pon_samples, 
                pon_cohorts=pon_cohorts, 
                refver=refver, 
                mode='wgs',
            )
        )

    return result


def get_preset_filter_somatic(
    with_pon=True, 
    pon_samples=None, 
    pon_cohorts=('BGI', 'PCAWG', 'SNULUNG'),
    refver='hg19',
    cutoff_diffbq=-5,
    cutoff_absbq=20,
    cutoff_diffmq=-15,
    cutoff_absmq=40,
    cutoff_clipovlp=2.0,
    cutoff_altcount=2,
    cutoff_otherratio=1.5,
    cutoff_totalcount=10,
    cutoff_unifpval=0.05,
):
    # sanity check
    if with_pon:
        if sum([
            pon_samples is not None,
            pon_cohorts is not None,
        ]) != 1:
            raise Exception(f'If "with_pon" is True, exactly 1 of "pon_samples" and "pon_cohorts" must be set.')

    result = FilterList(
        [
#            PopfreqFilter(
#                popnames=("GnomAD", "1000Genomes", "KOREAN", "Korea1K"), 
#                cutoff=0.01,
#            ),
            DiffMeanBQFilter(cutoff=cutoff_diffbq),
            AbsMeanBQFilter(cutoff=cutoff_absbq),
            DiffMeanMQFilter(cutoff=cutoff_diffmq),
            AbsMeanMQFilter(cutoff=cutoff_absmq),
            ClipoverlapFilter(cutoff=cutoff_clipovlp),
            ReadcountFilter(cutoff=cutoff_altcount),
            OthercountRatioFilter(cutoff=cutoff_otherratio),
            TotaldepthGTFilter(cutoff=cutoff_totalcount),
            VarposUniformFilter(cutoff=cutoff_unifpval),
        ]
    )
    if with_pon:
        result.append(
            get_ponfilter(
                pon_samples=pon_samples, 
                pon_cohorts=pon_cohorts, 
                refver=refver, 
                mode='wgs',
            )
        )

    return result


def get_preset_filter_panelseq(
    with_pon=True, 
    pon_samples=None, 
    pon_cohorts=None,
    refver='hg19',
):
    # sanity check
    if with_pon:
        if sum([
            pon_samples is not None,
            pon_cohorts is not None,
        ]) != 1:
            raise Exception(f'If "with_pon" is True, exactly 1 of "pon_samples" and "pon_cohorts" must be set.')

    result = FilterList(
        [
            PopfreqFilter(
                popnames=("GnomAD", "1000Genomes", "KOREAN", "Korea1K"), 
                cutoff=0.01,
            ),
            DiffMeanBQFilter(cutoff=-5),
            AbsMeanBQFilter(cutoff=20),
            DiffMeanMQFilter(cutoff=-15),
            AbsMeanMQFilter(cutoff=40),
            ClipoverlapFilter(cutoff=1),
            ReadcountFilter(cutoff=2),
            OthercountRatioFilter(cutoff=1.5, ref_length_cutoff=4),
        ]
    )
    if with_pon:
        result.append(
            get_ponfilter(pon_samples=pon_samples, pon_cohorts=pon_cohorts, refver=refver, mode='panel')
        )

    return result


# Filter Result classes

class FilterResultBase:
    @property
    def passed(self):
        return self['passed']

    @classmethod
    def init_blank(cls):
        result = cls()
        result['is_passed'] = None
        result['data'] = list()
        return result


class FilterResultInfo(FilterResultBase, AnnotItemInfoSingle):
    meta = {
        'ID': 'filter_result_info', 'Number': 1, 'Type': 'String',
        'Description': 'Result of variant filtering process. Only one item for a variant, not considering each sample.',
    }

    @classmethod
    def from_vr(cls, vr):
        return cls.from_vr_base(vr)

    def write(self, vr):
        self.write_base(vr)


class FilterResultFormatUnit(FilterResultBase, AnnotItemFormatSingle):
    pass


class FilterResultFormat(AnnotItemFormatSampledict):
    meta = {
        'ID': 'filter_result_format', 'Number': 1, 'Type': 'String',
        'Description': 'Result of variant filtering process. Separate result for each sample.',
    }
    unit_class = FilterResultFormatUnit

    @classmethod
    def from_vr(cls, vr):
        return cls.from_vr_base(vr)

    def write(self, vr):
        self.write_base(vr)


###################################################
# base classes for Filter Objects and FilterList

class FilterBase:
    """"Return value of "check" method: 
        True means the variant is valid
    """
    def __repr__(self):
        return f"<{self.__class__.__name__} ({self.params})>"


class SamplewiseFilter(FilterBase):
    """"check" method:
        - call signature: def check(self, vp, sampleid=None, allele_index=1)
    """
    @classmethod
    def add_filter(cls, vp, sampleid):
        vp.add_sample_filter(sampleid, cls.__name__)

    def apply(self, vp, sampleid):
        mask = self.check(vp, sampleid)
        if not mask:
            self.add_filter(vp, sampleid)

    def check_show(self, vp, sampleid, allele_index=1):
        mask = self.check(vp, sampleid, allele_index)
        print(mask, self, sep=(" " * 4))

    def _check_meanvalue_difference(
        self, readstats_key, vp, sampleid, allele_index, cutoff
    ):
        readstats = vp.readstats_dict[sampleid]
        target_value = readstats[readstats_key][allele_index]

        #other_allele_indexes = vp.get_other_allele_indexes(allele_index)
        other_allele_indexes = tuple(sorted(
            set(range(len(vp.vr.alleles))).difference({allele_index})
        ))
        other_value = readstats.get_allele_indexes_average(
            readstats_key, other_allele_indexes
        )

        if np.isnan([target_value, other_value]).any():
            return True
        else:
            return (target_value - other_value) >= cutoff

    def _check_absvalue(self, readstats_key, vp, sampleid, allele_index, cutoff):
        readstats = vp.readstats_dict[sampleid]
        target_value = readstats[readstats_key][allele_index]

        if np.isnan(target_value):
            return True
        else:
            return target_value >= cutoff


class NonSamplewiseFilter(FilterBase):
    """"check" method:
        - call signature: def check(self, vp, allele_index=1)
    """
    def check_show(self, vp, allele_index=1):
        mask = self.check(vp, allele_index)
        print(mask, self, sep=(" " * 4))


class FilterList(list):
    def __repr__(self):
        repr_string = [repr(x) for x in self]
        return f"<FilterList object ({repr_string})>"

    def get_check_results(self, vp, sampleid, allele_index):
        check_results = list()
        for fi in self:
            if isinstance(fi, NonSamplewiseFilter):
                check_results.append(
                    fi.check(vp, allele_index=allele_index)
                )
            elif isinstance(fi, SamplewiseFilter):
                check_results.append(
                    fi.check(vp, sampleid=sampleid, allele_index=allele_index)
                )
        return check_results

    def check(self, vp, sampleid, allele_index=1, write_result=False):
        check_result = all(self.get_check_results(vp=vp, sampleid=sampleid, allele_index=allele_index))
        if write_result:
            vp.filter['is_passed'] = check_result
        return check_result

    def check_show(self, vp, sampleid, allele_index=1):
        for fi in self:
            if isinstance(fi, NonSamplewiseFilter):
                fi.check_show(vp=vp, allele_index=allele_index)
            elif isinstance(fi, SamplewiseFilter):
                fi.check_show(vp=vp, sampleid=sampleid, allele_index=allele_index)

    def check_any(self, vp, sampleid, allele_index=1):
        return any(self.get_check_results(vp=vp, sampleid=sampleid, allele_index=allele_index))


######################################


class PopfreqFilter(NonSamplewiseFilter):
    def __init__(
        self, 
        popnames=("GnomAD", "1000Genomes", "KOREAN", "Korea1K"), 
        cutoff=0.01,
    ):
        self.params = {
            "popnames": popnames,
            "cutoff": cutoff,
        }

    def check(self, vp, allele_index=1):
        assert allele_index > 0
        return all(
            vp.popfreq[allele_index - 1].get_freq(popname) <= self.params["cutoff"]
            for popname in self.params["popnames"]
        )


class DiffMeanBQFilter(SamplewiseFilter):
    def __init__(self, cutoff=-5):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, allele_index=1):
        return self._check_meanvalue_difference(
            readstats_key="mean_BQs",
            vp=vp,
            sampleid=sampleid,
            allele_index=allele_index,
            cutoff=self.params["cutoff"],
        )


class AbsMeanBQFilter(SamplewiseFilter):
    def __init__(self, cutoff=20):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, allele_index=1):
        return self._check_absvalue(
            readstats_key="mean_BQs",
            vp=vp,
            sampleid=sampleid,
            allele_index=allele_index,
            cutoff=self.params["cutoff"],
        )


class DiffMeanMQFilter(SamplewiseFilter):
    def __init__(self, cutoff=-15):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, allele_index=1):
        return self._check_meanvalue_difference(
            readstats_key="mean_MQs",
            vp=vp,
            sampleid=sampleid,
            allele_index=allele_index,
            cutoff=self.params["cutoff"],
        )


class AbsMeanMQFilter(SamplewiseFilter):
    def __init__(self, cutoff=40):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, allele_index=1):
        return self._check_absvalue(
            readstats_key="mean_MQs",
            vp=vp,
            sampleid=sampleid,
            allele_index=allele_index,
            cutoff=self.params["cutoff"],
        )


class ClipoverlapFilter(SamplewiseFilter):
    def __init__(self, cutoff=1):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, allele_index=1):
        readstats = vp.readstats_dict[sampleid]
        clipoverlap_count = readstats["rppcounts"]["softclip_overlap"]
        target_count = readstats["rppcounts"][allele_index]

        if clipoverlap_count == 0:
            return True
        else:
            ratio = target_count / clipoverlap_count
            return ratio > self.params["cutoff"]


class CliplenFilter(SamplewiseFilter):
    def __init__(self, cutoff=20):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, allele_index=None):
        readstats = vp.readstats_dict[sampleid]
        avg_cliplen = vp.readstats_dict[sampleid].get_allele_indexes_average(
            key='mean_cliplens', allele_indexes=None,
        )
        return avg_cliplen < self.params['cutoff']


class VarposUniformFilter(SamplewiseFilter):
    def __init__(self, cutoff=0.05):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, allele_index=1):
        readstats = vp.readstats_dict[sampleid]
        pval = readstats["varpos_uniform_pvalues"][allele_index]
        if np.isnan(pval):
            return True
        else:
            return pval >= self.params["cutoff"]


class VarposValueFilter(SamplewiseFilter):
    def __init__(self, cutoff=0.3):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, allele_index=1):
        return self._check_absvalue(
            readstats_key="mean_varpos_fractions",
            vp=vp,
            sampleid=sampleid,
            allele_index=allele_index,
            cutoff=self.params["cutoff"],
        )


class ReadcountFilter(SamplewiseFilter):
    def __init__(self, cutoff=2):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, allele_index=1):
        return self._check_absvalue(
            readstats_key="rppcounts",
            vp=vp,
            sampleid=sampleid,
            allele_index=allele_index,
            cutoff=self.params["cutoff"],
        )


class OthercountRatioFilter(SamplewiseFilter):
    def __init__(self, cutoff=1.5, ref_length_cutoff=None):
        self.params = {
            "cutoff": cutoff,
            'ref_length_cutoff': ref_length_cutoff,
        }

    def check(self, vp, sampleid, allele_index=1):
        if self.params['ref_length_cutoff'] is not None:
            if len(vp.vcfspec.ref) >= self.params['ref_length_cutoff']:
                return True

        readstats = vp.readstats_dict[sampleid]
        target_count = readstats["rppcounts"][allele_index]
        other_count = readstats["rppcounts"][-1]

        if other_count == 0:
            return True
        else:
            ratio = target_count / other_count
            return ratio > self.params["cutoff"]


class TotaldepthGTFilter(SamplewiseFilter):
    def __init__(self, cutoff=20):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, allele_index=None, exclude_other=False):
        return vp.readstats_dict[sampleid].get_total_rppcount(
            exclude_other=exclude_other,
        ) >= self.params['cutoff']


class TotaldepthLTFilter(SamplewiseFilter):
    def __init__(self, cutoff=200):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, allele_index=None, exclude_other=False):
        return vp.readstats_dict[sampleid].get_total_rppcount(
            exclude_other=exclude_other,
        ) <= self.params['cutoff']


############################

# PON filter

@common.get_deco_arg_choices({'mode': ('wgs', 'panel')})
@common.get_deco_num_set_differently(('pon_samples', 'pon_cohorts'), 1)
def get_ponfilter(pon_samples=None, pon_cohorts=None, refver=None, mode='wgs', **kwargs):
    if pon_cohorts is not None:
        pon_samples = libponbams.get_pon_sample_names(pon_cohorts, refver)
    else:
        pon_samples = list(pon_samples)

    if mode == 'panel':
        return PonFilterPanelseq(samples=pon_samples, **kwargs)
    elif mode == 'wgs':
        return PonFilterWGS(samples=pon_samples, **kwargs)


class PonFilterBase(SamplewiseFilter):
    default_cache = {
        'samples_with_low_depth': None,

        'lower_sample_fraction': None,
        'lower_mean': None,
        'lower_stdev': None,
        'lower_cutoff': None,
        'lower_different': None,

        'upper_sample_fraction': None,
        'upper_mean': None,
        'upper_stdev': None,
        'upper_cutoff': None,
        'upper_different': None,

        'is_germline': None,
        'het_vaf_mean': None,
        'germline_fraction_among_upper': None,
        'germline_fraction_among_all': None,

        'is_global': None,

        'upper_result': None,
        
        'total_result': None,
    }

    def __repr__(self):
        infostr = ", ".join(
            [
                f"samples={self.samples}",
                f"params={self.params}",
            ]
        )
        return f"<{self.__class__.__name__} ({infostr})>"

    @staticmethod
    def get_vaf_array(allele_index, readstats_dict, exclude_other):
        """Returns:
        np.ndarray composed of VAFs of PON samples, where 
            "excluded_pon_samples" are excluded and 
            nan values are removed
        Result is sorted in ascending order.
        """
        pon_vafs = list()
        for sampleid, readstats in readstats_dict.items():
#            readstats = readstats_dict[sampleid]
#            denom = readstats.get_total_rppcount()
#            if denom == 0:
#                vaf = np.nan
#            else:
#                numer = readstats["rppcounts"][allele_index]
#                vaf = numer / denom

            vaf = readstats_dict[sampleid].get_vaf(alleleclass=allele_index, exclude_other=exclude_other)
            if not np.isnan(vaf):
                pon_vafs.append(vaf)

        pon_vafs.sort()
        pon_vaf_array = np.array(pon_vafs)

        return pon_vaf_array

    def get_rppcount_array(self, allele_index, readstats_dict):
        return np.array([
            readstats["rppcounts"][allele_index]
            for readstats in readstats_dict.values()
        ])

    def get_used_pon_samples(self, query_sampleid, excluded_pon_samples, exclude_query):
        if excluded_pon_samples is None:
            excluded_pon_samples = set()
        else:
            excluded_pon_samples = set(excluded_pon_samples)

        if exclude_query:
            excluded_pon_samples.add(query_sampleid)

        used_pon_samples = set(self.samples).difference(excluded_pon_samples)
#        if len(used_pon_samples) < self.params["valid_sample_num_cutoff"]:
#            self.logger.warning(
#                f'The number of valid PON samples ({len(used_pon_samples)}) is '
#                f'less than the parameter "valid_sample_num_cutoff" ({self.params["valid_sample_num_cutoff"]}).'
#            )

        return used_pon_samples

    def set_common_params(self, vp, allele_index, query_sampleid, excluded_pon_samples, exclude_query, pon_depth_cutoff, exclude_other):
        # init used_pon_samples
        used_pon_samples = self.get_used_pon_samples(query_sampleid, excluded_pon_samples, exclude_query)
        if not used_pon_samples.issubset(set(vp.readstats_dict.keys())):
            raise Exception(f'Not all of PON sample IDs are included in the VariantPlus object.')

        # additional PON sample removal based on total read depth 
        pon_samples_with_low_depth = set()
        for sampleid in used_pon_samples:
            if vp.readstats_dict[sampleid].get_total_rppcount(exclude_other=exclude_other) < pon_depth_cutoff:
                pon_samples_with_low_depth.add(sampleid)
        used_pon_samples.difference_update(pon_samples_with_low_depth)
        if self.params['save_cache']:
            self.cache['samples_with_low_depth'] = pon_samples_with_low_depth

        # make readstats_dict
        readstats_dict = {
            sampleid: vp.readstats_dict[sampleid] 
            for sampleid in used_pon_samples
        }

        # others
        pon_vaf_array = self.get_vaf_array(allele_index, readstats_dict, exclude_other)  # nan is removed from "pon_vaf_array"
        query_vaf = vp.get_vaf(query_sampleid, allele_index, exclude_other=exclude_other)
        #pon_rppcount_array = self.get_rppcount_array(allele_index, readstats_dict)
        #query_rppcount = vp.readstats_dict[sampleid]["rppcounts"][allele_index]

        return used_pon_samples, readstats_dict, pon_vaf_array, query_vaf

    def reset_cache(self):
        self.cache = self.__class__.default_cache.copy()

    def bisect_vafs(self, pon_vaf_array):
        lower_selector = pon_vaf_array <= self.params['bisect_cutoff']
        return pon_vaf_array[lower_selector], pon_vaf_array[~lower_selector]

    def check_greater_than_lower(self, query_vaf, pon_vafs_lower, pon_vaf_array):
        sample_fraction = len(pon_vafs_lower) / len(pon_vaf_array)
        if sample_fraction <= self.params['lower_sample_fraction']:
            mean = None
            stdev = None
            cutoff = None
            is_different = True
        else:
            mean = np.mean(pon_vafs_lower)
            stdev = np.std(pon_vafs_lower)
            cutoff = min(
                0.4,  # 221231
                max(
                    mean + (self.params['lower_stdev_factor'] * stdev),
                    mean * self.params['lower_ratio']
                )
            )
            is_different = query_vaf >= cutoff

        if self.params['save_cache']:
            self.cache['lower_sample_fraction'] = sample_fraction
            self.cache['lower_mean'] = mean
            self.cache['lower_stdev'] = stdev
            self.cache['lower_cutoff'] = cutoff
            self.cache['lower_different'] = is_different

        return mean, stdev, cutoff, is_different

    def check_different_to_upper(self, query_vaf, pon_vafs_upper, pon_vaf_array):
        sample_fraction = len(pon_vafs_upper) / len(pon_vaf_array)
        if sample_fraction <= self.params['upper_sample_fraction']:
            mean = None
            stdev = None
            cutoff = None
            is_different = True
        else:
            mean = np.mean(pon_vafs_upper)
            stdev = np.std(pon_vafs_upper)
            cutoff = max(
                mean + (self.params['upper_stdev_factor'] * stdev),
                mean * self.params['upper_ratio']
            )
            is_different = query_vaf >= cutoff

        if self.params['save_cache']:
            self.cache['upper_sample_fraction'] = sample_fraction
            self.cache['upper_mean'] = mean
            self.cache['upper_stdev'] = stdev
            self.cache['upper_cutoff'] = cutoff
            self.cache['upper_different'] = is_different

        return mean, stdev, cutoff, is_different

    def check_global(self, pon_vafs_upper, pon_vaf_array):
        result = (len(pon_vafs_upper) / len(pon_vaf_array)) > self.params['global_sample_ratio_cutoff']
        if self.params['save_cache']:
            self.cache['is_global'] = result

        return result


class PonFilterPanelseqOld(PonFilterBase):
    @common.get_deco_arg_choices({"mode": ("mean", "max", "median")})
    def __init__(
        self,
        samples,
        verbose=False,

        mode="mean",
        #valid_sample_num_cutoff=20,
        min_germline_vaf=0.2,
        germline_sample_ratio_cutoff=0.1,
        # max_noise_vaf=0.1,
        lowest_subset_fraction=0.7,
        lowest_subset_num_cutoff=10,
        lowest_subset_snr_cutoff=5,
        nearby_ratio=1.3,
        #nearby_subset_fraction_cutoff=0.3,
        nearby_subset_num_cutoff=5,
        nearby_subset_snr_cutoff=3,

        pon_depth_cutoff=30,
    ):
        # set logger
        self.logger = (PON_LOGGER_DEBUG if verbose else PON_LOGGER_INFO)

        # sanity check
        if nearby_ratio < 1:
            raise Exception(
                f'"nearby_ratio" argument must be equal to or greater than 1.'
            )

        # main
        self.samples = samples
        self.params = {
            #"valid_sample_num_cutoff": valid_sample_num_cutoff,
            "min_germline_vaf": min_germline_vaf,
            "germline_sample_ratio_cutoff": germline_sample_ratio_cutoff,
            # "max_noise_vaf": max_noise_vaf,
            "lowest_subset_fraction": lowest_subset_fraction,
            "lowest_subset_num_cutoff": lowest_subset_num_cutoff,
            "lowest_subset_snr_cutoff": lowest_subset_snr_cutoff,
            "nearby_ratio": nearby_ratio,
            #"nearby_subset_fraction_cutoff": nearby_subset_fraction_cutoff,
            "nearby_subset_num_cutoff": nearby_subset_num_cutoff,
            "nearby_subset_snr_cutoff": nearby_subset_snr_cutoff,
            "mode": mode,
            "pon_depth_cutoff": pon_depth_cutoff,
        }

    def check(
        self,
        vp,
        sampleid,
        allele_index=1,
        excluded_pon_samples=None,
        exclude_query=True,
        do_germlinepat=True,
        do_lowest=True,
        do_nearby=True,
        exclude_other=False,
    ):
        # set params
        used_pon_samples, readstats_dict, pon_vaf_array, query_vaf = self.set_common_params(
            vp, allele_index, sampleid, excluded_pon_samples, exclude_query, self.params["pon_depth_cutoff"], exclude_other,
        )
        # main
        if np.isnan(query_vaf):
            return True

        subtests = list()
        if do_germlinepat:
            subtests.append(
                not self.check_germline_pattern(pon_vaf_array)
            )
        if do_lowest:
            subtests.append(
                # is_gt_noise = self.check_greater_than_lowest_count(query_rppcount, pon_rppcount_array)
                self.check_greater_than_lowest_vaf(query_vaf, pon_vaf_array)
            )
        if do_nearby:
            subtests.append(self.check_greater_than_nearby(query_vaf, pon_vaf_array))

        return all(subtests)

    ###################

    def get_lowest_summary_count(self, pon_rppcount_array):
        rppcounts_nonzero = pon_rppcount_array[pon_rppcount_array > 0]
        subset_num = int(len(rppcounts_nonzero) * self.params["subset_fraction"])

        if subset_num >= self.params["subset_num_cutoff"]:
            subset = sorted(rppcounts_nonzero)[:subset_num]
            noise_summary = getattr(np, self.params["mode"])(subset)
        else:
            noise_summary = None

        return noise_summary

    def check_greater_than_lowest_count(self, query_rppcount, pon_rppcount_array):
        noise_summary = self.get_lowest_summary_count(pon_rppcount_array)
        if noise_summary is None:
            return True
        else:
            snr = query_rppcount / noise_summary
            return snr >= self.params["lowest_subset_snr_cutoff"]

    def get_lowest_summary_vaf(self, pon_vaf_array):
        subset_num = int(pon_vaf_array.shape[0] * self.params["lowest_subset_fraction"])
        if subset_num == 0:
            lowest_summary = None
        else:
            subset = pon_vaf_array[:subset_num]
            lowest_summary = getattr(np, self.params["mode"])(subset)

        return lowest_summary

    def check_greater_than_lowest_vaf(self, query_vaf, pon_vaf_array):
        #assert not np.isnan(query_vaf)
        lowest_summary = self.get_lowest_summary_vaf(pon_vaf_array)
        self.logger.debug(f"lowest_summary {lowest_summary}")

        if (lowest_summary is None) or (lowest_summary == 0):
            return True
        else:
            snr = query_vaf / lowest_summary
            self.logger.debug(f"snr {snr}")

            return snr >= self.params["lowest_subset_snr_cutoff"]

    def check_greater_than_nearby(self, query_vaf, pon_vaf_array):
        #assert not np.isnan(query_vaf)
        upper_bound = query_vaf * self.params["nearby_ratio"]
        lower_bound = query_vaf / self.params["nearby_ratio"]
        vafs_nearby_query = pon_vaf_array[
            np.logical_and(pon_vaf_array <= upper_bound, pon_vaf_array >= lower_bound)
        ]
        self.logger.debug(f"upper_bound, {upper_bound}")
        self.logger.debug(f"lower_bound, {lower_bound}")
        self.logger.debug(f"vafs_nearby_query {vafs_nearby_query}")

        if len(vafs_nearby_query) == 0:
            return True
        else:
            #if len(vafs_nearby_query) >= len(pon_vafs) * self.params["nearby_subset_fraction_cutoff"]:
            if len(vafs_nearby_query) >= self.params["nearby_subset_num_cutoff"]:
                nearby_summary = getattr(np, self.params["mode"])(vafs_nearby_query)
                self.logger.debug(f"nearby_summary {nearby_summary}")

                if nearby_summary == 0:
                    return True
                else:
                    snr = query_vaf / nearby_summary
                    self.logger.debug(f'snr {snr}')
                    return snr >= self.params["nearby_subset_snr_cutoff"]
            else:
                return True

    def check_germline_pattern(self, pon_vaf_array):
        if pon_vaf_array.shape[0] == 0:
            self.logger.debug(f"germline_sample_ratio {np.nan}")
            return True
        else:
            n_germline_vaf = (pon_vaf_array >= self.params["min_germline_vaf"]).sum()
            germline_sample_ratio = n_germline_vaf / pon_vaf_array.shape[0]
            self.logger.debug(f"germline_sample_ratio {germline_sample_ratio}")
            return germline_sample_ratio >= self.params["germline_sample_ratio_cutoff"]


class PonFilterWGS(PonFilterBase):
    default_cache = PonFilterBase.default_cache 

    def __init__(
        self,
        samples,
        verbose=False,
        save_cache=False,

        check_global=True,

        bisect_cutoff=0.15,
        lower_stdev_factor=3,
        lower_ratio=3,
        lower_sample_fraction=0.2,
        upper_stdev_factor=3,
        upper_ratio=3,
        upper_sample_fraction=0.2,

        germline_het_range=(0.3, 0.7),
        germline_hom_range=(0.9, 1),
        het_vaf_mean_range=(0.4, 0.6),
        germline_fraction_among_upper_cutoff=0.8,
        germline_fraction_among_all_cutoff=0.1,

        global_sample_ratio_cutoff=0.9,

        pon_depth_cutoff=10,
    ):
        # set logger
        self.logger = (PON_LOGGER_DEBUG if verbose else PON_LOGGER_INFO)

        # sanity check

        # main
        self.samples = samples
        self.params = {
            'save_cache': save_cache,
            'check_global': check_global,

            'bisect_cutoff': bisect_cutoff,
            'lower_stdev_factor': lower_stdev_factor,
            'lower_ratio': lower_ratio,
            'lower_sample_fraction': lower_sample_fraction,
            'upper_stdev_factor': upper_stdev_factor,
            'upper_ratio': upper_ratio,
            'upper_sample_fraction': upper_sample_fraction,

            'germline_het_range': germline_het_range,
            'germline_hom_range': germline_hom_range,
            'het_vaf_mean_range': het_vaf_mean_range,
            'germline_fraction_among_upper_cutoff': germline_fraction_among_upper_cutoff,
            'germline_fraction_among_all_cutoff': germline_fraction_among_all_cutoff,

            'global_sample_ratio_cutoff': global_sample_ratio_cutoff,

            'pon_depth_cutoff': pon_depth_cutoff,
        }

    def check(
        self,
        vp,
        sampleid,
        allele_index=1,
        excluded_pon_samples=None,
        exclude_query=True,
        exclude_other=False,
    ):
        # init cache
        if self.params['save_cache']:
            self.reset_cache()
        # set params
        used_pon_samples, readstats_dict, pon_vaf_array, query_vaf = self.set_common_params(
            vp, allele_index, sampleid, excluded_pon_samples, exclude_query, self.params["pon_depth_cutoff"], exclude_other,
        )
        # main
        if np.isnan(query_vaf):
            if self.params['save_cache']:
                self.cache['total_result'] = True
            return True
        elif len(pon_vaf_array) == 0:
            if self.params['save_cache']:
                self.cache['total_result'] = True
            return True
        else:
            pon_vafs_lower, pon_vafs_upper = self.bisect_vafs(pon_vaf_array)

            lower_mean, lower_stdev, lower_cutoff, lower_different = self.check_greater_than_lower(query_vaf, pon_vafs_lower, pon_vaf_array)
            upper_mean, upper_stdev, upper_cutoff, upper_different = self.check_different_to_upper(query_vaf, pon_vafs_upper, pon_vaf_array)
            is_germline = self.check_is_germline(pon_vafs_upper, pon_vaf_array)
            is_global = self.check_global(pon_vafs_upper, pon_vaf_array)

            if self.params['check_global']:
                upper_result = self.get_upper_result(upper_different, is_germline, is_global)
            else:
                upper_result = self.get_upper_result_without_global(upper_different, is_germline)

            total_result = self.get_total_result(lower_different, upper_result)

            return total_result

    ###########################

    def check_is_germline(self, pon_vafs_upper, pon_vaf_array):
        if len(pon_vafs_upper) == 0:
            het_vaf_mean = None
            germline_fraction_among_upper = None
            germline_fraction_among_all = None
            result = False
        else:
            upper_vaf_is_het = np.logical_and(
                pon_vafs_upper >= self.params['germline_het_range'][0],
                pon_vafs_upper <= self.params['germline_het_range'][1],
            )
            upper_vaf_is_hom = np.logical_and(
                pon_vafs_upper >= self.params['germline_hom_range'][0],
                pon_vafs_upper <= self.params['germline_hom_range'][1],
            )
            upper_vaf_is_germline = np.logical_or(upper_vaf_is_het, upper_vaf_is_hom)

            if not upper_vaf_is_het.any():
                het_vaf_mean = np.nan
            else:
                het_vaf_mean = pon_vafs_upper[upper_vaf_is_het].mean()

            germline_upper_sample_count = upper_vaf_is_germline.sum()
            germline_fraction_among_upper = germline_upper_sample_count / len(pon_vafs_upper)
            germline_fraction_among_all = germline_upper_sample_count / len(pon_vaf_array)

            result = (
                (
                    np.isnan(het_vaf_mean)
                    or (
                        het_vaf_mean >= self.params['het_vaf_mean_range'][0]
                        and het_vaf_mean <= self.params['het_vaf_mean_range'][1]
                    )
                )
                and germline_fraction_among_upper >= self.params['germline_fraction_among_upper_cutoff']
                and germline_fraction_among_all >= self.params['germline_fraction_among_all_cutoff']
            )

        if self.params['save_cache']:
            self.cache['het_vaf_mean'] = het_vaf_mean
            self.cache['germline_fraction_among_upper'] = germline_fraction_among_upper
            self.cache['germline_fraction_among_all'] = germline_fraction_among_all
            self.cache['is_germline'] = result

        return result

    def get_upper_result(self, upper_different, is_germline, is_global):
        if is_global:
            result = False
        elif is_germline:
            result = True
        else:
            result = upper_different

        if self.params['save_cache']:
            self.cache['upper_result'] = result

        return result

    def get_upper_result_without_global(self, upper_different, is_germline):
        if is_germline:
            result = True
        else:
            result = upper_different

        if self.params['save_cache']:
            self.cache['upper_result'] = result

        return result

    def get_total_result(self, lower_different, upper_result):
        result = (lower_different and upper_result)
        if self.params['save_cache']:
            self.cache['total_result'] = result

        return result


class PonFilterPanelseqGermline(PonFilterBase):
    default_cache = PonFilterBase.default_cache 

    def __init__(
        self,
        samples,
        verbose=False,
        save_cache=False,

        check_global=True,

        bisect_cutoff=0.2,
        lower_stdev_factor=3,
        lower_ratio=3,
        lower_sample_fraction=0.2,
        upper_stdev_factor=3,
        upper_ratio=3,
        upper_sample_fraction=0.2,

        #germline_het_range=(0.4, 0.6),
        #germline_hom_range=(0.9, 1),
        global_sample_ratio_cutoff=0.9,

        pon_depth_cutoff=50,
    ):
        # set logger
        self.logger = (PON_LOGGER_DEBUG if verbose else PON_LOGGER_INFO)

        # sanity check

        # main
        self.samples = samples
        self.params = {
            'save_cache': save_cache,
            'check_global': check_global,

            'bisect_cutoff': bisect_cutoff,
            'lower_stdev_factor': lower_stdev_factor,
            'lower_ratio': lower_ratio,
            'lower_sample_fraction': lower_sample_fraction,
            'upper_stdev_factor': upper_stdev_factor,
            'upper_ratio': upper_ratio,
            'upper_sample_fraction': upper_sample_fraction,

            #'germline_het_range': germline_het_range,
            #'germline_hom_range': germline_hom_range,
            "global_sample_ratio_cutoff": global_sample_ratio_cutoff,

            "pon_depth_cutoff": pon_depth_cutoff,
        }

    def check(
        self,
        vp,
        sampleid,
        allele_index=1,
        excluded_pon_samples=None,
        exclude_query=True,
        exclude_other=False,
    ):
        # init cache
        if self.params['save_cache']:
            self.reset_cache()
        # set params
        used_pon_samples, readstats_dict, pon_vaf_array, query_vaf = self.set_common_params(
            vp, allele_index, sampleid, excluded_pon_samples, exclude_query, self.params["pon_depth_cutoff"], exclude_other,
        )
        # main
        if np.isnan(query_vaf):
            if self.params['save_cache']:
                self.cache['total_result'] = True
            return True
        elif len(pon_vaf_array) == 0:
            if self.params['save_cache']:
                self.cache['total_result'] = True
            return True
        else:
            pon_vafs_lower, pon_vafs_upper = self.bisect_vafs(pon_vaf_array)

            lower_mean, lower_stdev, lower_cutoff, lower_different = self.check_greater_than_lower(query_vaf, pon_vafs_lower, pon_vaf_array)
            is_global = self.check_global(pon_vafs_upper, pon_vaf_array)

            if self.params['check_global']:
                total_result = self.get_total_result(lower_different, is_global)
            else:
                total_result = self.get_total_result_without_global(lower_different)

            return total_result

    ###########################

    def get_total_result(self, lower_different, is_global):
        result = (lower_different and (not is_global))
        if self.params['save_cache']:
            self.cache['total_result'] = result

        return result

    def get_total_result_without_global(self, lower_different):
        result = lower_different
        if self.params['save_cache']:
            self.cache['total_result'] = result

        return result


class PonFilterPanelseqSomatic(PonFilterBase):
    default_cache = PonFilterBase.default_cache 

    def __init__(
        self,
        samples,
        verbose=False,
        save_cache=False,

        check_global=True,

        bisect_cutoff=0.2,
        lower_stdev_factor=3,
        lower_ratio=3,
        lower_sample_fraction=0.2,
        upper_stdev_factor=3,
        upper_ratio=3,
        upper_sample_fraction=0.2,

        #germline_het_range=(0.4, 0.6),
        #germline_hom_range=(0.9, 1),
        germline_sample_ratio_range=(0.2, 1),
        global_sample_ratio_cutoff=0.9,

        pon_depth_cutoff=50,
    ):
        # set logger
        self.logger = (PON_LOGGER_DEBUG if verbose else PON_LOGGER_INFO)

        # sanity check

        # main
        self.samples = samples
        self.params = {
            'save_cache': save_cache,
            'check_global': check_global,

            'bisect_cutoff': bisect_cutoff,
            'lower_stdev_factor': lower_stdev_factor,
            'lower_ratio': lower_ratio,
            'lower_sample_fraction': lower_sample_fraction,
            'upper_stdev_factor': upper_stdev_factor,
            'upper_ratio': upper_ratio,
            'upper_sample_fraction': upper_sample_fraction,

            #'germline_het_range': germline_het_range,
            #'germline_hom_range': germline_hom_range,
            'germline_sample_ratio_range': germline_sample_ratio_range,
            'global_sample_ratio_cutoff': global_sample_ratio_cutoff,

            "pon_depth_cutoff": pon_depth_cutoff,
        }

    def check(
        self,
        vp,
        sampleid,
        allele_index=1,
        excluded_pon_samples=None,
        exclude_query=True,
        exclude_other=False,
    ):
        # init cache
        if self.params['save_cache']:
            self.reset_cache()
        # set params
        used_pon_samples, readstats_dict, pon_vaf_array, query_vaf = self.set_common_params(
            vp, allele_index, sampleid, excluded_pon_samples, exclude_query, self.params["pon_depth_cutoff"], exclude_other,
        )
        # main
        if np.isnan(query_vaf):
            if self.params['save_cache']:
                self.cache['total_result'] = True
            return True
        elif len(pon_vaf_array) == 0:
            if self.params['save_cache']:
                self.cache['total_result'] = True
            return True
        else:
            pon_vafs_lower, pon_vafs_upper = self.bisect_vafs(pon_vaf_array)

            lower_mean, lower_stdev, lower_cutoff, lower_different = self.check_greater_than_lower(query_vaf, pon_vafs_lower, pon_vaf_array)
            is_germline = self.check_is_germline(pon_vafs_upper, pon_vaf_array)
            is_global = self.check_global(pon_vafs_upper, pon_vaf_array)

            if self.params['check_global']:
                total_result = self.get_total_result(lower_different, is_germline, is_global)
            else:
                total_result = self.get_total_result_without_global(lower_different, is_germline)

            return total_result

    ###########################

    def check_is_germline(self, pon_vafs_upper, pon_vaf_array):
        upper_sample_ratio = len(pon_vafs_upper) / len(pon_vaf_array)
        result = (
            upper_sample_ratio >= self.params['germline_sample_ratio_range'][0]
            and upper_sample_ratio <= self.params['germline_sample_ratio_range'][1]
        )
        if self.params['save_cache']:
            self.cache['is_germline'] = result

        return result

    def get_total_result(self, lower_different, is_germline, is_global):
        result = (lower_different and (not is_germline) and (not is_global))
        if self.params['save_cache']:
            self.cache['total_result'] = result

        return result

    def get_total_result_without_global(self, lower_different, is_germline):
        result = (lower_different and (not is_germline))
        if self.params['save_cache']:
            self.cache['total_result'] = result

        return result



############################



def get_transcript_subtype_filter(key):
    def vpfilter(vp):
        return any(
            feature["transcript_subtype_flags"][key]
            for feature in vp.annotdb.transcript_canon_ovlp.values()
        )

    return vpfilter


CODING_GENE_INVOLVED = get_transcript_subtype_filter("is_coding")


def get_consequence_filter(key):
    def vpfilter(vp):
        return any(
            feature["consequence_flags"][key]
            for feature in vp.annotdb.transcript_canon_ovlp.values()
        )

    return vpfilter


PROTEIN_ALTERING = get_consequence_filter("is_protein_altering")
NOT_PROTEIN_ALTERING = get_consequence_filter("is_not_protein_altering")


def get_genename_filter(gene_list, canonical_only=True):
    gene_set = set(gene_list)
    if canonical_only:

        def vpfilter(vp):
            vp_genes = set(vp.get_gene_names(canonical_only=True))
            return not vp_genes.isdisjoint(gene_set)

    else:

        def vpfilter(vp):
            vp_genes = set(vp.get_gene_names(canonical_only=False))
            return not vp_genes.isdisjoint(gene_set)

    return vpfilter


## function versions of filters
#
#def filter_popfreq(
#    vp, allele_index=1,
#    popnames=("GnomAD", "1000Genomes", "KOREAN", "Korea1K"), 
#    cutoff=0.01,
#):
#    return all(
#        vp.popfreq[allele_index - 1].get_freq(popname) <= cutoff
#        for popname in popnames
#    )
#
#
#def filter_diffMeanBQ(
#    vp, sampleid, allele_index=1,
#    cutoff=-5,
#):
#    return _check_meanvalue_difference(
#        readstats_key="mean_BQs",
#        vp=vp,
#        sampleid=sampleid,
#        allele_index=allele_index,
#        cutoff=cutoff,
#    )
#
#
#def filter_absMeanBQ(
#    vp, sampleid, allele_index=1,
#    cutoff=20,
#):
#    return _check_absvalue(
#        readstats_key="mean_BQs",
#        vp=vp,
#        sampleid=sampleid,
#        allele_index=allele_index,
#        cutoff=cutoff,
#    )
#
#
#def filter_diffMeanMQ(
#    vp, sampleid, allele_index=1,
#    cutoff=-15,
#):
#    return _check_meanvalue_difference(
#        readstats_key="mean_MQs",
#        vp=vp,
#        sampleid=sampleid,
#        allele_index=allele_index,
#        cutoff=cutoff,
#    )
#
#
#def filter_absMeanMQ(
#    vp, sampleid, allele_index=1,
#    cutoff=40,
#):
#    return _check_absvalue(
#        readstats_key="mean_MQs",
#        vp=vp,
#        sampleid=sampleid,
#        allele_index=allele_index,
#        cutoff=cutoff,
#    )
#
#
#def filter_clipoverlap(
#    vp, sampleid, allele_index=1,
#    cutoff=1,
#):
#    readstats = vp.readstats_dict[sampleid]
#    clipoverlap_count = readstats["rppcounts"]["softclip_overlap"]
#    target_count = readstats["rppcounts"][allele_index]
#
#    if clipoverlap_count == 0:
#        return True
#    else:
#        return (target_count / clipoverlap_count) > cutoff
#
#
#def filter_varReadCount(
#    vp, sampleid, allele_index=1,
#    cutoff=2,
#):
#    return self._check_absvalue(
#        readstats_key="rppcounts",
#        vp=vp,
#        sampleid=sampleid,
#        allele_index=allele_index,
#        cutoff=cutoff,
#    )
#
#
#def filter_otherReadRatio(
#    vp, sampleid, allele_index=1,
#    cutoff=1.5,
#):
#    readstats = vp.readstats_dict[sampleid]
#    target_count = readstats["rppcounts"][allele_index]
#    other_count = readstats["rppcounts"][-1]
#
#    if other_count == 0:
#        return True
#    else:
#        return (target_count / other_count) > cutoff
#
#
#def filter_varposUniform(
#    vp, sampleid, allele_index=1,
#    cutoff=0.05,
#):
#    readstats = vp.readstats_dict[sampleid]
#    pval = readstats["varpos_uniform_pvalues"][allele_index]
#    if np.isnan(pval):
#        return True
#    else:
#        return pval >= cutoff
#
#
#def filter_varposValue(
#    vp, sampleid, allele_index=1,
#    cutoff=0.3,
#):
#    return self._check_absvalue(
#        readstats_key="mean_varpos_fractions",
#        vp=vp,
#        sampleid=sampleid,
#        allele_index=allele_index,
#        cutoff=cutoff,
#    )
#
#
## filter function helpers
#
#def _check_meanvalue_difference(readstats_key, vp, sampleid, allele_index, cutoff, how='ge'):
#    assert how in ('ge', 'le')
#
#    readstats = vp.readstats_dict[sampleid]
#    #other_allele_indexes = vp.get_other_allele_indexes(allele_index)
#    other_allele_indexes = tuple(sorted(
#        set(range(len(vp.vr.alleles))).difference({allele_index})
#    ))
#
#    target_value = readstats[readstats_key][allele_index]
#    other_value = readstats.get_allele_indexes_mean(
#        readstats_key, other_allele_indexes
#    )
#
#    if np.isnan([target_value, other_value]).any():
#        return True
#    else:
#        if how == 'ge':
#            return (target_value - other_value) >= cutoff
#        elif how == 'le':
#            return (target_value - other_value) <= cutoff
#
#
#def _check_absvalue(readstats_key, vp, sampleid, allele_index, cutoff, how='ge'):
#    assert how in ('ge', 'le')
#
#    readstats = vp.readstats_dict[sampleid]
#    target_value = readstats[readstats_key][allele_index]
#
#    if np.isnan(target_value):
#        return True
#    else:
#        if how == 'ge':
#            return target_value >= cutoff
#        elif how == 'le':
#            return target_value <= cutoff



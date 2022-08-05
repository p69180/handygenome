import re
import functools

import pysam
import numpy as np

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))


FORMAT_FILTER_META = {
    "ID": "filter",
    "Type": "String",
    "Number": ".",
    "Description": (
        f"Works as FILTER column. The variant is true one if only 'PASS' exists."
    ),
}


def add_format_filter_meta(vcfheader):
    vcfheader.add_meta(key="FORMAT", items=FORMAT_FILTER_META.items())


###################################################


class VariantPlusFilter:
    def __repr__(self):
        return f"<{self.__class__.__name__} ({self.params})>"

    def check(self, vp, sampleid=None, alleleindex=1):
        pass

    def show_result(self, vp, sampleid=None, alleleindex=1):
        mask = self.check(vp, sampleid, alleleindex)
        print(mask, self, sep=(" " * 4))


class SamplewiseFilter(VariantPlusFilter):
    @classmethod
    def add_filter(cls, vp, sampleid):
        vp.add_sample_filter(sampleid, cls.__name__)

    def apply(self, vp, sampleid):
        mask = self.check(vp, sampleid)
        if not mask:
            self.add_filter(vp, sampleid)

    def _check_meanvalue_difference(
        self, readstats_key, vp, sampleid, alleleindex, cutoff
    ):
        readstats = vp.readstats_dict[sampleid]
        other_alleleindexes = vp.get_other_alleleindexes(alleleindex)

        target_value = readstats[readstats_key][alleleindex]
        other_value = readstats.get_alleleindexes_mean(
            readstats_key, other_alleleindexes
        )

        if np.isnan([target_value, other_value]).any():
            return True
        else:
            return (target_value - other_value) >= cutoff

    def _check_absvalue(self, readstats_key, vp, sampleid, alleleindex, cutoff):
        readstats = vp.readstats_dict[sampleid]
        target_value = readstats[readstats_key][alleleindex]

        if np.isnan(target_value):
            return True
        else:
            return target_value >= cutoff


class NonSamplewiseFilter(VariantPlusFilter):
    pass


class FilterList(list):
    def __init__(self, filters=list()):
        self.extend(filters)

    def __repr__(self):
        repr_string = [repr(x) for x in self]
        return f"<FilterList object ({repr_string})>"

    def check_all(self, vp, sampleid, alleleindex=1):
        return all(
            vf.check(vp, sampleid=sampleid, alleleindex=alleleindex) for vf in self
        )

    def check_any(self, vp, sampleid, alleleindex=1):
        return any(
            vf.check(vp, sampleid=sampleid, alleleindex=alleleindex) for vf in self
        )


#######################################


class PopfreqFilter(NonSamplewiseFilter):
    def __init__(
        self, popnames=("GnomAD", "1000Genomes", "KOREAN", "Korea1K"), cutoff=0.01
    ):
        self.params = {
            "popnames": popnames,
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid=None, alleleindex=1):
        return all(
            vp.get_popfreq(popname) <= self.params["cutoff"]
            for popname in self.params["popnames"]
        )


class PonFilter(SamplewiseFilter):
    @common.get_deco_arg_choices({"mode": ("mean", "max", "median")})
    def __init__(
        self,
        vcfspec,
        readstats_dict,
        valid_sample_num_cutoff=20,
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
        mode="mean",
    ):
        # sanity check
        if nearby_ratio < 1:
            raise Exception(
                f'"nearby_ratio" argument must be equal to or greater than 1.'
            )
        if len(readstats_dict) < valid_sample_num_cutoff:
            raise Exception(f'The number of samples in "readstats_dict" is less than "valid_sample_num_cutoff" argument. Please use a PON cohort with a larger number of samples, or decrease "valid_sample_num_cutoff" parameter.')

        # main
        self.vcfspec = vcfspec
        self.readstats_dict = readstats_dict
        self.sampleids = tuple(sorted(self.readstats_dict.keys()))
        self.nsample = len(self.readstats_dict)
        self.params = {
            "valid_sample_num_cutoff": valid_sample_num_cutoff,
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
        }

    def __repr__(self):
        infostr = ", ".join(
            [
                f"{self.vcfspec}",
                f"sample number: {self.nsample}",
                f"params: {self.params}",
            ]
        )
        return f"<PonFilter ({infostr})>"

    # @functools.cache
    def get_vaf_dict(self, alleleindex):
        vaf_dict = dict()
        for sampleid in self.sampleids:
            readstats = self.readstats_dict[sampleid]
            denom = readstats.get_total_rppcount()
            if denom == 0:
                vaf = np.nan
            else:
                numer = readstats["rppcounts"][alleleindex]
                vaf = numer / denom
            vaf_dict[sampleid] = vaf

        return vaf_dict

    @functools.cache
    def get_vafs(self, alleleindex, excluded_pon_samples, verbose):
        """Returns:
        np.ndarray composed of VAFs of PON samples, where samples in
            "excluded_pon_samples" and samples with nan vafs are excluded.
        Result is sorted in ascending order.
        """
        vaf_dict = self.get_vaf_dict(alleleindex)
        vafs = list()
        for sampleid in self.sampleids:
            vaf = vaf_dict[sampleid]
            if (sampleid not in excluded_pon_samples) and (not np.isnan(vaf)):
                vafs.append(vaf)

        if verbose:
            if len(vafs) < self.params["valid_sample_num_cutoff"]:
                print('The number of valid PON samples ({len(vafs)}) is less than the parameter "valid_sample_num_cutoff" ({self.params["valid_sample_num_cutoff"]}).')

        vafs.sort()
        vafs = np.array(vafs)

        return vafs

    # @functools.cache
    def get_rppcount_dict(self, alleleindex):
        rppcount_dict = dict()
        for sampleid in self.sampleids:
            readstats = self.readstats_dict[sampleid]
            rppcount_dict[sampleid] = readstats["rppcounts"][alleleindex]

        return rppcount_dict

    # @functools.cache
    def get_rppcounts(self, alleleindex):
        rppcount_dict = self.get_rppcount_dict(alleleindex)
        rppcounts = np.array([rppcount_dict[key] for key in self.sampleids])

        return rppcounts

    # @functools.cache
    def get_lowest_summary_count(self, alleleindex):
        rppcounts = self.get_rppcounts(alleleindex)
        rppcounts_nonzero = rppcounts[rppcounts > 0]
        subset_num = int(len(rppcounts_nonzero) * self.params["subset_fraction"])

        if subset_num >= self.params["subset_num_cutoff"]:
            subset = sorted(rppcounts_nonzero)[:subset_num]
            noise_summary = getattr(np, self.params["mode"])(subset)
        else:
            noise_summary = None

        return noise_summary

    def check_greater_than_lowest_count(self, vp, sampleid, alleleindex):
        noise_summary = self.get_lowest_summary_count(alleleindex)
        if noise_summary is None:
            return True
        else:
            rppcount = vp.readstats_dict[sampleid]["rppcounts"][alleleindex]
            snr = rppcount / noise_summary
            return snr >= self.params["lowest_subset_snr_cutoff"]

    # @functools.cache
    def get_lowest_summary_vaf(self, alleleindex, excluded_pon_samples, verbose):
        vafs = self.get_vafs(alleleindex, excluded_pon_samples, verbose)
        subset_num = int(len(vafs) * self.params["lowest_subset_fraction"])
        if subset_num == 0:
            lowest_summary = None
        else:
            subset = vafs[:subset_num]
            lowest_summary = getattr(np, self.params["mode"])(subset)

#            if subset_num >= self.params["lowest_subset_num_cutoff"]:
#                subset = vafs[:subset_num]
#                lowest_summary = getattr(np, self.params["mode"])(subset)
#            else:
#                lowest_summary = None

        return lowest_summary

    #    @functools.cache
    #    def get_noise_summary_vaf_deprecated(self, alleleindex, excluded_pon_samples):
    #        vafs = self.get_vafs(alleleindex, excluded_pon_samples, verbose)
    #        vafs_nonzero = vafs[vafs > 0]
    #        subset_num = int(len(vafs_nonzero) * self.params["lowest_subset_fraction"])
    #        if subset_num >= self.params["lowest_subset_num_cutoff"]:
    #            subset = vafs_nonzero[:subset_num]
    #            noise_summary = getattr(np, self.params["mode"])(subset)
    #        else:
    #            noise_summary = None
    #
    #        return noise_summary

    def check_greater_than_lowest_vaf(
        self, vp, sampleid, alleleindex, excluded_pon_samples, verbose
    ):
        query_vaf = vp.get_vaf(sampleid, alleleindex)
        if np.isnan(query_vaf):
            return True
        else:
            lowest_summary = self.get_lowest_summary_vaf(
                alleleindex, excluded_pon_samples, verbose
            )
            if verbose:
                print("lowest_summary", lowest_summary)

            if (lowest_summary is None) or (lowest_summary == 0):
                return True
            else:
                snr = query_vaf / lowest_summary
                if verbose:
                    print("snr", snr)

                return snr >= self.params["lowest_subset_snr_cutoff"]

    def check_greater_than_nearby(
        self, vp, sampleid, alleleindex, excluded_pon_samples, verbose
    ):
        query_vaf = vp.get_vaf(sampleid, alleleindex)
        if np.isnan(query_vaf):
            return True
        else:
            vafs = self.get_vafs(alleleindex, excluded_pon_samples, verbose)
            upper_bound = query_vaf * self.params["nearby_ratio"]
            lower_bound = query_vaf / self.params["nearby_ratio"]
            vafs_nearby_query = vafs[
                np.logical_and(vafs <= upper_bound, vafs >= lower_bound)
            ]
            if verbose:
                print("upper_bound", upper_bound)
                print("lower_bound", lower_bound)
                print("vafs_nearby_query", vafs_nearby_query)

            if len(vafs_nearby_query) == 0:
                return True
            else:
                #if len(vafs_nearby_query) >= len(vafs) * self.params["nearby_subset_fraction_cutoff"]:
                if len(vafs_nearby_query) >= self.params["nearby_subset_num_cutoff"]:
                    nearby_summary = getattr(np, self.params["mode"])(vafs_nearby_query)
                    if verbose:
                        print("nearby_summary", nearby_summary)

                    if nearby_summary == 0:
                        return True
                    else:
                        snr = query_vaf / nearby_summary
                        if verbose:
                            print("snr", snr)

                        return snr >= self.params["nearby_subset_snr_cutoff"]
                else:
                    return True

    # @functools.cache
    def check_germline_pattern(self, alleleindex, excluded_pon_samples, verbose):
        vafs = self.get_vafs(alleleindex, excluded_pon_samples, verbose)

        n_germline_vaf = (vafs >= self.params["min_germline_vaf"]).sum()
        germline_sample_ratio = n_germline_vaf / len(vafs)

        if verbose:
            print("germline_sample_ratio", germline_sample_ratio)
        return germline_sample_ratio >= self.params["germline_sample_ratio_cutoff"]

    def check(
        self,
        vp,
        sampleid,
        alleleindex=1,
        excluded_pon_samples=None,
        exclude_target=False,
        verbose=False,
        do_germlinepat=True,
        do_lowest=True,
        do_nearby=True,
    ):
        assert (
            self.vcfspec == vp.vcfspec
        ), f"Vcfspecs are different between PON and input VariantPlus object."
        # subtract self if required
        if excluded_pon_samples is None:
            excluded_pon_samples = set()
        else:
            excluded_pon_samples = set(excluded_pon_samples)

        if exclude_target:
            excluded_pon_samples.add(sampleid)

        # turn "excluded_pon_samples" argument into tuple for caching
        excluded_pon_samples = tuple(sorted(excluded_pon_samples))

        # main
        subtests = list()
        if do_germlinepat:
            subtests.append(
                not self.check_germline_pattern(
                    alleleindex,
                    excluded_pon_samples,
                    verbose,
                )
            )
        if do_lowest:
            subtests.append(
                # is_gt_noise = self.check_greater_than_lowest_count(vp, sampleid, alleleindex)
                self.check_greater_than_lowest_vaf(
                    vp,
                    sampleid,
                    alleleindex,
                    excluded_pon_samples,
                    verbose,
                )
            )
        if do_nearby:
            subtests.append(
                self.check_greater_than_nearby(
                    vp,
                    sampleid,
                    alleleindex,
                    excluded_pon_samples,
                    verbose,
                )
            )

        return all(subtests)


class DiffMeanBQFilter(SamplewiseFilter):
    def __init__(self, cutoff=-5):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        return self._check_meanvalue_difference(
            readstats_key="mean_BQs",
            vp=vp,
            sampleid=sampleid,
            alleleindex=alleleindex,
            cutoff=self.params["cutoff"],
        )


class AbsMeanBQFilter(SamplewiseFilter):
    def __init__(self, cutoff=20):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        return self._check_absvalue(
            readstats_key="mean_BQs",
            vp=vp,
            sampleid=sampleid,
            alleleindex=alleleindex,
            cutoff=self.params["cutoff"],
        )


class DiffMeanMQFilter(SamplewiseFilter):
    def __init__(self, cutoff=-15):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        return self._check_meanvalue_difference(
            readstats_key="mean_MQs",
            vp=vp,
            sampleid=sampleid,
            alleleindex=alleleindex,
            cutoff=self.params["cutoff"],
        )


class AbsMeanMQFilter(SamplewiseFilter):
    def __init__(self, cutoff=40):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        return self._check_absvalue(
            readstats_key="mean_MQs",
            vp=vp,
            sampleid=sampleid,
            alleleindex=alleleindex,
            cutoff=self.params["cutoff"],
        )


class ClipoverlapFilter(SamplewiseFilter):
    def __init__(self, cutoff=1):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        readstats = vp.readstats_dict[sampleid]
        clipoverlap_count = readstats["rppcounts"]["softclip_overlap"]
        target_count = readstats["rppcounts"][alleleindex]

        if clipoverlap_count == 0:
            return True
        else:
            ratio = target_count / clipoverlap_count
            return ratio > self.params["cutoff"]


class VarposUniformFilter(SamplewiseFilter):
    def __init__(self, cutoff=0.05):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        readstats = vp.readstats_dict[sampleid]
        pval = readstats["varpos_uniform_pvalues"][alleleindex]
        if np.isnan(pval):
            return True
        else:
            return pval >= self.params["cutoff"]


class VarposValueFilter(SamplewiseFilter):
    def __init__(self, cutoff=0.3):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        return self._check_absvalue(
            readstats_key="mean_varpos_fractions",
            vp=vp,
            sampleid=sampleid,
            alleleindex=alleleindex,
            cutoff=self.params["cutoff"],
        )


class ReadcountFilter(SamplewiseFilter):
    def __init__(self, cutoff=2):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        return self._check_absvalue(
            readstats_key="rppcounts",
            vp=vp,
            sampleid=sampleid,
            alleleindex=alleleindex,
            cutoff=self.params["cutoff"],
        )


class OthercountRatioFilter(SamplewiseFilter):
    def __init__(self, cutoff=1.5):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        readstats = vp.readstats_dict[sampleid]
        target_count = readstats["rppcounts"][alleleindex]
        other_count = readstats["rppcounts"][-1]

        if other_count == 0:
            return True
        else:
            ratio = target_count / other_count
            return ratio > self.params["cutoff"]


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

import re
import pprint
import os
import itertools
import random

import pysam
import pyranges as pr
import numpy as np
import pandas as pd

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
workflow = importlib.import_module(".".join([top_package_name, "workflow"]))

libbreakends = importlib.import_module(
    ".".join([top_package_name, "sv", "breakends"])
)

varianthandler = importlib.import_module(
    ".".join([top_package_name, "variantplus", "varianthandler"])
)
variantviz = importlib.import_module(
    ".".join([top_package_name, "variantplus", "variantviz"])
)
vpfilter = importlib.import_module(
    ".".join([top_package_name, "variantplus", "vpfilter"])
)
infoformat = importlib.import_module(
    ".".join([top_package_name, "variantplus", "infoformat"])
)
initvcf = importlib.import_module(".".join([top_package_name, "vcfeditor", "initvcf"]))

#annotationdb = importlib.import_module(
#    ".".join([top_package_name, "annotation", "annotationdb"])
#)
ensembl_feature = importlib.import_module(
    ".".join([top_package_name, "annotation", "ensembl_feature"])
)
ensembl_parser = importlib.import_module(
    ".".join([top_package_name, "annotation", "ensembl_parser"])
)
libpopfreq = importlib.import_module(
    ".".join([top_package_name, "annotation", "popfreq"])
)
libcosmic = importlib.import_module(
    ".".join([top_package_name, "annotation", "cosmic"])
)

readplus = importlib.import_module(".".join([top_package_name, "readplus", "readplus"]))
libreadstats = importlib.import_module(
    ".".join([top_package_name, "annotation", "readstats"])
)
alleleinfosetup = importlib.import_module(
    ".".join([top_package_name, "readplus", "alleleinfosetup"])
)
alleleinfosetup_sv = importlib.import_module(
    ".".join([top_package_name, "readplus", "alleleinfosetup_sv"])
)


READCOUNT_FORMAT_KEY = "allele_readcounts"


class VariantPlus:
    """Attributes:
        vr
        refver
        fasta
        chromdict
        vcfspec
        is_sv
        is_bnd1
        bnds
        annotdb
        annotdb_bnd1
        annotdb_bnd2
    """
    # constructors
    @classmethod
    def from_vr(cls, vr, refver=None, fasta=None, chromdict=None, **kwargs):
        result = cls()

        result.vr = vr
        result.refver = common.infer_refver_vr(result.vr) if (refver is None) else refver
        result.fasta = (
            pysam.FastaFile(common.DEFAULT_FASTA_PATHS[result.refver])
            if fasta is None
            else fasta
        )
        result.chromdict = (
            common.ChromDict(fasta=result.fasta) if chromdict is None else chromdict
        )
        result.vcfspec = common.Vcfspec.from_vr(result.vr)

        result.init_common(**kwargs)

        return result

    @classmethod
    def from_vcfspec(cls, vcfspec, refver, fasta=None, chromdict=None, **kwargs):
        result = cls()

        result.vcfspec = vcfspec
        result.refver = refver
        result.fasta = (
            pysam.FastaFile(common.DEFAULT_FASTA_PATHS[result.refver])
            if fasta is None
            else fasta
        )
        result.chromdict = (
            common.ChromDict(fasta=result.fasta) if chromdict is None else chromdict
        )
        result.vr = initvcf.create_vr(chromdict=result.chromdict, vcfspec=result.vcfspec)

        result.init_common(**kwargs)

        return result

    @classmethod
    def from_bnds(cls, bnds, refver, fasta=None, chromdict=None, **kwargs):
        vcfspec = bnds.get_vcfspec_bnd1()
        return cls.from_vcfspec(vcfspec, refver, fasta, chromdict, **kwargs)

    def init_common(
        self, 
        init_popfreq=True, 
        init_cosmic=True,
        init_transcript=True,
        init_regulatory=True,
        init_motif=False,
        init_repeat=True,
        init_readstats=True,
        popfreq_metadata=None, 
        cosmic_metadata=None,
    ):
        # SV attrs
        self.is_sv = varianthandler.check_SV(self.vr)
        self.set_bnds_attributes()  # is_bnd1, bnds

        # popfreq & cosmic
        if init_popfreq:
            self._init_popfreq(metadata=popfreq_metadata)
        if init_cosmic:
            self._init_cosmic(metadata=cosmic_metadata)

        # transcript, regulatory, motif, repeat
        if init_transcript:
            self._init_transcript()
        if init_regulatory:
            self._init_regulatory()
        if init_motif:
            self._init_motif()
        if init_repeat:
            self._init_repeat()

        # readstats_data_dict, rpplists
        self.rpplist_dict = dict()
        self.readstats_data_dict = dict()
        if init_readstats:
            self._init_readstats()

        # sampleids
        self.sampleids = tuple(self.vr.header.samples)

#    @common.get_deco_num_set(("vr", "vcfspec", "bnds"), 1)
#    def __init__(
#        self,
#        vr=None,
#        vcfspec=None,
#        bnds=None,
#        refver=None,
#        fasta=None,
#        chromdict=None,
#        #add_annotdb_infometas=True,
#        #set_annotdb=True,
#        init_popfreq=True,
#        init_cosmic=True,
#        init_transcript=True,
#        init_regulatory=True,
#        init_motif=False,
#        init_repeat=True,
#        init_readstats=True,
#        #annotitem_septype=annotationdb.DEFAULT_SEPTYPE,
#        popfreq_metadata=None, 
#        cosmic_metadata=None,
#    ):
#        """Args:
#            vr: pysam.VariantRecord instance
#            fasta: pysam.FastaFile instance
#            chromdict: handygenome.common.Chromdict instance
#        How to initialize:
#            1) "vr"
#            2) "vcfspec" and "refver"
#            3) "bnds" and "refver"
#        """
#
#        def init_from_vr(self, vr, refver, fasta, chromdict):
#            self.vr = vr
#            self.refver = common.infer_refver_vr(self.vr) if refver is None else refver
#            self.fasta = (
#                pysam.FastaFile(common.DEFAULT_FASTA_PATHS[self.refver])
#                if fasta is None
#                else fasta
#            )
#            self.chromdict = (
#                common.ChromDict(fasta=self.fasta) if chromdict is None else chromdict
#            )
#            self.vcfspec = varianthandler.get_vcfspec(self.vr)
#
#        def init_from_vcfspec(self, vcfspec, refver, fasta, chromdict):
#            if refver is None:
#                raise Exception(f"If initializing with vcfspec, refver must be given.")
#
#            self.vcfspec = vcfspec
#            self.refver = refver
#
#            self.fasta = (
#                pysam.FastaFile(common.DEFAULT_FASTA_PATHS[self.refver])
#                if fasta is None
#                else fasta
#            )
#            self.chromdict = (
#                common.ChromDict(fasta=self.fasta) if chromdict is None else chromdict
#            )
#
#            self.vr = initvcf.create_vr(chromdict=self.chromdict, vcfspec=self.vcfspec)
#
#        def init_from_bnds(self, bnds, refver, fasta, chromdict):
#            if refver is None:
#                raise Exception(
#                    f"If initializing with Breakends, refver must be given."
#                )
#
#            vcfspec = bnds.get_vcfspec_bnd1()
#            init_from_vcfspec(self, vcfspec, refver, fasta, chromdict)
#
#        # MAIN
#        if vr is not None:
#            init_from_vr(self, vr, refver, fasta, chromdict)
#        elif vcfspec is not None:
#            init_from_vcfspec(self, vcfspec, refver, fasta, chromdict)
#        elif bnds is not None:
#            init_from_bnds(self, bnds, refver, fasta, chromdict)
#
#        # SV attrs
#        self.is_sv = varianthandler.check_SV(self.vr)
#        self.set_bnds_attributes()  # is_bnd1, bnds
#
#        # annotdb
##        self.annotitem_septype = annotitem_septype
##        if add_annotdb_infometas:
##            annotationdb.add_infometas(self.vr.header)
##        if set_annotdb:
##            self.set_annotdb()
#
#        # popfreq & cosmic
#        if init_popfreq:
#            self._init_popfreq(metadata=popfreq_metadata)
#        if init_cosmic:
#            self._init_cosmic(metadata=cosmic_metadata)
#
#        # transcript, regulatory, motif, repeat
#        if init_transcript:
#            self._init_transcript()
#        if init_regulatory:
#            self._init_regulatory()
#        if init_motif:
#            self._init_motif()
#        if init_repeat:
#            self._init_repeat()
#
#        # readstats_data_dict, rpplists
#        self.rpplist_dict = dict()
#        self.readstats_data_dict = dict()
#        if init_readstats:
#            self._init_readstats()
#
#        # sampleids
#        self.sampleids = tuple(self.vr.header.samples)

    def __repr__(self):
        vr_string = "\t".join(str(self.vr).split("\t")[:5])
        return f"<VariantPlus object ({vr_string})>"

    # popfreq & cosmic
    def _init_popfreq(self, metadata=None):
        self.popfreq = libpopfreq.PopfreqInfoALTlist.from_vr(self.vr, metadata=metadata)

    def _init_cosmic(self, metadata=None):
        self.cosmic = libcosmic.CosmicInfoALTlist.from_vr(self.vr, metadata=metadata)

    def create_popfreq(self):
        if self.refver == 'GRCh37':
            refver = 'hg19'
        else:
            refver = self.refver

        dbsnp_vcf = libpopfreq.POPFREQ_VCFS[refver]
        popfreqlist = libpopfreq.PopfreqInfoALTlist.from_vcfspec(
            vcfspec=self.vcfspec, 
            dbsnp_vcf=dbsnp_vcf,
            fasta=self.fasta,
            metadata=libpopfreq.PopfreqMetadata.from_vcfheader(dbsnp_vcf.header),
        )
        return popfreqlist

    def create_cosmic(self):
        if self.refver == 'GRCh37':
            refver = 'hg19'
        else:
            refver = self.refver

        cosmic_vcf = libcosmic.COSMIC_VCFS[refver]
        cosmiclist = libcosmic.CosmicInfoALTlist.from_vcfspec(
            vcfspec=self.vcfspec, 
            cosmic_vcf=cosmic_vcf,
            fasta=self.fasta,
            metadata=libcosmic.CosmicMetadata.from_vcfheader(cosmic_vcf.header),
        )
        return cosmiclist

    # ensembl feature annotations
    def set_annotdb(self):
        if self.is_sv:
            self.annotdb = None
            self.annotdb_bnd1 = annotationdb.AnnotDB(
                "bnd1",
                self.refver,
                self.fasta,
                self.chromdict,
                vr=self.vr,
                septype=self.annotitem_septype,
            )
            self.annotdb_bnd2 = annotationdb.AnnotDB(
                "bnd2",
                self.refver,
                self.fasta,
                self.chromdict,
                vr=self.vr,
                septype=self.annotitem_septype,
            )
        else:
            self.annotdb = annotationdb.AnnotDB(
                "plain",
                self.refver,
                self.fasta,
                self.chromdict,
                vr=self.vr,
                septype=self.annotitem_septype,
            )
            self.annotdb_bnd1 = None
            self.annotdb_bnd2 = None

    def init_features(self):
        self._init_transcript()
        self._init_regulatory()
        self._init_motif()
        self._init_repeat()

    def _init_transcript(self):
        self.transcript = ensembl_feature.TranscriptSetALTlist.from_vr(self.vr)

    def _init_regulatory(self):
        self.regulatory = ensembl_feature.RegulatorySet.from_vr(self.vr)

    def _init_motif(self):
        self.motif = ensembl_feature.MotifSet.from_vr(self.vr)

    def _init_repeat(self):
        self.repeat = ensembl_feature.RepeatSet.from_vr(self.vr)

    def update_vep(self, vep_vr, alt_idx, distance, tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats):
        self.update_cmdline_vep(vep_vr, alt_idx, overwrite=True, create_new=True)
        self._update_post_vep_transcript(tabixfile_geneset, distance, alt_idx)
        self._update_post_vep_regulatory(tabixfile_regulatory, distance, alt_idx)
        self._update_post_vep_repeat(tabixfile_repeats, distance, alt_idx)

    def update_cmdline_vep(self, vep_vr, alt_idx, overwrite=True, create_new=True):
        parsed = ensembl_parser.parse_cmdline_vep(vep_vr)
        self.transcript[alt_idx].update_other(
            parsed["transcript"], overwrite=overwrite, create_new=create_new
        )
        self.regulatory.update_other(
            parsed["regulatory"], overwrite=overwrite, create_new=create_new
        )
        self.motif.update_other(
            parsed["motif"], overwrite=overwrite, create_new=create_new
        )

    def _update_post_vep_transcript(self, tabixfile_geneset, distance, alt_idx):
        self.transcript[alt_idx].update_ensembl_gff(
            vcfspec=self.vcfspec,
            chromdict=self.chromdict,
            distance=distance,
            tabixfile_geneset=tabixfile_geneset,
            refver=self.refver,
            fill_missing_canonical=True,
        )
        self.transcript[alt_idx].set_missing_distances(
            vcfspec=self.vcfspec, alt_idx=alt_idx
        )

    def _update_post_vep_regulatory(self, tabixfile_regulatory, distance, alt_idx):
        self.regulatory.update_ensembl_gff(
            vcfspec=self.vcfspec,
            chromdict=self.chromdict,
            distance=distance,
            tabixfile_regulatory=tabixfile_regulatory,
        )
        self.regulatory.set_missing_distances(vcfspec=self.vcfspec, alt_idx=alt_idx)

    def _update_post_vep_repeat(self, tabixfile_repeats, distance, alt_idx):
        self.repeat.update_repeatmasker_bed(
            vcfspec=self.vcfspec,
            chromdict=self.chromdict,
            distance=distance,
            tabixfile_repeats=tabixfile_repeats,
        )
        self.repeat.set_missing_distances(vcfspec=self.vcfspec, alt_idx=alt_idx)

    # alleleindexes
    def get_other_alleleindexes(self, allele_index):
        """Returns all integer allele indexes other than the input"""

        all_alleleindexes = set(range(len(self.vr.alleles)))
        other_alleleindexes = all_alleleindexes.difference([allele_index])
        return tuple(sorted(other_alleleindexes))

    # miscellaneous
    def get_gr(self):
        return pr.from_dict(
            {
                "Chromosome": [self.vr.contig],
                "Start": [self.vr.pos - 1],
                "End": [self.vr.pos],
            }
        )

    def get_gene_names(self, canonical_only=True):
        if canonical_only:
            return [
                feature["gene_name"]
                for feature in self.annotdb.transcript_canon_ovlp.values()
            ]
        else:
            return [
                feature["gene_name"]
                for feature in self.annotdb.transcript_ovlp.values()
            ]

    def check_intergenic(self):
        return len(self.annotdb.transcript) == 0

    def get_info(self, key, collapse_tuple=True):
        return infoformat.get_value_info(self.vr, key, collapse_tuple=collapse_tuple)

    def get_format(self, sampleid, key, collapse_tuple=True):
        return infoformat.get_value_format(
            self.vr, sampleid, key, collapse_tuple=collapse_tuple
        )

    def set_info(self, key, val, typeconv=True):
        infoformat.set_info(self.vr, key, val, typeconv=typeconv)

    def set_format(self, sampleid, key, val, typeconv=True):
        infoformat.set_format(self.vr, sampleid, key, val, typeconv=typeconv)

    def check_NA_info(self, key):
        return infoformat.check_NA_info(self.vr, key)

    def check_NA_format(self, sampleid, key):
        return infoformat.check_NA_format(self.vr, sampleid, key)

    def show_info(self):
        infoformat.show_info(self.vr)

    def show_format(self):
        infoformat.show_format(self.vr)

    def _refine_vr_InfoFormatValues(self):
        infoformat.refine_vr_InfoFormatValue(self.vr)

    # breakends
    def get_vr_bnd2(self):
        assert self.bnds is not None, f'"bnds" attribute must be set.'

        vr_bnd2 = self.vr.header.new_record()
        vr_bnd2.id = self.bnds.get_id_bnd2()
        vcfspec_bnd2 = self.bnds.get_vcfspec_bnd2()
        varianthandler.apply_vcfspec(vr_bnd2, vcfspec_bnd2)
        vr_bnd2.info["MATEID"] = self.vr.id

        return vr_bnd2

    def set_bnds_attributes(self):
        if self.is_sv:
            vr_svinfo = libbreakends.get_vr_svinfo_standard_vr(
                self.vr,
                self.fasta,
                self.chromdict,
            )
            self.is_bnd1 = vr_svinfo["is_bnd1"]
            self.bnds = libbreakends.get_bnds_from_vr_svinfo(
                self.vr, vr_svinfo, self.fasta, self.chromdict
            )
        else:
            self.is_bnd1 = None
            self.bnds = None

    # rpplist
    def make_rpplist(
        self,
        bam,
        set_alleleinfo=True,
        rpplist_kwargs={
            'view': False,
            'no_matesearch': True,
        },
        alleleinfo_kwargs=dict(),
    ):
        if self.is_sv:
            rpplist = readplus.get_rpplist_sv(
                bam=bam,
                fasta=self.fasta,
                chromdict=self.chromdict,
                bnds=self.bnds,
                **rpplist_kwargs,
            )
            if set_alleleinfo:
                rpplist.update_alleleinfo_sv(
                    bnds=self.bnds,
                    **alleleinfo_kwargs,
                )
        else:
            rpplist = readplus.get_rpplist_nonsv(
                bam=bam,
                fasta=self.fasta,
                chromdict=self.chromdict,
                chrom=self.vcfspec.chrom,
                start0=self.vcfspec.pos0,
                end0=self.vcfspec.end0,
                **rpplist_kwargs,
            )
            if set_alleleinfo:
                rpplist.update_alleleinfo(
                    vcfspec=self.vcfspec, 
                    **alleleinfo_kwargs,
                )

        return rpplist

    def set_rpplist(
        self, sampleid, bam, 
        set_alleleinfo=True,
        rpplist_kwargs={
            'view': False,
            'no_matesearch': True,
        },
        alleleinfo_kwargs=dict(),
    ):
        self.rpplist_dict[sampleid] = self.make_rpplist(
            bam, set_alleleinfo=set_alleleinfo,
            rpplist_kwargs=rpplist_kwargs,
            alleleinfo_kwargs=alleleinfo_kwargs,
        )

#    @common.get_deco_num_set(('bam', 'sampleid'), 1)
#    def get_rpplist(self, sampleid=None, bam=None):
#        if sampleid is not None:
#            if sampleid not in self.rpplist_dict:
#                raise Exception(f"rpplist for the sampleid {sampleid} is not created. To make it, run VariantPlus.set_rpplist method.")
#            return self.rpplist_dict[sampleid]
#        elif bam is not None:
#            return self.make_rpplist(bam=bam)

    # readstats_data
    def make_readstats_data(self, bam, **kwargs):
        return libreadstats.get_readstats_data(
            self.vcfspec,
            bam,
            self.fasta,
            self.chromdict,
            **kwargs,
        )

    def set_readstats_data(
        self, 
        sampleid, 
        bam, 
        verbose=False, 
        rpplist_kwargs={
            'view': False,
            'no_matesearch': True,
        },
        alleleinfo_kwargs=dict(),
    ):
        if verbose:
            print(f"Getting readstats_data for the sample {sampleid}")

        if sampleid not in self.rpplist_dict.keys():
            self.set_rpplist(
                sampleid, 
                bam, 
                set_alleleinfo=True,
                rpplist_kwargs=rpplist_kwargs,
                alleleinfo_kwargs=alleleinfo_kwargs,
            )

        readstats_data = libreadstats.rpplist_to_readstats_data(
            self.rpplist_dict[sampleid],
            self.vcfspec,
        )
        self.readstats_data_dict[sampleid] = readstats_data

#    def get_readstats_data(self, sampleid, bam=None):
#        if sampleid not in self.readstats_data_dict:
#            if bam is None:
#                raise Exception(
#                    f"readstats_data for the sampleid {sampleid} is not "
#                    f"prepared. To make it, bam is required."
#                )
#            self.set_readstats_data(sampleid, bam, verbose=False)
#        return self.readstats_data_dict[sampleid]

    # readstats
    def _init_readstats(self):
        self.readstats_dict = dict()
        for sampleid in self.vr.samples.keys():
            if self.check_NA_format(sampleid, libreadstats.ReadStats.meta["ID"]):
                self.readstats_dict[sampleid] = None
            else:
                readstats = libreadstats.ReadStats.from_vr(self.vr, sampleid)
                readstats.postprocess()
                self.readstats_dict[sampleid] = readstats

    def make_readstats(self, bam, **kwargs):
        return libreadstats.get_readstats(self.vcfspec, bam, self.fasta, self.chromdict, **kwargs)

    def set_readstats(self, sampleid, bam):
        if sampleid not in self.readstats_data_dict.keys():
            self.set_readstats_data(sampleid, bam)

        readstats_data = self.readstats_data_dict[sampleid]
        self.readstats_dict[sampleid] = libreadstats.summarize_readstats_data(readstats_data)

    # readstats-related-others
    def calc_readcounts(self, bam, no_matesearch=True):
        """Designed only for non-sv"""

        rpplist = self.make_rpplist(
            bam, no_matesearch=no_matesearch, set_alleleinfo=True
        )
        readcounts = rpplist.get_readcounts(self.vcfspec)

        return readcounts

    def get_ponfilter(self, sampleids=None, **kwargs):
        if sampleids is None:
            readstats_dict = self.readstats_dict
        else:
            readstats_dict = {
                sampleid: self.readstats_dict[sampleid] for sampleid in sampleids
            }

        return vpfilter.PonFilter(self.vcfspec, readstats_dict, **kwargs)

    def get_total_rppcount(self, sampleid):
        """Sum of rpp counts for all alleleclasses except None"""
        readstats = self.readstats_dict[sampleid]
        return readstats.get_total_rppcount()

    def get_vaf(self, sampleid, allele_index=1, ndigits=None):
        readstats = self.readstats_dict[sampleid]
        vaf = readstats.get_vaf(allele_index)
        if np.isnan(vaf):
            return vaf
        else:
            if ndigits is None:
                return vaf
            else:
                return round(vaf, ndigits)

    # filters
    def reset_sample_filter(self, sampleid):
        self.set_format(sampleid, vpfilter.FORMAT_FILTER_META["ID"], None)

    def reset_sample_filter_all(self):
        for sampleid in self.vr.header.samples:
            self.reset_sample_filter(sampleid)

    def check_sample_filter(self, sampleid):
        sample_filter = self.get_format(
            sampleid, vpfilter.FORMAT_FILTER_META["ID"], collapse_tuple=False
        )
        return sample_filter == ("PASS",)

    def add_sample_filter(self, sampleid, value: str):
        if value == "PASS":
            self.set_format(sampleid, vpfilter.FORMAT_FILTER_META["ID"], ("PASS",))
        else:
            if self.check_NA_format(sampleid, vpfilter.FORMAT_FILTER_META["ID"]):
                new_val = (value,)
            else:
                old_val = self.get_format(
                    sampleid, vpfilter.FORMAT_FILTER_META["ID"], collapse_tuple=False
                )
                new_val = tuple(set(old_val + (value,)))

            self.set_format(sampleid, vpfilter.FORMAT_FILTER_META["ID"], new_val)

    def show_sample_filters(self):
        for sampleid in self.vr.header.samples:
            print(
                sampleid,
                self.get_format(
                    sampleid, vpfilter.FORMAT_FILTER_META["ID"], collapse_tuple=False
                ),
            )

    # visualizations
    @common.get_deco_num_set(('bam', 'sampleid'), 1)
    def show_igv(
        self, igv, bam=None, sampleid=None, tmpbam_dir=None,
        rpplist_kwargs={
            'view': True,
            'no_matesearch': True,
        },
        alleleinfo_kwargs=dict(),
    ):
        # get rpplist
        if sampleid is not None:
            if sampleid not in self.rpplist_dict.keys():
                raise Exception(f"rpplist for the sampleid {sampleid} is not created. To bypass this problem, 1) run this method with 'bam' argument, or 2) run VariantPlus.set_rpplist method first.")
            rpplist = self.rpplist_dict[sampleid]
        elif bam is not None:
            rpplist = self.make_rpplist(
                bam, set_alleleinfo=True,
                rpplist_kwargs=rpplist_kwargs,
                alleleinfo_kwargs=alleleinfo_kwargs,
            )

        # add alleleinfo tag to ReadPlus and ReadPlusPair objects
        if self.is_sv:
            rpplist.set_alleleinfo_tag_rp_sv(self.bnds)
            rpplist.set_alleleinfo_tag_rpp_sv(self.bnds)
        else:
            rpplist.set_alleleinfo_tag_rp(self.vcfspec)
            rpplist.set_alleleinfo_tag_rpp(self.vcfspec)

        # set tmpbam_path
        if tmpbam_dir is None:
            tmpbam_dir = os.getcwd()

        if sampleid is None:
            prefix = None
        else:
            prefix = f'{sampleid}_for_show_'
        tmpbam_path = workflow.get_tmpfile_path(prefix=prefix, suffix=".bam", where=tmpbam_dir)

        # main
        rpplist.write_bam(outfile_path=tmpbam_path)
        igv.load([tmpbam_path])

        if self.is_sv:
            (reads_range0_bnd1, reads_range0_bnd2) = rpplist.get_ref_range0_sv(
                self.bnds
            )
            pos_range0_bnd1 = self.bnds.get_pos_range0_bnd1()
            pos_range0_bnd2 = self.bnds.get_pos_range0_bnd2()

            width_bnd1 = max(
                pos_range0_bnd1.start - reads_range0_bnd1.start,
                reads_range0_bnd1.stop - pos_range0_bnd1.stop,
            )
            width_bnd2 = max(
                pos_range0_bnd2.start - reads_range0_bnd2.start,
                reads_range0_bnd2.stop - pos_range0_bnd2.stop,
            )
            width = max(width_bnd1, width_bnd2)

            igv.goto([self.bnds], width=width)
        else:
            reads_ref_range0 = rpplist.get_ref_range0(self.vcfspec.chrom)
            width = max(
                self.vcfspec.pos0 - reads_ref_range0.start,
                reads_ref_range0.stop - self.vcfspec.pos0,
            )

            igv.goto([self.vcfspec], width=width)

        igv.cmd("group")
        igv.viewaspairs()
        igv.cmd(f"colorBy TAG {readplus.ALLELEINFO_TAG_RPP}")

        os.remove(tmpbam_path)

    def show_readcounts(self, sampleids=None, **kwargs):
        if sampleids is None:
            sampleids = sorted(self.vr.samples.keys())

        variantviz.show_readcounts(
            self, sampleid_order=sampleids, title=str(self.vcfspec), **kwargs
        )

#    @common.get_deco_num_set(("sampleid", "bam"), 1)
#    def show_readstats_data(self, sampleid=None, bam=None, varpos_key="left"):
#        if sampleid is None:
#            readstats_data = libreadstats.get_readstats_data(
#                self.vcfspec, bam, self.fasta, self.chromdict
#            )
#        else:
#            readstats_data = self.get_readstats_data(sampleid, bam=bam)
#
#        variantviz.show_readstats_data(
#            readstats_data, title=sampleid, varpos_key=varpos_key
#        )

    ##################################################################

    # others
    def get_popfreq(self, popname):
        val = self.annotdb.popfreq[popname]
        if val is None:
            val = 0

        return val

    def get_cosmic_total_occur(self):
        val = self.annotdb.cosmic["total_occurrence"]
        if val is None:
            val = 0

        return val


class VariantPlusList(list):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._gr = None

    def sample(self, n=1):
        if n == 1:
            return random.choice(self)
        else:
            return random.sample(self, k=n)

    def get_df(self, vaf_sampleid=None, as_gr=False):
        if vaf_sampleid is None:
            sample_ids = list(self[0].vr.header.samples)
            if len(sample_ids) == 0:
                vaf_sampleid = None
            else:
                vaf_sampleid = sample_ids[0]

        chroms = list()
        poss = list()
        starts = list()
        ends = list()
        refs = list()
        alts = list()
        vafs = list()

        for vp in self:
            chroms.append(vp.vcfspec.chrom)
            poss.append(vp.vcfspec.pos)
            starts.append(vp.vcfspec.pos0)
            ends.append(vp.vcfspec.end0)
            refs.append(vp.vcfspec.ref)
            alts.append(vp.vcfspec.alts[0])

            if vaf_sampleid is None:
                vafs.append(np.nan)
            else:
                vaf = vp.get_vaf(vaf_sampleid, allele_index=1, ndigits=None)
                vafs.append(vaf)  # vaf can be np.nan

        if as_gr:
            return pr.from_dict(
                {
                    'Chromosome': chroms,
                    'Start': starts,
                    'End': ends,
                    'REF': refs,
                    'ALT': alts,
                    'VAF': vafs,
                },
                int64=False,
            )
        else:
            return pd.DataFrame.from_dict(
                {
                    'CHROM': chroms,
                    'POS': poss,
                    'REF': refs,
                    'ALT': alts,
                    'VAF': vafs,
                }
            )

    def get_gr(self, vaf_sampleid=None):
        return self.get_df(vaf_sampleid=vaf_sampleid, as_gr=True)

    def get_isec(self, gr):
        overlaps = self.get_gr().count_overlaps(gr, overlap_col="count")
        marks = overlaps.df.loc[:, "count"] > 0
        vplist_filtered = VariantPlusList()
        vplist_filtered.extend(itertools.compress(self, marks))
        return vplist_filtered

    def get_sigresult(self, catalogue_type="sbs96", **kwargs):
        libsig = importlib.import_module(
            ".".join([top_package_name, "signature", "signatureresult"])
        )

        vcfspec_iter = (vp.vcfspec for vp in self)
        refver = self[0].refver
        sigresult = libsig.get_sigresult_from_vcfspecs(
            vcfspec_iter, refver=refver, catalogue_type=catalogue_type, **kwargs
        )
        return sigresult


def get_vp_sortkey(chromdict):
    vr_sortkey = common.get_vr_sortkey(chromdict)

    def vp_sortkey(vp):
        return vr_sortkey(vp.vr)

    return vp_sortkey

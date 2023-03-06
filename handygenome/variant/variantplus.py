import re
import pprint
import os
import itertools
import functools
import random
import shutil
import logging
import uuid

import pysam
import pyranges as pr
import numpy as np
import pandas as pd

import handygenome.common as common
import handygenome.workflow as workflow
import handygenome.variant.varianthandler as varianthandler
import handygenome.variant.filter as libfilter
import handygenome.variant.infoformat as infoformat
import handygenome.sv.breakends as libbnd
import handygenome.plot.filterinfo as plot_filterinfo
import handygenome.variant.infoformat as infoformat
import handygenome.vcfeditor.initvcf as initvcf
import handygenome.vcfeditor.headerhandler as headerhandler
import handygenome.annotation.data as annotdata
import handygenome.annotation.ensembl_feature as ensembl_feature
import handygenome.annotation.ensembl_parser as ensembl_parser
import handygenome.annotation.ensembl_rest as ensembl_rest
import handygenome.annotation.popfreq as libpopfreq
import handygenome.annotation.cosmic as libcosmic
import handygenome.annotation.oncokb as liboncokb
import handygenome.read.readplus as readplus
import handygenome.annotation.readstats as libreadstats
import handygenome.variant.vcfspec as libvcfspec
import handygenome.vcfeditor.indexing as indexing
import handygenome.cnv.misc as cnvmisc
from handygenome.variant.filter import FilterResultInfo, FilterResultFormat
import handygenome.variant.ponbams as libponbams


#READCOUNT_FORMAT_KEY = "allele_readcounts"


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
    # CONSTRUCTORS
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
        result.vcfspec = libvcfspec.Vcfspec.from_vr(result.vr, refver=result.refver)

        result.init_common(**kwargs)

        return result

    @classmethod
    def from_vcfspec(cls, vcfspec, chromdict=None, **kwargs):
        result = cls()

        result.vcfspec = vcfspec
        result.refver = result.vcfspec.refver
        result.fasta = result.vcfspec.fasta
        result.chromdict = (
            #common.ChromDict(fasta=result.fasta) 
            common.DEFAULT_CHROMDICTS[result.refver]
            if chromdict is None else 
            chromdict
        )
        result.vr = initvcf.create_vr(chromdict=result.chromdict, vcfspec=result.vcfspec)

        result.init_common(**kwargs)

        return result

    @classmethod
    def from_bnds(cls, bnds, chromdict=None, **kwargs):
        vcfspec = bnds.get_vcfspec_bnd1()
        return cls.from_vcfspec(vcfspec, chromdict, **kwargs)

    def init_common(
        self, 
        init_popfreq=True, 
        init_cosmic=True,
        popfreq_metadata=None, 
        cosmic_metadata=None,

        init_transcript=True,
        init_regulatory=True,
        init_motif=True,
        init_repeat=True,

        init_readstats=True,

        init_oncokb=True,

        init_filterresult=True,

        sampleid_list=None,
    ):
        # SV attrs
        self.is_sv = varianthandler.check_SV(self.vr)
        self.set_bnds_attributes()  # is_bnd1, bnds

        # popfreq & cosmic
        self._init_popfreq(metadata=popfreq_metadata, init_blank=(not init_popfreq))
        self._init_cosmic(metadata=cosmic_metadata, init_blank=(not init_cosmic))

        # features
        self._init_transcript(init_blank=(not init_transcript))
        self._init_regulatory(init_blank=(not init_regulatory))
        self._init_motif(init_blank=(not init_motif))
        self._init_repeat(init_blank=(not init_repeat))

        # readstats, readstats_data, rpplists
        self.rpplist_dict = dict()
        self.readstats_data_dict = dict()
        self._init_readstats(init_blank=(not init_readstats), sampleid_list=sampleid_list)

        # oncokb
        self._init_oncokb(init_blank=(not init_oncokb))

        # filter result
        self._init_filterresult(init_blank=(not init_filterresult))

    def __repr__(self):
        if self.transcript is None:
            gene_string = str(None)
        else:
            gene_string_buffer = {alt_idx: set() for alt_idx in range(len(self.vcfspec.alts))}
            for alt_idx, tr_set in enumerate(self.transcript):
                for tr in tr_set.canon_ovlp.values():
                    if tr.get_hgvsp_genename(one_letter=True) is None:
                        if tr.hgvsc_genename is None:
                            gene_string_buffer[alt_idx].add(str(None))
                        else:
                            gene_string_buffer[alt_idx].add(tr.hgvsc_genename)
                    else:
                        gene_string_buffer[alt_idx].add(tr.get_hgvsp_genename(one_letter=True))

                if (len(gene_string_buffer[alt_idx]) > 1) and (str(None) in gene_string_buffer[alt_idx]):
                    gene_string_buffer[alt_idx].remove(str(None))
                
            gene_string = '; '.join(f'alt_index {key}: {", ".join(val)}' for key, val in gene_string_buffer.items())

        #return f'<VariantPlus(chrom={self.vcfspec.chrom}, pos={self.vcfspec.pos}, ref={self.vcfspec.ref}, alts={self.vcfspec.alts}, gene={gene_string})>'
        return f'<VariantPlus(vcfspec={self.vcfspec}, gene=({gene_string}))>'

    # BASICS
    @property
    def chrom(self):
        return self.vcfspec.chrom

    @property
    def pos(self):
        return self.vcfspec.pos

    @property
    def ref(self):
        return self.vcfspec.ref

    @property
    def alts(self):
        return self.vcfspec.alts

    @property
    def alleles(self):
        return self.vcfspec.alleles

    @property
    def start0(self):
        return self.vcfspec.start0

    @property
    def end0(self):
        return self.vcfspec.end0

    @property
    def sampleids(self):
        return tuple(self.vr.header.samples)

    # BREAKENDS
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
            vr_svinfo = libbnd.get_vr_svinfo_standard_vr(
                self.vr,
                self.fasta,
                self.chromdict,
            )
            self.is_bnd1 = vr_svinfo["is_bnd1"]
            self.bnds = libbnd.get_bnds_from_vr_svinfo(
                self.vr, vr_svinfo, self.fasta, self.chromdict
            )
        else:
            self.is_bnd1 = None
            self.bnds = None

    # ANNOTATION GENERAL
    def write_annots_to_vr(self, fill_missing_sample=True):
        # Vcfspec
        if not self.vcfspec.components.is_missing:
            self.vcfspec.components.write(vr)

        # INFO
        for key in (
            'popfreq',
            'cosmic',
            'transcript',
            'regulatory',
            'motif',
            'repeat',
            'oncokb',
            #'filter',
        ):
            annot = getattr(self, key)
            if not annot.is_missing:
                annot.write(self.vr)
        # FORMAT
        for key in (
            'readstats_dict',
            #'filter_bysample',
        ):
            annot = getattr(self, key)
            if not annot.is_missing:
                if fill_missing_sample:
                    self.vr = annot.get_sample_augmented_vr(self.vr)
                self.vr = annot.write(self.vr)
    
    # ANNOTATION - filterresult
    def _init_filterresult(self, init_blank=False):
        if init_blank:
            self.filter = FilterResultInfo.init_missing()
            self.filter_bysample = FilterResultFormat.init_missing(self.vr)
        else:
            self.filter = FilterResultInfo.from_vr(self.vr)
            self.filter_bysample = FilterResultFormat.from_vr(self.vr)

    # ANNOTATION - popfreq & cosmic
    def _init_popfreq(self, metadata=None, init_blank=False):
        if init_blank:
            self.popfreq = libpopfreq.PopfreqInfoALTlist.init_missing(self.vr)
        else:
            self.popfreq = libpopfreq.PopfreqInfoALTlist.from_vr(self.vr, metadata=metadata)

    def _init_cosmic(self, metadata=None, init_blank=False):
        if init_blank:
            self.cosmic = libcosmic.CosmicInfoALTlist.init_missing(self.vr)
        else:
            self.cosmic = libcosmic.CosmicInfoALTlist.from_vr(self.vr, metadata=metadata)

    def create_popfreq(self):
        dbsnp_vcf = annotdata.VCFS_DBSNP[self.refver]
        self.popfreq = libpopfreq.PopfreqInfoALTlist.from_vcfspec(
            vcfspec=self.vcfspec, 
            dbsnp_vcf=dbsnp_vcf,
            metadata=libpopfreq.PopfreqMetadata.from_vcfheader(dbsnp_vcf.header),
        )

    def create_cosmic(self):
        cosmic_vcf = annotdata.VCFS_COSMIC[self.refver]
        self.cosmic = libcosmic.CosmicInfoALTlist.from_vcfspec(
            vcfspec=self.vcfspec, 
            cosmic_vcf=cosmic_vcf,
            metadata=libcosmic.CosmicMetadata.from_vcfheader(cosmic_vcf.header),
        )

    # ANNOTATION - oncokb
    def _init_oncokb(self, init_blank=False):
        if init_blank:
            self.oncokb = liboncokb.OncoKBInfoALTlist.init_missing(self.vr)
        else:
            self.oncokb = liboncokb.OncoKBInfoALTlist.from_vr(self.vr)

    def create_oncokb(self, token):
        oncokb_altlist = liboncokb.OncoKBInfoALTlist()
        oncokb_altlist.extend(
            liboncokb.query_hgvsg_post(
                hgvsg_list=list(x.to_hgvsg() for x in self.vcfspec.iter_monoalts()), 
                token=token, tumor_type=None, evidence_types=None,
            )
        )
        self.oncokb = oncokb_altlist

    # ANNOTATION - ensembl features
    def _init_transcript(self, init_blank=False):
        if init_blank:
            self.transcript = ensembl_feature.TranscriptSetALTlist.init_missing(self.vr)
        else:
            self.transcript = ensembl_feature.TranscriptSetALTlist.from_vr(self.vr)

    def _init_regulatory(self, init_blank=False):
        if init_blank:
            self.regulatory = ensembl_feature.RegulatorySet.init_missing()
        else:
            self.regulatory = ensembl_feature.RegulatorySet.from_vr(self.vr)

    def _init_motif(self, init_blank=False):
        if init_blank:
            self.motif = ensembl_feature.MotifSet.init_missing()
        else:
            self.motif = ensembl_feature.MotifSet.from_vr(self.vr)

    def _init_repeat(self, init_blank=False):
        if init_blank:
            self.repeat = ensembl_feature.RepeatSet.init_missing()
        else:
            self.repeat = ensembl_feature.RepeatSet.from_vr(self.vr)

    def _init_features(self, init_blank=False):
        self._init_transcript(init_blank=init_blank)
        self._init_regulatory(init_blank=init_blank)
        self._init_motif(init_blank=init_blank)
        self._init_repeat(init_blank=init_blank)

    def create_features(
        self, vep_method='rest', distance=5000, 
        with_CADD=True, with_Phenotypes=False, with_canonical=True,
        with_mane=True, with_miRNA=False, with_numbers=True, 
        with_protein=True, with_ccds=True, with_hgvs=True,
    ):
        assert vep_method in ('rest', 'local')

        # run vep and parse the result
        rest_vep_results = ensembl_rest.vep_post(
            refver=self.refver,
            vcfspec_list=list(self.vcfspec.iter_monoalts()),
            distance=distance, 
            with_CADD=with_CADD, with_Phenotypes=with_Phenotypes, with_canonical=with_canonical,
            with_mane=with_mane, with_miRNA=with_miRNA, with_numbers=with_numbers, 
            with_protein=with_protein, with_ccds=with_ccds, with_hgvs=with_hgvs,
        )
        parsed = [ensembl_parser.parse_rest_vep(x) for x in rest_vep_results]
        transcript = ensembl_feature.TranscriptSetALTlist()
        for x in parsed:
            transcript.append(x['transcript'])
        regulatory = parsed[0]['regulatory']
        motif = parsed[0]['motif']
        # postprocess
        tabixfile_geneset = annotdata.TABIXFILES_GENESET[self.refver]
        tabixfile_regulatory = annotdata.TABIXFILES_REGULATORY[self.refver]
        tabixfile_repeats = annotdata.TABIXFILES_REPEATS[self.refver]
        # transcript
        for alt_idx, tr_set in enumerate(transcript):
            tr_set.update_ensembl_gff(
                vcfspec=self.vcfspec,
                chromdict=self.chromdict,
                distance=distance,
                tabixfile_geneset=tabixfile_geneset,
                refver=self.refver,
                fill_missing_canonical=True,
            )
            tr_set.set_missing_distances(vcfspec=self.vcfspec, alt_idx=alt_idx)
        # regulatory
        regulatory.update_ensembl_gff(
            vcfspec=self.vcfspec,
            chromdict=self.chromdict,
            distance=distance,
            tabixfile_regulatory=tabixfile_regulatory,
        )
        regulatory.set_missing_distances(vcfspec=self.vcfspec, alt_idx=0)
        # repeat
        repeat = ensembl_feature.RepeatSet.init_nonmissing()
        repeat.update_repeatmasker_bed(
            vcfspec=self.vcfspec,
            chromdict=self.chromdict,
            distance=distance,
            tabixfile_repeats=tabixfile_repeats,
        )
        repeat.set_missing_distances(vcfspec=self.vcfspec, alt_idx=0)

        # assign attributes
        self.transcript = transcript
        self.regulatory = regulatory
        self.motif = motif
        self.repeat = repeat

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
#    def get_other_alleleindexes(self, allele_index):
#        """Returns all integer allele indexes other than the input"""
#
#        all_alleleindexes = set(range(len(self.vr.alleles)))
#        other_alleleindexes = all_alleleindexes.difference([allele_index])
#        return tuple(sorted(other_alleleindexes))

    # ANNOTATION - readstats and related
    def _init_readstats(self, init_blank=False, sampleid_list=None):
        if init_blank:
            self.readstats_dict = libreadstats.ReadStatsSampledict.init_missing(self.vr)
        else:
            self.readstats_dict = libreadstats.ReadStatsSampledict.from_vr(self.vr, sampleid_list=sampleid_list)

    def create_readstats(
        self, bam_dict,
        rpplist_kwargs={
            'view': False,
            'no_matesearch': True,
        },
        alleleinfo_kwargs={},
    ):
        """Args:
            bam_dict: keys - sample ID, values - pysam.AlignmentFile
        Makes a ReadStatsSampledict object, keys of which are the same as the
            keys of "bam_dict" parameter.
        """

        self.readstats_dict = libreadstats.ReadStatsSampledict()
        for sampleid, bam in bam_dict.items():
            self.readstats_dict[sampleid] = libreadstats.ReadStats.from_bam(
                vcfspec=self.vcfspec, 
                bam=bam, 
                fasta=self.fasta, 
                chromdict=self.chromdict,
                rpplist_kwargs=rpplist_kwargs,
                alleleinfo_kwargs=alleleinfo_kwargs,
            )

    def update_readstats(
        self, bam_dict, 
        rpplist_kwargs={
            'view': False,
            'no_matesearch': True,
        },
        alleleinfo_kwargs={},
        countonly=False,
    ):
        self.readstats_dict.update_bams(
            bam_dict=bam_dict,
            rpplist_kwargs=rpplist_kwargs,
            alleleinfo_kwargs=alleleinfo_kwargs,
            countonly=countonly,
        )

    def make_rpplist(
        self,
        bam,
        set_alleleclass=True,
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
            if set_alleleclass:
                rpplist.update_alleleinfo_sv(bnds=self.bnds, **alleleinfo_kwargs)
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
            if set_alleleclass:
                rpplist.update_alleleclass(vcfspec=self.vcfspec, **alleleinfo_kwargs)

        return rpplist

    def set_rpplist(
        self, sampleid, bam, 
        set_alleleclass=True,
        rpplist_kwargs={
            'view': False,
            'no_matesearch': True,
        },
        alleleinfo_kwargs=dict(),
    ):
        self.rpplist_dict[sampleid] = self.make_rpplist(
            bam, set_alleleclass=set_alleleclass,
            rpplist_kwargs=rpplist_kwargs,
            alleleinfo_kwargs=alleleinfo_kwargs,
        )

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

        #if sampleid not in self.rpplist_dict.keys():
        self.set_rpplist(
            sampleid, 
            bam, 
            set_alleleclass=True,
            rpplist_kwargs=rpplist_kwargs,
            alleleinfo_kwargs=alleleinfo_kwargs,
        )

#        readstats_data = libreadstats.rpplist_to_readstats_data(
#            self.rpplist_dict[sampleid],
#            self.vcfspec,
#        )
        self.readstats_data_dict[sampleid] = libreadstats.rpplist_to_readstats_data(
            self.rpplist_dict[sampleid],
            self.vcfspec,
        )

#    def calc_readcounts(self, bam, no_matesearch=True):
#        """Designed only for non-sv"""
#
#        rpplist = self.make_rpplist(
#            bam, no_matesearch=no_matesearch, set_alleleclass=True,
#        )
#        readcounts = rpplist.get_readcounts(self.vcfspec)
#
#        return readcounts

    def get_total_rppcount(self, sampleid, exclude_other=False):
        """Sum of rpp counts for all alleleclasses except None"""
        readstats = self.readstats_dict[sampleid]
        return readstats.get_total_rppcount(exclude_other=exclude_other)

    def get_vaf(self, sampleid, allele_index=1, exclude_other=False, ndigits=None):
        readstats = self.readstats_dict[sampleid]
        vaf = readstats.get_vaf(alleleclass=allele_index, exclude_other=exclude_other)
        if np.isnan(vaf):
            return vaf
        else:
            if ndigits is None:
                return vaf
            else:
                return round(vaf, ndigits)

    def get_sorted_vafs(self, sampleid, exclude_other=True, reverse=True):
        return self.readstats_dict[sampleid].get_sorted_vafs(exclude_other=exclude_other, reverse=reverse)

    def show_vafs(self, sampleid, exclude_other=False):
        for allele_index, allele in enumerate(self.alleles):
            print(allele, self.get_vaf(sampleid, allele_index=allele_index, exclude_other=exclude_other))

    # CNV-RELATED
    def fetch_cnvinfo(self, segments_df):
        nearest = self.get_gr().nearest(pr.PyRanges(segments_df))
        seg_row = nearest.df.iloc[0, :]
        return int(seg_row.CNt), int(seg_row.A), int(seg_row.B)

    def get_ccf_CNm(self, segments_df, is_female, cellularity, sampleid, allele_index=1, likelihood='diff'):
        # set params
        CNt, A, B = self.fetch_cnvinfo(segments_df)
        if cnvmisc.check_haploid(is_female, self.chrom):
            CNn = 1
            CNB = None
        else:
            CNn = 2
            CNB = B

        if likelihood == 'binom':
            readstats = self.readstats_dict[sampleid]
            total_read = readstats.get_total_rppcount()
            var_read = readstats['rppcounts'][allele_index]
        elif likelihood == 'diff':
            total_read = None
            var_read = None

        return cnvmisc.get_ccf_CNm(
            vaf=self.get_vaf(sampleid, allele_index), 
            cellularity=cellularity, 
            CNt=CNt, 
            CNB=CNB, 
            likelihood=likelihood, 
            total_read=total_read, 
            var_read=var_read, 
            CNn=CNn,
        )

    # filters
    def reset_sample_filter(self, sampleid):
        self.set_format(sampleid, libfilter.FORMAT_FILTER_META["ID"], None)

    def reset_sample_filter_all(self):
        for sampleid in self.vr.header.samples:
            self.reset_sample_filter(sampleid)

    def check_sample_filter(self, sampleid):
        sample_filter = self.get_format(
            sampleid, libfilter.FORMAT_FILTER_META["ID"], collapse_tuple=False
        )
        return sample_filter == ("PASS",)

    def add_sample_filter(self, sampleid, value: str):
        if value == "PASS":
            self.set_format(sampleid, libfilter.FORMAT_FILTER_META["ID"], ("PASS",))
        else:
            if self.check_NA_format(sampleid, libfilter.FORMAT_FILTER_META["ID"]):
                new_val = (value,)
            else:
                old_val = self.get_format(
                    sampleid, libfilter.FORMAT_FILTER_META["ID"], collapse_tuple=False
                )
                new_val = tuple(set(old_val + (value,)))

            self.set_format(sampleid, libfilter.FORMAT_FILTER_META["ID"], new_val)

    def show_sample_filters(self):
        for sampleid in self.vr.header.samples:
            print(
                sampleid,
                self.get_format(
                    sampleid, libfilter.FORMAT_FILTER_META["ID"], collapse_tuple=False
                ),
            )

    # visualizations
    def show_igv(
        self, igv, bam_dict, 
        max_width_oneside=200,
        new=True,
        readheight='collapse',
        viewaspairs=True,
        rpplist_kwargs={
            'view': True,
            'no_matesearch': True,
        },
        alleleinfo_kwargs=dict(),
    ):
        if self.is_sv:
            raise Exception(f'SV is not supported yet')

        assert readheight in ('collapse', 'expand', 'squish')

        # make rpplists
        rpplist_dict = dict()
        for sampleid, bam in bam_dict.items():
            rpplist = self.make_rpplist(
                bam, 
                set_alleleclass=True,
                rpplist_kwargs=rpplist_kwargs,
                alleleinfo_kwargs=alleleinfo_kwargs,
            )
            if self.is_sv:
                rpplist.set_alleleinfo_tag_rp_sv(self.bnds)
                rpplist.set_alleleinfo_tag_rpp_sv(self.bnds)
            else:
                rpplist.set_alleleclass_tag_rp(self.vcfspec)
                rpplist.set_alleleclass_tag_rpp(self.vcfspec)

            rpplist_dict[sampleid] = rpplist

        # make temporary directory and write temp bam files
        tmpbam_dir = workflow.get_tmpfile_path(prefix='tmpdir_show_igv_', dir=os.getcwd(), delete=False, is_dir=True)
        tmpbam_paths = {
            sampleid: os.path.join(tmpbam_dir, f'{sampleid}_for_show.bam')
            for sampleid in bam_dict.keys()
        }
        for sampleid, rpplist in rpplist_dict.items():
            rpplist.write_bam(outfile_path=tmpbam_paths[sampleid])
        # igv manipulation
        if new:
            igv.cmd('new')
        igv.load(tuple(tmpbam_paths.values()))

        if not self.is_sv:
            mins = list()
            maxs = list()
            for rpplist in rpplist_dict.values():
                reads_ref_range0 = rpplist.get_ref_range0(self.vcfspec.chrom)
                mins.append(reads_ref_range0.start)
                maxs.append(reads_ref_range0.stop)

            show_range_start = max(self.vcfspec.start0 - max_width_oneside, min(mins))
            show_range_end = min(self.vcfspec.end0 + max_width_oneside, max(maxs))
            show_locus = (self.vcfspec.chrom, show_range_start, show_range_end)

            igv.goto([show_locus], width=0)
        else:
            pass

        if viewaspairs:
            igv.viewaspairs()
        else:
            igv.viewaspairs_off()

        igv.cmd(readheight)
        igv.cmd("group TAG a2")
        igv.cmd(f"colorBy TAG {readplus.ALLELECLASS_TAG_RPP}")
        shutil.rmtree(tmpbam_dir)

    def show_readstats_data(self, bam, alt_index=0, varpos_key='left'):
        readstats_data = self.make_readstats_data(bam)
        fig = plot_filterinfo.show_readstats_data(readstats_data, alt_index=alt_index, title=repr(self.vcfspec), varpos_key=varpos_key)
        fig.show()

    show_filterinfo = show_readstats_data  # alias

    def show_readcounts(self, allele_index=1, samples=None, mask=None):
        fig, ax_bar, ax_dot = plot_filterinfo.show_readcounts(
            self.readstats_dict, 
            allele_index=allele_index,
            samples=samples,
            mask=mask,
            title=repr(self.vcfspec),
        )
        fig.show()

    #@common.get_deco_num_set_differently(('pon_samples', 'pon_cohorts'), 1)
    def show_pon(self, query_sample, allele_index=1, exclude_query=True, pon_cohorts=None, pon_samples=None, **kwargs):
        fig = plot_filterinfo.show_pon(
            query_sample=query_sample, 
            vp=self, 
            allele_index=allele_index, 
            exclude_query=exclude_query, 
            pon_cohorts=pon_cohorts,
            pon_samples=pon_samples, 
            **kwargs,
        )
        fig.show()
        

#    #@common.get_deco_num_set(('bam', 'sampleid'), 1)
#    def show_igv_old(
#        self, igv, bam=None, sampleid=None, tmpbam_name=None,
#        rpplist_kwargs={
#            'view': True,
#            'no_matesearch': True,
#        },
#        alleleclass_kwargs=dict(),
#        viewaspairs=True,
#    ):
#        """Args:
#            tmpbam_name: Basename of temporary bam. It will be shown as track name on IGV. If not set, when "sampleid" argument is set, it will be "<sampleid>_for_show.bam".
#        """
#        # sanity check
#        if tmpbam_name is not None:
#            if not tmpbam_name.endswith('.bam'):
#                raise Exception(f'"tmpbam_name" must end with ".bam". Otherwise, IGV cannot appropriately load it.')
#
#        # get rpplist
#        if sampleid is not None:
#            if sampleid not in self.rpplist_dict.keys():
#                raise Exception(
#                    f"rpplist for the sampleid {sampleid} is not created. To bypass this problem, 1) run this method with 'bam' argument, or 2) run VariantPlus.set_rpplist method first."
#                )
#            rpplist = self.rpplist_dict[sampleid]
#        elif bam is not None:
#            rpplist = self.make_rpplist(
#                bam, set_alleleclass=True,
#                rpplist_kwargs=rpplist_kwargs,
#                alleleclass_kwargs=alleleclass_kwargs,
#            )
#
#        if len(rpplist) == 0:
#            return
#
#        # add alleleclass tag to ReadPlus and ReadPlusPair objects
#        if self.is_sv:
#            rpplist.set_alleleinfo_tag_rp_sv(self.bnds)
#            rpplist.set_alleleinfo_tag_rpp_sv(self.bnds)
#        else:
#            rpplist.set_alleleclass_tag_rp(self.vcfspec)
#            rpplist.set_alleleclass_tag_rpp(self.vcfspec)
#
#        # set tmpbam_path
#        tmpbam_dir = workflow.get_tmpfile_path(prefix='tmpdir_show_igv_', dir=os.getcwd(), delete=False, is_dir=True)
#        if tmpbam_name is None:  # tmpbam_name is basename of temporary bam path
#            if sampleid is None:
#                tmpbam_path = workflow.get_tmpfile_path(suffix='.bam', dir=tmpbam_dir, delete=False)
#            else:
#                tmpbam_path = os.path.join(tmpbam_dir, f'{sampleid}_for_show.bam')
#        else:
#            tmpbam_path = os.path.join(tmpbam_dir, tmpbam_name)
#
#        # main
#        rpplist.write_bam(outfile_path=tmpbam_path)
#        igv.load([tmpbam_path])
#
#        if self.is_sv:
#            (reads_range0_bnd1, reads_range0_bnd2) = rpplist.get_ref_range0_sv(
#                self.bnds
#            )
#            pos_range0_bnd1 = self.bnds.get_pos_range0_bnd1()
#            pos_range0_bnd2 = self.bnds.get_pos_range0_bnd2()
#
#            width_bnd1 = max(
#                pos_range0_bnd1.start - reads_range0_bnd1.start,
#                reads_range0_bnd1.stop - pos_range0_bnd1.stop,
#            )
#            width_bnd2 = max(
#                pos_range0_bnd2.start - reads_range0_bnd2.start,
#                reads_range0_bnd2.stop - pos_range0_bnd2.stop,
#            )
#            width = max(width_bnd1, width_bnd2)
#
#            igv.goto([self.bnds], width=width)
#        else:
#            reads_ref_range0 = rpplist.get_ref_range0(self.vcfspec.chrom)
#            width = max(
#                self.vcfspec.pos0 - reads_ref_range0.start,
#                reads_ref_range0.stop - self.vcfspec.pos0,
#            )
#
#            igv.goto([self.vcfspec], width=width)
#
#        igv.cmd("group")
#
#        if viewaspairs:
#            igv.viewaspairs()
#        else:
#            igv.viewaspairs_off()
#
#        igv.cmd(f"colorBy TAG {readplus.ALLELEINFO_TAG_RPP}")
#
#        shutil.rmtree(tmpbam_dir)

    #def show_readcounts(self, sampleids=None, **kwargs):
    #    pass
        #if sampleids is None:
        #    sampleids = sorted(self.vr.samples.keys())

        #variantviz.show_readcounts(
        #    self, sampleid_order=sampleids, title=str(self.vcfspec), **kwargs
        #)

    # miscellaneous
    def get_gr(self):
        return pr.from_dict(
            {
                "Chromosome": [self.vr.contig],
                "Start": [self.vr.pos - 1],
                "End": [self.vr.pos],
            }
        )

    def get_gene_names(self, canonical=True, overlap=True, coding=False):
        result = set()
        for tr_set in self.transcript:
            for tr in tr_set.values():
                if canonical:
                    if not tr['is_canonical']:
                        continue
                if overlap:
                    if tr['distance'] is not None:
                        continue
                if coding:
                    if not tr['subtype_flags']['coding']:
                        continue
                result.add(tr['gene_name'])

        return result
#        return set(
#            itertools.chain.from_iterable(
#                tr_set.get_gene_names(canonical=canonical, overlap=overlap)
#                for tr_set in self.transcript
#            )
#        )

    def get_protein_changes(self, one_letter=True):
        result = set()
        for tr_set in self.transcript:
            for tr in tr_set.canon_ovlp.values():
                hgvsp = tr.get_hgvsp_genename(one_letter=one_letter)
                if hgvsp is not None:
                    result.add(hgvsp)
                
        return result

    #def check_intergenic(self):
    #    return len(self.annotdb.transcript) == 0

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


class VariantPlusList(list):
    ################
    # initializers #
    ################

    def __init__(
        self, vcf_path=None, vcf=None, refver=None,
        fasta=None, chromdict=None,
        logging_lineno=1000, preview_lineno=30, 
        verbose=True, logger=None,
    ):
        # vcf_path, vcf
        self.vcf_path = vcf_path
        if vcf is None:
            if self.vcf_path is None:
                self.vcf = None
            else:
                self.vcf = pysam.VariantFile(self.vcf_path, "r")
        else:
            self.vcf = vcf

        # refver, fasta, chromdict
        if refver is None:
            if self.vcf is None:
                raise Exception(f'When "refver" is not set, "vcf_path" must be set.')
            self.refver = common.infer_refver_vcfheader(self.vcf.header)
        else:
            self.refver = refver

        if fasta is None:
            self.fasta = common.DEFAULT_FASTAS[self.refver]
        else:
            self.fasta = fasta

        if chromdict is None:
            self.chromdict = common.ChromDict(fasta=self.fasta)
        else:
            self.chromdict = chromdict

        # lineno, verbose, logger
        self.logging_lineno = logging_lineno
        self.preview_lineno = preview_lineno
        self.verbose = verbose
        if logger is None:
            self.logger = self.get_logger(self.verbose)
        else:
            self.logger = logger

    @classmethod
    def concat(cls, others, **kwargs):
        others = list(others)
        if len(set(x.refver for x in others)) != 1:
            raise Exception(f'Different refvers of input VariantPlusList objects')

        result = cls(refver=others[0].refver, **kwargs)
        for other in others:
            result.extend(other)

        return result

    @classmethod
    def from_vps(cls, vps, **kwargs):
        vp_iterator = iter(vps)
        vp = next(vp_iterator)
        result = cls(refver=vp.refver, **kwargs)
        result.append(vp)
        result.extend(vp_iterator)
        return result

    @classmethod
    def from_vcf(
        cls,
        vcf_path,
        init_all_attrs=True,
        init_popfreq=False,
        init_cosmic=False,
        init_transcript=False,
        init_regulatory=False,
        init_motif=False,
        init_repeat=False,
        init_readstats=False,
        init_oncokb=False,
        init_filterresult=False,
        sampleid_list=None,
        prop=None,
        **kwargs,
    ):
        """Args:
            prop: probability for random selection. If None, all entries are selected.
        """
        # initialize results
        result = cls(vcf_path=vcf_path, **kwargs)

        # set params
        if init_all_attrs:
            init_popfreq = True
            init_cosmic = True
            init_transcript = True
            init_regulatory = True
            init_motif = True
            init_repeat = True
            init_readstats = True
            init_oncokb = True
            init_filterresult = True

        if init_popfreq:
            popfreq_metadata = libpopfreq.PopfreqMetadata.from_vcfheader(result.vcf.header)
        else:
            popfreq_metadata = None

        if init_cosmic:
            cosmic_metadata = libcosmic.CosmicMetadata.from_vcfheader(result.vcf.header)
        else:
            cosmic_metadata = None

#        '''
#        Following code is a failed attempt to parallelize using multiprocessing.
#        Many custom classes cannot be pickled.
#
#        def vr_iterator(vcf):
#            for vr in vcf.fetch():
#                if varianthandler.check_SV(vr):
#                    vr_svinfo = libbnd.get_vr_svinfo_standard_vr(vr, result.fasta, result.chromdict)
#                    if vr_svinfo["is_bnd1"]:
#                        yield vr
#                else:
#                    yield vr
#
#        refver, fasta, chromdict = result.refver, result.fasta, result.chromdict
#        with multiprocessing.Pool(ncore) as p:
#            NR = 0
#            for vr_subiter in common.grouper(vr_iterator(result.vcf), result.logging_lineno):
#                vp_sublist = p.starmap(
#                    _init_helper_make_vp, 
#                    (
#                        (
#                            vr, 
#                            refver, 
#                            fasta, 
#                            chromdict, 
#                            init_popfreq, 
#                            init_cosmic,
#                            popfreq_metadata, 
#                            cosmic_metadata,
#                            init_transcript,
#                            init_regulatory,
#                            init_motif,
#                            init_repeat,
#                            init_readstats,
#                            init_oncokb,
#                        )
#                        for vr in vr_subiter
#                    )
#                )
#                NR += len(vp_sublist)
#                result.logger.info(f'Processing {NR:,}th line')
#                result.extend(vp_sublist)
#
#        '''

        # create vps
        if prop is None:
            vr_iterator = common.iter_lineno_logging(result.vcf.fetch(), result.logger, result.logging_lineno)
        else:
            vr_iterator = common.iter_lineno_logging(
                common.bernoulli_iterator(result.vcf.fetch(), p=prop, block_size=int(1e5)), 
                result.logger, 
                result.logging_lineno,
            )

        for vr in vr_iterator:
            if varianthandler.check_SV(vr):
                vr_svinfo = libbnd.get_vr_svinfo_standard_vr(vr, result.fasta, result.chromdict)
                if not vr_svinfo["is_bnd1"]:
                    continue

            vp = VariantPlus.from_vr(
                vr=vr,
                refver=result.refver,
                fasta=result.fasta,
                chromdict=result.chromdict,

                init_popfreq=init_popfreq,
                init_cosmic=init_cosmic,
                popfreq_metadata=popfreq_metadata, 
                cosmic_metadata=cosmic_metadata,

                init_transcript=init_transcript,
                init_regulatory=init_regulatory,
                init_motif=init_motif,
                init_repeat=init_repeat,
                init_readstats=init_readstats,
                init_oncokb=init_oncokb,
                init_filterresult=init_filterresult,

                sampleid_list=sampleid_list,
            )
            result.append(vp)

        return result

    ########
    # repr #
    ########

    def __repr__(self):
        return self.show(self.preview_lineno)

    def show(self, lineno):
        tmp = list()
        tmp.append(f'<VariantPlusList of length {len(self)} [')
        for idx, vp in enumerate(self):
            if idx == lineno:
                break
            tmp.append(f'\t{idx}\t{vp}')
        if len(self) > lineno:
            tmp.append(f'\t...')
        tmp.append(f']')
        return '\n'.join(tmp)

    @staticmethod
    def get_logger(verbose):
        formatter = logging.Formatter(
            fmt='[%(asctime)s.%(msecs)03d] VariantPlusList: %(message)s', 
            datefmt='%Z %Y-%m-%d %H:%M:%S'
        )
        return workflow.get_logger(
            name=str(uuid.uuid4()),
            level=('info' if verbose else 'error'),
            #level='info',
            formatter=formatter,
        )

    ##########
    # others #
    ##########

    @functools.cached_property
    def vcfspec_list(self):
        return [vp.vcfspec for vp in self]

    @functools.cached_property
    def index_gr(self):
        gr = self.get_df(vaf_sampleid=None, as_gr=True, omit_vaf=True)
        gr.Vcfspec = self.vcfspec_list
        gr.Index = list(range(gr.df.shape[0]))  # In this form of assignment, left-side term must be a list.
        return gr

    def fetch_by_vcfspec(self, vcfspec):
        selector = (self.index_gr.Vcfspec == vcfspec)
        if not selector.any():
            raise ValueError(f'Input vcfspec is not present in the VariantPlusList.')
        #idx = self.index_gr.Vcfspec.loc[selector].index[0]
        idx = self.index_gr.Index.loc[selector][0]
        return self[idx]

    def fetch_by_coord(self, chrom, start0=None, end0=None):
        if (start0 is None) and (end0 is None):
            gr_subset = self.index_gr[chrom]
        else:
            if start0 is None:
                start0 = 0
            if end0 is None:
                end0 = self.chromdict[chrom]
            gr_subset = self.index_gr[chrom, start0:end0]

        if gr_subset.df.shape[0] == 0:
            return list()
        else:
            sl = slice(gr_subset.Index.iloc[0], gr_subset.Index.iloc[-1] + 1)
            return self[sl]

    def sample(self, n=1):
        if len(self) == 0:
            raise Exception(f'The VariantPlusList is empty.')

        if n == 1:
            return random.choice(self)
        else:
            return random.sample(self, k=n)

    def get_df(
        self, 
        alt_index=0,
        vaf_sampleid=None, 
        as_gr=False, 
        omit_vaf=False, 
        get_vaf_kwargs={
            #'allele_index': 1, 
            'exclude_other': False, 
            'ndigits': None,
        },
    ):
        # parameter setups
        get_vaf_kwargs['allele_index'] = alt_index + 1

        if vaf_sampleid is None:
            sample_ids = list(self[0].vr.header.samples)
            if len(sample_ids) == 0:
                vaf_sampleid = None
                omit_vaf = True
            else:
                vaf_sampleid = sample_ids[0]

        # main
        chroms = list()
        pos1s = list()
        starts = list()
        ends = list()
        refs = list()
        alts = list()
        vafs = list()
        for vp in self:
            chroms.append(vp.vcfspec.chrom)
            pos1s.append(vp.vcfspec.pos)
            starts.append(vp.vcfspec.pos0)
            ends.append(vp.vcfspec.end0)
            refs.append(vp.vcfspec.ref)
            alts.append(vp.vcfspec.alts[alt_index])

            if not omit_vaf:
                if vaf_sampleid is None:
                    vafs.append(np.nan)
                else:
                    vaf = vp.get_vaf(vaf_sampleid, **get_vaf_kwargs)
                    vafs.append(vaf)  # vaf can be np.nan
        # make data
        data = {
            'Chromosome': chroms,
            'Start': starts,
            'End': ends,
            'POS': pos1s,
            'REF': refs,
            'ALT': alts,
        }
        if not omit_vaf:
            data.update({'vaf': vafs})
        # result
        if as_gr:
            return pr.from_dict(data)
        else:
            return pd.DataFrame.from_dict(data)

    def get_gr(self, vaf_sampleid=None):
        return self.get_df(vaf_sampleid=vaf_sampleid, as_gr=True)

    def spawn(self):
        kwargs = {
            key: getattr(self, key) for key in 
            (
                'vcf_path', 'vcf', 'refver', 
                'logging_lineno', 'preview_lineno', 'verbose', 
                'fasta', 'chromdict',
            )
        }
        result = self.__class__(**kwargs)
        return result

    def filter(self, key):
        result = self.spawn()
        for vp in self:
            if key(vp):
                result.append(vp)
        return result
    
    def isec(self, gr):
        overlaps = self.get_gr().count_overlaps(gr, overlap_col="count")
        marks = overlaps.df.loc[:, "count"] > 0
        result = self.spawn()
        result.extend(itertools.compress(self, marks))
        return result

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

    def select_by_mutation_type(self, muttype):
        result = VariantPlusList(refver=self.refver)
        for vp in self:
            if vp.vcfspec.get_mutation_type() == muttype:
                result.append(vp)
        return result

    def select_ins(self):
        return self.select_by_mutation_type('ins')

    def select_del(self):
        return self.select_by_mutation_type('del')

    def select_indel(self):
        result = VariantPlusList(refver=self.refver)
        result.extend(self.select_ins())
        result.extend(self.select_del())
        return result

    def get_vp_sortkey(self):
        vr_sortkey = common.get_vr_sortkey(self.chromdict)
        def vp_sortkey(vp):
            return vr_sortkey(vp.vr)
        return vp_sortkey

    def sort(self):
        super().sort(key=self.get_vp_sortkey())

    # writing related methods
    def get_output_header(self):
        if len(self) == 0:
            merged_header = initvcf.create_header(self.chromdict)
            conflicting_keys = None
            return merged_header, conflicting_keys
        else:
            return headerhandler.merge_vcfheaders(vp.vr.header for vp in self)

    @staticmethod
    def write_each_vr(vr, out_vcf, conflicting_keys):
        # reheader if sample ordering is different
        if tuple(vr.header.samples) != tuple(out_vcf.header.samples):
            out_vr = varianthandler.reheader(vr, out_vcf.header)
        else:
            out_vr = vr.copy()

        # remove conflicting keys
        for key in conflicting_keys['info']:
            if key in out_vr.info:
                del out_vr.info[key]
        for key in conflicting_keys['format']:
            for format_item in out_vr.samples.values():
                if key in format_item:
                    del format_item[key]

        # write
        out_vcf.write(out_vr)

    def write(self, outfile_path, mode_bcftools="z", mode_pysam=None, index=True):
        # write annotation items to vr
        for vp in self:
            vp.write_annots_to_vr(fill_missing_sample=True)
        # prepare header
        header, conflicting_keys = self.get_output_header()
        # main
        self.sort()
        mode_pysam = common.write_mode_arghandler(mode_bcftools, mode_pysam)
        with pysam.VariantFile(outfile_path, mode=mode_pysam, header=header) as out_vcf:
            #for vp in self:
            for vp in common.iter_lineno_logging(self, self.logger, self.logging_lineno):
                if vp.is_sv:
                    self.write_each_vr(vp.vr, out_vcf, conflicting_keys)
                    self.write_each_vr(vp.get_vr_bnd2(), out_vcf, conflicting_keys)
                else:
                    self.write_each_vr(vp.vr, out_vcf, conflicting_keys)
        if index:
            indexing.index_vcf(outfile_path)


# not used
def _init_helper_make_vp(
    vr, refver, fasta, chromdict, 

    init_popfreq, 
    init_cosmic,
    popfreq_metadata, 
    cosmic_metadata,

    init_transcript,
    init_regulatory,
    init_motif,
    init_repeat,
    init_readstats,
    init_oncokb,
):
    return VariantPlus.from_vr(
        vr=vr,
        refver=result.refver,
        fasta=result.fasta,
        chromdict=result.chromdict,

        init_popfreq=init_popfreq,
        init_cosmic=init_cosmic,
        popfreq_metadata=popfreq_metadata, 
        cosmic_metadata=cosmic_metadata,

        init_transcript=init_transcript,
        init_regulatory=init_regulatory,
        init_motif=init_motif,
        init_repeat=init_repeat,
        init_readstats=init_readstats,
        init_oncokb=init_oncokb,
    )


def load_vcf(vcf_path):
    return VariantPlusList.from_vcf(vcf_path)



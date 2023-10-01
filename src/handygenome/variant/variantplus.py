import re
import pprint
import os
import itertools
import functools
import collections
import random
import shutil
import logging
import uuid
import array
import multiprocessing
#import importlib

import pysam
import pyranges as pr
import numpy as np
import pandas as pd

import handygenome.tools as tools
import handygenome.deco as deco
import handygenome.refgenome.refgenome as refgenome
import handygenome.logutils as logutils
import handygenome.workflow as workflow
import handygenome.variant.varianthandler as varianthandler
import handygenome.variant.filter as libfilter
import handygenome.variant.infoformat as infoformat
import handygenome.sv.breakends as libbnd
import handygenome.plot.filterinfo as plot_filterinfo
import handygenome.variant.infoformat as infoformat
import handygenome.vcfeditor.initvcf as initvcf
import handygenome.vcfeditor.headerhandler as headerhandler
#import handygenome.annotation.data as annotdata
import handygenome.annotation.ensembl_feature as ensembl_feature
import handygenome.annotation.ensembl_parser as ensembl_parser
import handygenome.annotation.ensembl_rest as ensembl_rest
import handygenome.annotation.popfreq as libpopfreq
import handygenome.annotation.cosmic as libcosmic
import handygenome.annotation.oncokb as liboncokb
import handygenome.read.readplus as readplus
import handygenome.annotation.readstats as libreadstats
from handygenome.annotation.readstats import AlleleclassError
import handygenome.variant.vcfspec as libvcfspec
import handygenome.vcfeditor.indexing as indexing
import handygenome.cnv.misc as cnvmisc
from handygenome.variant.filter import FilterResultInfo, FilterResultFormat
import handygenome.variant.ponbams as libponbams
import handygenome.vcfeditor.misc as vcfmisc

from handygenome.genomedf.genomedf import GenomeDataFrame
import handygenome.genomedf.genomedf_utils as genomedf_utils


#READCOUNT_FORMAT_KEY = "allele_readcounts"


class VariantPlusInitParams(dict):
    def __init__(self):
        """prop: probability for random selection. If None, all entries are selected.
        """
        self.update({
            'init_popfreq': False,
            'init_cosmic': False,
            'init_transcript': False,
            'init_regulatory': False,
            'init_motif': False,
            'init_repeat': False,
            'init_readstats': False,
            'init_oncokb': False,
            'init_filterresult': False,

            'sampleid_list': None,
        })

    def update_args(self, init_all_attrs, vp_init_params_arg):
        """'init_all_attrs' takes precedence over 'vp_init_params_arg'"""

        assert set(vp_init_params_arg.keys()).issubset(self.keys()), (
            f'Unavailable keys for VariantPlus initation parameters. Valid keys are: '
            f'{set(self.keys())}'
        )

        self.update(vp_init_params_arg)
        if init_all_attrs:
            for key in (
                'init_popfreq',
                'init_cosmic',
                'init_transcript',
                'init_regulatory',
                'init_motif',
                'init_repeat',
                'init_readstats',
                'init_oncokb',
                'init_filterresult',
            ):
                self[key] = True


class VariantPlus:
    # CONSTRUCTORS
    @classmethod
    def from_vr(
        cls, vr, 
        refver=None, fasta=None, chromdict=None,
        init_all_attrs=False,
        vp_init_params=dict(),
        preset_vp_init_params=None,
        popfreq_metadata=None, 
        cosmic_metadata=None,
    ):
        # set params
        result = cls()

        result.vr = vr
        result.refver = refgenome.infer_refver_vr(result.vr) if (refver is None) else refver
        result.fasta = (
            refgenome.get_fasta(result.refver)
            if fasta is None
            else fasta
        )
        result.chromdict = (
            refgenome.ChromDict.from_fasta(result.fasta) if chromdict is None else chromdict
        )
        result.vcfspec = libvcfspec.Vcfspec.from_vr(result.vr, refver=result.refver)

        # set init params
        result.set_vp_init_params(init_all_attrs, vp_init_params, preset_vp_init_params)

        # init
        result.init_common(popfreq_metadata, cosmic_metadata)

        return result

    @classmethod
    def from_vcfspec(
        cls, vcfspec, 
        chromdict=None, 
        init_all_attrs=False,
        vp_init_params=dict(),
        preset_vp_init_params=None,
        popfreq_metadata=None, 
        cosmic_metadata=None,
    ):
        # set params
        result = cls()

        result.vcfspec = vcfspec
        result.refver = result.vcfspec.refver
        result.fasta = result.vcfspec.fasta
        result.chromdict = (
            refgenome.get_chromdict(result.refver)
            if chromdict is None else 
            chromdict
        )
        result.vr = initvcf.create_vr(chromdict=result.chromdict, vcfspec=result.vcfspec)

        # set init params
        result.set_vp_init_params(init_all_attrs, vp_init_params, preset_vp_init_params)

        # init
        result.init_common(popfreq_metadata, cosmic_metadata)

        return result

    @classmethod
    def from_bnds(cls, bnds, chromdict=None, **kwargs):
        vcfspec = bnds.get_vcfspec_bnd1()
        return cls.from_vcfspec(vcfspec, chromdict, **kwargs)

    def set_vp_init_params(self, init_all_attrs, vp_init_params, preset_vp_init_params):
        if preset_vp_init_params is None:
            self.vp_init_params = VariantPlusInitParams()
            self.vp_init_params.update_args(init_all_attrs, vp_init_params)
        else:
            self.vp_init_params = preset_vp_init_params

    def init_common(self, popfreq_metadata, cosmic_metadata):
        # SV attrs
        self.is_sv = varianthandler.check_SV(self.vr)
        self.set_bnds_attributes()  # is_bnd1, bnds

        # popfreq & cosmic
        self._init_popfreq(
            metadata=popfreq_metadata, 
            init_blank=(not self.vp_init_params['init_popfreq']),
        )
        self._init_cosmic(
            metadata=cosmic_metadata, 
            init_blank=(not self.vp_init_params['init_cosmic']),
        )

        # features
        self._init_transcript(init_blank=(not self.vp_init_params['init_transcript']))

        self._init_transcript(init_blank=(not self.vp_init_params['init_transcript']))
        self._init_regulatory(init_blank=(not self.vp_init_params['init_regulatory']))
        self._init_motif(init_blank=(not self.vp_init_params['init_motif']))
        self._init_repeat(init_blank=(not self.vp_init_params['init_repeat']))

        # readstats, readstats_data, rpplists
        self.rpplist_dict = dict()
        self.readstats_data_dict = dict()
        self._init_readstats(
            init_blank=(not self.vp_init_params['init_readstats']), 
            sampleid_list=self.vp_init_params['sampleid_list'],
        )

        # oncokb
        self._init_oncokb(init_blank=(not self.vp_init_params['init_oncokb']))

        # filter result
        self._init_filterresult(init_blank=(not self.vp_init_params['init_filterresult']))

    def __repr__(self):
        alteration_strings = self.get_alteration_strings(one_letter=True)
        if alteration_strings is None:
            string = str(None)
        else:
            tmp = list()
            for alt_idx, tups in alteration_strings.items():
                gene_alt_concat = ", ".join(
                    " ".join(x[:2]) for x in tups
                )
                tmp.append(f'alt_index {alt_idx}: {gene_alt_concat}')
            string = '; '.join(tmp)

        return f'<VariantPlus(vcfspec={self.vcfspec}, alteration=({string}))>'

    def get_alteration_strings(self, one_letter=True):
        if self.transcript is None:
            return None
        else:
            alteration_strings = {alt_idx: set() for alt_idx in range(len(self.vcfspec.alts))}
            for alt_idx, tr_set in enumerate(self.transcript):
                for tr in tr_set.canon_ovlp.values():
                    hgvsp_genename = tr.get_hgvsp_genename(one_letter=one_letter)
                    if hgvsp_genename is None:
                        hgvsc_genename = tr.get_hgvsc_genename()
                        if hgvsc_genename is None:
                            alteration_strings[alt_idx].add(str(None))
                        else:
                            alteration_strings[alt_idx].add((hgvsc_genename, 'noncoding'))
                    else:
                        alteration_strings[alt_idx].add((hgvsp_genename, 'coding'))

                if (
                    (len(alteration_strings[alt_idx]) > 1) 
                    and (str(None) in alteration_strings[alt_idx])
                ):
                    alteration_strings[alt_idx].remove(str(None))

            alteration_strings = {
                key: set(
                    tuple(x[0].split()) + (x[1],) 
                    for x in val
                )
                for key, val in alteration_strings.items()
            }
                
            return alteration_strings

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

    def get_allele(self, allele_index):
        if allele_index == 0:
            return self.ref
        else:
            alt_index = allele_index - 1
            try:
                return self.alts[alt_index]
            except IndexError:
                return pd.NA

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
            self.vcfspec.components.write(self.vr)

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
                hgvsg_list=list(x.to_hgvsg() for x in self.vcfspec.iter_annotation_forms()), 
                token=token, 
                tumor_type=None, 
                evidence_types=None,
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
            for val in self.readstats_dict.values():
                val.vcfspec = self.vcfspec
                val.chromdict = self.chromdict
                val.fasta = self.fasta

    @property
    def readstats(self):
        """Alias"""
        return self.readstats_dict

    def create_readstats(
        self, bam_dict,
        rpplist_kwargs=dict(),
        alleleinfo_kwargs=dict(),
        **kwargs,
    ):
        """Args:
            bam_dict: keys - sample ID, values - pysam.AlignmentFile
        Makes a ReadStatsSampledict object, keys of which are the same as the
            keys of "bam_dict" parameter.
        """

        self.readstats_dict = libreadstats.ReadStatsSampledict.from_bam_dict(
            bam_dict=bam_dict,
            vcfspec=self.vcfspec, 
            fasta=self.fasta, 
            chromdict=self.chromdict,
            rpplist_kwargs=rpplist_kwargs,
            alleleinfo_kwargs=alleleinfo_kwargs,
            **kwargs,
        )

#        self.readstats_dict = libreadstats.ReadStatsSampledict()
#        for sampleid, bam in bam_dict.items():
#            self.readstats_dict[sampleid] = libreadstats.ReadStats.from_bam(
#                vcfspec=self.vcfspec, 
#                bam=bam, 
#                fasta=self.fasta, 
#                chromdict=self.chromdict,
#                rpplist_kwargs=rpplist_kwargs,
#                alleleinfo_kwargs=alleleinfo_kwargs,
#            )

    def update_readstats(
        self, bam_dict, **kwargs, 
    ):
        if len(self.readstats_dict) == 0:
            self.create_readstats(bam_dict, **kwargs)
        else:
            self.readstats_dict.update_bams(
                bam_dict=bam_dict, **kwargs,
            )

    def make_readstats(
        self, 
        bam, 
        rpplist_kwargs={},
        alleleinfo_kwargs={},
    ):
        return libreadstats.ReadStats.from_bam(
            vcfspec=self.vcfspec, 
            bam=bam, 
            fasta=self.fasta, 
            chromdict=self.chromdict,
            rpplist_kwargs=rpplist_kwargs,
            alleleinfo_kwargs=alleleinfo_kwargs,
        )

    def make_rpplist(
        self,
        bam,
        set_alleleclass=True,
        rpplist_kwargs={},
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
        rpplist_kwargs={},
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
        rpplist_kwargs={},
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

    @deco.get_deco_atleast1d(['sampleids'])
    def get_vafs(self, sampleids, n_allele=2, exclude_other=False):
        result = {sid: list() for sid in sampleids}
        for sid in sampleids:
            readstats = self.readstats_dict[sid]
            for allele_idx in range(n_allele):
                try:
                    vaf = readstats.get_vaf(
                        alleleclass=allele_idx, 
                        exclude_other=exclude_other,
                    )
                except AlleleclassError:
                    vaf = np.nan
                result[sid].append(vaf)
        return result

    def get_vaf(self, sampleid, allele_index=1, exclude_other=False, ndigits=None):
        if isinstance(sampleid, (list, tuple)):
            return [
                self.get_vaf_singlesample(x, allele_index=allele_index, exclude_other=exclude_other, ndigits=ndigits)
                for x in sampleid
            ]
        else:
            return self.get_vaf_singlesample(sampleid, allele_index=allele_index, exclude_other=exclude_other, ndigits=ndigits)

    def get_vaf_singlesample(self, sampleid, allele_index=1, exclude_other=False, ndigits=None):
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
        CNt = int(seg_row.CNt)
        A = np.nan if np.isnan(seg_row.A) else int(seg_row.A)
        B = np.nan if np.isnan(seg_row.B) else int(seg_row.B)
        return CNt, A, B

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
        colorby='FIRST_OF_PAIR_STRAND',
        colorby_alleleclass=False,
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
        if colorby_alleleclass:
            igv.cmd(f"colorBy TAG {readplus.ALLELECLASS_TAG_RPP}")
        else:
            igv.cmd(f"colorBy {colorby}")

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

    # miscellaneous
    def get_gr(self):
        return pr.from_dict(
            {
                "Chromosome": [self.chrom],
                "Start": [self.start0],
                "End": [self.end0],
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

    def get_protein_changes(self, one_letter=True):
        result = set()
        for tr_set in self.transcript:
            for tr in tr_set.canon_ovlp.values():
                hgvsp = tr.get_hgvsp_genename(one_letter=one_letter)
                if hgvsp is not None:
                    result.add(hgvsp)
                
        return result

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
        self, 
        vcf_path=None, vcf=None, 
        refver=None,
        fasta=None, chromdict=None,
        logging_lineno=1000, 
        preview_lineno=15, 
        verbose=True, 
        logger=None,
        init_all_attrs=False,
        vp_init_params=dict(),
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
                raise Exception(f'When "refver" is not set, "vcf_path" or "vcf" must be set.')
            self.refver = refgenome.infer_refver_vcfheader(self.vcf.header)
        else:
            self.refver = refver

        if fasta is None:
            self.fasta = refgenome.get_fasta(self.refver)
        else:
            self.fasta = fasta

        if chromdict is None:
            self.chromdict = refgenome.ChromDict.from_fasta(self.fasta)
        else:
            self.chromdict = chromdict

        # vp initiation parameters
        self.vp_init_params = VariantPlusInitParams()
        self.vp_init_params.update_args(init_all_attrs, vp_init_params)

        # lineno, verbose, logger
        self.logging_lineno = logging_lineno
        self.preview_lineno = preview_lineno
        self.verbose = verbose
        if logger is None:
            self.logger = self.get_logger(self.verbose)
        else:
            self.logger = logger

        # others
        self.is_sorted = False
        self._index_gr = None

    @classmethod
    def concat(cls, vplist_list, **kwargs):
        """Args:
            vplist_list: An iterable of VariantPlusList objects
        """
        vplist_list = list(vplist_list)
        if len(set(x.refver for x in vplist_list)) != 1:
            raise Exception(f'Different refvers of input VariantPlusList objects')

        result = cls(refver=vplist_list[0].refver, **kwargs)
        for other in vplist_list:
            result.extend(other)

        return result

    @classmethod
    def from_vps(
        cls, vps, 
        **kwargs
    ):
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
        logging_lineno=1000, 
        preview_lineno=15, 
        verbose=True, 
        logger=None,
        init_all_attrs=False,
        vp_init_params=dict(),

        prop=None,
        chrom=None,
        start0=None,
        end0=None,
    ):
        """Args:
            prop: probability for random selection. If None, all entries are selected.
        """
        result = cls(
            vcf_path=vcf_path, 
            logging_lineno=logging_lineno, 
            preview_lineno=preview_lineno, 
            verbose=verbose, 
            logger=logger,
            init_all_attrs=init_all_attrs,
            vp_init_params=vp_init_params,
        )

        result.load_vps_from_vcf(
            prop=prop,
            chrom=chrom,
            start0=start0,
            end0=end0,
        )

        return result

    @classmethod
    def from_vcf_lazy(
        cls,
        vcf_path,
        logging_lineno=1000, 
        preview_lineno=15, 
        verbose=True, 
        logger=None,
        init_all_attrs=False,
        vp_init_params=dict(),
    ):
        result = cls(
            vcf_path=vcf_path, 
            logging_lineno=logging_lineno, 
            preview_lineno=preview_lineno, 
            verbose=verbose, 
            logger=logger,
            init_all_attrs=init_all_attrs,
            vp_init_params=vp_init_params,
        )
        return result

    def load_vps_from_vcf(
        self,
        prop=None,
        chrom=None,
        start0=None,
        end0=None,
        vpfilter=None,
    ):
        # set other params
        if self.vp_init_params['init_popfreq']:
            popfreq_metadata = libpopfreq.PopfreqMetadata.from_vcfheader(self.vcf.header)
        else:
            popfreq_metadata = None

        if self.vp_init_params['init_cosmic']:
            cosmic_metadata = libcosmic.CosmicMetadata.from_vcfheader(self.vcf.header)
        else:
            cosmic_metadata = None

        # make iterator 
        vp_iterator = self.get_vp_iter_from_vcf(
            chrom=chrom, 
            start0=start0, 
            end0=end0, 
            prop=prop,
        )

        # run iterator
        if vpfilter is None:
            self.extend(vp_iterator)
        else:
            self.extend(filter(vpfilter, vp_iterator))

        # set sorted flag
        if self.vcf is not None:
            if self.vcf.index is not None:
                self.is_sorted = True

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
#            for vr_subiter in tools.chunk_iter(vr_iterator(result.vcf), result.logging_lineno):
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

    #######################
    # initializer helpers #
    #######################

#    def _run_vcf_fetch(self, chrom, start0, end0, respect_refregion=False):
#        assert isinstance(self.vcf, pysam.VariantFile)
#        assert (start0 is None) == (end0 is None)
#        
#        if chrom is None:
#            fetcher = self.vcf.fetch()
#        else:
#            if ((start0 is None) or (end0 is None)):
#                fetcher = self.vcf.fetch(contig=chrom)
#            else:
#                fetcher = self.vcf.fetch(contig=chrom, start=start0, stop=end0)
#
#        if start0 is None:
#            def check_vr_inclusion(vr, chrom, start0, end0):
#                return vr.contig == chrom
#        else:
#            if respect_refregion:
#                def check_vr_inclusion(vr, chrom, start0, end0):
#                    return (
#                        (vr.contig == chrom)
#                        and ((vr.pos - 1) >= start0)
#                        and (vr.pos <= end0)
#                    )
#            else:
#                def check_vr_inclusion(vr, chrom, start0, end0):
#                    return (
#                        (vr.contig == chrom)
#                        and ((vr.pos - 1) >= start0)
#                        and (vr.pos <= end0)
#                    )
#
#        for vr in fetcher:
#            if check_vr_inclusion(vr, chrom, start0, end0):
#                yield vr

    def get_vr_iter_from_vcf(self, chrom=None, start0=None, end0=None, prop=None):
        assert isinstance(self.vcf, pysam.VariantFile)

        #fetcher = self._run_vcf_fetch(chrom, start0, end0)
        #fetcher = vcfmisc.get_vr_fetcher(
        fetcher = vcfmisc.get_vr_fetcher_new(
            self.vcf, self.refver, chrom, start0, end0,
        )

        if prop is None:
            vr_iterator = fetcher
        else:
            vr_iterator = tools.bernoulli_iterator(fetcher, p=prop, block_size=int(1e5))

        return vr_iterator

    def _get_vp_iter_from_vr_iter(self, vr_iterator):
        for vr in vr_iterator:
            if varianthandler.check_SV(vr):
                vr_svinfo = libbnd.get_vr_svinfo_standard_vr(vr, self.fasta, self.chromdict)
                if not vr_svinfo["is_bnd1"]:
                    continue

            vp = VariantPlus.from_vr(
                vr=vr,
                refver=self.refver,
                fasta=self.fasta,
                chromdict=self.chromdict,
                preset_vp_init_params=self.vp_init_params,
            )

            yield vp

    def get_vp_iter_from_vcf(
        self, 
        *,
        chrom=None, start0=None, end0=None, 
        prop=None,
        vpfilter=None,
    ):
        # "prop" is treated in VariantRecord iteration step for performance
        vr_iterator = self.get_vr_iter_from_vcf(
            chrom, start0, end0, prop
        )
        vp_iterator = self._get_vp_iter_from_vr_iter(vr_iterator)
        vp_iterator = filter(vpfilter, vp_iterator)
        if self.logging_lineno is not None:
            vp_iterator = logutils.iter_lineno_logging(vp_iterator, self.logger, self.logging_lineno)
        return vp_iterator

    def get_vp_iter_from_self(
        self, 
        chrom=None, start0=None, end0=None, 
        prop=None,
        vpfilter=None,
    ):
        # filter by coordinates
        if chrom is None:
            vp_iterator = iter(self)
        else:
            # self gets sorted
            coord_indexes = self.get_coord_indexes(chrom, start0, end0)
            if coord_indexes is None:
                vp_iterator = iter(())
            else:
                sl = slice(coord_indexes.iloc[0], coord_indexes.iloc[-1] + 1)
                vp_iterator = itertools.islice(self, sl.start, sl.stop, sl.step)
        # filter by prop
        if prop is not None:
            vp_iterator = tools.bernoulli_iterator(vp_iterator, p=prop, block_size=int(1e5))
        # filter by vpfilter
        vp_iterator = filter(vpfilter, vp_iterator)

        if self.logging_lineno is not None:
            vp_iterator = logutils.iter_lineno_logging(vp_iterator, self.logger, self.logging_lineno)

        return vp_iterator


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

    ###########################################
    # show readstats summary for all variants #
    ###########################################

    def get_all_vp_readstats(self, key, sampleid, allele_index=1):
        return np.fromiter(
            (
                vp.readstats_dict[sampleid][key][allele_index] 
                for vp in self
            ),
            dtype=float,
        )

    def get_all_vp_readstats_diff(self, key, sampleid, allele_indexes=(1, 0)):
        return np.fromiter(
            (
                (
                    vp.readstats_dict[sampleid][key][allele_indexes[0]] 
                    - vp.readstats_dict[sampleid][key][allele_indexes[1]] 
                )
                for vp in self
            ),
            dtype=float,
        )

    def get_all_vp_readstats_average(self, key, sampleid):
        return np.fromiter(
            (
                vp.readstats_dict[sampleid].get_allele_indexes_average(
                    key, allele_indexes=None,
                ) 
                for vp in self
            ),
            dtype=float,
        )

    def get_all_vp_total_depths(self, sampleid, exclude_other=False):
        return np.fromiter(
            (
                vp.readstats_dict[sampleid].get_total_rppcount(exclude_other=exclude_other)
                for vp in self
            ),
            dtype=float,
        )

    def get_all_vp_vafs(self, sampleid, allele_index=1, exclude_other=False):
        return np.fromiter(
            (
                vp.get_vaf(sampleid, allele_index=allele_index, exclude_other=exclude_other)
                for vp in self
            ),
            dtype=float,
        )

    get_vafs = get_all_vp_vafs

    #################
    # set operation #
    #################

    def difference(self, other):
        left_vcfspecs = set(x.vcfspec for x in self)
        right_vcfspecs = set(x.vcfspec for x in other)
        diff = left_vcfspecs.difference(right_vcfspecs)

        result = self.__class__(refver=self.refver)
        for vp in self:
            if vp.vcfspec in diff:
                result.append(vp)

        return result

    ##########
    # others #
    ##########

    def set_unsorted(self):
        self.is_sorted = False
        self._index_gr = None

    def append(self, vp):
        #super().append(vp)
        list.append(self, vp)
        self.set_unsorted()

    def extend(self, vp_iter):
        #super().extend(vp_iter)
        list.extend(self, vp_iter)
        self.set_unsorted()

    def get_index_gr(self):
        if not self.is_sorted:
            self.sort()

        gr = self.get_df(vaf_sampleid=None, as_gr=True)
        gr.Vcfspec = self
        gr.Index = list(range(gr.df.shape[0]))
            # In this form of assignment, right-side term must be a list.

        return gr

    #@functools.cached_property
    @property
    def index_gr(self):
        if self._index_gr is None:
            self._index_gr = self.get_index_gr()
        return self._index_gr

    def fetch_by_vcfspec(self, vcfspec):
        selector = (self.index_gr.Vcfspec == vcfspec)
        #if not selector.any():
        #    raise ValueError(f'Input vcfspec is not present in the VariantPlusList.')
        #idx = self.index_gr.Vcfspec.loc[selector].index[0]
        #idx = self.index_gr.Index.loc[selector][0]
        indexes = np.where(selector)[0]
        if len(indexes) == 0:
            raise ValueError(f'Input Vcfspec is not present in the VariantPlusList.')
        elif len(indexes) > 1:
            raise Exception(f'Duplicate Vcfspec in the VariantPlusList')

        return self[indexes[0]]

    def get_coord_indexes(self, chrom, start0, end0):
        if (start0 is None) and (end0 is None):
            gr_subset = self.index_gr[chrom]
        else:
            if start0 is None:
                start0 = 0
            if end0 is None:
                end0 = self.chromdict[chrom]
            gr_subset = self.index_gr[chrom, start0:end0]

        if gr_subset.df.shape[0] == 0:
            return None
        else:
            return gr_subset.Index

    def fetch_by_coord(self, chrom, start0=None, end0=None):
        coord_indexes = self.get_coord_indexes(chrom, start0, end0)
        if coord_indexes is None:
            #return self.spawn()
            return iter(tuple())
        else:
            #sl = slice(coord_indexes.iloc[0], coord_indexes.iloc[-1] + 1)
            #return self[sl]
            return itertools.slice(self, coord_indexes.iloc[0], coord_indexes.iloc[-1] + 1)

    def sample(self, n=1):
        if len(self) == 0:
            raise Exception(f'The VariantPlusList is empty.')

        if n == 1:
            return random.choice(self)
        else:
            return random.sample(self, k=n)

#    def get_df_multiprocessing(
#        self, 
#        alt_index=0,
#        vaf_sampleid=None, 
#        as_gr=False, 
#        get_vaf_kwargs={
#            'exclude_other': False, 
#            'ndigits': None,
#        },
#        lazy=False,
#        vcf_iter_kwargs=dict(),
#        nproc=1,
#    ):
#        # parameter setups
#        get_vaf_kwargs['allele_index'] = alt_index + 1
#
#        if vaf_sampleid is None:
#            if lazy:
#                first_vr = next(self.get_vr_iter_from_vcf(**vcf_iter_kwargs))
#                sample_ids = list(first_vr.header.samples)
#            else:
#                sample_ids = list(first_vp.vr.header.samples)
#
#            if len(sample_ids) == 0:
#                vaf_sampleid = None
#                #omit_vaf = True
#            else:
#                vaf_sampleid = sample_ids
#        else:
#            if not isinstance(vaf_sampleid, (tuple, list)):
#                vaf_sampleid = [vaf_sampleid]
#
#        # run multiprocessing
#        if lazy:
#            vp_iterator = self.get_vp_iter_from_vcf(**vcf_iter_kwargs)
#        else:
#            vp_iterator = iter(self)
#
#        with multiprocessing.Pool(nproc) as pool:
#            source_data = pool.starmap(
#                vplist_get_df_subproc, 
#                (
#                    (vp, alt_index, vaf_sampleid, get_vaf_kwargs) 
#                    for vp in vp_iterator
#                ),
#            )
#
#        # make data
#        columns = ['Chromosome', 'Start', 'End', 'REF', 'ALT']
#        if vaf_sampleid is not None:
#            columns.extend(f'vaf_{sid}' for sid in vaf_sampleid)
#        data = dict(zip(columns, zip(*source_data)))
#
#        # result
#        if as_gr:
#            return pr.from_dict(data)
#        else:
#            return pd.DataFrame.from_dict(data)

#    def get_df(
#        self, 
#        alt_index=0,
#        vaf_sampleid=None, 
#        as_gr=False, 
#        #omit_vaf=False, 
#        lazy=False,
#        get_vaf_kwargs=dict(),
#        vcf_iter_kwargs=dict(),
#    ):
#        # parameter setups
#        get_vaf_kwargs['allele_index'] = alt_index + 1
#
#        if vaf_sampleid is None:
#            if lazy:
#                first_vr = next(self.get_vr_iter_from_vcf(**vcf_iter_kwargs))
#                sample_ids = list(first_vr.header.samples)
#            else:
#                sample_ids = list(self[0].vr.header.samples)
#
#            if len(sample_ids) == 0:
#                vaf_sampleid = None
#                #omit_vaf = True
#            else:
#                vaf_sampleid = sample_ids
#        else:
#            if not isinstance(vaf_sampleid, (tuple, list)):
#                vaf_sampleid = [vaf_sampleid]
#
#        # main
#        chroms = list()
#        pos1s = list()
#        starts = list()
#        ends = list()
#        refs = list()
#        alts = list()
#
#        if vaf_sampleid is not None:
#            vafs = [list() for x in vaf_sampleid]
#
#        if lazy:
#            vp_iterator = self.get_vp_iter_from_vcf(**vcf_iter_kwargs)
#        else:
#            vp_iterator = iter(self)
#
#        for vp in vp_iterator:
#            chroms.append(vp.vcfspec.chrom)
#            pos1s.append(vp.vcfspec.pos)
#            starts.append(vp.vcfspec.pos0)
#            ends.append(vp.vcfspec.end0)
#            refs.append(vp.vcfspec.ref)
#            alts.append(vp.vcfspec.alts[alt_index])
#
#            if vaf_sampleid is not None:
#                for vaf_val, sublist in zip(
#                    vp.get_vaf(vaf_sampleid, **get_vaf_kwargs),
#                    vafs,
#                ):
#                    sublist.append(vaf_val)
#
#        # make data
#        data = {
#            'Chromosome': chroms,
#            'Start': starts,
#            'End': ends,
#            'POS': pos1s,
#            'REF': refs,
#            'ALT': alts,
#        }
#        if vaf_sampleid is not None:
#            for sid, vaf_sublist in zip(vaf_sampleid, vafs):
#                data.update({f'vaf_{sid}': vaf_sublist})
#
#        # result
#        if as_gr:
#            return pr.from_dict(data)
#        else:
#            return pd.DataFrame.from_dict(data)

    def get_df(
        self, 
        chrom=None, start0=None, end0=None,
        alt_index=0,
        vaf_sampleid=None, 
        as_gr=False, 
        lazy=False,
        get_vaf_kwargs=dict(),
        vcf_iter_kwargs=dict(),
    ):
        # parameter setups
        get_vaf_kwargs['allele_index'] = alt_index + 1

        if vaf_sampleid is None:
            vaf_sampleid = None
#            if lazy:
#                first_vr = next(self.get_vr_iter_from_vcf(**vcf_iter_kwargs))
#                sample_ids = list(first_vr.header.samples)
#            else:
#                sample_ids = list(self[0].vr.header.samples)
#
#            if len(sample_ids) == 0:
#                vaf_sampleid = None
#            else:
#                vaf_sampleid = sample_ids
        else:
            if not isinstance(vaf_sampleid, (tuple, list)):
                vaf_sampleid = [vaf_sampleid]

        # main
        chroms = list()
        pos1s = list()
        starts = list()
        ends = list()
        refs = list()
        alts = list()

        if vaf_sampleid is not None:
            vafs = [list() for x in vaf_sampleid]

        if lazy:
            vp_iterator = self.get_vp_iter_from_vcf(
                chrom=chrom, start0=start0, end0=end0, 
                **vcf_iter_kwargs,
            )
        else:
            if chrom is None:
                vp_iterator = iter(self)
            else:
                vp_iterator = self.fetch_by_coord(chrom, start0, end0)

        for vp in vp_iterator:
            chroms.append(vp.vcfspec.chrom)
            pos1s.append(vp.vcfspec.pos)
            starts.append(vp.vcfspec.pos0)
            ends.append(vp.vcfspec.end0)
            refs.append(vp.vcfspec.ref)
            alts.append(vp.vcfspec.alts[alt_index])

            if vaf_sampleid is not None:
                for vaf_val, sublist in zip(
                    vp.get_vaf(vaf_sampleid, **get_vaf_kwargs),
                    vafs,
                ):
                    sublist.append(vaf_val)

        # make data
        data = {
            'Chromosome': chroms,
            'Start': starts,
            'End': ends,
            'POS': pos1s,
            'REF': refs,
            'ALT': alts,
        }
        if vaf_sampleid is not None:
            for sid, vaf_sublist in zip(vaf_sampleid, vafs):
                data.update({f'vaf_{sid}': vaf_sublist})

        # result
        if as_gr:
            return pr.from_dict(data)
        else:
            return pd.DataFrame.from_dict(data).astype(
                {'Start': int, 'End': int, 'POS': int}
            )
    
    @staticmethod
    def vafdf_allele_columns(num):
        return ['REF'] + [f'ALT{x + 1}' for x in range(num - 1)]

    @deco.get_deco_atleast1d(['sampleids'])
    def get_vafdf(
        self, 
        sampleids, 
        chrom=None, start0=None, end0=None,
        n_allele=2,
        exclude_other=False,

        lazy=False,
        prop=None,
        vpfilter=None,

        multiindex=False,
        sampleid_index_name='sampleid',
    ):
        # set data component lists
        chroms = list()
        if multiindex:
            pos1s = list()
        starts = list()
        ends = list()

        alleles = [list() for x in range(n_allele)]
        vafs = {
            sid: [list() for x in range(n_allele)] 
            for sid in sampleids
        }

        # set vp iterator
        if lazy:
            vp_iterator = self.get_vp_iter_from_vcf(
                chrom=chrom, start0=start0, end0=end0, 
                prop=prop,
                vpfilter=vpfilter,
            )
        else:
            if chrom is None:
                vp_iterator = iter(self)
            else:
                vp_iterator = self.fetch_by_coord(chrom, start0, end0)

        # extract information for vps
        for vp in vp_iterator:
            this_vafs = vp.get_vafs(sampleids, n_allele=n_allele, exclude_other=exclude_other)
            # skip if vafs are invalid values
            vafs_are_valid = all(
                (not np.isnan(sublist).any())
                for sid, sublist in this_vafs.items()
            )
            if not vafs_are_valid:
                continue

            # populate data
            chroms.append(vp.vcfspec.chrom)
            if multiindex:
                pos1s.append(vp.vcfspec.pos)
            starts.append(vp.vcfspec.pos0)
            ends.append(vp.vcfspec.end0)

            for allele_idx, sublist in enumerate(alleles):
                sublist.append(vp.get_allele(allele_idx))

            for sid, this_vaf_list in this_vafs.items():
                for vaflist, vafval in zip(vafs[sid], this_vaf_list):
                    vaflist.append(vafval)

        # result
        if multiindex:
            data = list()
            columns_lv2 = list()

            data.append(chroms)
            columns_lv2.append('Chromosome')
            data.append(starts)
            columns_lv2.append('Start')
            data.append(ends)
            columns_lv2.append('End')
            data.append(pos1s)
            columns_lv2.append('POS')

            allele_colnames = self.vafdf_allele_columns(n_allele)
            for colname, sublist in zip(allele_colnames, alleles):
                data.append(sublist)
                columns_lv2.append(colname)

            for sid in sampleids:
                for allele_col, vaflist in zip(allele_colnames, vafs[sid]):
                    data.append(vaflist)
                    columns_lv2.append(f'{allele_col}_vaf')

            columns = pd.MultiIndex.from_arrays(
                [
                    np.concatenate([np.repeat(None, 4 + n_allele), np.repeat(sampleids, n_allele)]),
                    columns_lv2,
                ],
                names=[sampleid_index_name, None],
            )

            return pd.DataFrame(list(zip(*data)), columns=columns)
        else:
            data = dict()

            data['chroms'] = chroms
            data['start0s'] = starts
            data['end0s'] = ends

            allele_colnames = self.vafdf_allele_columns(n_allele)
            for colname, sublist in zip(allele_colnames, alleles):
                data[colname] = sublist
            for sid in sampleids:
                for allele_col, vaflist in zip(allele_colnames, vafs[sid]):
                    formatkey = f'{allele_col}_vaf'
                    colname = sid + VCFDataFrame.SAMPLEID_SEP + formatkey
                    data[colname] = vaflist

            return VCFDataFrame.from_data(refver=self.refver, **data)

    def get_gr(self, vaf_sampleid=None):
        return self.get_df(vaf_sampleid=vaf_sampleid, as_gr=True)

    def spawn(self):
        kwargs = {
            key: getattr(self, key) for key in 
            (
                'vcf_path', 'vcf', 
                'refver', 'fasta', 'chromdict',
                'logging_lineno', 'preview_lineno', 'verbose', 
            )
        }
        kwargs['vp_init_params'] = self.vp_init_params
        kwargs['init_all_attrs'] = False
        result = self.__class__(**kwargs)
        return result

    def filter(self, key):
        result = self.spawn()
        for vp in self:
            if key(vp):
                result.append(vp)
        return result
    
    def isec(self, gr):
        if len(self) == 0:
            return self
        else:
            overlaps = self.get_gr().count_overlaps(gr, overlap_col="count")
            marks = overlaps.df.loc[:, "count"] > 0
            result = self.spawn()
            result.extend(itertools.compress(self, marks))
            return result

    def get_sigresult(self, catalogue_type="sbs96", **kwargs):
        import handygenome.signature.signatureresult as libsig

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
        vr_sortkey = varianthandler.get_vr_sortkey(self.chromdict)
        def vp_sortkey(vp):
            return vr_sortkey(vp.vr)
        return vp_sortkey

    def sort(self):
        super().sort(key=self.get_vp_sortkey())
        self.is_sorted = True

    ###########################
    # writing-related methods #
    ###########################

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
        if not self.is_sorted:
            self.sort()
        mode_pysam = vcfmisc.write_mode_arghandler(mode_bcftools, mode_pysam)
        with pysam.VariantFile(outfile_path, mode=mode_pysam, header=header) as out_vcf:
            for vp in self.get_vp_iter_from_self():
                if vp.is_sv:
                    self.write_each_vr(vp.vr, out_vcf, conflicting_keys)
                    self.write_each_vr(vp.get_vr_bnd2(), out_vcf, conflicting_keys)
                else:
                    self.write_each_vr(vp.vr, out_vcf, conflicting_keys)
        if index:
            indexing.index_vcf(outfile_path)

    def write_with_filter(
        self, outfile_path, vpfilter, mode_bcftools="z", mode_pysam=None, index=True
    ):
        """Loads VariantRecords from VCF file, filters, and writes"""
        header = self.vcf.header.copy()
        mode_pysam = vcfmisc.write_mode_arghandler(mode_bcftools, mode_pysam)
        with pysam.VariantFile(outfile_path, mode=mode_pysam, header=header) as out_vcf:
            for vp in self.get_vp_iter_from_vcf(vpfilter=vpfilter):
                if vp.is_sv:
                    out_vcf.write(vp.vr)
                    out_vcf.write(vp.get_vr_bnd2())
                else:
                    out_vcf.write(vp.vr)

        if index:
            indexing.index_vcf(outfile_path)


#########
# vcfdf #
#########

class VCFDataFrame(GenomeDataFrame):
    COMMON_COLUMNS = (
        list(genomedf_utils.COMMON_COLUMNS) 
        + ['REF', 'ALT1']
    )
    SAMPLEID_SEP = ':::'

    ##############
    # properties #
    ##############

    @property
    def POS(self):
        return self['Start'] + 1

    @property
    def nonannot_columns(self):
        return self.colname_getter(self.colnamefilter_nonannot)

    @property
    def allele_columns(self):
        return [
            x for x in self.nonannot_columns 
            if x not in genomedf_utils.COMMON_COLUMNS
        ]

    @property
    def info_columns(self):
        return self.colname_getter(self.colnamefilter_info)

    @property
    def format_columns(self):
        return self.colname_getter(self.colnamefilter_format)

    @property
    def samples(self):
        return tools.unique_keeporder(
            [
                x.split(self.__class__.SAMPLEID_SEP)[0] 
                for x in self.format_columns
            ]
        )

    @property
    def format_keys(self):
        return tools.unique_keeporder(
            [
                x.split(self.__class__.SAMPLEID_SEP)[1] 
                for x in self.format_columns
            ]
        )

    #######################
    # column name filters #
    #######################

    def colname_getter(self, colnamefilter):
        return tools.unique_keeporder(filter(colnamefilter, self.columns))

    @classmethod
    def colnamefilter_nonannot(cls, x):
        return (
            (x in cls.COMMON_COLUMNS)
            or (re.fullmatch('ALT[0-9]+', x) is not None)
        )

    @classmethod
    def colnamefilter_info(cls, x):
        return (
            (cls.SAMPLEID_SEP not in x)
            and (not cls.colnamefilter_nonannot(x))
        )

    @classmethod
    def colnamefilter_format(cls, x):
        return (
            (cls.SAMPLEID_SEP in x)
            and (not cls.colnamefilter_nonannot(x))
        )

    ###############
    # sanitycheck #
    ###############
        
    @classmethod
    def sanitycheck_df(cls, df):
        super().sanitycheck_df(df)
        assert all(
            (cls.SAMPLEID_SEP not in x)
            for x in df.columns
            if cls.colnamefilter_nonannot(x)
        )
        assert all(
            (x.count(cls.SAMPLEID_SEP) == 1)
            for x in df.columns
            if cls.colnamefilter_format(x)
        )

    ###########
    # getters #
    ###########

    def get_format(self, samples, keys):
        samples = np.atleast_1d(samples)
        keys = np.atleast_1d(keys)
        selected_colnames = list()
        for s, k in itertools.product(samples, keys):
            colname = s + self.__class__.SAMPLEID_SEP + k
            selected_colnames.append(colname)

        return self.df.loc[:, selected_colnames]

    def get_sample_annots(self, sample):
        colnames = [
            x for x in self.format_columns 
            if x.split(self.__class__.SAMPLEID_SEP)[0] == sample
        ]
        return self.loc[:, colnames]


#########################################
# parallelized vaf dataframe generation #
#########################################

@deco.get_deco_atleast1d(['sampleids'])
def get_vafdf(
    vcf_path, 
    sampleids, 
    n_allele=2,
    nproc=1,
    exclude_other=False,
    prop=None,
    vpfilter=None,
    verbose=True,
):
    """Columns: MultiIndex with 2 levels
        level 1: (name:                                              None |                                 sid1 |                                 sid2
        level 2: Chromosome, Start, End, REF, [ALT1, [ALT2, ...]]   REF_vaf, [ALT1_vaf, [ALT2_vaf, ...]]   REF_vaf, [ALT1_vaf, [ALT2_vaf, ...]]
    """
    # get VCF fetch regions for each parallel job
    if verbose:
        logutils.log(f'Extracting vcf position information') 
    fetchregion_gdf_list = vcfmisc.get_vcf_fetchregions_new(vcf_path, n=nproc, nproc=nproc)

    # run multiprocess jobs
    args = (
        (
            fetchregion_gdf,
            vcf_path, 
            sampleids, 
            n_allele, 
            exclude_other,
            prop,
            vpfilter,
        )
        for fetchregion_gdf in fetchregion_gdf_list
    )
    with multiprocessing.Pool(nproc) as pool:
        if verbose:
            logutils.log(f'Running parallel jobs') 
        mp_result = pool.starmap(get_vafdf_targetfunc, args)
        if verbose:
            logutils.log(f'Concatenating split job dataframes') 
        result = VCFDataFrame.concat(itertools.chain.from_iterable(mp_result))

    return result


def get_vafdf_targetfunc(
    fetchregion_gdf,
    vcf_path, 
    sampleids, 
    n_allele, 
    exclude_other,
    prop,
    vpfilter,
):
    vafdf_list = list()
    for chrom, start0, end0 in fetchregion_gdf.iter_coords():
        vplist = VariantPlusList.from_vcf_lazy(
            vcf_path, 
            logging_lineno=None, 
            verbose=False,
            init_all_attrs=False,
            vp_init_params=dict(
                init_readstats=True,
                sampleid_list=sampleids,
            ),
        )
        vafdf = vplist.get_vafdf(
            sampleids=sampleids,
            chrom=chrom, start0=start0, end0=end0,
            n_allele=n_allele,
            exclude_other=exclude_other,
            lazy=True,
            prop=prop,
            vpfilter=vpfilter,
        )
        vafdf_list.append(vafdf)

    return vafdf_list


#def get_vafdf_targetfunc_old(
#    position_info, 
#    refver, 
#    vcf_path, 
#    sampleids, 
#    n_allele, 
#    exclude_other,
#    prop,
#    vpfilter,
#):
#    fetchargs_list = vcfmisc.get_fetchargs_from_vcf_positions(position_info, refver)
#    dflist = list()
#    for fetchargs in fetchargs_list:
#        df = get_vafdf_nonparallel(
#            vcf_path, 
#            sampleids, 
#            chrom=fetchargs[0], start0=fetchargs[1], end0=fetchargs[2],
#            n_allele=n_allele,
#
#            exclude_other=exclude_other,
#            prop=prop,
#            vpfilter=vpfilter,
#        )
#        dflist.append(df)
#
#    return dflist
#
#
#@deco.get_deco_atleast1d(['sampleids'])
#def get_vafdf_nonparallel(
#    vcf_path, 
#    sampleids, 
#    chrom=None, start0=None, end0=None,
#    n_allele=2, 
#    exclude_other=False,
#    prop=None,
#    vpfilter=None,
#):
#    vplist = VariantPlusList.from_vcf_lazy(
#        vcf_path, 
#        logging_lineno=None, 
#        verbose=False,
#        init_all_attrs=False,
#        vp_init_params=dict(
#            init_readstats=True,
#            sampleid_list=sampleids,
#        ),
#    )
#    vafdf = vplist.get_vafdf(
#        sampleids=sampleids,
#        chrom=chrom, start0=start0, end0=end0,
#        n_allele=n_allele,
#        exclude_other=exclude_other,
#        lazy=True,
#        prop=prop,
#        vpfilter=vpfilter,
#    )
#    return vafdf


#def get_vafdf_old(
#    vcf_path, 
#    sampleid, 
#    nproc=1,
#    alt_index=0, 
#    get_vaf_kwargs=dict(),
#    vcf_iter_kwargs=dict(),
#    verbose=True,
#):
#    # setup params
#    with pysam.VariantFile(vcf_path) as vcf:
#        refver = refgenome.infer_refver_vcfheader(vcf.header)
#
#    logutils.log(f'Extracting vcf position information') 
#    all_position_info = vcfmisc.get_vcf_positions(
#        vcf_path, as_iter=False, verbose=False,
#    )
#    split_position_info = [
#        x for x in np.array_split(all_position_info, nproc) 
#        if x.shape[0] != 0
#    ]
#
#    # run multiprocess jobs
#    logutils.log(f'Running parallel jobs') 
#    args = (
#        (position_info, refver, vcf_path, sampleid, alt_index, get_vaf_kwargs, vcf_iter_kwargs)
#        for position_info in split_position_info
#    )
#    with multiprocessing.Pool(nproc) as pool:
#        mp_result = pool.starmap(_get_vafdf_targetfunc_old, args)
#        logutils.log(f'Concatenating split job dataframes') 
#        result = pd.concat(itertools.chain.from_iterable(mp_result), axis=0)
#
#    result.reset_index(drop=True, inplace=True)
#
#    return result
#
#
#def _get_vafdf_targetfunc_old(
#    position_info, refver, vcf_path, sampleid, alt_index, get_vaf_kwargs, vcf_iter_kwargs,
#):
#    fetchargs_list = vcfmisc.get_fetchargs_from_vcf_positions(position_info, refver)
#    dflist = list()
#    for fetchargs in fetchargs_list:
#        df = get_vafdf_nonparallel_old(
#            vcf_path, sampleid, 
#            chrom=fetchargs[0], start0=fetchargs[1], end0=fetchargs[2],
#            alt_index=alt_index, 
#            as_gr=False, 
#            get_vaf_kwargs=get_vaf_kwargs,
#            vcf_iter_kwargs=vcf_iter_kwargs,
#            logging_lineno=1000, 
#            verbose=False,
#        )
#        dflist.append(df)
#
#    return dflist
#
#
#def get_vafdf_nonparallel_old(
#    vcf_path, sampleid, 
#    chrom=None, start0=None, end0=None,
#    alt_index=0, 
#    as_gr=False, 
#    get_vaf_kwargs=dict(),
#    vcf_iter_kwargs=dict(),
#    logging_lineno=1000, 
#    verbose=True,
#):
#    sampleid = list(np.atleast_1d(sampleid))
#    vplist = VariantPlusList.from_vcf_lazy(
#        vcf_path, 
#        logging_lineno=logging_lineno, 
#        verbose=verbose,
#        init_all_attrs=False,
#        vp_init_params=dict(
#            init_readstats=True,
#            sampleid_list=sampleid,
#        ),
#    )
#    vafdf = vplist.get_df(
#        chrom=chrom, start0=start0, end0=end0,
#        alt_index=alt_index, 
#        vaf_sampleid=sampleid,
#        as_gr=as_gr,
#        lazy=True,
#        get_vaf_kwargs=get_vaf_kwargs,
#        vcf_iter_kwargs=vcf_iter_kwargs,
#    )
#    return vafdf


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



import re
import pprint
import os
import itertools
import random

import pysam
import pyranges
import numpy as np

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
breakends_module = importlib.import_module('.'.join([top_package_name, 'variantplus', 'breakends']))
varianthandler = importlib.import_module('.'.join([top_package_name, 'variantplus', 'varianthandler']))
variantviz = importlib.import_module('.'.join([top_package_name, 'variantplus', 'variantviz']))
vpfilter = importlib.import_module('.'.join([top_package_name, 'variantplus', 'vpfilter']))
infoformat = importlib.import_module('.'.join([top_package_name, 'variantplus', 'infoformat']))
initvcf = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'initvcf']))
annotationdb = importlib.import_module('.'.join([top_package_name, 'annotation', 'annotationdb']))
readplus = importlib.import_module('.'.join([top_package_name, 'readplus', 'readplus']))
libreadstats = importlib.import_module('.'.join([top_package_name, 'annotation', 'readstats']))
alleleinfosetup = importlib.import_module('.'.join([top_package_name, 'readplus', 'alleleinfosetup']))
alleleinfosetup_sv = importlib.import_module('.'.join([top_package_name, 'readplus', 'alleleinfosetup_sv']))


READCOUNT_FORMAT_KEY = 'allele_readcounts'


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

    @common.get_deco_num_set(('vr', 'vcfspec', 'bnds'), 1)
    def __init__(self, vr=None, vcfspec=None, bnds=None,
                 refver=None, fasta=None, chromdict=None,
                 add_annotdb_infometas=True, set_annotdb=True,
                 set_readstats=True, 
                 annotitem_septype=annotationdb.DEFAULT_SEPTYPE):
        """Args:
            vr: pysam.VariantRecord instance
            fasta: pysam.FastaFile instance
            chromdict: handygenome.common.Chromdict instance
        """

        def init_from_vr(self, vr, refver, fasta, chromdict):
            self.vr = vr
            self.refver = (common.infer_refver_vr(self.vr)
                           if refver is None else
                           refver)
            self.fasta = (
                pysam.FastaFile(common.DEFAULT_FASTA_PATHS[self.refver])
                if fasta is None else
                fasta)
            self.chromdict = (common.ChromDict(fasta=self.fasta)
                              if chromdict is None else
                              chromdict)
            self.vcfspec = varianthandler.get_vcfspec(self.vr)

        def init_from_vcfspec(self, vcfspec, refver, fasta, chromdict):
            if refver is None:
                raise Exception(
                    f'If initializing with vcfspec, refver must be given.')

            self.vcfspec = vcfspec
            self.refver = refver

            self.fasta = (
                pysam.FastaFile(common.DEFAULT_FASTA_PATHS[self.refver])
                if fasta is None else
                fasta)
            self.chromdict = (common.ChromDict(fasta=self.fasta)
                              if chromdict is None else
                              chromdict)

            self.vr = initvcf.create_vr(chromdict=self.chromdict, 
                                        vcfspec=self.vcfspec)

        def init_from_bnds(self, bnds, refver, fasta, chromdict):
            if refver is None:
                raise Exception(
                    f'If initializing with Breakends, refver must be given.')

            vcfspec = bnds.get_vcfspec_bnd1()
            init_from_vcfspec(self, vcfspec, refver, fasta, chromdict)

        # MAIN
        #assert len(vr.alts) == 1, (
        #    f'Multiallelic variant record is not allowed:\n{vr}')
        if vr is not None:
            init_from_vr(self, vr, refver, fasta, chromdict)
        elif vcfspec is not None:
            init_from_vcfspec(self, vcfspec, refver, fasta, chromdict)
        elif bnds is not None:
            init_from_bnds(self, bnds, refver, fasta, chromdict)

        # SV attrs
        self.is_sv = varianthandler.check_SV(self.vr)
        self.set_bnds_attributes()  # is_bnd1, bnds

        # annotdb
        self.annotitem_septype = annotitem_septype
        if add_annotdb_infometas:
            annotationdb.add_infometas(self.vr.header)
        if set_annotdb:
            self.set_annotdb()

        # readstats_datas, rpplists
        self.rpplist_dict = dict()
        self.readstats_datas = dict()
        self.readstats_dict = dict()
        if set_readstats:
            self.load_readstats()

        # sampleids
        self.sampleids = tuple(self.vr.header.samples)

    def __repr__(self):
        vr_string = '\t'.join(str(self.vr).split('\t')[:5])
        return f'<VariantPlus object ({vr_string})>'

    def set_readcounts(self):
        pass

    def set_annotdb(self):
        if self.is_sv:
            self.annotdb = None
            self.annotdb_bnd1 = annotationdb.AnnotDB(
                'bnd1', self.refver, self.fasta, self.chromdict,
                vr=self.vr, septype=self.annotitem_septype)
            self.annotdb_bnd2 = annotationdb.AnnotDB(
                'bnd2', self.refver, self.fasta, self.chromdict,
                vr=self.vr, septype=self.annotitem_septype)
        else:
            self.annotdb = annotationdb.AnnotDB(
                'plain', self.refver, self.fasta, self.chromdict, vr=self.vr,
                septype=self.annotitem_septype)
            self.annotdb_bnd1 = None
            self.annotdb_bnd2 = None

    ##################################################################

    def get_other_alleleindexes(self, alleleindex):
        """Returns all integer allele indexes other than the input"""

        all_alleleindexes = set(range(len(self.vr.alleles)))
        other_alleleindexes = all_alleleindexes.difference([alleleindex])
        return tuple(sorted(other_alleleindexes))

    ##################################################################

    def get_gr(self):
        return pyranges.from_dict({'Chromosome': [self.vr.contig],
                                   'Start': [self.vr.pos - 1],
                                   'End': [self.vr.pos]})

    def get_gene_names(self, canonical_only=True):
        if canonical_only:
            return [feature['gene_name'] 
                    for feature in 
                    self.annotdb.transcript_canon_ovlp.values()]
        else:
            return [feature['gene_name'] 
                    for feature in 
                    self.annotdb.transcript_ovlp.values()]

    def check_intergenic(self):
        return len(self.annotdb.transcript) == 0

    def get_info(self, key, collapse_tuple=True):
        return infoformat.get_value_info(self.vr, key, 
                                         collapse_tuple=collapse_tuple)

    def get_format(self, sampleid, key, collapse_tuple=True):
        return infoformat.get_value_format(self.vr, sampleid, key,
                                           collapse_tuple=collapse_tuple)

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

    ##################################################################

    def get_vr_bnd2(self):
        assert self.bnds is not None, f'"bnds" attribute must be set.'

        vr_bnd2 = self.vr.header.new_record()
        vr_bnd2.id = self.bnds.get_id_bnd2()
        vcfspec_bnd2 = self.bnds.get_vcfspec_bnd2()
        varianthandler.apply_vcfspec(vr_bnd2, vcfspec_bnd2)
        vr_bnd2.info['MATEID'] = self.vr.id

        return vr_bnd2

    ##################################################################

    def set_bnds_attributes(self):
        if self.is_sv:
            vr_svinfo = breakends_module.get_vr_svinfo_standard_vr(
                    self.vr, self.fasta, self.chromdict,
                    )
            self.is_bnd1 = vr_svinfo['is_bnd1']
            self.bnds = breakends_module.get_bnds_from_vr_svinfo(
                self.vr, vr_svinfo, self.fasta, self.chromdict)
        else:
            self.is_bnd1 = None
            self.bnds = None

    def _refine_vr_InfoFormatValues(self):
        infoformat.refine_vr_InfoFormatValue(self.vr)

    ##################################################################

    def make_rpplist(
            self, bam,
            no_matesearch=False, set_alleleinfo=False,
            fetch_padding_common=readplus.FETCH_PADDING_COMMON,
            fetch_padding_sv=readplus.FETCH_PADDING_SV,
            new_fetch_padding=readplus.NEW_FETCH_PADDING,
            long_insert_threshold=readplus.LONG_INSERT_THRESHOLD,
            flanklen=alleleinfosetup.DEFAULT_FLANKLEN,
            flanklen_parside=alleleinfosetup_sv.DEFAULT_FLANKLEN_PARSIDE,
            flanklen_bndside=alleleinfosetup_sv.DEFAULT_FLANKLEN_BNDSIDE):
        if self.is_sv:
            rpplist = readplus.get_rpplist_sv(
                bam=bam, fasta=self.fasta, chromdict=self.chromdict, 
                bnds=self.bnds, view=False,
                fetch_padding_common=fetch_padding_common,
                fetch_padding_sv=fetch_padding_sv,
                new_fetch_padding=new_fetch_padding,
                long_insert_threshold=long_insert_threshold)
            if set_alleleinfo:
                rpplist.update_alleleinfo_sv(
                    bnds=self.bnds,
                    flanklen_parside=flanklen_parside,
                    flanklen_bndside=flanklen_bndside)
        else:
            rpplist = readplus.get_rpplist_nonsv(
                bam=bam, fasta=self.fasta, chromdict=self.chromdict, 
                chrom=self.vcfspec.chrom, start0=self.vcfspec.pos0, 
                end0=self.vcfspec.end0, view=False, 
                no_matesearch=no_matesearch,
                fetch_padding_common=fetch_padding_common,
                new_fetch_padding=new_fetch_padding,
                long_insert_threshold=long_insert_threshold)
            if set_alleleinfo:
                rpplist.update_alleleinfo(
                    vcfspec=self.vcfspec, flanklen=flanklen)

        return rpplist

    def set_rpplist(self, sampleid, bam,
                    no_matesearch=False, set_alleleinfo=False):
        self.rpplist_dict[sampleid] = self.make_rpplist(
            bam, no_matesearch=no_matesearch, 
            set_alleleinfo=set_alleleinfo)

    def get_rpplist(self, sampleid, bam=None):
        if sampleid not in self.rpplist_dict:
            if bam is None:
                raise Exception(
                    f'rpplist for the sampleid {sampleid} is not '
                    f'prepared. To make it, bam is required.')
            self.set_rpplist(sampleid, bam)
        return self.rpplist_dict[sampleid]

    def make_readstats_data(self, bam, no_matesearch=False, **kwargs):
        readstats_data = libreadstats.get_readstats_data(
            self.vcfspec, bam, self.fasta, self.chromdict,
            no_matesearch=no_matesearch, **kwargs)

        return readstats_data

    def set_readstats_data(self, sampleid, bam, verbose=False):
        if verbose:
            print(f'Getting readstats_data for the sample {sampleid}')
        self.readstats_datas[sampleid] = self.make_readstats_data(bam)

    def get_readstats_data(self, sampleid, bam=None):
        if sampleid not in self.readstats_datas:
            if bam is None:
                raise Exception(
                    f'readstats_data for the sampleid {sampleid} is not '
                    f'prepared. To make it, bam is required.')
            self.set_readstats_data(sampleid, bam, verbose=False)
        return self.readstats_datas[sampleid]

    def load_readstats(self):
        for sampleid in self.vr.samples.keys():
            if self.check_NA_format(sampleid, libreadstats.ReadStats.meta['ID']):
                self.readstats_dict[sampleid] = None
            else:
                readstats = libreadstats.ReadStats()
                readstats.load_vr_format(self.vr, sampleid)
                readstats.postprocess()
                self.readstats_dict[sampleid] = readstats

    ##################################################################

    def reset_sample_filter(self, sampleid):
        self.set_format(sampleid, vpfilter.FORMAT_FILTER_META['ID'], None)

    def reset_sample_filter_all(self):
        for sampleid in self.vr.header.samples:
            self.reset_sample_filter(sampleid)

    def check_sample_filter(self, sampleid):
        sample_filter = self.get_format(sampleid, 
                                        vpfilter.FORMAT_FILTER_META['ID'],
                                        collapse_tuple=False)
        return sample_filter == ('PASS',)

    def add_sample_filter(self, sampleid, value: str):
        if value == 'PASS':
            self.set_format(sampleid, vpfilter.FORMAT_FILTER_META['ID'], ('PASS',))
        else:
            if self.check_NA_format(sampleid, vpfilter.FORMAT_FILTER_META['ID']):
                new_val = (value,)
            else:
                old_val = self.get_format(sampleid, 
                                          vpfilter.FORMAT_FILTER_META['ID'], 
                                          collapse_tuple=False)
                new_val = tuple(set(old_val + (value,)))

            self.set_format(sampleid, vpfilter.FORMAT_FILTER_META['ID'], new_val)

    def show_sample_filters(self):
        for sampleid in self.vr.header.samples:
            print(sampleid, self.get_format(sampleid, 
                                            vpfilter.FORMAT_FILTER_META['ID'], 
                                            collapse_tuple=False))

    ##################################################################

    def show_igv(self, bam, igv, tmpbam_dir=None, no_matesearch=False):
        rpplist = self.make_rpplist(bam, no_matesearch=no_matesearch,
                                    set_alleleinfo=True)
        if self.is_sv:
            rpplist.set_alleleinfo_tag_rp_sv(self.bnds)
            rpplist.set_alleleinfo_tag_rpp_sv(self.bnds)
        else:
            rpplist.set_alleleinfo_tag_rp(self.vcfspec)
            rpplist.set_alleleinfo_tag_rpp(self.vcfspec)

        # set tmpbam_path
        if tmpbam_dir is None:
            tmpbam_dir = os.getcwd()
        tmpbam_path = workflow.get_tmpfile_path(suffix='.bam', 
                                                where=tmpbam_dir)

        # main
        rpplist.write_bam(outfile_path=tmpbam_path)
        igv.load([tmpbam_path])

        if self.is_sv:
            (reads_range0_bnd1, 
             reads_range0_bnd2) = rpplist.get_ref_range0_sv(self.bnds)
            pos_range0_bnd1 = self.bnds.get_pos_range0_bnd1()
            pos_range0_bnd2 = self.bnds.get_pos_range0_bnd2()

            width_bnd1 = max(
                pos_range0_bnd1.start - reads_range0_bnd1.start,
                reads_range0_bnd1.stop - pos_range0_bnd1.stop)
            width_bnd2 = max(
                pos_range0_bnd2.start - reads_range0_bnd2.start,
                reads_range0_bnd2.stop - pos_range0_bnd2.stop)
            width = max(width_bnd1, width_bnd2)

            igv.goto([self.bnds], width=width)
        else:
            reads_ref_range0 = rpplist.get_ref_range0(self.vcfspec.chrom)
            width = max(
                self.vcfspec.pos0 - reads_ref_range0.start,
                reads_ref_range0.stop - self.vcfspec.pos0)
                
            igv.goto([self.vcfspec], width=width)

        igv.cmd('group')
        igv.viewaspairs()
        igv.cmd(f'colorBy TAG {readplus.ALLELEINFO_TAG_RPP}')

        os.remove(tmpbam_path)

    def show_readcounts(self, sampleids=None, **kwargs):
        if sampleids is None:
            sampleids = sorted(self.vr.samples.keys())

        variantviz.show_readcounts(self, 
                                     sampleid_order=sampleids,
                                     title=str(self.vcfspec),
                                     **kwargs)

    @common.get_deco_num_set(('sampleid', 'bam'), 1)
    def show_readstats_data(self, sampleid=None, bam=None, varpos_key='left'):
        if sampleid is None:
            readstats_data = libreadstats.get_readstats_data(
                self.vcfspec, bam, self.fasta, self.chromdict)
        else:
            readstats_data = self.get_readstats_data(sampleid, bam=bam)

        variantviz.show_readstats_data(readstats_data, title=sampleid,
                                      varpos_key=varpos_key)

    ##################################################################

    def calc_readcounts(self, bam, no_matesearch=True):
        """Designed only for non-sv"""

        rpplist = self.make_rpplist(bam, no_matesearch=no_matesearch,
                                    set_alleleinfo=True)
        readcounts = rpplist.get_readcounts(self.vcfspec)
        
        return readcounts

    def get_ponfilter(self, sampleids=None, **kwargs):
        if sampleids is None:
            readstats_dict = self.readstats_dict
        else:
            readstats_dict = {sampleid: self.readstats_dict[sampleid]
                              for sampleid in sampleids}

        return vpfilter.PonFilter(self.vcfspec, readstats_dict, **kwargs)

    def get_total_rppcount(self, sampleid):
        """Sum of rpp counts for all alleleclasses except None"""
        readstats = self.readstats_dict[sampleid]
        return readstats.get_total_rppcount()

    def get_vaf(self, sampleid, alleleindex=1, ndigits=None):
        readstats = self.readstats_dict[sampleid]
        vaf = readstats.get_vaf(alleleindex)
        if np.isnan(vaf):
            return vaf
        else:
            if ndigits is None:
                return vaf
            else:
                return round(vaf, ndigits)

    def get_popfreq(self, popname):
        val = self.annotdb.popfreq[popname]
        if val is None:
            val = 0

        return val

    def get_cosmic_total_occur(self):
        val = self.annotdb.cosmic['total_occurrence']
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

    def set_gr(self):
        chroms = list()
        starts = list()
        ends = list()
        for vp in self:
            ref_range0 = vp.vcfspec.REF_range0
            chroms.append(vp.vr.contig)
            starts.append(ref_range0.start)
            ends.append(ref_range0.stop)
        self._gr = pyranges.from_dict({'Chromosome': chroms,
                                       'Start': starts,
                                       'End': ends})

    def get_gr(self):
        if self._gr is None:
            self.set_gr()
        return self._gr

    def get_isec(self, gr):
        overlaps = self.get_gr().count_overlaps(gr, overlap_col='count')
        marks = overlaps.df.loc[:, 'count'] > 0
        vplist_filtered = VariantPlusList()
        vplist_filtered.extend(itertools.compress(self, marks))
        return vplist_filtered

    def get_sigresult(self, catalogue_type='sbs96', **kwargs):
        signatureresult = importlib.import_module('.'.join([top_package_name, 'signature', 'signatureresult']))

        vcfspec_iter = (vp.vcfspec for vp in self)
        refver = self[0].refver
        sigresult = signatureresult.get_sigresult_from_vcfspecs(
            vcfspec_iter, refver=refver, catalogue_type=catalogue_type,
            **kwargs)
        return sigresult


def get_vp_sortkey(chromdict):
    vr_sortkey = common.get_vr_sortkey(chromdict)
    def vp_sortkey(vp):
        return vr_sortkey(vp.vr)

    return vp_sortkey

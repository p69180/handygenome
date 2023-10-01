import re
import sys
import os
import inspect
import itertools
import functools
import collections
import json
import pickle
import shutil

import pysam
import scipy
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import handygenome
import handygenome.deco as deco
import handygenome.tools as tools
import handygenome.peakutils as peakutils
import handygenome.logutils as logutils
import handygenome.refgenome.refgenome as refgenome
import handygenome.refgenome.refverfile as refverfile
from handygenome.refgenome.refverfile import PARUnavailableError

from handygenome.genomedf.genomedf import GenomeDataFrame as GDF
from handygenome.cnv.depth import DepthRawDataFrame as DepthRawDF
from handygenome.cnv.depth import DepthSegmentDataFrame as DepthSegDF
from handygenome.cnv.baf import BAFRawDataFrame as BAFRawDF
from handygenome.cnv.baf import BAFSegmentDataFrame as BAFSegDF
#from handygenome.cnv.segment import SegmentDataFrame as SegDF
from handygenome.cnv.cnvsegment import CNVSegmentDataFrame as CNVSegDF

import handygenome.cnv.cncall as libcncall
import handygenome.cnv.mosdepth as libmosdepth
import handygenome.cnv.baf as libbaf
import handygenome.cnv.bafcorrection as bafcorrection

import handygenome.plot.genomeplot as libgenomeplot
from handygenome.plot.genomeplot import GenomePlotter




#############
# CNVSample #
#############

class CNVSample:
    default_window = {'panel': 100, 'wgs': 1000}
    simple_attr_keys = (
        'bam_path',
        'vcf_path',
        'sampleid',
        'refver',
        'is_female',
        'mode',
        'nproc',
        'verbose',
        'ploidy',
        'genomeplotter_kwargs',
    )
    save_paths_basenames = {
        'simple_attrs': 'simple_attrs.json',

        # target region
        'raw_target_region': 'raw_target_region.tsv.gz',

        # raw data
        'depth_rawdata': 'depth_rawdata.tsv.gz',
        'baf_rawdata': 'baf_rawdata.tsv.gz',

        # segment
        'depth_segment': 'depth_segment.tsv.gz',
        'baf_noedit_segment': 'baf_noedit_segment.tsv.gz',
        'baf_edited_segment': 'baf_edited_segment.tsv.gz',
        'baf_rmzero_noedit_segment': 'baf_rmzero_noedit_segment.tsv.gz',
        'merged_segment': 'merged_segment.tsv.gz',
    }

    def __repr__(self):
        return f'{self.__class__.__name__} object (sampleid:{self.simple_attrs["sampleid"]})'

    def __getattr__(self, key):
        #if key in self.__class__.simple_attr_keys:
        if key in self.simple_attrs.keys():
            return self.simple_attrs[key]
        else:
            super().__getattr__(key)

    ##########
    # logger #
    ##########

    def log(self, msg, verbose=None):
        if verbose is None:
            verbose = self.simple_attrs["verbose"]
        if verbose:
            logutils.log(f'CNVSample {self.simple_attrs["sampleid"]} - {msg}')

    #############
    # decorator #
    #############

    @staticmethod
    def deco_sanitycheck_sampletype(func):
        sig = inspect.signature(func)
        assert 'sampletype' in sig.parameters.keys()

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()
            if ba.arguments['sampletype'] not in ('normal', 'tumor'):
                raise Exception(f'"sampletype" argument must be either "normal" or "tumor".')

            return func(*args, **kwargs)

        return wrapper

    @staticmethod
    def get_deco_logging(msg):
        def decorator(func):
            sig = inspect.signature(func)

            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                ba = sig.bind(*args, **kwargs)
                ba.apply_defaults()

                ba.arguments['self'].log('Beginning ' + msg)
                result = func(*args, **kwargs)
                ba.arguments['self'].log('Finished ' + msg)

                return result

            return wrapper

        return decorator

    @staticmethod
    def deco_nproc(func):
        sig = inspect.signature(func)
        if 'nproc' not in sig.parameters.keys():
            raise Exception(
                f'Decorated plotter method does not have '
                f'"nproc" parameter.'
            )

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()
            if ba.arguments['nproc'] == 0:
                ba.arguments['nproc'] = ba.arguments['self'].nproc

            return func(*ba.args, **ba.kwargs)

        return wrapper

    @staticmethod
    def deco_plotter(func):
        sig = inspect.signature(func)
        if 'genomeplotter_kwargs' not in sig.parameters.keys():
            raise Exception(
                f'Decorated plotter method does not have '
                f'"genomeplotter_kwargs" parameter.'
            )

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()

            if ba.arguments['genomeplotter_kwargs'] is not None:
                cnvsample_obj = ba.arguments['self']
                gplotter_kwargs = ba.arguments['genomeplotter_kwargs'].copy()

                gplotter_kwargs['refver'] = cnvsample_obj.simple_attrs["refver"]
                new_gplotter = GenomePlotter(**gplotter_kwargs)

                old_gplotter = cnvsample_obj.reset_genomeplotter(new_gplotter)
                result = func(*args, **kwargs)
                cnvsample_obj.reset_genomeplotter(old_gplotter)
            else:
                result = func(*args, **kwargs)

            return result

        return wrapper

    ###############
    # initializer #
    ###############

    #def __del__(self):
    #    del self.raw_target_region
    #    for datatype, sampletype in itertools.product(
    #        ['depth', 'baf'], ['normal', 'tumor'],
    #    ):
    #        del self.rawdata[datatype][sampletype]

    @classmethod
    @deco.get_deco_num_set_differently(
        ('target_region_path', 'target_region_gdf'), 2, 'lt'
    )
    def from_files(
        cls, 

        *,

        # required simple attrs
        sampleid, 
        refver,
        is_female,
        mode,

        # depth
        bam_path,

        # baf
        vcf_path=None,
        vcf_sampleid=None,
        bafdf=None,

        # optional simple attrs
        target_region_path=None,
        target_region_gdf=None,

        ploidy=2,
        nproc=1,
        genomeplotter_kwargs=dict(),
        verbose=True,

        # others
        norm_method='plain',
    ):
        cls.init_sanitycheck(
            mode, 
            target_region_path,
            target_region_gdf,
            vcf_path,
            vcf_sampleid,
            bafdf,
        )

        result = cls()

        # set simple attributes
        result.simple_attrs = {
            'bam_path': bam_path,
            'vcf_path': vcf_path,
            'vcf_sampleid': vcf_sampleid,
            'sampleid': sampleid,
            'refver': refver,
            'is_female': is_female,
            'mode': mode,
            'nproc': nproc,
            'verbose': verbose,
            'ploidy': ploidy,
            'genomeplotter_kwargs': genomeplotter_kwargs,
        }

        # initiate gdf data dict
        result.init_gdfs()

        # target region
        result.load_raw_target_region(target_region_path, target_region_gdf)

        # raw data
        result.load_depth(bam_path)
        #result.make_normalized_depth()
        result.load_baf(vcf_path, vcf_sampleid, bafdf)

        # genomeplotter
        result.init_genomeplotter(**result.simple_attrs['genomeplotter_kwargs'])

        # plotdata
        result.plotdata_cache = dict()

        return result

    def init_gdfs(self):
        self.gdfs = {
            'raw_target_region': None,
            'target_region': None,
            'excluded_region': None,
            'depth_rawdata': None,
            'baf_rawdata': None,
            'depth_segment': None,
            'baf_noedit_segment': {baf_idx: None for baf_idx in self.get_baf_indexes()},
            'baf_edited_segment': {baf_idx: None for baf_idx in self.get_baf_indexes()},
            'baf_rmzero_noedit_segment': {baf_idx: None for baf_idx in self.get_baf_indexes()},
            'merged_segment': None,
        }


    @staticmethod
    def init_sanitycheck(
        mode, 
        target_region_path,
        target_region_gdf,
        vcf_path,
        vcf_sampleid,
        bafdf,
    ):
        assert mode in ('wgs', 'panel')
        target_region_given = (
            (target_region_path is not None)
            or (target_region_gdf is not None)
        )
        if (mode == 'panel') and (not target_region_given):
            raise Exception(f'target_region must be given when "mode" is "panel"')
        elif (mode == 'wgs') and target_region_given:
            raise Exception(f'target_region must not be given when "mode" is "wgs"')

        # baf loading args
        if not (
            ((vcf_path is not None) and (vcf_sampleid is not None))
            ^ (bafdf is not None)
        ):
            raise Exception(f'BAF loading argument usage: 1) "vcf_path" and "vcf_sampleid" ; 2) "bafdf"')

    def load_raw_target_region(self, target_region_path, target_region_gdf):
        # target_region
        if target_region_path is not None:
            raw_target_region = GDF.read_tsv(target_region_path).drop_annots()
        elif target_region_gdf is not None:
            assert isinstance(target_region_gdf, GDF)
            raw_target_region = GDF.from_frame(
                target_region_gdf.drop_annots().df,
                refver=self.refver,
            )
        else:
            raw_target_region = None

        self.gdfs['raw_target_region'] = raw_target_region

    @deco_nproc
    def load_depth(self, bam_path, window=None, t=1, nproc=0):
        self.log(f'Beginning making depth rawdata from bam file')

        # prepare arguments
        if self.simple_attrs["mode"] == 'wgs':
            region_gdf = GDF.all_regions(self.refver, assembled_only=True)
        elif self.simple_attrs["mode"] == 'panel':
            region_gdf = self.gdfs['raw_target_region']

        if window is None:
            window = self.__class__.default_window[self.simple_attrs["mode"]]

        # main
        depth_gdf = libmosdepth.run_mosdepth(
            bam_path, 
            t=t, 
            nproc=nproc,
            use_median=False, 

            region_gdf=region_gdf,

            window=window,
            verbose=self.simple_attrs['verbose'],
        )
        self.gdfs['depth_rawdata'] = depth_gdf

        self.log(f'Finished making depth rawdata from bam file')

    @deco_nproc
    def load_baf(self, vcf_path, vcf_sampleid, bafdf, nproc=0):
        if bafdf is None:
            self.log(f'Beginning making BAF rawdata from vcf file')

            bafdfs_bysample = libbaf.get_bafdfs_from_vcf(
                vcf_path=vcf_path,
                sampleids=vcf_sampleid,
                n_allele=self.simple_attrs["ploidy"],
                nproc=nproc,
                exclude_other=False,
            )
            bafdf = bafdfs_bysample[vcf_sampleid]

            self.log(f'Finished making BAF rawdata from vcf file')

        bafdf = bafdf.remove_overlapping_rows()
        if self.mode == 'panel':
            bafdf = bafdf.intersect(self.get_raw_target_region())

        self.gdfs['baf_rawdata'] = bafdf

#        # set target_region
#        self.set_target_region(sampleid, mode, target_region)
#
#        # normal mean ploidy
#        self.set_normal_mean_ploidy(sampleid)
#
#        # load germline vcf
#        if germline_vcf_path is not None:
#            LOGGER_INFO.info('Loading tumor germline vcf')
#            self.load_germline_vcf(
#                sampleid=sampleid, 
#                vcf_path=germline_vcf_path, 
#                vcf_sampleid_tumor=vcf_sampleid_tumor,
#                vcf_sampleid_normal=vcf_sampleid_normal,
#                logging_lineno=50000,
#                nproc=vcfload_nproc,
#            )
#        elif germline_vafdf_path is not None:
#            LOGGER_INFO.info('Loading tumor germline vaf dataframe')
#            self.load_germline_vafdf(sampleid, germline_vafdf_path)
#
#        #self.postprocess_bafdf(sampleid)
#
#        # postprocess depths
#        self.postprocess_depth(sampleid, verbose=verbose, norm_method=norm_method)

    #################
    # genomeplotter #
    #################

    def init_genomeplotter(self, **kwargs):
        self.genomeplotter = GenomePlotter(refver=self.simple_attrs["refver"], **kwargs)

    def reset_genomeplotter(self, new_genomeplotter=None):
        old_genomeplotter = self.genomeplotter

        if new_genomeplotter is None:
            new_genomeplotter = GenomePlotter(self.simple_attrs["refver"])
        self.genomeplotter = new_genomeplotter

        return old_genomeplotter

    ##############
    # properties #
    ##############

    ##############
    # CN calling #
    ##############

    def add_clonal_solution_germline(self):
        """Only for germline sample"""
        cnv_seg_gdf = self.get_merged_segment()
        clonal_CNt, clonal_Bt = libcncall.clonal_solution_from_germline(
            cnv_seg_gdf, self.is_female, self.ploidy,
        )
        cnv_seg_gdf.assign_clonal_CN(clonal_CNt)
        cnv_seg_gdf.assign_clonal_B(clonal_Bt)

    #########################
    # target region setting #
    #########################

    def set_excluded_region(self):
        excluded_region_list = list()

        # baf
        for baf_seg_gdf in self.get_baf_edited_segment_dict().values():
            bafseg_selector = baf_seg_gdf.get_filter()
            baf_excl_region_gdf = baf_seg_gdf.drop_annots().loc[~bafseg_selector, :]
            excluded_region_list.append(baf_excl_region_gdf)

        # depth
        depth_seg_gdf = self.get_depth_segment() 
        depthseg_selector = depth_seg_gdf.get_filter()
        depth_excl_region_gdf = depth_seg_gdf.drop_annots().loc[~depthseg_selector, :]
        excluded_region_list.append(depth_excl_region_gdf)

        # N region
        Nregion_gdf = refverfile.get_Nregion_gdf(self.refver)
        excluded_region_list.append(Nregion_gdf)
        
        # male Y PAR
        if not self.is_female:
            try:
                par_gdf = refverfile.get_par_gdf(self.refver)
            except PARUnavailableError:
                pass
            else:
                y_chrom = self.chromdict.XY_names[1]
                y_par_gdf = par_gdf.subset_chroms(y_chrom)
                excluded_region_list.append(y_par_gdf)

        # get union of all
        excluded_region = functools.reduce(
            lambda x, y: x.set_union(y), 
            excluded_region_list,
        )
        self.gdfs['excluded_region'] = excluded_region

    def set_target_region(self):
        """Must be run after depth and baf segment processing"""

        # make allregion
        if self.mode == 'wgs':
            allregion = GDF.all_regions(refver=self.refver, assembled_only=True)
        elif mode == 'panel':
            allregion = self.get_raw_target_region()

        if (
            self.is_female
            and self.chromdict.check_has_XY()
        ):
            y_chrom = self.chromdict.XY_names[1]
            selected_chroms = set(allregion.chroms).difference([y_chrom])
            allregion = allregion.subset_chroms(selected_chroms)

        # subtract excluded region
        self.set_excluded_region()
        target_region = allregion.subtract(self.get_excluded_region())
        target_region.sort()

        self.gdfs['target_region'] = target_region

    def find_onecopy_depth(self):
        return self.get_depth_segment().find_onecopy_depth(
            self.is_female, self.ploidy,
        )

    def find_germline_intcopy_validbaf_region(self, factors=(0.3, 2)):
        assert self.simple_attrs["mode"] == 'wgs'
        assert self.get_depth_segment('normal') is not None
        assert all(
            (self.get_baf_segment('normal', baf_idx) is not None)
            for baf_idx in self.get_baf_indexes()
        )

        depth_seg_gdf = self.get_depth_segment('normal').copy()
        baf_seg_gdfs = self.get_baf_segment_dict('normal')

        # onecopy depth
        onecopy_depth = self.get_onecopy_depth(depth_seg_gdf)
        depth_seg_gdf['onecopy_depth_ratio'] = (
            depth_seg_gdf[DepthRawDF.norm_depth_colname]
            / onecopy_depth
        )

        print(depth_seg_gdf)


#        global_mean = np.average(
#            upscaled_depth_df['mean_depth'], 
#            weights=(upscaled_depth_df['End'] - upscaled_depth_df['Start']),
#        )
#        selector = segdf['depth_segment_mean'].between(
#            global_mean * factors[0], 
#            global_mean * factors[1], 
#        )
#        included_segments = segdf.loc[selector, :]
#
#        target_region_gr = pr.PyRanges(depth_df).drop().overlap(
#            pr.PyRanges(included_segments)
#        ).merge()
#
#        return target_region_gr


    ####################
    # pickle save/load #
    ####################

    def save_pickle(self, filepath):
        if not filepath.endswith('.pickle'):
            raise Exception(f'File name must end with ".pickle"')
        
        if os.path.exists(filepath):
            logutils.log(f'Overwriting an existing file', level='warning')

        with open(filepath, 'wb') as outfile:
            pickle.dump(self, outfile)

    @classmethod
    def load_pickle(cls, filepath):
        with open(filepath, 'rb') as infile:
            return pickle.load(infile)

    ######################
    # textfile save/load #
    ######################

    def save_simple_attrs(self, topdir):
        savepath = os.path.join(
            topdir, 
            self.__class__.save_paths_basenames['simple_attrs'],
        )
        with open(savepath, 'wt') as outfile:
            json.dump(self.simple_attrs, outfile)

    def load_simple_attrs(self, topdir):
        savepath = os.path.join(
            topdir, 
            self.__class__.save_paths_basenames['simple_attrs'],
        )
        with open(savepath, 'rt') as infile:
            self.simple_attrs = json.load(infile)

    ###
    
    def save_nonbafidx_gdf(self, gdf_key, topdir):
        gdf = self.gdfs[gdf_key]
        if gdf is not None:
            savepath = os.path.join(
                topdir, 
                self.__class__.save_paths_basenames[gdf_key],
            )
            gdf.write_tsv(savepath)

    def load_nonbafidx_gdf(self, gdf_key, topdir, GDF_class):
        savepath = os.path.join(
            topdir, 
            self.__class__.save_paths_basenames[gdf_key],
        )
        if os.path.exists(savepath):
            self.gdfs[gdf_key] = GDF_class.read_tsv(savepath, self.simple_attrs["refver"])

    ###

    def save_bafidx_gdf(self, gdf_key, topdir):
        for baf_idx in self.get_baf_indexes():
            gdf = self.gdfs[gdf_key][baf_idx]
            if gdf is not None:
                basename = re.sub(
                    '\.tsv\.gz$', 
                    f'_{baf_idx}.tsv.gz', 
                    self.__class__.save_paths_basenames[gdf_key],
                )
                savepath = os.path.join(topdir, basename)
                gdf.write_tsv(savepath)

    def load_bafidx_gdf(self, gdf_key, topdir, GDF_class):
        for baf_idx in self.get_baf_indexes():
            basename = re.sub(
                '\.tsv\.gz$', 
                f'_{baf_idx}.tsv.gz', 
                self.__class__.save_paths_basenames[gdf_key],
            )
            savepath = os.path.join(topdir, basename)
            if os.path.exists(savepath):
                self.gdfs[gdf_key][baf_idx] = GDF_class.read_tsv(savepath, self.simple_attrs["refver"])

    def save_tsv(self, topdir, force=False):
        if os.path.exists(topdir):
            if force:
                shutil.rmtree(topdir)
            else:
                raise Exception(f'Output directory must not exist in advance')
        os.mkdir(topdir)

        self.log(f'Beginning saving tsv files')
        
        # simple attrs
        self.save_simple_attrs(topdir)

        # target region
        self.save_nonbafidx_gdf('raw_target_region', topdir)

        # raw data
        self.save_nonbafidx_gdf('depth_rawdata', topdir)
        self.save_nonbafidx_gdf('baf_rawdata', topdir)

        # segment
        self.save_nonbafidx_gdf('depth_segment', topdir)
        self.save_bafidx_gdf('baf_noedit_segment', topdir)
        self.save_bafidx_gdf('baf_edited_segment', topdir)
        self.save_bafidx_gdf('baf_rmzero_noedit_segment', topdir)
        self.save_nonbafidx_gdf('merged_segment', topdir)

        self.log(f'Finished saving tsv files')

    @classmethod
    def load_tsv(cls, topdir):
        result = cls()

        # simple attrs
        result.load_simple_attrs(topdir)

        # initiate gdf data dict
        result.init_gdfs()

        # target region
        result.load_simple_attrs(topdir)

        # raw data
        result.load_nonbafidx_gdf('depth_rawdata', topdir, DepthRawDF)
        result.load_nonbafidx_gdf('baf_rawdata', topdir, BAFRawDF)

        # segment
        result.load_nonbafidx_gdf('depth_segment', topdir, DepthSegDF)
        result.load_bafidx_gdf('baf_noedit_segment', topdir, BAFSegDF)
        result.load_bafidx_gdf('baf_edited_segment', topdir, BAFSegDF)
        result.load_bafidx_gdf('baf_rmzero_noedit_segment', topdir, BAFSegDF)
        result.load_nonbafidx_gdf('merged_segment', topdir, CNVSegDF)

        # genomeplotter
        result.init_genomeplotter(**result.simple_attrs['genomeplotter_kwargs'])

        # plotdata
        result.plotdata_cache = dict()

        return result
            
    ##############
    # properties #
    ##############

    @property
    def chromdict(self):
        return refgenome.get_chromdict(self.simple_attrs["refver"])

    #################
    # data fetchers #
    #################

    def get_raw_target_region(self):
        return self.gdfs['raw_target_region']

    def get_target_region(self):
        return self.gdfs['target_region']

    def get_excluded_region(self):
        return self.gdfs['excluded_region']

    def get_depth_rawdata(self):
        return self.gdfs['depth_rawdata']

    def get_baf_rawdata(self, baf_idx, rmzero=False):
        baf_raw_gdf = self.gdfs['baf_rawdata']
        if baf_raw_gdf is None:
            return None
        else:
            baf_raw_gdf = baf_raw_gdf.choose_annots(baf_idx)
            if rmzero:
                selector = baf_raw_gdf[baf_idx] > 0
                baf_raw_gdf = baf_raw_gdf.loc[selector, :]

            return baf_raw_gdf

    def get_depth_segment(self):
        return self.gdfs['depth_segment']

    def get_baf_noedit_segment(self, baf_idx, rmzero=False):
        segment_dict_key = ('baf_rmzero_noedit_segment' if rmzero else 'baf_noedit_segment')
        result = self.gdfs[segment_dict_key][baf_idx]
        return result

    def get_baf_noedit_segment_dict(self, rmzero=False):
        return {
            baf_idx: self.get_baf_noedit_segment(baf_idx, rmzero=rmzero)
            for baf_idx in self.get_baf_indexes()
        }

    def get_baf_edited_segment(self, baf_idx, rmzero=False):
        segment_dict_key = 'baf_edited_segment'
        result = self.gdfs[segment_dict_key][baf_idx]
        return result

    def get_baf_edited_segment_dict(self, rmzero=False):
        return {
            baf_idx: self.get_baf_edited_segment(baf_idx, rmzero=rmzero)
            for baf_idx in self.get_baf_indexes()
        }

    def get_merged_segment(self):
        return self.gdfs['merged_segment']

    ##############################
    # depth rawdata modification #
    ##############################

    @get_deco_logging(f'depth segmentation')
    @deco_nproc
    def make_depth_segment(self, rawdepth=False, nproc=0, **kwargs):
        kwargs['nproc'] = nproc
        depth_raw_gdf = self.get_depth_rawdata()
        annot_colname = (
            DepthRawDF.depth_colname
            if rawdepth else
            DepthRawDF.norm_depth_colname
        )
        self.gdfs['depth_segment'] = depth_raw_gdf.get_segment(
            annot_colname=annot_colname,
            **kwargs,
        )

    @get_deco_logging(f'depth segment modification')
    @deco_nproc
    def edit_depth_segment(self, nproc=0, readlength=151, MQ_cutoff=(40, None)):
        self.make_normalized_depth()
        self.add_MQ_to_depth_rawdata(readlength=readlength, nproc=nproc)
        self.add_rawdata_to_depth_segment(rawdepth=False, nproc=nproc)
        self.add_filter_to_depth_segment(MQ_cutoff=MQ_cutoff)

    ###

    def make_normalized_depth(self):
        self.get_depth_rawdata().set_normalized_depth()

    @get_deco_logging(f'MQ annotation to depth rawdata')
    @deco_nproc
    def add_MQ_to_depth_rawdata(self, readlength=151, nproc=0):
        depth_rawdata = self.get_depth_rawdata()
        if not depth_rawdata.check_has_MQ():
            depth_rawdata.add_MQ(
                bam_path=self.bam_path,
                readlength=readlength,
                nproc=nproc,
                verbose=False,
            )

    @get_deco_logging(f'rawdata annotation to depth segment')
    @deco_nproc
    def add_rawdata_to_depth_segment(
        self, rawdepth=False, nproc=0,
    ):
        depth_seg_gdf = self.get_depth_segment()
        depth_raw_gdf = self.get_depth_rawdata()
        depth_seg_gdf.add_rawdata_info(
            depth_raw_gdf, 
            merge_methods=['mean', 'std'],
            rawdepth=rawdepth,
            nproc=nproc,
        )

    def add_filter_to_depth_segment(self, MQ_cutoff=(40, None)):
        self.get_depth_segment().add_filter(MQ_cutoff=MQ_cutoff)

    ############################
    # baf rawdata modification #
    ############################

    @get_deco_logging(f'BAF segmentation')
    @deco_nproc
    def make_baf_noedit_segment(self, nproc=0, rmzero=False, **kwargs):
        kwargs['nproc'] = nproc
        for baf_idx in self.get_baf_indexes():
            self.log(f'Beginning {baf_idx}, rmzero={rmzero}')

            baf_raw_gdf = self.get_baf_rawdata(baf_idx, rmzero=rmzero)
            baf_seg_gdf = baf_raw_gdf.get_segment(baf_idx, **kwargs)
            #baf_seg_gdf = baf_seg_gdf.fill_gaps(edit_first_last=(self.mode == 'wgs'))

            segment_dict_key = ('baf_rmzero_noedit_segment' if rmzero else 'baf_noedit_segment')
            self.gdfs[segment_dict_key][baf_idx] = baf_seg_gdf

            self.log(f'Finished {baf_idx}, rmzero={rmzero}')

    @deco_nproc
    def add_rawdata_to_bafseg_base(
        self, 
        baf_seg_gdf, 
        baf_raw_gdf, 
        baf_idx,
        distinfo_kwargs,
        nproc=0,
    ):
        baf_seg_gdf.add_rawdata_info(
            baf_raw_gdf, 
            baf_idx,
            distinfo_keys=['ndata', 'center', 'width', 'left_ips', 'right_ips'],
            distinfo_kwargs=distinfo_kwargs,
            merge_methods=['mean', 'std'],
            nproc=nproc,
        )
        #return annotated_baf_seg_gdf

    @get_deco_logging(f'rawdata annotation to baf segment')
    @deco_nproc
    def add_rawdata_to_baf_noedit_segment(self, nproc=0):
        for baf_idx in self.get_baf_indexes():
            self.add_rawdata_to_bafseg_base(
                baf_seg_gdf=self.get_baf_noedit_segment(baf_idx), 
                baf_raw_gdf=self.get_baf_rawdata(baf_idx), 
                baf_idx=baf_idx,
                distinfo_kwargs=dict(),
                nproc=nproc,
            )

    @deco_nproc
    @get_deco_logging(f'modifying BAF segment')
    def make_baf_edited_segment(
        self, 
        rmzero=False, 
        distinfo_kwargs=dict(),
        nproc=0,

        # low ndata merging
        ndata_cutoff=30, 
        merge_low_ndata=True,

        # corrected baf cutoff
        corrected_baf_cutoff=0.41,

        # filtering
        width_cutoff=None,
        #center_cutoff=(0.4, np.inf),
        center_cutoff=None,
    ):
        for baf_idx in self.get_baf_indexes():
            logutils.log(f'Beginning job for baf index {baf_idx}')

            baf_raw_gdf = self.get_baf_rawdata(baf_idx, rmzero=rmzero)

            # 1. annotate with bafdistinfo
            baf_edited_seg_gdf = self.get_baf_noedit_segment(baf_idx, rmzero=rmzero).copy()
            self.add_rawdata_to_bafseg_base(
                baf_seg_gdf=baf_edited_seg_gdf, 
                baf_raw_gdf=baf_raw_gdf, 
                baf_idx=baf_idx,
                distinfo_kwargs=distinfo_kwargs,
                nproc=nproc,
            )

            # 2. fill gap
            baf_edited_seg_gdf = baf_edited_seg_gdf.fill_gaps(edit_first_last=(self.mode == 'wgs'))

            # 3. merge low ndata segments
            if merge_low_ndata:
                logutils.log(f'Beginning merging low ndata segments')

                cycle = 0
                while True:
                    cycle += 1
                    logutils.log(f'Cycle {cycle}')

                    # merge between low segments
                    logutils.log(f'Beginning merging low ndata segments')
                    baf_edited_seg_gdf = baf_edited_seg_gdf.merge_low_ndata_segments(
                        cutoff=ndata_cutoff,
                        nproc=nproc,
                    )
                    self.add_rawdata_to_bafseg_base(
                        baf_seg_gdf=baf_edited_seg_gdf, 
                        baf_raw_gdf=baf_raw_gdf, 
                        baf_idx=baf_idx,
                        distinfo_kwargs=distinfo_kwargs,
                        nproc=nproc,
                    )
                    logutils.log(f'Finished merginig low ndata segments')

                    # incorporate low segments into adjacent high segments
                    logutils.log(f'Beginning incorporating into adjacent high ndata segments')
                    baf_edited_seg_gdf = baf_edited_seg_gdf.incoporate_low_ndata_segments(
                        cutoff=ndata_cutoff,
                        nproc=nproc,
                    )
                    self.add_rawdata_to_bafseg_base(
                        baf_seg_gdf=baf_edited_seg_gdf, 
                        baf_raw_gdf=baf_raw_gdf, 
                        baf_idx=baf_idx,
                        distinfo_kwargs=distinfo_kwargs,
                        nproc=nproc,
                    )
                    logutils.log(f'Finished incorporating into adjacent high ndata segments')

                    if not np.isnan(baf_edited_seg_gdf.get_dist_ndata()).any():
                        break

                logutils.log(f'Finished merging low ndata segments')

#            # 3. merge low ndata segments with themselves
#            if merge_low_ndata:
#                logutils.log(f'Beginning merging low ndata segments')
#                baf_edited_seg_gdf = baf_edited_seg_gdf.merge_low_ndata_segments(
#                    cutoff=ndata_cutoff,
#                    nproc=nproc,
#                )
#                self.add_rawdata_to_bafseg_base(
#                    baf_seg_gdf=baf_edited_seg_gdf, 
#                    baf_raw_gdf=baf_raw_gdf, 
#                    baf_idx=baf_idx,
#                    distinfo_kwargs=distinfo_kwargs,
#                    nproc=nproc,
#                )
#                logutils.log(f'Finished merginig low ndata segments')
#
#            # 4. merge remaining low ndata segments into adjacent high ndata segments
#            if merge_low_ndata:
#                logutils.log(f'Beginning incorporating into adjacent high ndata segments')
#                baf_edited_seg_gdf = baf_edited_seg_gdf.incoporate_low_ndata_segments(
#                    cutoff=ndata_cutoff,
#                    nproc=nproc,
#                )
#                self.add_rawdata_to_bafseg_base(
#                    baf_seg_gdf=baf_edited_seg_gdf, 
#                    baf_raw_gdf=baf_raw_gdf, 
#                    baf_idx=baf_idx,
#                    distinfo_kwargs=distinfo_kwargs,
#                    nproc=nproc,
#                )
#                logutils.log(f'Finished incorporating into adjacent high ndata segments')

            # 4. add corrected baf column
            baf_edited_seg_gdf.add_corrected_baf(cutoff=corrected_baf_cutoff)

            # 5. add filter column
            baf_edited_seg_gdf.add_filter(width_cutoff=width_cutoff, center_cutoff=center_cutoff)

            # result
            self.gdfs['baf_edited_segment'][baf_idx] = baf_edited_seg_gdf

    #########################
    # merged segment making #
    #########################

    @deco_nproc
    @get_deco_logging(f'merging depth and baf segments')
    def make_merged_segment(self, nproc=0):
        depth_segdf = self.get_depth_segment()
        baf_edited_segdf_dict = self.get_baf_edited_segment_dict()
        sub_seg_gdfs = [depth_segdf] + list(baf_edited_segdf_dict.values())

        merged_seg_gdf = CNVSegDF.from_segments(sub_seg_gdfs, nproc=nproc)
        merged_seg_gdf = merged_seg_gdf.intersect(self.get_target_region())

        cols_to_drop = (
            [depth_segdf.get_filter_colname()]
            + [x.get_filter_colname() for x in baf_edited_segdf_dict.values()]
        )
        merged_seg_gdf.drop_annots(cols_to_drop, inplace=True)

        self.gdfs['merged_segment'] = merged_seg_gdf

    ############
    # plotdata #
    ############

    def make_upscaled_depth(self, binsize):
        pass

    def make_plotdata_base(self, data, cache_key=None):
        #self.log(f'Beginning plotdata generation for {cache_key}')
        plotdata = self.genomeplotter.cconv.prepare_plot_data(data)
        #self.plotdata_cache[cache_key] = plotdata
        #self.log(f'Finished plotdata generation for {cache_key}')

        return plotdata

    @get_deco_logging(f'making BAF rawdata plotdata')
    def make_baf_rawdata_plotdata(self, rmzero=False):
        plotdata_dict = dict()
        for baf_idx in self.get_baf_indexes():
            data = self.get_baf_rawdata(baf_idx, rmzero=rmzero)
            #plotdata_key = f'{sampletype}_baf_raw_{baf_idx}'
            plotdata_dict[baf_idx] = self.make_plotdata_base(data)

        return plotdata_dict

    @get_deco_logging(f'making BAF segment plotdata')
    def make_baf_noedit_segment_plotdata(self, rmzero=False):
        segment_dict_key = ('baf_rmzero_noedit_segment' if rmzero else 'baf_noedit_segment')
        data_dict = self.gdfs[segment_dict_key]
        plotdata_dict = dict()
        for baf_idx, data in data_dict.items():
            #plotdata_key = f'{sampletype}_baf_segment_{baf_idx}'
            plotdata = self.make_plotdata_base(data)
            plotdata_dict[baf_idx] = plotdata

        return plotdata_dict

    @get_deco_logging(f'making BAF segment plotdata')
    def make_baf_edited_segment_plotdata(self, rmzero=False):
        data_dict = self.gdfs['baf_edited_segment']
        plotdata_dict = dict()
        for baf_idx, data in data_dict.items():
            #plotdata_key = f'{sampletype}_baf_segment_{baf_idx}'
            plotdata = self.make_plotdata_base(data)
            plotdata_dict[baf_idx] = plotdata

        return plotdata_dict

    @get_deco_logging(f'making depth rawdata plotdata')
    def make_depth_rawdata_plotdata(self):
        data = self.get_depth_rawdata()
        #plotdata_key = f'{sampletype}_depth_raw'
        return self.make_plotdata_base(data)

    @get_deco_logging(f'making depth segment plotdata')
    def make_depth_segment_plotdata(self):
        data = self.get_depth_segment()
        #plotdata_key = f'{sampletype}_depth_segment'
        return self.make_plotdata_base(data)

    @get_deco_logging(f'making merged segment plotdata')
    def make_merged_segment_plotdata(self):
        data = self.get_merged_segment()
        #plotdata_key = f'{sampletype}_depth_segment'
        return self.make_plotdata_base(data)

    #########################
    # HIGH level ax drawers #
    #########################

    @deco_plotter
    def draw_depth_and_baf(
        self,
        axd=None,

        # figure title - only works when figure object is not already created
        title=None,
        suptitle_kwargs=dict(),
        subplots_kwargs=dict(),

        # kwargs for low level drawers
        draw_baf_ax_kwargs=dict(),
        draw_depth_ax_kwargs=dict(),
        draw_common_kwargs=dict(n_xlabel=20),

        # drawing region configuration
        genomeplotter_kwargs=None,
    ):
        subplots_kwargs = (
            {
                'figsize': (30, 12),
                'gridspec_kw': {'hspace': 0.6},
            }
            | subplots_kwargs
        )

        if axd is None:
            mosaic = (
                [[baf_idx] for baf_idx in self.get_baf_indexes()]
                + [['depth']]
            )
            fig, axd = plt.subplot_mosaic(mosaic, **subplots_kwargs)
            if title is None:
                title = self.get_default_title()
            fig.suptitle(title, **suptitle_kwargs)
        else:
            fig = None

        self.draw_depth_ax(axd['depth'], draw_common_kwargs=draw_common_kwargs, **draw_depth_ax_kwargs)
        baf_axd_keys = set(axd.keys()).intersection(self.get_baf_indexes())
        baf_axd = {k: axd[k] for k in baf_axd_keys}
        self.draw_baf_ax(baf_axd, draw_common_kwargs=draw_common_kwargs, **draw_baf_ax_kwargs)

        return fig, axd

    ########################
    # LOW level ax drawers #
    ########################

    @deco_plotter
    def draw_baf_ax(
        self, 
        ax_dict=None, 

        # whether to remove zero-baf points
        rmzero=False,

        # fig generation params
        title=None,
        suptitle_kwargs=dict(),
        subplots_kwargs=dict(
            figsize=(30, 5),
        ),

        # drawing kwargs
        dot_kwargs=dict(),
        line_segmean_kwargs=dict(),
        line_corr_segmean_kwargs=dict(),
        line_predict_kwargs=dict(),
        line_predict_clonal_kwargs=dict(),

        # segment drawing
        draw_segment=True,
        edited_segment=True,
        segment_value_type='center',

        # axes setting
        ylabels=dict(),
        ylabel_kwargs=dict(),
        ymin=-0.1,
        ymax=0.6,
        draw_common_kwargs=dict(n_xlabel=20),

        # drawing region configuration
        genomeplotter_kwargs=None,
    ):
        """Args:
            ax_dict: dict (keys: baf0, baf1, ... ; values: ax)
            ylabels: dict (keys: baf0, baf1, ... ; values: label text)
        """
        if ax_dict is None:
            fig, axs = plt.subplots(self.simple_attrs["ploidy"] - 1, 1, **subplots_kwargs)

            ax_dict = dict()
            for baf_idx, ax in zip(
                self.get_baf_indexes(),
                np.atleast_1d(axs),
            ):
                ax_dict[baf_idx] = ax

            if title is None:
                title = self.get_default_title()
            fig.suptitle(title, **suptitle_kwargs)
        else:
            fig = None

        # prepare plotdata
        rawdata_plotdata_dict = self.make_baf_rawdata_plotdata(rmzero=rmzero)
        if draw_segment:
            if edited_segment:
                segment_dict = self.get_baf_edited_segment_dict()
                segment_plotdata_dict = self.make_baf_edited_segment_plotdata(rmzero=rmzero)
            else:
                segment_dict = self.get_baf_noedit_segment_dict()
                segment_plotdata_dict = self.make_baf_noedit_segment_plotdata(rmzero=rmzero)

        # handle kwargs
        max_nrow = max(gdf.nrow for gdf in rawdata_plotdata_dict.values())
        dot_kwargs = (
            {
                'color': 'black', 
                'markersize': 0.3, 
                'alpha': libgenomeplot.calc_dot_alpha_baf(max_nrow),
            }
            | dot_kwargs
        )
        line_segmean_kwargs = (
            {'color': 'tab:blue', 'linewidth': 2, 'alpha': 0.8}
            | line_segmean_kwargs
        )
        line_corr_segmean_kwargs = (
            {'color': 'tab:green', 'linewidth': 2, 'alpha': 0.8}
            | line_corr_segmean_kwargs
        )

        for baf_idx, ax in ax_dict.items():
            # draw raw data
            self.genomeplotter.draw_dots(
                ax, 
                plotdata=rawdata_plotdata_dict[baf_idx],
                y_colname=baf_idx, 
                plot_kwargs=dot_kwargs,
                draw_common=False, 
            )

            # draw segment
            if draw_segment:
                y_colname = segment_dict[baf_idx].get_distinfo_colname(segment_value_type)
                self.genomeplotter.draw_hlines(
                    ax, 
                    plotdata=segment_plotdata_dict[baf_idx],
                    y_colname=y_colname, 
                    plot_kwargs=line_segmean_kwargs,
                    draw_common=False, 
                )

            # set axes styles
            if baf_idx in ylabels.keys():
                ylabel = ylabels[baf_idx]
            else:
                ylabel = f'BAF ({baf_idx})'
            self.setup_yaxis(
                ax, 
                ylabel=ylabel, 
                ylabel_kwargs=ylabel_kwargs, 
                ymin=ymin, 
                ymax=ymax,
                yticks=np.round(np.arange(0, 0.6, 0.1), 1),
            )
            self.genomeplotter.draw_common(
                ax, 
                **draw_common_kwargs,
            )

            # draw excluded region
            self.draw_excluded_region(ax)

        return fig, ax_dict

    @deco_plotter
    def draw_depth_ax(
        self, 
        ax=None, 

        # fig generation params
        title=None,
        suptitle_kwargs=dict(),
        subplots_kwargs=dict(
            figsize=(30, 5),
        ),

        # normalized vs raw depth
        normalized=True,

        # drawing kwargs
        dot_kwargs=dict(),
        line_segmean_kwargs=dict(),

        # segment drawing
        draw_segment=True,

        # axes setting
        ylabel=None,
        ylabel_kwargs=dict(),
        ymax=None,
        ymin=None,
        draw_common_kwargs=dict(n_xlabel=20),

        # drawing region configuration
        genomeplotter_kwargs=None,
    ):
        """Args:
            ax_dict: dict (keys: index of baf (0, 1, ...) ; values: ax)
            ylabels: dict (keys: index of baf (0, 1, ...) ; values: label text)
        """
        if ax is None:
            fig, ax = plt.subplots(1, 1, **subplots_kwargs)
            if title is None:
                title = self.get_default_title()
            fig.suptitle(title, **suptitle_kwargs)
        else:
            fig = None

        # prepare plotdata
        rawdata_plotdata = self.make_depth_rawdata_plotdata()
        if draw_segment:
            segment_plotdata = self.make_depth_segment_plotdata()

        # set drawing parameters
        dot_kwargs = (
            {
                'color': 'black', 
                'markersize': 0.3, 
                'alpha': libgenomeplot.calc_dot_alpha_baf(rawdata_plotdata.nrow),
            }
            | dot_kwargs
        )
        line_segmean_kwargs = (
            {'color': 'tab:blue', 'linewidth': 2, 'alpha': 0.8}
            | line_segmean_kwargs
        )

        # draw raw data
        if normalized:
            raw_y_colname = DepthRawDF.norm_depth_colname
        else:
            raw_y_colname = DepthRawDF.depth_colname

        self.genomeplotter.draw_dots(
            ax, 
            plotdata=rawdata_plotdata,
            y_colname=raw_y_colname, 
            plot_kwargs=dot_kwargs,
            draw_common=False, 
        )

        # draw segment
        if draw_segment:
            segment_y_colname = DepthSegDF.norm_depth_mean_colname
            self.genomeplotter.draw_hlines(
                ax, 
                plotdata=segment_plotdata,
                y_colname=segment_y_colname, 
                plot_kwargs=line_segmean_kwargs,
                draw_common=False, 
            )

        # set axes styles
        if ylabel is None:
            prefix = 'Normalized' if normalized else 'Raw'
            ylabel = f'{prefix} Depth'
        if ymax is None:
            #ymax = np.nanmean(rawdata_plotdata[raw_y_colname]) * 2
            ymax = np.nanquantile(rawdata_plotdata[raw_y_colname], 0.999)
        if ymin is None:
            ymin = -ymax * 0.01
        yticks = np.round(
            np.linspace(0, ymax, 10), 
            (2 if normalized else 1),
        )

        self.setup_yaxis(
            ax, 
            ylabel=ylabel, 
            ylabel_kwargs=ylabel_kwargs, 
            ymin=ymin, 
            ymax=ymax,
            yticks=yticks,
        )
        self.genomeplotter.draw_common(
            ax, 
            **draw_common_kwargs,
        )

        # draw excluded region
        self.draw_excluded_region(ax)

        return fig, ax

    def draw_excluded_region(self, ax):
        self.genomeplotter.draw_bgcolors(
            ax,
            data=self.get_excluded_region(),
            colors='magenta',
            draw_common=False,
        )

    def draw_data_ax(self):
        pass

    def draw_solution_ax(self):
        pass

    ######################
    # ax drawing helpers #
    ######################

    @classmethod
    def setup_yaxis(
        cls,
        ax, 
        ylabel=None, 
        ylabel_kwargs=dict(), 
        ymin=None, 
        ymax=None,
        yticks=None,
    ):
        if ylabel is not None:
            ax.set_ylabel(ylabel, **ylabel_kwargs)
        ax.set_ylim(ymin, ymax)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks, size=cls.get_yticklabel_size(yticks))

    @staticmethod
    def get_yticklabel_size(yticks):
        return min((200 /len(yticks)), 10)

    def get_default_title(self):
        sampleid = self.simple_attrs["sampleid"]
        is_female = self.simple_attrs["is_female"]
        return f'sampleid={sampleid}, is_female={is_female}'

    def get_baf_indexes(self):
        return libbaf.get_baf_indexes(self.simple_attrs["ploidy"])











#########################################################################



######################
# CNV analyzer class #
######################

PLOTTER_DECORATOR_REGION_ARG_NAMES = (
    'region_chroms',
    'region_start0s',
    'region_end0s',
    'weights',
    'region_gaps',
)

def plotter_decorator(func):
    sig = inspect.signature(func)
    new_sig = sig.replace(
        parameters=itertools.chain(
            sig.parameters.values(),
            (
                inspect.Parameter(
                    x, inspect.Parameter.POSITIONAL_OR_KEYWORD, default=None,
                ) for x in PLOTTER_DECORATOR_REGION_ARG_NAMES
            ),
        )
    )

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        ba = new_sig.bind(*args, **kwargs)
        ba.apply_defaults()
        argdict = ba.arguments

        fig, axd = plotter_decorator_core(
            basemethod=func, 
            **argdict,
        )

        return fig, axd

    return wrapper


def plotter_decorator_core(
    basemethod, 

    region_chroms=None, 
    region_start0s=None, 
    region_end0s=None, 
    weights=None,
    region_gaps=None,

    **kwargs,
):
    if (
        (region_gaps is not None)
        or (region_chroms is not None)
    ):
        self = kwargs['self']

        # make new genomeplotter
        if region_chroms is None:
            new_region_gdf = self.genomeplotter.region_gdf.copy()
        else:
            new_region_gdf = libgenomeplot.make_new_region_gdf(
                refver=self.refver, 
                chroms=region_chroms, 
                start0s=region_start0s, 
                end0s=region_end0s, 
                weights=weights,
            )

        new_genomeplotter = GenomePlotter(
            self.refver, 
            region_gdf=new_region_gdf,
            region_gaps=region_gaps,
        )

        # switch and run
        old_genomeplotter = self.genomeplotter
        self.genomeplotter = new_genomeplotter

        fig, axd = basemethod(**kwargs)

        self.genomeplotter = old_genomeplotter
    else:
        fig, axd = basemethod(**kwargs)

    return fig, axd

class CNVPlotter:
    def __init__(
        self, 
        refver=handygenome.PARAMS['default_refver'], 
        region_df=None, 
        chroms=None, start0s=None, end0s=None, weights=None,
        region_gaps=None,
    ):
        self.refver = refgenome.standardize_refver(refver)

        region_df, region_gaps = libgenomeplot.handle_region_args(
            refver, 
            region_df=region_df, 
            chroms=chroms, 
            start0s=start0s, 
            end0s=end0s, 
            weights=weights,
            region_gaps=region_gaps,
        )

        self.genomeplotter = GenomePlotter(
            refver, region_df=region_df, region_gaps=region_gaps,
        )
        self.data = dict()
        self.default_binsize = 100

    def reset_genomeplotter(
        self, 
        region_df=None, 
        chroms=None, start0s=None, end0s=None, weights=None,
        region_gaps=None,
    ):
        region_df, region_gaps = libgenomeplot.handle_region_args(
            refver=self.refver, 
            region_df=region_df, 
            chroms=chroms, 
            start0s=start0s, 
            end0s=end0s, 
            weights=weights,
            region_gaps=region_gaps,
        )
        self.genomeplotter = GenomePlotter(
            self.refver, region_df=region_df, region_gaps=region_gaps,
        )

    def save_data(self, sampleid, outfile_path):
        with open(outfile_path, 'wb') as outfile:
            pickle.dump(self.data[sampleid], outfile)

    def load_data(self, sampleid, infile_path):
        with open(infile_path, 'rb') as infile:
            self.data[sampleid] = pickle.load(infile)

    ##############
    # mainstream #
    ##############

    @deco.get_deco_num_set_differently(
        ('normal_bam_path', 'normal_depth_path', 'normal_depth_df'), 2, 'lt',
    )
    @deco.get_deco_num_set_differently(
        ('tumor_bam_path', 'tumor_depth_path', 'tumor_depth_df'), 2, 'lt',
    )
    @deco.get_deco_num_set_differently(
        ('germline_vcf_path', 'germline_vafdf_path'), 2, 'lt',
    )
    def add_sample_file_new(
        self, 
        sampleid, 
        is_female,

        germline_vcf_path=None,
        germline_vafdf_path=None,
        vcf_sampleid_tumor=None,
        vcf_sampleid_normal=None,

        *,

        norm_method='plain',

        mode='wgs',
        target_region=None,

        normal_bam_path=None,
        normal_depth_path=None, 
        normal_depth_df=None, 

        tumor_bam_path=None,
        tumor_depth_path=None, 
        tumor_depth_df=None, 

        vcfload_nproc=1,
        #mosdepth_postprocess_kwargs=dict(),

        verbose=True,
    ):
        """Args:
            *_depth_path: mosdepth output file
            germline_vcf_path: germline variant vcf
            vcf_sampleids: tuple of (normal sample id, tumor sample id)
        """
        def depth_loading_helper(
            self, 
            bam_path, 
            depth_path, 
            depth_df, 
            sampleid, 
            sampletype,
        ):
            assert sampletype in ('tumor', 'normal')

            if bam_path is not None:
                LOGGER_INFO.info(f'Loading {sampletype} depth - running mosdepth')
                depth_df = libmosdepth.run_mosdepth(
                    bam_path, 
                    t=8, 
                    use_median=False, 
                    region_bed_path=None, 
                    region_gr=(
                        None
                        if self.data[sampleid]['mode'] == 'wgs' else
                        self.data[sampleid]['target_region']
                    ), 
                    window_size=self.default_binsize, 
                    donot_subset_bam=True,
                    as_gr=False, 
                    load_perbase=False,
                )
            elif depth_path is not None:
                LOGGER_INFO.info(f'Loading {sampletype} depth - reading mosdepth output file')
                depth_df = libmosdepth.load_mosdepth_output(
                    depth_path, depth_colname='mean_depth', as_gr=False,
                )
            elif depth_df is not None:
                LOGGER_INFO.info(f'Loading {sampletype} depth - using the given depth dataframe')
                assert isinstance(depth_df, pd.DataFrame)
                assert set(depth_df.columns) == {'Chromosome', 'Start', 'End', 'mean_depth'}

            if depth_df is not None:
                self.data[sampleid][f'{sampletype}_depth'] = depth_df

        def sanity_check(
            mode, target_region, vcf_sampleid_tumor, vcf_sampleid_normal, germline_vcf_path,
        ):
            # target_region
            assert mode in ('wgs', 'panel')
            if (mode == 'panel') and (target_region is None):
                raise Exception(f'"target_region" must be given when "mode" is "panel"')
            elif (mode == 'wgs') and (target_region is not None):
                raise Exception(f'"target_region" must not be given when "mode" is "wgs"')

            # germline VCF file arguments
            if germline_vcf_path is not None:
                if vcf_sampleid_tumor is None:
                    raise Exception(f'When "germline_vcf_path" is used, "vcf_sampleid_tumor" must be given.')
                if (mode == 'wgs') and (vcf_sampleid_normal is None):
                    raise Exception(f'"vcf_sampleid_normal" must be given when "mode" is "wgs"')

        # main
        sanity_check(mode, target_region, vcf_sampleid_tumor, vcf_sampleid_normal, germline_vcf_path)
        self.data[sampleid] = dict()
        self.data[sampleid]['is_female'] = is_female
        self.data[sampleid]['mode'] = mode

        # load normal
        depth_loading_helper(
            self, 
            normal_bam_path, 
            normal_depth_path, 
            normal_depth_df,
            sampleid, 
            'normal', 
        )

        # load tumor
        depth_loading_helper(
            self, 
            tumor_bam_path, 
            tumor_depth_path, 
            tumor_depth_df,
            sampleid, 
            'tumor', 
        )

        # set target_region
        self.set_target_region(sampleid, mode, target_region)

        # normal mean ploidy
        self.set_normal_mean_ploidy(sampleid)

        # load germline vcf
        if germline_vcf_path is not None:
            LOGGER_INFO.info('Loading tumor germline vcf')
            self.load_germline_vcf(
                sampleid=sampleid, 
                vcf_path=germline_vcf_path, 
                vcf_sampleid_tumor=vcf_sampleid_tumor,
                vcf_sampleid_normal=vcf_sampleid_normal,
                logging_lineno=50000,
                nproc=vcfload_nproc,
            )
        elif germline_vafdf_path is not None:
            LOGGER_INFO.info('Loading tumor germline vaf dataframe')
            self.load_germline_vafdf(sampleid, germline_vafdf_path)

        #self.postprocess_bafdf(sampleid)

        # postprocess depths
        self.postprocess_depth(sampleid, verbose=verbose, norm_method=norm_method)

    def load_germline_vafdf(self, sampleid, vafdf_path):
        self.data[sampleid]['original_baf'] = pd.read_csv(
            vafdf_path,
            sep='\t',
            dtype={
                'Chromosome': 'string',  
                'Start': int,   
                'End': int,     
                'vaf_raw_tumor': float,   
                'baf_raw_tumor': float,   
                'vaf_raw_normal': float,  
                'baf_raw_normal': float,
            },
        )

    def postprocess_bafdf(self, sampleid):
        modified_baf = self.data[sampleid]['original_baf'].copy()
        #self.data[sampleid]['baf'] = modified_baf.loc[modified_baf['baf_raw_tumor'] > 0, :]

    def add_bafpeak_to_segment(self, sampleid, bw=1):
        # join segment and raw bafs
        left = self.data[sampleid]['baf_segment']

        right = self.data[sampleid]['original_baf']
        right = right.loc[
            right['baf_raw_tumor'] > 0, 
            ['Chromosome', 'Start', 'End', 'baf_raw_tumor'],
        ]

        joined = pyranges_helper.join(
            left, right, how='left', merge=None, sort=True, refver=self.refver,
        )
        # find peaks
        groupkey = cnvmisc.genome_df_groupkey(joined, refver=self.refver)
        peaks = list()
        for k, v in joined.groupby(groupkey)[['baf_raw_tumor']]:
            baf_values = v['baf_raw_tumor'].to_numpy()
            peaks.append(
                bafcorrection.infer_baf_density(baf_values, bw=bw, rmzero=False)
            )
        # assign values
        assert len(peaks) == left.shape[0], (
            f'The number of groupby groups and segment df row number are different'
        )
        left['baf_segment_peak'] = peaks
        self.data[sampleid]['baf_segment'] = left

    def postprocess_depth(self, sampleid, norm_method='plain', verbose=False):
        logger = (LOGGER_DEBUG if verbose else LOGGER_INFO)

        # get GC df
        if norm_method == 'plain':
            gc_df = None
        else:
            if self.data[sampleid]['mode'] == 'wgs':
                logger.info(f'Getting gc fraction dataframe')
                gc_df = libgcfraction.get_gc_df(
                    self.refver, 
                    self.default_binsize, 
                    coords_as_index=True,
                )
            elif self.data[sampleid]['mode'] == 'panel':
                gc_df = None

        # main
        for key in ('normal', 'tumor'):
            if f'{key}_depth' not in self.data[sampleid].keys():
                continue

            logger.info(f'Beginning postprocess of {key} depth')
            refver_arg = (
                self.refver 
                if self.data[sampleid]['mode'] == 'panel' else 
                None
            )
            output_depth_df, gcbin_average_depths = cnvmisc.postprocess_depth_df(
                self.data[sampleid][f'{key}_depth'],
                refver=refver_arg,
                gc_df=gc_df,
                included_region=self.data[sampleid]['target_region'],
                as_gr=False,
                verbose=verbose,
                norm_method=norm_method,
                #add_norm_depth=True,
            )

            self.data[sampleid][f'{key}_depth'] = output_depth_df
            self.data[sampleid][f'{key}_gcdata'] = gcbin_average_depths

    def set_depthratio(self, sampleid):
        depthratio_df = cnvmisc.make_depth_ratio(
            self.data[sampleid]['tumor_depth'], 
            self.data[sampleid]['normal_depth'],
            #make_depthratio_mystyle=False,
            #make_depthratio_plain=True,
            as_gr=False,
        )
        depthratio_df.rename(
            columns={
                #'depth_ratio_sequenzastyle': 'depthratio_raw_seqzstyle',
                #'depth_ratio_plain': 'depthratio_raw_plain',
                'depthratio': 'depthratio_raw',
            }, 
            inplace=True,
        )
        self.data[sampleid]['depthratio'] = depthratio_df

    def upscale_preprocessing(self, input_df):
        result = input_df.copy()
        annot_cols = cnvmisc.get_genome_df_annotcols(input_df)
        result.loc[result['excluded'], annot_cols] = np.nan
        result.drop('excluded', axis=1, inplace=True)
        return result

    def upscale_depthratio(self, sampleid, binsize=1000):
        #input_df = self.upscale_preprocessing(self.data[sampleid]['depthratio'])
        input_df = self.data[sampleid]['depthratio']
        self.data[sampleid]['depthratio_upscaled'] = cnvmisc.upsize_depth_df_bin(
            input_df, 
            size=binsize, 
            refver=self.refver,
        )

    def upscale_depth(self, sampleid, binsize=1000, do_normal=True, do_tumor=True):
        if do_normal:
            self.data[sampleid]['normal_depth_upscaled'] = cnvmisc.upsize_depth_df_bin(
                self.data[sampleid]['normal_depth'], 
                size=binsize, 
                refver=self.refver,
            )

        if do_tumor:
            self.data[sampleid]['tumor_depth_upscaled'] = cnvmisc.upsize_depth_df_bin(
                self.data[sampleid]['tumor_depth'], 
                size=binsize, 
                refver=self.refver,
            )

    def make_segments(
        self,
        sampleid,
        winsorize=False,
        depthratio_gamma=None,
        depthratio_kmin=None,
        baf_gamma=100,
        baf_kmin=None,
        verbose=False,
        segment_baf_cutoff=0.1,

        bafcorrection_cutoff=None,
        bafcorrector=None,

        bw=1,
    ):
        depthratio_df = (
            self.data[sampleid]['depthratio_upscaled']
            if 'depthratio_upscaled' in self.data[sampleid] else
            self.data[sampleid]['depthratio']
        )
        depth_segment, baf_segment = _make_segments_main(
            depthratio_df=depthratio_df,
            mode=self.data[sampleid]['mode'],
            refver=self.refver,
            winsorize=winsorize,
            depthratio_gamma=depthratio_gamma,
            depthratio_kmin=depthratio_kmin,
            baf_gamma=baf_gamma,
            baf_kmin=baf_kmin,
            verbose=verbose,

            baf_df=self.data[sampleid]['original_baf'],
            target_region=self.data[sampleid]['target_region'],
            baf_cutoff=segment_baf_cutoff,
        )
        self.data[sampleid]['depthratio_segment'] = depth_segment
        self.data[sampleid]['baf_segment'] = baf_segment

        self.make_segments_postprocess(
            sampleid, 
            bw=bw,
            bafcorrection_cutoff=bafcorrection_cutoff,
            bafcorrector=bafcorrector,
        )

    def make_segments_postprocess(
        self, sampleid, bw=1,
        bafcorrection_cutoff=None,
        bafcorrector=None,
    ):
        # add baf peak
        self.add_bafpeak_to_segment(sampleid, bw=bw)

        # make merged segments
        self.data[sampleid]['merged_segment'] = self.make_merged_segment(
            sampleid,
            self.data[sampleid]['depthratio_segment'],
        )
        
        # add corrected baf
        self.add_corrected_baf(
            sampleid, 
            round_cutoff=bafcorrection_cutoff,
            bafcorrector=bafcorrector,
        )

    def add_corrected_baf(
        self, 
        sampleid, 
        round_cutoff=None,
        bafcorrector=None,
    ):
        if bafcorrector is None:
            bafcorrector = bafcorrection.load_bafcorrect_func(x_cutoff=round_cutoff)

        for df in (
            self.data[sampleid]['baf_segment'], 
            self.data[sampleid]['merged_segment'], 
        ):
            df['corrected_baf_segment_mean'] = bafcorrector(df['baf_segment_peak'].to_numpy())

    def get_cp_from_twodata(
        self, sampleid, depthratio1, CNt1, depthratio2, CNt2,
    ):
        return cnvmisc.get_cp_from_twodata(
            depthratio1, 
            CNt1, 
            depthratio2, 
            CNt2, 
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'], 
            CNn=2, 
            tumor_avg_depth_ratio=1, 
            normal_avg_depth_ratio=1,
        )

    def calculate_tumor_ploidy(self, sampleid, segment_df):
        if self.data[sampleid]['mode'] == 'wgs':
            weights = segment_df['End'] - segment_df['Start']
        else:
            stripped_segment_df = segment_df.loc[:, ['Chromosome', 'Start', 'End']].copy()
            all_indexes = list(range(stripped_segment_df.shape[0]))  # index of segments
            stripped_segment_df['index'] = all_indexes

            target_region = cnvmisc.arg_into_gr(self.data[sampleid]['target_region'])
            index_annotated_targetregion_df = cnvmisc.annotate_region_with_segment(  
                # each region is annotated with corresponding segment index
                target_region[[]],
                stripped_segment_df,
                as_gr=False,
            )

            index_annotated_targetregion_df['length'] = (
                index_annotated_targetregion_df['End']
                - index_annotated_targetregion_df['Start']
            )
            weights_dict = index_annotated_targetregion_df.loc[
                :, ['length', 'index']
            ].groupby('index').sum().to_dict()['length']
            weights = [
                (weights_dict[x] if x in weights_dict else 0)
                for x in all_indexes   
            ]

        return np.average(segment_df['CNt'], weights=weights)

    @plotter_decorator
    def plot_beforecp(
        self, 
        sampleid, 
        figsize=None, 
        hspace=None,
        draw_invalid_regions=False, 
        use_saved_plotdata=False,
        use_merged_segment=True,

        draw_depthratio_hist=True,
        #rm_haploid_from_hist=True,

        depthratio_hist_threshold=None,
        depthratio_hist_bw=None,

        draw_depth=False,
        is_rawdepth=True,
        depth_binsize=10000,
        depth_ymax=None,

        n_xlabel=None,
        depthratio_ymax=None,

        depthratio_dot_kwargs=dict(),
        depthratio_line_segmean_kwargs=dict(),
        #depthratio_line_predict_kwargs=dict(),

        depthratio_hist_annotate_kwargs=dict(),
        depthratio_hist_plot_kwargs=dict(),

        depth_dot_kwargs=dict(),

        baf_dot_kwargs=dict(),
        baf_line_segmean_kwargs=dict(),
        baf_line_corr_segmean_kwargs=dict(),
        #baf_line_predict_kwargs=dict(),

        draw_upscaled_depthratio=True,
        draw_upscaled_depth=True,
    ):
        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
        self.make_plotdata_basic(sampleid, draw_upscaled_depthratio)
        if draw_depth:
            self.make_plotdata_fordepth(sampleid, use_upscaled=draw_upscaled_depth, binsize=depth_binsize)
        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')

        fig, axd, bottom_axkey = self.make_axd(
            figsize=figsize, 
            hspace=hspace, 
            draw_depthratio_hist=draw_depthratio_hist, 
            draw_solution=False,
            draw_tumor_baf=True,
            draw_depthratio=True,
            draw_normal_baf=draw_depth,
            draw_tumor_depth=draw_depth,
            draw_normal_depth=draw_depth,
        )

        fig.suptitle(
            f'sample_id={sampleid}, is_female={self.data[sampleid]["is_female"]}',
            fontsize=20,
        )

        self.draw_baf_ax(
            sampleid, 
            axd['baf'], 

            use_merged_segment=True,

            draw_predicted=False,
            draw_corrected=True,
            draw_segmean=True,
            n_xlabel=(n_xlabel if 'baf' == bottom_axkey else None),

            is_tumor=True,

            mark_unfit_regions=False,

            dot_kwargs=baf_dot_kwargs,
            line_segmean_kwargs=baf_line_segmean_kwargs,
            line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,
        )
        self.draw_depthratio_ax(
            sampleid, 
            axd['depthratio'], 
            use_merged_segment=True,

            draw_predicted=False,
            draw_segmean=True,
            draw_deviation=False,
            draw_depthratio_peaks=False,

            n_xlabel=(n_xlabel if 'depthratio' == bottom_axkey else None),
            dot_kwargs=depthratio_dot_kwargs,
            line_segmean_kwargs=depthratio_line_segmean_kwargs,
            ymax=depthratio_ymax,

            use_upscaled=draw_upscaled_depthratio,
        )

        if draw_depthratio_hist:
            peak_depthratios = self.draw_depthratio_hist_ax(
                sampleid,
                axd['depthratio_hist'],
                use_merged_segment=True, 
                depth_ylim=axd['depthratio'].get_ylim(),
                #rm_haploid=rm_haploid_from_hist,
                peak_threshold=depthratio_hist_threshold,
                bw=depthratio_hist_bw,
                annotate_kwargs=depthratio_hist_annotate_kwargs,
                plot_kwargs=depthratio_hist_plot_kwargs,
            )

            for y in peak_depthratios:
                axd['depthratio'].axhline(y, color='orange', linewidth=1, alpha=0.6)

        if draw_depth:
            self.draw_depth_bundle(
                sampleid, axd, n_xlabel, is_rawdepth, 
                depth_dot_kwargs, depth_ymax, 
                baf_dot_kwargs, baf_line_segmean_kwargs, baf_line_corr_segmean_kwargs,
                bottom_axkey,
                draw_upscaled_depth,
            )

        return fig, axd

    @plotter_decorator
    def plot_aftercp_freeccf(
        self, 
        sampleid, 
        cellularity,
        ploidy,

        figsize=(30, 16), 
        hspace=None,

        n_xlabel=None,
        depthratio_ymax=None,
        CN_ymax=None,
        subCN_ymax=None,

        depthratio_std_factor=1,

        draw_depth=False,
        is_rawdepth=True,
        depth_binsize=10000,
        depth_ymax=None,

        depthratio_dot_kwargs=dict(),
        depthratio_line_segmean_kwargs=dict(),
        depthratio_line_predict_kwargs=dict(),
        depthratio_line_predict_clonal_kwargs=dict(),

        baf_dot_kwargs=dict(),
        baf_line_segmean_kwargs=dict(),
        baf_line_corr_segmean_kwargs=dict(),
        baf_line_predict_kwargs=dict(),
        baf_line_predict_clonal_kwargs=dict(),

        CN_line_CNt_kwargs=dict(),
        subCN_line_CNt_kwargs=dict(),
        ccf_bar_kwargs=dict(),

        depth_ratio_diff=None,
        baf_diff=0.05,

        Bn=1,

        min_N_CNt_candidates=5,
        N_CNt_candidates_fraction=0.5,

        limited_clonal=True,

        draw_upscaled_depthratio=True,
        draw_upscaled_depth=True,
    ):
        LOGGER_INFO.info(f'Beginning calculation of subclonal solution')
        self.make_CN_solution_freeccf(
            sampleid,
            cellularity,
            ploidy,
            depth_ratio_diff=depth_ratio_diff,
            baf_diff=baf_diff,
            min_N_CNt_candidates=min_N_CNt_candidates,
            N_CNt_candidates_fraction=N_CNt_candidates_fraction,
            limited_clonal=limited_clonal,
        )
        self.add_freeccf_solution_to_segment(sampleid)

        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
        self.make_plotdata_basic(sampleid, draw_upscaled_depthratio)
        if draw_depth:
            self.make_plotdata_fordepth(sampleid, use_upscaled=draw_upscaled_depth, binsize=depth_binsize)
        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')

        fig, axd, bottom_axkey = self.make_axd(
            figsize=figsize, 
            hspace=hspace, 
            draw_depthratio_hist=False, 
            draw_solution=True,
            draw_tumor_baf=True,
            draw_depthratio=True,
            draw_normal_baf=draw_depth,
            draw_tumor_depth=draw_depth,
            draw_normal_depth=draw_depth,
        )

        fig.suptitle(
            ', '.join([
                f'sample_id={sampleid}',
                f'is_female={self.data[sampleid]["is_female"]}',
                f'cellularity={round(cellularity, 3)}',
                f'ploidy={round(ploidy, 3)}',
            ]),
            fontsize=20,
        )

        # depth
        self.draw_depthratio_ax(
            sampleid, 
            axd['depthratio'], 
            use_merged_segment=True,

            draw_predicted=True,
            draw_segmean=True,
            draw_deviation=False,
            draw_depthratio_peaks=False,

            mark_unfit_regions=True,

            cellularity=cellularity,
            ploidy=ploidy,
            draw_integerCN_lines=True,

            std_factor=depthratio_std_factor,
            n_xlabel=n_xlabel,
            dot_kwargs=depthratio_dot_kwargs,
            line_segmean_kwargs=depthratio_line_segmean_kwargs,
            line_predict_kwargs=depthratio_line_predict_kwargs,
            line_predict_clonal_kwargs=depthratio_line_predict_clonal_kwargs,
            ymax=depthratio_ymax,

            use_upscaled=draw_upscaled_depthratio,
        )
        if draw_depth:
            self.draw_depth_bundle(
                sampleid, axd, n_xlabel, is_rawdepth, 
                depth_dot_kwargs, depth_ymax, 
                baf_dot_kwargs, baf_line_segmean_kwargs, baf_line_corr_segmean_kwargs,
                bottom_axkey,
                draw_upscaled_depth,
            )

        # baf
        self.draw_baf_ax(
            sampleid, 
            axd['baf'], 
            use_merged_segment=True,
            draw_predicted=True,
            draw_corrected=True,
            draw_segmean=True,
            n_xlabel=n_xlabel,

            is_tumor=True,
            mark_unfit_regions=True,

            dot_kwargs=baf_dot_kwargs,
            line_segmean_kwargs=baf_line_segmean_kwargs,
            line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,
            line_predict_kwargs=baf_line_predict_kwargs,
            line_predict_clonal_kwargs=baf_line_predict_clonal_kwargs,
        )

        # clonal CN
        self.draw_CN_ax(
            sampleid, 
            axd['clonal_CN'],
            n_xlabel=n_xlabel,
            line_CNt_kwargs=CN_line_CNt_kwargs,
            ymax=CN_ymax,
            draw_CNt=True,
            draw_A=False,
            draw_B=True,
        )
        # subclonal CN
        self.draw_subclonal_CN_ax(
            sampleid, 
            axd['subclonal_CN'],
            n_xlabel=n_xlabel,
            line_CNt_kwargs=subCN_line_CNt_kwargs,
            ymax=subCN_ymax,
            draw_CNt=True,
            draw_A=False,
            draw_B=True,
        )
        # ccf
        self.draw_ccf_ax(
            sampleid,
            axd['ccf'],
            n_xlabel=n_xlabel,
            bar_kwargs=ccf_bar_kwargs,
        )

        return fig, axd

    def show_ccfs(self, sampleid, bandwidth=0.1):
        self.select_fixed_ccfs(sampleid, bandwidth=bandwidth)
        ccf_plotdata = self.data[sampleid]['ccf_plotdata']
        self.show_ccfs_main(
            ccfs=ccf_plotdata['ccfs'],
            lengths=ccf_plotdata['lengths'],
            density=ccf_plotdata['density'],
            peak_values=ccf_plotdata['peak_values'],
        )

    @plotter_decorator
    def plot_aftercp_fixedccf(
        self, 
        sampleid, 
        cellularity,
        ploidy,

        figsize=(30, 16), 
        hspace=None,

        n_xlabel=None,
        depthratio_ymax=None,
        CN_ymax=None,
        subCN_ymax=None,

        depthratio_std_factor=1,

        draw_depth=False,
        is_rawdepth=True,
        depth_binsize=10000,
        depth_ymax=None,

        depthratio_dot_kwargs=dict(),
        depthratio_line_segmean_kwargs=dict(),
        depthratio_line_predict_kwargs=dict(),
        depthratio_line_predict_clonal_kwargs=dict(),

        baf_dot_kwargs=dict(),
        baf_line_segmean_kwargs=dict(),
        baf_line_corr_segmean_kwargs=dict(),
        baf_line_predict_kwargs=dict(),
        baf_line_predict_clonal_kwargs=dict(),

        CN_line_CNt_kwargs=dict(),
        subCN_line_CNt_kwargs=dict(),
        ccf_bar_kwargs=dict(),

        depth_ratio_diff=None,
        baf_diff=0.05,

        min_N_CNt_candidates=5,
        N_CNt_candidates_fraction=0.5,

        ccf_bw=0.1,

        update_plotdata=False,
        CNt_diff_factor=0.1,

        mark_unfit_regions=False,

        limited_clonal=True,

        draw_upscaled_depthratio=True,
        draw_upscaled_depth=True,

        merge_same_chroms=True,
    ):
        LOGGER_INFO.info(f'Beginning calculation of subclonal solution')
        self.make_CN_solution_after_ccfs(
            sampleid,
            cellularity,
            ploidy,
            min_N_CNt_candidates=min_N_CNt_candidates,
            N_CNt_candidates_fraction=N_CNt_candidates_fraction,
            CNt_diff_factor=CNt_diff_factor,
            limited_clonal=limited_clonal,
        )
        self.add_fixedccf_solution_to_segment(sampleid)
        LOGGER_INFO.info(f'Finished calculation of subclonal solution')

        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
        if update_plotdata:
            self.add_solution_to_plotdata(sampleid)
        else:
            self.make_plotdata_basic(sampleid, draw_upscaled_depthratio)
            if draw_depth:
                self.make_plotdata_fordepth(sampleid, use_upscaled=draw_upscaled_depth, binsize=depth_binsize)
        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')

        fig, axd, bottom_axkey = self.make_axd(
            figsize=figsize, 
            hspace=hspace, 
            draw_depthratio_hist=False, 
            draw_solution=True,
            draw_tumor_baf=True,
            draw_depthratio=True,
            draw_normal_baf=draw_depth,
            draw_tumor_depth=draw_depth,
            draw_normal_depth=draw_depth,
        )

        fig.suptitle(
            ', '.join([
                f'sample_id={sampleid}',
                f'is_female={self.data[sampleid]["is_female"]}',
                f'cellularity={round(cellularity, 3)}',
                f'ploidy={round(ploidy, 3)}',
            ]),
            fontsize=20,
        )

        # depth
        self.draw_depthratio_ax(
            sampleid, 
            axd['depthratio'], 
            use_merged_segment=True,

            draw_predicted=True,
            draw_segmean=True,
            draw_deviation=False,
            draw_depthratio_peaks=False,

            mark_unfit_regions=mark_unfit_regions,

            cellularity=cellularity,
            ploidy=ploidy,
            draw_integerCN_lines=True,

            std_factor=depthratio_std_factor,
            n_xlabel=(n_xlabel if 'depthratio' == bottom_axkey else None),
            dot_kwargs=depthratio_dot_kwargs,
            line_segmean_kwargs=depthratio_line_segmean_kwargs,
            line_predict_kwargs=depthratio_line_predict_kwargs,
            line_predict_clonal_kwargs=depthratio_line_predict_clonal_kwargs,
            ymax=depthratio_ymax,

            use_upscaled=draw_upscaled_depthratio,

            merge_same_chroms=merge_same_chroms,
        )
        if draw_depth:
            self.draw_depth_bundle(
                sampleid, axd, n_xlabel, is_rawdepth, 
                depth_dot_kwargs, depth_ymax, 
                baf_dot_kwargs, baf_line_segmean_kwargs, baf_line_corr_segmean_kwargs,
                bottom_axkey,
                draw_upscaled_depth,
            )

        # baf
        self.draw_baf_ax(
            sampleid, 
            axd['baf'], 
            use_merged_segment=True,
            draw_predicted=True,
            draw_corrected=True,
            draw_segmean=True,
            n_xlabel=(n_xlabel if 'baf' == bottom_axkey else None),

            is_tumor=True,
            mark_unfit_regions=mark_unfit_regions,

            dot_kwargs=baf_dot_kwargs,
            line_segmean_kwargs=baf_line_segmean_kwargs,
            line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,
            line_predict_kwargs=baf_line_predict_kwargs,
            line_predict_clonal_kwargs=baf_line_predict_clonal_kwargs,

            merge_same_chroms=merge_same_chroms,
        )

        # clonal CN
        self.draw_CN_ax(
            sampleid, 
            axd['clonal_CN'],
            n_xlabel=(n_xlabel if 'clonal_CN' == bottom_axkey else None),
            line_CNt_kwargs=CN_line_CNt_kwargs,
            ymax=CN_ymax,
            draw_CNt=True,
            draw_A=False,
            draw_B=True,

            merge_same_chroms=merge_same_chroms,
        )
        # subclonal CN
        self.draw_subclonal_CN_ax(
            sampleid, 
            axd['subclonal_CN'],
            n_xlabel=(n_xlabel if 'subclonal_CN' == bottom_axkey else None),
            line_CNt_kwargs=subCN_line_CNt_kwargs,
            ymax=subCN_ymax,
            draw_CNt=True,
            draw_A=False,
            draw_B=True,

            merge_same_chroms=merge_same_chroms,
        )
        # ccf
        self.draw_ccf_ax(
            sampleid,
            axd['ccf'],
            n_xlabel=(n_xlabel if 'ccf' == bottom_axkey else None),
            bar_kwargs=ccf_bar_kwargs,

            merge_same_chroms=merge_same_chroms,
        )
        fixed_ccfs = self.data[sampleid]['fixed_ccfs']
        for y in fixed_ccfs:
            axd['ccf'].axhline(y, color='red', linewidth=1)

        return fig, axd

    @plotter_decorator
    def plot_custom(
        self, 
        sampleid, 

        fig=None,
        axd=None,
        figsize=(30, 16), 
        hspace=None,
        height_ratios=None,
        title=None,
        title_size=20,
        title_y=0.95,

        draw_depthratio_hist=False, 
        draw_tumor_baf=False,
        draw_depthratio=False,
        draw_normal_baf=False,
        draw_normal_depth=False,
        draw_tumor_depth=False,
        draw_feature=False,

        tumor_baf_ylabel=None,
        depthratio_ylabel=None,
        normal_baf_ylabel=None,
        normal_depth_ylabel=None,
        tumor_depth_ylabel=None,
        feature_ylabel=None,
        tumor_baf_ylabel_kwargs=dict(),
        depthratio_ylabel_kwargs=dict(),
        normal_baf_ylabel_kwargs=dict(),
        normal_depth_ylabel_kwargs=dict(),
        tumor_depth_ylabel_kwargs=dict(),
        feature_ylabel_kwargs=dict(),

        n_xlabel=None,
        depthratio_ymax=None,
        CN_ymax=None,
        subCN_ymax=None,

        is_rawdepth=True,
        depth_binsize=1000,
        depth_ymax=None,
        use_upscaled_depthratio=True,
        use_upscaled_tumor_depth=True,
        use_upscaled_normal_depth=True,

        depthratio_dot_kwargs=dict(),
        depthratio_line_segmean_kwargs=dict(),
        depthratio_line_predict_kwargs=dict(),
        depthratio_line_predict_clonal_kwargs=dict(),

        baf_dot_kwargs=dict(),
        baf_line_segmean_kwargs=dict(),
        baf_line_corr_segmean_kwargs=dict(),
        baf_line_predict_kwargs=dict(),
        baf_line_predict_clonal_kwargs=dict(),

        feature_text_kwargs=dict(),
        feature_line_kwargs=dict(),

        CN_line_CNt_kwargs=dict(),
        subCN_line_CNt_kwargs=dict(),
        ccf_bar_kwargs=dict(),

        normal_depth_dot_kwargs=dict(),
        tumor_depth_dot_kwargs=dict(),

        feature_df=None,
        draw_feature_label=True,
        #feature_as_dot=False,

        #depth_ratio_diff=None,
        #baf_diff=0.05,

        #min_N_CNt_candidates=5,
        #N_CNt_candidates_fraction=0.5,

        #ccf_bw=0.1,

        #update_plotdata=False,
        #CNt_diff_factor=0.1,

        #mark_unfit_regions=False,

        split_spines=True,
        merge_same_chroms=True,
    ):
        # make plotdata
        LOGGER_INFO.info(f'Beginning conversion of data coordinates into plot coordinates')
        if draw_tumor_baf:
            self.make_tumor_baf_plotdata(sampleid)
        if draw_depthratio:
            self.make_depthratio_plotdata(sampleid, use_upscaled=use_upscaled_depthratio)
        if draw_normal_baf:
            self.make_normal_baf_plotdata(sampleid)
        if draw_normal_depth:
            self.make_depth_plotdata(
                sampleid, 
                is_tumor=False, 
                use_upscaled=use_upscaled_normal_depth, 
                binsize=depth_binsize,
            )
        if draw_tumor_depth:
            self.make_depth_plotdata(
                sampleid, 
                is_tumor=True, 
                use_upscaled=use_upscaled_tumor_depth, 
                binsize=depth_binsize,
            )
        LOGGER_INFO.info(f'Finished conversion of data coordinates into plot coordinates')

        # main
        assert not (
            (fig is None)
            and (axd is not None)
        )
        fig_not_given = (fig is None)
        axd_not_given = (axd is None)
        if axd is None:
            fig, axd, bottom_axkey = self.make_axd(
                figsize=figsize, 
                hspace=hspace, 
                height_ratios=height_ratios,
                draw_depthratio_hist=draw_depthratio_hist, 
                draw_solution=False,
                draw_tumor_baf=draw_tumor_baf,
                draw_depthratio=draw_depthratio,
                draw_normal_baf=draw_normal_baf,
                draw_tumor_depth=draw_tumor_depth,
                draw_normal_depth=draw_normal_depth,
                draw_feature=draw_feature,
                draw_xlabel=(n_xlabel is not None),
                fig=fig,
            )
            
        if fig_not_given:
            if title is None:
                title = f'sample_id={sampleid}, is_female={self.data[sampleid]["is_female"]}'
            fig.suptitle(title, fontsize=title_size, y=title_y)

        if draw_tumor_baf:
            self.draw_baf_ax(
                sampleid, 
                axd['baf'], 
                use_merged_segment=True,

                draw_predicted=False,
                draw_corrected=False,
                draw_segmean=False,

                n_xlabel=n_xlabel,

                is_tumor=True,

                mark_unfit_regions=False,

                dot_kwargs=baf_dot_kwargs,
                line_segmean_kwargs=baf_line_segmean_kwargs,
                line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,

                split_spines=split_spines,

                ylabel=tumor_baf_ylabel,
                ylabel_kwargs=tumor_baf_ylabel_kwargs,

                modify_ax=axd_not_given,
                merge_same_chroms=merge_same_chroms,
            )
        if draw_depthratio:
            self.draw_depthratio_ax(
                sampleid, 
                axd['depthratio'], 
                use_merged_segment=True,

                draw_predicted=False,
                draw_segmean=False,
                draw_deviation=False,
                draw_depthratio_peaks=False,

                n_xlabel=n_xlabel,
                dot_kwargs=depthratio_dot_kwargs,
                line_segmean_kwargs=depthratio_line_segmean_kwargs,
                ymax=depthratio_ymax,

                split_spines=split_spines,

                ylabel=depthratio_ylabel,
                ylabel_kwargs=depthratio_ylabel_kwargs,

                modify_ax=axd_not_given,
                merge_same_chroms=merge_same_chroms,

                use_upscaled=use_upscaled_depthratio,
            )
        if draw_normal_baf:
            self.draw_baf_ax(
                sampleid, 
                axd['normal_baf'], 

                use_merged_segment=True,

                draw_predicted=False,
                draw_corrected=False,
                draw_segmean=False,

                n_xlabel=n_xlabel,

                is_tumor=False,

                mark_unfit_regions=False,

                dot_kwargs=baf_dot_kwargs,
                line_segmean_kwargs=baf_line_segmean_kwargs,
                line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,

                split_spines=split_spines,

                ylabel=normal_baf_ylabel,
                ylabel_kwargs=normal_baf_ylabel_kwargs,

                modify_ax=axd_not_given,
                merge_same_chroms=merge_same_chroms,
            )
        if draw_normal_depth:
            self.draw_depth_ax(
                sampleid,
                axd['normal_depth'],
                n_xlabel=n_xlabel,
                is_tumor=False,
                is_rawdepth=is_rawdepth,
                dot_kwargs=normal_depth_dot_kwargs,
                ymax=depth_ymax,

                split_spines=split_spines,

                ylabel=normal_depth_ylabel,
                ylabel_kwargs=normal_depth_ylabel_kwargs,

                modify_ax=axd_not_given,
                merge_same_chroms=merge_same_chroms,
            )
        if draw_tumor_depth:
            self.draw_depth_ax(
                sampleid,
                axd['tumor_depth'],
                n_xlabel=n_xlabel,
                is_tumor=True,
                is_rawdepth=is_rawdepth,
                dot_kwargs=tumor_depth_dot_kwargs,
                ymax=depth_ymax,

                split_spines=split_spines,

                ylabel=tumor_depth_ylabel,
                ylabel_kwargs=tumor_depth_ylabel_kwargs,

                modify_ax=axd_not_given,
                merge_same_chroms=merge_same_chroms,
            )
        if draw_feature:
            self.draw_feature_ax(
                axd['feature'],
                feature_df=feature_df,

                n_xlabel=n_xlabel,

                ylabel=feature_ylabel,
                ylabel_kwargs=feature_ylabel_kwargs,

                text_kwargs=feature_text_kwargs,
                line_kwargs=feature_line_kwargs,

                split_spines=split_spines,
                merge_same_chroms=merge_same_chroms,

                #feature_as_dot=feature_as_dot,
                draw_label=draw_feature_label,
            )

        return fig, axd

    ############
    # plotting #
    ############

    def make_axd(
        self, 
        figsize=None, 
        hspace=None, 
        height_ratios=None,
        draw_depthratio_hist=False, 
        draw_solution=False,
        draw_tumor_baf=False,
        draw_depthratio=False,
        draw_normal_baf=False,
        draw_tumor_depth=False,
        draw_normal_depth=False,
        draw_feature=False,
        draw_xlabel=False,
        fig=None,
    ):
        # row names
        row_names = list()
        if draw_solution:
            row_names.extend(['ccf', 'subclonal_CN', 'clonal_CN'])
        if draw_tumor_baf:
            row_names.append('baf')
        if draw_depthratio:
            row_names.append('depthratio')
        if draw_normal_baf:
            row_names.append('normal_baf')
        if draw_normal_depth:
            row_names.append('normal_depth')
        if draw_tumor_depth:
            row_names.append('tumor_depth')
        if draw_feature:
            row_names.append('feature')

        # default figsize
        if figsize is None:
            figsize = (30, 5 * len(row_names))

        # mosaic
        if draw_depthratio_hist:
            depthratio_idx = row_names.index('depthratio')
            mosaic = [
                (
                    [name, 'empty_upper'] 
                    if idx < depthratio_idx else
                    (
                        [name, 'empty_lower'] 
                        if idx > depthratio_idx else
                        [name, 'depthratio_hist']
                    )
                )
                for idx, name in enumerate(row_names)
            ]
        else:
            mosaic = [[name,] for name in row_names]

        # gridspec_kw
        if hspace is None:
            hspace = (0.4 if draw_xlabel else 0.1)
        gridspec_kw = dict(
            hspace=hspace, 
            height_ratios=height_ratios,
        )
        if draw_depthratio_hist:
            gridspec_kw.update(dict(width_ratios=[1, 0.1], wspace=0.02))

        # result
        if fig is None:
            fig, axd = plt.subplot_mosaic(
                mosaic,
                figsize=figsize,
                gridspec_kw=gridspec_kw,
            )
        else:
            axd = fig.subplot_mosaic(
                mosaic,
                figsize=figsize,
                gridspec_kw=gridspec_kw,
            )

        if draw_depthratio_hist:
            if 'empty_upper' in axd:
                axd['empty_upper'].axis('off')
            if 'empty_lower' in axd:
                axd['empty_lower'].axis('off')

        bottom_axkey = row_names[-1]

        return fig, axd, bottom_axkey

    @staticmethod
    def get_yticklabel_size(yticks):
        return min((200 /len(yticks)), 10)

#    def draw_grids(self, ax, ys, line_params=dict(), merge_same_chroms=True):
#        line_params = (
#            dict(color='black', linewidth=0.2, alpha=0.5)
#            | line_params
#        )
#        chroms, start0s, end0s = zip(
#            *self.genomeplotter.cconv.get_chrom_borders(
#                merge_same_chroms=merge_same_chroms,
#            )
#        )
#        for y in ys:
#            ax.hlines(
#                np.repeat(y, len(start0s)), start0s, end0s, 
#                **line_params,
#            )

    #def draw_centromeres(self, ax):
    #    self.genomeplotter.draw_centromeres(ax)

    def draw_centromeres_type2(self, ax):
        self.genomeplotter.draw_centromeres_type2(ax)

    def make_plotdata_basic(self, sampleid, use_upscaled):
        self.make_depthratio_plotdata(sampleid, use_upscaled)
        self.make_tumor_baf_plotdata(sampleid)
        self.make_segment_plotdata(sampleid)

    def make_plotdata_fordepth(self, sampleid, use_upscaled=True, binsize=100000):
        self.make_normal_baf_plotdata(sampleid)

        LOGGER_INFO.info(f'Beginning tumor depth data processing')
        self.make_depth_plotdata(
            sampleid, is_tumor=True, use_upscaled=use_upscaled, binsize=binsize
        )
        LOGGER_INFO.info(f'Finished tumor depth data processing')

        LOGGER_INFO.info(f'Beginning normal depth data processing')
        self.make_depth_plotdata(
            sampleid, is_tumor=False, use_upscaled=use_upscaled, binsize=binsize,
        )
        LOGGER_INFO.info(f'Finished normal depth data processing')

    def make_depth_plotdata(self, sampleid, is_tumor, use_upscaled=True, binsize=1000):
        # set params
        datakey = ('tumor_depth' if is_tumor else 'normal_depth')
        if use_upscaled:
            datakey = datakey + '_upscaled'
        plotdata_key = ('tumor_depth_plotdata' if is_tumor else 'normal_depth_plotdata')
        #y_colname = ('mean_depth' if is_rawdepth else 'sequenza_style_norm_mean_depth')

        # upscale raw depth df
        relevant_chroms = [
            x for x in self.genomeplotter.cconv.totalregion_df['Chromosome']
            if not x.startswith('-')
        ]
        original_df = self.data[sampleid][datakey]
        relevant_chroms_df = original_df.loc[
            original_df['Chromosome'].isin(relevant_chroms), 
            :
        ]

#        if use_upscaled:
#            assert binsize is not None
#            input_df = cnvmisc.upsize_depth_df_bin(
#                relevant_chroms_df, 
#                size=binsize, 
#                refver=self.refver,
#            )
#        else:
#            input_df = relevant_chroms_df

        # turn into plotdata
        self.data[sampleid][plotdata_key] = self.genomeplotter.prepare_plot_data(
            relevant_chroms_df
        )

    def make_depthratio_plotdata(self, sampleid, use_upscaled=True):
        depthratio_df = (
            self.data[sampleid]['depthratio_upscaled']
            if use_upscaled else
            self.data[sampleid]['depthratio']
        )

        relevant_chroms = [
            x for x in self.genomeplotter.cconv.totalregion_df['Chromosome']
            if not x.startswith('-')
        ]
        relevant_chroms_df = depthratio_df.loc[
            depthratio_df['Chromosome'].isin(relevant_chroms), 
            :
        ]

        self.data[sampleid]['depthratio_raw_plotdata'] = (
            self.genomeplotter.prepare_plot_data(relevant_chroms_df)
        )

    def make_tumor_baf_plotdata(self, sampleid):
        bafdf = self.data[sampleid]['original_baf']
        tumor_baf_df = bafdf.loc[bafdf['baf_raw_tumor'] > 0, :]
        self.data[sampleid]['baf_raw_plotdata'] = (
            self.genomeplotter.prepare_plot_data(tumor_baf_df)
        )

    def make_normal_baf_plotdata(self, sampleid):
        bafdf = self.data[sampleid]['original_baf']
        normal_baf_df = bafdf.loc[bafdf['baf_raw_normal'] > 0, :]
        self.data[sampleid]['normal_baf_raw_plotdata'] = (
            self.genomeplotter.prepare_plot_data(normal_baf_df)
        )

    def make_segment_plotdata(self, sampleid):
        self.data[sampleid]['merged_segment_plotdata'] = (
            self.genomeplotter.prepare_plot_data(
                self.data[sampleid]['merged_segment']
            )
        )

    def make_plotdata_aftercp_wobaf(self, sampleid, use_saved_plotdata):
        raw_plot_map = {
            'depthratio_upscaled': 'depthratio_raw_plotdata',
            'merged_segment': 'merged_segment_plotdata',
        }
        def helper(raw_key):
            plot_key = raw_plot_map[raw_key]
            if (
                (not use_saved_plotdata) 
                or (plot_key not in self.data[sampleid])
            ):
                self.data[sampleid][plot_key] = self.genomeplotter.prepare_plot_data(
                    self.data[sampleid][raw_key]
                )

        helper('depthratio_upscaled')
        helper('merged_segment')

    def draw_feature_ax(
        self, 
        ax, 
        feature_df,

        n_xlabel=None,

        ylabel=None,
        ylabel_kwargs=dict(),

        #feature_as_dot=False,
        draw_label=True,

        text_kwargs=dict(),
        line_kwargs=dict(),

        split_spines=True,
        merge_same_chroms=True,

        chromlabel_kwargs=dict(), 
    ):
        if ylabel is None:
            ylabel = 'features'

        ax.set_ylim(0, 1)
        ax.set_ylabel(ylabel, **ylabel_kwargs)
        ax.set_yticks([])
        self.genomeplotter.draw_common(
            ax, 
            n_xlabel=n_xlabel, 
            split_spines=split_spines,
            merge_same_chroms=merge_same_chroms,
            chromlabel_kwargs=chromlabel_kwargs,
        )

        self.genomeplotter.draw_features(
            ax,
            df=feature_df,

            draw_label=draw_label,
            #feature_as_dot=feature_as_dot,

            y_features=None,
            y_labels=None,

            text_kwargs=text_kwargs,
            line_kwargs=line_kwargs,
        )

    def draw_depth_bundle(
        self, 
        sampleid, axd, n_xlabel, is_rawdepth, 
        depth_dot_kwargs, depth_ymax,
        baf_dot_kwargs, baf_line_segmean_kwargs, baf_line_corr_segmean_kwargs,
        bottom_axkey,
        use_upscaled,
    ):
        self.draw_depth_ax(
            sampleid,
            axd['tumor_depth'],
            n_xlabel=(n_xlabel if 'tumor_depth' == bottom_axkey else None),
            is_tumor=True,
            is_rawdepth=is_rawdepth,
            dot_kwargs=depth_dot_kwargs,
            ymax=depth_ymax,

            use_upscaled=use_upscaled,
        )
        self.draw_depth_ax(
            sampleid,
            axd['normal_depth'],
            n_xlabel=(n_xlabel if 'normal_depth' == bottom_axkey else None),
            is_tumor=False,
            is_rawdepth=is_rawdepth,
            dot_kwargs=depth_dot_kwargs,
            ymax=depth_ymax,

            use_upscaled=use_upscaled,
        )
        self.draw_baf_ax(
            sampleid, 
            axd['normal_baf'], 

            use_merged_segment=True,

            draw_predicted=False,
            draw_corrected=False,
            draw_segmean=False,

            n_xlabel=(n_xlabel if 'normal_baf' == bottom_axkey else None),

            is_tumor=False,

            mark_unfit_regions=False,

            dot_kwargs=baf_dot_kwargs,
            line_segmean_kwargs=baf_line_segmean_kwargs,
            line_corr_segmean_kwargs=baf_line_corr_segmean_kwargs,

            make_segment_plotdata=False,
        )

    def draw_depth_ax(
        self,
        sampleid,
        ax,
        n_xlabel,
        is_tumor,
        is_rawdepth,
        dot_kwargs=dict(),
        ymax=None,
        add_color=True,
        split_spines=True,
        ylabel=None,
        ylabel_kwargs=dict(),

        make_raw_plotdata=True,
        use_upscaled=True,

        modify_ax=True,
        merge_same_chroms=True,

        chromlabel_kwargs=dict(),
    ):
        # prepare plotdata
        if make_raw_plotdata:
            self.make_depth_plotdata(sampleid, is_tumor=is_tumor, use_upscaled=use_upscaled)

        # set params
        y_colname = (
            'mean_depth' 
            if is_rawdepth else 
            'norm_mean_depth'
        )
        plotdata_key = ('tumor_depth_plotdata' if is_tumor else 'normal_depth_plotdata')
        plotdata_df = self.data[sampleid][plotdata_key]

        # calc dot alpha
        n_dots = plotdata_df.shape[0]
        default_alpha = calc_dot_alpha_depth(n_dots)

        # handle kwargs
        dot_kwargs = (
            #{'color': 'black', 'markersize': 0.3, 'alpha': 0.05}
            {'color': 'black', 'markersize': 0.3, 'alpha': default_alpha}
            | dot_kwargs
        )

        # draw raw data
        if add_color:
            colors = np.repeat('blue', plotdata_df.shape[0])
            colors[np.where(plotdata_df['excluded'])[0]] = 'red'
            self.genomeplotter.draw_dots_scatter(
                ax, 
                df_plotdata=plotdata_df,
                y_colname=y_colname, 
                color_vals=colors,
                plot_kwargs=dot_kwargs,
            )
        else:
            self.genomeplotter.draw_dots(
                ax, 
                df_plotdata=plotdata_df,
                y_colname=y_colname, 
                plot_kwargs=dot_kwargs,
            )

        # set axes attributes
        if modify_ax:
            if ylabel is None:
                ylabel_1 = ('tumor' if is_tumor else 'normal')
                ylabel_2 = ('raw' if is_rawdepth else 'normalized')
                ylabel = f'{ylabel_1} sample {ylabel_2} depth'
            ax.set_ylabel(ylabel, **ylabel_kwargs)

            if ymax is None:
                ymax = np.nanmean(plotdata_df[y_colname]) * 2
            ax.set_ylim(-ymax * 0.1, ymax)
            roundnum = (1 if is_rawdepth else 2)
            yticks = np.round(np.linspace(0, ymax, 10), roundnum)

            ax.set_yticks(yticks)
            ax.set_yticklabels(yticks)

            self.genomeplotter.draw_common(
                ax, 
                n_xlabel=n_xlabel, 
                split_spines=split_spines, 
                merge_same_chroms=merge_same_chroms,
                chromlabel_kwargs=chromlabel_kwargs,
            )

    def draw_depthratio_ax(
        self, 
        sampleid, 
        ax, 

        use_merged_segment=True, 

        cellularity=None,
        ploidy=None,

        draw_integerCN_lines=False,

        draw_predicted=True,
        draw_segmean=True,
        draw_deviation=False,
        draw_depthratio_peaks=False,

        std_factor=1,

        mark_unfit_regions=False,

        n_xlabel=None,
        dot_kwargs=dict(),
        line_segmean_kwargs=dict(),
        line_predict_kwargs=dict(),
        line_predict_clonal_kwargs=dict(),
        ymax=None,

        split_spines=True,

        ylabel=None,
        ylabel_kwargs=dict(),

        make_raw_plotdata=True,
        make_segment_plotdata=True,
        use_upscaled=True,

        modify_ax=True,
        merge_same_chroms=True,

        chromlabel_kwargs=dict(),
    ):
        # determine if segment information is plotted
        draw_segment_info = any([
            draw_predicted,
            draw_segmean,
        ])

        # prepare plotdata
        if make_raw_plotdata:
            self.make_depthratio_plotdata(sampleid, use_upscaled=use_upscaled)
        if make_segment_plotdata:
            self.make_segment_plotdata(sampleid)

        # calc dot alpha
        raw_plotdata = self.data[sampleid]['depthratio_raw_plotdata']
        n_dots = raw_plotdata.shape[0]
        default_alpha = calc_dot_alpha_depth(n_dots)

        # handle kwargs
        dot_kwargs = (
            {'color': 'black', 'markersize': 0.3, 'alpha': default_alpha}
            | dot_kwargs
        )
        line_segmean_kwargs = (
            {'color': 'tab:blue', 'linewidth': 5, 'alpha': 0.6}
            | line_segmean_kwargs
        )
        line_predict_kwargs = (
            {'color': 'tab:red', 'linewidth': 1.5, 'alpha': 0.6}
            | line_predict_kwargs
        )
        line_predict_clonal_kwargs = (
            {'color': 'tab:orange', 'linewidth': 1.5, 'alpha': 0.6}
            | line_predict_clonal_kwargs
        )

        # draw raw data
        self.genomeplotter.draw_dots(
            ax, 
            y_colname='depthratio_raw', 
            df_plotdata=raw_plotdata, 
            plot_kwargs=dot_kwargs,
        )

        # draw segment information
        if draw_segment_info:
            if use_merged_segment:
                segment_plotdata = self.data[sampleid]['merged_segment_plotdata']
            else:
                segment_plotdata = self.data[sampleid]['depthratio_segment_plotdata']

            if draw_segmean:
                self.genomeplotter.draw_hlines(
                    ax, 
                    df_plotdata=segment_plotdata,
                    y_colname='depthratio_segment_mean', 
                    plot_kwargs=line_segmean_kwargs,
                )
            if draw_predicted:
                self.genomeplotter.draw_hlines(
                    ax, 
                    df_plotdata=segment_plotdata,
                    y_colname='depthratio_predicted', 
                    plot_kwargs=line_predict_kwargs,
                )
            if draw_depthratio_peaks:
                segments_gr = pr.PyRanges(self.data[sampleid]['merged_segment'])
                peak_xs = cnvmisc.get_depthratio_peaks(
                    segments_gr.depthratio_segment_mean, 
                    lengths=segments_gr.lengths(), 
                    limits=(0, 2), 
                    step=0.01, 
                    peak_cutoff=1e8,
                )
            if draw_deviation:
                self.genomeplotter.draw_bgcolors(
                    ax,
                    df_plotdata=segment_plotdata,
                    ymins=(
                        segment_plotdata['depthratio_segment_mean']
                        - std_factor * segment_plotdata['depthratio_segment_std']
                    ),
                    ymaxs=(
                        segment_plotdata['depthratio_segment_mean']
                        + std_factor * segment_plotdata['depthratio_segment_std']
                    ),
                    colors='yellow',
                    plot_kwargs=dict(alpha=0.4),
                )

            if mark_unfit_regions:
                self.draw_unfit_region(ax, segment_plotdata)

        # integer CN lines
        if draw_integerCN_lines:
            assert ((cellularity is not None) and (ploidy is not None))
            self.draw_depthratio_ax_integer_CNs(ax, cellularity, ploidy, sampleid)

        # modify axes
        if modify_ax:
            if ylabel is None:
                ylabel = 'depth ratio'
            ax.set_ylabel(ylabel, **ylabel_kwargs)

            if ymax is None:
                df = self.data[sampleid]['depthratio']
                y_values = df['depthratio_raw'].loc[~df['excluded']].dropna()
                ymax = y_values.quantile(0.999)

                ax.set_ylim(-ymax * 0.1, ymax)
                yticks = np.round(np.arange(0, ymax + 0.25, 0.25), 2)
            else:
                ax.set_ylim(-0.1, ymax)
                yticks = np.round(np.linspace(0, ymax, 8), 2)

            ax.set_yticks(yticks)
            ax.set_yticklabels(yticks, size=self.get_yticklabel_size(yticks))

            self.genomeplotter.draw_common(
                ax,
                n_xlabel=n_xlabel, 
                split_spines=split_spines,
                merge_same_chroms=merge_same_chroms,
                chromlabel_kwargs=chromlabel_kwargs,
            )

    def draw_depthratio_ax_integer_CNs(self, ax, cellularity, ploidy, sampleid):
        plotregion_df_withCN = self.add_CNn_to_segment(
            segment_df=self.genomeplotter.region_df,
            mode=self.data[sampleid]['mode'],
            refver=self.refver,
            is_female=self.data[sampleid]['is_female'],
            target_region=self.data[sampleid]['target_region'],
        )
        plotdata = self.genomeplotter.prepare_plot_data(plotregion_df_withCN)

        integer_depthratios = cnvmisc.theoretical_depth_ratio(
            CNt=np.arange(0, 10, 1)[:, np.newaxis],
            cellularity=cellularity,
            tumor_ploidy=ploidy,
            CNn=plotregion_df_withCN['CNn'].to_numpy()[np.newaxis, :],
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
        )  # ndim == 2
        for ys in integer_depthratios:
            plotdata['integer_depthratio'] = ys
            self.genomeplotter.draw_hlines(
                ax, 
                y_colname='integer_depthratio',
                df_plotdata=plotdata,
                plot_kwargs=dict(linewidth=0.7, zorder=0, color='black', alpha=0.7),
            )

    def draw_unfit_region(self, ax, plotdata):
        self.genomeplotter.draw_bgcolors(
            ax,
            df_plotdata=plotdata.loc[plotdata['polyploid_unfit'], :],
            colors='yellow',
            plot_kwargs=dict(alpha=0.2),
        )
        self.genomeplotter.draw_bgcolors(
            ax,
            df_plotdata=plotdata.loc[plotdata['polyploid_unfit_bafonly'], :],
            colors='green',
            plot_kwargs=dict(alpha=0.2),
        )
        self.genomeplotter.draw_bgcolors(
            ax,
            df_plotdata=plotdata.loc[plotdata['monoploid_unfit'], :],
            colors='blue',
            plot_kwargs=dict(alpha=0.2),
        )

    def draw_depthratio_hist_ax(
        self,
        sampleid,
        ax,
        use_merged_segment, 
        depth_ylim,
        #rm_haploid,
        kde=True,
        bw=None,
        peak_threshold=None,
        annotate_kwargs=dict(),
        plot_kwargs=dict(),
    ):
        """Histogram only includes CNn == 2 positions"""
        # handle kwargs
        annotate_kwargs = (
            dict(
                va='center', ha='left', size=8,
                arrowprops=dict(arrowstyle='-', color='black', linewidth=0.5),
            )
            | annotate_kwargs
        )
        plot_kwargs = (
            #dict(linewidth=0.8)
            dict(facecolor='tab:blue', alpha=1)
            | plot_kwargs
        )

        if use_merged_segment:
            segment_df = self.data[sampleid]['merged_segment']
        else:
            segment_df = self.data[sampleid]['depthratio_segment']

        segment_df = segment_df.loc[segment_df['CNn'] == 2, :]

        depthratio_list = segment_df['depthratio_segment_mean']
        weights = (segment_df['End'] - segment_df['Start'])

        # set ylim
        ax.set_ylim(*depth_ylim)

        if kde:
            peak_values, peak_densities, density = cnvmisc.get_density_peaks(
                depthratio_list, 
                weights=weights, 
                xs=None, 
                threshold=peak_threshold, 
                bw_method=bw, 
                invert=False,
            )
            peak_depthratios = peak_values
            ax_ymin, ax_ymax = ax.get_ylim()
            fill_ys = np.linspace(ax_ymin, ax_ymax, 100)
            fill_xs = density(fill_ys)
        else:
            histpeaks = cnvmisc.get_hist_peaks(
                depthratio_list, 
                weights=weights,
                bins=np.arange(0, depthratio_list.max(), 0.01),
                threshold=peak_threshold,
            )
            peak_depthratios = histpeaks['peak_values']
            fill_ys = histpeaks['bin_midpoints']
            fill_xs = histpeaks['hist']

        # set xlim
        ax.set_xlim(0, max(fill_xs))

        # draw data
        ax.fill_betweenx(y=fill_ys, x1=fill_xs, **plot_kwargs)
        ytext_list = np.linspace(*ax.get_ylim(), len(peak_depthratios))
        for y, ytext in zip(peak_depthratios, ytext_list):
            ax.axhline(y, color='black', linewidth=0.5)
            ax.annotate(
                round(y, 3), 
                (ax.get_xlim()[1], y),
                xytext=(ax.get_xlim()[1] * 1.1, ytext),
                annotation_clip=False,
                **annotate_kwargs,
            )
        ax.set_yticks([])

        return peak_depthratios

    def draw_baf_ax(
        self, 
        sampleid, 
        ax, 

        use_merged_segment=True, 

        draw_corrected=False,
        draw_predicted=False,
        draw_segmean=False,

        n_xlabel=None,

        is_tumor=True,

        mark_unfit_regions=False,

        dot_kwargs=dict(),
        line_segmean_kwargs=dict(),
        line_corr_segmean_kwargs=dict(),
        line_predict_kwargs=dict(),
        line_predict_clonal_kwargs=dict(),

        split_spines=True,

        ylabel=None,
        ylabel_kwargs=dict(),

        make_raw_plotdata=True,
        make_segment_plotdata=True,

        modify_ax=True,
        merge_same_chroms=True,

        chromlabel_kwargs=dict(),
    ):
        # prepare plotdata
        if make_raw_plotdata:
            if is_tumor:
                self.make_tumor_baf_plotdata(sampleid)
            else:
                self.make_normal_baf_plotdata(sampleid)
        if make_segment_plotdata:
            self.make_segment_plotdata(sampleid)

        # find plotdata
        if is_tumor:
            raw_plotdata = self.data[sampleid]['baf_raw_plotdata']
        else:
            raw_plotdata = self.data[sampleid]['normal_baf_raw_plotdata']

        # calc default alpha
        n_dots = raw_plotdata.shape[0]
        default_alpha = calc_dot_alpha_baf(n_dots)

        # handle kwargs
        dot_kwargs = (
            #{'color': 'black', 'markersize': 0.3, 'alpha': 0.01}
            {'color': 'black', 'markersize': 0.3, 'alpha': default_alpha}
            | dot_kwargs
        )
        line_segmean_kwargs = (
            {'color': 'tab:blue', 'linewidth': 2, 'alpha': 0.4}
            | line_segmean_kwargs
        )
        line_corr_segmean_kwargs = (
            {'color': 'tab:green', 'linewidth': 2, 'alpha': 0.4}
            | line_corr_segmean_kwargs
        )
        line_predict_kwargs = (
            {'color': 'tab:red', 'linewidth': 2, 'alpha': 0.4}
            | line_predict_kwargs
        )
        line_predict_clonal_kwargs = (
            {'color': 'tab:orange', 'linewidth': 2, 'alpha': 0.4}
            | line_predict_clonal_kwargs
        )

        # draw raw data
        if is_tumor:
            self.genomeplotter.draw_dots(
                ax, 
                df_plotdata=raw_plotdata, 
                y_colname='baf_raw_tumor', 
                plot_kwargs=dot_kwargs,
            )
        else:
            self.genomeplotter.draw_dots(
                ax, 
                df_plotdata=raw_plotdata, 
                y_colname='baf_raw_normal', 
                plot_kwargs=dot_kwargs,
            )

        # draw segment information
        draw_segment_info = any([
            draw_corrected,
            draw_predicted,
            draw_segmean,
        ])

        if draw_segment_info:
            if is_tumor:
                segment_plotdata = self.data[sampleid]['merged_segment_plotdata']
            else:
                segment_plotdata = self.data[sampleid]['normal_baf_segment_plotdata']

            if draw_segmean:
                self.genomeplotter.draw_hlines(
                    ax, 
                    df_plotdata=segment_plotdata,
                    y_colname='baf_segment_peak', 
                    #y_colname='baf_segment_mean', 
                    plot_kwargs=line_segmean_kwargs,
                )
            if draw_corrected:
                self.genomeplotter.draw_hlines(
                    ax, 
                    df_plotdata=segment_plotdata,
                    y_colname='corrected_baf_segment_mean', 
                    plot_kwargs=line_corr_segmean_kwargs,
                )

            if draw_predicted:
                self.genomeplotter.draw_hlines(
                    ax, 
                    df_plotdata=segment_plotdata,
                    y_colname='baf_predicted', 
                    plot_kwargs=line_predict_kwargs,
                )

            if mark_unfit_regions:
                self.draw_unfit_region(ax, segment_plotdata)

        # set axes attributes
        if ylabel is None:
            ylabel = 'tumor baf' if is_tumor else 'normal baf'
        ax.set_ylabel(ylabel, **ylabel_kwargs)
        ax.set_ylim(-0.6 * 0.1, 0.6)

        yticks = np.round(np.arange(0, 0.6, 0.1), 1)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks, size=self.get_yticklabel_size(yticks))

        if modify_ax:
            self.genomeplotter.draw_common(
                ax, 
                n_xlabel=n_xlabel, 
                split_spines=split_spines, 
                merge_same_chroms=merge_same_chroms,
                chromlabel_kwargs=chromlabel_kwargs,
            )

    def draw_CN_ax(
        self, 
        sampleid, 
        ax,
        n_xlabel=None,
        line_A_kwargs=dict(),
        line_B_kwargs=dict(),
        line_CNt_kwargs=dict(),
        ymax=None,
        draw_CNt=True,
        draw_A=True,
        draw_B=True,

        split_spines=True,

        ylabel=None,
        ylabel_kwargs=dict(),

        merge_same_chroms=True,
        chromlabel_kwargs=dict(),
    ):
        # handle kwargs
        line_CNt_kwargs = (
            {'color': 'black'}
            | line_CNt_kwargs
        )
        line_A_kwargs = (
            {'color': 'red'}
            | line_A_kwargs
        )
        line_B_kwargs = (
            {'color': 'blue'}
            | line_B_kwargs
        )

        # draw data
        if draw_CNt:
            self.genomeplotter.draw_hlines(
                ax, 
                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
                y_colname='clonal_CNt', 
                offset=0.1,
                plot_kwargs=line_CNt_kwargs,
            )
        if draw_A:
            self.genomeplotter.draw_hlines(
                ax, 
                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
                y_colname='clonal_A', 
                offset=0,
                plot_kwargs=line_A_kwargs,
            )
        if draw_B:
            self.genomeplotter.draw_hlines(
                ax, 
                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
                y_colname='clonal_B', 
                offset=-0.1,
                plot_kwargs=line_B_kwargs,
            )

        # set axes attributes
        if ylabel is None:
            ylabel = 'tumor clonal copy number'
        ax.set_ylabel(ylabel, **ylabel_kwargs)

        if ymax is None:
            df = self.data[sampleid]['merged_segment']
            ymax = df['clonal_CNt'].quantile(0.99)
        else:
            pass

        max_ticknum = 15
        step = np.ceil(ymax / 15).astype(int)
        yticks = np.arange(0, ymax, step).astype(int)

        ax.set_ylim(-0.5, ymax)
        ax.set_yticks(yticks, minor=False)
        #ax.set_yticks(np.arange(0, ymax, 1).astype(int), minor=True)
        ax.set_yticklabels(yticks, size=self.get_yticklabel_size(yticks))

        self.genomeplotter.draw_common(
            ax, 
            n_xlabel=n_xlabel, 
            split_spines=split_spines, 
            merge_same_chroms=merge_same_chroms,
            chromlabel_kwargs=chromlabel_kwargs,
        )

    def draw_subclonal_CN_ax(
        self, 
        sampleid, 
        ax,
        n_xlabel=None,
        line_A_kwargs=dict(),
        line_B_kwargs=dict(),
        line_CNt_kwargs=dict(),
        ymax=None,
        draw_CNt=True,
        draw_A=True,
        draw_B=True,

        split_spines=True,

        ylabel=None,
        ylabel_kwargs=dict(),

        merge_same_chroms=True,
        chromlabel_kwargs=dict(),
    ):
        # handle kwargs
        line_CNt_kwargs = (
            {'color': 'black'}
            | line_CNt_kwargs
        )
        line_A_kwargs = (
            {'color': 'red'}
            | line_A_kwargs
        )
        line_B_kwargs = (
            {'color': 'blue'}
            | line_B_kwargs
        )

        # draw CNt
        if draw_CNt:
            self.genomeplotter.draw_hlines(
                ax, 
                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
                y_colname='subclonal_CNt', 
                offset=0.1,
                plot_kwargs=line_CNt_kwargs,
            )
        if draw_A:
            self.genomeplotter.draw_hlines(
                ax, 
                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
                y_colname='subclonal_A', 
                offset=0,
                plot_kwargs=line_A_kwargs,
            )
        if draw_B:
            self.genomeplotter.draw_hlines(
                ax, 
                df_plotdata=self.data[sampleid]['merged_segment_plotdata'],
                y_colname='subclonal_B', 
                offset=-0.1,
                plot_kwargs=line_B_kwargs,
            )

        # set axes attributes
        if ylabel is None:
            ylabel = 'tumor subclonal copy number'
        ax.set_ylabel(ylabel, **ylabel_kwargs)

        if ymax is None:
            df = self.data[sampleid]['merged_segment']
            ymax = df['subclonal_CNt'].quantile(0.99)
        else:
            pass

        max_ticknum = 15
        step = np.ceil(ymax / 15).astype(int)
        yticks = np.arange(0, ymax, step).astype(int)

        ax.set_ylim(-0.5, ymax)
        ax.set_yticks(yticks, minor=False)
        #ax.set_yticks(np.arange(0, ymax, 1).astype(int), minor=True)
        ax.set_yticklabels(yticks, size=self.get_yticklabel_size(yticks))

        self.genomeplotter.draw_common(
            ax, 
            n_xlabel=n_xlabel, 
            split_spines=split_spines, 
            merge_same_chroms=merge_same_chroms,
            chromlabel_kwargs=chromlabel_kwargs,
        )

    def draw_ccf_ax(
        self, 
        sampleid, 
        ax,
        n_xlabel=None,
        bar_kwargs=dict(),

        split_spines=True,

        ylabel=None,
        ylabel_kwargs=dict(),

        merge_same_chroms=True,
        chromlabel_kwargs=dict(),

        mark_clonal_region=False,
    ):
        bar_kwargs = (
            dict()
            | bar_kwargs
        )

        plotdata = self.data[sampleid]['merged_segment_plotdata']

        self.genomeplotter.draw_bars(
            ax, 
            y_colname='ccf', 
            df_plotdata=plotdata,
            plot_kwargs=bar_kwargs,
        )

        if ylabel is None:
            ylabel = 'subclonal fraction'
        ax.set_ylabel(ylabel, **ylabel_kwargs)

        ax.set_ylim(-0.05, 1.05)
        yticks = np.round(np.arange(0, 1.1, 0.1), 1)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)

        if mark_clonal_region:
            self.genomeplotter.draw_bgcolors(
                ax,
                df_plotdata=plotdata.loc[plotdata['ccf'].isna(), :],
                colors='green',
                plot_kwargs=dict(alpha=0.1),
            )

        self.genomeplotter.draw_common(
            ax, 
            n_xlabel=n_xlabel, 
            split_spines=split_spines, 
            merge_same_chroms=merge_same_chroms,
            chromlabel_kwargs=chromlabel_kwargs,
        )

    def show_ccfs_main(
        self,
        ccfs,
        lengths,
        density,
        peak_values,
        draw_dots=False,
    ):
        fig, ax = plt.subplots(figsize=(30, 5))
        ax.hist(
            ccfs, 
            range=(0, 1), 
            bins=50, 
            weights=lengths,
            density=True,
        )

        density_xs = np.arange(0, 1, 0.01)
        density_ys = density(density_xs)
        ax.plot(density_xs, density_ys)

        for x in peak_values:
            ax.axvline(x, color='red')

        if draw_dots:
            ylim = ax.get_ylim()
            dot_xs = ccfs
            dot_ys = scipy.stats.uniform.rvs(loc=ylim[0], scale=(ylim[1] - ylim[0]), size=len(ccfs))
            colors = fit.labels_
            ax.scatter(dot_xs, dot_ys, c=colors, alpha=0.7, s=4)

        return fig, ax





    ###############################
    # helpers of segment creation #
    ###############################

#    def make_depth_segment(
#        self,
#        sampleid,
#        winsorize=False,
#        gamma=None,
#        kmin=None,
#        verbose=True,
#    ):
#        self.data[sampleid]['depthratio_segment'] = _make_depth_segment(
#            self.data[sampleid]['depthratio_upscaled'],
#            mode=self.data[sampleid]['mode'],
#            refver=self.refver,
#            winsorize=winsorize,
#            gamma=gamma,
#            kmin=kmin,
#            verbose=verbose,
#        )
#
#    def make_baf_segment(
#        self,
#        sampleid,
#        winsorize=False,
#        gamma=None,
#        kmin=None,
#        verbose=False,
#        bafcorrector=bafcorrection.load_bafcorrect_func(),
#    ):
#        self.data[sampleid]['baf_segment'] = _make_baf_segment(
#            baf_df=self.data[sampleid]['baf'],
#            target_region=self.data[sampleid]['target_region'],
#            mode=self.data[sampleid]['mode'],
#            refver=self.refver,
#            winsorize=winsorize,
#            gamma=gamma,
#            kmin=kmin,
#            verbose=verbose,
#            bafcorrector=bafcorrector,
#        )

    def make_merged_segment(self, sampleid, depthratio_segment):
        merged_segment = pyranges_helper.isec_union(
            depthratio_segment,
            self.data[sampleid]['baf_segment'],
        )
        merged_segment = pyranges_helper.join(
            merged_segment, 
            depthratio_segment,
            how='left',
            merge=None,
            find_nearest=True,
            sort=True,
            refver=self.refver,
        )
        merged_segment = pyranges_helper.join(
            merged_segment, 
            self.data[sampleid]['baf_segment'],
            how='left',
            merge=None,
            find_nearest=True,
            sort=True,
            refver=self.refver,
        )

        # reduce by merging adjacent segments with identical annotation values
        merged_segment = self.deduplicate_merged_segment(merged_segment)

        # add CNn
        merged_segment = self.add_CNn_to_segment(
            merged_segment,
            self.data[sampleid]['mode'],
            self.refver,
            self.data[sampleid]['is_female'],
            self.data[sampleid]['target_region'],
        )

        # add std and std/mean ratio
        #self.add_depthratio_std_to_segment(sampleid)

        # fit to target region
        #merged_segment = self.fit_segment_to_targetregion(sampleid, merged_segment)

        return merged_segment

        #self.data[sampleid]['merged_segment'] = merged_segment

    def fit_segment_to_targetregion(self, sampleid, segment_df):
        segment_gr = cnvmisc.arg_into_gr(segment_df)
        seg_subset = segment_gr.intersect(self.data[sampleid]['target_region'])
        target_diff_seg = self.data[sampleid]['target_region'].subtract(segment_gr)
        target_diff_seg = pyranges_helper.join(
            target_diff_seg, 
            segment_gr,
            how='left',
            merge=None,
            find_nearest=True,
            as_gr=True,
        )
        result = pr.concat([seg_subset, target_diff_seg]).df
        result = cnvmisc.sort_genome_df(result, self.refver)

        return result

    def deduplicate_merged_segment(self, merged_segment):
        merged_segment = cnvmisc.sort_genome_df(merged_segment, self.refver)

        chromdict = refgenome.get_chromdict(self.refver)
        annot_arr = np.concatenate(
            [
                merged_segment['Chromosome'].apply(chromdict.contigs.index).to_numpy()[:, np.newaxis],
                merged_segment.loc[:, ['depthratio_segment_mean', 'baf_segment_mean']].to_numpy(),
            ],
            axis=1,
        )
        values, counts, groupkey = tools.array_grouper(annot_arr, omit_values=True)
        groupby = merged_segment.groupby(groupkey)
        dedup_segdf = groupby.first()
        dedup_segdf['End'] = groupby['End'].last()

        return dedup_segdf

#    def postprocess_segment(self, sampleid, cellularity, ploidy):
#        # add CNt and B
#        cpinfo = self.data[sampleid]['cpscores'][(cellularity, ploidy)]
#        self.data[sampleid]['merged_segment']['CNt'] = cpinfo['CNt_list']
#        self.data[sampleid]['merged_segment']['B'] = cpinfo['B_list']
#        self.data[sampleid]['merged_segment']['A'] = (
#            self.data[sampleid]['merged_segment']['CNt']
#            - self.data[sampleid]['merged_segment']['B']
#        )
#
#        # add theoreticals
#        self.data[sampleid]['merged_segment'] = cnvmisc.add_theoreticals_to_segment(
#            self.data[sampleid]['merged_segment'], 
#            cellularity=cellularity, 
#            tumor_ploidy=ploidy, 
#            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'], 
#            is_female=self.data[sampleid]['is_female'],
#        )

    ####################
    # solution finding #
    ####################

    def make_CN_solution_freeccf(
        self,
        sampleid,
        cellularity,
        ploidy,
        depth_ratio_diff=None,
        baf_diff=None,
        min_N_CNt_candidates=5,
        N_CNt_candidates_fraction=0.5,
        limited_clonal=True,
    ):
        segdf = self.data[sampleid]['merged_segment']
        (
            clonal_solution, 
            flags, 
            freeccf_solution,
            calculated_depth_ratio, 
            calculated_baf,
            average_CNt,
        ) = cnvmisc.find_solution_before_ccfs(
            depth_ratio=segdf['depthratio_segment_mean'], 
            baf=segdf['corrected_baf_segment_mean'],
            CNn=segdf['CNn'],
            lengths=(segdf['End'] - segdf['Start']),
            cellularity=cellularity,
            tumor_ploidy=ploidy,
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
            depth_ratio_diff=depth_ratio_diff,
            baf_diff=baf_diff,
            Bn=1,
            min_N_CNt_candidates=min_N_CNt_candidates,
            N_CNt_candidates_fraction=N_CNt_candidates_fraction,

            limited_clonal=limited_clonal,
        )

        freeccf_result = {
            #'fixed_ccfs': fixed_ccfs, 
            #'ccf_plotdata': ccf_plotdata, 
            'clonal_solution': clonal_solution, 
            'flags': flags, 
            'freeccf_solution': freeccf_solution,
            'calculated_depth_ratio': calculated_depth_ratio, 
            'calculated_baf': calculated_baf,
        }
        self.data[sampleid]['freeccf_result'] = freeccf_result
        self.data[sampleid]['average_CNt'] = average_CNt

    def select_fixed_ccfs(self, sampleid, bandwidth=0.1):
        segdf = self.data[sampleid]['merged_segment']
        lengths = (segdf['End'] - segdf['Start']).to_numpy()
        fixed_ccfs, ccf_plotdata = cnvmisc.select_fixed_ccfs(
            freeccf_solution=self.data[sampleid]['freeccf_result']['freeccf_solution'], 
            lengths=lengths, 
            flags=self.data[sampleid]['freeccf_result']['flags'], 
            bandwidth=bandwidth,
        )
        self.data[sampleid]['fixed_ccfs'] = fixed_ccfs
        self.data[sampleid]['ccf_plotdata'] = ccf_plotdata

    def make_CN_solution_after_ccfs(
        self,
        sampleid,
        cellularity,
        ploidy,
        min_N_CNt_candidates=5,
        N_CNt_candidates_fraction=0.5,
        CNt_diff_factor=0.1,
        limited_clonal=True,
    ):
        segdf = self.data[sampleid]['merged_segment']
        solution = cnvmisc.find_solution_after_ccfs(
            depth_ratio=segdf['depthratio_segment_mean'].to_numpy(),
            baf=segdf['corrected_baf_segment_mean'].to_numpy(),
            CNn=segdf['CNn'].to_numpy(),
            lengths=(segdf['End'] - segdf['Start']).to_numpy(),
            cellularity=cellularity,
            tumor_ploidy=ploidy,
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
            average_CNt=self.data[sampleid]['average_CNt'],

            fixed_ccfs=self.data[sampleid]['fixed_ccfs'], 
            clonal_solution=self.data[sampleid]['freeccf_result']['clonal_solution'], 
            flags=self.data[sampleid]['freeccf_result']['flags'],

            Bn=1,
            min_N_CNt_candidates=min_N_CNt_candidates,
            N_CNt_candidates_fraction=N_CNt_candidates_fraction,

            CNt_diff_factor=CNt_diff_factor,

            limited_clonal=limited_clonal,
        )

        self.data[sampleid]['CN_solution'] = solution

#    def make_CN_solution_onestep(
#        self,
#        sampleid,
#        cellularity,
#        ploidy,
#        depth_ratio_diff=None,
#        baf_diff=0.05,
#        min_N_CNt_candidates=5,
#        N_CNt_candidates_fraction=0.5,
#        ccf_bw=0.1,
#    ):
#        self.make_CN_solution_freeccf(
#            sampleid,
#            cellularity,
#            ploidy,
#            depth_ratio_diff=depth_ratio_diff,
#            baf_diff=baf_diff,
#            min_N_CNt_candidates=min_N_CNt_candidates,
#            N_CNt_candidates_fraction=N_CNt_candidates_fraction,
#        )
#        self.select_fixed_ccfs(sampleid, bandwidth=ccf_bw)
#        self.make_CN_solution_after_ccfs(
#            sampleid,
#            cellularity,
#            ploidy,
#            min_N_CNt_candidates=min_N_CNt_candidates,
#            N_CNt_candidates_fraction=N_CNt_candidates_fraction,
#        )

    def add_freeccf_solution_to_segment(self, sampleid):
        #subclonal_theoretical_depthratio = self.data[sampleid]['freeccf_result']['calculated_depth_ratio']
        #subclonal_theoretical_baf = self.data[sampleid]['freeccf_result']['calculated_baf']
        cnvmisc.add_freeccf_solution_to_segment(
            segment_df=self.data[sampleid]['merged_segment'], 
            clonal_solution=self.data[sampleid]['freeccf_result']['clonal_solution'],
            freeccf_subclonal_solution=self.data[sampleid]['freeccf_result']['freeccf_solution'], 
            flags=self.data[sampleid]['freeccf_result']['flags'], 
            unfit_region_depth_ratio=self.data[sampleid]['freeccf_result']['calculated_depth_ratio'], 
            unfit_region_baf=self.data[sampleid]['freeccf_result']['calculated_baf'], 
        )

    def add_fixedccf_solution_to_segment(self, sampleid):
        cnvmisc.add_fixedccf_solution_to_segment(
            segment_df=self.data[sampleid]['merged_segment'], 
            fixedccf_solution=self.data[sampleid]['CN_solution'],
        )

    def add_solution_to_plotdata(self, sampleid):
        assert (
            self.data[sampleid]['merged_segment'].loc[:, ['Chromosome', 'Start', 'End']]
            == self.data[sampleid]['merged_segment_plotdata'].loc[:, ['Chromosome', 'Start', 'End']]
        ).all(axis=None)

        cnvmisc.add_fixedccf_solution_to_segment(
            segment_df=self.data[sampleid]['merged_segment_plotdata'], 
            fixedccf_solution=self.data[sampleid]['CN_solution'],
        )

    def add_CN_pred_to_segment(self, sampleid, cellularity, ploidy):
        if 'cpscores' not in self.data[sampleid]:
            self.data[sampleid]['cpscores'] = dict()

        # add CNt and B
        cpinfo = cnvmisc.calc_cp_score(
            segment_df=self.data[sampleid]['merged_segment'], 
            cellularity=cellularity, 
            ploidy=ploidy, 
            is_female=self.data[sampleid]['is_female'], 
            CNt_weight=1, 
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
        )
        self.data[sampleid]['cpscores'][(cellularity, ploidy)] = cpinfo

        self.data[sampleid]['merged_segment']['CNt'] = cpinfo['CNt_list']
        self.data[sampleid]['merged_segment']['B'] = cpinfo['B_list']
        self.data[sampleid]['merged_segment']['A'] = (
            self.data[sampleid]['merged_segment']['CNt']
            - self.data[sampleid]['merged_segment']['B']
        )

        # add theoreticals
        self.data[sampleid]['merged_segment'] = cnvmisc.add_theoreticals_to_segment(
            self.data[sampleid]['merged_segment'], 
            cellularity=cellularity, 
            tumor_ploidy=ploidy, 
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'], 
            is_female=self.data[sampleid]['is_female'],
        )

    def add_CN_pred_to_segment_new(self, sampleid, cellularity, ploidy):
        # calculate CNt values
        depth_ratios = self.data[sampleid]['merged_segment']['depthratio_segment_mean'].to_numpy()
        CNns = self.data[sampleid]['merged_segment']['CNn'].to_numpy()
        clonal_CNts = cnvmisc.calc_clonal_CNt(
            depth_ratio=depth_ratios,
            CNn=CNns,

            cellularity=cellularity,
            tumor_ploidy=ploidy,
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
        )

        subclonal_CNts, ccfs = cnvmisc.calc_subclonal_CNt(
            depth_ratio=depth_ratios,
            clonal_CNt=clonal_CNts,
            CNn=CNns,

            cellularity=cellularity,
            tumor_ploidy=ploidy,
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
            tumor_avg_depth_ratio=1,
            normal_avg_depth_ratio=1,
            only_max_ccf=True,
        )

        # calculate theoretical depths
        predicted_depth_ratios = cnvmisc.theoretical_depth_ratio_subclone(
            clonal_CNt=clonal_CNts, 
            subclonal_CNt=subclonal_CNts,
            ccf=ccfs,
            cellularity=cellularity, 
            tumor_ploidy=ploidy, 
            CNn=CNns, 
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
            tumor_avg_depth_ratio=1, 
            normal_avg_depth_ratio=1,
        )

        predicted_depth_ratios_clonal = cnvmisc.theoretical_depth_ratio(
            CNt=clonal_CNts, 
            cellularity=cellularity, 
            tumor_ploidy=ploidy, 
            CNn=CNns, 
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
            tumor_avg_depth_ratio=1, 
            normal_avg_depth_ratio=1,
        )

        # annotate segment dataframe
        merged_segdf = self.data[sampleid]['merged_segment']

        merged_segdf['clonal_CNt'] = clonal_CNts
        merged_segdf['subclonal_CNt'] = subclonal_CNts
        merged_segdf['ccf'] = ccfs
        merged_segdf['depthratio_predicted'] = predicted_depth_ratios
        merged_segdf['depthratio_predicted_clonal'] = predicted_depth_ratios_clonal

        self.data[sampleid]['merged_segment'] = merged_segdf

        # add standard deviation of depth ratios
        self.add_depthratio_std_to_segment(sampleid)

    def add_depthratio_std_to_segment(self, sampleid):
        self.data[sampleid]['merged_segment'] = add_depthratio_std(
            self.data[sampleid]['merged_segment'], 
            self.data[sampleid]['depthratio_upscaled'],
            self.refver,
        )

    def add_baf_std_to_segment(self, sampleid):
        merged_segdf = self.data[sampleid]['merged_segment']
        right_df = self.data[sampleid]['original_baf']
        right_df = right_df.loc[
            right_df['baf_raw_tumor'] > 0, 
            ['Chromosome', 'Start', 'End', 'baf_raw_tumor'],
        ]
        joined_segdf = pyranges_helper.join(
            merged_segdf,
            right_df,
            how='left', merge='mean', add_std=True, ddof=0,
            sort=True, refver=self.refver,
        )
        joined_segdf.drop('baf_raw_tumor', axis=1, inplace=True)
        joined_segdf.rename(
            columns={'baf_raw_tumor_std': 'baf_segment_std'}, 
            inplace=True,
        )

        self.data[sampleid]['merged_segment'] = joined_segdf

    ####################
    # cp value finding #
    ####################

    def get_valid_CNts(
        self, sampleid, depthratio1, depthratio2, CNt1_range=range(0, 6), CNt2_maxdiff=6,
    ):
        return cnvmisc.get_valid_CNts_from_depthratios(
            depthratio1, 
            depthratio2, 
            normal_mean_ploidy=self.data[sampleid]['normal_mean_ploidy'], 
            CNn=2, 
            tumor_avg_depth_ratio=1, 
            normal_avg_depth_ratio=1, 
            CNt1_range=CNt1_range, 
            CNt2_maxdiff=CNt2_maxdiff,
        )

    def calc_cpscore(
        self, 
        sampleid, 
        CNt_weight=libcncall.DEFAULT_CNT_WEIGHT,
        nproc=1,
    ):
        cpscore_dict = cnvmisc.get_cp_score_dict(
            self.data[sampleid]['merged_segment'], 
            refver=self.refver, 
            is_female=self.data[sampleid]['is_female'], 
            target_region_gr=self.data[sampleid]['target_region'], 
            CNt_weight=CNt_weight, 
            normal_ploidy=self.data[sampleid]['normal_mean_ploidy'],
            nproc=nproc,
        )
        peak_values, dfs = cnvmisc.get_peak_info(cpscore_dict)

        self.data[sampleid]['cpscores'] = cpscore_dict
        self.data[sampleid]['peak_values'] = peak_values
        self.data[sampleid]['peak_dfs'] = dfs

    def show_peaks(self, sampleid, figsize=(20, 20), **kwargs):
        fig, ax = cnvmisc.show_heatmap_peaks_new(
            self.data[sampleid]['peak_dfs'], figsize=figsize, **kwargs,
        )

    #########
    # depth #
    #########

    def load_bam_for_depth(self, sampleid, bam_path, sampletype):
        assert sampletype in ('normal', 'tumor')
        mosdepth_df = self._run_mosdepth(
            bam_path,
            binsize=self.default_binsize,
            region_df=(
                None 
                if self.data[sampleid]['mode'] == 'wgs' else 
                self.data[sampleid]['target_region']
            ),
        )

        self.data.setdefault(sampleid, dict())
        self.data[sampleid][f'{sampletype}_depth'] = mosdepth_df

    def load_mosdepth_df(
        self, sampleid, mosdepth_df, sampletype, use_default_gc=True, **kwargs,
    ):
        assert sampletype in ('normal', 'tumor')

        self.data.setdefault(sampleid, dict())
        depth_df, gcbin_average_depths = self._postprocess_mosdepth_df(
            mosdepth_df, 
            self.data[sampleid]['mode'],
            sampletype, 
            use_default_gc=use_default_gc, 
            **kwargs,
        )
        self.data[sampleid][f'{sampletype}_depth'] = depth_df
        self.data[sampleid][f'{sampletype}_gcdata'] = gcbin_average_depths

    def _run_mosdepth(self, bam_path, binsize, region_df):
        mosdepth_df, _ = libmosdepth.run_mosdepth(
            bam_path, 
            t=8, 
            use_median=False, 
            region_bed_path=None, 
            region_gr=region_df, 
            window_size=binsize, 
            donot_subset_bam=True,
            as_gr=False, 
            load_perbase=False,
        )

        return mosdepth_df

    def _postprocess_mosdepth_df(
        self, 
        mosdepth_df, 
        mode,
        sampletype, 
        use_default_gc=True, 
        **kwargs,
    ):
        assert sampletype in ('normal', 'tumor')

        # set preset_cutoffs
        #if mode == 'panel':
            #preset_cutoffs = 'panel'
        #elif mode == 'wgs':
            #if sampletype == 'normal':
                #preset_cutoffs = 'normal_wgs'
            #elif sampletype == 'tumor':
                #preset_cutoffs = 'wgs'
        #kwargs['preset_cutoffs'] = preset_cutoffs

        # set gc_df
        if use_default_gc:
            gc_df = libgcfraction.get_gc_df(
                self.refver, self.default_binsize, coords_as_index=True,
            )
            kwargs['gc_df'] = gc_df

        kwargs['as_gr'] = False
        depth_df, gcbin_average_depths = cnvmisc.postprocess_depth_df(
            mosdepth_df, 
            **kwargs,
        )
        return depth_df, gcbin_average_depths

    #######
    # baf #
    #######

    def load_germline_vcf(
        self, 
        sampleid, 
        vcf_path, 
        vcf_sampleid_tumor,
        vcf_sampleid_normal=None,
        nproc=1, 
        logging_lineno=50000,
    ):
        self.data.setdefault(sampleid, dict())

        if vcf_sampleid_normal is None:
            vcf_sampleids = [vcf_sampleid_tumor]
        else:
            vcf_sampleids = [vcf_sampleid_tumor, vcf_sampleid_normal]

        # load vafdf
        vaf_df = variantplus.get_vafdf(
            vcf_path, 
            sampleid=vcf_sampleids, 
            nproc=nproc,
        )
        # rename vaf columns
        vaf_df.rename(
            columns={f'vaf_{vcf_sampleid_tumor}': 'vaf_raw_tumor'}, inplace=True,
        )
        if vcf_sampleid_normal is not None:
            vaf_df.rename(
                columns={f'vaf_{vcf_sampleid_normal}': 'vaf_raw_normal'}, inplace=True,
            )

        # postprocess
        #vaf_df = vaf_df.loc[vaf_df['vaf_raw'].notna().to_numpy(), :]
        vaf_df.reset_index(drop=True, inplace=True)
        vaf_df['baf_raw_tumor'] = cnvmisc.get_bafs(vaf_df['vaf_raw_tumor'])
        if 'vaf_raw_normal' in vaf_df:
            vaf_df['baf_raw_normal'] = cnvmisc.get_bafs(vaf_df['vaf_raw_normal'])

        selected_cols = [
            'Chromosome', 
            'Start', 
            'End', 
            'vaf_raw_tumor', 
            'baf_raw_tumor', 
        ]
        if vcf_sampleid_normal is not None:
            selected_cols.append('vaf_raw_normal')
            selected_cols.append('baf_raw_normal')
        vaf_df = vaf_df.loc[:, selected_cols]
        self.data[sampleid]['original_baf'] = vaf_df

    #################
    # other helpers #
    #################

    def set_normal_mean_ploidy(self, sampleid):
        self.data[sampleid]['normal_mean_ploidy'] = cnvmisc.get_normal_mean_ploidy(
            self.refver, 
            self.data[sampleid]['is_female'], 
            self.data[sampleid]['target_region'],
        )

    @staticmethod
    def add_CNn_to_segment(
        segment_df,
        mode,
        refver,
        is_female,
        target_region=None,
    ):
        if mode == 'wgs':
            segment_df = rcopynumber.add_CNn_to_wgs_segment_gr(
                segment_df, refver, is_female, as_gr=False,
            )
        elif mode == 'panel':
            assert target_region is not None
            segment_df = rcopynumber.add_CNn_to_targetseq_segment_gr(
                segment_df, target_region, refver, is_female, as_gr=False,
            )
        return segment_df


# segmentation functions for running with multiprocessing

def handle_gamma_kmin_args(gamma, kmin, mode):
    if gamma is None:
        if mode == 'wgs':
            gamma = 40
        elif mode == 'panel':
            gamma = 30
    if kmin is None:
        if mode == 'wgs':
            kmin = 5
        elif mode == 'panel':
            kmin = 1

    return gamma, kmin


def _make_depth_segment(
    depthratio_df,
    mode,
    refver,
    winsorize=False,
    gamma=None,
    kmin=None,
    verbose=False,
):
    gamma, kmin = handle_gamma_kmin_args(gamma, kmin, mode)

    #if 'depthratio_raw' in depthratio_df.columns:
    #    input_df = depthratio_df.rename(columns={'depthratio_raw': 'depth_raw'})
    #else:
    #    input_df = depthratio_df
    input_df = depthratio_df.rename(columns={'depthratio_raw': 'depth_raw'})
    segment_df, _ = rcopynumber.run_rcopynumber_unified(
        depth_df=input_df,
        refver=refver,

        as_gr=False, 
        winsorize=winsorize,
        #compact=(mode == 'panel'), 
        compact=False,
        verbose=verbose,
        gamma=gamma,
        kmin=kmin,
    )

    segment_df.rename(columns={'depth_segment_mean': 'depthratio_segment_mean'}, inplace=True)

    segment_df = add_depthratio_std(segment_df, depthratio_df, refver)

    return segment_df


def add_depthratio_std(segment_df, raw_df, refver):
    # make left df
    left_df = segment_df.drop(
        segment_df.columns.intersection(
            ['depthratio_segment_std', 'depthratio_segment_std_mean_ratio'], 
        ),
        axis=1, 
        inplace=False,
    )
    # make right df
    right_df = raw_df.loc[
        :, ['Chromosome', 'Start', 'End', 'depthratio_raw']
    ].copy()
    # join
    joined_segdf = pyranges_helper.join(
        left_df,
        right_df,
        how='left', merge='mean', add_std=True, ddof=0,
        sort=True, refver=refver,
    )
    joined_segdf.drop('depthratio_raw', axis=1, inplace=True)
    joined_segdf.rename(
        columns={'depthratio_raw_std': 'depthratio_segment_std'}, 
        inplace=True,
    )

    # add std/mean ratio
    joined_segdf['depthratio_segment_std_mean_ratio'] = (
        joined_segdf['depthratio_segment_std']
        / joined_segdf['depthratio_segment_mean']
    ).to_numpy()

    return joined_segdf


def _make_baf_segment(
    baf_df,
    target_region,
    mode,
    refver,
    winsorize=False,
    gamma=None,
    kmin=None,
    verbose=False,
    baf_cutoff=0.1,
    #bafcorrector=bafcorrection.load_bafcorrect_func(),
):
    gamma, kmin = handle_gamma_kmin_args(gamma, kmin, mode)

    targetovlp_baf_df = pr.PyRanges(baf_df).overlap(target_region).df
    input_df = targetovlp_baf_df.loc[:, ['Chromosome', 'Start', 'End']]
    input_df['depth_raw'] = targetovlp_baf_df['baf_raw_tumor']
    input_df = input_df.loc[input_df['depth_raw'] > baf_cutoff, :]

    segment_df, _ = rcopynumber.run_rcopynumber_unified(
        depth_df=input_df,
        refver=refver,

        as_gr=False, 
        winsorize=winsorize,
        #compact=(mode == 'panel'), 
        compact=False,
        verbose=verbose,
        gamma=gamma,
        kmin=kmin,
    )

    segment_df.rename(columns={'depth_segment_mean': 'baf_segment_mean'}, inplace=True)

    return segment_df


def _make_segments_targetfunc_depth(shdict, shdict_key, **kwargs):
    shdict[shdict_key] = _make_depth_segment(**kwargs)


def _make_segments_targetfunc_baf(shdict, shdict_key, **kwargs):
    shdict[shdict_key] = _make_baf_segment(**kwargs)


def _make_segments_main(
    depthratio_df,
    mode,
    refver,
    winsorize,
    depthratio_gamma,
    depthratio_kmin,
    baf_gamma,
    baf_kmin,
    verbose,

    baf_df,
    target_region,
    baf_cutoff,
    #bafcorrector,
):
    with multiprocessing.Manager() as manager:
        shdict = manager.dict()
        subp1 = multiprocessing.Process(
            target=_make_segments_targetfunc_depth,
            args=(
                shdict,
                'depth_segment',
            ),
            kwargs=dict(
                depthratio_df=depthratio_df,
                mode=mode,
                refver=refver,
                winsorize=winsorize,
                gamma=depthratio_gamma,
                kmin=depthratio_kmin,
                verbose=verbose,
            ),
        )
        subp2 = multiprocessing.Process(
            target=_make_segments_targetfunc_baf,
            args=(
                shdict, 
                'baf_segment',
            ),
            kwargs=dict(
                baf_df=baf_df,
                target_region=target_region,
                mode=mode,
                refver=refver,
                winsorize=winsorize,
                gamma=baf_gamma,
                kmin=baf_kmin,
                verbose=verbose,
                baf_cutoff=baf_cutoff,
            ),
        )

        subp1.start()
        subp2.start()

        subp1.join()
        subp2.join()

        # result
        depth_segment = shdict['depth_segment']
        baf_segment = shdict['baf_segment']

    return depth_segment, baf_segment


###############################################################


def make_targetseq_cnvplot(data_df, draw_invalid_regions=False):
    # sanity check
    required_cols = {
        'CNt', 'B', 
        'baf_segment_mean', 'baf_predicted', 'baf_raw',
        'depthratio_segment_mean', 'depthratio_predicted', 'depthratio_raw',
    }
    assert required_cols.issubset(data_df.columns)

    # setup
    if isinstance(data_df, pr.PyRanges):
        data_df = data_df.df

    gplotter = GenomePlotter(data_df)
    fig, axd = plt.subplot_mosaic(
        [
            ['CN',], 
            ['baf',], 
            ['depth',],
        ],
        figsize=(30, 13),
    )
    for ax in axd.values():
        gplotter.set_xlim(ax)

    # CN
    axd['CN'].set_ylabel('CN')

    gplotter.draw_hlines(
        axd['CN'], df=data_df, y_colname='CNt', offset=0.1,
        plot_kwargs={'color': 'black'},
    )
    gplotter.draw_hlines(
        axd['CN'], df=data_df, y_colname='B', 
        plot_kwargs={'color': 'blue'},
    )
    gplotter.fit_spines_to_regions(axd['CN'])

    # baf
    axd['baf'].set_ylabel('baf')
    axd['baf'].set_ylim(0, 0.6)

    gplotter.draw_dots(
        axd['baf'], df=data_df, y_colname='baf_raw', 
        plot_kwargs={'color': 'gray', 'markersize': 1.5, 'alpha': 0.5},
    )
    gplotter.draw_hlines(
        axd['baf'], df=data_df, y_colname='baf_segment_mean', 
        plot_kwargs={'color': 'black', 'linewidth': 2, 'alpha': 0.3},
    )
    gplotter.draw_hlines(
        axd['baf'], df=data_df, y_colname='baf_predicted', 
        plot_kwargs={'color': 'green', 'linewidth': 2, 'alpha': 0.3},
    )
    gplotter.fit_spines_to_regions(axd['baf'])

    # depthratio
    axd['depth'].set_ylabel('depth ratio')
    axd['depth'].set_ylim(
        0, 
        data_df['depthratio_raw'].max() * 1.1
    )

    gplotter.draw_dots(
        axd['depth'], df=data_df, y_colname='depthratio_raw', 
        plot_kwargs={'color': 'gray', 'markersize': 1.5, 'alpha': 0.5},
    )
    gplotter.draw_hlines(
        axd['depth'], df=data_df, y_colname='depthratio_segment_mean', 
        plot_kwargs={'color': 'black', 'linewidth': 2, 'alpha': 0.5},
    )
    gplotter.draw_hlines(
        axd['depth'], df=data_df, y_colname='depthratio_predicted', 
        plot_kwargs={'color': 'green', 'linewidth': 2, 'alpha': 0.5},
    )

    if draw_invalid_regions:
        selector = np.logical_or(
            data_df['depthratio_raw'].isna().to_numpy,
            data_df['baf_raw'].isna().to_numpy,
        )
        gplotter.draw_bgcolors(
            axd['depth'], df=data_df.loc[selector, :], 
            plot_kwargs=dict(color='red', alpha=0.01)
        )

    gplotter.fit_spines_to_regions(axd['depth'])

    return fig, axd


def draw_targetseq_cnvplot_from_data_precp(
    tumor_depth_df, 
    normal_depth_df,
    germline_df,

    region_gr,
    refver, 
    is_female, 
):
    assert isinstance(tumor_depth_df, pd.DataFrame)
    assert isinstance(normal_depth_df, pd.DataFrame)
    assert isinstance(germline_df, pd.DataFrame)

    if 'baf_raw' not in germline_df.columns:
        germline_df['baf_raw'] = cnvmisc.get_bafs(germline_df['vaf'])
    germline_gr = pr.PyRanges(germline_df)

    # make depth ratio
    depthratio_gr = cnvmisc.make_depth_ratio(tumor_depth_df, normal_depth_df, as_gr=True)
    depthratio_df = depthratio_gr.df

    # run R copynumber
    segment_gr, depth_baf_gr = rcopynumber.run_rcopynumber_unified(
        depthratio_df.rename(columns={'depth_ratio_sequenzastyle': 'depth_raw'}),
        baf_df=germline_df, 
        refver=refver, 
        as_gr=True, 
        compact=True, 
        winsorize=False,
        gamma=30, 
        kmin=1,
    )
    segment_gr = segment_gr[[]]

    # postprocess segment df
    segment_gr = rcopynumber.add_CNn_to_targetseq_segment_gr(segment_gr, region_gr, refver=refver, is_female=is_female)
    segment_gr = pyranges_helper.join(
        segment_gr, 
        depthratio_gr[['depth_ratio_sequenzastyle']], 
        how='left', merge='mean', find_nearest=False, as_gr=True,
    )
    segment_gr = pyranges_helper.join(
        segment_gr, 
        germline_gr[['baf_raw']], 
        how='left', merge='mean', find_nearest=False, as_gr=True,
    )
    segment_df = segment_gr.df
    segment_df.rename(
        {'depth_ratio_sequenzastyle': 'depthratio_segment_mean', 'baf_raw': 'baf_segment_mean'},
        inplace=True,
        axis=1,
    )

    return segment_df, depth_baf_gr


def draw_targetseq_cnvplot_from_data_choosecp(
    segment_df, 
    depth_baf_gr,

    region_gr,
    refver,
    is_female,

    CNt_weight=libcncall.DEFAULT_CNT_WEIGHT,
    nproc=None,
):
    cpscore_dict = cnvmisc.get_cp_score_dict(segment_df, refver=refver, is_female=is_female, target_region_gr=region_gr, CNt_weight=CNt_weight, nproc=nproc)
    peak_values, dfs = cnvmisc.get_peak_info(cpscore_dict)

    cnvmisc.show_heatmap_peaks_new(dfs, figsize=(20, 20))

    return peak_values, dfs, cpscore_dict


def draw_targetseq_cnvplot_from_data_postcp(
    cellularity,
    ploidy,
    segment_df, 
    depth_baf_gr,

    region_gr,
    refver,
    is_female,

    draw_invalid_regions=False,
):
    normal_mean_ploidy = cnvmisc.get_normal_mean_ploidy(refver, is_female, target_region_gr=region_gr)

    # postprocess segment df
    segment_df = cnvmisc.add_CNt_to_segment(
        segment_df, 
        cellularity=cellularity, 
        tumor_ploidy=ploidy, 
        is_female=is_female, 
        normal_ploidy=normal_mean_ploidy,
    )
    segment_df = cnvmisc.add_theoreticals_to_segment(
        segment_df, 
        cellularity=cellularity, 
        tumor_ploidy=ploidy, 
        normal_ploidy=normal_mean_ploidy, 
        is_female=is_female,
    )

    # make a single merged plot data df
    depth_baf_gr.depthratio_raw = depth_baf_gr.depth_raw
    depth_baf_gr = cnvmisc.annotate_region_with_segment(depth_baf_gr, pr.PyRanges(segment_df), as_gr=True)

    # draw plot
    fig, axd = make_targetseq_cnvplot(depth_baf_gr, draw_invalid_regions=draw_invalid_regions)

    return fig, axd



import re
import sys
import os
import inspect
import itertools
import functools
import operator
import collections
import json
import pickle
import shutil
import warnings

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
from handygenome.peakutils import DensityPeaks
import handygenome.logutils as logutils
import handygenome.refgenome.refgenome as refgenome
import handygenome.refgenome.refverfile as refverfile
from handygenome.refgenome.refverfile import PARUnavailableError

import handygenome.genomedf.genomedf as libgdf
from handygenome.genomedf.genomedf import GenomeDataFrame as GDF
import handygenome.cnv.depth as libdepth
from handygenome.variant.vcfdataframe import VariantDataFrame

from handygenome.cnv.depth import DepthRawDataFrame as DepthRawDF
from handygenome.cnv.depth import DepthSegmentDataFrame as DepthSegDF

import handygenome.cnv.baf as libbaf
from handygenome.cnv.baf import BAFRawDataFrame as BAFRawDF
from handygenome.cnv.baf import BAFSegmentDataFrame as BAFSegDF
import handygenome.cnv.bafsimul as bafsimul

from handygenome.cnv.cnvsegment import CNVSegmentDataFrame as CNVSegDF

import handygenome.cnv.cnvcall as cnvcall
from handygenome.cnv.cnvcall import KCError
import handygenome.cnv.mosdepth as libmosdepth

import handygenome.plot.genomeplot as libgenomeplot
import handygenome.plot.misc as plotmisc
from handygenome.plot.genomeplot import GenomePlotter
import handygenome.genomedf.genomedf_draw as genomedf_draw
from handygenome.genomedf.genomedf_draw import GenomeDrawingFigureResult


##############
# decorators #
##############

def deco_nproc(func):
    sig = inspect.signature(func)

    req_params = set(['nproc', 'self'])
    if not set(req_params).issubset(sig.parameters.keys()):
        raise Exception(
            f'Decorated plotter method does not have required parameters: {req_params}'
        )

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        ba = sig.bind(*args, **kwargs)
        ba.apply_defaults()
        if ba.arguments['nproc'] == 0:
            ba.arguments['nproc'] = ba.arguments['self'].nproc

        return func(*ba.args, **ba.kwargs)

    return wrapper


def get_deco_logging(msg):
    def decorator(func):
        sig = inspect.signature(func)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()

            if hasattr(ba.arguments['self'], 'log'):
                logfunc = ba.arguments['self'].log
            else:
                logfunc = logutils.log

            logfunc('Beginning ' + msg)
            result = func(*args, **kwargs)
            logfunc('Finished ' + msg)

            return result

        return wrapper

    return decorator


def deco_corrector_id(func):
    sig = inspect.signature(func)
    req_params = set(['corrector_id'])
    if not set(req_params).issubset(sig.parameters.keys()):
        raise Exception(
            f'Decorated plotter method does not have required parameters: {req_params}'
        )

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        ba = sig.bind(*args, **kwargs)
        ba.apply_defaults()
        if ba.arguments['corrector_id'] is None:
            draw_result_corrector_id = ba.arguments['self'].draw_result_corrector_id
            if draw_result_corrector_id is None:
                raise Exception(f'"draw_result_corrector_id" is not set.')
            ba.arguments['corrector_id'] = draw_result_corrector_id

        return func(*ba.args, **ba.kwargs)

    return wrapper


##################
# common methods #
##################

def getattr_base(obj, key):
    if key in obj.simple_attrs.keys():
        return obj.simple_attrs[key]
    else:
        return super(obj.__class__, obj).__getattribute__(key)


#############
# CNVSample #
#############

class CNVSample:
    default_window = {'panel': 100, 'wgs': 1000}
    save_paths_basenames = {
        'simple_attrs': 'simple_attrs.json',
        'solution_attrs': 'solution_attrs.pickle',

        # target region
        'raw_target_region': 'raw_target_region.tsv.gz',
        'target_region': 'target_region.tsv.gz',
        'excluded_region': 'excluded_region.tsv.gz',

        # raw data
        'depth_rawdata': 'depth_rawdata.tsv.gz',
        'baf_rawdata': 'baf_rawdata.tsv.gz',
        #'baf_hetalt_rawdata': 'baf_hetalt_rawdata.tsv.gz',

        # segment
        'depth_segment': 'depth_segment.tsv.gz',

        'baf_segment': 'baf_segment.tsv.gz',
        'baf_noedit_segment': 'baf_noedit_segment.tsv.gz',
        'baf_edited_segment': 'baf_edited_segment.tsv.gz',
        'baf_rmzero_noedit_segment': 'baf_rmzero_noedit_segment.tsv.gz',

        'merged_segment': 'merged_segment.tsv.gz',
    }
    gdf_types = {
        'raw_target_region': GDF,
        'target_region': GDF,
        'excluded_region': GDF,

        'depth_rawdata': DepthRawDF,
        'baf_rawdata': BAFRawDF,
        'filtered_baf_rawdata': BAFRawDF,
        #'baf_hetalt_rawdata': BAFRawDF,

        'depth_segment': DepthSegDF,

        'baf_segment': BAFSegDF,
        'baf_noedit_segment': BAFSegDF,
        'baf_edited_segment': BAFSegDF,
        'baf_rmzero_noedit_segment': BAFSegDF,

        'merged_segment': CNVSegDF,
    }
    corrected_gdf_types = {
        'depth_rawdata': DepthRawDF,
        'depth_segment': DepthSegDF,
        'depth_excluded_region': GDF,

        'baf_rawdata': BAFRawDF,
        'filtered_baf_rawdata': BAFRawDF,

        'merged_segment': CNVSegDF,
    }
    corrected_gdf_types_bafidx = {
        'baf_segment': BAFSegDF,
    }

    def __repr__(self):
        return f'{self.__class__.__name__} object (sampleid: {self.simple_attrs["sampleid"]})'

    def __getattr__(self, key):
        return getattr_base(self, key)

#    def __getattr__(self, key):
#        if key in self.simple_attrs.keys():
#            return self.simple_attrs[key]
#        else:
#            try:
#                super().__getattr__(key)
#            except AttributeError:
#                raise AttributeError(f'{repr(self.__class__)} object has no attribute {repr(key)}')

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
        is_germline,
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
        nproc=min(10, os.cpu_count()),
        genomeplotter_kwargs=None,
        verbose=True,

        # others
        #norm_method='plain',
    ):
        bafinfo_given = cls.init_sanitycheck(
            mode, 
            target_region_path,
            target_region_gdf,
            vcf_path,
            vcf_sampleid,
            bafdf,
        )

        result = cls()

        # set simple attributes
        result.init_simple_attrs(
            bam_path=bam_path,
            vcf_path=vcf_path,
            vcf_sampleid=vcf_sampleid,
            sampleid=sampleid,
            refver=refgenome.standardize(refver),
            is_female=is_female,
            is_germline=is_germline,
            mode=mode,
            nproc=nproc,
            verbose=verbose,
            ploidy=ploidy,
            genomeplotter_kwargs=None,
        )

        # initiate gdf data dict
        result.init_gdfs()

        # initiate cache
        result.init_solution_attrs()
        result.init_cache()

        # target region
        result.load_raw_target_region(target_region_path, target_region_gdf)

        # genomeplotter
        result.set_genomeplotter_kwargs(genomeplotter_kwargs)
        result.init_genomeplotter(**result.simple_attrs['genomeplotter_kwargs'])

        # raw data
        result.load_depth(bam_path)
        if bafinfo_given:
            result.load_baf(
                vcf_path=vcf_path, 
                vcf_sampleid=vcf_sampleid, 
                bafdf=bafdf,
            )

        # plotdata
        result.plotdata_cache = dict()

        return result

    def init_simple_attrs(
        self,
        bam_path=None,
        vcf_path=None,
        vcf_sampleid=None,
        sampleid=None,
        refver=refgenome.standardize('hg19'),
        is_female=None,
        is_germline=None,
        mode=None,
        nproc=10,
        verbose=True,
        ploidy=2,
        genomeplotter_kwargs=dict(),
    ):
        self.simple_attrs = {
            'bam_path': bam_path,
            'vcf_path': vcf_path,
            'vcf_sampleid': vcf_sampleid,
            'sampleid': sampleid,
            'refver': refver,
            'is_female': is_female,
            'is_germline': is_germline,
            'mode': mode,
            'nproc': nproc,
            'verbose': verbose,
            'ploidy': ploidy,
            'genomeplotter_kwargs': genomeplotter_kwargs,
        }

    def init_gdfs(self):
        self.gdfs = {
            'raw_target_region': None,
            'target_region': None,
            'excluded_region': None,

            'depth_rawdata': None,
            'baf_rawdata': None,
            'filtered_baf_rawdata': None,

            'depth_segment': None,
            'baf_segment': {baf_idx: None for baf_idx in self.get_baf_indexes()},
            'baf_noedit_segment': {baf_idx: None for baf_idx in self.get_baf_indexes()},
            'baf_edited_segment': {baf_idx: None for baf_idx in self.get_baf_indexes()},
            'baf_rmzero_noedit_segment': {baf_idx: None for baf_idx in self.get_baf_indexes()},

            'merged_segment': None,
            'merged_merged_segment': None,
        }
        self.corrected_gdfs = dict()

    def init_solution_attrs(self):
        self.solution_attrs = dict()

    def init_solution_attrs_item(self, corrector_id):
        self.solution_attrs[corrector_id] = {
            'solution_candidates': None,
            'selected_solution': None,
            'selected_regions': None,
        }

    def init_cache(self):
        pass
        #self.draw_result = None
        #self.draw_result_corrector_id = None
        #self.mean_values = None

    def set_genomeplotter_kwargs(self, genomeplotter_kwargs):
        if genomeplotter_kwargs is None:
            if self.mode == 'wgs':
                genomeplotter_kwargs = dict()
            elif self.mode == 'panel':
                region_gdf = self.gdfs['raw_target_region'].copy().drop_annots()
                genomeplotter_kwargs = dict(
                    chroms=list(region_gdf.chroms),
                    start0s=list(int(x) for x in region_gdf.start0s),
                    end0s=list(int(x) for x in region_gdf.end0s),
                    weights=1,
                    region_gaps=0,
                )
        self.simple_attrs['genomeplotter_kwargs'] = genomeplotter_kwargs

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
        bafinfo_as_vcf = (vcf_path is not None) and (vcf_sampleid is not None)
        bafinfo_as_bafdf = (bafdf is not None)
        bafinfo_given = (bafinfo_as_vcf or bafinfo_as_bafdf)
        if bafinfo_as_vcf and bafinfo_as_bafdf:
            raise Exception(f'BAF loading argument usage: 1) "vcf_path" and "vcf_sampleid" ; 2) "bafdf"')
        if (mode == 'wgs') and (not bafinfo_given):
            raise Exception(f'BAF data must be given with a WGS sample.')

        return bafinfo_given

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
    def load_depth(self, bam_path, window=None, t=1, readlength=151, nproc=0):
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

        # postprocess
        self.get_depth_rawdata().set_normalized_depth()
        self.add_MQ_to_depth_rawdata(readlength=readlength, nproc=nproc)

        self.log(f'Finished making depth rawdata from bam file')

    @deco_nproc
    def load_baf(self, vcf_path=None, vcf_sampleid=None, bafdf=None, nproc=0):
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

        bafdf.remove_overlapping_rows()
        if self.mode == 'panel':
            bafdf = bafdf.intersect(self.get_raw_target_region())
            bafdf.sort()

        self.gdfs['baf_rawdata'] = bafdf

    def check_bafinfo_exists(self):
        return self.gdfs['baf_rawdata'] is not None

    def check_nobias(self):
        return (self.mode == 'wgs')

    def check_germline_solution_available(self):
        return self.check_nobias() and self.is_germline

    ##############
    # mainstream #
    ##############

    @deco_nproc
    def mainstream_processing(
        self, 
        nproc=0, 

        # depth filter
        MQ_cutoff_wgs=(45, None),
        norm_depth_cutoff_wgs=(0, None),
        MQ_cutoff_panel=(45, None),
        norm_depth_cutoff_panel=None,
        raw_depth_cutoff_panel=(50, None),

        # baf seg modification
        ndata_cutoff=30, 
        corrected_baf_cutoff=0.41,

        # skipping
        skip_baf_noedit_segment=False,
        skip_depth_segment=False,
    ):
        """For germline sample"""
        assert self.is_germline

        if self.check_nobias():
            assert self.check_bafinfo_exists()

        # depth
        if self.check_nobias():  # depth segmentation can be done only in this case
            if not skip_depth_segment:
                self.make_depth_segment(nproc=nproc)
                #self.add_rawdata_to_depth_segment(rawdepth=False, nproc=nproc)
            self.add_filter_to_depth_segment(
                MQ_cutoff=MQ_cutoff_wgs, 
                norm_depth_cutoff=norm_depth_cutoff_wgs,
            )
        else:
            self.get_depth_rawdata().add_filter(
                MQ_cutoff=MQ_cutoff_panel, 
                norm_depth_cutoff=norm_depth_cutoff_panel, 
                raw_depth_cutoff=raw_depth_cutoff_panel,
            )
#                self.add_filter_to_depth_rawdata(
#                    MQ_cutoff=MQ_cutoff_panel, 
#                    norm_depth_cutoff=norm_depth_cutoff_panel,
#                    raw_depth_cutoff=raw_depth_cutoff_panel,
#                )

        # baf
        if self.check_bafinfo_exists():
            self.make_baf_segment(nproc=nproc)

#            if not skip_baf_noedit_segment:
#                self.make_baf_noedit_segment(nproc=nproc)
#            self.make_baf_edited_segment(
#                nproc=nproc,
#                ndata_cutoff=ndata_cutoff, 
#                merge_low_ndata=None,
#                corrected_baf_cutoff=corrected_baf_cutoff,
#            )

        # cnv segment and germline solution
        if self.check_nobias():
            self.make_merged_segment(nproc=nproc)

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

    #########################
    # CN calling (solution) #
    #########################

    #def add_clonal_solution_to_germline_sample(self):
    #    """Only for non-biased germline sample"""
    #    assert self.is_germline
    #    assert self.check_nobias()
#
#        cnv_seg_gdf = self.get_merged_segment()
#        cnv_seg_gdf.add_clonal_solution_germline(self.is_female, self.ploidy)

    def add_clonal_solution_to_corrected_sample(self, *, corrector_id, cellularity, K):
        """merged segment for the corrector id must be annotated with germline CN
        with 'add_gCN_to_corrected_merged_segment' method of CNVManager object
        """
        merged_segment = self.get_merged_segment(corrector_id=corrector_id)
        merged_segment.add_clonal_solution_targetsample(cellularity=cellularity, K=K)

        # save cellularity, K, mean ploidy
        corrected_gdfs = self.corrected_gdfs[corrector_id]
        corrected_gdfs['cellularity'] = cellularity
        corrected_gdfs['K'] = K
        corrected_gdfs['ploidy'] = merged_segment.get_average_ploidy()

    def test_solution_fitnesses(self, corrector_id, *, lowest_depth, delta_depth, CNg, ncand=10):
        """Arguments:
            lowest_depth: Must be the depth value of the lowest depth cluster
            ncand: The number of candidate CNt at the level corresponding to 'lowest_depth' value (CNt candidates become 0, 1, 2, ..., ncand)
        """
        all_data = list()

        # calculate fitness for each CNt candidate
        merged_segment = self.get_merged_segment(corrector_id=corrector_id)
        for CNt in range(ncand):
            try:
                K, cellularity = cnvcall.find_Kc_withdelta(
                    depth=lowest_depth, 
                    CNt=CNt, 
                    CNg=CNg, 
                    delta_depth=delta_depth,
                )
            except KCError:
                continue
            else:
                merged_segment.add_clonal_solution_targetsample(cellularity=cellularity, K=K)
                fitness = merged_segment.get_Bt_fitness()
                all_data.append({'CNt': CNt, 'K': K, 'cellularity': cellularity, 'fitness': fitness})

        merged_segment.drop_solution_columns()

        # find best solution
        best_dict = min(all_data, key=(lambda x: x['fitness']))

        return best_dict, all_data
            
    def drop_solution(self, corrector_id):
        self.get_merged_segment(corrector_id=corrector_id).drop_solution_columns()

        corrected_gdfs = self.corrected_gdfs[corrector_id]
        corrected_gdfs['cellularity'] = None
        corrected_gdfs['K'] = None
        corrected_gdfs['ploidy'] = None

    #########################
    # target region setting #
    #########################

    def set_excluded_region_germline(self):
        """Must be run after running 'add_clonal_solution_to_germline_sample'"""

        excluded_region_list = list()

        # baf
        if self.check_bafinfo_exists():
            for baf_seg_gdf in self.get_baf_edited_segment_dict().values():
                bafseg_selector = baf_seg_gdf.get_filter()
                baf_excl_region_gdf = baf_seg_gdf.drop_annots().loc[~bafseg_selector, :]
                excluded_region_list.append(baf_excl_region_gdf)

        # depth
        if self.check_nobias():
            depth_seg_gdf = self.get_depth_segment() 
            selector = np.logical_not(depth_seg_gdf.get_filter())
            depth_excl_region_gdf = depth_seg_gdf.drop_annots().loc[selector, :]
        else:
            depth_rawdata = self.get_depth_rawdata()
            selector = np.logical_not(depth_rawdata.get_filter())
            depth_excl_region_gdf = depth_rawdata.drop_annots().loc[selector, :]

        if self.is_female:
            Y_chrom = self.get_chromdict().XY_names[1]
            depth_excl_region_gdf = depth_excl_region_gdf.subset_chroms(
                set(depth_excl_region_gdf.chroms).difference([Y_chrom])
            )

        excluded_region_list.append(depth_excl_region_gdf)

        # germline CN 0
        if self.check_nobias():
            merged_segment = self.get_merged_segment()
            selector = merged_segment.get_clonal_CN() == 0
            CN0_region_gdf = merged_segment.loc[selector, :].drop_annots()
        else:
            self.add_default_CNg_Bg_to_depth_rawdata()
            depth_rawdata = self.get_depth_rawdata()
            selector = depth_rawdata[cnvcall.DEFAULT_CNG_COLNAME] == 0
            CN0_region_gdf = depth_rawdata.loc[selector, :].drop_annots()
        excluded_region_list.append(CN0_region_gdf)

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
                chromdict = self.get_chromdict()
                y_chrom = chromdict.XY_names[1]
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
        elif self.mode == 'panel':
            allregion = self.get_raw_target_region()

        chromdict = self.get_chromdict()
        if (
            self.is_female
            and chromdict.check_has_XY()
        ):
            y_chrom = chromdict.XY_names[1]
            selected_chroms = set(allregion.chroms).difference([y_chrom])
            allregion = allregion.subset_chroms(selected_chroms)

        # subtract excluded region
        self.set_excluded_region_germline()
        target_region = allregion.subtract(self.get_excluded_region())
        target_region.sort()

        self.gdfs['target_region'] = target_region

    def find_onecopy_depth(self):
        return self.get_depth_segment().find_onecopy_depth(
            self.is_female, self.ploidy,
        )

#    def find_germline_intcopy_validbaf_region(self, factors=(0.3, 2)):
#        assert self.simple_attrs["mode"] == 'wgs'
#        assert self.get_depth_segment('normal') is not None
#        assert all(
#            (self.get_baf_segment('normal', baf_idx) is not None)
#            for baf_idx in self.get_baf_indexes()
#        )
#
#        depth_seg_gdf = self.get_depth_segment('normal').copy()
#        baf_seg_gdfs = self.get_baf_segment_dict('normal')
#
#        # onecopy depth
#        onecopy_depth = self.get_onecopy_depth(depth_seg_gdf)
#        depth_seg_gdf['onecopy_depth_ratio'] = (
#            depth_seg_gdf[DepthRawDF.norm_depth_colname]
#            / onecopy_depth
#        )
#


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

        self.log(f'Beginning saving as pickle')
        
        if os.path.exists(filepath):
            logutils.log(f'Overwriting an existing file', level='warning')

        with open(filepath, 'wb') as outfile:
            pickle.dump(self, outfile)

        self.log(f'Finished saving')

    @classmethod
    def load_pickle(cls, filepath):
        with open(filepath, 'rb') as infile:
            return pickle.load(infile)

    ###########################################
    # pickle save/load, each gdf as dataframe #
    ###########################################

    def save_nonbafidx_gdf_aspickle(self, gdf_key, topdir):
        gdf = self.gdfs[gdf_key]
        if gdf is not None:
            basename = self.__class__.save_paths_basenames[gdf_key]
            assert basename.endswith('.tsv.gz')
            pkl_basename = re.sub(r'\.tsv\.gz$', '.pickle', basename)

            savepath = os.path.join(topdir, pkl_basename)
            with open(savepath, 'wb') as outfile:
                pickle.dump(gdf.df, outfile)

    def save_bafidx_gdf_aspickle(self, gdf_key, topdir):
        if gdf_key in self.gdfs:
            for baf_idx in self.get_baf_indexes():
                gdf = self.gdfs[gdf_key][baf_idx]
                if gdf is not None:
                    basename = re.sub(
                        r'\.tsv\.gz$', 
                        f'_{baf_idx}.tsv.gz', 
                        self.__class__.save_paths_basenames[gdf_key],
                    )
                    assert basename.endswith('.tsv.gz')
                    pkl_basename = re.sub(r'\.tsv\.gz$', '.pickle', basename)

                    savepath = os.path.join(topdir, pkl_basename)
                    with open(savepath, 'wb') as outfile:
                        pickle.dump(gdf.df, outfile)

    def load_nonbafidx_gdf_aspickle(self, gdf_key, topdir):
        basename = self.__class__.save_paths_basenames[gdf_key]
        assert basename.endswith('.tsv.gz')
        pkl_basename = re.sub(r'\.tsv\.gz$', '.pickle', basename)
        savepath = os.path.join(topdir, pkl_basename)

        if os.path.exists(savepath):
            with open(savepath, 'rb') as infile:
                df = pickle.load(infile)
            gdf_type = self.__class__.gdf_types[gdf_key]
            self.gdfs[gdf_key] = gdf_type.from_frame(df, refver=self.refver)

    def load_bafidx_gdf_aspickle(self, gdf_key, topdir):
        for baf_idx in self.get_baf_indexes():
            basename = re.sub(
                r'\.tsv\.gz$', 
                f'_{baf_idx}.tsv.gz', 
                self.__class__.save_paths_basenames[gdf_key],
            )
            assert basename.endswith('.tsv.gz')
            pkl_basename = re.sub(r'\.tsv\.gz$', '.pickle', basename)

            savepath = os.path.join(topdir, pkl_basename)
            if os.path.exists(savepath):
                with open(savepath, 'rb') as infile:
                    df = pickle.load(infile)
                gdf_type = self.__class__.gdf_types[gdf_key]
                self.gdfs[gdf_key][baf_idx] = gdf_type.from_frame(df, refver=self.refver)

    def save_corrected_gdfs(self, corrector_id, savepath):
        saved_dict = self.corrected_gdfs[corrector_id].copy()

        for key in self.__class__.corrected_gdf_types.keys():
            if key not in saved_dict:
                continue
            if saved_dict[key] is None:
                continue

            saved_dict[key] = saved_dict[key].df

        for key in self.__class__.corrected_gdf_types_bafidx.keys():
            if key not in saved_dict:
                continue
            subdic = saved_dict[key]
            for subkey in tuple(subdic.keys()):
                try:
                    subdic[subkey] = subdic[subkey].df
                except Exception as exc:
                    raise Exception(f'key={key}, subkey={subkey}') from exc

        with open(savepath, 'wb') as outfile:
            pickle.dump(saved_dict, outfile)

    def load_corrected_gdfs(self, savepath):
        with open(savepath, 'rb') as infile:
            corrected_gdfs = pickle.load(infile)

        for key, gdf_type in self.__class__.corrected_gdf_types.items():
            if key not in corrected_gdfs:
                continue
            if corrected_gdfs[key] is None:
                continue
            if not isinstance(corrected_gdfs[key], GDF):
                corrected_gdfs[key] = gdf_type.from_frame(corrected_gdfs[key], refver=self.refver)

        for key, gdf_type in self.__class__.corrected_gdf_types_bafidx.items():
            if key not in corrected_gdfs:
                continue
            subdic = corrected_gdfs[key]
            for subkey in tuple(subdic.keys()):
                subdic[subkey] = gdf_type.from_frame(subdic[subkey], refver=self.refver)

        self.corrected_gdfs[corrected_gdfs['id']] = corrected_gdfs

    def save_pickle_asdf(self, topdir, force=False):
        if os.path.exists(topdir):
            if force:
                shutil.rmtree(topdir)
            else:
                raise Exception(f'Output directory must not exist in advance')
        os.mkdir(topdir)

        self.log(f'Beginning saving to a directory, with each gdf pickled as a DataFrame')

        # simple and solution attrs
        self.save_simple_attrs(topdir)
        self.save_solution_attrs(topdir)

        # target region
        self.save_nonbafidx_gdf_aspickle('raw_target_region', topdir)
        self.save_nonbafidx_gdf_aspickle('target_region', topdir)
        self.save_nonbafidx_gdf_aspickle('excluded_region', topdir)

        # raw data
        self.save_nonbafidx_gdf_aspickle('depth_rawdata', topdir)
        self.save_nonbafidx_gdf_aspickle('baf_rawdata', topdir)
        #self.save_nonbafidx_gdf_aspickle('baf_hetalt_rawdata', topdir)

        # segment
        self.save_nonbafidx_gdf_aspickle('depth_segment', topdir)
        self.save_bafidx_gdf_aspickle('baf_segment', topdir)
        self.save_bafidx_gdf_aspickle('baf_noedit_segment', topdir)
        self.save_bafidx_gdf_aspickle('baf_edited_segment', topdir)
        self.save_bafidx_gdf_aspickle('baf_rmzero_noedit_segment', topdir)
        self.save_nonbafidx_gdf_aspickle('merged_segment', topdir)

        # corrected_gdfs
        cgdfs_dir = os.path.join(topdir, 'corrected_gdfs')
        os.mkdir(cgdfs_dir)
        for sid in self.corrected_gdfs.keys():
            savepath = os.path.join(cgdfs_dir, sid)
            self.save_corrected_gdfs(sid, savepath)

        # solution

        self.log(f'Finished saving')

    @classmethod
    def load_pickle_asdf(
        cls, topdir,
        load_all=True,
        load_rawdata=False,
        load_segment=False,
        load_corrected=False,
    ):
        if not os.path.exists(topdir):
            raise Exception(f'Directory does not exist: {topdir}')

        if load_all:
            load_rawdata = True
            load_segment = True
            load_corrected = True

        result = cls()

        # simple attrs
        result.load_simple_attrs(topdir)
        result.load_solution_attrs(topdir)

        # initiate gdf data dict
        result.init_gdfs()

        # initiate cache
        result.init_cache()

        # target region
        result.load_nonbafidx_gdf_aspickle('raw_target_region', topdir)
        result.load_nonbafidx_gdf_aspickle('target_region', topdir)
        result.load_nonbafidx_gdf_aspickle('excluded_region', topdir)

        # raw data
        if load_rawdata:
            result.load_nonbafidx_gdf_aspickle('depth_rawdata', topdir)
            result.load_nonbafidx_gdf_aspickle('baf_rawdata', topdir)
            #result.load_nonbafidx_gdf_aspickle('baf_hetalt_rawdata', topdir)

        # segment
        if load_segment:
            result.load_nonbafidx_gdf_aspickle('depth_segment', topdir)
            result.load_bafidx_gdf_aspickle('baf_segment', topdir)
            result.load_bafidx_gdf_aspickle('baf_noedit_segment', topdir)
            result.load_bafidx_gdf_aspickle('baf_edited_segment', topdir)
            result.load_bafidx_gdf_aspickle('baf_rmzero_noedit_segment', topdir)

        result.load_nonbafidx_gdf_aspickle('merged_segment', topdir)

        # corrected_gdfs
        if load_corrected:
            cgdfs_dir = os.path.join(topdir, 'corrected_gdfs')
            for fname in os.listdir(cgdfs_dir):
                savepath = os.path.join(cgdfs_dir, fname)
                result.load_corrected_gdfs(savepath)

        # genomeplotter
        result.init_genomeplotter(**result.genomeplotter_kwargs)

        # plotdata
        #result.plotdata_cache = dict()

        return result

#    def save_pickle_asdf_notworking(self, filepath):
#        """This does not work"""
#        # sanitycheck
#        if not filepath.endswith('.pickle'):
#            raise Exception(f'File name must end with ".pickle"')
#
#        # make a dict of df from self.gdfs
#        self.log(f'Converting gdfs into dfs')
#
#        def converter(val):
#            if val is None:
#                return None
#            else:
#                return val.df
#
#        gdfs_as_df = dict()
#        for key, val in self.gdfs.items():
#            if isinstance(val, dict):
#                gdfs_as_df[key] = dict()
#                for subkey, subval in val.items():
#                    gdfs_as_df[key][subkey] = converter(subval)
#            else:
#                gdfs_as_df[key] = converter(val)
#
#        original_gdfs = self.gdfs
#        self.gdfs = gdfs_as_df
#
#        # main
#        self.log(f'Beginning saving as pickle')
#        
#        if os.path.exists(filepath):
#            logutils.log(f'Overwriting an existing file', level='warning')
#
#        with open(filepath, 'wb') as outfile:
#            pickle.dump(self, outfile)
#
#        self.log(f'Finished saving')
#
#        # revert to dict of gdfs
#        self.gdfs = original_gdfs

    ######################
    # textfile save/load #
    ######################

    def get_simple_attrs_savepath(self, topdir):
        return os.path.join(
            topdir, 
            self.__class__.save_paths_basenames['simple_attrs'],
        )

    def save_simple_attrs(self, topdir):
        with open(self.get_simple_attrs_savepath(topdir), 'wt') as outfile:
            json.dump(self.simple_attrs, outfile)

    def load_simple_attrs(self, topdir):
        savepath = self.get_simple_attrs_savepath(topdir)
        if os.path.exists(savepath):
            with open(savepath, 'rt') as infile:
                self.simple_attrs = json.load(infile)
        else:
            self.init_simple_attrs()

    ###

    def get_solution_attrs_savepath(self, topdir):
        return os.path.join(
            topdir, 
            self.__class__.save_paths_basenames['solution_attrs'],
        )

    def save_solution_attrs(self, topdir):
        with open(self.get_solution_attrs_savepath(topdir), 'wb') as outfile:
            pickle.dump(self.solution_attrs, outfile)

    def load_solution_attrs(self, topdir):
        savepath = self.get_solution_attrs_savepath(topdir)
        if not os.path.exists(savepath):
            self.init_solution_attrs()
        else:
            with open(savepath, 'rb') as infile:
                self.solution_attrs = pickle.load(infile)

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
                    r'\.tsv\.gz$', 
                    f'_{baf_idx}.tsv.gz', 
                    self.__class__.save_paths_basenames[gdf_key],
                )
                savepath = os.path.join(topdir, basename)
                gdf.write_tsv(savepath)

    def load_bafidx_gdf(self, gdf_key, topdir, GDF_class):
        for baf_idx in self.get_baf_indexes():
            basename = re.sub(
                r'\.tsv\.gz$', 
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
        self.save_solution_attrs(topdir)

        # target region
        self.save_nonbafidx_gdf('raw_target_region', topdir)
        self.save_nonbafidx_gdf('target_region', topdir)
        self.save_nonbafidx_gdf('excluded_region', topdir)

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

        # simple and solution attrs
        result.load_simple_attrs(topdir)
        result.load_solution_attrs(topdir)

        # initiate gdf data dict
        result.init_gdfs()

        # initiate cache
        result.init_cache()

        # target region
        result.load_nonbafidx_gdf('raw_target_region', topdir, GDF)
        result.load_nonbafidx_gdf('target_region', topdir, GDF)
        result.load_nonbafidx_gdf('excluded_region', topdir, GDF)

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

    #@property
    #def chromdict(self):
    #    return refgenome.get_chromdict(self.simple_attrs["refver"])

    def get_chromdict(self):
        return refgenome.get_chromdict(self.refver)

    #################
    # data fetchers #
    #################

    def get_raw_target_region(self):
        return self.gdfs['raw_target_region']

    def get_target_region(self):
        return self.gdfs['target_region']

    #def get_excluded_region(self):
    #    return self.gdfs['excluded_region']

    def get_depth_rawdata(self, corrector_id=None):
        if corrector_id is not None:
            return self.corrected_gdfs[corrector_id]['depth_rawdata']
        else:
            return self.gdfs['depth_rawdata']

    def get_baf_rawdata(self, corrector_id=None, baf_idx=None, rmzero=False):
        if corrector_id is None:
            baf_raw_gdf = self.gdfs['baf_rawdata']
        else:
            baf_raw_gdf = self.corrected_gdfs[corrector_id]['baf_rawdata']

        if baf_raw_gdf is None:
            return None
        else:
            if baf_idx is None:
                return baf_raw_gdf
            else:
                baf_raw_gdf = baf_raw_gdf.choose_annots(baf_idx)
                if rmzero:
                    selector = baf_raw_gdf[baf_idx] > 0
                    baf_raw_gdf = baf_raw_gdf.loc[selector, :]

                return baf_raw_gdf

    def get_filtered_baf_rawdata(self, corrector_id=None):
        if corrector_id is None:
            return self.gdfs['filtered_baf_rawdata']
        else:
            return self.corrected_gdfs[corrector_id]['filtered_baf_rawdata']

    def get_depth_segment(self, corrector_id=None):
        if corrector_id is not None:
            return self.corrected_gdfs[corrector_id]['depth_segment']
        else:
            return self.gdfs['depth_segment']

    def get_baf_segment_dict(self, corrector_id=None):
        if corrector_id is None:
            return self.gdfs['baf_segment']
        else:
            return self.corrected_gdfs[corrector_id]['baf_segment']

    def get_baf_segment(self, baf_idx, corrector_id=None):
        return self.get_baf_segment_dict(corrector_id)[baf_idx]

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

    def get_merged_segment(self, corrector_id=None):
        if corrector_id is None:
            return self.gdfs['merged_segment']
        else:
            return self.corrected_gdfs[corrector_id]['merged_segment']

    def get_merged_cnv_segment(self, corrector_id=None):
        if corrector_id is None:
            return self.gdfs['merged_merged_segment']
        else:
            return self.corrected_gdfs[corrector_id]['merged_merged_segment']

    def get_K(self, corrector_id):
        return self.corrected_gdfs[corrector_id]['K']

    def get_cellularity(self, corrector_id):
        return self.corrected_gdfs[corrector_id]['cellularity']

    def get_ploidy(self, corrector_id):
        return self.corrected_gdfs[corrector_id]['ploidy']

    def get_solution_params(self, corrector_id):
        return {
            key: self.corrected_gdfs[corrector_id][key]
            for key in ['K', 'cellularity', 'ploidy']
        }

    ###########################
    # depth data modification #
    ###########################

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

    @get_deco_logging(f'default germline CN/B annotation to depth rawdata')
    @deco_nproc
    def add_default_CNg_Bg_to_depth_rawdata(self, nproc=0):
        cnvcall.add_default_CNg_Bg(
            gdf=self.get_depth_rawdata(),
            is_female=self.is_female,
            inplace=True,
            nproc=nproc,
        )

#    def make_depth_segment_base(self, depth_rawdata, annot_colname, nproc, **kwargs):
#        depth_segment = depth_rawdata.get_segment(
#            annot_colname=annot_colname,
#            drop_annots=True,
#            nproc=nproc,
#            **kwargs,
#        )
#      
#        # fill gaps (when annotation values contain nan, segments may be gapped)
#        depth_segment = depth_segment.fill_gaps(edit_first_last=True)
#
#        # trim into target region
#        raw_target_region = self.get_raw_target_region()
#        if raw_target_region is not None:
#            depth_segment = depth_segment.intersect(raw_target_region, nproc=nproc)
#
#        return depth_segment

    @get_deco_logging(f'depth segmentation')
    @deco_nproc
    def make_depth_segment(
        self, 
        nproc=0, 
        segment_kwargs=dict(),
        winsorize=libdepth.DEFAULT_WINSORIZE,
    ):
        depth_rawdata = self.get_depth_rawdata()
        depth_segment = depth_rawdata.get_segment(
            target_region=self.get_raw_target_region(),
            rawdepth=False,
            nproc=nproc,
            **segment_kwargs,
        )
        depth_segment.add_rawdata_info(
            depth_rawdata, 
            merge_methods=['mean', 'std'],
            rawdepth=False,
            nproc=nproc,
            winsorize=winsorize,
        )
        self.gdfs['depth_segment'] = depth_segment

#        depth_rawdata = self.get_depth_rawdata()
#        annot_colname = (
#            DepthRawDF.depth_colname
#            if rawdepth else
#            DepthRawDF.norm_depth_colname
#        )
#        depth_segment = self.make_depth_segment_base(
#            depth_rawdata=depth_rawdata, 
#            annot_colname=annot_colname, 
#            nproc=nproc, 
#            **kwargs,
#        )
#        self.gdfs['depth_segment'] = depth_segment

#    @get_deco_logging(f'depth segment modification')
#    @deco_nproc
#    def edit_depth_segment(self, nproc=0, readlength=151, MQ_cutoff=(40, None)):
#        self.make_normalized_depth()
#        self.add_MQ_to_depth_rawdata(readlength=readlength, nproc=nproc)
#        self.add_rawdata_to_depth_segment(rawdepth=False, nproc=nproc)
#        self.add_filter_to_depth_segment(MQ_cutoff=MQ_cutoff)

    ###

    #def make_normalized_depth(self):
    #    self.get_depth_rawdata().set_normalized_depth()

#    @get_deco_logging(f'rawdata annotation to depth segment')
#    @deco_nproc
#    def add_rawdata_to_depth_segment(
#        self, rawdepth=False, nproc=0,
#    ):
#        depth_seg_gdf = self.get_depth_segment()
#        depth_raw_gdf = self.get_depth_rawdata()
#        depth_seg_gdf.add_rawdata_info(
#            depth_raw_gdf, 
#            merge_methods=['mean', 'std'],
#            rawdepth=rawdepth,
#            nproc=nproc,
#        )

#    def add_filter_to_depth_rawdata_base(self, depth_rawdata, MQ_cutoff, norm_depth_cutoff, raw_depth_cutoff):
#        depth_rawdata.add_filter(
#            MQ_cutoff=MQ_cutoff, 
#            norm_depth_cutoff=norm_depth_cutoff, 
#            raw_depth_cutoff=raw_depth_cutoff,
#        )

#    @get_deco_logging(f'filter addition to depth rawdata')
#    def add_filter_to_depth_rawdata(self, MQ_cutoff=(40, None), norm_depth_cutoff=(0, None), raw_depth_cutoff=(0, None)):
#        self.add_filter_to_depth_rawdata_base(
#            depth_rawdata=self.get_depth_rawdata(), 
#            MQ_cutoff=MQ_cutoff, 
#            norm_depth_cutoff=norm_depth_cutoff, 
#            raw_depth_cutoff=raw_depth_cutoff,
#        )

    @get_deco_logging(f'filter addition to depth segment')
    def add_filter_to_depth_segment(self, MQ_cutoff=(40, None), norm_depth_cutoff=(0, None)):
        self.get_depth_segment().add_filter(MQ_cutoff=MQ_cutoff, norm_depth_cutoff=norm_depth_cutoff)

    def make_upscaled_depth(self, binsize):
        pass

    #########################
    # baf data modification #
    #########################

    @get_deco_logging(f'BAF segmentation')
    @deco_nproc
    def make_baf_segment(self, nproc=0, segment_kwargs=dict()):
        assert self.is_germline
        self.gdfs['baf_segment'], self.gdfs['filtered_baf_rawdata'] = self.get_baf_rawdata().get_segment_dict(
            target_region=self.get_raw_target_region(), 

            germline_baf_rawdata=None,

            germline_CN_gdf=None, 
            is_female=None,
            skip_hetalt_filter=False,

            return_filtered_rawdata=True,

            nproc=nproc,

            **segment_kwargs,
        )

    @get_deco_logging(f'BAF segmentation')
    @deco_nproc
    def make_baf_noedit_segment(self, nproc=0, rmzero=False):
        for baf_idx in self.get_baf_indexes():
            self.log(f'Beginning {baf_idx}, rmzero={rmzero}')

            # 1. make segment
            baf_raw_gdf = self.get_baf_rawdata(baf_idx=baf_idx, rmzero=rmzero)
            baf_seg_gdf = baf_raw_gdf.get_segment(baf_idx, drop_annots=True, nproc=nproc)

            # 2. add rawdata
            self.add_rawdata_to_bafseg_base(
                baf_seg_gdf=baf_seg_gdf, 
                baf_raw_gdf=baf_raw_gdf, 
                baf_idx=baf_idx,
                distinfo_kwargs=dict(),
                nproc=nproc,
            )

            # 3. isec with raw target region (panel)
            raw_target_region = self.get_raw_target_region()
            if raw_target_region is not None:
                baf_seg_gdf = baf_seg_gdf.intersect(raw_target_region)

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

#    @get_deco_logging(f'rawdata annotation to baf segment')
#    @deco_nproc
#    def add_rawdata_to_baf_noedit_segment(self, nproc=0):
#        for baf_idx in self.get_baf_indexes():
#            self.add_rawdata_to_bafseg_base(
#                baf_seg_gdf=self.get_baf_noedit_segment(baf_idx), 
#                baf_raw_gdf=self.get_baf_rawdata(baf_idx), 
#                baf_idx=baf_idx,
#                distinfo_kwargs=dict(),
#                nproc=nproc,
#            )

    @deco_nproc
    @get_deco_logging(f'modifying BAF segment')
    def make_baf_edited_segment(
        self, 
        rmzero=False, 
        distinfo_kwargs=dict(),
        nproc=0,

        # low ndata merging
        ndata_cutoff=30, 
        merge_low_ndata=None,

        # corrected baf cutoff
        corrected_baf_cutoff=0.41,

        # filtering
        width_cutoff=None,
        #center_cutoff=(0.4, np.inf),
        center_cutoff=None,
    ):
        if merge_low_ndata is None:
            if self.mode == 'wgs':
                merge_low_ndata = True
            elif self.mode == 'panel':
                merge_low_ndata = False

        for baf_idx in self.get_baf_indexes():
            logutils.log(f'Beginning job for baf index {baf_idx}')

            baf_raw_gdf = self.get_baf_rawdata(baf_idx=baf_idx, rmzero=rmzero)

            # 2. fill gap
            baf_edited_seg_gdf = self.get_baf_noedit_segment(baf_idx, rmzero=rmzero).copy()  # this is annotated with rawdata
            baf_edited_seg_gdf = baf_edited_seg_gdf.fill_gaps(edit_first_last=True)
            raw_target_region = self.get_raw_target_region()
            if raw_target_region is not None:
                baf_edited_seg_gdf = baf_edited_seg_gdf.intersect(raw_target_region)

            # 3. merge low ndata segments
            if merge_low_ndata:
                logutils.log(f'Beginning handling of low ndata segments')

                cycle = 0
                while True:
                    cycle += 1
                    logutils.log(f'Cycle {cycle}')

                    # merge between low segments
                    logutils.log(f'Beginning merging with themselves')
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
                    logutils.log(f'Finished merging with themselves')

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

                logutils.log(f'Finished handling of low ndata segments')

            # 4. add corrected baf column
            baf_edited_seg_gdf.add_corrected_baf_simple(cutoff=corrected_baf_cutoff)

            # 5. add filter column
            baf_edited_seg_gdf.add_filter(width_cutoff=width_cutoff, center_cutoff=center_cutoff)

            # result
            self.gdfs['baf_edited_segment'][baf_idx] = baf_edited_seg_gdf

    ###########################
    # merged segment handling #
    ###########################

    @deco_nproc
    @get_deco_logging(f'merging depth and baf segments')
    def make_merged_segment(self, nproc=0):
        assert self.is_germline
        
        if not self.check_bafinfo_exists():
            logutils.log(f'Making merged segment without BAF information', level='warning')

        depth_segdf = self.get_depth_segment()
        if self.check_bafinfo_exists():
            #baf_edited_segdf_dict = self.get_baf_edited_segment_dict()
            #sub_seg_gdfs = [depth_segdf] + list(baf_edited_segdf_dict.values())
            sub_seg_gdfs = [depth_segdf] + list(self.get_baf_segment_dict().values())
        else:
            sub_seg_gdfs = [depth_segdf]

        merged_seg_gdf = CNVSegDF.from_segments(sub_seg_gdfs, nproc=nproc)

        if self.check_bafinfo_exists():
            cols_to_drop = (
                [DepthSegDF.filter_colname]
                #+ [x.get_filter_colname() for x in baf_edited_segdf_dict.values()]
                + [x.get_filter_colname() for x in self.get_baf_segment_dict().values()]
            )
        else:
            cols_to_drop = [DepthSegDF.filter_colname]
        merged_seg_gdf.drop_annots(cols_to_drop, inplace=True)

        # trim into target region
        raw_target_region = self.get_raw_target_region()
        if raw_target_region is not None:
            merged_seg_gdf = merged_seg_gdf.intersect(raw_target_region, nproc=nproc)

        # add germline solution
        merged_seg_gdf.add_clonal_solution_germline(self.is_female, self.ploidy)
        #self.add_clonal_solution_to_germline_sample()

        self.gdfs['merged_segment'] = merged_seg_gdf

    #@get_deco_logging(f'default germline CN/B annotation to merged segment')
    #@deco_nproc
    #def add_default_CNg_Bg_to_merged_segment(self, nproc=0):
    #    cnvcall.add_default_CNg_Bg(
    #        gdf=self.get_merged_segment(),
    #        is_female=self.is_female,
    #        inplace=True,
    #        nproc=nproc,
    #    )

    def make_merged_cnv_segment(self, corrector_id=None):
        # make
        nonmerged_cnv_gdf = self.get_merged_segment(corrector_id=corrector_id)
        merged_cnv_gdf = nonmerged_cnv_gdf.merge_by_CN()

        # intersect with target region
        raw_target = self.get_raw_target_region()
        if raw_target is not None:
            merged_cnv_gdf = merged_cnv_gdf.intersect(raw_target)
        excl_region = self.get_corrected_depth_excluded_region(corrector_id=corrector_id)
        merged_cnv_gdf = merged_cnv_gdf.subtract(excl_region)

        # annotate with germline CN values
        merged_cnv_gdf = merged_cnv_gdf.join(
            nonmerged_cnv_gdf, 
            (
                [nonmerged_cnv_gdf.get_clonal_CN_colname(germline=True)]
                + [
                    nonmerged_cnv_gdf.get_clonal_B_colname(baf_index=x, germline=True)
                    for x in nonmerged_cnv_gdf.get_baf_indexes()
                ]
            ), 
            merge='longest',
            suffixes={'longest': ''},
        )

        # assign
        if corrector_id is None:
            self.gdfs['merged_merged_segment'] = merged_cnv_gdf
        else:
            self.corrected_gdfs[corrector_id]['merged_merged_segment'] = merged_cnv_gdf

    ##########################
    # corrected gdfs - maker #
    ##########################

    def make_corrected_gdfs_init(self, corrector_id):
        self.corrected_gdfs[corrector_id] = {
            'id': corrector_id,
            'solution': None,
            #'cellularity': None,
            #'K': None,
            #'ploidy': None,

            'depth_rawdata': None,
            'depth_segment': None,
            'depth_excluded_region': None,

            'baf_rawdata': None,
            'filtered_baf_rawdata': None,
            'baf_segment': None,

            'merged_segment': None,
            'merged_merged_segment': None,
        }

    @deco_nproc
    def make_corrected_gdfs_depth(
        self, 
        corrector_id, 
        corrector_gdf, 
        MQ_cutoff=(40, None), 
        norm_depth_cutoff=None,
        raw_depth_cutoff=(0, None),
        nproc=0,
        segment_kwargs=dict(),
        winsorize=libdepth.DEFAULT_WINSORIZE,
    ):
        assert raw_depth_cutoff[0] is not None
        corrected_depth_rawdata = self.get_depth_rawdata().copy()
        corrected_depth_rawdata.sort()
        corrector_gdf.sort()

        # 1. sanitycheck
        assert corrected_depth_rawdata.drop_annots() == corrector_gdf.drop_annots()

        # 2. add corrected norm depth
        self.log(f'Corrected data "{corrector_id}" - Calculating corrected norm depth')
        corrector_depth = corrector_gdf.norm_depth
        corrector_depth[corrector_depth == 0] = np.nan  # corrector depth may contain zeros (e.g. leading N region where CN call is not 0)
        corrected_norm_depth = corrected_depth_rawdata.norm_depth / corrector_depth
        corrected_norm_depth[
            corrected_depth_rawdata.raw_depth < raw_depth_cutoff[0]
        ] = np.nan
            # Where raw panel seq depth is too low, results in artefactual spikes in depth segment data
        corrected_depth_rawdata[DepthRawDF.norm_depth_colname] = corrected_norm_depth

        # 3. make segment from corrected norm depth
        self.log(f'Corrected data "{corrector_id}" - Making depth segment')
        corrected_depth_segment = corrected_depth_rawdata.get_segment(
            target_region=self.get_raw_target_region(),
            rawdepth=False,
            nproc=nproc,
            **segment_kwargs,
        )

        # 4. add rawdata information to segment
        self.log(f'Corrected data "{corrector_id}" - Adding rawdata to depth segment')
        corrected_depth_segment.add_rawdata_info(
            corrected_depth_rawdata, 
            merge_methods=['mean', 'std'],
            nproc=nproc,
            winsorize=winsorize,
        )

        # 5. extract excluded region
        self.log(f'Corrected data "{corrector_id}" - Creating excluded region from depth rawdata')
        corrected_depth_rawdata.add_filter(
            MQ_cutoff=MQ_cutoff, 
            norm_depth_cutoff=norm_depth_cutoff, 
            raw_depth_cutoff=raw_depth_cutoff,
        )

        #self.add_filter_to_depth_rawdata_base(
        #    depth_rawdata=corrected_depth_rawdata, 
        #    MQ_cutoff=MQ_cutoff, 
        #    norm_depth_cutoff=norm_depth_cutoff, 
        #    raw_depth_cutoff=raw_depth_cutoff,
        #)

        excl_selector = np.logical_not(corrected_depth_rawdata.get_filter())
        excluded_region = corrected_depth_rawdata.drop_annots().loc[excl_selector, :]
        excluded_region = excluded_region.merge()

        # 6. assign
        self.corrected_gdfs.setdefault(corrector_id, dict())
        self.corrected_gdfs[corrector_id]['depth_rawdata'] = corrected_depth_rawdata
        self.corrected_gdfs[corrector_id]['depth_segment'] = corrected_depth_segment
        self.corrected_gdfs[corrector_id]['depth_excluded_region'] = excluded_region

    @deco_nproc
    def make_corrected_gdfs_baf_germlinesample(
        self, 
        corrector_id, 
        corrector_gdf, 
        nproc=0,
    ):
        assert self.is_germline

        self.log(f'Applying corrector to VAF rawdata')

        # do join
        joined_baf_rawdata = self.get_baf_rawdata().variant_join(corrector_gdf, how='left', nproc=nproc)

        # collect corrector values, excluding self id
        corrector_values = dict()
        selected_cols_list = list()
        for pf in joined_baf_rawdata.allele_columns:
            selected_cols = [
                x for x in joined_baf_rawdata.columns
                if (x.startswith(f'{pf}_vaf_') and (x != f'{pf}_vaf_{self.sampleid}'))
            ]
            selected_cols.sort()
            selected_cols_list.append(
                [re.sub(r'^.+_vaf_', '', x) for x in selected_cols]
            )

            corrector_values[pf] = joined_baf_rawdata[selected_cols]

        assert all((x == y) for x, y in tools.pairwise(selected_cols_list))

        # choose a random not-nan index from each row
        notna_col_indexes = list()
        all_notna = ~(np.isnan(corrector_values['REF']))
        rng = np.random.default_rng()
        for notna_row in all_notna:
            not_na_indexes = np.nonzero(notna_row)[0]
            if len(not_na_indexes) == 0:
                notna_col_indexes.append(0)
            else:
                notna_col_indexes.append(rng.choice(not_na_indexes))
        notna_col_indexes = np.asarray(notna_col_indexes)

        # select corrector values according to the not-nan indexes
        selected_corrector_values = dict()
        row_indexes = np.arange(corrector_values['REF'].shape[0])
        for key, val in corrector_values.items():
            selected_corrector_values[key] = val[(row_indexes, notna_col_indexes)]

        # drop unneccessary columns
        chosen_cols = [x + '_vaf' for x in joined_baf_rawdata.allele_columns]
        corrected_baf_rawdata = joined_baf_rawdata.choose_annots(chosen_cols)

        # do division
        for key, val in selected_corrector_values.items():
            corrected_baf_rawdata[f'{key}_vaf'] = corrected_baf_rawdata[f'{key}_vaf'] / val

        # remove rows where corrector could not be found
        selector = ~np.isnan(corrected_baf_rawdata['REF_vaf'])
        corrected_baf_rawdata = corrected_baf_rawdata.loc[selector, :]

        # normalize VAF values and make BAF values
        corrected_baf_rawdata.add_baf(fit_self_vafs=True)

        # assign
        if corrector_id not in self.corrected_gdfs:
            self.make_corrected_gdfs_init(corrector_id)
        self.corrected_gdfs[corrector_id]['baf_rawdata'] = corrected_baf_rawdata

    def make_corrected_gdfs_baf_pairedgermline(
        self, 
        corrector_id, 
        germline_baf_rawdata,
    ):
        # 1. synchronize coordinates between tumor and germline baf rawdata
        corrected_baf_rawdata = self.get_baf_rawdata().copy()
        corrected_baf_rawdata.sort()
        germline_baf_rawdata.sort()
        assert (
            corrected_baf_rawdata[corrected_baf_rawdata.nonannot_columns] 
            == germline_baf_rawdata[germline_baf_rawdata.nonannot_columns] 
        ).all()

        # 2. keep only hetalt positions determined from paired germline data
        germline_baf_rawdata.add_baf_hetalt_flag()
        ishet = germline_baf_rawdata.get_hetalt_selector()
        germline_baf_rawdata = germline_baf_rawdata.loc[ishet, :]
        corrected_baf_rawdata = corrected_baf_rawdata.loc[ishet, :]

        # 2. do division
        corrected_baf_rawdata[corrected_baf_rawdata.vaf_columns] = (
            corrected_baf_rawdata[corrected_baf_rawdata.vaf_columns]
            / germline_baf_rawdata[germline_baf_rawdata.vaf_columns]
        )

        # 3. normalize VAF values and make BAF values
        corrected_baf_rawdata.add_baf(fit_self_vafs=True)

        # 4. assign
        if corrector_id not in self.corrected_gdfs:
            self.make_corrected_gdfs_init(corrector_id)
        self.corrected_gdfs[corrector_id]['baf_rawdata'] = corrected_baf_rawdata

    @deco_nproc
    def make_corrected_gdfs_baf(
        self, 
        corrector_id, 
        corrector_gdf, 
        germline_CN_gdf=None,
        germline_baf_rawdata=None,
        nproc=0,
        segment_kwargs=dict(),
    ):
        if self.is_germline:
            assert germline_baf_rawdata is None
            self.make_corrected_gdfs_baf_germlinesample(
                corrector_id=corrector_id, 
                corrector_gdf=corrector_gdf, 
                nproc=nproc,
            )
        else:
            assert germline_baf_rawdata is not None

            # HETALT FILTERING IS DONE WITHIN THIS METHOD
            self.make_corrected_gdfs_baf_pairedgermline(
                corrector_id=corrector_id, 
                germline_baf_rawdata=germline_baf_rawdata,
            )

        # segment
        self.log(f'Making BAF segment from the corrected data')

        corrected_baf_rawdata = self.corrected_gdfs[corrector_id]['baf_rawdata']
        skip_hetalt_filter = (not self.is_germline)  # In order to skip hetalt filtering which is already done in case of non-germline sample
        corrected_baf_segment_dict, filtered_baf_rawdata = corrected_baf_rawdata.get_segment_dict(
            target_region=self.get_raw_target_region(),

            germline_baf_rawdata=germline_baf_rawdata,

            germline_CN_gdf=germline_CN_gdf, 
            is_female=self.is_female,
            skip_hetalt_filter=skip_hetalt_filter,

            return_filtered_rawdata=True,

            nproc=nproc, 

            **segment_kwargs,
        )
        self.corrected_gdfs[corrector_id]['baf_segment'] = corrected_baf_segment_dict
        self.corrected_gdfs[corrector_id]['filtered_baf_rawdata'] = filtered_baf_rawdata

    @deco_nproc
    def make_corrected_gdfs_baf_old(
        self, 
        corrector_id, 
        corrector_gdf, 
        germline_CN_gdf=None,
        n_sample_cutoff=5,
        nproc=0,
    ):
        # rawdata
        self.log(f'Applying corrector to BAF rawdata')
        corrected_baf_rawdata = self.get_baf_rawdata().variant_join(
            corrector_gdf,
            how='left',
        )
        n_sample_unselector = (corrected_baf_rawdata['n_samples'] < n_sample_cutoff)  # nan is removed here
        for colname in corrected_baf_rawdata.vaf_columns:
            numer = corrected_baf_rawdata[colname]
            denom = corrected_baf_rawdata[colname + '_corrector']
            denom[
                np.logical_or(n_sample_unselector, np.isnan(denom))
            ] = 1
            corrected_baf_rawdata[colname] = numer / denom

        corrected_baf_rawdata.add_baf(fit_self_vafs=True)

        self.corrected_gdfs.setdefault(corrector_id, dict())
        self.corrected_gdfs[corrector_id]['baf_rawdata'] = corrected_baf_rawdata

        # segment
        self.log(f'Making BAF segment from the corrected data')
        corrected_baf_segment_dict = corrected_baf_rawdata.get_segment_dict(
            target_region=self.get_raw_target_region(),
            germline_CN_gdf=germline_CN_gdf, 
            is_female=self.is_female,
            nproc=nproc, 
        )
        self.corrected_gdfs[corrector_id]['baf_segment'] = corrected_baf_segment_dict

    @deco_nproc
    def make_corrected_gdfs_merged_segment(
        self, 
        corrector_id, 
        germline_CN_gdf=None,
        nproc=0,
    ):
        self.log(f'Corrected data "{corrector_id}" - Making merged segment')

        # 1. make merged segment
        corrected_gdfs = self.corrected_gdfs[corrector_id]
        sub_seg_gdfs = (
            [corrected_gdfs['depth_segment']] 
            + list(corrected_gdfs['baf_segment'].values())
        )
        merged_seg_gdf = CNVSegDF.from_segments(sub_seg_gdfs, nproc=nproc)

        # 2. drop filter columns
        cols_to_drop = (
            [DepthSegDF.filter_colname]
            + [x.get_filter_colname() for x in corrected_gdfs['baf_segment'].values()]
        )
        merged_seg_gdf.drop_annots(cols_to_drop, inplace=True)

        # 3. subtract excluded region
        merged_seg_gdf = merged_seg_gdf.subtract(corrected_gdfs['depth_excluded_region'])

        # 4. intersect with target region
        raw_target_region = self.get_raw_target_region()
        if raw_target_region is not None:
            merged_seg_gdf = merged_seg_gdf.intersect(raw_target_region)

        # 5. add germline CN
        if self.is_germline:
            merged_seg_gdf.add_clonal_solution_germline(self.is_female, self.ploidy)
        else:
            if germline_CN_gdf is None:
                self.make_corrected_gdfs_add_germline_CN_without_pairednormal(
                    merged_seg_gdf,
                    nproc=nproc,
                )
            else:
                merged_seg_gdf = self.make_corrected_gdfs_add_germline_CN(
                    cnv_seg_gdf=merged_seg_gdf,
                    paired_germline_CN_gdf=germline_CN_gdf,
                    nproc=nproc,
                )

        # result
        corrected_gdfs['merged_segment'] = merged_seg_gdf

    @deco_nproc
    def make_corrected_gdfs_add_germline_CN(
        self, 
        cnv_seg_gdf,
        paired_germline_CN_gdf,
        nproc=0,
    ):
        #self.log(f'Corrected data "{corrector_id}" - Adding germline CN to merged segment')
        #cnv_seg_gdf = self.corrected_gdfs[corrector_id]['merged_segment']

        # add CN of paired germline sample
        gCN_colname = paired_germline_CN_gdf.get_clonal_CN_colname(germline=False)
        gB_colname_dict = paired_germline_CN_gdf.get_clonal_B_colname_dict()
        added_cols = [gCN_colname] + list(gB_colname_dict.values())
        cnv_seg_gdf = cnv_seg_gdf.join(
            paired_germline_CN_gdf,
            how='left',
            right_gdf_cols=added_cols,
            merge='longest',
            overlapping_length=True,
            omit_N=True,
            suffixes={'longest': ''},
            nproc=nproc,
        )
        cnv_seg_gdf.assign_clonal_CN(
            data=cnv_seg_gdf[gCN_colname],
            germline=True,
        )
        for baf_index, colname in gB_colname_dict.items():
            cnv_seg_gdf.assign_clonal_B(
                data=cnv_seg_gdf[colname],
                baf_index=baf_index,
                germline=True,
            )

        # remove unused columns
        germline_CN_columns = (
            [cnv_seg_gdf.get_clonal_CN_colname(germline=True)]
            + [
                cnv_seg_gdf.get_clonal_B_colname(baf_index, germline=True)
                for baf_index in cnv_seg_gdf.get_baf_indexes()
            ]
        )
        cnv_seg_gdf.drop_annots(
            set(added_cols).difference(germline_CN_columns), 
            inplace=True,
        )

        return cnv_seg_gdf

        #self.corrected_gdfs[corrector_id]['merged_segment'] = cnv_seg_gdf

    @deco_nproc
    def make_corrected_gdfs_add_germline_CN_without_pairednormal(
        self, 
        cnv_seg_gdf,
        #corrector_id, 
        nproc=0,
    ):
        #self.log(f'Corrected data "{corrector_id}" - Adding germline CN to merged segment')

        #cnv_seg_gdf = self.corrected_gdfs[corrector_id]['merged_segment']

        cnvcall.add_default_CNg_Bg(
            gdf=cnv_seg_gdf,
            is_female=self.is_female,
            inplace=True,
            nproc=nproc,
        )
        cnv_seg_gdf.assign_clonal_CN(
            cnv_seg_gdf[cnvcall.DEFAULT_CNG_COLNAME],
            germline=True,
        )
        cnv_seg_gdf.assign_clonal_B(
            cnv_seg_gdf[cnvcall.DEFAULT_BG_COLNAME],
            baf_index='baf0',
            germline=True,
        )
        cnv_seg_gdf.drop_annots(
            [cnvcall.DEFAULT_CNG_COLNAME, cnvcall.DEFAULT_BG_COLNAME],
            inplace=True,
        )

        #self.corrected_gdfs[corrector_id]['merged_segment'] = cnv_seg_gdf

    ############################
    # corrected gdfs - fetcher #
    ############################

    #def get_corrected_depth_rawdata(self, corrector_id):
    #    return self.corrected_gdfs[corrector_id]['depth_rawdata']

    #def get_corrected_depth_segment(self, corrector_id):
    #    return self.corrected_gdfs[corrector_id]['depth_segment']

    def get_corrected_depth_excluded_region(self, corrector_id):
        return self.corrected_gdfs[corrector_id]['depth_excluded_region']

    #def get_corrected_depth_merged_segment(self, corrector_id):
    #    return self.corrected_gdfs[corrector_id]['merged_segment']

    ############
    # plotdata #
    ############

#    def make_plotdata_base(self, data, cache_key=None):
#        #self.log(f'Beginning plotdata generation for {cache_key}')
#        plotdata = self.genomeplotter.cconv.prepare_plot_data(data)
#        #self.plotdata_cache[cache_key] = plotdata
#        #self.log(f'Finished plotdata generation for {cache_key}')
#
#        return plotdata
#
#    @get_deco_logging(f'making BAF rawdata plotdata')
#    def make_baf_rawdata_plotdata(self, rmzero=False):
#        plotdata_dict = dict()
#        for baf_idx in self.get_baf_indexes():
#            data = self.get_baf_rawdata(baf_idx, rmzero=rmzero)
#            #plotdata_key = f'{sampletype}_baf_raw_{baf_idx}'
#            plotdata_dict[baf_idx] = self.make_plotdata_base(data)
#
#        return plotdata_dict
#
#    @get_deco_logging(f'making BAF segment plotdata')
#    def make_baf_noedit_segment_plotdata(self, rmzero=False):
#        segment_dict_key = ('baf_rmzero_noedit_segment' if rmzero else 'baf_noedit_segment')
#        data_dict = self.gdfs[segment_dict_key]
#        plotdata_dict = dict()
#        for baf_idx, data in data_dict.items():
#            #plotdata_key = f'{sampletype}_baf_segment_{baf_idx}'
#            plotdata = self.make_plotdata_base(data)
#            plotdata_dict[baf_idx] = plotdata
#
#        return plotdata_dict
#
#    @get_deco_logging(f'making BAF segment plotdata')
#    def make_baf_edited_segment_plotdata(self, rmzero=False):
#        data_dict = self.gdfs['baf_edited_segment']
#        plotdata_dict = dict()
#        for baf_idx, data in data_dict.items():
#            #plotdata_key = f'{sampletype}_baf_segment_{baf_idx}'
#            plotdata = self.make_plotdata_base(data)
#            plotdata_dict[baf_idx] = plotdata
#
#        return plotdata_dict
#
#    @get_deco_logging(f'making depth rawdata plotdata')
#    def make_depth_rawdata_plotdata(self):
#        data = self.get_depth_rawdata()
#        #plotdata_key = f'{sampletype}_depth_raw'
#        return self.make_plotdata_base(data)
#
#    @get_deco_logging(f'making depth segment plotdata')
#    def make_depth_segment_plotdata(self):
#        data = self.get_depth_segment()
#        #plotdata_key = f'{sampletype}_depth_segment'
#        return self.make_plotdata_base(data)
#
#    @get_deco_logging(f'making merged segment plotdata')
#    def make_merged_segment_plotdata(self, corrected=False, corrector_id=None):
#        data = self.get_merged_segment(corrected=corrected, corrector_id=corrector_id)
#        #plotdata_key = f'{sampletype}_depth_segment'
#        return self.make_plotdata_base(data)

    #########################
    # HIGH level ax drawers #
    #########################

    def make_mosaic(
        self,
        draw_depth=False,
        draw_baf=False,
        draw_MQ=False,
        draw_CN=False,
        draw_mut=False,
        draw_depth_peaks=False,
        with_legend=True,
    ):
        mosaic = list()
        if draw_CN:
            mosaic.append(['CN'])
        if draw_baf:
            mosaic.extend([[baf_idx] for baf_idx in self.get_baf_indexes()])
        if draw_depth:
            mosaic.append(['depth'])
        if draw_MQ:
            mosaic.append(['MQ'])
        if draw_mut:
            mosaic.append(['mut'])

        if with_legend:
            for x in mosaic:
                x.append(x[0] + '_legend')

        if draw_depth_peaks:
            for idx, x in enumerate(mosaic):
                if x[0] == 'depth':
                    x.append('depth_peaks')
                else:
                    x.append(f'blank{idx}')

        return mosaic

    def make_axd(
        self,
        draw_depth=False,
        draw_baf=False,
        draw_MQ=False,
        draw_CN=False,
        draw_mut=False,
        draw_depth_peaks=False,
        subplots_kwargs=dict(),
        figsize=None,
        with_legend=True,
    ):
        mosaic = self.make_mosaic(
            draw_depth=draw_depth,
            draw_baf=draw_baf,
            draw_MQ=draw_MQ,
            draw_CN=draw_CN,
            draw_mut=draw_mut,
            draw_depth_peaks=draw_depth_peaks,
            with_legend=with_legend,
        )

        default_subplots_kwargs = {
            'figsize': (12, 4 * len(mosaic)),
            'gridspec_kw': (
                {'hspace': 0.6, 'width_ratios': [1, 0.2]}
                if with_legend else
                {'hspace': 0.6}
            ),
        }
        if draw_depth_peaks:
            default_subplots_kwargs['gridspec_kw'].update(
                {'width_ratios': [1, 0.1], 'wspace': 0.02}
            )
        subplots_kwargs = default_subplots_kwargs | subplots_kwargs
        if figsize is not None:
            subplots_kwargs.update({'figsize': figsize})

        fig, axd = plt.subplot_mosaic(mosaic, **subplots_kwargs)
        for key, ax in axd.items():
            if key.endswith('_legend'):
                plotmisc.clear_ax(ax)

        for key, ax in axd.items():
            if key.startswith('blank'):
                ax.axis('off')

        return fig, axd

    @deco_nproc
    def draw(
        self,
        axd=None,
        genomeplotter=None,

        draw_all=True,
        draw_depth=False,
        draw_baf=False,
        draw_MQ=False,
        draw_mut=False,
        draw_CN=False,

        corrector_id=None,
        depthtype=None,
        chromwise_peaks=False,
        chromwise_peaks_density_kwargs=dict(),
        chromwise_peaks_line_kwargs=dict(),
        draw_excluded=False,

        vplist=None,
        #vcf_sampleid=None,

        #draw_depth_peaks=False,
        #draw_chromwise_peak=True,

        frac=0.1,
        draw_depth_rawdata=True,
        draw_baf_rawdata=True,
        filtered_baf=True,
        draw_MQ_rawdata=True,

        # figure title - only works when figure object is not already created
        figsize=None,
        subplots_kwargs=dict(),

        # kwargs for low level drawers
        draw_common_kwargs=dict(),

        depth_rawdata_kwargs=dict(),
        depth_segment_kwargs=dict(),
        #depth_peaks_kwargs=dict(),
        baf_rawdata_kwargs=dict(),
        baf_segment_kwargs=dict(),
        MQ_kwargs=dict(),
        CN_kwargs=dict(),
        excluded_kwargs=dict(),

        setup_axes=True,
        with_legend=True,

        #chromwise_peak_bw_method=1,
        #chromwise_peak_limit=None,

        # axes setting
        ylabel_prefix='',
        title=False,
        suptitle_kwargs=dict(),
        omit_title_mode=False,

        # multiprocessing
        nproc=0,
        verbose=True,
    ):
        if draw_all:
            draw_depth = True
            draw_baf = True
            draw_mut = True
            draw_CN = True

        # sanitycheck
        if not any([draw_depth, draw_baf, draw_MQ, draw_CN]):
            raise Exception(f'At least one of drawing flags must be set.')

        #assert not ((not draw_depth) and draw_depth_peaks)

        if genomeplotter is None:
            genomeplotter = self.genomeplotter

        # gdraw result
        gdraw_axresult_list = list()

        # make axd
        if axd is None:
            fig, axd = self.make_axd(
                draw_depth=draw_depth,
                draw_baf=draw_baf,
                draw_MQ=draw_MQ,
                draw_CN=draw_CN,
                draw_mut=draw_mut,
                draw_depth_peaks=False,
                subplots_kwargs=subplots_kwargs,
                figsize=figsize,
                with_legend=with_legend,
            )
        else:
            fig = next(iter(axd.values())).figure

        # prepare CN plotdata
        CN_gdf = self.get_merged_segment(corrector_id=corrector_id)
        CN_exists = (CN_gdf is not None)
        #CN_has_solution = (CN_exists and CN_gdf.check_has_CN())
        if CN_exists:
            CN_plotdata = genomeplotter.make_plotdata(
                CN_gdf, 
                log_suffix=' (Copy number)', 
                nproc=nproc,
                verbose=verbose,
            )

        # draw CN
        if draw_CN:
            axd['CN'].set_xticks([])
            axd['CN'].set_yticks([])

        # draw BAF
        if draw_baf:
            # segment
            #baf_segment_gdf_dict = self.get_baf_edited_segment_dict()
            baf_segment_gdf_dict = self.get_baf_segment_dict(corrector_id=corrector_id)
            baf_axd_keys = set(axd.keys()).intersection(self.get_baf_indexes())
            baf_axd = {k: axd[k] for k in baf_axd_keys}
            for baf_idx, baf_ax in baf_axd.items():
                baf_segment_gdf = baf_segment_gdf_dict[baf_idx]
                assert baf_segment_gdf is not None
                #if baf_segment_gdf is not None:
                gdraw_result = baf_segment_gdf.draw_baf(
                    ax=baf_ax, 
                    genomeplotter=genomeplotter,
                    setup_axes=setup_axes,
                    title=None,
                    ylabel_prefix=ylabel_prefix,
                    nproc=nproc,
                    verbose=verbose,
                    draw_common_kwargs=draw_common_kwargs,
                    **baf_segment_kwargs,
                )
                gdraw_axresult_list.append(gdraw_result)

            # segment - corrected baf
            if CN_exists:
                for baf_idx, baf_ax in baf_axd.items():
                    legend_ax = (
                        axd[baf_idx + '_legend']
                        if with_legend else
                        None
                    )
                    gdraw_result = CN_gdf.draw_corrected_baf(
                        baf_index=baf_idx,
                        ax=baf_ax,
                        legend_ax=legend_ax,
                        genomeplotter=genomeplotter,
                        plotdata=CN_plotdata,

                        nproc=nproc,
                        verbose=verbose,
                    )
                    gdraw_axresult_list.append(gdraw_result)

            # baf rawdata
            if draw_baf_rawdata:
                #baf_rawdata_gdf = self.gdfs['baf_rawdata']
                if filtered_baf:
                    baf_rawdata_gdf = self.get_filtered_baf_rawdata(corrector_id=corrector_id)
                else:
                    baf_rawdata_gdf = self.get_baf_rawdata(corrector_id=corrector_id)

                if baf_rawdata_gdf is not None:
                    plotdata_data = (
                        baf_rawdata_gdf
                        if frac is None else
                        baf_rawdata_gdf.sample(frac=frac)
                    )
                    baf_rawdata_plotdata = genomeplotter.make_plotdata(
                        plotdata_data, 
                        log_suffix=' (BAF raw data)', 
                        nproc=nproc,
                        verbose=verbose,
                    )
                    for baf_index, baf_ax in baf_axd.items():
                        gdraw_result = baf_rawdata_gdf.draw_baf(
                            baf_index, 
                            ax=baf_ax,
                            genomeplotter=genomeplotter,
                            plotdata=baf_rawdata_plotdata,
                            ylabel_prefix=ylabel_prefix,
                            setup_axes=setup_axes,
                            verbose=verbose,
                            **baf_rawdata_kwargs,
                        )
                        gdraw_axresult_list.append(gdraw_result)

        # prepare plotdata for depth and MQ
        if draw_depth or draw_MQ:
            use_corrected = (corrector_id is not None)

            if depthtype is None:
                #depthtype = ('corr' if use_corrected else 'norm')
                depthtype = 'norm'

            if use_corrected:
                excl_region = self.get_corrected_depth_excluded_region(corrector_id=corrector_id)
            else:
                excl_region = GDF.init_empty(refver=self.refver)

            if draw_depth_rawdata:
                depth_rawdata = self.get_depth_rawdata(corrector_id=corrector_id)
                plotdata_data = (
                    depth_rawdata
                    if frac is None else
                    depth_rawdata.sample(frac=frac)
                )
                depth_raw_plotdata = genomeplotter.make_plotdata(
                    plotdata_data, 
                    log_suffix=' (Depth raw data)', 
                    nproc=nproc,
                    verbose=verbose,
                )

            depth_segment = self.get_depth_segment(corrector_id=corrector_id)
            depth_segment_exists = (depth_segment is not None)
            if depth_segment_exists:
                depth_seg_plotdata = genomeplotter.make_plotdata(
                    depth_segment, 
                    log_suffix=' (Depth segment)', 
                    nproc=nproc,
                    verbose=verbose,
                )

        # draw depth
        if draw_depth:
            default_ylabel, y_colname = depth_segment.get_colname_ylabel(depthtype)

            if depth_segment_exists:
                gdraw_result = depth_segment.draw_depth(
                    ax=axd['depth'],
                    genomeplotter=genomeplotter,
                    depthtype=depthtype,
                    plotdata=depth_seg_plotdata,
                    chromwise_peaks=chromwise_peaks,
                    chromwise_peaks_density_kwargs=chromwise_peaks_density_kwargs,
                    chromwise_peaks_line_kwargs=chromwise_peaks_line_kwargs,
                    setup_axes=setup_axes,
                    draw_common_kwargs=draw_common_kwargs,
                    ylabel_prefix=ylabel_prefix,
                    verbose=verbose,
                    **depth_segment_kwargs,
                )
                gdraw_axresult_list.append(gdraw_result)

            if draw_depth_rawdata:
                gdraw_result = depth_rawdata.draw_depth(
                    ax=axd['depth'],
                    genomeplotter=genomeplotter,
                    depthtype=depthtype,
                    plotdata=depth_raw_plotdata,
                    setup_axes=(not depth_segment_exists),
                    ylabel_prefix=ylabel_prefix,
                    verbose=verbose,
                    **depth_rawdata_kwargs,
                )
                gdraw_axresult_list.append(gdraw_result)

        # draw MQ
        if draw_MQ:
            if draw_MQ_rawdata:
                gdraw_result = depth_rawdata.draw_MQ(
                    ax=axd['MQ'],
                    genomeplotter=genomeplotter,
                    plotdata=depth_raw_plotdata,
                    setup_axes=(not depth_segment_exists),
                    ylabel_prefix=ylabel_prefix,
                    verbose=verbose,
                    **MQ_kwargs,
                )
                gdraw_axresult_list.append(gdraw_result)
            if depth_segment_exists:
                gdraw_result = depth_segment.draw_MQ(
                    ax=axd['MQ'],
                    genomeplotter=genomeplotter,
                    plotdata=depth_seg_plotdata,
                    setup_axes=setup_axes,
                    draw_common_kwargs=draw_common_kwargs,
                    ylabel_prefix=ylabel_prefix,
                    verbose=verbose,
                    **MQ_kwargs,
                )
                gdraw_axresult_list.append(gdraw_result)

        # draw mutation
        if (draw_mut and (vplist is not None)):
            #if vcf_sampleid is None:
            #    vcf_sampleid = self.sampleid
            vp_gdf = vplist.get_vafdf(self.sampleid, add_depth=True, exclude_other=False)
            vp_gdf.draw_vaf(
                ax=axd['mut'],
                genomeplotter=genomeplotter,
                setup_axes=setup_axes,
                draw_common_kwargs=draw_common_kwargs,
                ylabel_prefix=ylabel_prefix,
                verbose=verbose,
            )

        # draw excluded
        if draw_excluded and (corrector_id is not None):
            excluded_region = self.get_corrected_depth_excluded_region(corrector_id)
            excluded_region = excluded_region.assign(colors='red')
            excluded_region_plotdata = genomeplotter.make_plotdata(
                excluded_region, 
                log_suffix=' (excluded region)', 
                nproc=nproc,
                verbose=verbose,
            )
            excluded_kwargs = (
                dict(
                    alpha=0.01,
                ) | excluded_kwargs
            )
            for ax in [val for key, val in axd.items() if not key.startswith('blank')]:
                excluded_region.draw_boxes(
                    ax=ax, 
                    color_colname='colors', 
                    plot_kwargs=excluded_kwargs,
                    plotdata=excluded_region_plotdata,
                    setup_axes=setup_axes,
                )

        draw_result = CNVDrawingResult(
            cnvsample=self,
            corrector_id=corrector_id,
            fig=fig, 
            axd=axd,
            genomeplotter=genomeplotter,
            vplist=vplist,
        )
        #self.draw_result_corrector_id = corrector_id

        ###

        selected_sol = self.get_selected_solution(corrector_id=corrector_id)
        if selected_sol is None:
            title = self.make_title(selected_solution=None, omit_mode=omit_title_mode)
            plotmisc.draw_suptitle(draw_result.fig, title, **suptitle_kwargs)
        else:
            try:
                if draw_CN:
                    draw_result.draw_solution(
                        index=selected_sol['index'], 
                        verbose=verbose, 
                        omit_title_mode=omit_title_mode,
                        **CN_kwargs,
                    )
            except:
                raise

        return draw_result

    #def disconnect_plot(self):
    #    self.draw_result.disconnect()

    ########################
    # LOW level ax drawers #
    ########################

    def draw_data_ax(self):
        pass

    def draw_solution_ax(self):
        pass

    ######################
    # ax drawing helpers #
    ######################

    def make_title(self, selected_solution=None, omit_mode=False):
        basic_keys = (
            ['sampleid', 'is_female']
            if omit_mode else
            ['sampleid', 'is_female', 'mode']
        )
        result = ', '.join(
            f'{key}={getattr(self, key)}' 
            for key in basic_keys
        )
        if selected_solution is not None:
            cellularity = selected_solution['cellularity']
            ploidy = selected_solution['ploidy']
            result = (
                result 
                + '\n' 
                + f'cellularity={cellularity:.3}, ploidy={ploidy:.3}'
            )
        return result

    def get_baf_indexes(self):
        return libbaf.get_baf_indexes(self.simple_attrs["ploidy"])

    ####################
    # solution related #
    ####################

#    def make_mean_values(self):
#        result = self.draw_result.get_mean_values()
#        depth_rawdata = self.get_depth_rawdata(corrector_id=self.draw_result_corrector_id)
#        result['raw_depth_rawdata'] = genomedf_draw.get_region_mean_values(
#            plotdata=depth_rawdata,
#            y_colname=depth_rawdata.__class__.depth_colname, 
#            region_gdfs=self.draw_result.get_region_gdfs(),
#        )
#        result['corrected_baf_rawdata'] = bafsimul.predict_true_baf(
#            np.stack(
#                [result['baf_rawdata'], result['raw_depth_rawdata']],
#                axis=1,
#            )
#        )
#
#        return result
#
#    def set_mean_values(self):
#        self.mean_values = self.make_mean_values()
#
#    def get_mean_values(self):
#        return self.mean_values

    def get_solution_candidates(self, corrector_id):
        return self.solution_attrs[corrector_id]['solution_candidates']

    def set_selected_solution(self, corrector_id, index):
        self.solution_attrs[corrector_id]['selected_solution'] = self.pick_solution(index, corrector_id)

    def get_selected_solution(self, corrector_id):
        if corrector_id in self.solution_attrs:
            return self.solution_attrs[corrector_id]['selected_solution']
        else:
            return None

    def pick_solution(self, index, corrector_id):
        index = tuple(index)
        result = dict()
        sol_cands = self.get_solution_candidates(corrector_id=corrector_id)
        for key in ('CNt0', 'onecopy_depth', 'num_division', 'cellularity', 'K', 'targetval', 'ploidy'):
            result[key] = sol_cands[key][index]
        for key in ('CNt', 'Bt'):
            result[key] = sol_cands[key][np.index_exp[:,] + index]
        result['index'] = index

        return result

    def list_bestfit_solutions(self, q, corrector_id, allvalue=True):
        sol_cands = self.get_solution_candidates(corrector_id=corrector_id)

        #if allvalue:
        #    targetvals = sol_cands['all_targetval']
        #else:
        targetvals = sol_cands['targetval']

        quantile = np.nanquantile(targetvals, q)
        indexes = np.nonzero(targetvals < quantile)
        CNt0s = sol_cands['CNt0'][indexes]
        ndivs = sol_cands['num_division'][indexes]
        return indexes, CNt0s, ndivs

    def show_solution(self, index):
        mean_values = self.get_mean_values()
        nseg = len(mean_values['depths'])
        solution = self.pick_solution(index)

        fig, axd = plt.subplot_mosaic(
            [['CN'], ['baf'], ['depth']], 
            gridspec_kw=dict(hspace=0.4),
            figsize=(8, 10),
        )

        xs = np.arange(nseg)
        width = 0.8

        axd['CN'].set_xlim(-1, nseg)
        axd['CN'].hlines(solution['CNt'], xs - width/2, xs + width/2, alpha=0.4, color='black')
        axd['CN'].hlines(solution['Bt'], xs - width/2, xs + width/2, alpha=0.4, color='blue')
        axd['CN'].set_title('CN')

        axd['baf'].set_xlim(-1, nseg)
        axd['baf'].set_ylim(-0.02, 0.52)
        axd['baf'].hlines(mean_values['baf_rawdata'], xs - width/2, xs + width/2, alpha=0.4)
        axd['baf'].hlines(solution['predicted_baf'], xs - width/2, xs + width/2, alpha=0.4, color='tab:red')
        axd['baf'].set_title('BAF')

        axd['depth'].set_xlim(-1, nseg)
        axd['depth'].set_ylim(0, np.floor(mean_values['depth_rawdata'].max()) + 1)
        axd['depth'].hlines(mean_values['depth_rawdata'], xs - width/2, xs + width/2, alpha=0.4)
        axd['depth'].set_title('depth')

        depth0 = mean_values['depth_rawdata'].min()
        for y in depth0 + solution['onecopy_depth'] * np.arange(20):
            axd['depth'].axhline(y, color='tab:red', alpha=0.4, linestyle='dotted')

        return fig, axd



#########################################################################
#########################################################################
#########################################################################


class CNVManager:
    @classmethod
    def init_empty(cls, refver, nproc=os.cpu_count()):
        result = cls()
        # simple attrs
        result.simple_attrs = {
            'refver': refgenome.standardize(refver),
            'nproc': nproc,
        }

        # data dicts
        result.init_datadicts()

        return result

    def init_datadicts(self):
        self.cnvsamples = dict()
        #self.depth_correctors = dict()
        self.correctors = dict()

    def __getattr__(self, key):
        return getattr_base(self, key)

#    def __getattr__(self, key):
#        if key in self.simple_attrs.keys():
#            return self.simple_attrs[key]
#        else:
#            try:
#                super().__getattr__(key)
#            except AttributeError:
#                raise AttributeError(f'{repr(self.__class__)} object has no attribute {repr(key)}')

    #####

    @property
    def samples(self):
        return self.cnvsamples

    def get_cnvsample_dir(self, topdir):
        return os.path.join(topdir, 'cnvsample')

    def add_cnvsample(self, cnvsample, force=False):
        if not force:
            assert cnvsample.sampleid not in self.cnvsamples
        self.cnvsamples[cnvsample.sampleid] = cnvsample

    def load_cnvsample_pickle(
        self, cnvsample_topdir,
        load_all=True,
        load_rawdata=False,
        load_segment=False,
        load_corrected=False,
    ):
        logutils.log(f'Loading CNVSample from directory {cnvsample_topdir}')
        cnvsample = CNVSample.load_pickle_asdf(
            cnvsample_topdir,
            load_all=load_all,
            load_rawdata=load_rawdata,
            load_segment=load_segment,
            load_corrected=load_corrected,
        )
        self.add_cnvsample(cnvsample)

    def save_cnvsample(self, cnvsample_id, all_topdir, force=False):
        logutils.log(f'Saving CNVSample {cnvsample_id}')
        cnvsample = self.cnvsamples[cnvsample_id]

        cnvsample_topdir = os.path.join(self.get_cnvsample_dir(all_topdir), cnvsample_id)
        cnvsample.save_pickle_asdf(cnvsample_topdir, force=force)

    #####

#    def save_depth_corrector(self, corrector_id, savepath):
#        logutils.log(f'Saving depth corrector {corrector_id}')
#        saved_dict = self.depth_correctors[corrector_id].copy()
#
#        saved_dict['depth_rawdata'] = saved_dict['depth_rawdata'].df
#        if saved_dict['raw_target_region'] is not None:
#            saved_dict['raw_target_region'] = saved_dict['raw_target_region'].df
#
#        with open(savepath, mode='wb') as outfile:
#            pickle.dump(saved_dict, outfile)
#
#    def save_all_depth_correctors(self, topdir):
#        depthcorrector_dir = os.path.join(topdir, 'depth_correctors')
#        os.makedirs(depthcorrector_dir, exist_ok=True)
#        for corrector_id in self.depth_correctors.keys():
#            savepath = os.path.join(depthcorrector_dir, corrector_id + '.pickle')
#            self.save_depth_corrector(corrector_id, savepath)
#
#    def load_depth_corrector(self, savepath):
#        logutils.log(f'Loading depth corrector from directory {savepath}')
#        with open(savepath, 'rb') as infile:
#            corrector_dict = pickle.load(infile)
#
#        corrector_dict['depth_rawdata'] = DepthRawDF.from_frame(corrector_dict['depth_rawdata'], refver=self.refver)
#        if corrector_dict['raw_target_region'] is not None:
#            corrector_dict['raw_target_region'] = GDF.from_frame(corrector_dict['raw_target_region'], refver=self.refver)
#
#        self.depth_correctors[corrector_dict['id']] = corrector_dict

    #####

    def get_corrector_dir(self, topdir):
        return os.path.join(topdir, 'correctors')

    def save_corrector(self, corrector_id, savepath):
        logutils.log(f'Saving corrector {corrector_id}')
        saved_dict = self.correctors[corrector_id].copy()

        for key in (
            'raw_target_region', 
            'depth_rawdata', 
            'vaf_rawdata', 
            #'vaf_rawdata_nomerge',
        ):
            if key == 'raw_target_region':
                if saved_dict[key] is None:
                    continue
                else:
                    saved_dict[key] = saved_dict[key].df

            elif key == 'vaf_rawdata':
                if key in saved_dict:
                    saved_dict[key] = saved_dict[key].df

            elif key == 'depth_rawdata':
                saved_dict[key] = saved_dict[key].df

        with open(savepath, mode='wb') as outfile:
            pickle.dump(saved_dict, outfile)

    def save_all_correctors(self, topdir):
        corrector_dir = self.get_corrector_dir(topdir)
        os.makedirs(corrector_dir, exist_ok=True)
        for corrector_id in self.correctors.keys():
            savepath = os.path.join(corrector_dir, corrector_id + '.pickle')
            self.save_corrector(corrector_id, savepath)

    def load_corrector(self, savepath):
        logutils.log(f'Loading corrector from directory {savepath}')
        with open(savepath, 'rb') as infile:
            corrector_dict = pickle.load(infile)

        corrector_dict['depth_rawdata'] = DepthRawDF.from_frame(
            corrector_dict['depth_rawdata'], refver=self.refver,
        )

        if 'vaf_rawdata' in corrector_dict:
            corrector_dict['vaf_rawdata'] = VariantDataFrame.from_frame(
                corrector_dict['vaf_rawdata'], refver=self.refver,
            )

        #if 'vaf_rawdata_nomerge' in corrector_dict:
        #    corrector_dict['vaf_rawdata_nomerge'] = VariantDataFrame.from_frame(
        #        corrector_dict['vaf_rawdata_nomerge'], refver=self.refver,
        #    )
        if corrector_dict['raw_target_region'] is not None:
            corrector_dict['raw_target_region'] = GDF.from_frame(
                corrector_dict['raw_target_region'], refver=self.refver,
            )

        self.correctors[corrector_dict['id']] = corrector_dict

    #####

    def get_simple_attrs_path(self, topdir):
        return os.path.join(topdir, 'simple_attrs.json')

    def save_pickle(self, topdir, save_correctors=True, sampleids=None):
        os.makedirs(topdir, exist_ok=True)

        # simple attrs
        with open(self.get_simple_attrs_path(topdir), 'wt') as outfile:
            json.dump(self.simple_attrs, outfile)

        # save cnvsamples
        os.makedirs(self.get_cnvsample_dir(topdir), exist_ok=True)
        if sampleids is not None:
            assert set(sampleids).issubset(self.cnvsamples.keys())
            sids_to_save = sampleids
        else:
            sids_to_save = tuple(self.cnvsamples.keys())
        for sid in sids_to_save:
            self.save_cnvsample(sid, topdir, force=True)

        # save depth correctors
        if save_correctors:
            self.save_all_correctors(topdir)

    @classmethod
    def load_pickle(
        cls, topdir, 
        cnvsamples=None, 
        load_corrector=False,
        correctors=None,

        ###

        load_all=True,
        load_rawdata=False,
        load_segment=False,
        load_corrected=False,
    ):
        result = cls()
        result.init_datadicts()

        # simple attrs
        with open(result.get_simple_attrs_path(topdir), 'rt') as infile:
            result.simple_attrs = json.load(infile)

        # cnvsamples
        logutils.log(f'Beginning loading CNVSamples')
        cnvsample_dir = result.get_cnvsample_dir(topdir)

        all_fnames = sorted(os.listdir(cnvsample_dir))
        if cnvsamples is None:
            cnvsamples = all_fnames
        else:
            cnvsamples = list(cnvsamples)
        if not set(cnvsamples).issubset(all_fnames):
            raise Exception(
                f'Some of the specified CNVSample IDs are not present. '
                f'Existing CNVSample IDs: {all_fnames}'
            )

        for fname in cnvsamples:
            result.load_cnvsample_pickle(
                os.path.join(cnvsample_dir, fname),
                load_all=load_all,
                load_rawdata=load_rawdata,
                load_segment=load_segment,
                load_corrected=load_corrected,
            )

        logutils.log(f'Finished loading CNVSamples')

        # correctors
        if load_corrector:
            logutils.log(f'Beginning loading correctors')
            corrector_dir = result.get_corrector_dir(topdir)

            for fname in sorted(os.listdir(corrector_dir)):

                corrector_id = re.sub(r'\.pickle$', '', fname)
                if correctors is not None:
                    if corrector_id not in correctors:
                        continue

                result.load_corrector(os.path.join(corrector_dir, fname))

            logutils.log(f'Finished loading correctors')

        return result

    ##############################
    # depth corrector generation #
    ##############################

    @staticmethod
    def make_corrector_sanitycheck(cnvsample_list):
        # refver
        refvers = set(refgenome.standardize(x.refver) for x in cnvsample_list)
        if len(refvers) != 1:
            raise Exception(f'Reference versions differ between CNVSample objects.')
        refver = refvers.pop()

        # ploidy
        ploidies = set(x.ploidy for x in cnvsample_list)
        if len(ploidies) != 1:
            raise Exception(f'Ploidy values differ between CNVSample objects.')
        ploidy = ploidies.pop()

        # mode
        modes = set(x.mode for x in cnvsample_list)
        if len(modes) != 1:
            raise Exception(f'Mode values differ between CNVSample objects.')
        mode = modes.pop()

        # raw target region
        raw_target_list = [x.get_raw_target_region() for x in cnvsample_list]
        isnone = set((x is None) for x in raw_target_list)
        if len(isnone) != 1:
            raise Exception(f'Some of raw_target_region are None at some are not.')
        isnone = isnone.pop()
        if not isnone:
            all_same = all((x == y) for (x, y) in tools.pairwise(raw_target_list))
            if not all_same:
                raise Exception(f'Raw target regions differ between CNVsample objects.')
        raw_target_region = raw_target_list[0]

        return refver, ploidy, mode, raw_target_region

    @staticmethod
    def make_corrector_depth_onesample(cnvsample, nproc, apply_filter=False):
        if cnvsample.mode == 'wgs':
            depth_rawdata_gdf = cnvsample.get_depth_rawdata()
            depth_segment = cnvsample.get_depth_segment()
            merged_segment = cnvsample.get_merged_segment()
            assert merged_segment.check_has_CN()

            # 1-1. add calculated germline CN to depth rawdata
            depth_rawdata_gdf = depth_rawdata_gdf.join(
                merged_segment,
                how='left',
                right_gdf_cols=CNVSegDF.clonal_CN_colname,
                merge='longest',
                overlapping_length=True,
                suffixes={'longest': ''},
                omit_N=True,
                find_nearest=True,
                nproc=nproc,
            )

            # 1-2. add segment-based filters
            if apply_filter:
                depth_rawdata_gdf = depth_rawdata_gdf.join(
                    depth_segment,
                    how='left',
                    right_gdf_cols=DepthSegDF.filter_colname,
                    merge='longest',
                    overlapping_length=True,
                    suffixes={'longest': ''},
                    omit_N=True,
                    find_nearest=True,
                    nproc=nproc,
                )

            depth_rawdata_gdf.sort()

            # 2. prepare parameters
            CN = depth_rawdata_gdf[CNVSegDF.clonal_CN_colname].astype(float)
            raw_depth = depth_rawdata_gdf.raw_depth
            lengths = depth_rawdata_gdf.get_lengths()
            if apply_filter:
                filters = depth_rawdata_gdf[DepthSegDF.filter_colname]
            else:
                filters = np.full(CN.shape, True)

        elif cnvsample.mode == 'panel':
            # 1. add default CNg to depth rawdata
            cnvsample.add_default_CNg_Bg_to_depth_rawdata(nproc)
            depth_rawdata_gdf = cnvsample.get_depth_rawdata()
            depth_rawdata_gdf.sort()

            # 2. prepare parameters
            CN = depth_rawdata_gdf[cnvcall.DEFAULT_CNG_COLNAME].astype(float)
            raw_depth = depth_rawdata_gdf.raw_depth
            lengths = depth_rawdata_gdf.get_lengths()
            if apply_filter:
                filters = depth_rawdata_gdf.get_filter()
            else:
                filters = np.full(CN.shape, True)

        # 3. set CN == 0 regions as invalid
        CNzero_selector = (CN == 0)
        excl_selector = np.logical_or(CNzero_selector, ~filters)
        CN[CNzero_selector] = np.nan
        raw_depth[excl_selector] = np.nan

        # 4. correct raw depth by dividing with CN
        corrected_raw_depth = raw_depth * (cnvsample.ploidy / CN)  # where CN == 0, divided by np.nan

        # 5. make normalized depth from corrected_raw_depth
        corrector_depth = libdepth.make_normalized_depth(corrected_raw_depth, lengths)

        # 6. make depth gdf coordinate columns
        #depth_rawdata_coords = depth_rawdata_gdf.get_coordinate_array()

        return corrector_depth, depth_rawdata_gdf

    def make_depth_corrector(self, cnvsample_list, nproc):
        # 1. merge excluded regions
#        logutils.log(f'Beginning merging excluded regions')
#        excl_region_list = list()
#        for cnvsample in cnvsample_list:
#            excl_region_list.append(cnvsample.get_excluded_region())
#        merged_excl_region = libgdf.union(excl_region_list)
#        logutils.log(f'Finished merging excluded regions')

        # 1. make corrector depth for each cnvsample
        corrector_depth_list = list()
        depth_rawdata_gdf_list = list()
        for cnvsample in cnvsample_list:
            logutils.log(f'Making corrector depth for CNVSample {cnvsample.sampleid}')
            corrector_depth, depth_rawdata_gdf = self.make_corrector_depth_onesample(cnvsample, nproc)
            corrector_depth_list.append(corrector_depth)
            depth_rawdata_gdf_list.append(depth_rawdata_gdf)

        # 2. sanitycheck - assure all depth rawdata gdfs have identical coords
        logutils.log(f'Confirming identity of depth rawdata coordinates')
        coords_array_list = [x.get_coordinate_array() for x in depth_rawdata_gdf_list]
        all_same_coords = all((x == y).all() for (x, y) in tools.pairwise(coords_array_list))
        if not all_same_coords:
            raise Exception(f'Coordinates differ between depth rawdata of different CNVSamples')

        # 3. make average of all corrector depths
        logutils.log(f'Calculating average of corrector depths')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            avg_corrector_depth = np.nanmean(np.stack(corrector_depth_list, axis=1), axis=1)
                # in the stacked array, axis 0 means genomic positions, axis 1 means samples
                # invalid values are replaced with np.nan
                # averages are calculated ignoring invalid slots (useful for e.g. Y chromosome or other sample-specific 0-copy regions)
            
        # 4. make corrector depth gdf
        depth_corrector_gdf = depth_rawdata_gdf_list[0].drop_annots()
        depth_corrector_gdf[DepthRawDF.norm_depth_colname] = avg_corrector_depth

        return depth_corrector_gdf

    ############################
    # baf corrector generation #
    ############################

    def modify_with_predicted_baf(self, baf_gdf_copy, cnvsample, nproc):
        # make array of predicted bafs
        if cnvsample.check_nobias():
            logutils.log(f'Adding predicted baf values to baf rawdata gdf')
            germline_CN_gdf = cnvsample.get_merged_segment()
            for baf_index in cnvsample.get_baf_indexes():
                germline_CN_gdf.add_predicted_baf_germline(baf_index)

            added_cols = [
                germline_CN_gdf.get_predicted_baf_colname(baf_index)
                for baf_index in cnvsample.get_baf_indexes()
            ]
            baf_gdf_copy = baf_gdf_copy.join(
                germline_CN_gdf,
                how='left',
                right_gdf_cols=added_cols,
                merge='longest',
                overlapping_length=True,
                omit_N=True,
                suffixes={'longest': ''},
                nproc=nproc,
            )
            
            predicted_baf_list = [
                baf_gdf_copy[colname] for colname in added_cols
            ]
            predicted_bafs = np.stack(predicted_baf_list, axis=1)
        else:
            predicted_baf_list = list()
            for baf_index in cnvsample.get_baf_indexes():
                predicted_baf = cnvcall.get_predicted_baf(
                    CNt=cnvcall.get_default_CNg(baf_gdf_copy, cnvsample.is_female, nproc),
                    Bt=cnvcall.get_default_Bg(baf_gdf_copy, cnvsample.is_female, nproc),
                    cellularity=1,
                )
                predicted_baf_list.append(predicted_baf)
            predicted_bafs = np.stack(predicted_baf_list, axis=1)
        
        lastcol = 1 - predicted_bafs.sum(axis=1)
        predicted_vafs = np.concatenate([predicted_bafs, lastcol[:, np.newaxis]], axis=1)

        # apply correction to VAF values by dividing with corresponding predicted bafs
        logutils.log(f'Dividing raw VAF values with predicted values')
        vaf_values = baf_gdf_copy[baf_gdf_copy.vaf_columns]
        orders = vaf_values.argsort(axis=1).argsort(axis=1)
        predicted_vafs_reordered = np.stack(
            [predicted_vafs[idx, orders[idx, :]] for idx in range(predicted_vafs.shape[0])],
            axis=0,
        )
        corrected_vaf_values = vaf_values / predicted_vafs_reordered
        baf_gdf_copy.loc[:, baf_gdf_copy.vaf_columns] = corrected_vaf_values

        return baf_gdf_copy

    @deco_nproc
    def make_modified_vaf_gdf(self, cnvsample, divide_with_REF=True, nproc=0):
        assert cnvsample.is_germline

        logutils.log(f'Making copy of baf rawdata')
        if divide_with_REF:
            baf_gdf_copy = cnvsample.get_baf_rawdata().copy()
        else:
            baf_gdf_copy = cnvsample.get_baf_rawdata()

        logutils.log(f'Adding hetalt flag')
        baf_gdf_copy.add_baf_hetalt_flag()

        logutils.log(f'Removing invalid baf data')
        if cnvsample.check_nobias():
            baf_gdf_copy = baf_gdf_copy.filter_valid_bafs(
                germline_CN_gdf=cnvsample.get_merged_segment(), 
                remove_loh=False,
                nproc=nproc,
                verbose=False,
            )
        else:
            baf_gdf_copy = baf_gdf_copy.filter_valid_bafs(
                is_female=cnvsample.is_female, 
                remove_loh=False,
                nproc=nproc,
                verbose=False,
            )

        # divide with predicted baf
        logutils.log(f'Removing regions with imbalanced allele copy numbers')
        #baf_gdf_copy = self.modify_with_predicted_baf(baf_gdf_copy, cnvsample, nproc)
        if cnvsample.check_nobias():
            baf_gdf_copy = baf_gdf_copy.keep_equal_germline_allelecopy_region(
                germline_CN_gdf=cnvsample.get_merged_segment(), 
                nproc=nproc,
            )
        else:
            baf_gdf_copy = baf_gdf_copy.keep_equal_germline_allelecopy_region(
                is_female=cnvsample.is_female, 
                nproc=nproc,
            )

        # make raw VAF values into odds ratios (by dividing with REF_vaf)
        if divide_with_REF:
            logutils.log(f'Making VAFs into odds ratios')
            baf_gdf_copy.divide_with_REF_vaf()

        # final
        baf_gdf_copy = baf_gdf_copy.choose_annots(baf_gdf_copy.vaf_columns)
        if divide_with_REF:
            baf_gdf_copy.add_suffix(f'_{cnvsample.sampleid}', inplace=True)
        else:
            baf_gdf_copy = baf_gdf_copy.add_suffix(f'_{cnvsample.sampleid}', inplace=False)

        return baf_gdf_copy

    @deco_nproc
    @get_deco_logging(f'making vaf corrector')
    def make_vaf_corrector(self, cnvsample_list, nomerge=True, nproc=0):
        modified_vafdfs = list()
        for cnvsample in cnvsample_list:
            logutils.log(f'Processing {cnvsample.sampleid}')
            modified_vafdfs.append(
                self.make_modified_vaf_gdf(cnvsample, divide_with_REF=(not nomerge), nproc=nproc)
            )

        logutils.log(f'Merging vafdfs')
        merged_vaf_gdf = VariantDataFrame.variant_union(modified_vafdfs, nproc=nproc)
        if nomerge:
            vaf_corrector_gdf = merged_vaf_gdf
        else:
            subdf_arrays = dict()
            for prefix in merged_vaf_gdf.allele_columns:
                selected_cols = [x for x in merged_vaf_gdf.annot_columns if x.startswith(prefix)]
                subdf_arrays[prefix] = merged_vaf_gdf[selected_cols]

            vaf_corrector_gdf = merged_vaf_gdf.drop_annots()
            vaf_corrector_gdf['n_samples'] = (
                np.logical_not(
                    np.isnan(subdf_arrays['REF'])
                )
                .sum(axis=1)
            )
            for prefix, arr in subdf_arrays.items():
                vaf_corrector_gdf[f'{prefix}_vaf_corrector'] = np.nanmean(arr, axis=1)

        return vaf_corrector_gdf

    @deco_nproc
    @deco.get_deco_atleast1d(['cnvsample_idlist'])
    #def make_depth_corrector(self, sid, cnvsample_idlist, nproc=0, force=False):
    def make_corrector(self, corrector_id, cnvsample_idlist, nproc=0, force=False):
        # sanitycheck
        if not force:
            assert corrector_id not in self.correctors
        cnvsample_list = [self.cnvsamples[x] for x in cnvsample_idlist]
        refver, ploidy, mode, raw_target_region = self.make_corrector_sanitycheck(cnvsample_list)
        
        # depth corrector
        depth_corrector_gdf = self.make_depth_corrector(cnvsample_list, nproc)

        # vaf corrector
        vaf_corrector_gdf = self.make_vaf_corrector(cnvsample_list, nproc=nproc, nomerge=True)

        # save to dict
        self.correctors[corrector_id] = {
            'id': corrector_id,
            'cnvsample_ids': tuple(cnvsample_idlist),
            #'excluded_region': merged_excl_region,
            'mode': mode,
            'refver': refver,
            'ploidy': ploidy,

            'raw_target_region': raw_target_region,
            'depth_rawdata': depth_corrector_gdf,
            'vaf_rawdata': vaf_corrector_gdf,
            #'vaf_rawdata_nomerge': merged_vaf_gdf,
        }

    ####################################
    # add corrected depth to cnvsample #
    ####################################

    @deco_nproc
    def make_corrected_gdfs(
        self, 
        *,
        corrector_id, 
        cnvsample_id,
        paired_germline_id=None,
        do_vaf_correction=False,

        skip_hetalt_filter=False,  # works only when "do_vaf_correction" is False ; indicates not performing hetalt filtering with tumor-only sample
        force=False,
        nproc=0,

        segment_kwargs=dict(),
        winsorize=libdepth.DEFAULT_WINSORIZE,
    ):
        cnvsample = self.samples[cnvsample_id]
        corrector_dict = self.correctors[corrector_id]

        # sanitycheck
        if not force:
            assert corrector_id not in cnvsample.corrected_gdfs.keys()
        assert corrector_dict['mode'] == cnvsample.mode
        assert refgenome.compare_refvers(corrector_dict['refver'], cnvsample.refver)
        assert corrector_dict['ploidy'] == cnvsample.ploidy

        logutils.log(f'Beginning corrected data generation (cnvsample id: {cnvsample_id}, corrector id: {corrector_id})')

        # init entity
        cnvsample.make_corrected_gdfs_init(corrector_id)

        # make corrected depth
        MQ_cutoff = (45, None)
        norm_depth_cutoff = None
        if cnvsample.mode == 'wgs':
            raw_depth_cutoff = (0, None) 
        elif cnvsample.mode == 'panel':
            raw_depth_cutoff = (50, None) 

        cnvsample.make_corrected_gdfs_depth(
            corrector_id=corrector_id, 
            corrector_gdf=corrector_dict['depth_rawdata'],
            MQ_cutoff=MQ_cutoff, 
            norm_depth_cutoff=norm_depth_cutoff,
            raw_depth_cutoff=raw_depth_cutoff,
            nproc=nproc,
            segment_kwargs=segment_kwargs,
            winsorize=winsorize,
        )

        # make corrected vaf and baf
        if paired_germline_id is None:
            germline_CN_gdf = None
            germline_baf_rawdata = None
        else:
            germline_CN_gdf = self.samples[paired_germline_id].get_merged_segment()
            germline_baf_rawdata = self.samples[paired_germline_id].get_baf_rawdata()

        if do_vaf_correction:
            cnvsample.make_corrected_gdfs_baf(
                corrector_id=corrector_id, 
                corrector_gdf=corrector_dict['vaf_rawdata'], 
                germline_CN_gdf=germline_CN_gdf,
                germline_baf_rawdata=germline_baf_rawdata,
                nproc=nproc,
                segment_kwargs=segment_kwargs,
            )
        else:
            baf_segment_dict, filtered_baf_rawdata = cnvsample.get_baf_rawdata().get_segment_dict(
                target_region=cnvsample.get_raw_target_region(),

                germline_baf_rawdata=germline_baf_rawdata,

                germline_CN_gdf=germline_CN_gdf, 
                is_female=cnvsample.is_female,
                skip_hetalt_filter=skip_hetalt_filter,

                return_filtered_rawdata=True,

                nproc=nproc, 
                **segment_kwargs,
            )
            cnvsample.corrected_gdfs[corrector_id]['baf_segment'] = baf_segment_dict
            cnvsample.corrected_gdfs[corrector_id]['filtered_baf_rawdata'] = filtered_baf_rawdata
            cnvsample.corrected_gdfs[corrector_id]['baf_rawdata'] = cnvsample.get_baf_rawdata()

        # make merged segment
        cnvsample.make_corrected_gdfs_merged_segment(
            corrector_id=corrector_id,
            germline_CN_gdf=germline_CN_gdf,
            nproc=nproc,
        )

    ###########
    # drawing #
    ###########

    def make_axd(
        self,

        tumor_cnvsample,
        normal_cnvsample,

        draw_tumor_all,
        draw_tumor_depth,
        draw_tumor_baf,
        draw_tumor_CN,
        draw_tumor_MQ,

        draw_normal_all,
        draw_normal_depth,
        draw_normal_baf,
        draw_normal_CN,
        draw_normal_MQ,

        subplots_kwargs,
        figsize,
    ):
        # make mosaic
        tumor_mosaic = tumor_cnvsample.make_mosaic(
            draw_depth=draw_tumor_depth,
            draw_baf=draw_tumor_baf,
            draw_MQ=draw_tumor_MQ,
            draw_CN=draw_tumor_CN,
            draw_depth_peaks=False,
        )
        tumor_mosaic = [
            ['tumor_' + y for y in x]
            for x in tumor_mosaic
        ]

        normal_mosaic = normal_cnvsample.make_mosaic(
            draw_depth=draw_normal_depth,
            draw_baf=draw_normal_baf,
            draw_MQ=draw_normal_MQ,
            draw_CN=draw_normal_CN,
            draw_depth_peaks=False,
        )
        normal_mosaic = [
            ['normal_' + y for y in x]
            for x in normal_mosaic
        ]

        mosaic = list()
        mosaic.extend(tumor_mosaic)
        mosaic.extend(normal_mosaic)

        # make fig, axd
        subplots_kwargs = (
            {
                'figsize': (30, 6 * len(mosaic)),
                'gridspec_kw': {'hspace': 0.6},
            } | subplots_kwargs
        )
        if figsize is not None:
            subplots_kwargs.update({'figsize': figsize})

        fig, axd = plt.subplot_mosaic(mosaic, **subplots_kwargs)

        return fig, axd

    @deco_nproc
    def draw_tumor_normal(
        self, 
        tumor_id, 
        normal_id,

        tumor_corrector_id=None,
        normal_corrector_id=None,
        genomeplotter=None,

        draw_tumor_all=False,
        draw_tumor_depth=True,
        draw_tumor_baf=True,
        draw_tumor_CN=True,
        draw_tumor_MQ=False,

        tumor_chromwise_peaks=True,

        draw_normal_all=False,
        draw_normal_depth=True,
        draw_normal_baf=True,
        draw_normal_CN=True,
        draw_normal_MQ=False,

        draw_common_kwargs=dict(),
        tumor_draw_kwargs=dict(),
        normal_draw_kwargs=dict(),

        subplots_kwargs=dict(),
        figsize=None,

        nproc=0,
        verbose=True,
    ):
        # prepare axd
        tumor_cnvsample = self.samples[tumor_id]
        normal_cnvsample = self.samples[normal_id]

        fig, axd = self.make_axd(
            tumor_cnvsample=tumor_cnvsample,
            normal_cnvsample=normal_cnvsample,

            draw_tumor_all=draw_tumor_all,
            draw_tumor_depth=draw_tumor_depth,
            draw_tumor_baf=draw_tumor_baf,
            draw_tumor_CN=draw_tumor_CN,
            draw_tumor_MQ=draw_tumor_MQ,

            draw_normal_all=draw_normal_all,
            draw_normal_depth=draw_normal_depth,
            draw_normal_baf=draw_normal_baf,
            draw_normal_CN=draw_normal_CN,
            draw_normal_MQ=draw_normal_MQ,

            subplots_kwargs=subplots_kwargs,
            figsize=figsize,
        )

        tumor_axd = dict()
        normal_axd = dict()
        for key, ax in axd.items():
            if key.startswith('tumor_'):
                tumor_axd[re.sub('^tumor_', '', key)] = ax
            if key.startswith('normal_'):
                normal_axd[re.sub('^normal_', '', key)] = ax

        # draw
        if genomeplotter is None:
            genomeplotter = tumor_cnvsample.genomeplotter

        fig, _, genomeplotter = tumor_cnvsample.draw(
            axd=tumor_axd,
            genomeplotter=genomeplotter,

            draw_all=draw_tumor_all,
            draw_depth=draw_tumor_depth,
            draw_baf=draw_tumor_baf,
            draw_CN=draw_tumor_CN,
            draw_MQ=draw_tumor_MQ,

            corrector_id=tumor_corrector_id,
            chromwise_peaks=tumor_chromwise_peaks,

            ylabel_prefix='tumor ',

            nproc=nproc,
            verbose=verbose,

            draw_common_kwargs=draw_common_kwargs,
            **tumor_draw_kwargs,
        )

        fig, _, genomeplotter = normal_cnvsample.draw(
            axd=normal_axd,
            genomeplotter=genomeplotter,

            draw_all=draw_normal_all,
            draw_depth=draw_normal_depth,
            draw_baf=draw_normal_baf,
            draw_CN=draw_normal_CN,
            draw_MQ=draw_normal_MQ,

            corrector_id=normal_corrector_id,
            chromwise_peaks=False,

            ylabel_prefix='normal ',

            nproc=nproc,
            verbose=verbose,

            title=None,

            draw_common_kwargs=draw_common_kwargs,
            **normal_draw_kwargs,
        )

        if 'CN' in normal_axd:
            normal_axd['CN'].set_ylabel('germline copy number')

        return fig, axd, genomeplotter


class CNVDrawingResult:
    default_fig_width = 12
    depth_picker_cursor_style = dict(
        lw=1, 
        color='red',
    )
    depth_picker_pickedval_style = dict(
        lw=1, 
        color='purple',
    )
    distrib_ax_label = 'distrib'
    solpick_common_title = f'x axis: number of divisions of selected depth range'

    def __init__(
        self,
        cnvsample,
        corrector_id,
        fig, 
        axd,
        genomeplotter,
        vplist=None,
    ):
        self.cnvsample = cnvsample
        self.corrector_id = corrector_id
        self.fig = fig
        self.axd = axd
        self.genomeplotter = genomeplotter
        self.vplist = vplist

        #self.predict_lines = None

        self.init_depthpick_attrs()
        self.init_solutionpick_attrs()

    def init_depthpick_attrs(self):
        self.depthpick_attrs = {
            'fig': None,
            'axd': None,
            'valid_axd': None,

            'cursor_lines': dict(),
            'depthax_cursor_line': None,

            'cid': None,
            'picked_values': list(),
            'picked_lines': list(),
            'depthax_picked_lines': list(),
        }

    def init_solutionpick_attrs(self):
        self.solutionpick_attrs = {
            'data_line2ds': None,
            'pickmark': None,

            'lastpick_axlabel': None,
            'lastpick_index': None,

            'bestmark': None,
            'axd': None,

            'cid': None,
            'fig': None,

            'best_solution_text': None,
            #'events': list(),
            'draw_on_pick': None,

            'predict_lines': None,
            'ccf_dots': None,
        }

    ################
    # depth picker #
    ################

    def draw_depth_picker( 
        self, 
        chroms=None, 
        ncol=4, 
        mq_limit=50,
        bins=None,
        empty_prefix='__EMPTY__', 
        all_axlabel='__ALL__', 
        figsize=None, 
        gridspec_kw=dict(),
    ):
        assert 'depth' in self.axd.keys()

        if bins is None:
            depth_ax_ylim = self.axd['depth'].get_ylim()
            bins = np.linspace(depth_ax_ylim[0], depth_ax_ylim[1], 50)

        # remove existing picked lines
        self.remove_previous_depthpick()

        # initiate depthpick_attrs
        self.init_depthpick_attrs()

        # chroms handling
        if chroms is None:
            all_chroms = tools.unique_keeporder(self.cnvsample.genomeplotter.region_gdf.chroms)
            chroms = [
                x for x in all_chroms
                if refgenome.normalize_chrom(x, strip_chr=True).isdecimal()
            ]
        else:
            chroms = np.atleast_1d(chroms)
            assert chroms.ndim == 1
            chroms = list(chroms)

        assert not any(x.startswith(empty_prefix) for x in chroms)
        assert all_axlabel not in chroms
        chroms.insert(0, all_axlabel)
            
        # make mosaic
        nrow = np.ceil(len(chroms) / ncol).astype(int)
        mosaic = np.empty((nrow, ncol), dtype=object)
        mosaic.flat[:len(chroms)] = chroms
        mosaic.flat[len(chroms):] = [f'{empty_prefix}{x}' for x in range(mosaic.size - len(chroms))] 

        # make axd
        if figsize is None:
            figsize = (12, 2 * nrow)
        gridspec_kw = (dict(hspace=0.5) | gridspec_kw)
        fig, axd = plt.subplot_mosaic(mosaic, figsize=figsize, gridspec_kw=gridspec_kw)
        fig.suptitle(
            '\n'.join(
                [
                    f'DEPTH DISTRIBUTIONS',
                    f'Axes title: chromosome name',
                    f'y axis: density',
                    f'x axis: normalized depth',
                ]
            )
        )
        valid_axd = {key: val for key, val in axd.items() if (not key.startswith(empty_prefix))}
        for ax in [ax for (key, ax) in axd.items() if key.startswith(empty_prefix)]:
            plotmisc.clear_ax(ax)

        # set xlim
        bins_min = bins.min()
        bins_max = bins.max()
        ptp = bins.ptp()
        xmin = bins_min - ptp * 0.05
        xmax = bins_max + ptp * 0.05
        for ax in valid_axd.values():
                ax.set_xlim(xmin, xmax)

        # draw hists
        depth_seg = self.cnvsample.get_depth_segment(corrector_id=self.corrector_id)
        depth_seg_subset = depth_seg.loc[
            depth_seg[depth_seg.__class__.MQ_mean_colname] > mq_limit
        ]
        depth_seg_bychrom = depth_seg_subset.group_bychrom(sort=False)
        for chrom, ax in valid_axd.items():
            if chrom in depth_seg_bychrom:
                _ = ax.hist(
                    depth_seg_bychrom[chrom].norm_depth_mean, 
                    bins=bins,
                    weights=depth_seg_bychrom[chrom].lengths,
                    density=True,
                )
                ax.set_title(chrom)
        
        allax_gdf = depth_seg_subset.choose_chroms(chroms)
        _ = valid_axd[all_axlabel].hist(
            allax_gdf.norm_depth_mean, 
            bins=bins,
            weights=allax_gdf.lengths,
            density=True,
        )
        valid_axd[all_axlabel].set_title('ALL')

        # keep fig and axd
        self.depthpick_attrs['fig'] = fig
        self.depthpick_attrs['axd'] = axd
        self.depthpick_attrs['valid_axd'] = valid_axd

        # make cursor lines
        self.depthpick_attrs['cursor_lines'] = {
            chrom: ax.axvline(**self.__class__.depth_picker_cursor_style)
            for chrom, ax in self.depthpick_attrs['valid_axd'].items()
        }
        self.depthpick_attrs['depthax_cursor_line'] = self.axd['depth'].axhline(
            **self.__class__.depth_picker_cursor_style,
        )

        # connect event handlers
        self.depthpick_attrs['cid'] = {
            'move': self.depthpick_attrs['fig'].canvas.mpl_connect('motion_notify_event', self.depth_picker_on_move),
            'click': self.depthpick_attrs['fig'].canvas.mpl_connect('button_press_event', self.depth_picker_on_click),
        }

    def depth_picker_on_move(self, event):
        for linecol in self.depthpick_attrs['cursor_lines'].values():
            linecol.set_xdata(event.xdata)
        self.depthpick_attrs['depthax_cursor_line'].set_ydata(event.xdata)

    def depth_picker_on_click(self, event):
        self.depthpick_attrs['picked_values'].append(event.xdata)

        for ax in self.depthpick_attrs['valid_axd'].values():
            self.depthpick_attrs['picked_lines'].append(
                ax.axvline(event.xdata, **self.__class__.depth_picker_pickedval_style)
            )
        self.depthpick_attrs['depthax_picked_lines'].append(
            self.axd['depth'].axhline(event.xdata, **self.__class__.depth_picker_pickedval_style)
        )

    def disconnect_depth_picker(self):
        self.depthpick_attrs['fig'].canvas.mpl_disconnect(
            self.depthpick_attrs['cid']['move']
        )
        self.depthpick_attrs['fig'].canvas.mpl_disconnect(
            self.depthpick_attrs['cid']['click']
        )

    def remove_previous_depthpick(self):
        # remove existing picked lines
        artists_to_remove = list()
        artists_to_remove.extend(
            itertools.chain(
                self.depthpick_attrs['picked_lines'],
                self.depthpick_attrs['depthax_picked_lines'],
                self.depthpick_attrs['cursor_lines'].values(),
            )
        )
        if self.depthpick_attrs['depthax_cursor_line'] is not None:
            artists_to_remove.append(
                self.depthpick_attrs['depthax_cursor_line']
            )

        for x in artists_to_remove:
            try:
                x.remove()
            except ValueError:
                pass

    ###########################
    # set solution candidates #
    ###########################

    def set_solution_candidates(self, **kwargs):
        kwargs = dict(
            num_ndiv_cand=20,
            num_CNt0_cand=10,
            #CNt_handicap_factor=0.003,
            ndivs_offset_step=None,
            depthfit_cutoff=0.3,
            baffit_cutoff=0.3,
        ) | kwargs

        CN_gdf = self.cnvsample.get_merged_segment(corrector_id=self.corrector_id)
        sol_candidates = cnvcall.solution_from_depth_limits(
            picked_depths=self.depthpick_attrs['picked_values'],
            depths=CN_gdf.norm_depth_mean,
            bafs=CN_gdf.get_corrected_baf(baf_index='baf0'),
            depth_stds=CN_gdf.norm_depth_std,
            baf_stds=CN_gdf.get_baf_std(baf_index='baf0'),

            CNg=CN_gdf.get_clonal_CN(germline=True),
            Bg=CN_gdf.get_clonal_B(baf_index='baf0', germline=True),
            lengths=CN_gdf.lengths,
            **kwargs,
        )

        self.cnvsample.init_solution_attrs_item(self.corrector_id)
        #self.solution_attrs[corrector_id]['selected_regions'] = self.draw_result.get_region_gdfs()
        self.cnvsample.solution_attrs[self.corrector_id]['solution_candidates'] = sol_candidates

    def set_ccf_candidates(
        self, 
        vcf_sampleid=None, 
        exclude_other=False,
        clonal_pval=cnvcall.DEFAULT_CLONAL_PVALUE,
    ):
        if vcf_sampleid is None:
            vcf_sampleid = self.cnvsample.sampleid

        assert self.vplist is not None
        vp_gdf = self.vplist.get_vafdf(vcf_sampleid, add_depth=True, exclude_other=exclude_other)

        # get indexes of CNV segments corresponding to input mutations
        cnv_gdf_cp = self.cnvsample.get_merged_segment(self.corrector_id).copy()  # must be non-merged cnv segment
        gCN_colname = cnv_gdf_cp.get_clonal_CN_colname(germline=True)
        cnv_gdf_cp['index'] = np.arange(cnv_gdf_cp.nrow)
        vp_gdf_joined = vp_gdf.join(
            cnv_gdf_cp, 
            ['index', gCN_colname],
            merge='longest', 
            how='left', 
            suffixes={'longest': ''},
        )
        seg_indexes = vp_gdf_joined['index']
        index_isnan = np.isnan(seg_indexes)
        selector = np.where(index_isnan, 0, seg_indexes).astype(int)

        # select solution candidate elements
        solution_cands = self.cnvsample.get_solution_candidates(self.corrector_id)
        CNt = solution_cands['CNt'].take(selector, axis=0)
        CNt[index_isnan, :] = np.nan
        Bt = solution_cands['Bt'].take(selector, axis=0)
        Bt[index_isnan, :] = np.nan
        CNg = vp_gdf_joined[gCN_colname]
        CNg[index_isnan] = np.nan

        # make variant-related parameters
        vafs = vp_gdf_joined.get_format(vcf_sampleid, 'ALT1_vaf')
        total_depths = vp_gdf_joined['total_depth']
        alt_depths = vp_gdf_joined.get_format(vcf_sampleid, 'ALT1_depth').astype(float)
        alt_depths[alt_depths == 0] = np.nan

        # prepare cellularity
        cellularity = self.cnvsample.get_solution_candidates(self.corrector_id)['cellularity']

        # main
        CNm, ccf = cnvcall.find_ccf(
            vaf=vafs,
            total_depth=total_depths,
            alt_depth=alt_depths,
            CNg=CNg,
            CNt=CNt,
            Bt=Bt,
            cellularity=cellularity,
            clonal_pval=clonal_pval,
        )

        solution_cands['CNm'] = CNm
        solution_cands['ccf'] = ccf
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)  # mean of empty slices
            solution_cands['ccf_mean'] = np.nanmean(ccf, axis=0)
            solution_cands['clonal_fraction'] = np.nanmean((ccf == 1), axis=0)

    #####################################
    # show cellularity / ploidy distrib #
    #####################################

    def show_ploidy(self):
        ploidies = self.cnvsample.get_solution_candidates(self.corrector_id)['ploidy']
        minval = np.floor(np.nanmin(ploidies)) 
        maxval = np.ceil(np.nanmax(ploidies))
        xticks = np.arange(0, ploidies.shape[1] + 5, 5)

        fig, axs = plt.subplots(ploidies.shape[0], 1, figsize=(12, 7), sharex=True)
        fig.suptitle('ploidy distribution')
        fig.supxlabel('number of division')
        for idx, ax in enumerate(axs):
            ax.set_ylim(minval, maxval)
            ax.plot(ploidies[idx, :], '-o', lw=1)
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticks + 1)
            ax.axhline(7, color='tab:red')

    def show_cellularity(self):
        cellularities = self.cnvsample.get_solution_candidates(self.corrector_id)['cellularity']
        xticks = np.arange(0, cellularities.shape[1] + 5, 5)

        fig, axs = plt.subplots(cellularities.shape[0], 1, figsize=(12, 7), sharex=True)
        fig.suptitle('cellularity distribution')
        fig.supxlabel('number of division')
        for idx, ax in enumerate(axs):
            ax.set_ylim(0, 1)
            ax.plot(cellularities[idx, :], '-o', lw=1)
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticks + 1)

    ###################
    # solution picker #
    ###################

    def draw_solution_picker(self, figsize=(12, 10), gridspec_kw=dict(), all_targetval_offset=5, draw_on_pick=False):
        """Must be run after set_solution_candidates"""

        assert self.corrector_id in self.cnvsample.solution_attrs
        assert self.cnvsample.solution_attrs[self.corrector_id]['solution_candidates'] is not None

        self.init_solutionpick_attrs()
        self.solutionpick_attrs['draw_on_pick'] = draw_on_pick

        sol_cands = self.cnvsample.get_solution_candidates(corrector_id=self.corrector_id)
        assert sol_cands is not None

        #all_targetval_exist = (sol_cands['all_targetval'] is not None)

        # set data
        targetvals = sol_cands['targetval']
        depth_targetvals = sol_cands['depth_targetval']
        baf_targetvals = sol_cands['baf_targetval']
        #clonal_fractions = sol_cands['clonal_fraction']
        CNt0s = sol_cands['CNt0']

        if 'ccf_mean' in sol_cands:
            ccf_means = sol_cands['ccf_mean']
        else:
            ccf_means = None

        if 'ploidy' in sol_cands:
            ploidies = sol_cands['ploidy']
        else:
            ploidies = None

        # set axes ylim
        ymin = 0
        ymax = 1
        ywidth = ymax - ymin
        ylim = (ymin - 0.2 * ywidth, ymax + 0.2 * ywidth)
        #ploidy_ylim = (np.nanmin(ploidies), np.nanmax(ploidies))
        ploidy_ylim = (0, 8)
        ploidy_yticks = np.arange(0, ploidy_ylim[1] + 2, 2)

        # normalize data
        targetvals = tools.fit_data_to_range(targetvals, ymin, ymax)
        depth_targetvals = tools.fit_data_to_range(depth_targetvals, ymin, ymax)
        baf_targetvals = tools.fit_data_to_range(baf_targetvals, ymin, ymax)

        #if all_targetval_exist:
        #    all_targetvals = sol_cands['all_targetval'] + all_targetval_offset
        #    all_depth_targetvals = sol_cands['all_depth_targetval'] + all_targetval_offset
        #    all_baf_targetvals = sol_cands['all_baf_targetval'] + all_targetval_offset

        # set axes settings
        xmin = -2
        #xmax = targetvals.shape[1]
        max_valid_idx1 = self.get_max_valid_idx1()
        xmax = max_valid_idx1 + 3

        xticks = np.arange(0, max_valid_idx1 + 5, 5)
        xticks = xticks[xticks <= max_valid_idx1]
        xticklabels = sol_cands['num_division'][0, xticks]
        minor_xticks = np.arange(0, max_valid_idx1 + 1, 1)

        #if all_targetval_exist:
        #    arrs_to_concat = [
        #        targetvals.ravel(), depth_targetvals.ravel(), baf_targetvals.ravel(),
        #        all_targetvals.ravel(), all_depth_targetvals.ravel(), all_baf_targetvals.ravel(),
        #    ]
        #else:
        #arrs_to_concat = [
        #    targetvals.ravel(), depth_targetvals.ravel(), baf_targetvals.ravel(),
        #]
        #targetvals_concat = np.concatenate(arrs_to_concat)
        #ymin = np.nanmin(targetvals_concat)
        #ymax = np.nanmax(targetvals_concat)

        # make mosaic
        mosaic = list()
        for idx, CNt0 in enumerate(CNt0s[:, 0]):
            label = str(idx)
            mosaic.append(label)
        mosaic = [[x] for x in mosaic]
        mosaic.append([self.__class__.distrib_ax_label])

        # draw
        color_sum = 'black'
        color_depth = 'tab:orange'
        color_baf = 'tab:blue'
        color_ccf = 'tab:green'
        color_ploidy = 'tab:red'
        markersize = 3
        pickmark_lw = 5

        self.solutionpick_attrs['data_line2ds'] = dict()
        self.solutionpick_attrs['pickmark'] = dict()
        gridspec_kw = (
            {
                'height_ratios': ([1] * (len(mosaic) - 1)) + [3],
                'hspace': 0.8,
            }
            | gridspec_kw
        )
        fig, axd = plt.subplot_mosaic(mosaic, figsize=figsize, gridspec_kw=gridspec_kw)
        fig.suptitle(self.__class__.solpick_common_title)
        for axidx, (label, ax) in enumerate(axd.items()):
            if label == self.__class__.distrib_ax_label:
                continue

            idx = int(label)

            #############
            # main axes #
            #############

            _ = ax.plot(
                depth_targetvals[idx, :], '-o', 
                color=color_depth, markersize=markersize, alpha=0.5, linewidth=1, 
            )
            _ = ax.plot(
                baf_targetvals[idx, :], '-o', 
                color=color_baf, markersize=markersize, alpha=0.5, linewidth=1, 
            )
            _ = ax.plot(
                targetvals[idx, :], '-o', 
                color=color_sum, markersize=markersize, alpha=0.5, linewidth=1,
            )

            if ccf_means is not None:
                _ = ax.plot(
                    ccf_means[idx, :], '-o', 
                    color=color_ccf, markersize=markersize, alpha=0.5, linewidth=1,
                )

            picker_line2d = ax.axvline(0, color='tab:red', alpha=0.5, linewidth=pickmark_lw, visible=False)
            self.solutionpick_attrs['pickmark'][label] = picker_line2d

            #bestmark_linecol = ax.vlines([], *ax.get_ylim(), alpha=0.3, color='tab:green', linewidth=4, visible=False)
            #self.solutionpick_attrs['bestmark'][label] = bestmark_linecol

            #if axidx == len(axd) - 2:
            #    ax.set_xlabel('number of division', size='small')

            ax.set_xlim(xmin, xmax)
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
            ax.set_xticks(minor_xticks, minor=True)

            ax.set_ylim(ylim)
            ax.set_yticks([ymin, ymax])
            ax.set_ylabel(f'CNt0={CNt0s[idx, 0]}', size=10, rotation=0, labelpad=35, va='center')

            ax.tick_params(which='both', length=0, labelsize='small')
            ax.grid(which='both', axis='x')
            #ax.grid(which='major', axis='y')

            #########################
            # twinx axes for ploidy #
            #########################

            ax2 = ax.twinx()
            ax2.set_label(label)
            ax2.set_ylim(ploidy_ylim)
            if ploidies is not None:
                ax2.plot(
                    ploidies[idx, :], '-o', 
                    color=color_ploidy, markersize=markersize, alpha=0.5, linewidth=1,
                )
            ax2.set_yticks(ploidy_yticks)
            ax2.grid(which='major', axis='y')

        # initialize distrib axes
        axd[self.__class__.distrib_ax_label].clear()

        # draw legend
        handles = plotmisc.LegendHandles()
        handles.add_line(marker='o', linewidth=0, color=color_depth, label='depth')
        handles.add_line(marker='o', linewidth=0, color=color_baf, label='baf')
        handles.add_line(marker='o', linewidth=0, color=color_sum, label='depth + baf')

        if ccf_means is not None:
            handles.add_line(marker='o', linewidth=0, color=color_ccf, label='mean ccf')
        if ploidies is not None:
            handles.add_line(marker='o', linewidth=0, color=color_ploidy, label='ploidy')

        fig.legend(handles=handles, loc='upper right', bbox_to_anchor=(0.9, 0.99))

        # connect
        cid = fig.canvas.mpl_connect('button_press_event', self.solution_picker_on_press)
        self.solutionpick_attrs['cid'] = cid
        self.solutionpick_attrs['fig'] = fig

        # keep axd
        self.solutionpick_attrs['axd'] = axd

    def solution_picker_on_press(self, event):
        lastpick_label = self.solutionpick_attrs['lastpick_axlabel']
        if lastpick_label is not None:
            self.solutionpick_attrs['pickmark'][lastpick_label].set_visible(False)

        #self.solutionpick_attrs['events'].append(event)

        if event.inaxes is None:
            return

        label = event.inaxes.get_label()
        idx0 = int(label)
        idx1 = int(np.rint(event.xdata))
        picked_index = (idx0, idx1)

        if picked_index not in self.get_valid_solution_indexes():
            return

        # make pickmark visible
        pickmark = self.solutionpick_attrs['pickmark'][label]
        #pickmark.set_data([line_xs[idx1]], [line_ys[idx1]])
        pickmark.set_data([idx1, idx1], pickmark.get_data()[1])
        pickmark.set_visible(True)

        # change suptitle
        sol_cands = self.cnvsample.get_solution_candidates(self.corrector_id)
        picked_CNt0 = sol_cands['CNt0'][idx0, idx1].astype(int)
        picked_dd = sol_cands['onecopy_depth'][idx0, idx1]
        picked_ndivs = sol_cands['num_division'][idx0, idx1]
        picked_cellularity = sol_cands['cellularity'][idx0, idx1]
        picked_K = sol_cands['K'][idx0, idx1]

        title = '\n'.join([
            self.__class__.solpick_common_title + '\n',
            f'index={picked_index}',
            f'CNt0={picked_CNt0}, onecopy_depth={picked_dd:.3}, num_division={picked_ndivs}',
            f'cellularity={picked_cellularity:.3}, K={picked_K:.3}',
        ])
        self.solutionpick_attrs['fig'].suptitle(title, size=10)

        # assign
        self.solutionpick_attrs['lastpick_axlabel'] = label
        self.solutionpick_attrs['lastpick_index'] = picked_index
        self.cnvsample.set_selected_solution(self.corrector_id, picked_index)

        # draw targetval distribution
        self.draw_targetval_distribution(
            self.solutionpick_attrs['axd'][self.__class__.distrib_ax_label]
        )

        # draw solution on main figure
        if self.solutionpick_attrs['draw_on_pick']:
            self.draw_solution()

    def draw_targetval_distribution(self, ax=None):
        targetvals = self.cnvsample.get_solution_candidates(self.corrector_id)['targetval']
        chosen_index = self.cnvsample.get_selected_solution(self.corrector_id)['index']

        if ax is None:
            fig, ax = plt.subplots()

        ax.clear()

        ax.hist(targetvals.ravel(), bins=30)
        ax.axvline(targetvals[chosen_index], color='black')
        ax.annotate(
            'chosen solution', 
            (targetvals[chosen_index], ax.get_ylim()[1] * 0.8), 
            (30, 0),
            textcoords='offset points',
            arrowprops=dict(width=0.1, headwidth=3, color='black'),
            ha='left', va='center',
        )
        ax.set_title('Distribution of target values', y=1)

    def mark_best_solutions(self, q):
        if self.solutionpick_attrs['best_solution_text'] is not None:
            self.solutionpick_attrs['best_solution_text'].remove()

        if self.solutionpick_attrs['bestmark'] is not None:
            for x in self.solutionpick_attrs['bestmark']:
                x.remove()

        indexes, CNt0s, ndivs = self.cnvsample.list_bestfit_solutions(q=q, corrector_id=self.corrector_id)
        markers = list()
        for idx in zip(*indexes):
            idx0, idx1 = idx
            label = str(idx0)
            ax = self.solutionpick_attrs['axd'][label]
            line = ax.axvline(idx1, alpha=0.2, color='tab:green', linewidth=10)
            markers.append(line)
        self.solutionpick_attrs['bestmark'] = markers

        self.solutionpick_attrs['best_solution_text'] = self.solutionpick_attrs['fig'].text(
            0.05, 0.95,
            f'green boxes: best {q * 100:.3} % solutions',
            ha='left', va='center',
        )

    def disconnect_solution_picker(self):
        self.solutionpick_attrs['fig'].canvas.mpl_disconnect(
            self.solutionpick_attrs['cid']
        )

    def get_max_valid_idx1(self):
        sol_cands = self.cnvsample.get_solution_candidates(corrector_id=self.corrector_id)
        assert sol_cands is not None

        targetvals = sol_cands['targetval']
        depth_targetvals = sol_cands['depth_targetval']
        baf_targetvals = sol_cands['baf_targetval']

        max_valid_idx1 = max(
            (~np.isnan(targetvals)).nonzero()[1].max(),
            (~np.isnan(depth_targetvals)).nonzero()[1].max(),
            (~np.isnan(baf_targetvals)).nonzero()[1].max(),
        )

        return max_valid_idx1

    def get_valid_solution_indexes(self):
        sol_cands = self.cnvsample.get_solution_candidates(corrector_id=self.corrector_id)
        assert sol_cands is not None
        targetvals = sol_cands['targetval']
        result = list(zip(*np.nonzero(~np.isnan(targetvals))))
        return result

    ########
    # draw #
    ########

    def make_ccf_added_vp_gdf(self, vplist):
        vp_gdf = vplist.get_vafdf(self.cnvsample.sampleid, add_depth=True, exclude_other=False)

        solution = self.cnvsample.get_selected_solution(corrector_id=self.corrector_id)
        CN_gdf = self.cnvsample.get_merged_segment(corrector_id=self.corrector_id)
        vp_gdf.add_ccf(
            cnv_gdf=CN_gdf, 
            cellularity=solution['cellularity'], 
            vcf_sampleid=self.cnvsample.sampleid, 
            clonal_pval=cnvcall.DEFAULT_CLONAL_PVALUE,
        )
        return vp_gdf

    def draw_vplist(self, vplist, nproc=1):
        vp_gdf = self.make_ccf_added_vp_gdf(vplist)

        gdraw_result = vp_gdf.draw_vaf(
            ax=axd['mut'],
            genomeplotter=self.genomeplotter,
            setup_axes=False,
            title=None,
            nproc=nproc,
        )
        gdraw_result = vp_gdf.draw_ccf(
            ax=self.axd['mut'],
            genomeplotter=self.genomeplotter,
            setup_axes=False,
            title=None,
            nproc=nproc,
            vaf_legend=True,
        )

    def draw_solution(self, index=None, verbose=True, nproc=1, omit_title_mode=False, **kwargs):
        """self.cnvsample must be posessing solution candidates with the given corrector_id"""
        assert self.corrector_id in self.cnvsample.solution_attrs
        assert self.cnvsample.solution_attrs[self.corrector_id]['solution_candidates'] is not None

        nproc = 1

        if index is None:
            index = self.solutionpick_attrs['lastpick_index']
            if index is None:
                raise Exception(f'Solution index is not available')

        # remove previous solution drawing
        if self.solutionpick_attrs['predict_lines'] is not None:
            for x in self.solutionpick_attrs['predict_lines']:
                try:
                    x.remove()
                except ValueError:
                    pass
        if self.solutionpick_attrs['ccf_dots'] is not None:
            for x in self.solutionpick_attrs['ccf_dots']:
                try:
                    x.remove()
                except ValueError:
                    pass

        # init predict lines cache
        self.solutionpick_attrs['predict_lines'] = list()
        self.solutionpick_attrs['ccf_dots'] = list()

        # clear CN axes
        plotmisc.clear_ax(self.axd['CN'])

        # add solution to CN gdf
        solution = self.cnvsample.pick_solution(index, corrector_id=self.corrector_id)
        CN_gdf = self.cnvsample.get_merged_segment(corrector_id=self.corrector_id)
        CN_gdf.add_clonal_solution_targetsample(
            cellularity=solution['cellularity'], 
            K=solution['K'],
        )

        # draw CN
        if 'CN_legend' in self.axd:
            legend_ax = self.axd['CN_legend']
        else:
            legend_ax = None
        CN_gdf.draw_CN(
            ax=self.axd['CN'],
            legend_ax=legend_ax,
            genomeplotter=self.genomeplotter,
            setup_axes=True,
            title=None,
            nproc=nproc,
            verbose=verbose,
            **kwargs,
        )

        # draw predicted values
        gdraw_result = CN_gdf.draw_predicted_depth(
            ax=self.axd['depth'],
            genomeplotter=self.genomeplotter,
            nproc=nproc,
        )
        self.solutionpick_attrs['predict_lines'].append(gdraw_result.artist)

        baf_axd = {
            k: v for (k, v) in self.axd.items()
            if libbaf.check_valid_bafindex(k)
        }
        for baf_index, ax in baf_axd.items():
            gdraw_result = CN_gdf.draw_predicted_baf(
                baf_index=baf_index,
                ax=ax,
                genomeplotter=self.genomeplotter,
                nproc=nproc,
            )
            self.solutionpick_attrs['predict_lines'].append(gdraw_result.artist)

        # draw ccf
        if self.vplist is not None:
            vp_gdf = self.make_ccf_added_vp_gdf(self.vplist)
            gdraw_result = vp_gdf.draw_ccf(
                ax=self.axd['mut'],
                genomeplotter=self.genomeplotter,
                setup_axes=False,
                title=None,
                nproc=nproc,
                vaf_legend=True,
            )
            self.solutionpick_attrs['ccf_dots'].append(gdraw_result.artist)

        # set title
        title = self.cnvsample.make_title(
            selected_solution=self.cnvsample.pick_solution(
                index=index, corrector_id=self.corrector_id,
            ),
            omit_mode=omit_title_mode,
        )
        plotmisc.draw_suptitle(self.fig, title)

    # redraw main figure

    def redraw_mainfig(self, save=False, **kwargs):
        draw_result = self.cnvsample.draw(corrector_id=self.corrector_id, **kwargs)
        if save:
            self.fig = draw_result.fig
            self.axd = {ax.get_label(): ax for ax in self.fig.get_axes()}

    ###########
    # savefig #
    ###########

    def savefig(self, topdir, force=False):
        if not force:
            assert not os.path.exists(topdir)
        os.makedirs(topdir, exist_ok=True)

        main_path = os.path.join(topdir, f'cnv.pdf')
        self.fig.savefig(main_path)

        depth_picker_path =  os.path.join(topdir, f'depth_picker.pdf')
        self.depthpick_attrs['fig'].savefig(depth_picker_path)

        solution_picker_path =  os.path.join(topdir, f'solution_picker.pdf')
        self.solutionpick_attrs['fig'].savefig(solution_picker_path)


class CNVDrawingResult_old(GenomeDrawingFigureResult):
    def set_ylim(self, key, *, ymin=None, ymax=None):
        ax = self.axd[key]
        old_ymin, old_ymax = ax.get_ylim()
        if ymin is None:
            ymin = old_ymin
        if ymax is None:
            ymax = old_ymax
        ax.set_ylim(ymin, ymax)


class SolutionDrawingResult(GenomeDrawingFigureResult):
    def set_CN_ylim(self, *, ymin=None, ymax=None, yticks=None):
        for x in self.axresults['CN'].draw_common_artists:
            try:
                x.remove()
            except:
                pass
        if hasattr(self, 'old_artists'):
            for x in self.old_artists:
                try:
                    x.remove()
                except:
                    pass

        ax = self.axd['CN']
        old_ymin, old_ymax = ax.get_ylim()
        if ymin is None:
            ymin = old_ymin
        if ymax is None:
            ymax = old_ymax
        artists = genomedf_draw.draw_axessetup(
            ax=self.axd['CN'],
            genomeplotter=self.axresults['CN'].genomeplotter,
            ymin=ymin,
            ymax=ymax,
            yticks=yticks,
        )
        self.old_artists = artists


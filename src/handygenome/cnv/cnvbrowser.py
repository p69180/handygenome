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
from handygenome.variant.variantplus import VariantDataFrame
from handygenome.cnv.depth import DepthRawDataFrame as DepthRawDF
from handygenome.cnv.depth import DepthSegmentDataFrame as DepthSegDF
from handygenome.cnv.baf import BAFRawDataFrame as BAFRawDF
from handygenome.cnv.baf import BAFSegmentDataFrame as BAFSegDF
from handygenome.cnv.cnvsegment import CNVSegmentDataFrame as CNVSegDF

import handygenome.cnv.cnvcall as cnvcall
from handygenome.cnv.cnvcall import KCError
import handygenome.cnv.mosdepth as libmosdepth
import handygenome.cnv.baf as libbaf
import handygenome.cnv.bafcorrection as bafcorrection

import handygenome.plot.genomeplot as libgenomeplot
import handygenome.plot.misc as plotmisc
from handygenome.plot.genomeplot import GenomePlotter


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
        nproc=os.cpu_count(),
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
        result.simple_attrs = {
            'bam_path': bam_path,
            'vcf_path': vcf_path,
            'vcf_sampleid': vcf_sampleid,
            'sampleid': sampleid,
            'refver': refgenome.standardize(refver),
            'is_female': is_female,
            'is_germline': is_germline,
            'mode': mode,
            'nproc': nproc,
            'verbose': verbose,
            'ploidy': ploidy,
            #'genomeplotter_kwargs': genomeplotter_kwargs,
        }

        # initiate gdf data dict
        result.init_gdfs()

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
        }
        self.corrected_gdfs = dict()

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
        merged_segment.add_clonal_solution_corrected_depth(cellularity=cellularity, K=K)

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
                merged_segment.add_clonal_solution_corrected_depth(cellularity=cellularity, K=K)
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
                subdic[subkey] = subdic[subkey].df

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

        # simple attrs
        self.save_simple_attrs(topdir)

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

        self.log(f'Finished saving')

    @classmethod
    def load_pickle_asdf(cls, topdir):
        if not os.path.exists(topdir):
            raise Exception(f'Directory does not exist: {topdir}')

        result = cls()

        # simple attrs
        result.load_simple_attrs(topdir)

        # initiate gdf data dict
        result.init_gdfs()

        # target region
        result.load_nonbafidx_gdf_aspickle('raw_target_region', topdir)
        result.load_nonbafidx_gdf_aspickle('target_region', topdir)
        result.load_nonbafidx_gdf_aspickle('excluded_region', topdir)

        # raw data
        result.load_nonbafidx_gdf_aspickle('depth_rawdata', topdir)
        result.load_nonbafidx_gdf_aspickle('baf_rawdata', topdir)
        #result.load_nonbafidx_gdf_aspickle('baf_hetalt_rawdata', topdir)

        # segment
        result.load_nonbafidx_gdf_aspickle('depth_segment', topdir)
        result.load_bafidx_gdf_aspickle('baf_segment', topdir)
        result.load_bafidx_gdf_aspickle('baf_noedit_segment', topdir)
        result.load_bafidx_gdf_aspickle('baf_edited_segment', topdir)
        result.load_bafidx_gdf_aspickle('baf_rmzero_noedit_segment', topdir)
        result.load_nonbafidx_gdf_aspickle('merged_segment', topdir)

        # corrected_gdfs
        cgdfs_dir = os.path.join(topdir, 'corrected_gdfs')
        for fname in os.listdir(cgdfs_dir):
            savepath = os.path.join(cgdfs_dir, fname)
            result.load_corrected_gdfs(savepath)

        # genomeplotter
        result.init_genomeplotter(**result.genomeplotter_kwargs)

        # plotdata
        result.plotdata_cache = dict()

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

        # simple attrs
        result.load_simple_attrs(topdir)

        # initiate gdf data dict
        result.init_gdfs()

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
    def make_depth_segment(self, nproc=0, **kwargs):
        depth_rawdata = self.get_depth_rawdata()
        depth_segment = depth_rawdata.get_segment(
            target_region=self.get_raw_target_region(),
            rawdepth=False,
            nproc=nproc,
        )
        depth_segment.add_rawdata_info(
            depth_rawdata, 
            merge_methods=['mean', 'std'],
            rawdepth=False,
            nproc=nproc,
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
    def make_baf_segment(self, nproc=0):
        assert self.is_germline
        self.gdfs['baf_segment'] = self.get_baf_rawdata().get_segment_dict(
            target_region=self.get_raw_target_region(), 

            germline_baf_rawdata=None,

            germline_CN_gdf=None, 
            is_female=None,

            nproc=nproc,
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

    ##########################
    # corrected gdfs - maker #
    ##########################

    def make_corrected_gdfs_init(self, corrector_id):
        self.corrected_gdfs[corrector_id] = {
            'id': corrector_id,
            'cellularity': None,
            'K': None,
            'ploidy': None,

            'depth_rawdata': None,
            'depth_segment': None,
            'depth_excluded_region': None,

            'baf_rawdata': None,
            'filtered_baf_rawdata': None,
            'baf_segment': None,
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
    ):
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
        #corrected_depth_rawdata[DepthRawDF.corrected_norm_depth_colname] = corrected_norm_depth
        corrected_depth_rawdata[DepthRawDF.norm_depth_colname] = corrected_norm_depth

        # 3. make segment from corrected norm depth
        self.log(f'Corrected data "{corrector_id}" - Making depth segment')
        corrected_depth_segment = corrected_depth_rawdata.get_segment(
            target_region=self.get_raw_target_region(),
            rawdepth=False,
            nproc=nproc,
        )

        # 4. add rawdata information to segment
        self.log(f'Corrected data "{corrector_id}" - Adding rawdata to depth segment')
        corrected_depth_segment.add_rawdata_info(
            corrected_depth_rawdata, 
            merge_methods=['mean', 'std'],
            nproc=nproc,
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
            self.make_corrected_gdfs_baf_pairedgermline(
                corrector_id=corrector_id, 
                germline_baf_rawdata=germline_baf_rawdata,
            )

        # segment
        self.log(f'Making BAF segment from the corrected data')

        corrected_baf_rawdata = self.corrected_gdfs[corrector_id]['baf_rawdata']
        skip_hetalt_filter = (not self.is_germline)
        corrected_baf_segment_dict, filtered_baf_rawdata = corrected_baf_rawdata.get_segment_dict(
            target_region=self.get_raw_target_region(),

            germline_baf_rawdata=germline_baf_rawdata,

            germline_CN_gdf=germline_CN_gdf, 
            is_female=self.is_female,
            skip_hetalt_filter=skip_hetalt_filter,

            return_filtered_rawdata=True,

            nproc=nproc, 
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
        nproc=0,
    ):
        self.log(f'Corrected data "{corrector_id}" - Making merged segment')

        corrected_gdfs = self.corrected_gdfs[corrector_id]

        sub_seg_gdfs = (
            [corrected_gdfs['depth_segment']] 
            + list(corrected_gdfs['baf_segment'].values())
        )
        merged_seg_gdf = CNVSegDF.from_segments(sub_seg_gdfs, nproc=nproc)

        cols_to_drop = (
            [DepthSegDF.filter_colname]
            + [x.get_filter_colname() for x in corrected_gdfs['baf_segment'].values()]
        )
        merged_seg_gdf.drop_annots(cols_to_drop, inplace=True)

        # subtract excluded region
        merged_seg_gdf = merged_seg_gdf.subtract(corrected_gdfs['depth_excluded_region'])

        # intersect with target region
        raw_target_region = self.get_raw_target_region()
        if raw_target_region is not None:
            merged_seg_gdf = merged_seg_gdf.intersect(raw_target_region)

        # add germline solution
        if self.is_germline:
            merged_seg_gdf.add_clonal_solution_germline(self.is_female, self.ploidy)

        # result
        corrected_gdfs['merged_segment'] = merged_seg_gdf

    @deco_nproc
    def make_corrected_gdfs_add_germline_CN(
        self, 
        corrector_id, 
        paired_germline_CN_gdf,
        nproc=0,
    ):
        self.log(f'Corrected data "{corrector_id}" - Adding germline CN to merged segment')

        corrected_merged_seg = self.corrected_gdfs[corrector_id]['merged_segment']

        gCN_colname = paired_germline_CN_gdf.get_clonal_CN_colname(germline=False)
        gB_colname_dict = paired_germline_CN_gdf.get_clonal_B_colname_dict()
        added_cols = [gCN_colname] + list(gB_colname_dict.values())
        corrected_merged_seg = corrected_merged_seg.join(
            paired_germline_CN_gdf,
            how='left',
            right_gdf_cols=added_cols,
            merge='longest',
            overlapping_length=True,
            omit_N=True,
            suffixes={'longest': ''},
            nproc=nproc,
        )
        corrected_merged_seg.assign_clonal_CN(
            #corrected_merged_seg.get_clonal_CN(germline=False),
            data=corrected_merged_seg[gCN_colname],
            germline=True,
        )
        for baf_index, colname in gB_colname_dict.items():
            corrected_merged_seg.assign_clonal_B(
                #corrected_merged_seg.get_clonal_B(baf_index, germline=False),
                data=corrected_merged_seg[colname],
                baf_index=baf_index,
                germline=True,
            )

        germline_CN_columns = (
            [corrected_merged_seg.get_clonal_CN_colname(germline=True)]
            + [
                corrected_merged_seg.get_clonal_B_colname(baf_index, germline=True)
                for baf_index in corrected_merged_seg.get_baf_indexes()
            ]
        )
        corrected_merged_seg.drop_annots(
            set(added_cols).difference(germline_CN_columns), 
            inplace=True,
        )

        self.corrected_gdfs[corrector_id]['merged_segment'] = corrected_merged_seg

    @deco_nproc
    def make_corrected_gdfs_add_germline_CN_without_pairednormal(
        self, 
        corrector_id, 
        nproc=0,
    ):
        self.log(f'Corrected data "{corrector_id}" - Adding germline CN to merged segment')

        corrected_merged_seg = self.corrected_gdfs[corrector_id]['merged_segment']

        cnvcall.add_default_CNg_Bg(
            gdf=corrected_merged_seg,
            is_female=self.is_female,
            inplace=True,
            nproc=nproc,
        )
        corrected_merged_seg.assign_clonal_CN(
            corrected_merged_seg[cnvcall.DEFAULT_CNG_COLNAME],
            germline=True,
        )
        corrected_merged_seg.assign_clonal_B(
            corrected_merged_seg[cnvcall.DEFAULT_BG_COLNAME],
            baf_index='baf0',
            germline=True,
        )
        corrected_merged_seg.drop_annots(
            [cnvcall.DEFAULT_CNG_COLNAME, cnvcall.DEFAULT_BG_COLNAME],
            inplace=True,
        )

        self.corrected_gdfs[corrector_id]['merged_segment'] = corrected_merged_seg

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
        draw_depth_peaks=False,
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
        draw_depth_peaks=False,
        subplots_kwargs=dict(),
        figsize=None,
    ):
        mosaic = self.make_mosaic(
            draw_depth=draw_depth,
            draw_baf=draw_baf,
            draw_MQ=draw_MQ,
            draw_CN=draw_CN,
            draw_depth_peaks=draw_depth_peaks,
        )

        default_subplots_kwargs = {
            'figsize': (30, 6 * len(mosaic)),
            'gridspec_kw': {'hspace': 0.6},
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
        draw_CN=False,

        corrector_id=None,
        depthtype=None,
        chromwise_peaks=False,
        chromwise_peaks_density_kwargs=dict(),
        chromwise_peaks_line_kwargs=dict(),
        draw_excluded=False,

        #draw_depth_peaks=False,
        #draw_chromwise_peak=True,

        frac=0.1,
        draw_depth_rawdata=True,
        draw_baf_rawdata=True,
        draw_MQ_rawdata=True,

        # figure title - only works when figure object is not already created
        figsize=None,
        subplots_kwargs=dict(),

        # kwargs for low level drawers
        draw_common_kwargs=dict(),
        depth_kwargs=dict(),
        #depth_peaks_kwargs=dict(),
        baf_kwargs=dict(),
        MQ_kwargs=dict(),
        CN_kwargs=dict(),
        excluded_kwargs=dict(),
        #chromwise_peak_bw_method=1,
        #chromwise_peak_limit=None,

        # axes setting
        ylabel_prefix='',
#        setup_axes=True,
#        ylabel=None,
#        ylabel_kwargs=dict(),
#        ymax=None,
#        ymin=None,
#        yticks=None,
#        draw_common_kwargs=dict(),
#        rotate_chromlabel=None,
        title=False,
        suptitle_kwargs=dict(),

        # multiprocessing
        nproc=0,
        verbose=True,
    ):
        if draw_all:
            draw_depth = True
            draw_baf = True
            draw_MQ = True
            draw_CN = True

        # sanitycheck
        if not any([draw_depth, draw_baf, draw_MQ, draw_CN]):
            raise Exception(f'At least one of drawing flags must be set.')

        #assert not ((not draw_depth) and draw_depth_peaks)

        if genomeplotter is None:
            genomeplotter = self.genomeplotter

        # make axd
        if axd is None:
            fig, axd = self.make_axd(
                draw_depth=draw_depth,
                draw_baf=draw_baf,
                draw_MQ=draw_MQ,
                draw_CN=draw_CN,
                draw_depth_peaks=False,
                subplots_kwargs=subplots_kwargs,
                figsize=figsize,
            )
        else:
            fig = next(iter(axd.values())).figure

        # prepare CN plotdata
        CN_gdf = self.get_merged_segment(corrector_id=corrector_id)
        #CN_exists = (CN_gdf is not None) and CN_gdf.check_has_CN()
        CN_exists = (CN_gdf is not None)
        if CN_exists:
            CN_plotdata = genomeplotter.make_plotdata(
                CN_gdf, 
                log_suffix=' (Copy number)', 
                nproc=nproc,
                verbose=verbose,
            )

        # draw CN
        if draw_CN:
            if CN_exists and CN_gdf.check_has_CN():
                CN_gdf.draw_CN(
                    ax=axd['CN'],
                    plotdata=CN_plotdata,
                    genomeplotter=genomeplotter,
                    setup_axes=True,
                    title=None,
                    ylabel_prefix=ylabel_prefix,
                    nproc=nproc,
                    verbose=verbose,
                    draw_common_kwargs=draw_common_kwargs,
                    **CN_kwargs,
                )

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
                baf_segment_gdf.draw_baf(
                    ax=baf_ax, 
                    genomeplotter=genomeplotter,
                    setup_axes=True,
                    title=None,
                    ylabel_prefix=ylabel_prefix,
                    nproc=nproc,
                    verbose=verbose,
                    draw_common_kwargs=draw_common_kwargs,
                    **baf_kwargs,
                )

            # segment - corrected baf
            if CN_exists:
                for baf_idx, baf_ax in baf_axd.items():
                    CN_gdf.draw_corrected_baf(
                        baf_index=baf_idx,
                        ax=baf_ax,
                        genomeplotter=genomeplotter,
                        plotdata=CN_plotdata,

                        nproc=nproc,
                        verbose=verbose,
                    )

            # baf rawdata
            if draw_baf_rawdata:
                #baf_rawdata_gdf = self.gdfs['baf_rawdata']
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
                        baf_rawdata_gdf.draw_baf(
                            baf_index, 
                            ax=baf_ax,
                            genomeplotter=genomeplotter,
                            plotdata=baf_rawdata_plotdata,
                            ylabel_prefix=ylabel_prefix,
                            setup_axes=False,
                            verbose=verbose,
                            **baf_kwargs,
                        )

            # predicted value
            if CN_exists:
                for baf_idx, baf_ax in baf_axd.items():
                    y_colname = CN_gdf.get_predicted_baf_colname(baf_index=baf_idx)
                    if y_colname in CN_gdf.columns:
                        CN_gdf.draw_hlines(
                            ax=baf_ax,
                            y_colname=y_colname,
                            plotdata=CN_plotdata,
                            setup_axes=False,
                            plot_kwargs=dict(
                                color='red',
                                linewidth=2,
                                alpha=1,
                            ),
                            verbose=verbose,
                        )

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
            if draw_depth_rawdata:
                depth_rawdata.draw_depth(
                    ax=axd['depth'],
                    genomeplotter=genomeplotter,
                    depthtype=depthtype,
                    plotdata=depth_raw_plotdata,
                    setup_axes=(not depth_segment_exists),
                    ylabel_prefix=ylabel_prefix,
                    verbose=verbose,
                    **depth_kwargs,
                )
            if depth_segment_exists:
                depth_segment.draw_depth(
                    ax=axd['depth'],
                    genomeplotter=genomeplotter,
                    depthtype=depthtype,
                    plotdata=depth_seg_plotdata,
                    chromwise_peaks=chromwise_peaks,
                    chromwise_peaks_density_kwargs=chromwise_peaks_density_kwargs,
                    chromwise_peaks_line_kwargs=chromwise_peaks_line_kwargs,
                    setup_axes=True,
                    draw_common_kwargs=draw_common_kwargs,
                    ylabel_prefix=ylabel_prefix,
                    verbose=verbose,
                    **depth_kwargs,
                )

            # predicted value
            if CN_exists:
                y_colname = CN_gdf.get_predicted_depth_colname()
                if y_colname in CN_gdf.columns:
                    CN_gdf.draw_hlines(
                        ax=axd['depth'],
                        y_colname=y_colname,
                        plotdata=CN_plotdata,
                        setup_axes=False,
                        verbose=verbose,
                        plot_kwargs=dict(
                            color='red',
                            linewidth=2,
                            alpha=1,
                        ),
                    )

        # draw MQ
        if draw_MQ:
            if draw_MQ_rawdata:
                depth_rawdata.draw_MQ(
                    ax=axd['MQ'],
                    genomeplotter=genomeplotter,
                    plotdata=depth_raw_plotdata,
                    setup_axes=(not depth_segment_exists),
                    ylabel_prefix=ylabel_prefix,
                    verbose=verbose,
                    **MQ_kwargs,
                )
            if depth_segment_exists:
                depth_segment.draw_MQ(
                    ax=axd['MQ'],
                    genomeplotter=genomeplotter,
                    plotdata=depth_seg_plotdata,
                    setup_axes=True,
                    draw_common_kwargs=draw_common_kwargs,
                    ylabel_prefix=ylabel_prefix,
                    verbose=verbose,
                    **MQ_kwargs,
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
                    setup_axes=False,
                )

        # title
        if title is False:
            if corrector_id is None:
                cellularity = None
                ploidy = None
            else:
                cellularity = self.corrected_gdfs[corrector_id]['cellularity']
                ploidy = self.corrected_gdfs[corrector_id]['ploidy']
            title = self.get_default_title(cellularity=cellularity, ploidy=ploidy)

        if title is not None:
            plotmisc.draw_suptitle(fig, title, **suptitle_kwargs)

        return fig, axd, genomeplotter

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

    def get_default_title(self, cellularity=None, ploidy=None):
        result = ', '.join(
            f'{key}={getattr(self, key)}' 
            for key in ['sampleid', 'is_female', 'mode']
        )
        if (cellularity is not None) or (ploidy is not None):
            result = (
                result 
                + '\n' 
                + f'cellularity={round(cellularity, 3)}, ploidy={round(ploidy, 3)}'
            )
        return result

    def get_baf_indexes(self):
        return libbaf.get_baf_indexes(self.simple_attrs["ploidy"])


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

    def load_cnvsample_pickle(self, cnvsample_topdir):
        logutils.log(f'Loading CNVSample from directory {cnvsample_topdir}')
        cnvsample = CNVSample.load_pickle_asdf(cnvsample_topdir)
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

    def save_pickle(self, topdir):
        os.makedirs(topdir, exist_ok=True)

        # simple attrs
        with open(self.get_simple_attrs_path(topdir), 'wt') as outfile:
            json.dump(self.simple_attrs, outfile)

        # save cnvsamples
        os.makedirs(self.get_cnvsample_dir(topdir), exist_ok=True)
        for sid in self.cnvsamples.keys():
            self.save_cnvsample(sid, topdir, force=True)

        # save depth correctors
        self.save_all_correctors(topdir)

    @classmethod
    def load_pickle(cls, topdir, cnvsamples=None, load_corrector=False):
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
            result.load_cnvsample_pickle(os.path.join(cnvsample_dir, fname))

        logutils.log(f'Finished loading CNVSamples')

        # correctors
        if load_corrector:
            logutils.log(f'Beginning loading correctors')
            corrector_dir = result.get_corrector_dir(topdir)
            for fname in sorted(os.listdir(corrector_dir)):
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
        force=False,
        nproc=0,
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
            )
        else:
            baf_segment_dict, filtered_baf_rawdata = cnvsample.get_baf_rawdata().get_segment_dict(
                target_region=cnvsample.get_raw_target_region(),

                germline_baf_rawdata=germline_baf_rawdata,

                germline_CN_gdf=germline_CN_gdf, 
                is_female=cnvsample.is_female,

                return_filtered_rawdata=True,

                nproc=nproc, 
            )
            cnvsample.corrected_gdfs[corrector_id]['baf_segment'] = baf_segment_dict
            cnvsample.corrected_gdfs[corrector_id]['filtered_baf_rawdata'] = filtered_baf_rawdata

        # make merged segment
        cnvsample.make_corrected_gdfs_merged_segment(
            corrector_id=corrector_id,
            nproc=nproc,
        )

    @deco_nproc
    def add_gCN_to_corrected_merged_segment(
        self, *, sid, corrector_id, paired_germline_id=None, nproc=0,
    ):
        cnvsample = self.cnvsamples[sid]
        if (cnvsample.mode == 'wgs') and (paired_germline_id is None):
            raise Exception(f'With WGS sample, paired germline sample ID must be given.')

        #corrected_gdfs = cnvsample.corrected_gdfs[corrector_id]
        #corrected_merged_seg = cnvsample.corrected_gdfs[corrector_id]['merged_segment']

        if cnvsample.mode == 'wgs':
            paired_germline_CN_gdf = self.cnvsamples[paired_germline_id].get_merged_segment()
            cnvsample.make_corrected_gdfs_add_germline_CN(
                corrector_id=corrector_id, 
                paired_germline_CN_gdf=paired_germline_CN_gdf, 
                nproc=nproc,
            )
        elif cnvsample.mode == 'panel':
            cnvsample.make_corrected_gdfs_add_germline_CN_without_pairednormal(
                corrector_id=corrector_id, 
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


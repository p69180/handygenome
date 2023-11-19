import re
import os
import pickle
import itertools
import collections
import functools
import operator

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize
import scipy.interpolate
import sklearn.mixture

import handygenome
import handygenome.tools as tools
import handygenome.deco as deco
import handygenome.logutils as logutils
#import handygenome.tools as tools
import handygenome.peakutils as peakutils
from handygenome.peakutils import HistPeaks, DensityPeaks, DensityGenerationFailure, NoPeakError
import handygenome.refgenome.refgenome as refgenome
import handygenome.variant.variantplus as libvp
from handygenome.variant.variantplus import VCFDataFrame, VariantDataFrame
from handygenome.genomedf.genomedf import GenomeDataFrame, SegmentDataFrame
import handygenome.genomedf.genomedf_draw as genomedf_draw
import handygenome.plot.misc as plotmisc
import handygenome.cnv.cnvcall as cnvcall
#import handygenome.cnv.bafsimul as bafsimul


BAFINDEX_PAT_STRING = '(?P<bafindex_prefix>baf)(?P<num>[0-9]+)'
BAFINDEX_PAT = re.compile(BAFINDEX_PAT_STRING)


##############
# exceptions #
##############

class FailedNoTroughMethod(Exception):
    pass


##########
# basics #
##########

def get_baf_from_vaf(vafs):
    """Args:
        vafs: 
            - 1st dimension means genomic positions, 2nd dimension means allelic vafs
            - Desired format contains vaf values for all alleles. 
            (ndim == 2, shape[1] >= 2)
            - If values for only one allele is present, 
            (ndim == 1 or (ndim == 2 and shape[1] == 1)) 
            a transformed copy of input is generated which has values for 
            all alleles.
    """
    # preprocessing "vafs" into a 2d array
    vafs = np.atleast_1d(vafs)
    if vafs.ndim == 1:
        vafs = np.stack([vafs, 1 - vafs], axis=1)
    elif vafs.ndim == 2:
        if vafs.shape[1] == 1:
            vafs = np.stack([vafs[:, 0], 1 - vafs[:, 0]], axis=1)
    else:
        raise Exception(f'Input array must not have ndim larger than 2.')
        
    # main
    vafs = np.sort(vafs, axis=1)
    bafs = vafs[:, :-1]
    return bafs


def get_baf_indexes(ploidy):
    assert ploidy >= 2
    result = [f'baf{x}' for x in range(ploidy - 1)]
    assert all((BAFINDEX_PAT.fullmatch(x) is not None) for x in result)
    return result


def get_bafindex_number(baf_index):
    return int(BAFINDEX_PAT.fullmatch(baf_index).group('num'))


def check_valid_bafindex(baf_index):
    return (BAFINDEX_PAT.fullmatch(baf_index) is not None)


def decal_baf(bafs):
    return np.concatenate([bafs, 1 - bafs])


###################
# BAFRawDataFrame #
###################

class BAFRawDataFrame(VariantDataFrame):
    is_het_colname = 'is_het'
    drawresult_name = 'baf_rawdata'

    default_plot_kwargs = dict(
        zorder=0,
    )

    def get_baf_indexes(self):
        return get_baf_indexes(self.num_alleles)

    ###############
    # sanitycheck #
    ###############

    def remove_overlapping_rows_slow(self):
        new_df = (
            self.cluster().df
            .groupby('Cluster', sort=False)
            .first()
            .reset_index(drop=True)
        )
        return self.spawn(new_df)

    def remove_overlapping_rows(self):
        """Intended to remove overlapping variants (e.g. pos 1 A>G & pos 1 AC>GT)
        """
        filtered_subgdfs = list()
        for chrom, subgdf in self.group_bychrom().items():
            while True:
                overlap_selector = np.concatenate(
                    [[False], (subgdf.starts[1:] <= subgdf.ends[:-1])], 
                    axis=0,
                )
                if not overlap_selector.any():
                    break
                else:
                    subgdf = subgdf.loc[~overlap_selector, :]

            filtered_subgdfs.append(subgdf)

        new_gdf = GenomeDataFrame.concat(filtered_subgdfs)
        self.assign_df(new_gdf.df)
        self.sort()

    ########################
    # normalize VAF values #
    ########################

    def fit_vaf_sum_to_one(self, inplace=True):
        vaf_array = self[self.vaf_columns]
        vaf_sums = vaf_array.sum(axis=1)
        edited_vafs = vaf_array / vaf_sums[:, np.newaxis]
        if inplace:
            self[self.vaf_columns] = edited_vafs
        return edited_vafs

    def divide_with_REF_vaf(self):
        vaf_array = self[self.vaf_columns]
        oddsratio_array = vaf_array / vaf_array[:, 0][:, np.newaxis]
        self[self.vaf_columns] = oddsratio_array

    #######################################
    # create baf columns from vaf columns #
    #######################################

    def add_baf(self, fit_self_vafs=False):
        vafs = self.fit_vaf_sum_to_one(inplace=fit_self_vafs)
        bafs = get_baf_from_vaf(vafs)
        for col_idx, baf_idx in enumerate(self.get_baf_indexes()):
            self[baf_idx] = bafs[:, col_idx]

    ##################
    # hetalt marking #
    ##################

    def add_baf_hetalt_flag(self, cutoff=0.15):
        """Args:
            cutoff: If maximum VAF value is greater than (1 - cutoff) value, regarded as homalt
        """
        vaf_array = self.fit_vaf_sum_to_one(inplace=False)
        max_vafs = vaf_array.max(axis=1)
        is_het = (max_vafs < (1 - cutoff))
        self[self.__class__.is_het_colname] = is_het

    def add_baf_hetalt_flag_pairednormal(self, germline_baf_rawdata):
        self.sort()
        germline_baf_rawdata.sort()
        assert (
            self[self.nonannot_columns] 
            == germline_baf_rawdata[germline_baf_rawdata.nonannot_columns] 
            #self.get_coordinate_array()
            #== germline_baf_rawdata.get_coordinate_array()
        ).all()
        
        germline_baf_rawdata.add_baf_hetalt_flag()
        ishet_colname = self.__class__.is_het_colname
        self[ishet_colname] = germline_baf_rawdata[ishet_colname]

    def get_hetalt_selector(self):
        return self[self.__class__.is_het_colname]

    def subset_hetalt(self):
        """Args:
            cutoff: If maximum VAF value is greater than (1 - cutoff) value, regarded as homalt
        """
        return self.loc[self.get_hetalt_selector(), :]

    ##########
    # others #
    ##########

    @property
    def baf_columns(self):
        """baf0 comes first, the baf1, then baf2, ..."""
        result = [x for x in self.columns if re.fullmatch(r'baf[0-9]+', x)]
        result.sort(
            key=(
                lambda x: int(re.sub('^baf', '', x))
            )
        )  # baf0, baf1, baf2, ...
        return result

    def get_baf_array(self):
        """baf0 is on the left"""
        return self.df.loc[:, self.baf_columns].to_numpy()

    ###########################################
    # filtering rows which are invalid as baf #
    ###########################################

    def remove_invalid_vaf_rows(self):
        vaf_array = self[self.vaf_columns]
        invalid_selector = functools.reduce(
            np.logical_or,
            [
                (vaf_array.sum(axis=1) < 0.9), 
                (vaf_array == 0).any(axis=1),
            ],
        )
        return self.loc[~invalid_selector, :]

    def remove_germline_loh_region(self, germline_CN_gdf=None, is_female=None, nproc=1):
        if germline_CN_gdf is None:
            assert is_female is not None
            assert len(self.get_baf_indexes()) == 1  # ploidy == 2
            invalid_selector = (
                cnvcall.get_default_Bg(self, is_female, nproc) 
                == 0
            )
            result = self.loc[~invalid_selector, :]
        else:
            clonal_B_list = [
                germline_CN_gdf.get_clonal_B(baf_index, germline=False) 
                for baf_index in self.get_baf_indexes()
            ]
            invalid_selector = functools.reduce(np.logical_or, [(x == 0) for x in clonal_B_list])
            invalid_CN_gdf = germline_CN_gdf.loc[invalid_selector, :]
            result = self.subtract(invalid_CN_gdf)

        return result

    def keep_equal_germline_allelecopy_region(self, germline_CN_gdf=None, is_female=None, nproc=1):
        if germline_CN_gdf is None:
            assert is_female is not None
            assert len(self.get_baf_indexes()) == 1  # ploidy == 2
            default_Bg = cnvcall.get_default_Bg(self, is_female, nproc)
            default_Ag = cnvcall.get_default_CNg(self, is_female, nproc) - default_Bg
            valid_selector = np.logical_and(
                (default_Bg == default_Ag),
                (default_Bg != 0),
            )
            result = self.loc[valid_selector, :]
        else:
            clonal_CN = germline_CN_gdf.get_clonal_CN(germline=False)
            clonal_B_list = [
                germline_CN_gdf.get_clonal_B(baf_index, germline=False) 
                for baf_index in self.get_baf_indexes()
            ]
            clonal_A = clonal_CN - functools.reduce(operator.add, clonal_B_list)

            valid_selector = functools.reduce(
                np.logical_and, 
                [(clonal_A == x) for x in clonal_B_list],
            )
            valid_selector = np.logical_and(
                valid_selector,
                (clonal_A != 0),
            )
            valid_CN_gdf = germline_CN_gdf.loc[valid_selector, :]
            result = self.intersect(valid_CN_gdf, nproc=nproc)

        return result

    def filter_valid_bafs(
        self, 
        germline_CN_gdf=None, 
        is_female=None, 
        remove_loh=True, 
        skip_hetalt_filter=False,
        nproc=1,
        verbose=False,
    ):
        """Keeps only rows adequate for making BAF segment"""
        #1
        if verbose:
            logutils.log(f'Removing rows where sum of VAFs are too low OR any of VAF is 0')
        result = self.remove_invalid_vaf_rows()
        #2
        if remove_loh:
            if verbose:
                logutils.log(f'Removing regions where any germline B copy number is 0')
            result = result.remove_germline_loh_region(
                germline_CN_gdf=germline_CN_gdf, 
                is_female=is_female, 
                nproc=nproc,
            )
        #3
        if not skip_hetalt_filter:
            if verbose:
                logutils.log(f'Keeping only hetalt positions')
            result = result.subset_hetalt()

        return result

    def get_segment_dict(
        self, 
        target_region=None, 

        germline_baf_rawdata=None,

        germline_CN_gdf=None, 
        is_female=None,
        skip_hetalt_filter=False,

        return_filtered_rawdata=False,

        nproc=1, 
        **kwargs,
    ):
        """* Use "germline_CN_gdf" with germline sample
        * Only use hetalt positions
        """
        result = dict()

        germline_is_given = (germline_baf_rawdata is not None)
        remove_loh = (
            (germline_CN_gdf is not None)
            or (is_female is not None)
        )

        # add hetalt flag
        if not skip_hetalt_filter:
            if germline_is_given:
                self.add_baf_hetalt_flag_pairednormal(germline_baf_rawdata)
            else:
                self.add_baf_hetalt_flag()

        # filter
        filtered_self = self.filter_valid_bafs(
            germline_CN_gdf=germline_CN_gdf, 
            is_female=is_female, 
            remove_loh=remove_loh,
            skip_hetalt_filter=skip_hetalt_filter,
            nproc=nproc,
            verbose=True,
        )

        # main
        for baf_index in get_baf_indexes(self.num_alleles):
            seg_gdf = filtered_self.get_segment(
                annot_colname=baf_index,
                drop_annots=True,
                nproc=nproc,
                **kwargs,
            )
            seg_gdf = BAFSegmentDataFrame.from_frame(seg_gdf.df, refver=self.refver)
            seg_gdf = seg_gdf.fill_gaps(edit_first_last=True)
            seg_gdf.add_rawdata_info_simple(
                filtered_self, 
                baf_index,
                merge_methods=['mean', 'std'],
                nproc=nproc,
            )

            if target_region is not None:
                seg_gdf = seg_gdf.intersect(target_region, nproc=nproc)

            result[baf_index] = seg_gdf

        if return_filtered_rawdata:
            return result, filtered_self
        else:
            return result

    def draw_baf(
        self,
        baf_idx,
        ax=None,
        genomeplotter=None,

        frac=None,
        plotdata=None,

        # fig generation params
        title=None,
        suptitle_kwargs=dict(),
        subplots_kwargs=dict(),

        # drawing kwargs
        plot_kwargs=dict(),

        # axes setting
        setup_axes=True,
        ylabel=False,
        ylabel_prefix='',
        ylabel_kwargs=dict(),
        ymax=False,
        ymin=False,
        yticks=False,
        draw_common_kwargs=dict(),
        rotate_chromlabel=None,

        # multicore plotdata generation
        nproc=1,
        verbose=True,
    ):
        if ylabel is False:
            ylabel = f'BAF ({baf_idx})'
        if ymax is False:
            ymax = 0.52
        if ymin is False:
            ymin = -0.02
        if yticks is False:
            yticks = np.round(np.arange(0, 0.6, 0.1), 1)

        gdraw_result = self.draw_dots(
            y_colname=baf_idx,
            ax=ax,
            genomeplotter=genomeplotter,
            frac=frac,
            plotdata=plotdata,

            title=title,
            suptitle_kwargs=suptitle_kwargs,
            subplots_kwargs=subplots_kwargs,

            plot_kwargs=(self.__class__.default_plot_kwargs | plot_kwargs),

            setup_axes=setup_axes,
            ylabel=ylabel,
            ylabel_prefix=ylabel_prefix,
            ylabel_kwargs=ylabel_kwargs,
            ymax=ymax,
            ymin=ymin,
            yticks=yticks,
            draw_common_kwargs=draw_common_kwargs,
            rotate_chromlabel=rotate_chromlabel,

            nproc=nproc,
            log_suffix=' (BAF raw data)',
            verbose=verbose,
        )
        gdraw_result.set_name(self.__class__.drawresult_name)

        return gdraw_result


@deco.get_deco_atleast1d(['sampleids'])
def get_bafdfs_from_vcf(
    vcf_path, 
    sampleids, 
    n_allele=2,
    nproc=1,
    exclude_other=False,
    prop=None,
    vpfilter=None,
    verbose=True,
):
    vafdf = libvp.get_vafdf(
        vcf_path, 
        sampleids, 
        n_allele=n_allele,
        nproc=nproc,
        exclude_other=exclude_other,
        prop=prop,
        vpfilter=vpfilter,
        verbose=verbose,
    )
    return get_bafdfs_from_vafdf(vafdf)


def get_bafdfs_from_vafdf(vafdf):
    #assert isinstance(vafdf, VCFDataFrame), f'Invalid input object type: {type(vafdf)}'

    vaf_formatkeys = [f'{x}_vaf' for x in vafdf.allele_columns]
    result = dict()
    for sid in vafdf.samples:
        vafs = vafdf.get_format(sid, vaf_formatkeys).to_numpy()
        bafs = get_baf_from_vaf(vafs)

        src_data = dict()
        src_data['chroms'] = vafdf.chromosomes
        src_data['start0s'] = vafdf.starts
        src_data['end0s'] = vafdf.ends
        for colname in vafdf.allele_columns:
            src_data[colname] = vafdf[colname]

        # assign vaf columns
        vafs = vafdf.get_format(sid, vaf_formatkeys).to_numpy()
        for colname, values in zip(vaf_formatkeys, np.swapaxes(vafs, 0, 1)):
            src_data[colname] = values

        # make BAFRawDataFrame
        baf_gdf = BAFRawDataFrame.from_data(refver=vafdf.refver, **src_data)
        # make baf columns
        baf_gdf.add_baf(fit_self_vafs=False)

        # result
        result[sid] = baf_gdf

    return result


#######################
# BAFSegmentDataFrame #
#######################

class BAFSegmentDataFrame(SegmentDataFrame):
    distinfo_suffix = '_distinfo'
    corrected_baf_suffix = '_corrected'
    drawresult_name = 'baf_segment'
    default_plot_kwargs = dict(
        color='tab:blue',
        alpha=genomedf_draw.DEFAULT_SEGMENT_ALPHA,
        #alpha=0.9,
        linewidth=2,
        zorder=5,
    )

    #@property
    #def baf_colname(self):
    #    result = [x for x in self.columns if re.fullmatch(r'baf[0-9]+', x)]
    #    if len(result) != 1:
    #        raise Exception(f'Invalid column names: does not have a unique baf-index-named column.')
    #    return result[0]

    #################################
    # baf_index and colname getters #
    #################################

    def get_baf_index(self):
        baf_strings = set()
        for x in self.annot_cols:
            mat = re.match(r'baf[0-9]+', x)
            if mat is not None:
                baf_strings.add(mat.group(0))
        if len(baf_strings) != 1:
            raise Exception(f'Could not identify baf index from column names.')

        return baf_strings.pop()

    def get_baf_mean_colname(self, baf_index=None):
        if baf_index is None:
            baf_index = self.get_baf_index()
        return baf_index + '_mean'

    def get_baf_mean(self, baf_index=None):
        return self[self.get_baf_mean_colname(baf_index=baf_index)]

    def get_distinfo_colname(self, key, baf_index=None):
        assert key in (
            'ndata', 
            'center', 
            'width',
            'peak_index',
            'height',
            'left_ips',
            'right_ips',
        )

        if baf_index is None:
            baf_index = self.get_baf_index()
        prefix = baf_index + self.__class__.distinfo_suffix
        return prefix + '_' + key

    #def drop_distinfo(self):
    #    cols_to_drop = [x for x in self.columns if self.__class__.distinfo_suffix in x]
    #    return self.drop_annots(cols_to_drop)

    def get_corrected_baf_colname(self, baf_index=None):
        if baf_index is None:
            baf_index = self.get_baf_index()
        return baf_index + self.__class__.corrected_baf_suffix

    ########################
    # column value getters #
    ########################

    def get_dist_center(self, baf_index=None):
        return self[self.get_distinfo_colname('center', baf_index=baf_index)]

    def get_dist_ndata(self, baf_index=None):
        return self[self.get_distinfo_colname('ndata', baf_index=baf_index)]

    def get_corrected_baf(self, baf_index=None):
        return self[self.get_corrected_baf_colname(baf_index=baf_index)]

    #############
    # modifiers #
    #############

    def add_rawdata_info_simple(
        self, 
        baf_rawdata_gdf, 
        baf_index,
        merge_methods=['mean', 'std'],
        nproc=1,
    ):
        joined_gdf = self.drop_annots().join(
            baf_rawdata_gdf,
            right_gdf_cols=baf_index,
            how='left',
            merge=merge_methods,
            nproc=nproc,
        )
        self.assign_frame(joined_gdf.df)

    @deco.get_deco_atleast1d(['distinfo_keys'])
    def add_rawdata_info(
        self, 
        baf_rawdata_gdf, 
        baf_idx,
        distinfo_keys=['ndata', 'center', 'width', 'left_ips', 'right_ips'],
        distinfo_kwargs=dict(),
        merge_methods=['mean', 'std'],
        nproc=1,
    ):
        if distinfo_keys is not None:
            merge_methods = merge_methods + [bafdistinfo_aggfunc]

        joined_gdf = self.drop_annots().join(
            baf_rawdata_gdf,
            right_gdf_cols=baf_idx,
            how='left',
            merge=merge_methods,
            suffixes={
                'bafdistinfo_aggfunc': self.__class__.distinfo_suffix,
            },
            merge_kwargs=distinfo_kwargs,
            nproc=nproc,
        )

        joined_df = joined_gdf.df
        bafdistinfo_key = str(baf_idx) + self.__class__.distinfo_suffix
        bafdistinfo_list = joined_df[bafdistinfo_key]
        for key in distinfo_keys:
            joined_df[f'{bafdistinfo_key}_{key}'] = [
                (
                    x['hetalt'][key] 
                    if x['hetalt'] is not None else
                    np.nan
                )
                for x in bafdistinfo_list
            ]
        joined_df.drop(bafdistinfo_key, axis=1, inplace=True)

        self.assign_frame(joined_df)

    def add_corrected_baf_simple(self, cutoff=0.41):
        centers = self.get_dist_center()
        centers[centers > cutoff] = 0.5
        self[self.get_corrected_baf_colname(baf_index=None)] = centers

    def merge_low_ndata_segments(self, cutoff=30, nproc=1):
        ndata_colname = self.get_distinfo_colname('ndata', baf_index=None)
        return super().merge_low_ndata_segments(
            ndata_colname=ndata_colname, 
            cutoff=cutoff, 
            nproc=nproc,
        )

    def incoporate_low_ndata_segments(self, cutoff=30, nproc=1):
        ndata_colname = self.get_distinfo_colname('ndata', baf_index=None)
        return super().incoporate_low_ndata_segments(
            ndata_colname=ndata_colname, 
            cutoff=cutoff,
            nproc=nproc,
        )

    def get_filter_colname(self):
        baf_index = self.get_baf_index()
        return f'{baf_index}_selected'

    def add_filter(self, width_cutoff=None, center_cutoff=(0.4, None)):
        width_colname = self.get_distinfo_colname('width', baf_index=None)
        center_colname = self.get_distinfo_colname('center', baf_index=None)

        width_values = self[width_colname]
        width_selector = self.filter_helper(width_values, width_cutoff)

        center_values = self[center_colname]
        center_selector = self.filter_helper(center_values, center_cutoff)

        selector = np.logical_and(width_selector, center_selector)
        baf_index = self.get_baf_index()
        self[self.get_filter_colname()] = selector

    def get_filter(self):
        return self[self.get_filter_colname()]

    def draw_baf(
        self,
        ax=None,
        genomeplotter=None,
        plotdata=None,

        # fig generation params
        subplots_kwargs=dict(),

        # drawing kwargs
        plot_kwargs=dict(),

        # axes setting
        setup_axes=True,
        ylabel=False,
        ylabel_prefix='',
        ylabel_kwargs=dict(),
        ymax=False,
        ymin=False,
        yticks=False,
        draw_common_kwargs=dict(),
        rotate_chromlabel=None,
        title=None,
        suptitle_kwargs=dict(),

        # multicore plotdata generation
        nproc=1,
        verbose=True,
    ):
        if ylabel is False:
            ylabel = f'BAF mean ({self.get_baf_index()})'
            #ylabel = f'Corrected BAF ({self.get_baf_index()})'
        if ymax is False:
            ymax = 0.52
        if ymin is False:
            ymin = -0.02
        if yticks is False:
            yticks = np.round(np.arange(0, 0.6, 0.1), 1)

        y_colname_mean = self.get_baf_mean_colname()
        gdraw_result = self.draw_hlines(
            y_colname=y_colname_mean,
            ax=ax,
            genomeplotter=genomeplotter,
            plotdata=plotdata,

            subplots_kwargs=subplots_kwargs,

            plot_kwargs=(self.__class__.default_plot_kwargs | plot_kwargs),

            setup_axes=setup_axes,
            ylabel=ylabel,
            ylabel_prefix=ylabel_prefix,
            ylabel_kwargs=ylabel_kwargs,
            ymax=ymax,
            ymin=ymin,
            yticks=yticks,
            draw_common_kwargs=draw_common_kwargs,
            rotate_chromlabel=rotate_chromlabel,
            title=title,
            suptitle_kwargs=suptitle_kwargs,

            nproc=nproc,
            log_suffix=' (BAF segment)',
            verbose=verbose,
        )
        gdraw_result.set_name(self.__class__.drawresult_name)

        return gdraw_result


#############################################################
# baf distribution analyzer - invalid region identification #
#############################################################

def get_bafdistinfo_maxpeakidx_wotrough(peakresult):
    return np.argmax(peakresult['peak_properties']['peak_heights'])


def get_bafdistinfo_maxpeakidx(peakresult, trough_x):
    relevant_selector = peakresult['peak_xs'] > trough_x
    relevant_peakheights = peakresult['peak_properties']['peak_heights'][relevant_selector]
    max_peakheight = np.max(relevant_peakheights)
    max_peak_idxs = np.nonzero(
        peakresult['peak_properties']['peak_heights'] == max_peakheight
    )[0]
    if len(max_peak_idxs) != 1:
        raise Exception(f'More than one maximal peaks right to the trough point')
    return max_peak_idxs[0]


def get_bafdistinfo_findtrough_old(peakresult):
    return peakresult['xs'][np.argmin(peakresult['ys'])]


def get_bafdistinfo_findtrough(peakresult):
    negative_ys = np.max(peakresult['ys']) - peakresult['ys']
    negative_peaks, negative_peak_properties = scipy.signal.find_peaks(negative_ys)
    if len(negative_peaks) == 0:
        return None
    else:
        return peakresult['xs'][np.min(negative_peaks)]

#BAFPEAKINFO_KEYS = (
#    'peak_index',
#    'center',
#    'height',
#    'left_ips',
#    'right_ips',
#    'ndata',
#)

#class BAFPeakInfo(collections.namedtuple('BAFPeakInfo', BAFPEAKINFO_KEYS)):
#    pass


def get_bafdistinfo_make_peakinfo(peakresult, peak_array_index, reltoheight=None):
    peakinfo_src = dict()

    peakinfo_src['peak_index'] = peak_array_index
    peakinfo_src['center'] = peakresult['peak_xs'][peak_array_index]
    peakinfo_src['height'] = peakresult['peak_ys'][peak_array_index]

    if reltoheight is None:
        left_key = 'left_ips'
        right_key = 'right_ips'
    else:
        left_key = 'left_ips_reltoheight'
        right_key = 'right_ips_reltoheight'

    peakinfo_src['left_ips'] = peakresult['xcoord_values'][left_key][peak_array_index]
    peakinfo_src['right_ips'] = peakresult['xcoord_values'][right_key][peak_array_index]
    peakinfo_src['width'] = peakinfo_src['right_ips'] - peakinfo_src['left_ips']

    peakinfo_src['ndata'] = peakresult['ndata'][peak_array_index]

    #return BAFPeakInfo(**peakinfo_src)
    return peakinfo_src


def check_homalt(peakinfo):
    return (peakinfo['left_ips'] < 0)


def make_result_new(peakresult, reltoheight):
    peakinfo_list = [
        get_bafdistinfo_make_peakinfo(peakresult, idx, reltoheight)
        for idx in range(peakresult.npeaks)
    ]
    homalt_selector = np.fromiter(map(check_homalt, peakinfo_list), dtype=bool)
    homalt_indexes = np.nonzero(homalt_selector)[0]
    hetalt_indexes = np.nonzero(~homalt_selector)[0]

    # homalt
    if len(homalt_indexes) == 0:
        #result['homalt'] = None
        homalt_peakinfo = None
    else:
        homalt_peakinfo = peakinfo_list[homalt_indexes[0]]

    # hetalt
    if len(hetalt_indexes) == 0:
        hetalt_peakinfo = None
    else:
        hetalt_candidates = [peakinfo_list[x] for x in hetalt_indexes]
        maxheight_arg = np.argmax([x['height'] for x in hetalt_candidates])
        hetalt_peakinfo = hetalt_candidates[maxheight_arg]

    return {'homalt': homalt_peakinfo, 'hetalt': hetalt_peakinfo}


def get_bafdistinfo(
    bafs, 

    decal=False,

    #bw_method=None,
    bw_method=0.05,
    xs=None,
    bins=None,
    reltoheight=None, 

    hist=False,
    return_peakresult=False,

    **kwargs,
):
    if decal:
        bafs = decal_baf(bafs)

    if xs is None:
        if decal:
            xs = np.arange(-0.05, 1.05, 0.001)
        else:
            xs = np.arange(-0.05, 0.55, 0.001)
    if bins is None:
        step = 0.01
        if decal:
            bins = np.arange(0, 1 + step, step)
        else:
            bins = np.arange(0, 0.5 + step, step)

    if hist:
        peakresult = peakutils.find_hist_peaks(
            bafs, 
            reltoheight=reltoheight, 
            bins=bins,
            **kwargs,
        )
    else:
        peakresult = peakutils.find_density_peaks(
            bafs, 
            reltoheight=reltoheight, 
            xs=xs,
            bw_method=bw_method,
            **kwargs,
        )

    if decal:
        peakresult.undecal(0.5)
    #peakresult = peakresult.remove_contained_peaks()
    bafdistinfo = make_result_new(peakresult, reltoheight)

    if return_peakresult:
        return bafdistinfo, peakresult
    else:
        return bafdistinfo


def bafdistinfo_aggfunc(
    bafs, 
    #keys=['ndata', 'center'],
    **kwargs,
):
    #assert set(keys).issubset(BAFPEAKINFO_KEYS)

    #assert bafs.ndim == 2
    #assert bafs.shape[1] == 1
    #bafs = np.asarray(bafs)[:, 0]
    bafs = np.atleast_1d(np.squeeze(np.asarray(bafs)))

    try:
        bafdistinfo = get_bafdistinfo(
            bafs, hist=False, return_peakresult=False, **kwargs
        )
    except (DensityGenerationFailure, NoPeakError):
        try:
            bafdistinfo = get_bafdistinfo(
                bafs, hist=True, return_peakresult=False, **kwargs
            )
        except NoPeakError as exc:
            raise Exception(f'Even histogram could not detect a peak from the baf data.') from exc

    return bafdistinfo


def draw_bafarray_breadth(
    bafs,

    ax=None,
    distinfo_kwargs=dict(), 
    draw_peak_kwargs=dict(),

    hist=False,

    return_data=False,
):
    # prepare peakresult
    hist_bafdistinfo, hist_peakresult = get_bafdistinfo(
        bafs, 
        hist=True,
        return_peakresult=True,
        **distinfo_kwargs,
    )
    try:
        density_bafdistinfo, density_peakresult = get_bafdistinfo(
            bafs, 
            hist=False,
            return_peakresult=True,
            **distinfo_kwargs,
        )
    except DensityGenerationFailure:
        logutils.log(f'Density generation failed')
        density_bafdistinfo = None
        density_peakresult = None

    # prepare distinfo and returned peakresult
    if hist:
        bafdistinfo = hist_bafdistinfo 
        peakresult = hist_peakresult
    else:
        if density_bafdistinfo is None:
            bafdistinfo = hist_bafdistinfo 
            peakresult = hist_peakresult
        else:
            bafdistinfo = density_bafdistinfo 
            peakresult = density_peakresult

    # draw plot
    fig, ax = peakutils.draw_peaks(
        bafs, 
        ax=ax,
        histpeaks=hist_peakresult, 
        densitypeaks=density_peakresult, 
        omit_density=(density_peakresult is None),
        **draw_peak_kwargs,
    )

    if bafdistinfo['homalt'] is not None:
        ax.axvspan(
            bafdistinfo['homalt']['left_ips'], 
            bafdistinfo['homalt']['right_ips'], 
            color='magenta',
            alpha=0.2,
        )
        ax.text(bafdistinfo['homalt']['center'], ax.get_ylim()[1] * 1.05, 'hom')
    if bafdistinfo['hetalt'] is not None:
        ax.axvspan(
            bafdistinfo['hetalt']['left_ips'], 
            bafdistinfo['hetalt']['right_ips'], 
            color='cyan',
            alpha=0.2,
        )
        ax.text(bafdistinfo['hetalt']['center'], ax.get_ylim()[1] * 1.05, 'hom')

    ax.set_xlabel('BAF')
    ax.set_ylabel('density')

    if return_data:
        return fig, ax, bafdistinfo, peakresult
    else:
        return fig, ax


# unused ones

def get_bafdistinfo_base_usetrough(
    bafs, 

    #homalt_cutoff=None, 
    bw_method,
    xs,
    bins,
    reltoheight,

    hist=False,
    return_peakresult=False,
    use_trough=False,

    #homaltpeak_cutoff=0.05,
    **kwargs,
):
    #bafs = np.asarray(bafs)
    #if homalt_cutoff is not None:
    #    bafs = bafs[bafs > homalt_cutoff]

    if hist:
        peakresult = peakutils.find_hist_peaks(
            bafs, 
            reltoheight=reltoheight, 
            bins=bins,
            **kwargs,
        )
    else:
        peakresult = peakutils.find_density_peaks(
            bafs, 
            reltoheight=reltoheight, 
            xs=xs,
            bw_method=bw_method,
            **kwargs,
        )

    # remove contained peaks
    peakresult = peakresult.remove_contained_peaks()

    trough_x = get_bafdistinfo_findtrough(peakresult)
    if (trough_x is not None) and use_trough:
        result = make_result_using_trough(peakresult, reltoheight, trough_x)
    else:
        result = make_result_without_trough(peakresult, reltoheight)
    
    if return_peakresult:
        return result, peakresult
    else:
        return result


def make_result_using_trough(peakresult, reltoheight, trough_x):
    peaks_lt_trough = peakresult.subset(peakresult['peak_xs'] < trough_x)
    peaks_gt_trough = peakresult.subset(peakresult['peak_xs'] > trough_x)
    if (peaks_lt_trough.npeaks == 0) or (peaks_gt_trough.npeaks == 0):
        raise Exception(f'No peaks to the left or right of trough point')

    left_max_peak_idx = peaks_lt_trough.get_highest_peak_index()
    left_max_peak_array_idx = list(peakresult['peak_indexes']).index(left_max_peak_idx)
    left_max_peakinfo = get_bafdistinfo_make_peakinfo(peakresult, left_max_peak_array_idx, reltoheight)

    right_max_peak_idx = peaks_gt_trough.get_highest_peak_index()
    right_max_peak_array_idx = list(peakresult['peak_indexes']).index(right_max_peak_idx)
    right_max_peakinfo = get_bafdistinfo_make_peakinfo(peakresult, right_max_peak_array_idx, reltoheight)

    result = {
        'homalt': left_max_peakinfo,
        'hetalt': right_max_peakinfo,
    }
    return result


def make_result_without_trough(peakresult, reltoheight):
    result = dict()

    if peakresult.npeaks == 1:
        peakinfo = get_bafdistinfo_make_peakinfo(peakresult, 0, reltoheight)
        if check_homalt(peakinfo):
            result['homalt'] = peakinfo
            result['hetalt'] = None
        else:
            result['homalt'] = None
            result['hetalt'] = peakinfo
    else:
        max_peak_indexes = (
            np.argsort(peakresult['peak_properties']['peak_heights'])
        )[::-1][:2]  # (largest_idx, second_largest_idx)
        first_peakinfo = get_bafdistinfo_make_peakinfo(peakresult, max_peak_indexes[0], reltoheight)
        second_peakinfo = get_bafdistinfo_make_peakinfo(peakresult, max_peak_indexes[1], reltoheight)
        first_ishomalt = check_homalt(first_peakinfo)
        second_ishomalt = check_homalt(second_peakinfo)
        if first_ishomalt == second_ishomalt:
            raise FailedNoTroughMethod(
                f'Two maximal peaks are both hom or both het'
            )
        if first_ishomalt:
            result['homalt'] = first_peakinfo
            result['hetalt'] = second_peakinfo
        else:
            result['homalt'] = second_peakinfo
            result['hetalt'] = first_peakinfo

    return result

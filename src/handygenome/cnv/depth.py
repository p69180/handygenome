import multiprocessing
import functools

import pandas as pd
import numpy as np
import pyranges as pr
import matplotlib as mpl
import matplotlib.pyplot as plt
import pysam

import handygenome.peakutils as peakutils
import handygenome.logutils as logutils
from handygenome.genomedf.genomedf import GenomeDataFrame, SegmentDataFrame
import handygenome.cnv.cncall as cncall
from handygenome.cnv.cncall import DefaultCNgUnavailableError
import handygenome.plot.misc as plotmisc
from handygenome.peakutils import NoPeakError


def make_normalized_depth(depths, lengths):
    valid_idxs = ~np.isnan(depths)
    valid_depths = depths[valid_idxs]
    valid_lenghts = lengths[valid_idxs]
    average_depth = np.average(valid_depths, weights=valid_lenghts)
    result = depths / average_depth
    return result


class DepthRawDataFrame(GenomeDataFrame):
    depth_colname = 'depth'
    norm_depth_colname = 'normalized_depth'
    corrected_norm_depth_colname = 'corrected_normalized_depth'
    MQ_colname = 'MQ'
    filter_colname = 'depth_selected'

    DEFAULT_DTYPES = {
        'Chromosome': 'string', 
        'Start': int, 
        'End': int,
        depth_colname: float,
        norm_depth_colname: float,
    }

    def init_sanitycheck(self):
        pass
        #assert self.__class__.depth_colname in self.annot_cols

    def spawn(self, frame):
        return super().spawn(frame, use_median=self.use_median)

    def __init__(self, refver, use_median=False):
        super().__init__(refver)
        self.use_median = use_median

    @classmethod
    def from_frame(cls, frame, refver, dtype=dict(), use_median=False):
        result = cls(refver, use_median=use_median)
        result.assign_frame(frame, dtype)
        result.init_sanitycheck()
        return result

    @classmethod
    def read_tsv(cls, filename, refver, use_median=False):
        """Column names must be set appropriately"""
        frame = cls.make_frame_from_tsv(filename)
        return cls.from_frame(frame, refver, use_median=use_median)

    @classmethod
    def load_mosdepth(cls, filename, refver, use_median=False):
        assert filename.endswith('.regions.bed.gz')
        frame = cls.make_frame_from_tsv(filename, annot_cols=[cls.depth_colname])
        return cls.from_frame(frame, refver, use_median=use_median)

    def get_colname_ylabel(self, depthtype):
        assert depthtype in ('norm', 'normalized', 'raw', 'corr', 'corrected')
        if depthtype in ('norm', 'normalized'):
            default_ylabel = f'normalized depth'
            y_colname = self.__class__.norm_depth_colname
        elif depthtype in ('raw',):
            default_ylabel = f'raw depth'
            y_colname = self.__class__.depth_colname
        elif depthtype in ('corr', 'corrected'):
            default_ylabel = f'corrected depth'
            y_colname = self.__class__.corrected_norm_depth_colname

        return default_ylabel, y_colname

    #################
    # normalization #
    #################

    def set_normalized_depth(self):
        self[self.__class__.norm_depth_colname] = make_normalized_depth(self.depth, self.lengths)

    ################
    # segmentation #
    ################

    def get_segment(self, *args, **kwargs):
        seg = super().get_segment(*args, **kwargs)
        return DepthSegmentDataFrame.from_frame(seg.df, refver=self.refver)

    ##############
    # properties #
    ##############

    @property
    def depth(self):
        return self[self.__class__.depth_colname]

    raw_depth = depth

    @property
    def normalized_depth(self):
        return self[self.__class__.norm_depth_colname]

    norm_depth = normalized_depth

    @property
    def corrected_depth(self):
        return self[self.__class__.corrected_norm_depth_colname]

    corr_depth = corrected_depth

    @property
    def MQ(self):
        return self[self.__class__.MQ_colname]

    ###############
    # MQ addition #
    ###############

    def get_MQ_fetchregion_gdf(self, readlength=151):
        self.sort()

        new_chroms = list()
        new_start0s = list()
        new_end0s = list()
        halflen = int(readlength / 2)

        for chrom, subgdf in self.group_bychrom(sort=False).items():
            start0s = subgdf.start0s
            end0s = subgdf.end0s
            lengths = subgdf.lengths

            long_length_selector = (lengths > readlength)
            start0s[long_length_selector] = start0s[long_length_selector] + halflen
            end0s[long_length_selector] = end0s[long_length_selector] - halflen

            new_start0s.append(start0s)
            new_end0s.append(end0s)
            new_chroms.append(np.repeat(chrom, len(start0s)))

        result_chroms = np.concatenate(new_chroms)
        result_start0s = np.concatenate(new_start0s)
        result_end0s = np.concatenate(new_end0s)

        return self.spawn_from_data(
            chroms=result_chroms, 
            start0s=result_start0s, 
            end0s=result_end0s,
        )

    def add_MQ(self, bam_path, readlength=151, nproc=1, n_split=1000, verbose=False):
        # get MQ values
        fetchregion_gdf = self.get_MQ_fetchregion_gdf(readlength)  # self is sorted in this step
        split_fetchregions = fetchregion_gdf.equal_nrow_split(n=n_split)
        args = (
            (bam_path, gdf.chroms, gdf.start0s, gdf.end0s, verbose)
            for gdf in split_fetchregions
        )
        with multiprocessing.Pool(nproc) as pool:
            map_result = pool.starmap(make_MQ_array, args)

        # apply results
        self[self.__class__.MQ_colname] = np.concatenate(map_result)

    def check_has_MQ(self):
        return self.__class__.MQ_colname in self.columns

    #############
    # filtering #
    #############

    def add_filter(self, MQ_cutoff=(50, None), norm_depth_cutoff=(0, None), raw_depth_cutoff=None):
        MQ_selector = self.filter_helper(
            self[self.__class__.MQ_colname], 
            MQ_cutoff,
            include=False,
        )
        norm_depth_selector = self.filter_helper(
            self[self.__class__.norm_depth_colname], 
            norm_depth_cutoff,
            include=False,
        )
        raw_depth_selector = self.filter_helper(
            self[self.__class__.depth_colname], 
            raw_depth_cutoff,
            include=False,
        )
        self[self.__class__.filter_colname] = functools.reduce(
            np.logical_and,
            [MQ_selector, norm_depth_selector, raw_depth_selector]
        )

    def get_filter(self):
        return self[self.__class__.filter_colname]

    ###########
    # drawing #
    ###########

    def draw_depth(
        self,
        ax=None,
        genomeplotter=None,

        depthtype='norm',
        frac=None,
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
        log_suffix=None,
    ):
        # parameter setup by depth type 
        default_ylabel, y_colname = self.get_colname_ylabel(depthtype)

        # y axis params
        if ylabel is False:
            ylabel = default_ylabel
        if (ymax is False) and (depthtype != 'raw'):
            ymax = 2
        if (ymin is False) and (depthtype != 'raw'):
            ymin = 0

        fig, ax, genomeplotter = self.draw_dots(
            y_colname=y_colname,
            ax=ax,
            genomeplotter=genomeplotter,
            frac=frac,

            title=title,
            suptitle_kwargs=suptitle_kwargs,
            subplots_kwargs=subplots_kwargs,

            plot_kwargs=plot_kwargs,

            setup_axes=setup_axes,
            ylabel=ylabel,
            ylabel_prefix=ylabel_prefix,
            ylabel_kwargs=ylabel_kwargs,
            ymax=ymax,
            ymin=ymin,
            yticks=yticks,
            draw_common_kwargs=draw_common_kwargs,
            rotate_chromlabel=rotate_chromlabel,

            plotdata=plotdata,
            nproc=nproc,
            log_suffix=' (Depth raw data)',
        )

        return fig, ax

    def draw_MQ(
        self,
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
        log_suffix=None,
    ):
        y_colname = self.__class__.MQ_colname

        # y axis params
        if ylabel is False:
            ylabel = 'mapping quality'
        if ymax is False:
            ymax = 72
        if ymin is False:
            ymin = -2
        if yticks is False:
            yticks = np.arange(0, 80, 10).astype(int)

        fig, ax, genomeplotter = self.draw_dots(
            y_colname=y_colname,
            ax=ax,
            genomeplotter=genomeplotter,
            frac=frac,

            title=title,
            suptitle_kwargs=suptitle_kwargs,
            subplots_kwargs=subplots_kwargs,

            plot_kwargs=plot_kwargs,

            setup_axes=setup_axes,
            ylabel=ylabel,
            ylabel_prefix=ylabel_prefix,
            ylabel_kwargs=ylabel_kwargs,
            ymax=ymax,
            ymin=ymin,
            yticks=yticks,
            draw_common_kwargs=draw_common_kwargs,
            rotate_chromlabel=rotate_chromlabel,

            plotdata=plotdata,
            nproc=nproc,
            log_suffix=' (Mapping quality raw data)',
        )

        return fig, ax


def make_MQ_array(bam_path, chroms, start0s, end0s, verbose=False):
    if verbose:
        logutils.log(f'Beginning region chrom={chroms[0]}, start0={start0s[0]}, end0={end0s[0]}')
    with pysam.AlignmentFile(bam_path) as bam:
        def iterator():
            for c, s, e in zip(chroms, start0s, end0s):
                MQs = tuple(read.mapping_quality for read in bam.fetch(c, s, e))
                if len(MQs) == 0:
                    yield np.nan
                else:
                    yield np.mean(MQs)

        return np.fromiter(iterator(), dtype=float)


class DepthSegmentDataFrame(SegmentDataFrame):
    norm_depth_mean_colname = DepthRawDataFrame.norm_depth_colname + '_mean'
    norm_depth_std_colname = DepthRawDataFrame.norm_depth_colname + '_std'

    corrected_norm_depth_mean_colname = DepthRawDataFrame.corrected_norm_depth_colname + '_mean'
    corrected_norm_depth_std_colname = DepthRawDataFrame.corrected_norm_depth_colname + '_std'

    MQ_mean_colname = DepthRawDataFrame.MQ_colname + '_mean'
    filter_colname = 'depth_selected'

    @property
    def norm_depth_mean(self):
        return self[self.__class__.norm_depth_mean_colname]

    @property
    def norm_depth_std(self):
        return self[self.__class__.norm_depth_std_colname]

    @property
    def corrected_depth_mean(self):
        return self[self.__class__.corrected_norm_depth_mean_colname]

    def find_onecopy_depth(self, is_female, ploidy):
        try:
            default_CNg = cncall.get_default_CNg(self, is_female)
        except DefaultCNgUnavailableError:
            data = self.norm_depth_mean
            weights = self.lengths
        else:
            valid_selector = np.logical_and(
                (default_CNg != 0),
                ~np.isnan(default_CNg),
            )

            correction_factor = ploidy / default_CNg[valid_selector]
            data = self.norm_depth_mean[valid_selector] * correction_factor

            weights = self.lengths[valid_selector]

        # select copy-neutral depth segments
        peakresult = peakutils.find_density_peaks(
            data,
            weights=weights,
            xs=np.arange(0, 2, 0.01),
            bw_method=0.005,
        )
        modal_depth = peakresult['peak_xs'][np.argmax(peakresult['peak_ys'])]
        onecopy_depth = modal_depth / ploidy

        return onecopy_depth

    def add_rawdata_info(
        self, 
        depth_rawdata_gdf, 
        merge_methods=['mean', 'std'],
        rawdepth=False,
        nproc=1,
    ):
        #depth_colname = (
        #    DepthRawDataFrame.depth_colname
        #    if rawdepth else
        #    DepthRawDataFrame.norm_depth_colname
        #)
        #right_gdf_cols = [depth_colname]

        #MQ_colname = DepthRawDataFrame.MQ_colname
        #if MQ_colname in depth_rawdata_gdf.columns:
        #    right_gdf_cols.append(MQ_colname)

        joined_gdf = self.drop_annots().join(
            depth_rawdata_gdf,
            #right_gdf_cols=right_gdf_cols,
            how='left',
            merge=merge_methods,
            winsorize=(0.05, 0.05),
            nproc=nproc,
        )
        self.assign_frame(joined_gdf.df)

    def add_filter(self, MQ_cutoff=(50, None), norm_depth_cutoff=(0, None)):
        MQ_selector = self.filter_helper(
            self[self.__class__.MQ_mean_colname], 
            MQ_cutoff,
            include=False,
        )
        depth_selector = self.filter_helper(
            self[self.__class__.norm_depth_mean_colname], 
            norm_depth_cutoff,
            include=False,
        )
        self[self.__class__.filter_colname] = np.logical_and(MQ_selector, depth_selector)

    def get_filter(self):
        return self[self.__class__.filter_colname]

    ########
    # draw #
    ########

    def get_colname_ylabel(self, depthtype):
        assert depthtype in ('norm', 'normalized', 'corr', 'corrected')
        if depthtype in ('norm', 'normalized'):
            default_ylabel = f'normalized depth'
            y_colname = self.__class__.norm_depth_mean_colname
        elif depthtype in ('corr', 'corrected'):
            default_ylabel = f'corrected depth'
            y_colname = self.__class__.corrected_norm_depth_mean_colname

        return default_ylabel, y_colname

    def draw_depth(
        self,
        ax=None,
        genomeplotter=None,

        depthtype='norm',
        plotdata=None,

        chromwise_peaks=False,
        chromwise_peaks_kwargs=dict(),
        chromwise_peaks_fullspan=True,

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
        log_suffix=None,
    ):
        # parameter setup by depth type 
        default_ylabel, y_colname = self.get_colname_ylabel(depthtype)

        # y axis params
        if ylabel is False:
            ylabel = default_ylabel
        if (ymax is False) and (depthtype != 'raw'):
            ymax = 2
        if (ymin is False) and (depthtype != 'raw'):
            ymin = 0

        #default_ymin, default_ymax = plotmisc.get_boxplot_range(plotdata[y_colname])
        #if ymax is False:
        #    ymax = default_ymax
            #ymax = np.nanquantile(self[y_colname], 0.99)
        #if ymin is False:
        #    ymin = default_ymin
            #ymin = -ymax * 0.01

        #if yticks is False:
        #    yticks = np.round(np.linspace(0, ymax, 10), 2)

        fig, ax, genomeplotter = self.draw_hlines(
            y_colname=y_colname,
            ax=ax,
            genomeplotter=genomeplotter,

            title=title,
            suptitle_kwargs=suptitle_kwargs,
            subplots_kwargs=subplots_kwargs,

            plot_kwargs=plot_kwargs,

            setup_axes=setup_axes,
            ylabel=ylabel,
            ylabel_prefix=ylabel_prefix,
            ylabel_kwargs=ylabel_kwargs,
            ymax=ymax,
            ymin=ymin,
            yticks=yticks,
            draw_common_kwargs=draw_common_kwargs,
            rotate_chromlabel=rotate_chromlabel,

            plotdata=plotdata,
            nproc=nproc,
            log_suffix=' (Depth segment)',
        )

        if chromwise_peaks:
            densitypeak_kwargs = (
                dict(
                    limit=ax.get_ylim(),
                    prominence=(1, None), 
                    bw_method=0.05,
                    xs=np.arange(0, ax.get_ylim()[1], 0.01),
                ) | chromwise_peaks_kwargs
            )
            plotregion_chroms = set(genomeplotter.region_gdf.chroms)
            peakline_kwargs = dict(
                alpha=0.5,
                linewidth=0.8,
                color='orange',
            )

            # calculate peakresults
            peakresult_dict = dict()
            for chrom, sub_self in self.group_bychrom(sort=False).items():
                if chrom not in plotregion_chroms:
                    continue
                if sub_self.nrow == 0:
                    continue

                try:
                    peakresult = peakutils.find_density_peaks(
                        data=sub_self[y_colname],
                        weights=sub_self.lengths,
                        **densitypeak_kwargs,
                    )
                except NoPeakError:
                    continue
                else:
                    peakresult_dict[chrom] = peakresult

            for chrom, peakresult in peakresult_dict.items():
                chrom_plotregion = genomeplotter.region_gdf.subset_chroms(chrom)
                xmin = genomeplotter.genomic_to_plot(chrom, chrom_plotregion.start0s.min())
                xmax = genomeplotter.genomic_to_plot(chrom, chrom_plotregion.end0s.max() - 1)
                xmid = (xmin + xmax) / 2

                if chromwise_peaks_fullspan:
                    for y in peakresult['peak_xs']:
                        ax.axhline(y, **peakline_kwargs)
                else:
                    ax.hlines(peakresult['peak_xs'], xmin, xmax, **peakline_kwargs)

                for y in peakresult['peak_xs']:
                    ax.text(
                        xmid, y, str(round(y, 3)), 
                        ha='center', 
                        va='center', 
                        size=10, 
                        alpha=1,
                        color='aqua',
                        fontweight='bold',
                    )

        return fig, ax

    def draw_MQ(
        self,
        ax=None,
        genomeplotter=None,

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
        log_suffix=None,
    ):
        y_colname = self.__class__.MQ_mean_colname

        # y axis params
        if ylabel is False:
            ylabel = 'mapping quality'
        if ymax is False:
            ymax = 72
        if ymin is False:
            ymin = -2
        if yticks is False:
            yticks = np.arange(0, 80, 10).astype(int)

        fig, ax, genomeplotter = self.draw_hlines(
            y_colname=y_colname,
            ax=ax,
            genomeplotter=genomeplotter,

            title=title,
            suptitle_kwargs=suptitle_kwargs,
            subplots_kwargs=subplots_kwargs,

            plot_kwargs=plot_kwargs,

            setup_axes=setup_axes,
            ylabel=ylabel,
            ylabel_prefix=ylabel_prefix,
            ylabel_kwargs=ylabel_kwargs,
            ymax=ymax,
            ymin=ymin,
            yticks=yticks,
            draw_common_kwargs=draw_common_kwargs,
            rotate_chromlabel=rotate_chromlabel,

            plotdata=plotdata,
            nproc=nproc,
            log_suffix=' (Mapping quality segment)',
        )

        return fig, ax

    def draw_peaks(
        self,
        ax=None,
        depthtype='norm',

        genomeplotter=None,

        density_kwargs=dict(),
        hist_kwargs=dict(),
        density_peakline=True,

        return_peaks=False,
    ):
        if ax is None:
            fig, ax = plt.subplots()
        else:   
            fig = ax.figure

        ax_ylims = ax.get_ylim()

        # parameter setup by depth type 
        default_ylabel, y_colname = self.get_colname_ylabel(depthtype)

        # make common_kwargs
        if genomeplotter is None:
            sub_self = self
        else:
            sub_self = genomeplotter.get_intersect(self)
        common_kwargs = dict(
            data=sub_self[y_colname],
            weights=sub_self.lengths,
            limit=ax_ylims,
        )
        # density & hist kwargs
        density_kwargs = (
            dict(
                height=(0.1, None), 
                #prominence=(0.5, None), 
                bw_method=None,
                xs=np.arange(0, ax_ylims[1], 0.01),
            ) | density_kwargs
        )
        hist_kwargs = (
            dict(
                bins=np.arange(0, ax_ylims[1], 0.01),
            ) | hist_kwargs
        )

        densitypeaks = peakutils.find_density_peaks(
            **(common_kwargs | density_kwargs)
        )
        histpeaks = peakutils.find_hist_peaks(
            **(common_kwargs | hist_kwargs)
        )

        peakutils.draw_peaks(
            ax=ax,
            histpeaks=histpeaks,
            densitypeaks=densitypeaks,
            omit_hist_peakline=density_peakline,
            omit_density_peakline=(not density_peakline),
            lying=True,
        )

        if return_peaks:
            return densitypeaks, histpeaks

    def draw_depth_with_peaks(
        self,
        axd=None,
        depthtype='norm',
        genomeplotter=None,
        nproc=1,
        draw_hlines=True,
        ymin=False,
        ymax=False,
        depth_kwargs=dict(),
        peaks_kwargs=dict(),
    ):
        if axd is None:
            fig, axd = plt.subplot_mosaic(
                [['depth', 'depth_peaks']],
                figsize=(30, 6),
                gridspec_kw={
                    'width_ratios': (1, 0.1),
                    'wspace': 0.02,
                },
            )
        else:
            fig = next(iter(axd.values())).figure

        self.draw_depth(
            ax=axd['depth'],
            depthtype=depthtype,
            genomeplotter=genomeplotter,
            nproc=nproc,
            **(
                depth_kwargs
                | dict(
                    ymin=ymin,
                    ymax=ymax,
                )
            )
        )
        axd['depth_peaks'].set_ylim(axd['depth'].get_ylim())
        densitypeaks, histpeaks = self.draw_peaks(
            ax=axd['depth_peaks'],
            depthtype=depthtype,
            genomeplotter=genomeplotter,
            return_peaks=True,
            **peaks_kwargs,
        )
        axd['depth_peaks'].set_yticks([])
        if draw_hlines:
            for y in densitypeaks['peak_xs']:
                axd['depth'].axhline(
                    y, 
                    color='orange',
                    linewidth=1,
                )

        return fig, axd
        
    

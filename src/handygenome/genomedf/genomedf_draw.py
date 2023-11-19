import re
import os
import datetime
import inspect
import itertools
import contextlib
import multiprocessing
import functools
import warnings
import collections

import pandas as pd
import numpy as np
import pyranges as pr
import matplotlib as mpl
import matplotlib.pyplot as plt

#import handygenome.deco as deco
import handygenome.refgenome.refgenome as refgenome
import handygenome.tools as tools
#import handygenome.logutils as logutils

from handygenome.genomedf.genomedf_base import GenomeDataFrameBase

import handygenome.plot.misc as plotmisc
from handygenome.plot.gui import IntervalSelector
import handygenome.plot.genomeplot as libgenomeplot
from handygenome.plot.genomeplot import GenomePlotter


DEFAULT_SEGMENT_LINEWIDTH = 1.5
DEFAULT_SEGMENT_ALPHA = 0.9


def draw_axessetup(
    ax,
    genomeplotter,

    # y axes setting
    ylabel=None,
    ylabel_prefix=None,
    ylabel_kwargs=dict(),
    ymax=None,
    ymin=None,
    yticks=None,

    # draw common
    draw_common_kwargs=dict(),
    rotate_chromlabel=None,

    # figure suptitle
    fig=None,
    title=None,
    suptitle_kwargs=dict(),
):
    # ylabel
    if ylabel is not None:
        ylabel = ylabel_prefix + ylabel
        ax.set_ylabel(ylabel, **ylabel_kwargs)

    # yticks
    if yticks is False:
        yticks = None
    if yticks is not None:
        ax.set_yticks(yticks)
        ax.set_yticklabels(
            yticks, 
            size=plotmisc.get_yticklabel_size(len(yticks)),
        )

    # ymin, ymax
    if ymin is False:
        ymin = None
    if ymax is False:
        ymax = None
    ax.set_ylim(ymin, ymax)

    # draw_common
    draw_common_kwargs = (
        {'n_xlabel': None}
        | draw_common_kwargs
    )

    if rotate_chromlabel is not None:
        draw_common_kwargs.setdefault('chromlabel_kwargs', dict())
        draw_common_kwargs['chromlabel_kwargs'] |= {'rotation': rotate_chromlabel}

    draw_common_artists = genomeplotter.draw_common(ax, **draw_common_kwargs)

    # figure suptitle
    if title is not None:
        plotmisc.draw_suptitle(fig, title, **suptitle_kwargs)

    return draw_common_artists


class GenomeDataFrameDrawingBase:
    def get_default_genomeplotter(self):
        self_bychrom = self.group_bychrom(sort=False)
        kwargs = {'chroms': list(), 'start0s': list(), 'end0s': list(), 'refver': self.refver}
        for chrom, subgdf in self_bychrom.items():
            kwargs['chroms'].append(chrom)
            kwargs['start0s'].append(subgdf.start0s.min())
            kwargs['end0s'].append(subgdf.end0s.max())
        return GenomePlotter(**kwargs)

    ###################
    # drawing helpers #
    ###################
    
    def draw_preprocess(
        self, 
        y_colname,
        ax,
        genomeplotter,
        frac,

        # fig generation params
        subplots_kwargs,

        # logging
        log_suffix,
        nproc,
        verbose,

        # pre-generated plotdata
        plotdata,
    ):

        # fig/ax handling
        ax_is_given = (ax is not None)
        if ax_is_given:
            fig = ax.figure
        else:
            subplots_kwargs = (
                {'figsize': (20, 8)}
                | subplots_kwargs
            )
            fig, ax = plt.subplots(1, 1, **subplots_kwargs)

        # prepare genomeplotter
        if genomeplotter is None:
            genomeplotter = self.get_default_genomeplotter()

        assert refgenome.compare_refvers(self.refver, genomeplotter.refver)

        # prepare plotdata
        if plotdata is None:
            if y_colname is None:
                data = self
            else:
                data = self.choose_annots(y_colname)

            if frac is not None:
                data = data.sample(frac=frac)
            plotdata = genomeplotter.make_plotdata(
                data, 
                log_suffix=log_suffix,
                nproc=nproc,
                verbose=verbose,
            )

        return fig, ax, genomeplotter, plotdata

#    def draw_axessetup(
#        self, 
#        ax,
#        genomeplotter,
#
#        # y axes setting
#        ylabel,
#        ylabel_prefix,
#        ylabel_kwargs,
#        ymax,
#        ymin,
#        yticks,
#
#        # draw common
#        draw_common_kwargs,
#        rotate_chromlabel,
#
#        # figure suptitle
#        fig,
#        title,
#        suptitle_kwargs,
#    ):
#        # ylabel
#        if ylabel is not None:
#            ylabel = ylabel_prefix + ylabel
#            ax.set_ylabel(ylabel, **ylabel_kwargs)
#
#        # yticks
#        if yticks is False:
#            yticks = None
#        if yticks is not None:
#            ax.set_yticks(yticks)
#            ax.set_yticklabels(
#                yticks, 
#                size=plotmisc.get_yticklabel_size(len(yticks)),
#            )
#
#        # ymin, ymax
#        if ymin is False:
#            ymin = None
#        if ymax is False:
#            ymax = None
#        ax.set_ylim(ymin, ymax)
#
#        # draw_common
#        draw_common_kwargs = (
#            {'n_xlabel': None}
#            | draw_common_kwargs
#        )
#
#        if rotate_chromlabel is not None:
#            draw_common_kwargs.setdefault('chromlabel_kwargs', dict())
#            draw_common_kwargs['chromlabel_kwargs'] |= {'rotation': rotate_chromlabel}
#
#        draw_common_artists = genomeplotter.draw_common(ax, **draw_common_kwargs)
#
#        # figure suptitle
#        if title is not None:
#            plotmisc.draw_suptitle(fig, title, **suptitle_kwargs)
#
#        return draw_common_artists

    ################
    # main drawers #
    ################

    def draw_dots(
        self, 
        y_colname,
        ax=None, 
        genomeplotter=None,
        frac=None,

        # fig generation params
        title=None,
        suptitle_kwargs=dict(),
        subplots_kwargs=dict(),

        # drawing kwargs
        plot_kwargs=dict(),

        # axes setting
        setup_axes=True,
        ylabel=None,
        ylabel_prefix='',
        ylabel_kwargs=dict(),
        ymax=False,
        ymin=False,
        yticks=None,
        draw_common_kwargs=dict(),
        rotate_chromlabel=None,

        # pre-generated plotdata
        plotdata=None,

        # multicore plotdata generation
        nproc=1,
        log_suffix='',
        verbose=True,
    ):
        fig, ax, genomeplotter, plotdata = self.draw_preprocess(
            y_colname=y_colname,
            ax=ax,
            genomeplotter=genomeplotter,
            frac=frac,
            subplots_kwargs=subplots_kwargs,
            log_suffix=log_suffix,
            nproc=nproc,
            verbose=verbose,
            plotdata=plotdata,
        )

        # drawing parameters
        plot_kwargs = (
            {
                'color': 'black', 
                'markersize': 0.3, 
                'alpha': libgenomeplot.calc_dot_alpha_baf(plotdata.nrow),
            }
            | plot_kwargs
        )

        # do drawing
        line2d, _ = genomeplotter.draw_dots(
            ax,
            plotdata=plotdata,
            y_colname=y_colname,
            plot_kwargs=plot_kwargs,
            draw_common=False,
        )

        if setup_axes:
            default_ymin, default_ymax = plotmisc.get_boxplot_range(plotdata[y_colname])
            if ymin is False:
                ymin = default_ymin
            if ymax is False:
                ymax = default_ymax
            draw_common_artists = draw_axessetup(
                ax=ax,
                genomeplotter=genomeplotter,
                ylabel=ylabel,
                ylabel_prefix=ylabel_prefix,
                ylabel_kwargs=ylabel_kwargs,
                ymax=ymax,
                ymin=ymin,
                yticks=yticks,
                draw_common_kwargs=draw_common_kwargs,
                rotate_chromlabel=rotate_chromlabel,
                fig=fig,
                title=title,
                suptitle_kwargs=suptitle_kwargs,
            )
        else:
            draw_common_artists = None

        #return fig, ax, genomeplotter, plotdata
        #plothandle = GenomePlotHandle(fig, genomeplotter)
        return GenomeDrawingAxesResult(
            ax=ax, 
            genomeplotter=genomeplotter, 
            artist=line2d,
            draw_common_artists=draw_common_artists,
            y_colname=y_colname,
            plotdata=plotdata,
        )

    def draw_hlines(
        self, 
        y_colname,
        ax=None, 
        genomeplotter=None,
        offset=None,

        # fig generation params
        subplots_kwargs=dict(),

        # drawing kwargs
        plot_kwargs=dict(),

        # axes setting
        setup_axes=True,
        ylabel=None,
        ylabel_prefix='',
        ylabel_kwargs=dict(),
        ymax=False,
        ymin=False,
        yticks=None,
        draw_common_kwargs=dict(),
        rotate_chromlabel=None,
        title=None,
        suptitle_kwargs=dict(),

        # pre-generated plotdata
        plotdata=None,

        # multicore plotdata generation
        nproc=1,
        log_suffix='',
        verbose=True,
    ):
        fig, ax, genomeplotter, plotdata = self.draw_preprocess(
            y_colname=y_colname,
            ax=ax,
            genomeplotter=genomeplotter,
            frac=None,
            subplots_kwargs=subplots_kwargs,
            log_suffix=log_suffix,
            nproc=nproc,
            verbose=verbose,
            plotdata=plotdata,
        )

        # drawing parameters
        plot_kwargs = (
            {'color': 'tab:blue', 'linewidth': 2, 'alpha': 1}
            | plot_kwargs
        )

        # do drawing
        linecol, _ = genomeplotter.draw_hlines(
            ax,
            offset=offset,
            plotdata=plotdata,
            y_colname=y_colname,
            plot_kwargs=plot_kwargs,
            draw_common=False,
        )

        if setup_axes:
            default_ymin, default_ymax = plotmisc.get_boxplot_range(plotdata[y_colname])
            if ymin is False:
                ymin = default_ymin
            if ymax is False:
                ymax = default_ymax
            draw_common_artists = draw_axessetup(
                ax=ax,
                genomeplotter=genomeplotter,
                ylabel=ylabel,
                ylabel_prefix=ylabel_prefix,
                ylabel_kwargs=ylabel_kwargs,
                ymax=ymax,
                ymin=ymin,
                yticks=yticks,
                draw_common_kwargs=draw_common_kwargs,
                rotate_chromlabel=rotate_chromlabel,
                fig=fig,
                title=title,
                suptitle_kwargs=suptitle_kwargs,
            )
        else:
            draw_common_artists = None

        #return fig, ax, genomeplotter, plotdata
        #plothandle = GenomePlotHandle(fig, genomeplotter)
        return GenomeDrawingAxesResult(
            ax=ax, 
            genomeplotter=genomeplotter, 
            artist=linecol,
            draw_common_artists=draw_common_artists,
            y_colname=y_colname,
            plotdata=plotdata,
        )

    def draw_boxes(
        self, 

        ymin_colname=None,
        ymax_colname=None,
        color_colname=None,

        ax=None, 
        genomeplotter=None,

        # fig generation params
        title=None,
        suptitle_kwargs=dict(),
        subplots_kwargs=dict(),

        # drawing kwargs
        plot_kwargs=dict(),

        # axes setting
        setup_axes=True,
        ylabel=None,
        ylabel_prefix='',
        ylabel_kwargs=dict(),
        ymax=False,
        ymin=False,
        yticks=None,
        draw_common_kwargs=dict(),
        rotate_chromlabel=None,

        # pre-generated plotdata
        plotdata=None,

        # multicore plotdata generation
        nproc=1,
        log_suffix='',
        verbose=True,
    ):
        y_colnames = [
            x for x in [ymin_colname, ymax_colname, color_colname] 
            if x is not None
        ]
        if len(y_colnames) == 0:
            y_colnames = None

        fig, ax, genomeplotter, plotdata = self.draw_preprocess(
            y_colname=y_colnames,
            ax=ax,
            genomeplotter=genomeplotter,
            frac=None,
            subplots_kwargs=subplots_kwargs,
            log_suffix=log_suffix,
            nproc=nproc,
            verbose=verbose,
            plotdata=plotdata,
        )

        # drawing parameters
        plot_kwargs = (
            {'alpha': 0.1, 'zorder': 0}
            | plot_kwargs
        )
        if ymin_colname is None:
            ymins = None
        else:
            ymins = plotdata[ymin_colname]
        if ymax_colname is None:
            ymaxs = None
        else:
            ymaxs = plotdata[ymax_colname]
        if color_colname is None:
            colors = None
        else:
            colors = plotdata[color_colname]

        # do drawing
        rects = genomeplotter.draw_boxes(
            ax,
            plotdata=plotdata,

            ymins=ymins,
            ymaxs=ymaxs,
            colors=colors,

            rect_kwargs=plot_kwargs,
            draw_common=False,
        )

        if setup_axes:
            draw_common_artists = draw_axessetup(
                ax=ax,
                genomeplotter=genomeplotter,
                ylabel=ylabel,
                ylabel_prefix=ylabel_prefix,
                ylabel_kwargs=ylabel_kwargs,
                ymax=ymax,
                ymin=ymin,
                yticks=yticks,
                draw_common_kwargs=draw_common_kwargs,
                rotate_chromlabel=rotate_chromlabel,
                fig=fig,
                title=title,
                suptitle_kwargs=suptitle_kwargs,
            )
        else:
            draw_common_artists = None

        #return fig, ax, genomeplotter, plotdata
        #plothandle = GenomePlotHandle(fig, genomeplotter)
        return GenomeDrawingAxesResult(
            ax=ax, 
            genomeplotter=genomeplotter, 
            artist=rects,
            draw_common_artists=draw_common_artists,
            plotdata=plotdata,
        )


class GenomeDrawingResultBase:
    def connect(self):
        self.intv_selector.connect()

    def disconnect(self):
        self.intv_selector.disconnect()

    def reset(self):
        self.intv_selector.reset()

    def get_region_gdfs(self):
        intv_plotcoords = self.intv_selector.get_intervals()
        left_chroms, left_pos0s = self.genomeplotter.plot_to_genomic(intv_plotcoords[:, 0])
        right_chroms, right_pos0s = self.genomeplotter.plot_to_genomic(intv_plotcoords[:, 1])
        right_end0s = right_pos0s + 1

        results = list()
        for chrom1, position1, chrom2, position2 in zip(
            left_chroms, left_pos0s, right_chroms, right_end0s,
        ):
            if chrom1 != chrom2:
                raise Exception(f'On-press chrom and on-release chrom differ.')
                
            results.append(
                GenomeDataFrameBase.from_margins(
                    refver=self.genomeplotter.refver, 
                    chrom_left=chrom1, 
                    start0_left=position1, 
                    chrom_right=chrom2, 
                    end0_right=position2,
                )
            )

        return results

    def get_region_lengths(self):
        return [
            self.genomeplotter.region_gdf.intersect(x).lengths.sum()
            for x in self.get_region_gdfs()
        ]


def get_region_mean_values(plotdata, y_colname, region_gdfs):
    assert plotdata is not None
    assert y_colname is not None

    result = list()
    for gdf in region_gdfs:
        isec = plotdata.intersect(gdf)
        if isec.is_empty:
            result.append(np.nan)
        else:
            result.append(isec.get_lwavg(y_colname))
    return np.atleast_1d(result)


class GenomeDrawingAxesResult(GenomeDrawingResultBase):
    """variable name: gdraw_result"""
    def __init__(
        self, 
        ax, 
        genomeplotter, 
        artist,
        draw_common_artists,

        y_colname=None,
        plotdata=None,
        name=None,

    ):
        #self.fig = fig
        self.genomeplotter = genomeplotter
        self.ax = ax
        self.artist = artist
        self.draw_common_artists = draw_common_artists
        self.y_colname = y_colname
        self.plotdata = plotdata
        self.name = name

        self.intv_selector = IntervalSelector.from_ax(self.ax)

    @property
    def fig(self):
        return self.ax.figure

    def set_name(self, x):
        self.name = x

    def get_mean_values(self, region_gdfs=None):
        if region_gdfs is None:
            region_gdfs = self.get_region_gdfs()
        return get_region_mean_values(self.plotdata, self.y_colname, region_gdfs)


class GenomeDrawingFigureResult(GenomeDrawingResultBase):
    ########
    # init #
    ########

    def __init__(
        self, 
        axresult_list,
        omit_ax_labels=list(),
    ):
        self.axresults = {x.name: x for x in axresult_list}

        # sanitycheck
        self.init_sanitycheck()

        # derived attributes
        self.intv_selector = IntervalSelector.from_fig(self.fig, omit_ax_labels=omit_ax_labels)

    def init_sanitycheck(self):
        # 1. compare figures
        assert len(set(x.ax.figure for x in self.axresults.values())) == 1, (
            f'Figure object of input {GenomeDrawingAxesResult.__name__} objects differ' 
        )

        # 2. compare GenomePlotter
        gplotter_list = tuple(x.genomeplotter for x in self.axresults.values())
        assert all((x == y) for (x, y) in tools.pairwise(gplotter_list))

    #############
    # utilities #
    #############

    @property
    def first_axresult(self):
        return next(iter(self.axresults.values()))

    @property
    def fig(self):
        return self.first_axresult.ax.figure

    @property
    def axd(self):
        return {ax.get_label(): ax for ax in self.fig.get_axes()}

    @property
    def genomeplotter(self):
        return self.first_axresult.genomeplotter

    def get_mean_values(self, names=None):
        if names is None:
            names = tuple(self.axresults.keys())

        region_gdfs = self.get_region_gdfs()

        result = dict()
        for name in names:
            result[name] = self.axresults[name].get_mean_values(
                region_gdfs=region_gdfs,
            )

        return result



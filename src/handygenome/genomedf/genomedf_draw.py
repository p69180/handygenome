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
#import handygenome.tools as tools
#import handygenome.logutils as logutils

from handygenome.genomedf.genomedf_base import GenomeDataFrameBase

import handygenome.plot.misc as plotmisc
from handygenome.plot.gui import IntervalSelector
import handygenome.plot.genomeplot as libgenomeplot
from handygenome.plot.genomeplot import GenomePlotter


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

    def draw_axessetup(
        self, 
        ax,
        genomeplotter,

        # y axes setting
        ylabel,
        ylabel_prefix,
        ylabel_kwargs,
        ymax,
        ymin,
        yticks,

        # draw common
        draw_common_kwargs,
        rotate_chromlabel,

        # figure suptitle
        fig,
        title,
        suptitle_kwargs,
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

        genomeplotter.draw_common(ax, **draw_common_kwargs)

        # figure suptitle
        if title is not None:
            plotmisc.draw_suptitle(fig, title, **suptitle_kwargs)

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
        genomeplotter.draw_dots(
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
            self.draw_axessetup(
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

        #return fig, ax, genomeplotter, plotdata
        #plothandle = GenomePlotHandle(fig, genomeplotter)
        return GenomeDrawingResult(
            fig=fig, 
            ax=ax, 
            genomeplotter=genomeplotter, 
            #plothandle=plothandle, 
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
        genomeplotter.draw_hlines(
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
            self.draw_axessetup(
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

        #return fig, ax, genomeplotter, plotdata
        #plothandle = GenomePlotHandle(fig, genomeplotter)
        return GenomeDrawingResult(
            fig=fig, 
            ax=ax, 
            genomeplotter=genomeplotter, 
            #plothandle=plothandle, 
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
        genomeplotter.draw_boxes(
            ax,
            plotdata=plotdata,

            ymins=ymins,
            ymaxs=ymaxs,
            colors=colors,

            rect_kwargs=plot_kwargs,
            draw_common=False,
        )

        if setup_axes:
            self.draw_axessetup(
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

        #return fig, ax, genomeplotter, plotdata
        #plothandle = GenomePlotHandle(fig, genomeplotter)
        return GenomeDrawingResult(
            fig=fig, 
            ax=ax, 
            genomeplotter=genomeplotter, 
            #plothandle=plothandle, 
            plotdata=plotdata,
        )


#class GenomePlotHandle:
#    def __init__(self, fig, genomeplotter, omit_ax_labels=list()):
#        self.fig = fig
#        self.gplotter = genomeplotter
#        self.intv_selector = IntervalSelector(fig, omit_ax_labels=omit_ax_labels)
#
#    def connect(self):
#        self.intv_selector.connect()
#
#    def disconnect(self):
#        self.intv_selector.disconnect()
#
#    def get_region_gdfs(self):
#        intv_plotcoords = self.intv_selector.get_intervals()
#        left_chroms, left_pos0s = self.gplotter.plot_to_genomic(intv_plotcoords[:, 0])
#        right_chroms, right_pos0s = self.gplotter.plot_to_genomic(intv_plotcoords[:, 1])
#        right_end0s = right_pos0s + 1
#
#        results = list()
#        for chrom1, position1, chrom2, position2 in zip(
#            left_chroms, left_pos0s, right_chroms, right_end0s,
#        ):
#            if chrom1 != chrom2:
#                raise Exception(f'On-press chrom and on-release chrom differ.')
#                
#            results.append(
#                GenomeDataFrameBase.from_margins(
#                    self.gplotter.refver, 
#                    chrom1, 
#                    position1, 
#                    chrom2, 
#                    position2,
#                )
#            )
#
#        return results
#

class GenomeDrawingResult:
    """variable name: gdraw_result"""
    def __init__(
        self, fig, genomeplotter, 

        ax=None, 
        #plothandle=None, 
        y_colname=None,
        plotdata=None,
        name=None,

        omit_ax_labels=list(),
    ):
        self.fig = fig
        self.genomeplotter = genomeplotter
        self.ax = ax
        #self.plothandle = plothandle
        self.y_colname = y_colname
        self.plotdata = plotdata
        self.name = name

        self.intv_selector = IntervalSelector(self.fig, omit_ax_labels=omit_ax_labels)

    def set_name(self, x):
        self.name = x

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

    def get_mean_values(self):
        assert self.plotdata is not None
        assert self.y_colname is not None

        result = list()
        for gdf in self.get_region_gdfs():
            isec = self.plotdata.intersect(gdf)
            if isec.is_empty:
                result.append(np.nan)
            else:
                result.append(isec.get_lwavg(self.y_colname))
        return result

class GenomeDrawingResultGroup:
    def __init__(self, *args):
        self.gdraw_results = args




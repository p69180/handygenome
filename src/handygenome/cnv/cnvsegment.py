import multiprocessing
import functools
import re

import pandas as pd
import numpy as np
import pyranges as pr
import matplotlib as mpl
import matplotlib.pyplot as plt
import pysam

import handygenome.peakutils as peakutils
import handygenome.logutils as logutils
from handygenome.genomedf.genomedf import GenomeDataFrame
import handygenome.cnv.cncall as cncall
from handygenome.cnv.depth import DepthSegmentDataFrame
from handygenome.cnv.baf import BAFSegmentDataFrame
import handygenome.cnv.baf as libbaf
import handygenome.plot.misc as plotmisc


class CNVSegmentDataFrame(DepthSegmentDataFrame, BAFSegmentDataFrame):
    clonal_CN_colname = 'CN'
    clonal_B_colname_prefix = 'B'
    subclonal_CN_colname = 'subCN'
    subclonal_B_colname_prefix = 'subB'
    ccf_colname = 'ccf'
    predicted_depth_colname = 'predicted_depth'
    predicted_baf_colname_prefix = 'predicted_baf'

    clonal_B_colname_pat = re.compile(
        f'{clonal_B_colname_prefix}_(?P<baf_index>{libbaf.BAFINDEX_PAT_STRING})'
    )
    subclonal_B_colname_pat = re.compile(
        f'{subclonal_B_colname_prefix}_(?P<baf_index>{libbaf.BAFINDEX_PAT_STRING})'
    )


    @classmethod
    def from_segments(cls, seg_gdfs, nproc=1):
        assert len(seg_gdfs) > 0

        refvers = set(x.refver for x in seg_gdfs)
        if len(refvers) != 1:
            raise Exception(f'refver differs between segment gdfs: {refvers}')
        refver = refvers.pop()

        merged_seg_gdf = functools.reduce(
            lambda seg1, seg2: (
                seg1.drop_annots('N', inplace=False).isec_union(
                    seg2.drop_annots('N', inplace=False),
                    nproc=nproc,
                    drop_annots=False,
                    sort=True,
                )
            ),
            seg_gdfs,
        )

        return cls.from_frame(merged_seg_gdf.df, refver=refver)

    def check_has_CN(self):
        return self.__class__.clonal_CN_colname in self.columns

    def check_has_B_columns(self):
        return any(
            self.__class__.clonal_B_colname_pat.fullmatch(x)
            for x in self.columns
        )

    ###################
    # colname getters #
    ###################

    def get_clonal_CN_colname(self, germline=False):
        result = self.__class__.clonal_CN_colname
        if germline:
            result = 'g' + result
        return result

    def get_clonal_B_colname(self, baf_index, germline=False):
        assert libbaf.check_valid_bafindex(baf_index)
        result = f'{self.__class__.clonal_B_colname_prefix}_{baf_index}'
        if germline:
            result = 'g' + result
        return result

    def get_clonal_B_colname_dict(self):
        result = [x for x in self.columns if self.__class__.clonal_B_colname_pat.fullmatch(x)]
        result.sort(
            key=(
                lambda x: int(
                    self.__class__.clonal_B_colname_pat
                    .fullmatch(x)
                    .group('num')
                )
            )
        )
        result = dict(
            (self.__class__.clonal_B_colname_pat.fullmatch(x).group('baf_index'), x)
            for x in result
        )
        if len(result) == 0:
            raise Exception(f'Does not have any clonal B columns')
        return result

    def get_subclonal_CN_colname(self, germline=False):
        result = self.__class__.subclonal_CN_colname
        if germline:
            result = 'g' + result
        return result

    def get_subclonal_B_colname(self, baf_index, germline=False):
        assert libbaf.check_valid_bafindex(baf_index)
        result = f'{self.__class__.subclonal_B_colname_prefix}_{baf_index}'
        if germline:
            result = 'g' + result
        return result

    def get_all_subclonal_B_colnames(self):
        result = [x for x in self.columns if self.__class__.subclonal_B_colname_pat.fullmatch(x)]
        if len(result) == 0:
            raise Exception(f'Does not have any subclonal B columns')
        return result

    def get_predicted_depth_colname(self):
        return self.__class__.predicted_depth_colname

    def get_predicted_baf_colname(self, baf_index):
        return f'{self.__class__.predicted_baf_colname_prefix}_{baf_index}'

    ##########
    # assign #
    ##########

    def assign_clonal_CN(self, data, germline=False):
        self[self.get_clonal_CN_colname(germline=germline)] = data

    def assign_clonal_B(self, data, baf_index, germline=False):
        self[self.get_clonal_B_colname(baf_index, germline=germline)] = data

    def assign_subclonal_CN(self, data, germline=False):
        self[self.get_subclonal_CN_colname(germline=germline)] = data

    def assign_subclonal_B(self, data, baf_index, germline=False):
        self[self.get_subclonal_B_colname(baf_index, germline=germline)] = data

    def assign_ccf(self, data):
        self[self.__class__.ccf_colname] = data

    def assign_predicted_depth(self, data):
        self[self.get_predicted_depth_colname()] = data

    def assign_predicted_baf(self, data, baf_index):
        self[self.get_predicted_baf_colname(baf_index)] = data

    #########
    # fetch #
    #########

    def get_clonal_CN(self, germline=False):
        return self[self.get_clonal_CN_colname(germline=germline)]

    def get_clonal_B(self, baf_index, germline=False):
        return self[self.get_clonal_B_colname(baf_index, germline=germline)]

    def get_subclonal_CN(self, germline=False):
        return self[self.get_subclonal_CN_colname(germline=germline)]

    def get_subclonal_B(self, baf_index, germline=False):
        return self[self.get_subclonal_B_colname(baf_index, germline=germline)]

    def get_ccf(self):
        return self[self.__class__.ccf_colname]

    def get_predicted_depth(self):
        return self[self.get_predicted_depth_colname()]

    def get_predicted_baf(self, baf_index):
        return self[self.get_predicted_baf_colname(baf_index)]

    ########
    # draw #
    ########

    def draw(
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
    ):
        # y axis parameters
        if ylabel is False:
            ylabel = f'clonal copy number'
        if ymax is False:
            ymax = np.nanquantile(self.get_clonal_CN(), 0.99)
            #q25, q75 = np.nanquantile(self.get_clonal_CN(), [0.25, 0.75])
            #iqr = q75 - q25
            #ymax = q75 + 1.5 * iqr
        if ymin is False:
            ymin = -0.4
        if yticks is False:
            yticks = plotmisc.get_integer_yticks(
                ymin=0, 
                ymax=ymax, 
                num_ticks=15,
            )

        # prepare single plotdata
        if genomeplotter is None:
            genomeplotter = self.get_default_genomeplotter()
        if plotdata is None:
            plotdata = genomeplotter.make_plotdata(
                self, 
                log_suffix=' (clonal copy number)',
                nproc=nproc,
            )

        # plot kwargs
        plot_kwargs_base = (
            {'linewidth': 2, 'alpha': 1}
            | plot_kwargs
        )

        # draw CN
        CN_line_color = 'black'
        plot_kwargs = plot_kwargs_base | {'color': CN_line_color}
        offset = 0.1
        fig, ax, genomeplotter = self.draw_hlines(
            y_colname=self.get_clonal_CN_colname(),
            ax=ax,
            genomeplotter=genomeplotter,
            offset=offset,

            plotdata=plotdata,
            plot_kwargs=plot_kwargs,

            setup_axes=False,
            subplots_kwargs=subplots_kwargs,
        )

        # draw B
        B_colnames = self.get_clonal_B_colname_dict()
        colors = mpl.cm.cool(np.linspace(0, 1, len(B_colnames), endpoint=True))
        B_line_colors = dict(zip(iter(B_colnames.keys()), colors))
        for idx, (baf_idx, colname) in enumerate(B_colnames.items()):
            offset = idx * -0.1
            plot_kwargs = plot_kwargs_base | {'color': B_line_colors[baf_idx]}
            fig, ax, genomeplotter = self.draw_hlines(
                y_colname=colname,
                ax=ax,
                genomeplotter=genomeplotter,
                offset=offset,

                plotdata=plotdata,
                plot_kwargs=plot_kwargs,

                setup_axes=False,
            )

        # axes setup
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

        # make legend - intentionally placed after axes setup to use "bbox_to_anchor"
        handles = plotmisc.LegendHandles()
        handles.add_line(
            marker=None, 
            linewidth=4, 
            color=CN_line_color, 
            label='total copy number'
        )
        for baf_idx, color in B_line_colors.items():
            handles.add_line(
                marker=None, 
                linewidth=4, 
                color=color, 
                label=f'{baf_idx} B allele copy number',
            )
        ax.legend(handles=handles, loc='upper right', bbox_to_anchor=(1, 1.3))





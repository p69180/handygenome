import multiprocessing
import functools
import re

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pysam

import handygenome.tools as tools
import handygenome.peakutils as peakutils
import handygenome.logutils as logutils
from handygenome.genomedf.genomedf import GenomeDataFrame
import handygenome.cnv.cnvcall as cnvcall
from handygenome.cnv.depth import DepthSegmentDataFrame
from handygenome.cnv.baf import BAFSegmentDataFrame
import handygenome.cnv.baf as libbaf
import handygenome.plot.misc as plotmisc
import handygenome.cnv.bafsimul as bafsimul
import handygenome.genomedf.genomedf_draw as genomedf_draw
from handygenome.genomedf.genomedf_draw import GenomeDrawingFigureResult


BAF_COLNAME_PAT = re.compile(libbaf.BAFINDEX_PAT.pattern + '.*')


class CNVSegmentDataFrame(DepthSegmentDataFrame, BAFSegmentDataFrame):
    clonal_CN_colname = 'CN'
    clonal_B_colname_prefix = 'B'
    subclonal_CN_colname = 'subCN'
    subclonal_B_colname_prefix = 'subB'
    ccf_colname = 'ccf'

    predicted_depth_colname = 'predicted_depth'
    predicted_baf_colname_prefix = 'predicted_baf'
    predicted_depth_drawresult_name = 'predicted_depth'
    predicted_baf_drawresult_name_prefix = 'predicted_baf_'
    predicted_default_kwargs = dict(
        color='tab:green',
        alpha=genomedf_draw.DEFAULT_SEGMENT_ALPHA,
        linewidth=genomedf_draw.DEFAULT_SEGMENT_LINEWIDTH,
    )

    clonal_B_colname_pat = re.compile(
        f'{clonal_B_colname_prefix}_(?P<baf_index>{libbaf.BAFINDEX_PAT_STRING})'
    )
    subclonal_B_colname_pat = re.compile(
        f'{subclonal_B_colname_prefix}_(?P<baf_index>{libbaf.BAFINDEX_PAT_STRING})'
    )

    #corrected_baf_color = 'tab:red'
    corrected_baf_default_kwargs = dict(
        color='tab:red',
        alpha=genomedf_draw.DEFAULT_SEGMENT_ALPHA,
        linewidth=genomedf_draw.DEFAULT_SEGMENT_LINEWIDTH,
    )
    baf_drawresult_name = 'corrected_baf_segment'

    CN_drawresult_name = 'CN'
    B_drawresult_name = 'B'

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
        merged_seg_gdf = cls.from_frame(merged_seg_gdf.df, refver=refver)

        # add corrected baf
        for baf_index in merged_seg_gdf.get_baf_indexes():
            merged_seg_gdf.add_corrected_baf_interp(baf_index=baf_index)

        return merged_seg_gdf

    def check_has_CN(self):
        return self.__class__.clonal_CN_colname in self.columns

    def check_has_B_columns(self):
        return any(
            self.__class__.clonal_B_colname_pat.fullmatch(x)
            for x in self.columns
        )

    def get_baf_indexes(self):
        """baf0 comes first, the baf1, then baf2, ..."""
        matches = [BAF_COLNAME_PAT.fullmatch(x) for x in self.columns]
        matches = [x for x in matches if x is not None]
        matches.sort(key=(lambda x: int(x.group('num'))))
        result = tools.unique_keeporder(
            (x.group('bafindex_prefix') + x.group('num'))
            for x in matches
        )
        return result

    #####################
    # add corrected baf #
    #####################

    def add_corrected_baf_interp(self, baf_index):
        baf_mean_colname = self.get_baf_mean_colname(baf_index=baf_index)
        depth_colname = self.depth_mean_colname
        input_data = self[[baf_mean_colname, depth_colname]]
        corrected_bafs = bafsimul.predict_true_baf(input_data)
        self[self.get_corrected_baf_colname(baf_index=baf_index)] = corrected_bafs

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

    ###################################
    # predicted value column creation #
    ###################################

    def add_predicted_depth_germline(self, is_female, ploidy):
        onecopy_depth = self.find_onecopy_depth(is_female, ploidy)
        K = cnvcall.get_K_from_onecopy_depth_germline(onecopy_depth)

        predicted_depth = cnvcall.get_predicted_depth(
            CNt=self.get_clonal_CN(germline=False),
            CNg=None,
            cellularity=1,
            K=K,
        )
        self[self.get_predicted_depth_colname()] = predicted_depth

    def add_predicted_baf_germline(self, baf_index):
        predicted_baf = cnvcall.get_predicted_baf(
            CNt=self.get_clonal_CN(germline=False),
            Bt=self.get_clonal_B(baf_index=baf_index, germline=False),
            cellularity=1,
        )
        self[self.get_predicted_baf_colname(baf_index)] = predicted_baf

    def add_predicted_depth_nongermline(self, cellularity, K):
        predicted_depth = cnvcall.get_predicted_depth(
            CNt=self.get_clonal_CN(germline=False),
            CNg=self.get_clonal_CN(germline=True),
            cellularity=cellularity,
            K=K,
        )
        self[self.get_predicted_depth_colname()] = predicted_depth

    def add_predicted_baf_nongermline(self, baf_index, cellularity):
        predicted_baf = cnvcall.get_predicted_baf(
            CNt=self.get_clonal_CN(germline=False),
            Bt=self.get_clonal_B(baf_index=baf_index, germline=False),
            CNg=self.get_clonal_CN(germline=True),
            Bg=self.get_clonal_B(baf_index=baf_index, germline=True),
            cellularity=cellularity,
        )
        self[self.get_predicted_baf_colname(baf_index)] = predicted_baf

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

    #########################
    # solution-related ones #
    #########################

    def add_clonal_solution_germline(self, is_female, ploidy):
        assert isinstance(ploidy, int)
        clonal_CNt, clonal_Bt = cnvcall.clonal_solution_from_germline(
            self, is_female, ploidy,
        )
        self.assign_clonal_CN(clonal_CNt)
        for idx, baf_index in enumerate(self.get_baf_indexes()):
            self.assign_clonal_B(clonal_Bt[:, idx], baf_index)

        # predicted depth and baf
        self.add_predicted_depth_germline(is_female, ploidy)
        for baf_index in self.get_baf_indexes():
            self.add_predicted_baf_germline(baf_index)

    def add_clonal_solution_targetsample(self, cellularity, K):
        """depth and baf of this sample must be corrected"""
        #assert self.check_has_corrected_depth()

        # CNt
        CNt = cnvcall.find_clonal_CNt(
            #corrected_depth=self.corrected_depth_mean,
            corrected_depth=self.norm_depth_mean,
            cellularity=cellularity,
            K=K,
            CNg=self.get_clonal_CN(germline=True),
        )
        self.assign_clonal_CN(data=CNt, germline=False)

        # Bt
        for baf_index in self.get_baf_indexes():
            Bt = cnvcall.find_clonal_Bt(
                baf=self.get_corrected_baf(baf_index),
                CNt=self.get_clonal_CN(germline=False),
                cellularity=cellularity,
                CNg=self.get_clonal_CN(germline=True),
                Bg=self.get_clonal_B(baf_index, germline=True),
            )
            self.assign_clonal_B(data=Bt, baf_index=baf_index, germline=False)

        # predicted depth and baf
        self.add_predicted_depth_nongermline(cellularity, K)

        for baf_index in self.get_baf_indexes():
            self.add_predicted_baf_nongermline(baf_index, cellularity)

    def drop_solution_columns(self):
        cols_to_drop = set(self.columns).intersection(
            [
                self.get_clonal_CN_colname(germline=False),
                self.get_predicted_depth_colname(),
            ]
            + list(
                self.get_clonal_B_colname(baf_index, germline=False)
                for baf_index in self.get_baf_indexes()
            )
            + list(
                self.get_predicted_baf_colname(baf_index)
                for baf_index in self.get_baf_indexes()
            )
        )
        self.drop_annots(cols_to_drop, inplace=True)

    def get_average_ploidy(self):
        return cnvcall.get_average_ploidy(self.get_clonal_CN(), self.lengths)

    @staticmethod
    def get_fitness_helper(values, lengths):
        notnan_selector = np.logical_not(np.isnan(values))
        values = values[notnan_selector]
        lengths = lengths[notnan_selector]
        return np.sum(np.abs(values) * lengths)

    def get_CNt_fitness(self):
        return self.get_fitness_helper(
            (self.get_predicted_depth() - self.corrected_depth_mean),
            self.lengths,
        )

    def get_Bt_fitness_one_baf_index(self, baf_index):
        return self.get_fitness_helper(
            (self.get_predicted_baf(baf_index) - self.get_corrected_baf(baf_index)),
            self.lengths,
        )

    def get_Bt_fitness(self):
        return sum(self.get_Bt_fitness_one_baf_index(baf_index) for baf_index in self.get_baf_indexes())

    def get_solution_fitness(self):
        return self.get_CNt_fitness() + self.get_Bt_fitness()

    ########
    # draw #
    ########

    def draw_CN(
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
        verbose=True,
    ):
        # y axis parameters
        if ylabel is False:
            ylabel = f'clonal copy number'
        if ymax is False:
            #ymax = np.nanmax(self.get_clonal_CN()) + 0.4
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
                verbose=verbose,
            )

        # plot kwargs
        plot_kwargs_base = (
            {'linewidth': 2, 'alpha': 1}
            | plot_kwargs
        )

        # prepare gdraw list
        gdraw_result_list = list()

        # draw CN
        CN_line_color = 'black'
        plot_kwargs = plot_kwargs_base | {'color': CN_line_color}
        offset = 0.1
        #fig, ax, genomeplotter, plotdata = self.draw_hlines(
        gdraw_result_CN = self.draw_hlines(
            y_colname=self.get_clonal_CN_colname(),
            ax=ax,
            genomeplotter=genomeplotter,
            offset=offset,

            plotdata=plotdata,
            plot_kwargs=plot_kwargs,

            setup_axes=False,
            subplots_kwargs=subplots_kwargs,
        )

        ax = gdraw_result_CN.ax
        fig = gdraw_result_CN.fig
        genomeplotter = gdraw_result_CN.genomeplotter

        gdraw_result_CN.set_name(self.__class__.CN_drawresult_name)
        gdraw_result_list.append(gdraw_result_CN)

        # draw B
        B_colnames = self.get_clonal_B_colname_dict()
        colors = mpl.cm.cool(np.linspace(0, 1, len(B_colnames), endpoint=True))
        B_line_colors = dict(zip(iter(B_colnames.keys()), colors))
        #gdraw_result_B_
        for idx, (baf_idx, colname) in enumerate(B_colnames.items()):
            offset = idx * -0.1
            plot_kwargs = plot_kwargs_base | {'color': B_line_colors[baf_idx]}
            #fig, ax, genomeplotter, plotdata = self.draw_hlines(
            gdraw_result_B = self.draw_hlines(
                y_colname=colname,
                ax=ax,
                genomeplotter=genomeplotter,
                offset=offset,

                plotdata=plotdata,
                plot_kwargs=plot_kwargs,

                setup_axes=False,
            )
            gdraw_result_B.set_name(self.__class__.B_drawresult_name + '_' + baf_idx)
            gdraw_result_list.append(gdraw_result_B)

        # axes setup
        if setup_axes:
            genomedf_draw.draw_axessetup(
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

        gdraw_figresult = GenomeDrawingFigureResult(gdraw_result_list)
        return gdraw_figresult

    def draw_corrected_baf(
        self,
        baf_index,
        ax=None,
        genomeplotter=None,
        plotdata=None,

        # drawing kwargs
        plot_kwargs=dict(),

        # multicore plotdata generation
        nproc=1,
        verbose=True,
    ):
        y_colname_corr = self.get_corrected_baf_colname(baf_index=baf_index)
        #fig, ax, genomeplotter, plotdata = self.draw_hlines(
        gdraw_result = self.draw_hlines(
            y_colname=y_colname_corr,
            ax=ax,
            genomeplotter=genomeplotter,
            plotdata=plotdata,

            plot_kwargs=(
                self.__class__.corrected_baf_default_kwargs
                | plot_kwargs
            ),

            setup_axes=False,

            nproc=nproc,
            log_suffix=' (BAF segment - corrected baf)',
            verbose=verbose,
        )

        # make legend - intentionally placed after axes setup to use "bbox_to_anchor"
        handles = plotmisc.LegendHandles()
        handles.add_line(
            marker=None, 
            linewidth=4, 
            color=self.__class__.corrected_baf_default_kwargs['color'], 
            label='corrected BAF'
        )
        handles.add_line(
            marker=None, 
            linewidth=4, 
            color=BAFSegmentDataFrame.default_plot_kwargs['color'], 
            label='BAF data mean'
        )
        gdraw_result.ax.legend(handles=handles, loc='upper right', bbox_to_anchor=(1, 1.3))

        gdraw_result.set_name(self.__class__.baf_drawresult_name)
        return gdraw_result

    def draw_predicted_baf(
        self,
        baf_index,
        ax,

        genomeplotter=None,
        plotdata=None,
        plot_kwargs=dict(),

        nproc=1,
        verbose=True,
    ):
        plot_kwargs = (self.__class__.predicted_default_kwargs | plot_kwargs)
        gdraw_result = self.draw_hlines(
            ax=ax,
            y_colname=self.get_predicted_baf_colname(baf_index=baf_index),
            genomeplotter=genomeplotter,
            plotdata=plotdata,
            setup_axes=False,
            plot_kwargs=plot_kwargs,
            nproc=nproc,
            verbose=verbose,
        )
        gdraw_result.set_name(
            self.__class__.predicted_baf_drawresult_name_prefix + str(baf_index)
        )

        return gdraw_result

    def draw_predicted_depth(
        self,
        ax,

        genomeplotter=None,
        plotdata=None,
        plot_kwargs=dict(),

        nproc=1,
        verbose=True,
    ):
        plot_kwargs = (self.__class__.predicted_default_kwargs | plot_kwargs)
        gdraw_result = self.draw_hlines(
            ax=ax,
            y_colname=self.get_predicted_depth_colname(),
            genomeplotter=genomeplotter,
            plotdata=plotdata,
            setup_axes=False,
            plot_kwargs=plot_kwargs,
            nproc=nproc,
            verbose=verbose,
        )
        gdraw_result.set_name(self.__class__.predicted_depth_drawresult_name)

        return gdraw_result




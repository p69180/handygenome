import handygenome.plot.misc as plotmisc
from handygenome.plot.genomeplot import GenomePlotter


#######################
# CNV complex drawing #
#######################

def draw_two_cnvsample_overlay(
    cnvsample1,
    corrector_id1,
    cnvsample2,
    corrector_id2,
    roi_gdf,
    roi_expand_factor=100,
    roi_alpha=0.5,
    depth_ymax=None,
    CN_ymax=False,
):
    assert roi_gdf.nrow == 1

    width = roi_gdf.lengths[0]
    chrom = roi_gdf.chroms[0]
    view_start0 = roi_gdf.start0s[0] - width * roi_expand_factor
    view_end0 = roi_gdf.end0s[0] + width * roi_expand_factor

    gene = roi_gdf['gene'][0]

    tps_depth_kwargs = dict(color='tab:red', markersize=3, alpha=0.05)
    tps_baf_kwargs = dict(color='tab:red', markersize=3, alpha=0.3)
    tps_CN_kwargs = dict(color='tab:red', markersize=3, alpha=0.7, marker='o')
    tps_B_kwargs = dict(color='tab:green', markersize=3, alpha=0.7, marker='o')
    roi_kwargs = dict(linewidth=2, color='yellow', alpha=1, marker=None, label=f'{gene} gene')

    gplotter = GenomePlotter(refver=cnvsample1.refver, chroms=chrom, start0s=view_start0, end0s=view_end0)
    draw_result = cnvsample1.draw(
        genomeplotter=gplotter,
        corrector_id=corrector_id1,
        draw_all=False,
        draw_CN=True,
        draw_baf=True,
        draw_depth=True,
        draw_MQ=False,
        nproc=10,
        verbose=False,
        depth_segment_kwargs=dict(ymax=depth_ymax),
        CN_kwargs=dict(ymax=CN_ymax),
        draw_common_kwargs=dict(n_xlabel=10),
        subplots_kwargs=dict(gridspec_kw=dict(hspace=0.2, width_ratios=[1, 0.2], wspace=0)),
        omit_title_mode=False,
    )

    for key, ax in draw_result.axd.items():
        if key.endswith('_legend'):
            continue
        dres_box = roi_gdf.draw_boxes(
            ax=ax,
            genomeplotter=gplotter,
            plot_kwargs=dict(alpha=roi_alpha),
            verbose=False,
            setup_axes=False,
        )

    dres_tpsdepth = cnvsample2.get_depth_rawdata(corrector_id2).draw_depth(
        ax=draw_result.axd['depth'],
        plot_kwargs=tps_depth_kwargs,
        genomeplotter=gplotter,
        setup_axes=False,
        verbose=False,
    )
    dres_tpsbaf = cnvsample2.get_baf_rawdata(corrector_id=corrector_id2).draw_baf(
        baf_idx='baf0',
        ax=draw_result.axd['baf0'],
        plot_kwargs=tps_baf_kwargs,
        genomeplotter=gplotter,
        setup_axes=False,
        verbose=False,
    )
    dres_tpsCN = cnvsample2.get_merged_cnv_segment(corrector_id2).draw_dots(
        y_colname='CN',
        ax=draw_result.axd['CN'],

        plot_kwargs=tps_CN_kwargs,

        genomeplotter=gplotter,
        setup_axes=False,
        verbose=False,
    )
    dres_tpsB = cnvsample2.get_merged_cnv_segment(corrector_id2).draw_dots(
        y_colname='B_baf0',
        ax=draw_result.axd['CN'],

        plot_kwargs=tps_B_kwargs,

        genomeplotter=gplotter,
        setup_axes=False,
        verbose=False,
    )


    leghandles_depth = plotmisc.LegendHandles()
    leghandles_depth.add_line(
        color='black',
        marker='o',
        linewidth=0,
        markersize=1,
        label='WGS depth',
    )
    leghandles_depth.add_line(
        color=tps_depth_kwargs['color'],
        marker='o',
        linewidth=0,
        markersize=tps_depth_kwargs['markersize'],
        label='TPS depth',
    )
    leghandles_depth.add_line(**roi_kwargs)


    leghandles_baf = plotmisc.LegendHandles()
    leghandles_baf.add_line(
        color='black',
        marker='o',
        linewidth=0,
        markersize=1,
        label='WGS B allele freq',
    )
    leghandles_baf.add_line(
        color=tps_baf_kwargs['color'],
        marker='o',
        linewidth=0,
        markersize=tps_baf_kwargs['markersize'],
        label='TPS B allele freq',
    )
    leghandles_baf.add_line(**roi_kwargs)

    leghandles_CN = plotmisc.LegendHandles()
    leghandles_CN.add_line(
        color='black',
        marker=None,
        linewidth=2,
        label='WGS total CN',
    )
    leghandles_CN.add_line(
        color='cyan',
        marker=None,
        linewidth=2,
        label='WGS B allele CN',
    )
    leghandles_CN.add_line(
        linewidth=0,
        **tps_CN_kwargs,
        label='TPS total CN',
    )
    leghandles_CN.add_line(
        linewidth=0,
        **tps_B_kwargs,
        label='TPS B allele CN',
    )
    leghandles_CN.add_line(**roi_kwargs)


    legend_common_kwargs = dict(loc='center left', bbox_to_anchor=(0, 0.5), fontsize='small')
    draw_result.axd['depth_legend'].legend(handles=leghandles_depth, **legend_common_kwargs)
    draw_result.axd['baf0_legend'].legend(handles=leghandles_baf, **legend_common_kwargs)
    draw_result.axd['CN_legend'].legend(handles=leghandles_CN, **legend_common_kwargs)

    return draw_result

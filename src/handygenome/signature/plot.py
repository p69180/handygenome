import re
import collections
import functools
import operator

import IPython
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import handygenome.signature.misc as sigmisc
import handygenome.plot.misc as plotmisc


IPYTHON = IPython.get_ipython()
SBS96_KEY_PAT = re.compile(r'^(.*)\[((.*)>(.*))\](.*)$')


def make_sbs6_dict_from_sbs96(catalogue_sbs96: pd.Series):
    result = dict()
    for key, val in catalogue_sbs96.to_dict().items():
        sbs6_key = SBS96_KEY_PAT.sub('\\2', key)  # e.g. A[C>T]G -> C>T
        result.setdefault(sbs6_key, 0)
        result[sbs6_key] += val

    return result


def draw_bars_sbs(ax, catalogue_dict, title, ylim, catalogue_keys_sbs96, 
                  add_cat_labels=True):
    ax.tick_params('x', bottom=False, labelbottom=False)
    ax.tick_params('y', labelsize=20)
    ax.set_title(title, size=30, pad=10)
    ax.set_xlim(-1, 96)
    ax.set_ylim(*ylim)

    for idx, key in enumerate(catalogue_keys_sbs96):
        mat = SBS96_KEY_PAT.match(key)
        pre = mat.group(1)
        ref = mat.group(3)
        alt = mat.group(4)
        post = mat.group(5)
        catalogue_key = f'{pre}[{ref}>{alt}]{post}'
        #label = f'{pre}[{ref}>{alt}]{post}'

        x = idx
        y = catalogue_dict[catalogue_key]
        color = sigmisc.COLORS_SBS6[f'{ref}>{alt}']
        ax.bar(x, y, color=color, align='center', width=0.4, zorder=1000)
        if add_cat_labels:
            ax.text(x, 0, catalogue_key, size=15, rotation='vertical', 
                    ha='center', va='top', family='monospace', weight='bold')


def draw_bars_exposure(ax, exposure_dict, textsize=30):
    exposure_dict_nz = sorted(
        [(k, v) for (k, v) in exposure_dict.items() if v > 0],
        key=(lambda x: x[1]), reverse=True)
    x = [i[0] for i in exposure_dict_nz]
    y = [i[1] for i in exposure_dict_nz]

    ax.set_title('Exposures', size=textsize, pad=20)
    ax.tick_params('x', labelrotation=45, labelsize=15)
    ax.tick_params('y', labelsize=30)
    ax.barh(x, y)


def draw_pie_sbs6(ax, catalogue_dict_sbs6, catalogue_keys_sbs6, textsize=20):
    def autopct_func(pct):
        count = int(total_count * pct * 0.01)
        pct_str = str(round(pct, 1)) + '%'
        return f'{count:,}\n({pct_str})'

    labels = catalogue_keys_sbs6
    counts = [catalogue_dict_sbs6[k] for k in catalogue_keys_sbs6]
    total_count = sum(counts)
    colors = [sigmisc.COLORS_SBS6[k] for k in catalogue_keys_sbs6]

    ax.set_xticks([])
    ax.set_yticks([])
    ax.pie(counts, labels=labels, colors=colors, autopct=autopct_func, 
           pctdistance=1.3, labeldistance=1.6, startangle=90, 
           counterclock=False, textprops={'size': textsize})


def write_texts(ax, cossim):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(0.5, 0.5, f'cosine similarity = {cossim}', 
            size=20, transform=ax.transAxes, ha='center', va='center')


def draw_onesample(sigresult, title=None):
    # set ipython magic
    IPYTHON.magic('matplotlib inline')

    # get catalogue keys
    catalogue_keys_sbs96 = sigmisc.create_catalogue_keys_sbs96(as_tuple=False)
    catalogue_keys_sbs6 = sigmisc.create_catalogue_keys_sbs6(as_tuple=False)

    # sanity check
    if set(catalogue_keys_sbs96) != set(sigresult.catalogue.index):
        raise Exception(f'Input sigresult mutation type is not SBS96.')

    # prepare data
    catalogue_original = sigresult.catalogue
    exposure = sigresult.exposure
    cossim = sigresult.cossim

    catalogue_reconst = np.matmul(sigresult.sigdata, sigresult.exposure)
    catalogue_residual = catalogue_original - catalogue_reconst
    catalogue_dict_sbs6 = make_sbs6_dict_from_sbs96(catalogue_original)
    maxcount = max(map(max, [catalogue_original, 
                             catalogue_reconst,
                             catalogue_residual]))
    ylim = (0, maxcount)

    # plotting
    plt.rc('axes', edgecolor='lightgray', linewidth=2)
    fig = plt.figure(figsize=(40, 20), constrained_layout=True)
    if title is not None:
        fig.suptitle(title, size=40)

    subfigs = fig.subfigures(1, 2, width_ratios=[1, 3])

    axs_left = subfigs[0].subplots(3, 1, 
                                   gridspec_kw={'height_ratios': [5, 5, 1]})
    draw_pie_sbs6(axs_left[0], catalogue_dict_sbs6, catalogue_keys_sbs6)
    draw_bars_exposure(axs_left[1], exposure.to_dict())
    write_texts(axs_left[2], cossim)

    subfigs[1].supylabel('Mutation count', size=30)
    axs_right = subfigs[1].subplots(3, 1, gridspec_kw={'hspace': 0.1})
    draw_bars_sbs(axs_right[0], catalogue_original.to_dict(), 
                  title='Original', ylim=ylim, 
                  catalogue_keys_sbs96=catalogue_keys_sbs96)
    draw_bars_sbs(axs_right[1], catalogue_reconst.to_dict(), 
                  title='Reconstructed', ylim=ylim,
                  catalogue_keys_sbs96=catalogue_keys_sbs96)
    draw_bars_sbs(axs_right[2], catalogue_residual.to_dict(), 
                  title='Residual (original - reconstructed)', ylim=ylim,
                  catalogue_keys_sbs96=catalogue_keys_sbs96)

    #plt.show()
    return fig


def draw_onesample_new(sigresult, title=None, fig=None):
    # set ipython magic
    IPYTHON.magic('matplotlib inline')

    # get catalogue keys
    catalogue_keys_sbs96 = sigmisc.create_catalogue_keys_sbs96(as_tuple=False)
    catalogue_keys_sbs6 = sigmisc.create_catalogue_keys_sbs6(as_tuple=False)

    # sanity check
    if set(catalogue_keys_sbs96) != set(sigresult.catalogue.index):
        raise Exception(f'Input sigresult mutation type is not SBS96.')

    # prepare data
    catalogue_original = sigresult.catalogue
    exposure = sigresult.exposure
    cossim = sigresult.cossim

    catalogue_reconst = np.matmul(sigresult.sigdata, sigresult.exposure)
    catalogue_residual = catalogue_original - catalogue_reconst
    catalogue_dict_sbs6 = make_sbs6_dict_from_sbs96(catalogue_original)
    maxcount = np.concatenate(
        [catalogue_original, catalogue_reconst, catalogue_residual],
    ).max()
    ylim = (0, maxcount)

    # plotting
    #plt.rc('axes', edgecolor='lightgray', linewidth=2)
    if fig is None:
        fig = plt.figure(figsize=(40, 20), constrained_layout=True)

    subfigs = fig.subfigures(1, 2, width_ratios=[1, 3])

    axs_left = subfigs[0].subplots(
        3, 1, 
        gridspec_kw={'height_ratios': [5, 5, 1]},
    )
    #subfigs[1].supylabel('Mutation count', size=30)
    axs_right = subfigs[1].subplots(
        3, 1, 
        gridspec_kw={'hspace': 0.1},
    )

    axes_dict['pie'] = axs_left[0]
    axes_dict['exposures'] = axs_left[1]
    axes_dict['cossim'] = axs_left[2]
    axes_dict['sbs96_original'] = axs_right[0]
    axes_dict['sbs96_reconst'] = axs_right[1]
    axes_dict['sbs96_residual'] = axs_right[2]

    if title is not None:
        fig.suptitle(title, size=40)

    draw_pie_sbs6(axes_dict['pie'], catalogue_dict_sbs6, catalogue_keys_sbs6)
    draw_bars_exposure(axes_dict['exposures'], exposure.to_dict())
    write_texts(axes_dict['cossim'], cossim)

    draw_bars_sbs(axes_dict['sbs96_original'], catalogue_original.to_dict(), 
                  title='Original', ylim=ylim, 
                  catalogue_keys_sbs96=catalogue_keys_sbs96)
    draw_bars_sbs(axes_dict['sbs96_reconst'], catalogue_reconst.to_dict(), 
                  title='Reconstructed', ylim=ylim,
                  catalogue_keys_sbs96=catalogue_keys_sbs96)
    draw_bars_sbs(axes_dict['sbs96_residual'], catalogue_residual.to_dict(), 
                  title='Residual (original - reconstructed)', ylim=ylim,
                  catalogue_keys_sbs96=catalogue_keys_sbs96)

    #plt.show()
    return fig



##########################################
# drawing signatures of multiple samples #
##########################################



def stacked_barplot(
    exposure_list, 

    # merged signatures
    sigs_to_merge=list(),

    # signatures order
    sigorder_plotting=None,
    sigorder_legend=None,

    # maps
    colors=None,  # dict: signame -> color
    legend_labels=None,  # dict: signame -> legend label

    # x axis label
    xlabels=None,

    # others
    title=None,
    figsize=None,
    width_ratios=[1, 0.5],
    legend_kwargs=dict(),
    omit_legend=False,

    axd=None,
    xs=None,
    bar_width=None,
    align='center',
    rotation=0,
    zorder=3,
    dont_set_xticks=False,
):
    """
        Args examples:
            sigs_to_merge = [
                ('SBS2', 'SBS13'),
                ('SBS1', 'SBS5'),
            ]
            sigorder_legend = None #[f'SBS{x}' for x in ['1+5', '2+13', '4', '7a', '18']]
            legend_labels = {
                'SBS1+5': 'SBS1+5 (clock-like)',
                'SBS4': 'SBS4 (smoking)',
                'SBS18': 'SBS18 (ROS)',
                'SBS2+13': 'SBS2+13 (APOBEC)',
                'SBS36': 'SBS36 (MUTYH)',
                'SBS7a': 'SBS7a (UV)',
            }
    """

    # sanitycheck 1
    assert all(isinstance(x, pd.Series) for x in exposure_list)
    if xlabels is not None:
        assert len(xlabels) == len(exposure_list)

    # merge sig exp values
    merged_exposure_list = list()
    for exp in exposure_list:
        merged_exp = sigmisc.merge_signature(exp, sigs_to_merge)
        merged_exposure_list.append(merged_exp)

    # make sum of all samples
    all_sig_exps = functools.reduce(operator.add, merged_exposure_list)
    all_sig_exps.sort_values(ascending=True, inplace=True)
    all_sig_exps = all_sig_exps[all_sig_exps > 0]

    all_signames = set(all_sig_exps.index)

    # sanitycheck 2
    if sigorder_plotting is not None:
        #assert set(sigorder_plotting) == all_signames
        sigorder_plotting = list(sigorder_plotting)
    if sigorder_legend is not None:
        assert set(sigorder_legend) == all_signames
    if colors is not None:
        assert all_signames.issubset(colors.keys())
    #if legend_labels is not None:
    #    assert set(legend_labels.keys()) == all_signames

    # signature orders
    if sigorder_plotting is None:
        sigorder_plotting = list(all_sig_exps.index)  # signames with lower total amount comes at bottom
    if sigorder_legend is None:
        sigorder_legend = sorted(all_signames, key=sigmisc.signame_sortkey)

    # maps
    if colors is None:
        colors = dict(
            zip(
                sigorder_legend, 
                mpl.cm.tab10(np.linspace(0, 1, len(sigorder_legend))),
            )
        )
    if legend_labels is None:
        legend_labels = {name: name for name in all_signames}
    else:
        for key in all_signames.difference(legend_labels.keys()):
            legend_labels[key] = key

    # data for plotting
    def plotting_signame_sortkey(signame):
        if signame in sigorder_plotting:
            return sigorder_plotting.index(signame)
        else:
            return max(len(all_signames), len(sigorder_plotting))

    real_sigorder_plotting = sorted(all_signames, key=plotting_signame_sortkey)
    data_src = list()
    for signame in real_sigorder_plotting:
        data_src.append(
            [x[signame] for x in merged_exposure_list]
        )
    data = np.stack(data_src, axis=0)

    # plotting
    if xlabels is None:
        xlabels = np.arange(len(merged_exposure_list)).astype(str)

    if axd is None:
        if not omit_legend:
            fig, axd = plt.subplot_mosaic(
                [['main', 'legend']],
                figsize=figsize, 
                gridspec_kw=dict(width_ratios=width_ratios, wspace=0.2),
            )
        else:
            fig, axd = plt.subplot_mosaic(
                [['main']],
                figsize=figsize, 
                gridspec_kw=dict(),
            )
    else:
        fig = next(iter(axd.values())).figure

    # draw main
    _ = plotmisc.vstacked_barplot(
        data, 
        xlabels=xlabels, 
        colors=[colors[x] for x in real_sigorder_plotting],
        ax=axd['main'],
        xs=xs,
        bar_width=bar_width,
        rotation=rotation,
        align=align,
        zorder=zorder,
        dont_set_xticks=dont_set_xticks,
    )

    axd['main'].set_ylabel('proportion of signature')
    if title is None:
        title = f'Compositions of mutational signatures'
    fig.suptitle(title)

    # draw legend
    if not omit_legend:
        plotmisc.clear_ax(axd['legend'])
        leghandles = plotmisc.LegendHandles()
        for signame in sigorder_legend:
            leghandles.add_patch(color=colors[signame], label=legend_labels[signame])
        axd['legend'].legend(
            handles=leghandles, 
            bbox_to_anchor=(0.5, 0.5),
            loc='center',
            **legend_kwargs,
        )

    return fig, axd, sigorder_plotting



import re

import IPython
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
signature_misc = importlib.import_module('.'.join([top_package_name, 'signature', 'misc']))


IPYTHON = IPython.get_ipython()
SBS96_KEY_PAT = re.compile('^(.*)\[((.*)>(.*))\](.*)$')


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
        color = signature_misc.COLORS_SBS6[f'{ref}>{alt}']
        ax.bar(x, y, color=color, align='center', width=0.4, zorder=1000)
        if add_cat_labels:
            ax.text(x, 0, catalogue_key, size=15, rotation='vertical', 
                    ha='center', va='top', family='monospace', weight='bold')


def draw_bars_exposure(ax, exposure_dict):
    exposure_dict_nz = sorted(
        [(k, v) for (k, v) in exposure_dict.items() if v > 0],
        key=(lambda x: x[1]), reverse=True)
    x = [i[0] for i in exposure_dict_nz]
    y = [i[1] for i in exposure_dict_nz]

    ax.set_title('Exposures', size=30, pad=20)
    ax.tick_params('x', labelrotation=45, labelsize=15)
    ax.tick_params('y', labelsize=30)
    ax.barh(x, y)


def draw_pie_sbs6(ax, catalogue_dict_sbs6, catalogue_keys_sbs6):
    def autopct_func(pct):
        count = int(total_count * pct * 0.01)
        pct_str = str(round(pct, 1)) + '%'
        return f'{count:,}\n({pct_str})'

    labels = catalogue_keys_sbs6
    counts = [catalogue_dict_sbs6[k] for k in catalogue_keys_sbs6]
    total_count = sum(counts)
    colors = [signature_misc.COLORS_SBS6[k] for k in catalogue_keys_sbs6]

    ax.set_xticks([])
    ax.set_yticks([])
    ax.pie(counts, labels=labels, colors=colors, autopct=autopct_func, 
           pctdistance=1.3, labeldistance=1.6, startangle=90, 
           counterclock=False, textprops={'size': 20})


def write_texts(ax, cossim):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(0.5, 0.5, f'cosine similarity = {cossim}', 
            size=20, transform=ax.transAxes, ha='center', va='center')


def main(sigresult, sampleid):
    # set ipython magic
    IPYTHON.magic('matplotlib inline')

    # get catalogue keys
    catalogue_keys_sbs96 = signature_misc.create_catalogue_keys_sbs96(as_tuple=False)
    catalogue_keys_sbs6 = signature_misc.create_catalogue_keys_sbs6(as_tuple=False)

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
    fig.suptitle(sampleid, size=40)

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

    plt.show()


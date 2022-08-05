import collections

import numpy as np
import matplotlib.pyplot as plt

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))


def show_readcounts(vp, alleleindex=1, sampleid_order=None, sort_by_value=False, title=None,
                    mask=None, figsize=(20, 5),
                    count_ylim=None, vaf_ylim=None,
                    logscale_count=False, nan_value=-0.1):
    # sampleid_order argument sanity check
    if sampleid_order is not None:
        if not set(sampleid_order).issubset(vp.readstats_dict.keys()):
            raise Exception(f'"sampleid_order" must be a subset of the '
                            f'keys of "readstats_dict".')
    # set sampleid order
    if sampleid_order is None:
        sampleid_order = sorted(vp.readstats_dict.keys())
    if sort_by_value:
        items = [(key, vp.readstats_dict[key]['rppcounts'][alleleindex]) for key in sampleid_order]
        sampleid_order = [x[0] for x in sorted(items, key=(lambda x: x[1]))]

    # prepare data
    altcnts = np.array([
        vp.readstats_dict[sampleid]['rppcounts'][alleleindex]
        for sampleid in sampleid_order])
    vaf = np.array([
        vp.get_vaf(sampleid, alleleindex=alleleindex)
        for sampleid in sampleid_order])  # may include nan
    vaf[np.isnan(vaf)] = nan_value  # nan is substituted with a negative value
    
    # prepare mask
    if mask is None:
        mask = {key: True for key in sampleid_order}
    not_mask = {key: (not val) for key, val in mask.items()}

    mask_array = np.array([mask[key] for key in sampleid_order])
    not_mask_array = np.array([not_mask[key] for key in sampleid_order])
    
    x = np.arange(len(sampleid_order))
    y_maps = {
        'count': {
            'value': altcnts, 
            'label': 'ALT read count',
            },
        'vaf': {
            'value': vaf, 
            'label': 'vaf',
            },
        }
    y1_map = y_maps['count']
    y2_map = y_maps['vaf']
    
    y1 = y1_map['value']
    ylabel1 = y1_map['label']
    y2 = y2_map['value']
    ylabel2 = y2_map['label']    

    # plot parameters
    col_ax1 = 'tab:blue'
    col_good_bar = 'tab:blue'
    col_bad_bar = 'tab:red'
    col_ax2 = 'tab:red'
    col_dot = 'tab:orange'

    # ax1
    fig, ax1 = plt.subplots(figsize=figsize)

    ax1.bar(x[mask_array], y1[mask_array], color=col_good_bar, alpha=0.7)
    ax1.bar(x[not_mask_array], y1[not_mask_array], color=col_bad_bar, alpha=0.7)

    ax1.set(xticks=np.arange(len(sampleid_order)), 
            xticklabels=list(sampleid_order))
    ax1.set_xlabel('sample ID')
    ax1.tick_params(axis='x', direction='out', labelrotation=90)
    ax1.set_ylabel(f'{ylabel1} (BARPLOT)', color=col_ax1)
    ax1.tick_params(axis='y', color=col_ax1, labelcolor=col_ax1)
    if title is not None:
        ax1.set(title=title)
    if count_ylim is not None:
        ax1.set(ylim=count_ylim)
    if logscale_count:
        ax1.set_yscale('log')
        
    # ax2
    ax2 = ax1.twinx()
    ax2.scatter(x, y2, color=col_dot)
    if vaf_ylim is not None:
        ax2.set(ylim=vaf_ylim)
    ax2.set_ylabel(f'{ylabel2} (SCATTERPLOT)', color=col_ax2)
    ax2.tick_params(axis='y', color=col_ax2, labelcolor=col_ax2)

    plt.show()
    
    
###########################################

def show_readstats_data(readstats_data, title=None, varpos_key='left'):
    fig, axs = plt.subplots(2, 6, figsize=(40, 8))
    
    # title
    if title is not None:
        fig.suptitle(title)
        
    # MQ
    plotter_MQ(ax=axs[0, 0], data=readstats_data['MQ'][1], 
               title='MQ - ALT reads')
    plotter_MQ(ax=axs[1, 0], data=readstats_data['MQ'][0], 
               title='MQ - REF reads')    
    # BQ
    plotter_BQ(ax=axs[0, 1], data=readstats_data['BQ'][1], 
               title='BQ - ALT reads')
    plotter_BQ(ax=axs[1, 1], data=readstats_data['BQ'][0], 
               title='BQ - REF reads')
    # pariorient
    plotter_pairorient(ax=axs[0, 2], data=readstats_data['pairorient'][1], 
                       title='pairorient - ALT reads')
    plotter_pairorient(ax=axs[1, 2], data=readstats_data['pairorient'][0], 
                       title='pairorient - REF reads')    
    # cliplen
    plotter_cliplen(ax=axs[0, 3], data=readstats_data['cliplen'][1], 
                    title='Softclip length - ALT reads')
    plotter_cliplen(ax=axs[1, 3], data=readstats_data['cliplen'][0], 
                    title='Softclip length - REF reads')
    # variant position
    fraction = varpos_key.endswith('_fraction')
    plotter_varpos(
        ax=axs[0, 4],
        data=readstats_data[f'pos0_{varpos_key}'][1],
        title=f'variant position from {varpos_key} end - ALT reads',
        fraction=fraction)
    plotter_varpos(
        ax=axs[1, 4],
        data=readstats_data[f'pos0_{varpos_key}'][0],
        title=f'variant position from {varpos_key} end - REF reads',
        fraction=fraction)
    # read counts
    plotter_counts(ax=axs[0, 5], data=readstats_data['count'])
    plotter_counts_zoom(ax=axs[1, 5], data=readstats_data['count'])

    plt.show()
    
    
# show_readstats helpers
    
def plotter_hist(ax, data, title, xlim=None, binwidth=None, fontsize=18, 
                 text_mean=True, text_median=True):
    ax.set_title(title)
    
    if xlim is not None:
        ax.set_xlim(xlim[0], xlim[1])
        
    if (xlim is not None) and (binwidth is not None):
        bins = np.arange(xlim[0], xlim[1] + binwidth, binwidth)    
    else:
        bins = None
    ax.hist(data, bins=bins)
        
    if text_mean:
        mean = np.mean(data)
        ax.text(0.05, 0.9, f'mean: {mean}', transform=ax.transAxes, 
                fontsize=fontsize)
    if text_median:
        median = np.median(data)
        ax.text(0.05, 0.8, f'median: {median}', transform=ax.transAxes, 
                fontsize=fontsize)
        
def plotter_bar(ax, data, title, indexes, ylim=None, text=None, fontsize=18, 
                xlabelrotation=0):
    ax.set_title(title)
    ax.tick_params(axis='x', labelrotation=xlabelrotation)
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])
        
    ax.bar([str(x) for x in indexes], [data[x] for x in indexes])
    if text is not None:
        ax.text(0.05, 0.9, text, transform=ax.transAxes, fontsize=fontsize)
        
###
        
def plotter_pairorient(ax, data, title):
    data_counter = collections.Counter(data)
    if data_counter['f2r1'] == 0:
        ratio = np.nan
    else:
        ratio = np.round(data_counter['f1r2'] / data_counter['f2r1'], 3)
    text = f'f1r2/f2r1 ratio: {ratio}'
    plotter_bar(ax, data_counter, title, indexes=['f1r2', 'f2r1'], text=text)
        
def plotter_MQ(ax, data, title):
    plotter_hist(ax, data, title, xlim=(0, 70), binwidth=1, 
                 text_mean=True, text_median=True)
    
def plotter_BQ(ax, data, title):
    plotter_hist(ax, data, title, xlim=(0, 40), binwidth=1, 
                 text_mean=True, text_median=True)
    
def plotter_cliplen(ax, data, title):
    plotter_hist(ax, data, title, xlim=(0, max(data) + 10), binwidth=None, 
                 text_mean=True, text_median=True)
    
def plotter_varpos(ax, data, title, fraction=False):
    if fraction:
        xlim = (0, 1)
        binwidth = 0.05
    else:
        xlim = (0, 155)
        binwidth = 5
    plotter_hist(ax, data, title, xlim=xlim, binwidth=binwidth, 
                 text_mean=True, text_median=True)
    
def plotter_counts(ax, data):
    data_counter = collections.Counter(data)
    total_count = get_total_rppcount(data_counter)
    data_counter['total'] = total_count
    title = 'Read counts per alleleclass'
    plotter_bar(ax, data_counter, title, 
                indexes=['total', 0, 1, -1, None, 'softclip_overlap'], 
                ylim=None, text=None, xlabelrotation=45)
    
def plotter_counts_zoom(ax, data):
    data_counter = collections.Counter(data)
    total_count = get_total_rppcount(data_counter)
    data_counter['total'] = total_count
    title = 'Read counts per alleleclass, zoom on ALT and softclip-overlaps'
    
    ylim_max = int(max(data_counter[1], data_counter['softclip_overlap']) * 2)
    ylim = (0, ylim_max)
    
    ratio = round(data_counter['softclip_overlap'] / data_counter[1], 3)
    text = f'softclip_overlap / ALT read count ratio: {ratio}'
    
    plotter_bar(ax, data_counter, title, 
                indexes=['total', 0, 1, -1, None, 'softclip_overlap'], 
                ylim=ylim, text=text, xlabelrotation=45)

def get_total_rppcount(rppcount_data_counter):
    return sum(val for (key, val) in rppcount_data_counter.items()
               if isinstance(key, int))


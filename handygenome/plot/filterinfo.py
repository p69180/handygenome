import collections

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import handygenome.common as common
import handygenome.workflow as workflow
import handygenome.variant.ponbams as libponbams
import handygenome.variant.filter as libfilter
import handygenome.deco as deco


def show_readcounts(
    readstats_dict,
    allele_index=1, 
    samples=None, 
    sample_label_colors=None,
    sort_by_value=False, 
    exclude_other=False,
    mask=None, 
    xlabel='sample ID',
    title=None,
    title_pad=None,

    bar_data_type='count',
    dot_data_type='vaf',

    figsize=(20, 7),
    bar_ylim=None, 
    dot_ylim=(0, 1),
    logscale_count=False, 
    nan_value=0,
):
    """Args:
        samples: A sequence of sample IDs, which determines plotting order
        sort_by_value: If True, samples are further sorted by ALT rppcount value, in increasing order
        mask: A dictionary (keys: sample IDs, values: True or False); True bars are colored blue, False bars are red
    """
    # samples argument sanity check
    if samples is not None:
        if not set(samples).issubset(readstats_dict.keys()):
            raise Exception(f'"samples" must be a subset of the keys of "readstats_dict".')

    # handle samples argument
    if samples is None:
        samples = list()
        for x in sorted(readstats_dict.keys()):
            if not readstats_dict[x].is_missing:
                samples.append(x)
    if sort_by_value:
        samples = sorted(
            samples, 
            key=(lambda x: readstats_dict[x]['rppcounts'][allele_index])
        )
    # handle sample_label_colors argument
    if sample_label_colors is None:
        sample_label_colors = np.repeat('black', len(samples))

    # plot parameters
    col_bar = 'tab:blue'
    col_bar_label = col_bar #'tab:blue'
    col_bar_true = col_bar #'tab:blue'
    col_bar_false = 'tab:red'

    col_dot = 'tab:purple'
    col_dot_label = col_dot #'tab:red'
    col_dot_nan = 'tab:red'

    # prepare mask
    if mask is None:
        mask = {sid: True for sid in samples}
    not_mask = {key: (not val) for key, val in mask.items()}

    mask_array = np.array([mask[key] for key in samples])
    not_mask_array = np.array([not_mask[key] for key in samples])

    # prepare data
    target_rpp_counts = np.array(
        [
            readstats_dict[sid]['rppcounts'][allele_index]
            for sid in samples
        ]
    )
    total_rpp_counts = np.array(
        [
            readstats_dict[sid].get_total_rppcount(exclude_other=exclude_other)
            for sid in samples
        ]
    )
    nontarget_rpp_counts = total_rpp_counts - target_rpp_counts
    vafs = np.array(
        [
            readstats_dict[sampleid].get_vaf(
                alleleclass=allele_index, 
                exclude_other=exclude_other,
            )
            for sampleid in samples
        ]
    )
    #vafs[np.isnan(vafs)] = nan_value  # nan is substituted with a negative value

    y_candidates = {
        'count': {
            'value': target_rpp_counts, 
            'label': 'ALT read count',
        },
        'vaf': {
            'value': vafs, 
            'label': 'vaf',
        },
    }
    
    # prepare x, y axis values
    bar_data_type = 'count'
    dot_data_type = 'vaf'

    x = np.arange(len(samples))

    #y_bar = y_candidates[bar_data_type]['value']
    #y_bar_label = y_candidates[bar_data_type]['label']
    y_bar_target = target_rpp_counts
    y_bar_nontarget = nontarget_rpp_counts
    y_bar_label = f'Read pair count (BAR)\n(counts with relevant allele index is colored)'

    #y_dot = y_candidates[dot_data_type]['value']
    #y_dot_label = y_candidates[dot_data_type]['label']
    y_dot = vafs
    y_dot_label = f'VAF (SCATTER)'

    # init figure
    fig, ax_bar = plt.subplots(figsize=figsize)

    # ax_bar
    ax_bar.set(xticks=np.arange(len(samples)), xticklabels=list(samples))
    for label, color in zip(ax_bar.get_xticklabels(), sample_label_colors):
        label.set_color(color)

    # true bars
    width = 0.6
    ax_bar.bar(x[mask_array], y_bar_nontarget[mask_array], edgecolor='black', linewidth=1, fill=False, alpha=0.7, width=width)
    ax_bar.bar(x[mask_array], y_bar_target[mask_array], bottom=y_bar_nontarget[mask_array], edgecolor='black', linewidth=1, color=col_bar_true, alpha=0.7, width=width)
    # false bars
    ax_bar.bar(x[not_mask_array], y_bar_nontarget[not_mask_array], edgecolor='black', linewidth=1, fill=False, alpha=0.7, width=width)
    ax_bar.bar(x[not_mask_array], y_bar_target[not_mask_array], bottom=y_bar_nontarget[not_mask_array], edgecolor='black', linewidth=1, color=col_bar_false, alpha=0.7, width=width)
    #ax_bar.bar(x[not_mask_array], y_bar[not_mask_array], color=col_bar_false, alpha=0.7)

    ax_bar.set_xlabel(xlabel, size=15)
    ax_bar.tick_params(axis='x', direction='out', labelrotation=90)
    #ax_bar.set_ylabel(f'{y_bar_label} (BAR)', color=col_bar_label)
    ax_bar.set_ylabel(y_bar_label, color=col_bar_label)
    ax_bar.tick_params(axis='y', color=col_bar_label, labelcolor=col_bar_label)

    if title is not None:
        ax_bar.set_title(title, pad=title_pad)

    if bar_ylim is not None:
        ax_bar.set(ylim=bar_ylim)
    else:
        ax_bar.set(
            ylim=(0, np.ceil(max(total_rpp_counts)))
        )

    if logscale_count:
        ax_bar.set_yscale('log')
        
    # ax_dot
    ax_dot = ax_bar.twinx()
    y_isnan = np.isnan(y_dot)
    y_for_nan = np.repeat(nan_value, y_isnan.sum())

    markersize = 8
    ax_dot.plot(
        x[~y_isnan], 
        y_dot[~y_isnan], 
        color=col_dot,
        marker='o',
        markersize=markersize,
        linestyle='',
    )
    ax_dot.plot(
        x[y_isnan], 
        y_for_nan,
        color=col_dot_nan,
        marker='x',
        markersize=markersize,
        linestyle='',
    )

    if dot_ylim is not None:
        ax_dot.set(ylim=dot_ylim)
    else:
        ax_dot.set(ylim=(-0.03, max(y_dot) * 1.2))

    #ax_dot.set_ylabel(f'{y_dot_label} (SCATTER)', color=col_dot_label)
    ax_dot.set_ylabel(y_dot_label, color=col_dot_label)
    ax_dot.tick_params(axis='y', color=col_dot_label, labelcolor=col_dot_label)

    return fig, ax_bar, ax_dot


@deco.get_deco_num_set_differently(('pon_samples', 'pon_cohorts'), 1)
def show_pon(
    query_sample, 
    vp, 
    allele_index=1, 
    exclude_query=True, 
    exclude_other=False,
    pon_cohorts=None, 
    pon_samples=None, 
    ponfilter_mode='WGS',
    **kwargs,
):
    assert ponfilter_mode in ('WGS', 'PanelGermline', 'PanelSomatic')
    assert (
        (pon_cohorts is None and pon_samples is not None) or 
        (pon_cohorts is not None and pon_samples is None)
    )
    assert isinstance(pon_cohorts, (list, tuple, type(None))), f'"pon_cohorts" must be a list or a tuple.'
    assert isinstance(pon_samples, (list, tuple, type(None))), f'"pon_samples" must be a list or a tuple.'

    # set pon_samples from cohort name
    if pon_cohorts is not None:
        pon_samples = libponbams.get_pon_sample_names(pon_cohorts, vp.refver)
    else:
        pon_samples = list(pon_samples)

    # remove query sample from pon samples
    if exclude_query:
        if query_sample in pon_samples:
            pon_samples.remove(query_sample)

    # sanity check
    if query_sample not in vp.readstats_dict.keys():
        raise Exception(f'"query_sample" is not included in readstats_dict')
    if not set(pon_samples).issubset(vp.readstats_dict.keys()):
        raise Exception(f'"pon_samples" is not included in readstats_dict')

    # main
    # run ponfilter
    if ponfilter_mode == 'WGS':
        ponfilter = libfilter.PonFilterWGS(samples=pon_samples, save_cache=True, verbose=False, **kwargs)
    elif ponfilter_mode == 'PanelGermline':
        ponfilter = libfilter.PonFilterPanelseqGermline(samples=pon_samples, save_cache=True, verbose=False, **kwargs)
    elif ponfilter_mode == 'PanelSomatic':
        ponfilter = libfilter.PonFilterPanelseqSomatic(samples=pon_samples, save_cache=True, verbose=False, **kwargs)

    ponfilter.check(vp, query_sample, allele_index=allele_index, exclude_query=exclude_query)
    
    # prepare params for show_readcounts
    title = '\n'.join([
        #f'{repr(vp.vcfspec).strip()}',
        f'{vp.vcfspec}',
        ' ; '.join([
            f'ponfilter_mode: {ponfilter_mode}',
            f'bisect_cutoff: {ponfilter.params["bisect_cutoff"]}',
            f'allele_index: {allele_index}',
        ]),
        ' ; '.join([
            f'lower_different: {ponfilter.cache["lower_different"]}',
            f'upper_different: {ponfilter.cache["upper_different"]}',
            f'lower_sample_fraction: {str(ponfilter.cache["lower_sample_fraction"])[:4]}',
            f'upper_sample_fraction: {str(ponfilter.cache["upper_sample_fraction"])[:4]}',
        ]),
        ' ; '.join([
            f'is_germline: {ponfilter.cache["is_germline"]}',
            f'het_vaf_mean: {str(ponfilter.cache["het_vaf_mean"])[:4]}',
            f'germline_fraction_among_upper: {str(ponfilter.cache["germline_fraction_among_upper"])[:4]}',
            f'germline_fraction_among_all: {str(ponfilter.cache["germline_fraction_among_all"])[:4]}',
        ]),
#        ' ; '.join([
#            f'is_global: {ponfilter.cache["is_global"]}',
#        ]),
        ' ; '.join([
            f'total_result: {ponfilter.cache["total_result"]}',
        ]),
    ])
    all_samples = [query_sample] + pon_samples
    sample_label_colors = [
        'red' if x in ponfilter.cache['samples_with_low_depth'] else 'black'
        for x in all_samples
    ]

    if ponfilter.cache['upper_cutoff'] is not None:
        dot_ylim = (0, max(ponfilter.cache['upper_cutoff'] + 0.1, 1))
    else:
        dot_ylim = (0, 1)

    fig, ax_bar, ax_dot = show_readcounts(
        vp.readstats_dict, 
        allele_index=allele_index,
        samples=all_samples,
        sample_label_colors=sample_label_colors,
        exclude_other=exclude_other,
        xlabel='sample ID (low-depth samples in red)',
        title=title,
        title_pad=30,
        dot_ylim=dot_ylim,
    )

    # Partition between query and PON samples 
    ax_bar.axvline(x=0.5, color='black')
    #offset = 0.15
    ax_bar.annotate(
        'query sample', 
        xy=(0, ax_bar.get_ylim()[1] * 1.01), 
        size=10,
        annotation_clip=False,
        ha='right',
        va='bottom',
        bbox=dict(boxstyle='round', facecolor='white', alpha=1),
    )
    ax_bar.annotate(
        'PON samples', 
        #xy=(0.5 + offset, ax_bar.get_ylim()[1] * 0.95), 
        xy=((1 + len(all_samples)) / 2, ax_bar.get_ylim()[1] * 1.01), 
        size=10,
        annotation_clip=False,
        ha='center',
        va='bottom',
        bbox=dict(boxstyle='round', facecolor='white', alpha=1),
    )

    # cutoff hlines
    ax_dot.axhline(y=ponfilter.params['bisect_cutoff'], color='black', linestyle='--')
    if ponfilter.cache['lower_mean'] is not None:
        ax_dot.axhline(y=ponfilter.cache['lower_mean'], color='orangered')
        ax_dot.axhline(y=ponfilter.cache['lower_cutoff'], color='lime')
    if ponfilter.cache['upper_mean'] is not None:
        ax_dot.axhline(y=ponfilter.cache['upper_mean'], color='orangered')
        ax_dot.axhline(y=ponfilter.cache['upper_cutoff'], color='lime')

    for text, y in [
        ('lower vafs mean', ponfilter.cache['lower_mean']),
        ('lower vafs cutoff', ponfilter.cache['lower_cutoff']),
        ('upper vafs mean', ponfilter.cache['upper_mean']),
        ('upper vafs cutoff', ponfilter.cache['upper_cutoff']),
    ]:
        if y is None:
            continue

        ax_dot.annotate(
            text, 
            xy=(ax_dot.get_xlim()[1] * 1.05, y), 
            size=10,
            annotation_clip=False,
            ha='left',
            va='center',
            bbox=dict(boxstyle='round', facecolor='white', alpha=1),
        )

    ax_dot_xlim = ax_dot.get_xlim()
    rectangles = list()
    if 'germline_het_range' in ponfilter.params:
        rectangles.append(
            mpl.patches.Rectangle(
                (ax_dot_xlim[0], ponfilter.params['germline_het_range'][0]),
                ax_dot_xlim[1] - ax_dot_xlim[0],
                ponfilter.params['germline_het_range'][1] - ponfilter.params['germline_het_range'][0],
            ),
        )
    if 'germline_hom_range' in ponfilter.params:
        rectangles.append(
            mpl.patches.Rectangle(
                (ax_dot_xlim[0], ponfilter.params['germline_hom_range'][0]),
                ax_dot_xlim[1] - ax_dot_xlim[0],
                ponfilter.params['germline_hom_range'][1] - ponfilter.params['germline_hom_range'][0],
            ),
        )
    if len(rectangles) > 0:
        ax_dot.add_collection(
            mpl.collections.PatchCollection(
                rectangles,
                alpha=0.1,
                facecolor='red',
                linewidth=0,
            )
        )

    return fig
    
    
###########################################

def show_readstats_data(readstats_data, alt_index=0, title=None, varpos_key='left'):
    fig, axs = plt.subplots(
        2, 6, 
        figsize=(40, 15),
        gridspec_kw={'hspace': 0.4},
    )
    #fig.tight_layout()
    
    # title
    if title is not None:
        fig.suptitle(title)
        
    # MQ
    plotter_MQ(ax=axs[0, 0], data=readstats_data['MQ'][alt_index + 1], 
               title='MQ - ALT reads')
    plotter_MQ(ax=axs[1, 0], data=readstats_data['MQ'][0], 
               title='MQ - REF reads')    
    # BQ
    plotter_BQ(ax=axs[0, 1], data=readstats_data['BQ'][alt_index + 1], 
               title='BQ - ALT reads')
    plotter_BQ(ax=axs[1, 1], data=readstats_data['BQ'][0], 
               title='BQ - REF reads')
    # pariorient
    plotter_pairorient(ax=axs[0, 2], data=readstats_data['pairorient'][alt_index + 1], 
                       title='pairorient - ALT reads')
    plotter_pairorient(ax=axs[1, 2], data=readstats_data['pairorient'][0], 
                       title='pairorient - REF reads')    
    # cliplen
    plotter_cliplen(ax=axs[0, 3], data=readstats_data['cliplen'][alt_index + 1], 
                    title='Softclip length - ALT reads')
    plotter_cliplen(ax=axs[1, 3], data=readstats_data['cliplen'][0], 
                    title='Softclip length - REF reads')
    # variant position
    fraction = varpos_key.endswith('_fraction')
    plotter_varpos(
        ax=axs[0, 4],
        data=readstats_data[f'pos0_{varpos_key}'][alt_index + 1],
        title=f'variant position from {varpos_key} end\nALT reads',
        fraction=fraction)
    plotter_varpos(
        ax=axs[1, 4],
        data=readstats_data[f'pos0_{varpos_key}'][0],
        title=f'variant position from {varpos_key} end\nREF reads',
        fraction=fraction)
    # read counts
    plotter_counts(ax=axs[0, 5], data=readstats_data['count'])
    plotter_counts_zoom(ax=axs[1, 5], data=readstats_data['count'])

    #plt.show()
    return fig

    
    
# show_readstats helpers
    
def plotter_hist(
    ax, data, title, xlim=None, binwidth=None, 
    fontsize_annot=18, fontsize_title=20,
    text_mean=True, text_median=True,
):
    ax.set_title(title, size=fontsize_title, pad=10)
    
    if xlim is not None:
        ax.set_xlim(xlim[0], xlim[1])
        
    if (xlim is not None) and (binwidth is not None):
        bins = np.arange(xlim[0], xlim[1] + binwidth, binwidth)    
    else:
        bins = None
    ax.hist(data, bins=bins)
        
    if text_mean:
        mean = np.mean(data)
        ax.text(
            0.05, 0.9, f'mean: {mean}', transform=ax.transAxes, 
            fontsize=fontsize_annot,
        )
    if text_median:
        median = np.median(data)
        ax.text(
            0.05, 0.8, f'median: {median}', transform=ax.transAxes, 
            fontsize=fontsize_annot,
        )
        
def plotter_bar(
    ax, data, title, indexes, ylim=None, text=None,
    fontsize_annot=18, fontsize_title=20,
    xlabelrotation=0,
):
    ax.set_title(title, size=fontsize_title, pad=10)
    ax.tick_params(axis='x', labelrotation=xlabelrotation)
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])
        
    ax.bar([str(x) for x in indexes], [data[x] for x in indexes])
    if text is not None:
        ax.text(0.05, 0.9, text, transform=ax.transAxes, fontsize=fontsize_annot)
        
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
    title = 'Read counts per alleleclass,\nzoom on ALT and softclip-overlaps'
    
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


import numpy as np
#import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


def make_quadmesh(len_x, len_y, width_x, width_y):
    xs = np.tile(np.arange(0, (len_x + 1) * width_x, width_x), (len_y + 1, 1))
    ys = np.tile(np.arange(0, (len_y + 1) * width_y, width_y), (len_x + 1, 1)).T
    coords = np.dstack((xs, ys))
    return mpl.collections.QuadMesh(coords, edgecolors='black', facecolors='none', linewidth=0.5)


def setup_main(ax, samples, variants):
    ax.set_aspect(1)
    
    ax.set(xlim=(0, 2 * len(samples)), ylim=(0, len(variants)))
    ax.add_collection(make_quadmesh(2 * len(samples), len(variants), 1, 1))

    ax.set_xticks(np.arange(1, 2 * len(samples), 2), samples, minor=False)
    ax.tick_params(axis='x', which='major', bottom=False, rotation=90, pad=10)
    ax.set_xticks(np.arange(0.5, 2 * len(samples), 1), np.tile(['w', 'p'], len(samples)), minor=True)
    ax.tick_params(axis='x', which='minor', bottom=False, rotation=0, pad=2, labelsize=8)

    ax.set_yticks(np.arange(0.5, len(variants), 1), variant_names, va='center')
    ax.tick_params(axis='y', left=False)
    

def setup_bottom(ax, samples, title):
    ax.set_aspect(1)
    
    ax.set(xlim=(0, 2 * len(samples)), ylim=(0, 1))
    ax.add_collection(make_quadmesh(len(samples), 1, 2, 1))
    ax.tick_params(axis='both', which='both', labelbottom=False, labelleft=False, bottom=False, left=False)
    ax.set_ylabel(title, rotation=0, va='center', ha='right')

    
def setup_left(ax, variants, title):
    ax.set_aspect(1)
    
    ax.set(xlim=(0, 1), ylim=(0, len(variants)))
    ax.add_collection(make_quadmesh(1, len_y, 1, 1))
    ax.tick_params(axis='both', which='both', labelbottom=False, labelleft=False, bottom=False, left=False)
    ax.set_xlabel(title, rotation=90, va='top', ha='center')


def write_bottom(fig, ax_values, ax_cb, values, cbar_min, cbar_max):
    mappable = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=cbar_min, vmax=cbar_max), cmap=mpl.cm.cool)
    fig.colorbar(mappable, cax=ax_cb, orientation='horizontal')
    
    for t, x in zip([str(k) for k in values], np.arange(1, 2 * len(values), 2)):
        ax_values.annotate(t, (x, 0.5), ha='center', va='center', size=8)
    
    colors = mappable.to_rgba(values)
    pcol = mpl.collections.PatchCollection(
        [
            mpl.patches.Rectangle((x, 0), 2, 1, linewidth=0, color=c)
            for x, c in zip(np.arange(0, 2 * len(values), 2), colors)
        ],
        match_original=True,
    ) 
    ax_values.add_collection(pcol)
    

def get_var_id(vp):
    if not hasattr(vp, 'oncokb'):
        vp.oncokb = liboncokb.query_hgvsg(vp.vcfspec.to_hgvsg(), ONCOKB_TOKEN)
    var_id = re.sub('^.+the (.+) mutation.+$', '\\1', vp.oncokb['variantSummary'])
    return var_id

    
def fill_main(ax, wgs_driver_vplists, panel_driver_vplists, SVinfo,  samples, variants):
    patches = list()
    var_id_pat = re.compile('^.+the (.+) mutation.+$')
    
    for sampleid, x in zip(samples, np.arange(0.5, len(samples) * 2, 2)):
        x_mid_wgs = x
        x_mid_panel = x + 1
        for var_name, y_bot in zip(variant_names, np.arange(0, len(variants), 1)):
            # FUSION
            if '-' in var_name:
                gene1, gene2 = var_name.split('-')
                for dic in SVinfo['wgs'][sampleid]:
                    if gene1 in dic['genes'] and gene2 in dic['genes'] and dic['strand_aligned']:
                        patches.append(mpl.patches.Rectangle((x_mid_wgs - 0.5, y_bot), 1, 1, color='r', alpha=0.5, zorder=-1))
                        break
                for dic in SVinfo['panel'][sampleid]:
                    if gene1 in dic['genes'] and gene2 in dic['genes'] and dic['strand_aligned']:
                        patches.append(mpl.patches.Rectangle((x_mid_panel - 0.5, y_bot), 1, 1, color='r', alpha=0.5, zorder=-1))
                        break
            else:
                if var_name.startswith('EGFR'):
                    if var_name == 'EGFR others':
                        def vpmatcher(vp):
                            var_id = get_var_id(vp)
                            return ('EGFR' in var_id) and (var_id not in [    
                                'EGFR L858R',
                                'EGFR E746_A750del',
                                'EGFR A767_V769dup',])
                    else:
                        def vpmatcher(vp):
                            var_id = get_var_id(vp)
                            return var_id == var_name
                elif var_name.startswith('KRAS'):
                    if var_name == 'KRAS others':
                        def vpmatcher(vp):
                            var_id = get_var_id(vp)
                            return ('KRAS' in var_id) and (var_id not in [   
                                'KRAS G12C',
                                'KRAS G12D',
                                'KRAS G12V',                                
                            ])
                    else:
                        def vpmatcher(vp):
                            var_id = get_var_id(vp)
                            return var_id == var_name
                elif var_name.startswith('BRAF'):
                    if var_name == 'BRAF others':
                        def vpmatcher(vp):
                            var_id = get_var_id(vp)
                            return ('BRAF' in var_id) and (var_id not in [   
                                'BRAF G466V',
                            ])
                    else:
                        def vpmatcher(vp):
                            var_id = get_var_id(vp)
                            return var_id == var_name                   
                else:
                    if var_name == 'others':
                        def vpmatcher(vp):
                            var_id = get_var_id(vp)
                            return (
                                not any(x in var_id for x in ['EGFR', 'KRAS', 'BRAF',                             'PIK3CA', 
                            'PTEN',
                            'TP53', 
                            'KEAP1',
                            'ATM',]) 
                            )
                    else:
                        def vpmatcher(vp):
                            var_id = get_var_id(vp)
                            return var_name in var_id
                        
                hit = False
                for vp in wgs_driver_vplists[sampleid]:
                    if vpmatcher(vp):
                        hit = True
                        break
                if hit:
                    vaf = vp.get_vaf(f'{sampleid}_tumor')
                    if not hasattr(vp, 'ccfinfo'):
                        vp.ccfinfo = vp.get_ccf_CNm(segments_df=CNV_SEGMENTS[sampleid], is_female=IS_FEMALE[sampleid], cellularity=PREVIOUS_PURITY[sampleid], sampleid=f'{sampleid}_tumor')
                    ccf = vp.ccfinfo.ccf
                    
                    # print(vaf, ccf)
                    
                    ax.hlines(y=(y_bot + vaf), xmin=(x_mid_wgs - 0.5), xmax=(x_mid_wgs + 0.5), color='blue')
                    ax.hlines(y=(y_bot + ccf), xmin=(x_mid_wgs - 0.5), xmax=(x_mid_wgs + 0.5), color='orange')
                    
                hit = False
                for vp in panel_driver_vplists[sampleid]:
                    if vpmatcher(vp):
                        hit = True
                        break
                if hit:
                    vaf = vp.get_vaf(f'{sampleid}_panel')
                    if not hasattr(vp, 'ccfinfo'):
                        vp.ccfinfo = vp.get_ccf_CNm(segments_df=CNV_SEGMENTS[sampleid], is_female=IS_FEMALE[sampleid], cellularity=PREVIOUS_PURITY[sampleid], sampleid=f'{sampleid}_panel')
                    ccf = vp.ccfinfo.ccf
                    
                    # print(vaf, ccf)
                    
                    ax.hlines(y=(y_bot + vaf), xmin=(x_mid_panel - 0.5), xmax=(x_mid_panel + 0.5), color='blue')
                    ax.hlines(y=(y_bot + ccf), xmin=(x_mid_panel - 0.5), xmax=(x_mid_panel + 0.5), color='orange')

          
    ax.add_collection(mpl.collections.PatchCollection(patches, match_original=True))


def make_heatmap(variants, samples, purities, ploidies):
    fig = plt.figure(
        figsize=(20, 5), 
        layout='constrained',
    )

    axd = fig.subplot_mosaic(
        [
            ['lleft', 'left', 'main'],
            ['purity_cb', 'purity_cb', 'purity'],
            # ['.', '.', 'purity_cb'],
            ['ploidy_cb', 'ploidy_cb', 'ploidy'],
            # ['.', '.', 'ploidy_cb']
        ],
        width_ratios=(1, 1, len_x),
        height_ratios=(len_y, 1, 1),
    )

    setup_main(axd['main'], samples, variant_names)
    setup_bottom(axd['purity'], samples, 'purity')
    setup_bottom(axd['ploidy'], samples, 'ploidy')
    setup_left(axd['left'], variant_names, 'left1')
    setup_left(axd['lleft'], variant_names, 'left2')

    write_bottom(fig, axd['purity'], axd['purity_cb'], [purities[x] for x in samples], 0, 1)
    write_bottom(fig, axd['ploidy'], axd['ploidy_cb'], [ploidies[x] for x in samples], 1, 8)

    fill_main(axd['main'], wgs_driver_vplists, panel_driver_vplists, SVinfo, samples, variants)



    

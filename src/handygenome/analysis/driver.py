import pyranges as pr

import handygenome.annotation.data as annotdata
import handygenome.annotation.annotitem as libannotitem
import handygenome.analysis.feature as anal_feature


def gene_cosmic_occurrences(gene_names, refver):
    result = dict()
    cosmic_vcf = annotdata.VCFS_COSMIC[refver]
    # print(f'getting coordinates for each gene')
    gene_coords = anal_feature.get_gene_coords(gene_names, refver)
    for gene_name, interval in gene_coords.items():
        # print(f'counting occurrence for gene {gene_name}')
        occurrence = 0
        for vr in cosmic_vcf.fetch(interval.chrom, interval.start0, interval.end0):
            cosmic_info = libannotitem._decode(vr.info['cosmic_info'][0])
            occurrence += sum(cosmic_info['occurrence'].values())
        result[gene_name] = occurrence
    
    return result

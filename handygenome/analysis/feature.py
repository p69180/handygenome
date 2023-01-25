import pyranges as pr

from handygenome.common import Interval
import handygenome.annotation.ensembl_rest as ensembl_rest
import handygenome.annotation.data as annotdata


def get_gene_coords(gene_names, refver):
    result = dict()
    rest_result = ensembl_rest.lookup_symbol_post(gene_names, refver, expand=False)
    for x in rest_result:
        gene_name = x['display_name']
        chrom = x['seq_region_name']
        start1 = x['start']
        end1 = x['end']
        result[gene_name] = Interval(chrom=chrom, start1=start1, end1=end1)

    return result


def get_gene_exon_coords(gene_names, refver):
    result = dict()
    rest_result = ensembl_rest.lookup_symbol_post(gene_names, refver, expand=True)
    for x in rest_result:
        gene_name = x['display_name']
        chrom = x['seq_region_name']
        start0 = x['start'] - 1
        end0 = x['end']

        result[gene_name] = dict()
        result[gene_name]['chrom'] = chrom
        result[gene_name]['gene'] = range(start0, end0)

        transcript_info = parse_transcript_info_from_lookup_symbol(x)
        for key, val in transcript_info.items():
            if val['is_canonical']:
                result[gene_name]['exons'] = val['exons'].copy()

    return result


def parse_transcript_info_from_lookup_symbol(lookup_result):
    result = dict()
    for x in lookup_result['Transcript']:
        transcript_id = x['id']
        result[transcript_id] = dict()

        result[transcript_id]['is_canonical'] = (x['is_canonical'] == 1)
        result[transcript_id]['display_name'] = x['display_name']
        result[transcript_id]['exons'] = list()
        for y in x['Exon']:
            start0 = y['start'] - 1
            end0 = y['end']
            result[transcript_id]['exons'].append(range(start0, end0))

    return result


def get_geneset_gr(refver):
    return pr.read_gff3(annotdata.GENESET_PATH[refver])



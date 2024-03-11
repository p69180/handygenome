import functools

import numpy as np
import pyranges as pr
import pandas as pd
import Bio.Seq

import handygenome.deco as deco
import handygenome.annotation.ensembl_rest as ensembl_rest
import handygenome.annotation.data as annotdata
import handygenome.refgenome.refgenome as refgenome
from handygenome.genomedf.genomedf import GenomeDataFrame


class TranscriptData:
    def __repr__(self):
        repr_string = ', '.join(
            f'{key}={getattr(self, key)}' 
            for key in (
                'gene_name', 'enst',
            )
        )
        return f'<{self.__class__.__name__} ({repr_string})>'

    @classmethod
    def from_ensembl_rest_result(cls, rest_result, refver):
        result = cls()
        result.refver = refver

        # non-transcript information
        result.chrom = rest_result['seq_region_name']
        result.gene_name = rest_result['display_name']
        result.ensg = rest_result['id']
        assert rest_result['strand'] in [1, -1]
        result.is_forward = (rest_result['strand'] == 1)

        # pick canonical transcript subresult
        result.canonical_enst = rest_result['canonical_transcript'].split('.')[0]
        canon_subresult = [x for x in rest_result['Transcript'] if x['id'] == result.canonical_enst]
        assert len(canon_subresult) == 1
        canon_subresult = canon_subresult[0]

        # add canon tr information
        result.enst = canon_subresult['id']
        result.ensp = canon_subresult['Translation']['id']
        result.cds_start0 = canon_subresult['Translation']['start'] - 1
        result.cds_end0 = canon_subresult['Translation']['end']

        start0s = list()
        end0s =  list()
        ense = list()
        for x in canon_subresult['Exon']:
            ense.append(x['id'])
            start0s.append(x['start'] - 1)
            end0s.append(x['end'])
        result.exon_gdf = GenomeDataFrame.from_data(
            refver=refver,
            chroms=result.chrom,
            start0s=start0s,
            end0s=end0s,
            ense=ense,
        )

        result.exon_gdf.sort()
        exon_numbers = np.arange(result.exon_gdf.nrow) + 1
        if not result.is_forward:
            exon_numbers = exon_numbers[::-1]
        result.exon_gdf['number'] = exon_numbers

        return result

    #############
    # utilities #
    #############

    @property
    def fasta(self):
        return refgenome.get_fasta(self.refver)

    @functools.cached_property
    def intron_gdf(self):
        start0s = self.exon_gdf.end0s[:-1]
        end0s = self.exon_gdf.start0s[1:]
        numbers = np.arange(self.exon_gdf.nrow - 1) + 1
        if not self.is_forward:
            numbers = numbers[::-1]

        gdf = GenomeDataFrame.from_data(
            refver=self.refver,
            chroms=self.chrom,
            start0s=start0s,
            end0s=end0s,
            number=numbers,
        )
        return gdf

    @functools.cached_property
    def cds_margin_gdf(self):
        return GenomeDataFrame.from_data(
            refver=self.refver,
            chroms=self.chrom,
            start0s=self.cds_start0,
            end0s=self.cds_end0,
        )

    def get_cds_gdf(self, exons=(None, None)):
        assert isinstance(exons, tuple)
        assert len(exons) == 2
        if (exons[0] is not None) and (exons[1] is not None):
            assert exons[1] >= exons[0]

        all_cds_gdf = self.exon_gdf.intersect(self.cds_margin_gdf)
        exons_start = (all_cds_gdf['number'].min() if (exons[0] is None) else exons[0])
        exons_end = (all_cds_gdf['number'].max() if (exons[1] is None) else exons[1])

        ###

        selector = np.logical_and(
            (all_cds_gdf['number'] >= exons_start),
            (all_cds_gdf['number'] <= exons_end),
        )
        return all_cds_gdf.loc[selector, :]

    def get_cds(self, as_ref=False, exons=(None, None)):
        cds_gdf = self.get_cds_gdf(exons=exons)
        seq_list = [
            self.fasta.fetch(chrom, start0, end0) 
            for (chrom, start0, end0) in zip(cds_gdf.chroms, cds_gdf.start0s, cds_gdf.end0s)
        ]
        cds_as_ref = ''.join(seq_list)
        if as_ref:
            return cds_as_ref
        else:
            if self.is_forward:
                return cds_as_ref
            else:
                return Bio.Seq.reverse_complement(cds_as_ref)

    def get_protein(self):
        cds = self.get_cds(as_ref=False)
        return Bio.Seq.translate(cds, to_stop=False, cds=True)


@deco.get_deco_atleast1d(['gene_names'])
def get_gene_coords(gene_names, refver):
    result = dict()
    rest_result = ensembl_rest.lookup_symbol_post(gene_names, refver, expand=False)
    for x in rest_result:
        gene_name = x['display_name']
        chrom = x['seq_region_name']
        start1 = x['start']
        end1 = x['end']
        result[gene_name] = {
            'chrom': chrom,
            'start0': start1 - 1,
            'end0': end1,
            'gene': gene_name,
        }

    return result


@deco.get_deco_atleast1d(['gene_names'])
def get_gene_coords_gdf(gene_names, refver):
    """Genes not searchable from "Ensembl REST lookup symbol" are discarded
    """
    data = get_gene_coords(gene_names, refver)
    chroms = list()
    starts = list()
    ends = list()
    result_gene_names = list()

    for gene in gene_names:
        if gene in data.keys():
            subdic = data[gene]
            chroms.append(subdic['chrom'])
            starts.append(subdic['start0'])
            ends.append(subdic['end0'])
            result_gene_names.append(gene)
        else:
            chroms.append(pd.NA)
            starts.append(pd.NA)
            ends.append(pd.NA)
            result_gene_names.append(pd.NA)

    return GenomeDataFrame.from_data(
        refver=refver,
        chroms=chroms,
        start0s=starts,
        end0s=ends,
        gene=result_gene_names,
    )


@functools.cache
def get_all_gene_coord_df_fromlocal(refver):
    geneset_gr = get_geneset_gr(refver)

    row_selector = geneset_gr.ID.notna().to_numpy()
    subgr = geneset_gr[row_selector]

    selector2 = np.char.startswith(np.asarray(subgr.ID, dtype=str), 'gene:')
    subgr2 = subgr[selector2]

    result_df = pd.DataFrame.from_dict(
        {'Chromosome': subgr2.Chromosome, 'Start': subgr2.Start, 'End': subgr2.End, 'gene': subgr2.Name}
    )
    result_df.drop_duplicates('gene', inplace=True)

    return result_df


@deco.get_deco_atleast1d(['gene_names'])
def get_gene_coords_gdf_fromlocal(gene_names, refver):
    all_gene_coord_df = get_all_gene_coord_df_fromlocal(refver)
    left_df = pd.DataFrame.from_dict({'gene': gene_names})
    joined = left_df.join(all_gene_coord_df.set_index('gene'), on='gene', how='left')
    result = GenomeDataFrame.from_frame(frame=joined, refver=refver)
    return result


def get_ensid_name_map(self, refver):
    geneset_gr = self.get_geneset_gr(refver)
    selector = geneset_gr.gene_id.notna()
    all_ensid = geneset_gr.gene_id.loc[selector]
    all_genename = geneset_gr.Name.loc[selector]
    genename_map = dict(zip(all_ensid, all_genename))
    return genename_map


######


@deco.get_deco_atleast1d(['gene_names'])
def get_gene_exon_coords(gene_names, refver):
    result = dict()
    rest_result = ensembl_rest.lookup_symbol_post(gene_names, refver, expand=True)
    for x in rest_result:
        gene_name = x['display_name']
        result[gene_name] = dict()
        result[gene_name]['chrom'] = x['seq_region_name']
        result[gene_name]['gene'] = range(x['start'] - 1, x['end'])
        if x['strand'] not in (1, -1):
            raise Exception(f'"strand" value from REST lookup result is neither 1 nor -1.')
        result[gene_name]['is_forward'] = (x['strand'] == 1)

        transcript_info = parse_transcript_info_from_lookup_symbol(x)
        for key, val in transcript_info.items():
            if val['is_canonical']:
                result[gene_name]['exon'] = val['exon'].copy()

    return result


def parse_transcript_info_from_lookup_symbol(lookup_result):
    result = dict()
    for x in lookup_result['Transcript']:
        transcript_id = x['id']
        result[transcript_id] = dict()

        result[transcript_id]['is_canonical'] = (x['is_canonical'] == 1)
        result[transcript_id]['display_name'] = x['display_name']
        result[transcript_id]['exon'] = list()
        for y in x['Exon']:
            start0 = y['start'] - 1
            end0 = y['end']
            result[transcript_id]['exon'].append(range(start0, end0))

    return result


@refgenome.deco_refseq_refver
@functools.cache
def get_geneset_gr(refver):
    return pr.read_gff3(annotdata.GENESET_PATH[refver])


def get_gene_exon_coords_fromgff(gene_names, geneset_gr, refver=None, rest_result=None):
    """Returns CDS and UTRs as well as exons"""
    if rest_result is None:
        rest_result = ensembl_rest.lookup_symbol_post(gene_names, refver, expand=True)

    result = dict()
    for gene, rest_data in zip(gene_names, rest_result):
        subresult = dict()
        # chrom, gene, is_forward
        subresult = dict()
        subresult['chrom'] = rest_data['seq_region_name']
        subresult['gene'] = range(rest_data['start'] - 1, rest_data['end'])
        if rest_data['strand'] not in (1, -1):
            raise Exception(f'"strand" value from REST lookup result is neither 1 nor -1.')
        subresult['is_forward'] = (rest_data['strand'] == 1)
        subresult['canonical_transcript'] = rest_data['canonical_transcript'].split('.')[0]

        # exon, CDS, UTR
        subgr = geneset_gr[subresult['chrom'], subresult['gene'].start:subresult['gene'].stop]
        subgr = subgr[
            subgr.Parent == f'transcript:{subresult["canonical_transcript"]}'
        ]
        sources = subgr.Source.unique()
        if len(sources) > 1:
            source = max(sources, key=(lambda x: subgr[subgr.Source == x].df.shape[0]))
        else:
            source = sources[0]

        subgr = subgr[subgr.Source == source]
        for key, subdf in subgr.df.groupby('Feature'):
            subresult[key] = list(
                subdf.loc[:, ['Start', 'End']].apply(lambda x: range(*x), axis=1)
            )

        result[gene] = subresult

    return result



        



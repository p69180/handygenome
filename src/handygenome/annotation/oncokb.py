import itertools
import functools

import pandas as pd
import numpy as np

import handygenome.tools as tools
import handygenome.deco as deco
import handygenome.network as network
import handygenome.annotation.annotitem as annotitem
import handygenome.analysis.feature as libfeature
import handygenome.refgenome.refgenome as refgenome
from handygenome.variant.vcfspec import Vcfspec
from handygenome.genomedf.genomedf import GenomeDataFrame


URL_allCuratedGenes = 'https://www.oncokb.org/api/v1/utils/allCuratedGenes'
URL_cancerGeneList = 'https://www.oncokb.org/api/v1/utils/cancerGeneList'
URL_HGVSG = 'https://www.oncokb.org/api/v1/annotate/mutations/byHGVSg'
URL_byProteinChange = 'https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange'
URL_CNV = 'https://www.oncokb.org/api/v1/annotate/copyNumberAlterations'
URL_SV = 'https://www.oncokb.org/api/v1/annotate/structuralVariants'


class OncoKBInfo(annotitem.AnnotItemInfoSingle):
    def check_oncogenic(self):
        return self['oncogenic'] in ('Oncogenic', 'Likely Oncogenic')
        
    def check_targetable(self):
        return self['highestSensitiveLevel'] is not None

    def check_implication(self):
        return any(
            (self[key] is not None)
            for key in (
                'highestSensitiveLevel', 
                'highestResistanceLevel', 
                'highestDiagnosticImplicationLevel', 
                'highestPrognosticImplicationLevel',
            )
        )
    def check_meaningful(self):
        return self.check_implication() or self.check_oncogenic()


class OncoKBInfoALTlist(annotitem.AnnotItemInfoALTlist):
    meta = {
        "ID": "oncokb",
        "Number": "A",
        "Type": "String",
        "Description": "OncoKB data annotation, acquired throught web API(https://www.oncokb.org/swagger-ui/index.html).",
    }
    unit_class = OncoKBInfo

    @classmethod
    def from_vr(cls, vr):
        return cls.from_vr_base(vr)

    def write(self, vr):
        self.write_base(vr)


@functools.cache
def get_allCuratedGenes_cache(token, version=None):
    params= {'includeEvidence': True}
    if version is not None:
        params['version'] = version

    df = pd.DataFrame.from_dict(
        network.http_get(
            url=URL_allCuratedGenes,
            headers={
                'Authorization': f'Bearer {token}',
                'accept': 'application/json',
            },
            params=params,
        )
    )
    return df


@functools.cache
def get_allCuratedGenes(token, refver, version=None, drop_nocoord=True):
    refseq_refver = refgenome.get_refseq_refver(refver)
    assert refseq_refver in ['GRCh37', 'GRCh38']

    df = get_allCuratedGenes_cache(token, version=version)
    gene_names = df['hugoSymbol']
    gene_coord_gdf = libfeature.get_gene_coords_gdf_fromlocal(gene_names, refseq_refver)

    result_df = df.copy()
    result_df['Chromosome'] = gene_coord_gdf.df['Chromosome']  # may contain NA
    result_df['Start'] = gene_coord_gdf.df['Start']
    result_df['End'] = gene_coord_gdf.df['End']
    result = GenomeDataFrame.from_frame(frame=result_df, refver=refver)

    if drop_nocoord:
        result = result.loc[pd.notna(result.chroms), :]

    return result


@functools.cache
def get_cancerGeneList_cache(token, version=None):
    if version is None:
        params = None
    else:
        params = {'version': version}
    return pd.DataFrame.from_dict(
        network.http_get(
            url=URL_cancerGeneList,
            headers={
                'Authorization': f'Bearer {token}',
                'accept': 'application/json',
            },
            params=params,
        )
    )


@functools.cache
def get_cancerGeneList(token, refver, version=None, drop_nocoord=True):
    refseq_refver = refgenome.get_refseq_refver(refver)
    assert refseq_refver in ['GRCh37', 'GRCh38']

    df = get_cancerGeneList_cache(token, version=version)
    gene_names = df['hugoSymbol']
    gene_coord_gdf = libfeature.get_gene_coords_gdf_fromlocal(gene_names, refseq_refver)

    result_df = df.copy()
    result_df['Chromosome'] = gene_coord_gdf.df['Chromosome']
    result_df['Start'] = gene_coord_gdf.df['Start']
    result_df['End'] = gene_coord_gdf.df['End']
    result = GenomeDataFrame.from_frame(frame=result_df, refver=refver)

    if drop_nocoord:
        result = result.loc[pd.notna(result.chroms), :]
    return result


def query_proteinchange_get(hugosymbol, alteration, token, GRCh38=False, tumor_type=None, evidence_types=None):
    params={
        'hugoSymbol': hugosymbol,
        'alteration': alteration,
    }
    if tumor_type is not None:
        params['tumorType'] = tumor_type
    if evidence_types is not None:
        params['evidenceType'] = ','.join(evidence_types)

    result = OncoKBInfo.init_nonmissing()
    result.update_dict(
        network.http_get(
            url=URL_byProteinChange,
            params=params,
            headers={
                'Authorization': f'Bearer {token}',
                'accept': 'application/json',
            },
        )
    )
    return result


def query_hgvsg_get(hgvsg, token, tumor_type=None, evidence_types=None):
    """Args:
        hgvsg: e.g. '7:g.140453136A>T'
    """
    params={
        'hgvsg': hgvsg,
    }
    if tumor_type is not None:
        params['tumorType'] = tumor_type
    if evidence_types is not None:
        params['evidenceType'] = ','.join(evidence_types)

    result = OncoKBInfo.init_nonmissing()
    result.update_dict(
        network.http_get(
            url=URL_HGVSG,
            params=params,
            headers={
                'Authorization': f'Bearer {token}',
                'accept': 'application/json',
            },
        )
    )
    return result


def query_hgvsg_post(
    hgvsg_list, 
    token, 
    tumor_type=None, 
    evidence_types=None,
    GRCh38=False, 
    chunk_size=50,
):
    result = list()
    #NR = 0
    for hgvsg_chunk in tools.chunk_iter(hgvsg_list, chunk_size):
        #NR += 1
        #print(f'{NR * chunk_size} entries being processed')
        data = list()
        for hgvsg in hgvsg_chunk:
            dic = {'hgvsg': hgvsg}
            if tumor_type is not None:
                dic['tumorType'] = tumor_type
            if evidence_types is not None:
                dic['evidenceTypes'] = evidence_types
            if GRCh38:
                dic['referenceGenome'] = 'GRCh38'
            data.append(dic)

        for item in network.http_post(
            url=URL_HGVSG,
            data=data,
            headers={
                'Authorization': f'Bearer {token}',
                'accept': 'application/json',
                'Content-Type': 'application/json',
            },
        ):
            unit_oncokb = OncoKBInfo.init_nonmissing()
            unit_oncokb.update_dict(item)
            result.append(unit_oncokb)

    return result


@deco.get_deco_atleast1d(('gene_names', 'cnv_types', 'tumor_types'), keep_none=False)
def query_cnv_post(gene_names, cnv_types, token, tumor_types=None, GRCh38=False, chunk_size=50):
    """Args:
        gene_names: Hugo symbols
    """
    assert all(x in ('AMPLIFICATION', 'DELETION', 'GAIN', 'LOSS') for x in cnv_types)
    gene_names, cnv_types, tumor_types = np.broadcast_arrays(gene_names, cnv_types, tumor_types)

    indexes = iter(range(len(gene_names)))

    result = list()
    for idx_chunk in tools.chunk_iter(indexes, chunk_size):
        data = list()
        for idx in idx_chunk:
            dic = {
                'copyNameAlterationType': cnv_types[idx],
                "gene": {
                    #"entrezGeneId": 0,
                    "hugoSymbol": gene_names[idx],
                },
            }
            if tumor_types is not None:
                dic['tumorType'] = tumor_types[idx]
            if GRCh38:
                dic['referenceGenome'] = 'GRCh38'
            data.append(dic)

        for item in network.http_post(
            url=URL_CNV,
            data=data,
            headers={
                'Authorization': f'Bearer {token}',
                'accept': 'application/json',
                'Content-Type': 'application/json',
            },
        ):
            unit_oncokb = OncoKBInfo.init_nonmissing()
            unit_oncokb.update_dict(item)
            result.append(unit_oncokb)

    return result


@deco.get_deco_atleast1d(
    ('gene1_names', 'gene2_names', 'is_functional', 'sv_types', 'tumor_types'), 
    keep_none=False,
)
def query_sv_post(gene1_names, gene2_names, is_functional, token, sv_types='UNKNOWN', tumor_types=None, GRCh38=False, chunk_size=50):
    """Args:
        'None' in 'gene1_names' or 'gene2_names' means that gene is absent
    
    When 'is_functional' is True:
        - It is interpretated as fusion 
        - regardless of values of 'sv_types' value
        - even when one of genes are 'None'
    When 'is_functional' is False:
        - It is interpretated as truncation event
        - Only consequence for the 'gene1' is returned.
            e.g. gene1='ARID1A' & gene2='TP53' & is_functional=False : 
                Only consequence of ARID1A truncation is returned
    """

    assert all(
        x in ('DELETION', 'TRANSLOCATION', 'DUPLICATION', 'INSERTION', 'INVERSION', 'FUSION', 'UNKNOWN')
        for x in sv_types
    )
    gene1_names, gene2_names, is_functional, sv_types, tumor_types = np.broadcast_arrays(
        gene1_names, gene2_names, is_functional, sv_types, tumor_types
    )
    is_functional = [bool(x) for x in is_functional]

    indexes = iter(range(len(gene1_names)))

    result = list()
    for idx_chunk in tools.chunk_iter(indexes, chunk_size):
        data = list()
        for idx in idx_chunk:
            dic = {
                "functionalFusion": is_functional[idx],
                "geneA": {
                    "hugoSymbol": gene1_names[idx],
                },
                "geneB": {
                    "hugoSymbol": gene2_names[idx],
                },
                "structuralVariantType": sv_types[idx],
            }
            if tumor_types is not None:
                dic['tumorType'] = tumor_types[idx]
            if GRCh38:
                dic['referenceGenome'] = 'GRCh38'
            data.append(dic)

        for item in network.http_post(
            url=URL_SV,
            data=data,
            headers={
                'Authorization': f'Bearer {token}',
                'accept': 'application/json',
                'Content-Type': 'application/json',
            },
        ):
            unit_oncokb = OncoKBInfo.init_nonmissing()
            unit_oncokb.update_dict(item)
            result.append(unit_oncokb)

    return result


def add_oncokb_info(vp_list, token, **kwargs):
    vp_list = list(vp_list)
    hgvsg_list = list()
    vp_indexes = list()
    for vp_idx, vp in enumerate(vp_list):
        for sub_vcfspec in vp.vcfspec.iter_monoalts():
            hgvsg_list.append(sub_vcfspec.to_hgvsg())
            vp_indexes.append(vp_idx)

    query_result = query_hgvsg_post(hgvsg_list=hgvsg_list, token=token, **kwargs)
    for key, subiter in itertools.groupby(
        zip(vp_indexes, query_result),
        key=(lambda x: x[0]),
    ):
        oncokb_ALTlist = OncoKBInfoALTlist()
        for vp_idx, oncokb_item in subiter:
            oncokb_ALTlist.append(oncokb_item)
        vp_list[vp_idx].oncokb = oncokb_ALTlist


def make_OncoKBInfoALTlist_list(vr_iterator, token, **kwargs):
    hgvsg_list = list()
    vr_indexes = list()
    for vr_idx, vr in enumerate(vr_iterator):
        vcfspec = Vcfspec.from_vr(vr)
        for sub_vcfspec in vcfspec.iter_annotation_forms():
            hgvsg_list.append(sub_vcfspec.to_hgvsg())
            vr_indexes.append(vr_idx)

    result = list()
    query_result = query_hgvsg_post(hgvsg_list=hgvsg_list, token=token, **kwargs)
    for key, subiter in itertools.groupby(
        zip(vr_indexes, query_result),
        key=(lambda x: x[0]),
    ):
        oncokb_ALTlist = OncoKBInfoALTlist()
        for vr_idx, oncokb_item in subiter:
            oncokb_ALTlist.append(oncokb_item)
        result.append(oncokb_ALTlist)

    return result


#######
# CNV #
#######

def make_cancergene_cnv_annotation_data(token, refver=None, cancer_gene_list=None):
    if cancer_gene_list is None:
        oncokb_cancer_genes = get_cancerGeneList(token, refver=refver)
        oncokb_cancer_genes = oncokb_cancer_genes.loc[oncokb_cancer_genes.df['Chromosome'].notna(), :]
        cancer_gene_list = oncokb_cancer_genes['hugoSymbol']

    cancergene_del_oncokb_annots = query_cnv_post(cancer_gene_list, cnv_types='DELETION', token=token)
    cancergene_amp_oncokb_annots = query_cnv_post(cancer_gene_list, cnv_types='AMPLIFICATION', token=token)
    result = dict()
    for gene, delinfo, ampinfo in zip(cancer_gene_list, cancergene_del_oncokb_annots, cancergene_amp_oncokb_annots):
        result[gene] = {'del': delinfo, 'amp': ampinfo}

    return result


def make_cnv_meaningful_cancergene_gdf(refver, cancergene_cnv_annotation_data=None, token=None):
    if cancergene_cnv_annotation_data is None:
        cancergene_cnv_annotation_data = make_cancergene_cnv_annotation_data(token=token, refver=refver, cancer_gene_list=None)

    meaningful_cancergenes = {'amp': list(), 'del': list(), 'amp_target': list(), 'del_target': list()}
    for gene, subdic in cancergene_cnv_annotation_data.items():
        if subdic['amp'].check_oncogenic(): 
            meaningful_cancergenes['amp'].append(gene)
        if subdic['amp'].check_targetable():
            meaningful_cancergenes['amp_target'].append(gene)
        if subdic['del'].check_oncogenic():
            meaningful_cancergenes['del'].append(gene)
        if subdic['del'].check_targetable():
            meaningful_cancergenes['del_target'].append(gene)

    meaningful_cancergenes['any'] = list(set(meaningful_cancergenes['amp']).union(meaningful_cancergenes['del']))

    meaningful_cancergene_gdf = libfeature.get_gene_coords_gdf_fromlocal(meaningful_cancergenes['any'], refver)
    meaningful_cancergene_gdf['amp_oncogenic']  = [(x in meaningful_cancergenes['amp']) for x in meaningful_cancergene_gdf['gene']]
    meaningful_cancergene_gdf['del_oncogenic'] = [(x in meaningful_cancergenes['del']) for x in meaningful_cancergene_gdf['gene']]
    meaningful_cancergene_gdf['amp_target']  = [(x in meaningful_cancergenes['amp_target']) for x in meaningful_cancergene_gdf['gene']]
    meaningful_cancergene_gdf['del_target'] = [(x in meaningful_cancergenes['del_target']) for x in meaningful_cancergene_gdf['gene']]
    meaningful_cancergene_gdf = meaningful_cancergene_gdf.loc[pd.notna(meaningful_cancergene_gdf.chroms), :]

    return meaningful_cancergene_gdf, cancergene_cnv_annotation_data



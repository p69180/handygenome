import os
import pickle
import itertools

import numpy as np
import pandas as pd
import scipy.stats
import scipy.optimize
import scipy.interpolate
import sklearn.mixture

import handygenome
import handygenome.deco as deco
import handygenome.refgenome.refgenome as refgenome
import handygenome.variant.variantplus as libvp
from handygenome.variant.variantplus import VariantPlusList
from handygenome.genomedf import GenomeDataFrame


class BAFDataFrame(GenomeDataFrame):
    baf_colname = 'baf'

    def init_sanitycheck(self):
        assert self.__class__.baf_colname in self.annot_cols

    @classmethod
    def from_vcf(cls, vcf_path):
        pass

        #########


    def load_germline_vcf(
        self, 
        sampleid, 
        vcf_path, 
        vcf_sampleid_tumor,
        vcf_sampleid_normal=None,
        nproc=1, 
        logging_lineno=50000,
    ):
        self.data.setdefault(sampleid, dict())

        if vcf_sampleid_normal is None:
            vcf_sampleids = [vcf_sampleid_tumor]
        else:
            vcf_sampleids = [vcf_sampleid_tumor, vcf_sampleid_normal]

        # load vafdf
        vaf_df = variantplus.get_vafdf(
            vcf_path, 
            sampleid=vcf_sampleids, 
            nproc=nproc,
        )
        # rename vaf columns
        vaf_df.rename(
            columns={f'vaf_{vcf_sampleid_tumor}': 'vaf_raw_tumor'}, inplace=True,
        )
        if vcf_sampleid_normal is not None:
            vaf_df.rename(
                columns={f'vaf_{vcf_sampleid_normal}': 'vaf_raw_normal'}, inplace=True,
            )

        # postprocess
        #vaf_df = vaf_df.loc[vaf_df['vaf_raw'].notna().to_numpy(), :]
        vaf_df.reset_index(drop=True, inplace=True)
        vaf_df['baf_raw_tumor'] = cnvmisc.get_bafs(vaf_df['vaf_raw_tumor'])
        if 'vaf_raw_normal' in vaf_df:
            vaf_df['baf_raw_normal'] = cnvmisc.get_bafs(vaf_df['vaf_raw_normal'])

        selected_cols = [
            'Chromosome', 
            'Start', 
            'End', 
            'vaf_raw_tumor', 
            'baf_raw_tumor', 
        ]
        if vcf_sampleid_normal is not None:
            selected_cols.append('vaf_raw_normal')
            selected_cols.append('baf_raw_normal')
        vaf_df = vaf_df.loc[:, selected_cols]
        self.data[sampleid]['original_baf'] = vaf_df


def get_baf_from_vaf(vafs):
    """Args:
        vafs: 1st dimension means genomic positions, 2nd dimension means allelic vafs
    """
    # preprocessing "vafs" into a 2d array
    vafs = np.atleast_1d(vafs)
    if vafs.ndim == 1:
        vafs = np.stack([vafs, 1 - vafs], axis=1)
    elif vafs.ndim == 2:
        if vafs.shape[1] == 1:
            vafs = np.stack([vafs[:, 0], 1 - vafs[:, 0]], axis=1)
    else:
        raise Exception(f'Input array must not have ndim larger than 2.')
        
    # main
    vafs = np.sort(vafs, axis=1)
    bafs = vafs[:, :-1]
    return bafs


def get_bafdf_from_vafdf(vafdf, refver):
    sampleids = vafdf.columns.get_level_values('sampleid').unique()
    sampleids = sampleids[sampleids.notna()]

    result = dict()
    for sid in sampleids:
        vafs = vafdf.loc[:, pd.IndexSlice[sid, :]].to_numpy()
        bafs = get_baf_from_vaf(vafs)
        src_data = dict()
        src_data['refver'] = refver
        src_data['chroms'] = np.squeeze(vafdf.loc[:, (slice(None), 'Chromosome')].to_numpy())
        src_data['starts'] = np.squeeze(vafdf.loc[:, (slice(None), 'Start')].to_numpy())
        src_data['ends'] = np.squeeze(vafdf.loc[:, (slice(None), 'End')].to_numpy())
        for idx in range(bafs.shape[1]):
            baf_colname = f'baf{idx}'
            src_data[baf_colname] = bafs[:, idx]

        gdf = BAFDataFrame.from_data(**src_data)
        result[sid] = gdf

    return result


@deco.get_deco_atleast1d(['sampleids'])
def get_bafdf_from_vcf(
    vcf_path, 
    sampleids, 
    n_allele=2,
    nproc=1,
    exclude_other=False,
    prop=None,
    vpfilter=None,
    verbose=True,
):
    refver = refgenome.infer_refver_vcfpath(vcf_path)
    vafdf = libvp.get_vafdf(
        vcf_path, 
        sampleids, 
        n_allele=n_allele,
        nproc=nproc,
        exclude_other=exclude_other,
        prop=prop,
        vpfilter=vpfilter,
        verbose=verbose,
    )
    return get_bafdf_from_vafdf(vafdf, refver)


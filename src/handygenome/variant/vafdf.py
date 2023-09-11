import itertools
import multiprocessing

import numpy as np
import pandas as pd

import handygenome.deco as deco
import handygenome.refgenome.refgenome as refgenome
import handygenome.logutils as logutils
from handygenome.variant.variantplus import VariantPlusList
import handygenome.vcfeditor.misc as vcfmisc


def get_vafdf(
    vcf_path, 
    sampleids, 
    n_allele=2,
    nproc=1,
    exclude_other=False,
    prop=None,
    vpfilter=None,
    verbose=True,
):
    # setup params
    refver = refgenome.infer_refver_vcfpath(vcf_path)

    if verbose:
        logutils.log(f'Extracting vcf position information') 

    all_position_info = vcfmisc.get_vcf_positions(
        vcf_path, as_iter=False, verbose=False,
    )
    split_position_info = [
        x for x in np.array_split(all_position_info, nproc) 
        if x.shape[0] != 0
    ]

    # run multiprocess jobs
    if verbose:
        logutils.log(f'Running parallel jobs') 

    args = (
        (
            position_info, 
            refver, 
            vcf_path, 
            sampleids, 
            n_allele, 
            exclude_other,
            prop,
            vpfilter,
        )
        for position_info in split_position_info
    )
    with multiprocessing.Pool(nproc) as pool:
        mp_result = pool.starmap(get_vafdf_targetfunc, args)
        if verbose:
            logutils.log(f'Concatenating split job dataframes') 
        result = pd.concat(itertools.chain.from_iterable(mp_result), axis=0)

    result.reset_index(drop=True, inplace=True)
    return result


def get_vafdf_targetfunc(
    position_info, 
    refver, 
    vcf_path, 
    sampleids, 
    n_allele, 
    exclude_other,
    prop,
    vpfilter,
):
    fetchargs_list = vcfmisc.get_fetchargs_from_vcf_positions(position_info, refver)
    dflist = list()
    for fetchargs in fetchargs_list:
        df = get_vafdf_nonparallel(
            vcf_path, 
            sampleids, 
            chrom=fetchargs[0], start0=fetchargs[1], end0=fetchargs[2],
            n_allele=n_allele,

            exclude_other=exclude_other,
            prop=prop,
            vpfilter=vpfilter,
        )
        dflist.append(df)

    return dflist


@deco.get_deco_atleast1d(['sampleids'])
def get_vafdf_nonparallel(
    vcf_path, 
    sampleids, 
    chrom=None, start0=None, end0=None,
    n_allele=2, 
    exclude_other=False,
    prop=None,
    vpfilter=None,
):
    vplist = VariantPlusList.from_vcf_lazy(
        vcf_path, 
        logging_lineno=None, 
        verbose=False,
        init_all_attrs=False,
        vp_init_params=dict(
            init_readstats=True,
            sampleid_list=sampleids,
        ),
    )
    vafdf = vplist.get_vafdf(
        sampleids=sampleids,
        chrom=chrom, start0=start0, end0=end0,
        n_allele=n_allele,
        exclude_other=exclude_other,
        lazy=True,
        prop=prop,
        vpfilter=vpfilter,
    )
    return vafdf



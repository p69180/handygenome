import os
import collections
import functools

import handygenome.refgenome.refgenome as refgenome
import handygenome.refgenome.refverfile as refverfile
import handygenome.genomedf as genomedf
from handygenome.genomedf import GenomeDataFrame as GDF


MALE_HAPLOID_CHROMS = ('X', 'Y', 'chrX', 'chrY')
DEFAULT_CNT_WEIGHT = 5


class CCFInfo(
    collections.namedtuple('CCFInfo', ('CNm', 'ccf', 'mutated_allele'))
):
    pass


class CPPair(
    collections.namedtuple('CPPair', ('cellularity', 'ploidy'))
):
    pass


################
# miscellanous #
################

def check_haploid(is_female, chrom):
    return (not is_female) and (chrom in MALE_HAPLOID_CHROMS)


#@functools.cache
def get_default_CNn_diploid(refver, is_female):
    """Only accepts PAR-available reference versions"""
    chromdict = refgenome.get_chromdict(refver)
    autosomal_chroms = [
        x for x in chromdict.contigs 
        if (
            (refgenome.PAT_ASSEMBLED_CHROM.fullmatch(x) is not None)
            and (refgenome.normalize_chrom(x) not in MALE_HAPLOID_CHROMS)
        )
    ]  # 1, 2, ...,
    X_chrom, Y_chrom = chromdict.XY_names
    par_gdf = refverfile.get_par_gdf(refver)

    # autosomal
    all_chrom_gdf = GDF.all_regions(refver, assembled_only=True)
    autosomal_gdf = all_chrom_gdf.subset_chroms(autosomal_chroms)
    autosomal_gdf['CNn'] = 2

    # sex
    if is_female:
        sex_gdf = all_chrom_gdf.subset_chroms([X_chrom, Y_chrom])
        sex_gdf['CNn'] = 2
    else:
        all_sex_gdf = all_chrom_gdf.subset_chroms([X_chrom, Y_chrom])
        nonpar_sex_gdf = all_sex_gdf.subtract(par_gdf)
        nonpar_sex_gdf['CNn'] = 1

        par_sex_gdf = par_gdf.copy()
        is_X = (par_sex_gdf['Chromosome'] == X_chrom)
        par_sex_gdf.loc[is_X, 'CNn'] = 2
        par_sex_gdf.loc[~is_X, 'CNn'] = 0

        sex_gdf = genomedf.concat([nonpar_sex_gdf, par_sex_gdf])

    # result
    result = genomedf.concat([autosomal_gdf, sex_gdf])
    result.sort()
    return result




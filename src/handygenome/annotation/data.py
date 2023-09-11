import os

import pysam

import handygenome
import handygenome.refgenome.refgenome as refgenome


def get_dbsnp_path(refver):
    standard = refgenome.REFVERINFO.standardize(refver)
    refveritem = refgenome.REFVERINFO.access_refveritem(standard)
    if refveritem['dbsnp_path'] is None:
        return os.path.join(handygenome.DIRS['data'], 'popfreq', f'dbSNP_{standard}.vcf.gz')
    else:
        return refveritem['dbsnp_path']


REPEATS_PATH = refgenome.RefverDict({
    'GRCh37': '/home/users/pjh/References/ucsc_data_files/modified_files/grch37/repeatmasker.bed.gz',
    'GRCh38': '/home/users/pjh/References/ucsc_data_files/modified_files/grch38/repeatmasker.bed.gz',
})

GENESET_PATH = refgenome.RefverDict({
    'GRCh37': '/home/users/pjh/References/ensembl_data_files/modified_files/grch37/geneset_gff3_sorted.gz',
    'GRCh38': '/home/users/pjh/References/ensembl_data_files/modified_files/grch38/geneset_gff3_sorted.gz',
})

REGULATORY_PATH = refgenome.RefverDict({
    'GRCh37': '/home/users/pjh/References/ensembl_data_files/modified_files/grch37/ensembl_regulatory_sorted.grch37.gff.gz',
    'GRCh38': '/home/users/pjh/References/ensembl_data_files/modified_files/grch38/ensembl_regulatory_sorted.grch38.gff.gz',
})

DBSNP_PATHS = refgenome.RefverDict({
    'GRCh37': f'{handygenome.DIRS["data"]}/popfreq/dbSNP_b155_GRCh37.p13.vcf.gz',
    'GRCh38': f'{handygenome.DIRS["data"]}/popfreq/dbSNP_b156_GRCh38.p13.vcf.gz',
})

COSMIC_PATHS = refgenome.RefverDict({
    'GRCh37': f'{handygenome.DIRS["data"]}/cosmic/v96/grch37/nonSV.vcf.gz',
    'GRCh38': None,
})

###

TABIXFILES_REPEATS = refgenome.RefverDict({
    'GRCh37': pysam.TabixFile(REPEATS_PATH['GRCh37'], parser=pysam.asBed()),
    'GRCh38': pysam.TabixFile(REPEATS_PATH['GRCh38'], parser=pysam.asBed()),
})

TABIXFILES_GENESET = refgenome.RefverDict({
    'GRCh37': pysam.TabixFile(GENESET_PATH['GRCh37'], parser=pysam.asGTF()),
    'GRCh38': pysam.TabixFile(GENESET_PATH['GRCh38'], parser=pysam.asGTF()),
})

TABIXFILES_REGULATORY = refgenome.RefverDict({
    'GRCh37': pysam.TabixFile(REGULATORY_PATH['GRCh37'], parser=pysam.asGTF()),
    'GRCh38': pysam.TabixFile(REGULATORY_PATH['GRCh38'], parser=pysam.asGTF()),
})

VCFS_DBSNP = refgenome.RefverDict({
    'GRCh37': pysam.VariantFile(DBSNP_PATHS['GRCh37']),
    'GRCh38': pysam.VariantFile(DBSNP_PATHS['GRCh38']),
})

VCFS_COSMIC = refgenome.RefverDict({
    'GRCh37': pysam.VariantFile(COSMIC_PATHS['GRCh37']),
    'GRCh38': None,
})

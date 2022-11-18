import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


REPEATS_PATH = common.RefverDict({
    'GRCh37': '/home/users/pjh/References/ucsc_data_files/modified_files/grch37/repeatmasker.bed.gz',
    'GRCh38': '/home/users/pjh/References/ucsc_data_files/modified_files/grch38/repeatmasker.bed.gz',
})

GENESET_PATH = common.RefverDict({
    'GRCh37': '/home/users/pjh/References/ensembl_data_files/modified_files/grch37/geneset_gff3_sorted.gz',
    'GRCh38': '/home/users/pjh/References/ensembl_data_files/modified_files/grch38/geneset_gff3_sorted.gz',
})

REGULATORY_PATH = common.RefverDict({
    'GRCh37': '/home/users/pjh/References/ensembl_data_files/modified_files/grch37/ensembl_regulatory_sorted.grch37.gff.gz',
    'GRCh38': '/home/users/pjh/References/ensembl_data_files/modified_files/grch38/ensembl_regulatory_sorted.grch38.gff.gz',
})

DBSNP_PATHS = common.RefverDict({
    'GRCh37': '/home/users/pjh/scripts/python_genome_packages/data/popfreq/dbSNP_b155_GRCh37.p13.vcf.gz',
    'GRCh38': '/home/users/pjh/References/dbSNP38/modified_files/dbSNP_b155_GRCh38.p13.vcf.gz',
})

#COSMIC_PATHS = {
#    'hg19': {
#        'coding': '/home/users/pjh/References/COSMIC/hg19/modified_files/v95/CosmicMutantExport.vcf.gz',
#        'noncoding': '/home/users/pjh/References/COSMIC/hg19/modified_files/v95/CosmicNCV.vcf.gz',
#        },
#    'GRCh38': {
#        'coding': '/home/users/pjh/References/COSMIC/GRCh38/modified_files/v95/CosmicMutantExport.vcf.gz',
#        'noncoding': '/home/users/pjh/References/COSMIC/GRCh38/modified_files/v95/CosmicNCV.vcf.gz',
#        },
#    }

COSMIC_PATHS = common.RefverDict({
    'GRCh37': '/home/users/pjh/scripts/python_genome_packages/data/cosmic/v96/grch37/nonSV.vcf.gz',
    'GRCh38': None,
})

###

TABIXFILES_REPEATS = common.RefverDict({
    'GRCh37': pysam.TabixFile(REPEATS_PATH['GRCh37'], parser=pysam.asBed()),
    'GRCh38': pysam.TabixFile(REPEATS_PATH['GRCh38'], parser=pysam.asBed()),
})

TABIXFILES_GENESET = common.RefverDict({
    'GRCh37': pysam.TabixFile(GENESET_PATH['GRCh37'], parser=pysam.asGTF()),
    'GRCh38': pysam.TabixFile(GENESET_PATH['GRCh38'], parser=pysam.asGTF()),
})

TABIXFILES_REGULATORY = common.RefverDict({
    'GRCh37': pysam.TabixFile(REGULATORY_PATH['GRCh37'], parser=pysam.asGTF()),
    'GRCh38': pysam.TabixFile(REGULATORY_PATH['GRCh38'], parser=pysam.asGTF()),
})

VCFS_DBSNP = common.RefverDict({
    'GRCh37': pysam.VariantFile(DBSNP_PATHS['GRCh37']),
    'GRCh38': pysam.VariantFile(DBSNP_PATHS['GRCh38']),
})

#VCFS_COSMIC = {
#    'GRCh37': {
#        'coding': pysam.VariantFile(COSMIC_PATHS['GRCh37']['coding']),
#        'noncoding': pysam.VariantFile(COSMIC_PATHS['GRCh37']['noncoding']),
#        },
#    'GRCh38': {
#        'coding': pysam.VariantFile(COSMIC_PATHS['GRCh38']['coding']),
#        'noncoding': pysam.VariantFile(COSMIC_PATHS['GRCh38']['noncoding']),
#        },
#    }
VCFS_COSMIC = common.RefverDict({
    'GRCh37': pysam.VariantFile(COSMIC_PATHS['GRCh37']),
    'GRCh38': None,
})

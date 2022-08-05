import pysam


REPEATS_PATH = {
    'hg19': '/home/users/pjh/References/ucsc_data_files/modified_files/grch37/repeatmasker.bed.gz',
    'hg38': '/home/users/pjh/References/ucsc_data_files/modified_files/grch38/repeatmasker.bed.gz',
    }

GENESET_PATH = {
    'hg19': '/home/users/pjh/References/ensembl_data_files/modified_files/grch37/geneset_gff3_sorted.gz',
    'hg38': '/home/users/pjh/References/ensembl_data_files/modified_files/grch38/geneset_gff3_sorted.gz',
    }

REGULATORY_PATH = {
    'hg19': '/home/users/pjh/References/ensembl_data_files/modified_files/grch37/ensembl_regulatory_sorted.grch37.gff.gz',
    'hg38': '/home/users/pjh/References/ensembl_data_files/modified_files/grch38/ensembl_regulatory_sorted.grch38.gff.gz',
    }

DBSNP_PATHS = {
    'hg19': '/home/users/pjh/References/dbSNP37/modified_files_220602/dbSNP_b155_GRCh37.p13.vcf.gz',
    'hg38': '/home/users/pjh/References/dbSNP38/modified_files/dbSNP_b155_GRCh38.p13.vcf.gz',
    }

COSMIC_PATHS = {
    'hg19': {
        'coding': '/home/users/pjh/References/COSMIC/hg19/modified_files/v95/CosmicMutantExport.vcf.gz',
        'noncoding': '/home/users/pjh/References/COSMIC/hg19/modified_files/v95/CosmicNCV.vcf.gz',
        },
    'hg38': {
        'coding': '/home/users/pjh/References/COSMIC/hg38/modified_files/v95/CosmicMutantExport.vcf.gz',
        'noncoding': '/home/users/pjh/References/COSMIC/hg38/modified_files/v95/CosmicNCV.vcf.gz',
        },
    }

###

TABIXFILES_REPEATS = {
    'hg19': pysam.TabixFile(REPEATS_PATH['hg19'], parser=pysam.asBed()),
    'hg38': pysam.TabixFile(REPEATS_PATH['hg38'], parser=pysam.asBed()),
    }

TABIXFILES_GENESET = {
    'hg19': pysam.TabixFile(GENESET_PATH['hg19'], parser=pysam.asGTF()),
    'hg38': pysam.TabixFile(GENESET_PATH['hg38'], parser=pysam.asGTF()),
    }

TABIXFILES_REGULATORY = {
    'hg19': pysam.TabixFile(REGULATORY_PATH['hg19'], parser=pysam.asGTF()),
    'hg38': pysam.TabixFile(REGULATORY_PATH['hg38'], parser=pysam.asGTF()),
    }

VCFS_DBSNP = {
    'hg19': pysam.VariantFile(DBSNP_PATHS['hg19']),
    'hg38': pysam.VariantFile(DBSNP_PATHS['hg38']),
    }

VCFS_COSMIC = {
    'hg19': {
        'coding': pysam.VariantFile(COSMIC_PATHS['hg19']['coding']),
        'noncoding': pysam.VariantFile(COSMIC_PATHS['hg19']['noncoding']),
        },
    'hg38': {
        'coding': pysam.VariantFile(COSMIC_PATHS['hg38']['coding']),
        'noncoding': pysam.VariantFile(COSMIC_PATHS['hg38']['noncoding']),
        },
    }

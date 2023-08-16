import pysam


def index_vcf(vcf_path, overwrite=True):
    pysam.tabix_index(vcf_path, preset='vcf', index=(vcf_path + '.csi'),
                      csi=True, force=overwrite)


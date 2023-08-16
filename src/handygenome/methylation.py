def add_meta(vcfheader):
    # ALT metainfo
    vcfheader.add_meta(
        key='ALT',
        items=[('ID', 'CPGMET'), ('Description', 'CpG site methylation')])


import textwrap

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
infoformat = importlib.import_module('.'.join([top_package_name, 'variantplus', 'infoformat']))


class SomaticIndex:
    meta = {'ID': 'SomaticIndex', 'Type': 'Integer', 'Number': 1,
            'Description': (
                'The index of the allele corresponding to somatic mutation. '
                '0 means REF, 1 means the first ALT allele, 2 means the '
                'second ALT allele, and so on. It is 1 in most cases.')}
    type = 'INFO'

    def __init__(self, index=1):
        self.index = index

    def __repr__(self):
        return f'<SomaticIndex={self.index}>'

    @classmethod
    def from_vr(cls, vr):
        val = infoformat.get_info(vr, cls.meta['ID'], collapse_tuple=False)
        return cls(val)

    @classmethod
    def write_meta(cls, vcfheader):
        vcfheader.add_meta(key=cls.type, items=cls.meta.items())
        
    def write(self, vr):
        infoformat.set_info(vr, self.__class__.meta['ID'], self.index,
                            typeconv=False)
        

class GermlineIndexes:
    meta = {'ID': 'GermlineIndexes', 'Type': 'Integer', 'Number': '.',
            'Description': (
                'The indexes of the germline alleles. '
                '0 means REF, 1 means the first ALT allele, 2 means the '
                'second ALT allele, and so on. It is 0,0 in most cases.')}
    type = 'INFO'

    def __init__(self, indexes=(0, 0)):
        self.indexes = tuple(indexes)
        
    def __repr__(self):
        return f'<GermlineIndexes={self.indexes}>'

    @classmethod
    def from_vr(cls, vr):
        val = infoformat.get_info(vr, cls.meta['ID'], collapse_tuple=False)
        return cls(val)
    
    @classmethod
    def write_meta(cls, vcfheader):
        vcfheader.add_meta(key=cls.type, items=cls.meta.items())

    def write(self, vr):
        infoformat.set_info(vr, self.__class__.meta['ID'], self.indexes,
                            typeconv=False)

    def write_GT(self, vr, sampleid):
        vr.samples[sampleid].allele_indices = self.indexes
        vr.samples[sampleid].phased = False

import re

import pysam
import numpy as np

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
variantplus = importlib.import_module('.'.join([top_package_name, 'variantplus', 'variantplus']))


class PanelOfNormal:
    @common.get_deco_num_set(('vcfspec', 'vp'), 1)
    def __init__(self, data, vp=None, vcfspec=None, refver=None):
        """Args:
            data: Dict(<sampleid>: <sample bam>)

        """
        self.data = data
        if vcfspec is not None:
            if refver is None:
                raise Exception(
                    f'When initializing from vcfspec, refver must be given.')
            vp = variantplus.VariantPlus(vcfspec=vcfspec, refver=refver)
            
        for sampleid, bam in self.data.items():
            if bam is None:
                raise Exception(
                    f'When initializing from vcfspec, sample bam '
                    f'objects must be given.')
            readstats = 

    def 
        

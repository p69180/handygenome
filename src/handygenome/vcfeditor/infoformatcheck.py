import os
import itertools
import warnings
import textwrap

import pysam

import handygenome.refgenome.refgenome as refgenome
import handygenome.workflow as workflow
import handygenome.variant.varianthandler as varianthandler
import handygenome.vcfeditor.headerhandler as headerhandler
import handygenome.vcfeditor.initvcf as initvcf
import handygenome.vcfeditor.misc as vcfmisc
import handygenome.variant.vcfspec as libvcfspec


def parse_header(vcfheader):
    header_data = {'info': dict(), 'format': dict()}
    for key, val in vcfheader.info.items():
        metadata = {
            'name': val.name,
            'number': val.number,
        }
        header_data['info'][metadata['name']] = metadata
    for key, val in vcfheader.formats.items():
        metadata = {
            'name': val.name,
            'number': val.number,
        }
        header_data['format'][metadata['name']] = metadata
    return header_data


def main(vcf_path):
    vcf = pysam.VariantFile(vcf_path)
    header_data = parse_header(vcf.header)




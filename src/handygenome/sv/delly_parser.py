import sys
import os
import re

import pysam

import handygenome.deco as deco
import handygenome.refgenome.refgenome as refgenome
import handygenome.variant.varianthandler as varianthandler
from handygenome.variant.vcfspec import Vcfspec
from handygenome.sv.breakends import Breakends
from handygenome.variant.variantplus import VariantPlus, VariantPlusList


def get_bnds_from_delly_vr(vr):
    CT = vr.info['CT'].split('to')
    is5prime1 = (CT[0] == '5')
    is5prime2 = (CT[1] == '5')

    chrom1 = vr.chrom
    pos1 = vr.pos
    chrom2 = vr.info['CHR2']
    pos2 = varianthandler.get_END(vr)

    return Breakends(
        (chrom1, pos1 - 1, is5prime1),
        (chrom2, pos2 - 1, is5prime2),
        inserted_seq='',
        refver=refgenome.infer_refver_vcfheader(vr.header),
    )
     

@deco.get_deco_atleast1d(['vcf_path_list'])
def load_delly_vcf(vcf_path_list, passonly=True):
    refver = refgenome.infer_refver_vcfpath(vcf_path_list[0])
    vplist = VariantPlusList(refver=refver)

    bnds_list = list()

    for vcfpath in vcf_path_list:
        with pysam.VariantFile(vcfpath) as in_vcf:
            for vr in in_vcf.fetch():
                if not varianthandler.check_filter_pass(vr):
                    continue

                bnds = get_bnds_from_delly_vr(vr)
                if any((bnds == x) for x in bnds_list):
                    continue
                else:
                    bnds_list.append(bnds)

                vp = VariantPlus.from_bnds(bnds)
                vplist.append(vp)
            
    return vplist


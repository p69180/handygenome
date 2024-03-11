import pysam
import Bio.Seq
import pandas as pd

import handygenome.signature.misc as signature_misc
import handygenome.deco as deco
import handygenome.refgenome.refgenome as refgenome


#@deco.get_deco_arg_choices({'refver': common.AVAILABLE_REFVERS_PLUSNONE})
@deco.get_deco_num_set(('fasta', 'refver'), 1)
def get_sbs_context(
    chrom, pos, fasta=None, refver=None, pre=1, post=1, to_pyrimidine=False,
):
    if fasta is None:
        fasta = refgenome.get_fasta(refver)

    start = pos - 1 - pre
    end = pos + post
    seq = fasta.fetch(chrom, start, end).upper()
    center_base = seq[pre]
    if (center_base in 'AG') and to_pyrimidine:
        seq = Bio.Seq.reverse_complement(seq)

    return seq


def get_sbs96_catalogue_vcfspecs(vcfspec_iter, refver):
    """Multialleleic records and non-snv records are ignored."""

    def skip_vcfspec(vcfspec):
        """skipped if True"""
        if len(vcfspec.alts) != 1:  # multialleleic record
            return True
        else:
            if not (len(vcfspec.ref) == len(vcfspec.alts[0]) == 1):  
                # not a snv
                return True
            else:
                return False
                
    fasta = refgenome.get_fasta(refver)
    catalogue_keys = signature_misc.get_catalogue_keys('sbs96')
    data = dict((x, 0) for x in catalogue_keys)

    for vcfspec in vcfspec_iter:
        if skip_vcfspec(vcfspec):
            continue
        # get sequences
        context = get_sbs_context(vcfspec.chrom, vcfspec.pos, fasta=fasta, 
                                  pre=1, post=1, to_pyrimidine=False)
        alt = vcfspec.alts[0]
        # get reverse complements if REF is purine
        if vcfspec.ref in 'AG':
            context = Bio.Seq.reverse_complement(context)
            alt = Bio.Seq.reverse_complement(alt)
        # get catalogue key and update data
        cat_key = f'{context[0]}[{context[1]}>{alt}]{context[2]}'
        data[cat_key] += 1

    result = pd.Series(data=data, index=catalogue_keys)

    return result


import re

import pysam
import Bio.Seq

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
ensembl_rest = importlib.import_module('.'.join([top_package_name, 'annotation', 'ensembl_rest']))
ensembl_parser = importlib.import_module('.'.join([top_package_name, 'annotation', 'ensembl_parser']))
libvcfspec = importlib.import_module('.'.join([top_package_name, 'variant', 'vcfspec']))


HGVSG_PATTERNS = {
    'snv': re.compile('(?P<contig>[^:]+):[gm]\.(?P<pos1>[0-9]+)(?P<seq1>[ACGTN])>(?P<seq2>[ACGTN])'),
    'del': re.compile('(?P<contig>[^:]+):[gm]\.(?P<pos1>[0-9]+)(_(?P<pos2>[0-9]+))?del'),
    'dup': re.compile('(?P<contig>[^:]+):[gm]\.(?P<pos1>[0-9]+)(_(?P<pos2>[0-9]+))?dup'),
    'ins': re.compile('(?P<contig>[^:]+):[gm]\.(?P<pos1>[0-9]+)_(?P<pos2>[0-9]+)ins(?P<seq1>[ACGTN]+)'),
    'delins': re.compile('(?P<contig>[^:]+):[gm]\.(?P<pos1>[0-9]+)(_(?P<pos2>[0-9]+))?delins(?P<seq1>[ACGTN]+)'),
    'inv': re.compile('(?P<contig>[^:]+):[gm]\.(?P<pos1>[0-9]+)_(?P<pos2>[0-9]+)?inv'),
    }


def vcfspec_to_hgvsg(vcfspec):
    common.check_vcfspec_monoallele(vcfspec)

    chrom = vcfspec.chrom
    pos = vcfspec.pos
    ref = vcfspec.ref
    alt = vcfspec.alts[0]
    mttype = vcfspec.get_mutation_type(0)

    if mttype == 'snv':
        result = f'{chrom}:g.{pos}{ref}>{alt}'
    elif mttype == 'mnv':
        pos2 = pos1 + (len(ref) - 1)
        result = f'{chrom}:g.{pos}_{pos2}delins{alt}'
    elif mttype == 'ins':
        inserted_seq = alt[1:]
        result = f'{chrom}:g.{pos}_{pos+1}ins{inserted_seq}'
    elif mttype == 'del':
        pos1 = pos + 1
        pos2 = pos + (len(ref) - 1)
        result = f'{chrom}:g.{pos1}_{pos2}del'
    elif mttype == 'delins':
        pos2 = pos + (len(ref) - 1)
        result = f'{chrom}:g.{pos}_{pos2}delins{alt}'

    return result


def pattern_matcher(hgvsg):
    matched_patterns = list()
    for hgvs_mttype, pat in HGVSG_PATTERNS.items():
        mat = pat.fullmatch(hgvsg)
        if mat is not None:
            matched_patterns.append((hgvs_mttype, mat))

    if len(matched_patterns) == 1:
        return matched_patterns[0]
    elif len(matched_patterns) == 0:
        raise Exception(f'hgvsg string with an unknown pattern: {hgvsg}')
    elif len(matched_patterns) > 1:
        raise Exception(f'More than one pattern match for this hgvsg '
                        f'string: {hgvsg}')


def modify_chrom(chrom, fasta):
    """
    When input 'chrom' is not included in the reference contig 
        names(given by 'fasta'), guesses an appropriate contig name and 
        returns it.
    For example, if 'fasta' contig names include 'chr1', 'chr2', ... but 
        input 'chrom' is '4', returns 'chr4'.
    """

    if chrom in fasta.references:
        return chrom
    else:
        chrom_mat = common.RE_PATS['assembled_chromosome'].fullmatch(chrom)
        if chrom_mat is not None:
            if 'chr1' in fasta.references:
                modified_chrom = 'chr' + chrom_mat.group(2)
            elif '1' in fasta.references:
                modified_chrom = chrom_mat.group(2)
            else:
                raise Exception(f'Input fasta contig names do not include '
                                f'"1" nor "chr1".')

        elif chrom in ('chrM', 'MT'):
            if 'chrM' in fasta.references:
                modified_chrom = 'chrM'
            elif 'MT' in fasta.references:
                modified_chrom = 'MT'
            else:
                raise Exception(
                    f'Input "chrom" is mitochondrial but input fasta object '
                    f'does not include a mitonchondrial contig.')

        else:
            raise Exception(
                f'Input "chrom" ({chrom}) does not fit to the pattern of '
                f'numbered chromosomes or mitochondria.')

        return modified_chrom


def hgvsg_to_vcfspec(hgvsg, fasta, leftmost=True):
    def sanitycheck_pos1_lt_pos2(mat, hgvsg):
        if mat.group('pos2') is not None:
            if int(mat.group('pos1')) >= int(mat.group('pos2')):
                raise Exception(f'pos1 is not less than pos2 in the input '
                                f'hgvsg string: {hgvsg}')

    hgvs_mttype, mat = pattern_matcher(hgvsg)
    chrom = modify_chrom(mat.group('contig'), fasta) 

    #############

    if hgvs_mttype == 'snv':
        pos = int(mat.group('pos1'))
        ref = mat.group('seq1')
        alt = mat.group('seq2')

    elif hgvs_mttype == 'del':
        sanitycheck_pos1_lt_pos2(mat, hgvsg)

        pos = int(mat.group('pos1')) - 1
        fetch_end = (pos + 1 
                     if (mat.group('pos2') is None) else 
                     int(mat.group('pos2')))
        ref = fasta.fetch(chrom, pos - 1, fetch_end)
        alt = ref[0]

    elif hgvs_mttype == 'dup':
        sanitycheck_pos1_lt_pos2(mat, hgvsg)

        pos = int(mat.group('pos1')) - 1
        fetch_end = (pos + 1 
                     if (mat.group('pos2') is None) else 
                     int(mat.group('pos2')))
        alt = fasta.fetch(chrom, pos - 1, fetch_end)
        ref = alt[0]

    elif hgvs_mttype == 'ins':
        if abs( int(mat.group('pos1')) - int(mat.group('pos2')) ) != 1: 
            raise Exception(
                f'Input hgvsg pattern is assumed to be insertion but '
                f'difference between pos1 and pos2 is not 1 : {hgvsg}')
        
        pos = min( int(mat.group('pos1')), int(mat.group('pos2')) )
            # In COSMIC v95 GRCh38 CosmicNCV.tsv.gz file, 
            #   HGVSG for COSV58999368 is 18:g.12817321_12817320insA
        ref = fasta.fetch(chrom, pos - 1, pos)
        alt = ref + mat.group('seq1')

    elif hgvs_mttype == 'delins':
        sanitycheck_pos1_lt_pos2(mat, hgvsg)

        pos = int(mat.group('pos1')) - 1
        fetch_end = (pos + 1 
                     if (mat.group('pos2') is None) else 
                     int(mat.group('pos2')))
        ref = fasta.fetch(chrom, pos - 1, fetch_end)
        alt = ref[0] + mat.group('seq1')

    elif hgvs_mttype == 'inv':
        sanitycheck_pos1_lt_pos2(mat, hgvsg)

        pos = int(mat.group('pos1')) - 1
        fetch_end = int(mat.group('pos2'))
        ref = fasta.fetch(chrom, pos - 1, fetch_end)
        alt = ref[0] + Bio.Seq.reverse_complement(ref[1:])

    #############

    vcfspec = libvcfspec.Vcfspec(chrom, pos, ref, [alt])

    # when ref or alt contains 'n' or 'N', functions below do not 
    #   effectively work
    if leftmost:
        vcfspec = vcfspec.leftmost()

    return vcfspec
        

def hgvsg_to_vcfspec_refver(hgvsg, refver):
    assert refver in ('hg19', 'hg38')

    return hgvsg_to_vcfspec_fastapath(hgvsg, 
                                      common.DEFAULT_FASTA_PATHS[refver])


def hgvsg_to_vcfspec_fastapath(hgvsg, fasta_path):
    with pysam.FastaFile(fasta_path) as fasta:
        result = hgvsg_to_vcfspec(hgvsg, fasta)

    return result


def hgvsc_to_hgvsg(hgvsc, hg19):
    """
    Runs ensembl rest variantrecoder
    """

    assemblyspec = importlib.import_module(
        '.'.join([top_package_name, 'assemblyspec']))

    if hg19:
        raw_result = ensembl_rest.get_url_contents(
            f'http://grch37.rest.ensembl.org/variant_recoder/human/{hgvsc}')
    else:
        raw_result = ensembl_rest.get_url_contents(
            f'http://rest.ensembl.org/variant_recoder/human/{hgvsc}')

    if len(raw_result) != 1:
        raise Exception(f'Unexpected variantrecoder raw result '
                        f'pattern: {raw_result}')

    tmp1 = raw_result[0] # a 1-length dict
    if len(tmp1) != 1:
        raise Exception(f'Unexpected variantrecoder raw result '
                        f'pattern: {raw_result}')

    tmp2 = next(tmp1.values().__iter__())
    hgvsg_results = tmp2['hgvsg']
    hgvsg_results_filtered = [x for x in hgvsg_results 
                              if not x.startswith('LRG')]
    if len(hgvsg_results_filtered) != 1:
        raise Exception(f'Unexpected variantrecoder hgvsg results '
                        f'pattern: {hgvsg_results}')

    raw_hgvsg = hgvsg_results_filtered[0]
    raw_hgvsg_split = raw_hgvsg.split(':')
    assert len(raw_hgvsg_split) == 2
    raw_contig = raw_hgvsg_split[0]

    namemap = (assemblyspec.SPECS['grch37'] if hg19 else 
               assemblyspec.SPECS['grch38'])
    if raw_contig not in namemap.data['refseq']:
        raise Exception(f'contig name of variantrecoder result is not '
                        f'refseq: {raw_hgvsg}')
    new_contig = namemap.dicts[('refseq'), ('ucsc')][raw_contig]

    if hg19:
        new_contig = re.sub('^chr', '', new_contig)

    result = new_contig + ':' + raw_hgvsg_split[1]

    return result


def hgvsc_to_vcfspec(hgvsc, hg19, fasta, leftmost=True):
    hgvsg = hgvsc_to_hgvsg(hgvsc, hg19)
    vcfspec = hgvsg_to_vcfspec(hgvsg, fasta, leftmost=leftmost)

    return vcfspec


#def _hgvsc_to_hgvsg(hgvsc, hg19):
#    def get_rest_map_result(ID, start1, end1, hg19):
#        raw_result = ensembl_rest.map(ID, start1, end1, mode = 'cds', hg19 = hg19)
#        chrom, start1, end1, is_forward = ensembl_parser.parse_rest_map(raw_result, adjust_chrom = False)
#
#        return chrom, start1, end1, is_forward
#
#    hgvs_mttype, mat = pattern_matcher(hgvsc)
#
#    contig = mat.group('contig')
#    ID = contig.split('.')[0]
#    pos1 = int(mat.group('pos1'))
#    if 'pos2' in mat.groupdict():
#        pos2 = pos1 if (mat.group('pos2') is None) else int(mat.group('pos2'))
#    else:
#        pos2 = pos1
#    chrom, new_pos1, new_pos2, is_forward = get_rest_map_result(ID, pos1, pos2, hg19)
#
#    if hgvs_mttype == 'snv':
#        if is_forward:
#            seq1, seq2 = mat.group('seq1'), mat.group('seq2')
#        else:
#            seq1, seq2 = Bio.Seq.complement(mat.group('seq1')), Bio.Seq.complement(mat.group('seq2'))
#
#        hgvsg = f'{chrom}:g.{new_pos1}{seq1}>{seq2}'
#
#    elif hgvs_mttype == 'del':
#        if mat.group('pos2') is None:
#            hgvsg = f'{chrom}:g.{new_pos1}del'
#        else:
#            hgvsg = f'{chrom}:g.{new_pos1}_{new_pos2}del'
#
#    elif hgvs_mttype == 'dup':
#        if mat.group('pos2') is None:
#            hgvsg = f'{chrom}:g.{new_pos1}dup'
#        else:
#            hgvsg = f'{chrom}:g.{new_pos1}_{new_pos2}dup'
#
#    elif hgvs_mttype == 'ins':
#        if is_forward:
#            seq1 = mat.group('seq1')
#        else:
#            seq1 = Bio.Seq.reverse_complement(mat.group('seq1'))
#
#        hgvsg = f'{chrom}:g.{new_pos1}_{new_pos2}ins{seq1}'
#
#    elif hgvs_mttype == 'delins':
#        if is_forward:
#            seq1 = mat.group('seq1')
#        else:
#            seq1 = Bio.Seq.reverse_complement(mat.group('seq1'))
#
#        if mat.group('pos2') is None:
#            hgvsg = f'{chrom}:g.{new_pos1}delins{seq1}'
#        else:
#            hgvsg = f'{chrom}:g.{new_pos1}_{new_pos2}delins{seq1}'
#
#    elif hgvs_mttype == 'inv':
#        hgvsg = f'{chrom}:g.{new_pos1}_{new_pos2}inv'
#
#    return hgvsg

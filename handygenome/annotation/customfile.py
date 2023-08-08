import os
import warnings
import re
import logging

import pysam

import handygenome.annotation.ensembl_feature as ensembl_feature
import handygenome.variant.vcfspec as libvcfspec


DEFAULT_FETCHWIDTH = 10

PAT_GENESET_ID = re.compile('ID=(?P<type>[^:]+):(?P<id>[^;]+)')
    #   type: CDS, gene, transcript (obtained from: 
    #        zcat geneset_gff3_sorted.gz | 
    #            awk '
    #            {
    #                if ($9 ~ /ID=([^:]+):/) {
    #                    print gensub(/^.*ID=([^:]+):.*$/, "\\1", "g", $9)
    #                }
    #            }' | 
    #            sort | uniq)
PAT_GENESET_PARENT = re.compile('Parent=[^:]+:(?P<id>[^;]+)')
PAT_REGULATORY_ID = re.compile('(.+;)?id=([^;]+)(;.+)?')


# VCF

class NoVcfIndexError(Exception):
    pass


def fetch_relevant_vr_multialt(vcfspec, vcf, search_equivs=True, raise_with_multihit=True):
    result = list()
    for vcfspec_monoalt in vcfspec.iter_monoalts():
        result.append(
            fetch_relevant_vr(
                vcfspec_monoalt, 
                vcf, 
                search_equivs=search_equivs, 
                raise_with_multihit=raise_with_multihit,
            )
        )

    return result


def fetch_relevant_vr(vcfspec, vcf, search_equivs=True, raise_with_multihit=True):
    """Args:
        vcfspec: Must be with one ALT allele.
        vcf: pysam.VariantFile object
        search_equivs: If True, all equivalent forms are fetched. If False, 
            only vrs with identical chrom, pos, ref, alt are fetched.
    """

    def matcher_exact(query_vcfspec, target_vcfspec):
        return query_vcfspec == target_vcfspec

    def matcher_any(query_vcfspec, target_vcfspec):
        return any((query_vcfspec == x) for x in target_vcfspec.iter_monoalts())

    def search_equivalents(vcfspec, vcf, matcher):
        if vcfspec.chrom in vcf.header.contigs:
            equivs = vcfspec.get_equivalents()
            poslist = [x.pos for x in equivs]
            fetchresult = vcf.fetch(contig=vcfspec.chrom, 
                                    start=(min(poslist) - 1),
                                    end=max(poslist))

            relevant_vr_candidates = list()
            for vr in fetchresult:
                vr_vcfspec = libvcfspec.make_parsimonious(
                    libvcfspec.Vcfspec.from_vr(vr)
                )
                if any(matcher(query, vr_vcfspec) for query in equivs):
                    relevant_vr_candidates.append(vr)
        else:
            relevant_vr_candidates = list()

        return relevant_vr_candidates

    def search_identical(vcfspec, vcf, matcher):
        if vcfspec.chrom in vcf.header.contigs:
            fetchresult = vcf.fetch(contig=vcfspec.chrom,
                                    start=(vcfspec.pos - 1), end=vcfspec.pos)

            relevant_vr_candidates = list()
            for vr in fetchresult:
                vr_vcfspec = libvcfspec.Vcfspec.from_vr(vr)
                if matcher(vcfspec, vr_vcfspec):
                    relevant_vr_candidates.append(vr)
        else:
            relevant_vr_candidates = list()
    
        return relevant_vr_candidates

    def get_relevant_vr(relevant_vr_candidates, raise_with_multihit):
        if len(relevant_vr_candidates) == 0:
            relevant_vr = None
        elif len(relevant_vr_candidates) == 1:
            relevant_vr = relevant_vr_candidates[0]
        else:
            if raise_with_multihit:
                e_msg = list()
                e_msg.append(f'There are more than one relevant variant records:')
                for vr in relevant_vr_candidates:
                    e_msg.append(str(vr))
                raise Exception('\n'.join(e_msg))
            else:
                relevant_vr = relevant_vr_candidates[0]

        return relevant_vr

    # main
    vcfspec.check_monoalt(raise_with_false=True)

    if vcf.index is None:
        raise NoVcfIndexError(f'Input vcf is not indexed.')

    matcher = matcher_exact

    if search_equivs:
        relevant_vr_candidates = search_equivalents(vcfspec, vcf, matcher)
    else:
        relevant_vr_candidates = search_identical(vcfspec, vcf, matcher)

    relevant_vr = get_relevant_vr(relevant_vr_candidates, raise_with_multihit)

    return relevant_vr


def fetch_relevant_vr_vcfpath(vcfspec, vcf_path, search_equivs=True):
    with pysam.VariantFile(vcf_path, 'r') as vcf:
        try:
            relevant_vr = fetch_relevant_vr(vcfspec, vcf, search_equivs=search_equivs)
        except NoVcfIndexError:
            print(f'Input vcf path: {vcf_path}')
            raise

    return relevant_vr


# cosmic


# geneset gff3

def parse_transcript_tabixline(tabixline):
    transcript = ensembl_feature.Transcript(is_missing=False)

    attrs_parsed = parse_gff3_attrs(tabixline.attributes)
    transcript['id'] = attrs_parsed['ID'].split(':')[1]
    transcript['biotype'] = attrs_parsed['biotype']

    if tabixline.strand == '+':
        transcript['is_forward'] = True
    elif tabixline.strand == '-':
        transcript['is_forward'] = False
    else:
        raise Exception(f'Unexpected transcript tabixline strand value: '
                        f'"{tabixline.strand}"')

    transcript['chrom'] = tabixline.contig
    transcript['start0'] = tabixline.start
    transcript['start1'] = transcript['start0'] + 1
    transcript['end0'] = tabixline.end
    transcript['end1'] = transcript['end0']

    if 'Name' in attrs_parsed:
        transcript['transcript_name'] = attrs_parsed['Name']
    else:
        transcript['transcript_name'] = None

    transcript['gene_id'] = attrs_parsed['Parent'].split(':')[1]

    return transcript


def fetch_transcript(chrom, start0, end0, tabixfile_geneset):
    """Returns an empty TranscriptSet object if the tabixfile does not 
    contain chrom.
    """
    transcript_set = ensembl_feature.TranscriptSet(is_missing=False)
    if chrom in tabixfile_geneset.contigs:
        try:
            fetchresult = tabixfile_geneset.fetch(chrom, start0, end0)
        except Exception as e:
            raise Exception(f'tabixfile fetch failed: chrom: {chrom}, '
                            f'start0: {start0}, end0: {end0}') from e
            
        for tabixline in fetchresult:
            if check_tabixline_is_transcript(tabixline):
                transcript = parse_transcript_tabixline(tabixline)
                transcript_set[transcript['id']] = transcript

    return transcript_set


def fetch_transcript_tabixline(chrom, start0, end0, transcript_id_list, 
                               tabixfile_geneset):
    if chrom not in tabixfile_geneset.contigs:
        raise Exception(f'chrom is not contained in the tabixfile.')

    candidates = list()
    fetched = tabixfile_geneset.fetch(chrom, start0, end0)
    for tabixline in fetched:
        mat = PAT_GENESET_ID.search(tabixline.attributes)
        if mat is not None:
            if mat.group('type') == 'transcript':
                if mat.group('id') in transcript_id_list:
                    candidates.append(tabixline)

    if len(candidates) != len(transcript_id_list):
        raise Exception(f'Unsuccessful transcript fetch: '
                        f'transcript_id_list: {transcript_id_list}, '
                        f'coords: ({chrom}, {start0}, {end0})')

    return candidates


def check_geneset_tabixline_type(tabixline, type):
    assert type in ('gene', 'transcript', 'CDS', 'exon', 'UTR'), (
        f'Allowed "type" values are: '
        f'"gene", "transcript", "CDS", "exon", "UTR"')

    if type == 'exon':
        return tabixline.feature == 'exon'
    elif type == 'UTR':
        return tabixline.feature in ('five_prime_UTR', 'three_prime_UTR')
    else:
        mat = PAT_GENESET_ID.search(tabixline.attributes)
        if mat is None:
            return False
        else:
            return (mat.group('type') == type)


def check_tabixline_is_transcript(tabixline):
    mat = PAT_GENESET_ID.search(tabixline.attributes)
    if mat is None:
        return False
    else:
        return (mat.group('type') == 'transcript')


def get_parent(tabixline):
    mat = PAT_GENESET_PARENT.search(tabixline.attributes)
    if mat is None:
        return None
    else:
        return mat.group('id')


# repeats bed

def parse_repeat_tabixline(tabixline):
    repeat = ensembl_feature.Repeat(is_missing=False)

    repeat['chrom'] = tabixline.contig
    repeat['start0'] = tabixline.start
    repeat['end0'] = tabixline.end
    repeat['start1'] = repeat['start0'] + 1
    repeat['end1'] = repeat['end0']

    raw_split = tabixline.name.split(':')

    clsfml = raw_split[0].split('/')
    if len(clsfml) == 1:
        repeat['class'] = clsfml[0]
        repeat['family'] = clsfml[0]
    elif len(clsfml) == 2:
        repeat['class'] = clsfml[0]
        repeat['family'] = clsfml[1]
    else:
        raise Exception(f'Unexpected repeat name pattern: {tabixline.name}')

    repeat['name'] = raw_split[1]

    return repeat


def fetch_repeat(chrom, start0, end0, tabixfile_repeats):
    """Returns an empty EnsemblFeatureSet object if the tabixfile does not 
    contain chrom.
    """
    repeat_set = ensembl_feature.RepeatSet(is_missing=False)
    if chrom in tabixfile_repeats.contigs:
        try:
            fetchresult = tabixfile_repeats.fetch(chrom, start0, end0)
        except Exception as e:
            raise Exception(f'tabixfile fetch failed: chrom: {chrom}, '
                            f'start0: {start0}, end0: {end0}') from e

        for tabixline in fetchresult:
            repeat = parse_repeat_tabixline(tabixline)
            repeat_set[repeat['id']] = repeat

    return repeat_set


# regulatory gff3

def parse_regulatory_tabixline(tabixline):
    regulatory = ensembl_feature.Regulatory(is_missing=False)

    regulatory['chrom'] = tabixline.contig
    regulatory['start0'] = tabixline.start
    regulatory['end0'] = tabixline.end
    regulatory['start1'] = regulatory['start0'] + 1
    regulatory['end1'] = regulatory['end0']
    regulatory['biotype'] = tabixline.feature

    attrs_parsed = parse_gff3_attrs(tabixline.attributes)
    regulatory['id'] = attrs_parsed['id']
    regulatory['activity'] = dict(
        x.split(':') for x in attrs_parsed['activity'].split(',')
    )

    regulatory['bound_start1'] = int(attrs_parsed['bound_start'])
    regulatory['bound_end1'] = int(attrs_parsed['bound_end'])
    regulatory['bound_start0'] = regulatory['bound_start1'] - 1
    regulatory['bound_end0'] = regulatory['bound_end1']

    return regulatory


def fetch_regulatory_tabixline(chrom, start0, end0, regulatory_id_list, 
                               tabixfile_regulatory):
    if chrom not in tabixfile_regulatory.contigs:
        raise Exception(f'chrom is not contained in the tabixfile.')

    candidates = list()
    fetched = tabixfile_regulatory.fetch(chrom, start0, end0)
    for tabixline in fetched:
        ID = get_ID_regulatory_tabixline(tabixline)
        if ID in regulatory_id_list:
            candidates.append(tabixline)

    if len(candidates) != len(regulatory_id_list):
        raise Exception(f'Unsuccessful regulatory fetch: '
                        f'regulatory_id_list: {regulatory_id_list}, '
                        f'coords: ({chrom}, {start0}, {end0})')

    return candidates


def fetch_regulatory(chrom, start0, end0, tabixfile_regulatory):
    """Returns an empty AnnotItemDict object if the tabixfile does not 
    contain chrom.
    """

    regulatory_set = ensembl_feature.RegulatorySet(is_missing=False)
    if chrom in tabixfile_regulatory.contigs:
        try:
            fetchresult = tabixfile_regulatory.fetch(chrom, start0, end0)
        except Exception as e:
            raise Exception(f'tabixfile fetch failed: chrom: {chrom}, '
                            f'start0: {start0}, end0: {end0}') from e

        for tabixline in fetchresult:
            regulatory = parse_regulatory_tabixline(tabixline)
            regulatory_set.add_feature(regulatory)

    return regulatory_set


def get_ID_regulatory_tabixline(tabixline):
    return PAT_REGULATORY_ID.search(tabixline.attributes).group(2)

# misc.

def parse_gff3_attrs(raw_string):
    return dict(x.split('=') for x in raw_string.split(';'))


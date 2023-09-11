import re

import handygenome.tools as tools
import handygenome.network as network
import handygenome.deco as deco
import handygenome.refgenome.refgenome as refgenome


PREFIX_LOOKUP_ID = refgenome.RefverDict({
    'GRCh37': 'http://grch37.rest.ensembl.org/lookup/id',
    'GRCh38': 'http://rest.ensembl.org/lookup/id',
})

PREFIX_LOOKUP_SYMBOL = refgenome.RefverDict({
    'GRCh37': 'http://grch37.rest.ensembl.org/lookup/symbol',
    'GRCh38': 'http://rest.ensembl.org/lookup/symbol',
})

PREFIX_REGULATORY = refgenome.RefverDict({
    'GRCh37': 'http://grch37.rest.ensembl.org/regulatory/species/homo_sapiens/id',
    'GRCh38': 'http://rest.ensembl.org/regulatory/species/homo_sapiens/id',
})

PREFIX_VEP = refgenome.RefverDict({
    'GRCh37': 'http://grch37.rest.ensembl.org/vep/human/hgvs',
    'GRCh38': 'http://rest.ensembl.org/vep/human/hgvs',
    'GRCm39': 'http://rest.ensembl.org/vep/mouse/hgvs',
})

PREFIX_OVERLAP = refgenome.RefverDict({
    'GRCh37': 'http://grch37.rest.ensembl.org/overlap/region',
    'GRCh38': 'http://rest.ensembl.org/overlap/region',
})

PREFIX_MAP = refgenome.RefverDict({
    'GRCh37': 'http://grch37.rest.ensembl.org/map',
    'GRCh38': 'http://rest.ensembl.org/map',
})

PREFIX_SEQUENCE = refgenome.RefverDict({
    'GRCh37': 'http://grch37.rest.ensembl.org/sequence/id',
    'GRCh38': 'http://rest.ensembl.org/sequence/id',
    'GRCm39': 'http://rest.ensembl.org/sequence/id',
})

HTTP_HEADERS_POST = {
    'Content-Type': 'application/json',
    'Accept': 'application/json',
}

HTTP_HEADERS_GET = {
    'Content-Type': 'application/json',
    'Accept': 'application/json',
}

MAX_POST_SIZE = {
    'lookup_id': 1000,
    'lookup_symbol': 1000,
    'vep': 200,
    'sequence': 50,
}


def lookup_id(ID, refver, expand=False):
    """start, end are 1-based inclusive coordinates"""
    prefix = PREFIX_LOOKUP_ID[refver]
    url = '/'.join([prefix, ID])
    params = {'expand': int(expand)}

    return network.http_get(url, params=params, headers=HTTP_HEADERS_GET)


def lookup_id_post(ID_list, refver, expand=False):
    """start, end are 1-based inclusive coordinates"""
    result = list()

    url = PREFIX_LOOKUP_ID[refver]
    params = {'expand' : int(expand)}
    for sublist in tools.chunk_iter(ID_list, MAX_POST_SIZE['lookup_id']):
        data = {'ids' : sublist}
        subresult = network.http_post(url, data, params=params, headers=HTTP_HEADERS_POST)
        result.extend(subresult.values())

    return result


def lookup_symbol(symbol, refver, expand=False):
    prefix = PREFIX_LOOKUP_SYMBOL[refver]
    species = refgenome.RefverDict.find_species(refver)
    url = '/'.join([prefix, species, symbol])
    params = {'expand' : int(expand)}

    return network.http_get(url, params=params, headers=HTTP_HEADERS_GET)


def lookup_symbol_post(symbols, refver, expand=False):
    prefix = PREFIX_LOOKUP_SYMBOL[refver]
    species = refgenome.RefverDict.find_species(refver)
    url = '/'.join([prefix, species])
    params = {'expand' : int(expand)}

    result = list()
    for sublist in tools.chunk_iter(symbols, MAX_POST_SIZE['lookup_symbol']):
        data = {'symbols' : sublist}
        subresult = network.http_post(url, data, params=params, headers=HTTP_HEADERS_POST)
        result.extend(subresult.values())

    return result


def regulatory(ID, refver, activity=True):
    prefix = PREFIX_REGULATORY[refver]
    url = '/'.join([prefix, ID])
    params = {'activity' : int(activity)}

    result = network.http_get(url, params=params, headers=HTTP_HEADERS_GET)
    if isinstance(result, list) and len(result) == 1:
        result = result[0]
    else:
        raise Exception(f'REST regulatory result for ID {ID} is not a list with length of 1.')

    return result


def vep(
    refver, hgvsg=None, vcfspec=None, distance=5000, 
    with_CADD=True, with_Phenotypes=False, with_canonical=True, 
    with_mane=True, with_miRNA=False, with_numbers=True, 
    with_protein=True, with_ccds=True, with_hgvs=True,
):
    if hgvsg is None:
        if vcfspec is None: 
            raise Exception(f'When "hgvsg" is not set, "vcfspec" must be set.')
        vcfspec.check_monoalt(raise_with_false=True)
        hgvsg = vcfspec.to_hgvsg(alt_index=0)

    prefix = PREFIX_VEP[refver]
    url = '/'.join([prefix, hgvsg])

    params = {
        'distance' : distance,
        'CADD' : int(with_CADD),
        'canonical' : int(with_canonical),
        'mane' : int(with_mane),
        'miRNA' : int(with_miRNA),
        'numbers' : int(with_numbers),
        'protein' : int(with_protein),
        'ccds' : int(with_ccds),
        'hgvs' : int(with_hgvs),
    }
    if with_Phenotypes:
        # Adding 'Phenotypes' to params results in phenotype annotation
            # even when it is set to False
        params.update({'Phenotypes' : int(with_Phenotypes)})

    return network.http_get(url, params=params, headers=HTTP_HEADERS_GET)


def vep_post(
    refver, vcfspec_list=None, hgvsg_list=None, distance=5000,
    with_CADD=True, with_Phenotypes=False, with_canonical=True,
    with_mane=True, with_miRNA=False, with_numbers=True, 
    with_protein=True, with_ccds=True, with_hgvs=True,
):
    url = PREFIX_VEP[refver]
    if hgvsg_list is None:
        if vcfspec_list is None: 
            raise Exception(
                f'When "hgvsg_list" is not set, "vcfspec_list" must be set.'
            )
        if not all(
            vcfspec.check_monoalt(raise_with_false=False)
            for vcfspec in vcfspec_list
        ):
            raise Exception(f'Components of "vcfspec_list" must be all monoalt.')
        hgvsg_list = [vcfspec.to_hgvsg(alt_index=0) for vcfspec in vcfspec_list]

    params = {
        'distance': distance,
        'CADD': int(with_CADD),
        'canonical': int(with_canonical),
        'mane': int(with_mane),
        'miRNA': int(with_miRNA),
        'numbers': int(with_numbers),
        'protein': int(with_protein),
        'ccds': int(with_ccds),
        'hgvs': int(with_hgvs),
    }
    if with_Phenotypes:
        params.update({'Phenotypes': int(with_Phenotypes)})

    result = list()
    for sublist in tools.chunk_iter(hgvsg_list, MAX_POST_SIZE['vep']):
        data = {'hgvs_notations': sublist}
        result.extend(network.http_post(url, data, params=params, headers=HTTP_HEADERS_POST))

    return result


def overlap(
    chrom, start1, end1, refver,
    transcript=True, regulatory=True, motif=False, repeat=False,
):
    """start1, end1: 1-based closed system"""

    prefix = PREFIX_OVERLAP[refver]
    species = refgenome.RefverDict.find_species(refver)

    suffix_list = list()
    if transcript: 
        suffix_list.append('feature=transcript')
    if regulatory: 
        suffix_list.append('feature=regulatory')
    if motif: 
        suffix_list.append('feature=motif')
    if repeat: 
        suffix_list.append('feature=repeat')

    url = '/'.join([prefix, species, f'{chrom}:{start1}-{end1}'])
    if len(suffix_list) > 0:
        url += '?' + ';'.join(suffix_list)

    return network.http_get(url, headers=HTTP_HEADERS_GET)


@deco.get_deco_arg_choices({'mode': ('cdna', 'cds', 'translation')})
def map(ID, start1, end1, mode, refver):
    prefix = PREFIX_MAP[refver]
    url = '/'.join([prefix, mode, ID, f'{start1}..{end1}'])
    return network.http_get(url, headers=HTTP_HEADERS_GET)
    

def sequence_get(ID, refver, type='genomic'):
    prefix = PREFIX_SEQUENCE[refver]
    url = '/'.join([prefix, ID])
    params = {
        'type': type,
    }

    return network.http_get(
        url, 
        params=params, 
        headers={'Content-Type': 'text/plain'},
        text=True,
    )


def sequence_post(ID_list, refver, type='genomic'):
    url = PREFIX_SEQUENCE[refver]
    params = {
        'type': type,
    }

    result = list()
    for sublist in tools.chunk_iter(ID_list, MAX_POST_SIZE['sequence']):
        data = {'ids': sublist}
        result.extend(
            network.http_post(
                url, 
                data, 
                params=params, 
                headers=HTTP_HEADERS_POST,
            )
        )

    return result




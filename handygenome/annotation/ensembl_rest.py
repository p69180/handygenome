import re

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
hgvs_module = importlib.import_module('.'.join([top_package_name, 'hgvs']))


REFVERS_ALLOWED = ('hg19', 'hg38', 'mm39')

SPECIES = {
    'hg19': 'homo_sapiens',
    'hg38': 'homo_sapiens',
    'mm39': 'mus_musculus',
    }


def into_37(s):
    return re.sub('http://', 'http://grch37.', s)

PREFIX_ID = 'http://rest.ensembl.org/lookup/id'
PREFIX_ID_37 = into_37(PREFIX_ID)

PREFIX_SYMBOL = 'http://rest.ensembl.org/lookup/symbol'
PREFIX_SYMBOL_37 = into_37(PREFIX_SYMBOL)

PREFIX_REGULATORY = 'http://rest.ensembl.org/regulatory/species/homo_sapiens/id'
PREFIX_REGULATORY_37 = into_37(PREFIX_REGULATORY)

PREFIX_VEP = 'http://rest.ensembl.org/vep/human/hgvs'
PREFIX_VEP_37 = into_37(PREFIX_VEP)
PREFIX_VEP_MOUSE = 'http://rest.ensembl.org/vep/mouse/hgvs'

PREFIX_OVERLAP = 'http://rest.ensembl.org/overlap/region'
PREFIX_OVERLAP_37 = into_37(PREFIX_OVERLAP)

PREFIX_MAP = 'https://rest.ensembl.org/map'
PREFIX_MAP_37 = into_37(PREFIX_MAP)

HTTP_HEADERS_POST = {'Content-Type': 'application/json', 
                     'Accept': 'application/json'}
HTTP_HEADERS_GET = {'Content-Type': 'application/json', 
                    'Accept': 'application/json'}


@common.get_deco_arg_choices({'refver': REFVERS_ALLOWED})
def lookup_id(ID, refver, expand=True):
    """start, end are 1-based inclusive coordinates"""

    if refver == 'hg19':
        prefix = PREFIX_ID_37
    else:
        prefix = PREFIX_ID

    url = '/'.join([prefix, ID])
    params = {'expand': int(expand)}

    return common.http_get(url, params, HTTP_HEADERS_GET)


@common.get_deco_arg_choices({'refver': REFVERS_ALLOWED})
def lookup_id_post(IDs, refver, expand=True):
    """
    start, end are 1-based inclusive coordinates
    """

    if refver == 'hg19':
        url = PREFIX_ID_37
    else:
        url = PREFIX_ID

    params = { 'expand' : int(expand) }
    data = { 'ids' : IDs }

    return common.http_post(url, data, params, HTTP_HEADERS_POST)


@common.get_deco_arg_choices({'refver': REFVERS_ALLOWED})
def lookup_symbol(symbol, refver, expand=False):
    prefix = PREFIX_SYMBOL_37 if (refver == 'hg19') else PREFIX_SYMBOL
    species = SPECIES[refver]
    url = '/'.join( [ prefix, species, symbol ] )
    params = { 'expand' : int(expand) }

    return common.http_get(url, params, HTTP_HEADERS_GET)


@common.get_deco_arg_choices({'refver': REFVERS_ALLOWED})
def lookup_symbol_post(symbols, refver, expand=False):
    species = SPECIES[refver]
    if refver == 'hg19':
        url = '/'.join([PREFIX_SYMBOL_37, species])
    else:
        url = '/'.join([PREFIX_SYMBOL, species])

    params = { 'expand' : int(expand) }
    data = { 'symbols' : list(symbols) }

    return common.http_post(url, data, params, HTTP_HEADERS_POST)


# unavailable for mouse
@common.get_deco_arg_choices({'refver': ('hg19', 'hg38')})
def regulatory(ID, refver):
    prefix = PREFIX_REGULATORY_37 if (refver == 'hg19') else PREFIX_REGULATORY
    url = '/'.join([prefix, ID])
    params = { 'activity' : 1 }

    result = common.http_get(url, params, HTTP_HEADERS_GET)
    if isinstance(result, list) and len(result) == 1:
        result = result[0]
    else:
        raise Exception(f'REST regulatory result for ID {ID} is not a list with length of 1.')

    return result

@common.get_deco_arg_choices({'refver': REFVERS_ALLOWED})
#@common.get_deco_num_set_differently(('vcfspec', 'hgvsg'), 1)
def vep(refver, hgvsg=None, vcfspec=None, distance=5000, 
        with_CADD=True, with_Phenotypes=False, with_canonical=True, 
        with_mane=True, with_miRNA=False, with_numbers=True, 
        with_protein=True, with_ccds=True, with_hgvs=True):
    if refver == 'hg19':
        prefix = PREFIX_VEP_37
    elif refver == 'hg38':
        prefix = PREFIX_VEP
    elif refver == 'mm39':
        prefix = PREFIX_VEP_MOUSE
    
    if hgvsg is None:
        hgvsg = hgvs_module.vcfspec_to_hgvsg(vcfspec)

    url = '/'.join([prefix, hgvsg])

    params = {
        'distance' : distance,
        'CADD' : int(with_CADD),
        #'Phenotypes' : int(with_Phenotypes),
        'canonical' : int(with_canonical),
        'mane' : int(with_mane),
        'miRNA' : int(with_miRNA),
        'numbers' : int(with_numbers),
        'protein' : int(with_protein),
        'ccds' : int(with_ccds),
        'hgvs' : int(with_hgvs),
    }
    if with_Phenotypes:
        params.update( { 'Phenotypes' : int(with_Phenotypes) } )

    return common.http_get(url, params, HTTP_HEADERS_GET)


@common.get_deco_arg_choices({'refver': REFVERS_ALLOWED})
@common.get_deco_num_set_differently(('vcfspec_list', 'hgvsg_list'), 1)
def vep_post(refver, vcfspec_list=None, hgvsg_list=None, distance=5000,
             with_CADD=True, with_Phenotypes=False, with_canonical=True,
             with_mane=True, with_miRNA=False, with_numbers=True, 
             with_protein=True, with_ccds=True, with_hgvs=True):
    if refver == 'hg19':
        url = PREFIX_VEP_37
    elif refver == 'hg38':
        url = PREFIX_VEP
    elif refver == 'mm39':
        url = PREFIX_VEP_MOUSE

    if hgvsg_list is None:
        hgvsg_list = [hgvs_module.vcfspec_to_hgvsg(vcfspec) 
                      for vcfspec in vcfspec_list]

    data = {'hgvs_notations': hgvsg_list}

    params = {
        'distance': distance,
        'CADD': int(with_CADD),
        #'Phenotypes': int(with_Phenotypes),
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

    return common.http_post(url, data, params, HTTP_HEADERS_POST)


@common.get_deco_arg_choices({'refver': REFVERS_ALLOWED})
def overlap(chrom, start1, end1, refver,
            transcript=True, regulatory=True, motif=False, repeat=False):
    """start1, end1: 1-based closed system"""

    prefix = PREFIX_OVERLAP_37 if (refver == 'hg19') else PREFIX_OVERLAP
    species = SPECIES[refver]

    suffix_list = list()
    if transcript: 
        suffix_list.append('feature=transcript')
    if regulatory: 
        suffix_list.append('feature=regulatory')
    if motif: 
        suffix_list.append('feature=motif')
    if repeat: 
        suffix_list.append('feature=repeat')
    suffix = '?' + ';'.join(suffix_list)

    url = '/'.join([prefix, species, f'{chrom}:{start1}-{end1}']) + suffix

    return common.http_get(url, headers=HTTP_HEADERS_GET)


@common.get_deco_arg_choices({'refver': REFVERS_ALLOWED, 
                              'mode': ('cdna', 'cds', 'translation')})
def map(ID, start1, end1, mode, refver):
    if refver == 'hg19':
        prefix = '/'.join([PREFIX_MAP_37, mode])
    else:
        prefix = '/'.join([PREFIX_MAP, mode])

    url = '/'.join([prefix, ID, f'{start1}..{end1}'])

    return common.http_get(url, headers=HTTP_HEADERS_GET)
    

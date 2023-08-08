import os
import functools
import re
import itertools

import pandas as pd
import pyranges as pr

import handygenome
import handygenome.network as network
import handygenome.cnv.misc as cnvmisc
import handygenome.refgenome as refgenome


URL_BASE = 'https://api.genome.ucsc.edu'

CYTOBAND_URLS = refgenome.RefverDict({
    'GRCh37': 'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz', 
    'GRCh38': 'http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz', 
})

CYTOBAND_DIR = os.path.join(handygenome.DIRS['data'], 'cytoband')
os.makedirs(CYTOBAND_DIR, exist_ok=True)


def tohexcode(tup):
    return '#' + ''.join(hex(x)[2:] for x in tup)

# https://www.wikiwand.com/en/Chromosome_1
CYTOBAND_COLORMAP = {
    #'gvar': tohexcode((224, 224, 224)),
    'gvar': 'tab:pink',
    'stalk': tohexcode((112, 128, 144)),
    #'acen': tohexcode((110, 127, 143)),
    'acen': 'tab:red',
    'gneg': tohexcode((255, 255, 255)),
    'gpos25': tohexcode((217, 217, 217)),
    'gpos50': tohexcode((151, 151, 151)),
    'gpos75': tohexcode((99, 99, 99)),
    'gpos100': tohexcode((0, 0, 0)),
}


def get_avaliable_genomes():
    url = f'{URL_BASE}/list/ucscGenomes'
    raw = network.http_get(url)
    return raw['ucscGenomes']


@functools.cache
def get_refver_aliases():
    result = dict()
    genome_data = get_avaliable_genomes()
    for key, val in genome_data.items():
        if val['organism'] in ('Mouse', 'Human'):
            aliases = re.match('.*\((\S+ )?(\S+)\)', val['description']).group(2).split('/')
            result[key] = aliases

    return result


def standardize_refver(refver):
    for standard, aliases in get_refver_aliases().items():
        if refver in aliases:
            return standard
    raise Exception(f'Unknown reference version.')


def list_tracks(refver):
    refver = standardize_refver(refver)
    raw = network.http_get(f'{URL_BASE}/list/tracks?genome={refver}')
    return sorted(raw[refver].keys())


def rename_hg19_cytoband(cytoband_df):
    result = cytoband_df.copy()
    assemblyspec = refgenome.SPECS['hg19']
    new_chroms = result['Chromosome'].apply(
        lambda x: assemblyspec.convert(x, 'nochr_plus_genbank')
    )
    result['Chromosome'] = new_chroms
    result = result.loc[result['Chromosome'].notna(), :]

    return result


def get_cytoband_from_ucsc_api(refver, rename_hg19=True, as_gr=True):
    refver = standardize_refver(refver)
    if 'cytoBand' not in list_tracks(refver):
        raise Exception(f'"cytoBand" is not available for this reference genome.')

    raw = network.http_get(
        f'{URL_BASE}/getData/track?genome={refver};track=cytoBand'
    )
    data = list()
    for dic in itertools.chain.from_iterable(raw['cytoBand'].values()):
        data.append(
            (dic['chrom'], dic['chromStart'], dic['chromEnd'], dic['name'], dic['gieStain'])
        )

    chroms, starts, ends, names, stains = zip(*data)
    result = pd.DataFrame.from_dict(
        {
            'Chromosome': chroms,
            'Start': starts,
            'End': ends,
            'Name': names,
            'Stain': stains,
        }
    )
    if (refver == 'hg19') and rename_hg19:
        result = rename_hg19_cytoband(result)

    result = cnvmisc.sort_genome_df(result, refver=refver)
    if as_gr:
        result = pr.PyRanges(result)

    return result


def get_cytoband_path(refver):
    return os.path.join(CYTOBAND_DIR, f'{refver}.tsv.gz')


def write_cytoband(cytoband_df, refver):
    cytoband_df.to_csv(
        get_cytoband_path(refver),
        sep='\t',
        header=True,
        index=False,
    )


def get_cytoband(refver, as_gr=True):
    cytoband_path = get_cytoband_path(refver)
    if not os.path.exists(cytoband_path):
        cytoband_df = get_cytoband_from_ucsc_api(refver, rename_hg19=True, as_gr=False)
        write_cytoband(cytoband_df, refver)
    else:
        cytoband_df = pd.read_csv(
            cytoband_path,
            sep='\t',
            header=0,
            dtype={'Chromosome': str, 'Start': int, 'End': int, 'Name': str, 'Stain': str},
        )

    result = cytoband_df
    if as_gr:
        result = pr.PyRanges(result)

    return result
        

get_cytoband_gr = get_cytoband

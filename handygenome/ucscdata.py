import functools
import re
import itertools

import pandas as pd
import pyranges as pr

import handygenome.common as common
import handygenome.assemblyspec as libassemblyspec
import handygenome.cnv.misc as cnvmisc


URL_BASE = 'https://api.genome.ucsc.edu'


CYTOBAND_URLS = common.RefverDict({
    'GRCh37': 'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz', 
    'GRCh38': 'http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz', 
})


def get_avaliable_genomes():
    url = f'https://api.genome.ucsc.edu/list/ucscGenomes'
    raw = common.http_get(url)
    return raw['ucscGenomes']


def get_refver_aliases():
    result = dict()
    genome_data = get_avaliable_genomes()
    for key, val in genome_data.items():
        if val['organism'] in ('Mouse', 'Human'):
            aliases = re.match('.*\((\S+ )?(\S+)\)', val['description']).group(2).split('/')
            result[key] = aliases

    return result


REFVER_ALIASES = get_refver_aliases()


def standardize_refver(refver):
    for standard, aliases in REFVER_ALIASES.items():
        if refver in aliases:
            return standard
    raise Exception(f'Unknown reference version.')


def list_tracks(refver):
    refver = standardize_refver(refver)
    raw = common.http_get(f'https://api.genome.ucsc.edu/list/tracks?genome={refver}')
    return sorted(raw[refver].keys())


@functools.cache
def get_cytoband_gr(refver, rename_hg19=True, as_gr=True):
    refver = standardize_refver(refver)
    #if 'cytoBand' not in list_tracks(refver):
    #    raise Exception(f'"cytoBand" is not available for this reference genome.')
    raw = common.http_get(f'https://api.genome.ucsc.edu/getData/track?genome={refver};track=cytoBand')
    data = list()
    for dic in itertools.chain.from_iterable(raw['cytoBand'].values()):
        data.append(
            (dic['chrom'], dic['chromStart'], dic['chromEnd'], dic['name'], dic['gieStain'])
        )

    if (refver == 'hg19') and rename_hg19:
        new_data = list()
        assemblyspec = libassemblyspec.SPECS['hg19']
        for item in data:
            new_chrom = assemblyspec.convert(item[0], 'nochr_plus_genbank')
            if new_chrom is not None:
                new_data.append((new_chrom,) + item[1:])
        data = new_data

    chroms, starts, ends, names, stains = zip(*data)
    if as_gr:
        gr = pr.from_dict(
            {'Chromosome': chroms, 'Start': starts, 'End': ends, 'Name': names, 'Stain': stains}, 
            int64=False,
        )
        gr.sort()

        return gr
    else:
        df = pd.DataFrame.from_dict({
            'Chromosome': chroms,
            'Start': starts,
            'End': ends,
            'Name': names,
            'Stain': stains,
        })
        return cnvmisc.sort_genome_df(df, refver=refver)
        


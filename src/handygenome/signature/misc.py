import os
import itertools
import functools

import numpy as np
import pandas as pd

import handygenome
import handygenome.network as network
import handygenome.deco as deco


COLORS_SBS6 = {
    'C>A': np.array([3, 189, 239]) / 256,
    'C>G': np.array([1, 1, 1]) / 256,
    'C>T': np.array([228, 41, 38]) / 256,
    'T>A': np.array([203, 202, 202]) / 256,
    'T>C': np.array([162, 207, 99]) / 256, 
    'T>G': np.array([236, 199, 197]) / 256,
    'other': 'y',
}

AVAILABLE_REFVERS = ('GRCh38', 'GRCh37', 'mm10')
AVAILABLE_CAT_TYPES = ('sbs96', 'id83', 'dbs78', 'cn48')

DATA_DIR = os.path.join(handygenome.DIRS['data'], 'signature_data')
if not os.path.exists(DATA_DIR):
    os.mkdir(DATA_DIR)

# SIGNATURE_DATA_URLS
SIGNATURE_DATA_URLS = dict()
SIGNATURE_DATA_URLS['id83'] = 'https://cancer.sanger.ac.uk/signatures/documents/1907/COSMIC_v3.3_ID_GRCh37.txt'
SIGNATURE_DATA_URLS['cn48'] = 'https://cancer.sanger.ac.uk/signatures/documents/2044/COSMIC_v3.3_CN_GRCh37.txt'
SIGNATURE_DATA_URLS['GRCh37'] = {
    'sbs96': 'https://cancer.sanger.ac.uk/signatures/documents/1908/COSMIC_v3.3_SBS_GRCh37.txt',
    'id83': SIGNATURE_DATA_URLS['id83'],
    'dbs78': 'https://cancer.sanger.ac.uk/signatures/documents/1902/COSMIC_v3.3_DBS_GRCh37.txt',
    'cn48': SIGNATURE_DATA_URLS['cn48'],
    }
SIGNATURE_DATA_URLS['GRCh38'] = {
    'sbs96': 'https://cancer.sanger.ac.uk/signatures/documents/1909/COSMIC_v3.3_SBS_GRCh38.txt',
    'id83': SIGNATURE_DATA_URLS['id83'],
    'dbs78': 'https://cancer.sanger.ac.uk/signatures/documents/1903/COSMIC_v3.3_DBS_GRCh38.txt',
    'cn48': SIGNATURE_DATA_URLS['cn48'],
    }
SIGNATURE_DATA_URLS['mm10'] = {
    'sbs96': 'https://cancer.sanger.ac.uk/signatures/documents/1911/COSMIC_v3.3_SBS_mm10.txt',
    'id83': SIGNATURE_DATA_URLS['id83'],
    'dbs78': 'https://cancer.sanger.ac.uk/signatures/documents/1905/COSMIC_v3.3_DBS_mm10.txt',
    'cn48': SIGNATURE_DATA_URLS['cn48'],
    }
###

# SIGNATURE_DATA_FILE_PATHS
SIGNATURE_DATA_FILE_PATHS = dict()
for refver in AVAILABLE_REFVERS:
    subdict = dict()
    for catalogue_type in AVAILABLE_CAT_TYPES:
        url = SIGNATURE_DATA_URLS[refver][catalogue_type]
        path = os.path.join(DATA_DIR, os.path.basename(url))
        subdict[catalogue_type] = path
    SIGNATURE_DATA_FILE_PATHS[refver] = subdict
###


def get_cossim(v1, v2):
    if not (np.ndim(v1) == np.ndim(v2) == 1 and
            np.shape(v1) == np.shape(v2)):
        raise Exception(f'Each input argument must be equivalent to '
                        f'a one-dimensional array with the same length.')

    v1_2norm = np.linalg.norm(v1, ord=2, axis=None)
    v2_2norm = np.linalg.norm(v2, ord=2, axis=None)
    if v1_2norm == 0 or v2_2norm == 0:
        return np.nan
    else:
        return min(1.0, np.dot(v1, v2) / (v1_2norm * v2_2norm))


def create_catalogue_keys_sbs6(as_tuple=False):
    result = list()
    for ref in 'CT':
        alts = sorted(set('ACGT') - set(ref))
        for alt in alts:
            if as_tuple:
                item = (ref, alt)
            else:
                item = f'{ref}>{alt}'
            result.append(item)

    return tuple(result)


def create_catalogue_keys_sbs96(as_tuple=False):
    result = list()
    for ref in 'CT':
        alts = sorted(set('ACGT') - set(ref))
        for alt in alts:
            for pre, post in itertools.product('ACGT', repeat=2):
                if as_tuple:
                    item = (pre, ref, alt, post)
                else:
                    item = f'{pre}[{ref}>{alt}]{post}'
                result.append(item)

    return tuple(result)


@deco.get_deco_arg_choices({'refver': AVAILABLE_REFVERS})
@deco.get_deco_arg_choices({'catalogue_type': AVAILABLE_CAT_TYPES})
def download_signature_data(refver, catalogue_type):
    """Args:
        catalogue_type: Mutation type (e.g. sbs, id, dbs, cn)
    """
    data_url = SIGNATURE_DATA_URLS[refver][catalogue_type]
    data_file_path = SIGNATURE_DATA_FILE_PATHS[refver][catalogue_type]
    network.download(data_url, data_file_path)
    

@deco.get_deco_arg_choices({'refver': AVAILABLE_REFVERS})
@deco.get_deco_arg_choices({'catalogue_type': AVAILABLE_CAT_TYPES})
@functools.cache
def load_signature_data(refver, catalogue_type):
    """num_rows:
        sbs: 96
        id: 83
        dbs: 78
        cn: 48
    """
    data_file_path = SIGNATURE_DATA_FILE_PATHS[refver][catalogue_type]
    if not os.path.exists(data_file_path):
        download_signature_data(refver, catalogue_type)
    sigdata = pd.read_csv(data_file_path, sep='\t', index_col=0)

    return sigdata


@functools.cache
def get_catalogue_keys(catalogue_type):
    if catalogue_type in ('sbs96', 'id83', 'dbs78', 'cn48'):
        sigdata = load_signature_data('GRCh37', catalogue_type)
        return tuple(sigdata.index)
    elif catalogue_type == 'sbs6':
        return create_catalogue_keys_sbs6(as_tuple=False)
    else:
        raise Exception(f'Unavailable catalogue type')


def get_catalogue_type(pd_index):
    pd_index_set = set(pd_index)
    okay_cattypes = list()
    for catalogue_type in ('sbs96', 'id83', 'dbs78', 'cn48', 'sbs6'):
        if pd_index_set == set(get_catalogue_keys(catalogue_type)):
            okay_cattypes.append(catalogue_type)

    if len(okay_cattypes) == 1:
        return okay_cattypes[0]
    else:
        raise Exception(f'Unknown catalogue type')
        

import base64

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


def get_cosmic_auth_string(email, password):
    auth_string = base64.b64encode(f'{email}:{password}'.encode('utf-8')).decode()
    return auth_string


@common.get_deco_arg_choices({'refver': ('GRCh37', 'GRCh38')})
def get_cosmic_file_url(auth_string, refver, cosmicver, cosmicfile):
    open_url = f'https://cancer.sanger.ac.uk/cosmic/file_download/{refver}/cosmic/{cosmicver}/{cosmicfile}'
    url = common.http_get(open_url, headers={'Authorization' : f'Basic {auth_string}'})['url']
    return url


def download_cosmic_file(outfile_path, email, password, cosmicfile, cosmicver, refver):
    """Args:
        refver: "GRCh37" or "GRCh38"
        cosmicver: e.g. "v95"
        cosmicfile: e.g. "CosmicNCV.tsv.gz"
    """
    auth_string = get_cosmic_auth_string(email, password)
    cosmicfile_url = get_cosmic_file_url(auth_string, refver, cosmicver, cosmicfile)
    common.download(cosmicfile_url, outfile_path)

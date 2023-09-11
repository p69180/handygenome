import json

import handygenome
import handygenome.network as network
import handygenome.logutils as logutils
import handygenome.fastahandler as fastahandler


URL_AUTHORITY = 'ftp.ensembl.org'
URL_PATH_GENESET = '/pub/current_gff3'
URL_PATH_REGULATION = '/pub/current_regulation'
URL_PATH_GRCH37_GENESET = '/pub/grch37/current/gff3/homo_sapiens'
URL_PATH_GRCH37_REGULATION = '/pub/grch37/current/regulation/homo_sapiens'

class EnsemblPaths:

    toppaths_cache = os.path.join(
        handygenome.DATA_DIR, 
        'ensembl_toppaths.json',
    )

    def remove_toppaths_cache(self):
        os.remove(self.__class__.toppaths_cache)

    def load_toppaths_cache(self):
        with open(self.__class__.toppaths_cache, 'rt') as f:
            self.toppaths = json.load(f)

    def save_toppaths_cache(self):
        with open(self.__class__.toppaths_cache, 'wt') as f:
            json.dump(self.toppaths, f)

    ###

    def __init__(self):
        self.initialize(force_update=False)

    def initialize(self, force_update=False):
        if force_update:
            self.remove_toppaths_cache()
        self.init_toppaths()

    def init_toppaths(self):
        if os.path.exists(self.__class__.toppaths_cache):
            self.load_toppaths_cache()
        else:
            self.toppaths = dict(
                'gff3': dict(),
                'regulation': dict(),
            )

    def set_gff3_toppaths(self):

    #####




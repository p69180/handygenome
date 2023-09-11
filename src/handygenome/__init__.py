import sys
import os
import subprocess
import pathlib
import importlib
import importlib.resources
import importlib.metadata

import yaml

import handygenome


__version__ = importlib.metadata.version(handygenome.__name__)


#####################################################################################################
# set numpy threads (https://stackoverflow.com/questions/30791550/limit-number-of-threads-in-numpy) #
#####################################################################################################

os.environ["OMP_NUM_THREADS"] = "1"  # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1"  # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "1"  # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"  # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "1"  # export NUMEXPR_NUM_THREADS=6


############################
# USERDIR (~/.handygenome) #
############################

USERDIR = pathlib.Path.home() / '.handygenome'

USERDATA_DIR = USERDIR / 'data'
CONFIGPATH = USERDIR / 'config.yaml'


#####################
# CONFIG and PARAMS #
#####################

with CONFIGPATH.open() as f:
    OPTION = yaml.safe_load(f)


def get_exec_path(key, replace=None, whicharg=None):
    if OPTION[key] is None:
        if replace is None:
            if whicharg is None:
                whicharg = key

            p = subprocess.run(['which', whicharg], capture_output=True, text=True)
            if p.returncode == 0:
                return p.stdout.strip()
            else:
                return None
        else:
            return replace
    else:
        return OPTION[key]


PARAMS = dict()
PARAMS['python'] = get_exec_path('python', replace=sys.executable)
PARAMS['bash'] = get_exec_path('bash')
PARAMS['perl'] = get_exec_path('perl')
PARAMS['mosdepth'] = get_exec_path('mosdepth')
PARAMS['default_refver'] = OPTION['default_refver']


#################################
# data files within top package #
#################################

PKGPATH = importlib.resources.files(handygenome)
DATA_DIR = PKGPATH.joinpath('data')

DIRS = {
    #'data': os.path.join(PROJECT_PATH, 'data'),
    #'externals': os.path.join(PROJECT_PATH, 'externals'),
    'R': PKGPATH.joinpath('Rscripts'),
    'scripts': PKGPATH.joinpath('scripts'),
    #'package': os.path.join(PROJECT_PATH, 'handygenome'),
}

import sys
import os
import pathlib
import importlib
import importlib.resources
import importlib.metadata

import yaml

import handygenome


__version__ = importlib.metadata.version(handygenome.__name__)


# set numpy threads (https://stackoverflow.com/questions/30791550/limit-number-of-threads-in-numpy)
os.environ["OMP_NUM_THREADS"] = "1"  # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1"  # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "1"  # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"  # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "1"  # export NUMEXPR_NUM_THREADS=6


# USERDIR (~/.handygenome)
USERDIR = pathlib.Path.home() / '.handygenome'

CONFIGPATH = USERDIR / 'config.yaml'
with CONFIGPATH.open() as f:
    OPTION = yaml.safe_load(f)
if OPTION['python'] is None:
    OPTION['python'] = sys.executable

DATADIR = USERDIR / 'data'


#PROJECT_PATH = os.path.dirname(os.path.dirname(handygenome.__file__))  
#    # /home/users/pjh/scripts/python_genome_package_dev

# data files within top package
PKGPATH = importlib.resources.files(handygenome)
DIRS = {
    #'data': os.path.join(PROJECT_PATH, 'data'),
    #'externals': os.path.join(PROJECT_PATH, 'externals'),
    'R': PKGPATH.joinpath('Rscripts'),
    'scripts': PKGPATH.joinpath('scripts'),
    #'package': os.path.join(PROJECT_PATH, 'handygenome'),
}

#DIRS = {
#    'data': os.path.join(PROJECT_PATH, 'data'),
#    'externals': os.path.join(PROJECT_PATH, 'externals'),
#    'R': os.path.join(PROJECT_PATH, 'R'),
#    'scripts': os.path.join(PROJECT_PATH, 'scripts'),
#    'package': os.path.join(PROJECT_PATH, 'handygenome'),
#}



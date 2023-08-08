import sys
import os
import importlib

import yaml

import handygenome


# set numpy threads (https://stackoverflow.com/questions/30791550/limit-number-of-threads-in-numpy)
os.environ["OMP_NUM_THREADS"] = "1"  # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1"  # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "1"  # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"  # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "1"  # export NUMEXPR_NUM_THREADS=6


# constants
PROJECT_PATH = os.path.dirname(os.path.dirname(handygenome.__file__))  
    # /home/users/pjh/scripts/python_genome_package_dev
DIRS = {
    'data': os.path.join(PROJECT_PATH, 'data'),
    'externals': os.path.join(PROJECT_PATH, 'externals'),
    'R': os.path.join(PROJECT_PATH, 'R'),
    'scripts': os.path.join(PROJECT_PATH, 'scripts'),
    'package': os.path.join(PROJECT_PATH, 'handygenome'),
}


# load configuration
CONFIG_PATH = os.path.join(PROJECT_PATH, 'config.yaml')
with open(CONFIG_PATH) as f:
    OPTION = yaml.safe_load(f)
if OPTION['python'] is None:
    OPTION['python'] = sys.executable
    
    

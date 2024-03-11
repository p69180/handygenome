import sys
import os
import argparse
import tempfile

import handygenome.shell_wrapper as shell_wrapper


def argument_parsing(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--conda-prefix', 
        required=False,
        dest='conda_prefix',
        default=None,
        help=f'If not used, the conda environment used for package installation will be used',
    )
    parser.add_argument(
        '--conda-tmpfile-dir', 
        required=False,
        dest='conda_tmpfile_dir',
        default=shell_wrapper.DEFAULT_TMPFILE_DIR,
        help=f'Directory where temporary files will be created. Default is {shell_wrapper.DEFAULT_TMPFILE_DIR}',
    )
    parser.add_argument(
        'cmdargs',
        nargs='+',
        help='Commands that will be executed in a conda environment',
    )

    args = parser.parse_args(args)
    return args


def main():
    args = argument_parsing()
    shell_wrapper.run_with_condaenv(
        args.cmdargs,
        conda_prefix=args.conda_prefix,
        tmpfile_dir=args.conda_tmpfile_dir,
    )
    
    

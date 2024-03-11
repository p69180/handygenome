import sys
import os
import argparse

import handygenome.shell_wrapper as shell_wrapper


def argument_parsing(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--arg1', required=True)
    parser.add_argument('--arg2', required=True)

    args = parser.parse_args(args)
    return args


def main_core(args):
    args = argument_parsing(args)
    print(os.environ['CONDA_PREFIX'])
    print(args)


def main():
    shell_wrapper.run_func_with_script(
        main_core,
        args=[sys.argv[1:]],
        use_condaenv=True,
    )
    
    

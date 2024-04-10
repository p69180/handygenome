import os
import re
import shlex
import itertools

import numpy as np


def handle_cmdargs_arg(cmdargs, to_list=True):
    assert isinstance(cmdargs, (str, list, tuple))
    if to_list:
        if isinstance(cmdargs, (list, tuple)):
            return list(cmdargs)
        elif isinstance(cmdargs, str):
            return shlex.split(cmdargs)
    else:
        if isinstance(cmdargs, (list, tuple)):
            return shlex.join(cmdargs)
        elif isinstance(cmdargs, str):
            return cmdargs
            

def rename_logfile(logpath, success):
    assert logpath.endswith('.log')
    if success:
        os.rename(logpath, re.sub(r'\.log$', '.success', logpath))
    else:
        os.rename(logpath, re.sub(r'\.log$', '.failure', logpath))


def get_shebang(filepath):
    with open(filepath, 'rt') as infile:
        first_line = None
        for line in infile:
            line = line.strip()
            if len(line) > 0:
                first_line = line
                break

    if first_line is None:
        return None
    else:
        if first_line.startswith('#!'):
            return first_line
        else:
            return None
        

def check_has_shebang(filepath):
    return get_shebang(filepath) is not None


def broadcast_args(*args, nargs=None, check_length=False):
    # step1
    newargs = list()
    for x in args:
        if isinstance(x, list):
            newargs.append(x)
        else:
            newargs.append([x])

    # step2
    if nargs is None:
        nargs = max(len(x) for x in newargs)
    if check_length:
        if not all((len(x) in (1, nargs)) for x in newargs):
            raise Exception(f'Length of some of input args are not 1 or nargs.')

    result_args = list()
    for x in newargs:
        rpt_iter = itertools.chain.from_iterable(itertools.repeat(x))
        result_args.append(
            list(itertools.islice(rpt_iter, nargs))
        )

    return result_args
        

def prepare_kwargs(nargs, /, **kwargs):
    assert len(kwargs) > 0
    keys, vals = zip(*kwargs.items())
    bc_vals = broadcast_args(*vals, nargs=nargs, check_length=True)
    result = list()
    for subvals in zip(*bc_vals):
        result.append(dict(zip(keys, subvals)))
    return result



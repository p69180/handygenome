import os
import re
import shlex
import itertools
from collections import UserList

import numpy as np


class MultiArgsList(UserList):
    def __repr__(self):
        return f'{self.__class__.__name__}({super().__repr__()})'


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


def broadcast_args(*args, nargs=None, nargs_idx=None, check_length=False):
    # sanitycheck
    assert sum([(nargs is not None), (nargs_idx is not None)]) <= 1, (
        f'No more than one of "nargs" or "nargs_idx" should be set.'
    )

    # step1: modify input args into MultiArgsList objects
    newargs = list()
    for x in args:
        if isinstance(x, MultiArgsList):
            newargs.append(x)
        else:
            newargs.append(MultiArgsList([x]))

    # step2: set nargs
    if (nargs is None) and (nargs_idx is None):
        set_nargs = max(len(x) for x in newargs)
    else:
        if nargs is not None:
            set_nargs = nargs
        elif nargs_idx is not None:
            set_nargs = len(newargs[nargs_idx])

    # step3: check_length
    if check_length:
        if not all((len(x) in (1, set_nargs)) for x in newargs):
            errmsg = (
                f'Length of some of input args are not 1 or nargs.\n'
                f'args: {args}\n'
                f'nargs: {nargs}\n'
                f'nargs_idx: {nargs_idx}\n'
                f'check_length: {check_length}'
            )
            raise Exception(errmsg)

    # step4: make result
    result_args = list()
    for x in newargs:
        rpt_iter = itertools.chain.from_iterable(itertools.repeat(x))
        result_args.append(
            MultiArgsList(itertools.islice(rpt_iter, set_nargs))
        )

    return result_args
        

def prepare_kwargs(nargs, /, **kwargs):
    assert len(kwargs) > 0
    keys, vals = zip(*kwargs.items())
    bc_vals = broadcast_args(*vals, nargs=nargs, check_length=True)
    result = MultiArgsList()
    for subvals in zip(*bc_vals):
        result.append(dict(zip(keys, subvals)))
    return result



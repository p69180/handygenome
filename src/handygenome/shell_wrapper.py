import sys
import os
import re
import shutil
import subprocess
import tempfile
import textwrap
import shlex
import inspect
import pickle
import importlib
import stat
import builtins
import marshal

import handygenome.tools as tools
import handygenome.workflow as workflow


DEFAULT_TMPFILE_DIR = tempfile.gettempdir()
DEFAULT_PYTHON = sys.executable
CONDAENV_PATSTRING = r'(?P<conda_prefix>(?P<top_prefix>.+)/+envs/+(?P<envname>[^/]+))'
PYTHON_PATSTRING = r'(?P<python>/+bin/+python([0-9]+(\.[0-9]+)?)?)?'
CONDA_PREFIX_PAT = re.compile(CONDAENV_PATSTRING + PYTHON_PATSTRING)


#################
# conda wrapper #
#################

def get_conda_params(conda_prefix=None):
    if conda_prefix is None:
        python_path = sys.executable
        mat = CONDA_PREFIX_PAT.fullmatch(python_path)
        assert mat is not None, f'Invalid python executable path ({repr(python_path)})'
        assert mat.group('python') != '', f'Invalid python executable path'
        conda_prefix = mat.group('conda_prefix')
    else:
        conda_prefix = re.sub(r'/*$', '', conda_prefix)
        mat = CONDA_PREFIX_PAT.fullmatch(conda_prefix)
        assert mat is not None, f'Invalid conda prefix'
        assert mat.group('python') is None, f'Invalid conda prefix'
        python_path = conda_prefix + '/bin/python'
        assert os.path.exists(python_path)

    # get conda executable path
    conda_base_prefix = mat.group('top_prefix')
    conda_envname = mat.group('envname')
    conda_exec = os.path.join(conda_base_prefix, 'bin', 'conda')

    return {
        'prefix': conda_prefix,
        'base_prefix': conda_base_prefix,
        'envname': conda_envname,
        'exec': conda_exec,
        'python_path': python_path,
    }


def run_with_condaenv_prepare(cmd, conda_prefix=None, tmpfile_dir=DEFAULT_TMPFILE_DIR):
    """Args:
        cmd: Can be either str or a sequence.
        conda_prefix: e.g. <leading paths>/miniconda3/envs/<envname>
    """
    bash_path = tools.which('bash')
    env_path = tools.which('env')

    conda_params = get_conda_params(conda_prefix=conda_prefix)
    ld_library_path = os.path.join(conda_params['prefix'], 'lib')
    new_env = {
        'HOME': os.environ['HOME'],
        'LD_LIBRARY_PATH': ld_library_path,
    }

    if isinstance(cmd, (list, tuple)):
        cmd = shlex.join(cmd)

    with tempfile.NamedTemporaryFile(
        mode='wt', 
        prefix='conda_wrapper_script_', 
        suffix='.sh', 
        dir=tmpfile_dir,
        delete=False,
    ) as tmpfile:
        tmpfile.write(
            textwrap.dedent(f'''\
                eval "$({conda_params["exec"]} shell.bash hook)"
                conda activate {conda_params["envname"]}

                set -eu
                {cmd}
            ''')
        )

    run_cmd = [env_path, '-i']
    for key, val in new_env.items():
        run_cmd.append(f'{key}={val}')
    script_path = tmpfile.name
    run_cmd.extend([bash_path, '--noprofile', '--norc', script_path])

    return run_cmd, script_path


def run_with_condaenv(
    cmd, 
    conda_prefix=None, 
    tmpfile_dir=DEFAULT_TMPFILE_DIR, 

    use_ssh=False, 
    hostname=None,

    use_slurm=False,
    slurm_kwargs=dict(),
):
    """Args:
        cmd: Can be either str or a sequence.
        conda_prefix: e.g. <leading paths>/miniconda3/envs/<envname>
    """
    run_cmd, script_path = run_with_condaenv_prepare(
        cmd, 
        conda_prefix=conda_prefix, 
        tmpfile_dir=tmpfile_dir,
    )

    result = workflow.run_subprocess(
        run_cmd,
        use_ssh=use_ssh, 
        hostname=hostname,
        use_slurm=use_slurm,
        slurm_kwargs=slurm_kwargs,
    )

    os.remove(script_path)

    return result
    

##################
# python wrapper #
##################

def check_module_is_loadable(module, python=DEFAULT_PYTHON):
    if module.__spec__ is None:
        return False

    pycmd = textwrap.dedent(
        f'''\
            import importlib.util
            print((importlib.util.find_spec({repr(module.__name__)}) is not None), end="")
        '''
    )
    stdout = subprocess.check_output([python, '-c', pycmd], text=True)
    return eval(stdout)


def make_funcrun_script(
    func, 
    *, 
    args=tuple(), 
    kwargs=dict(), 
    python=DEFAULT_PYTHON, 
    tmpfile_dir=DEFAULT_TMPFILE_DIR,
):
    ##################################
    # make tmpdir and tmp file paths #
    ##################################
    tmpdir = tempfile.mkdtemp(prefix='python_func_running_tmpdir_', dir=tmpfile_dir)
    func_pklpath = os.path.join(tmpdir, 'function.pickle')
    args_pklpath = os.path.join(tmpdir, 'args.pickle')
    kwargs_pklpath = os.path.join(tmpdir, 'kwargs.pickle')
    result_pklpath = os.path.join(tmpdir, 'result.pickle')
    pyscript_path = os.path.join(tmpdir, 'run.py')

    with open(args_pklpath, 'wb') as f:
        pickle.dump(args, f)
    with open(kwargs_pklpath, 'wb') as f:
        pickle.dump(kwargs, f)

    ###################################
    # make python script key contents #
    ###################################
    funcspec = tools.get_callable_spec(func)
    # sanitycheck
    if (not funcspec['isfunction']) and (funcspec['isclsmthd'] is False):  
        # (not funcspec['isclsmthd']) is not appropriate because 
        # funcspec['isclsmthd'] is None in case of builtin constructor methods
        raise Exception(f'If "func" argument is a method, it must be a classmethod.')

    # If function cannot be loaded from a fresh python instance, 
    # make the function into pickle
    module_is_loadable = check_module_is_loadable(funcspec['module'], python=python)
    if not module_is_loadable:
        #raise Exception(f'Only a function defined in a module which can be loaded from a new python instance can be used')
        with open(func_pklpath, 'wb') as outfile:
            marshal.dump(func.__code__, outfile)

    # make func loading lines
    if module_is_loadable:
        funcload_lines = list()
        funcload_lines.append(f'module = importlib.import_module({repr(funcspec["module"].__spec__.name)})')
        if funcspec['isfunction']:
            funcload_lines.append(f'func = getattr(module, {repr(func.__name__)})')
        else:
            funcload_lines.append(f'bound_class = getattr(module, {repr(funcspec["bound_class"].__name__)})')
            funcload_lines.append(f'func = getattr(bound_class, {repr(func.__name__)})')
        funcload_lines = '\n'.join(funcload_lines) + '\n'
    else:
        funcload_lines = textwrap.dedent(f'''\
            with open({repr(func_pklpath)}, 'rb') as infile:
                funccode = marshal.load(infile)
            func = types.FunctionType(funccode, globals(), {repr(func.__name__)}) 
        ''')

    #######################
    # write python script #
    #######################
    with open(pyscript_path, 'wt') as outfile:
        outfile.write(
            textwrap.dedent(f'''\
                #!{python}
                import importlib
                import pickle
                import builtins
                import marshal
                import types

            ''')
        )
        outfile.write(funcload_lines)
        outfile.write(
            textwrap.dedent(f'''\
                with open({repr(args_pklpath)}, 'rb') as f:
                    args = pickle.load(f)
                with open({repr(kwargs_pklpath)}, 'rb') as f:
                    kwargs = pickle.load(f)
                result = func.__call__(*args, **kwargs)
                with open({repr(result_pklpath)}, 'wb') as f:
                    pickle.dump(result, f)
            ''')
        )
    os.chmod(pyscript_path, stat.S_IRUSR | stat.S_IXUSR)

    # result
    result = {
        'topdir': tmpdir,
        'args': args_pklpath,
        'kwargs': kwargs_pklpath,
        'result': result_pklpath,
        'python_script': pyscript_path,
    }
    return result


def run_func_with_script(
    func, 

    *, 

    args=tuple(), 
    kwargs=dict(), 
    python=DEFAULT_PYTHON, 
    python_tmpfile_dir=DEFAULT_TMPFILE_DIR,

    use_condaenv=False, 
    conda_prefix=None, 
    conda_tmpfile_dir=DEFAULT_TMPFILE_DIR, 

    use_ssh=False, 
    hostname=None,

    use_slurm=False,
    slurm_kwargs=dict(),
):
    assert sum([use_ssh, use_slurm]) <= 1, f'No more than one of "use_ssh", "use_slurm" may be set True'

    tmppaths = make_funcrun_script(func, tmpfile_dir=python_tmpfile_dir, python=python, args=args, kwargs=kwargs)
    cmdargs = [python, tmppaths['python_script']]

    if use_condaenv:
        tmpresult = run_with_condaenv(
            cmdargs, 
            conda_prefix=conda_prefix,
            tmpfile_dir=conda_tmpfile_dir,
            use_ssh=use_ssh, 
            hostname=hostname,
            use_slurm=use_slurm,
            slurm_kwargs=slurm_kwargs,
        )
    else:
        tmpresult = workflow.run_subprocess(
            cmdargs, 
            use_ssh=use_ssh, 
            hostname=hostname,
            use_slurm=use_slurm,
            slurm_kwargs=slurm_kwargs,
        )

    with open(tmppaths['result'], 'rb') as infile:
        result = pickle.load(infile)
    shutil.rmtree(tmppaths['topdir'])

    return result



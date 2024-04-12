import sys
import os
import re
import shutil
import subprocess
import multiprocessing
import tempfile
import textwrap
import shlex
import inspect
import pickle
import stat
import builtins
import marshal

import handygenome.tools as tools
import handygenome.deco as deco
#import handygenome.workflow as workflow
import handygenome.utils.workflow_utils as worfklow_utils
from handygenome.utils.workflow_utils import MultiArgsList
import handygenome.workflow.slurm as libslurm


#DEFAULT_TMPFILE_DIR = tempfile.gettempdir()
DEFAULT_TMPFILE_DIR = os.getcwd()
DEFAULT_PYTHON = sys.executable
CONDAENV_PATSTRING = r'(?P<conda_prefix>(?P<top_prefix>.+)/+envs/+(?P<envname>[^/]+))'
PYTHON_PATSTRING = r'(?P<python>/+bin/+python([0-9]+(\.[0-9]+)?)?)?'
CONDA_PREFIX_PAT = re.compile(CONDAENV_PATSTRING + PYTHON_PATSTRING)


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


##################################
# running command-line arguments #
##################################

def run_cmdargs_targetfunc(cmdargs_item, logpath_item):
    run_kwargs = dict(
        check=False,
        shell=True,
    )
    if logpath_item is not None:
        log = open(logpath_item, 'wt')
        run_kwargs.update(text=True, stdout=log, stderr=log)

    p = subprocess.run(cmdargs_item, **run_kwargs)

    if logpath_item is not None:
        log.close()
        worfklow_utils.rename_logfile(logpath_item, (p.returncode == 0))

    return p


@deco.get_deco_broadcast(
    ['cmdargs', 'logpath', 'hostname'],
    nargs_name='cmdargs',
    check_length=True,
)
def run_cmdargs_base(
    cmdargs, 
    logpath=None,
    max_run=None, # only works with non-slurm mode

    use_ssh=False, 
    hostname=None,
    set_wd=True,

    use_slurm=False,
    slurm_kwargs=dict(),  # nproc, jobname, intv_submit, intv_check, max_submit, etc.
):
    """- Run command-line arguments within a subshell, by either slurm,
        on another host over ssh, or on the current host.
    - If logpath is None, error message will be printed to stderr.
        Else, error messages will be written to the log file.

    Args:
        cmdargs: A single str object or a MultiArgsList of str objects
        set_wd: Only relevant when "use_ssh" is True. If "set_wd" is True, 
            move to current working directory before running commands
    """

    # sanity check
    assert sum([use_ssh, use_slurm]) <= 1, (
        f'No more than one of "use_ssh", "use_slurm" may be set True'
    )
    if any(
        ((x is not None) and (not x.endswith('.log')))
        for x in logpath
    ):
        raise Exception(f'All non-None "logpath" must end with ".log".')
    if use_ssh and any((x is None) for x in hostname):
        raise Exception(f'When "use_ssh" is True, "hostname" must be set.')

    # main
    if use_slurm:
        # slurm automatically sets wd
        success_list, exitcode_list, script_path_list = libslurm.run_cmdargs_with_slurm(
            cmdargs, 
            logpath=logpath,
            **slurm_kwargs,
        )
        subp_result = {
            'success': success_list,
            'exitcode': exitcode_list,
            'jobscript_path': script_path_list,
            'subprocess_obj': None,
        }
    else:
        if use_ssh:
            ssh_path = tools.which('ssh')
            new_cmdargs = MultiArgsList()
            for cmdargs_item, hostname_item in zip(cmdargs, hostname):
                if set_wd:
                    new_item = f'{ssh_path} {hostname_item} "cd {os.getcwd()} ; {cmdargs_item}"'
                else:
                    new_item = f'{ssh_path} {hostname_item} "{cmdargs_item}"'
                new_cmdargs.append(new_item)
            cmdargs = new_cmdargs

        with multiprocessing.Pool(max_run) as pool:
            plist = pool.starmap(
                run_cmdargs_targetfunc, 
                zip(cmdargs, logpath),
            )

        subp_result = {
            'success': [(x.returncode == 0) for x in plist],
            'exitcode': [x.returncode for x in plist],
            'jobscript_path': None,
            'subprocess_obj': plist,
        }

    subp_result['all_success'] = all(subp_result['success'])

    return subp_result


def prepare_conda_wrapping(
    cmdargs_item, 
    conda_prefix=None, 
    condawrap_script_dir=DEFAULT_TMPFILE_DIR,
    condawrap_script_prefix='conda_wrapping_script_',
    shell="bash",
):
    """Args:
        cmdargs_item: Must be a str object.
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

    with tempfile.NamedTemporaryFile(
        mode='wt', 
        prefix=condawrap_script_prefix, 
        suffix='.sh', 
        dir=condawrap_script_dir,
        delete=False,
    ) as tmpfile:
        tmpfile.write(
            textwrap.dedent(f'''\
                eval "$({conda_params["exec"]} shell.{shell} hook)"
                conda activate {conda_params["envname"]}

                set -eu
                {cmdargs_item}
            ''')
        )

    new_cmdargs = [env_path, '-i']
    for key, val in new_env.items():
        new_cmdargs.append(f'{key}={val}')
    script_path = tmpfile.name
    new_cmdargs.extend([bash_path, '--noprofile', '--norc', script_path])

    new_cmdargs = shlex.join(new_cmdargs)

    return new_cmdargs, script_path


@deco.get_deco_broadcast(
    ['cmdargs', 'logpath', 'hostname'],
    nargs_name='cmdargs',
    check_length=True,
)
def run_cmdargs(
    cmdargs, 
    logpath=None,
    max_run=None, # only works with non-slurm mode

    use_condaenv=False,
    conda_prefix=None, 
    condawrap_script_dir=DEFAULT_TMPFILE_DIR, 
    remove_condawrap_script_dir=True,

    use_ssh=False, 
    hostname=None,

    use_slurm=False,
    slurm_kwargs=dict(),

    raise_with_failure=True,
):
    """
    ####################
    # For end-user use #
    ####################
    Args:
        cmdargs: Must be str object or a list of them
        conda_prefix: e.g. <leading paths>/miniconda3/envs/<envname>

    Automatically remove conda-wrapping-related scripts
    """
    if use_condaenv:
        new_cmdargs = MultiArgsList()
        condawrap_script_list = list()
        for zidx, cmdargs_item in tools.zenumerate(cmdargs):
            new_cmdargs_item, condawrap_script = prepare_conda_wrapping(
                cmdargs_item, 
                conda_prefix=conda_prefix, 
                condawrap_script_dir=condawrap_script_dir,
                condawrap_script_prefix=f'conda_wrapping_script_{zidx}_',
            )
            new_cmdargs.append(new_cmdargs_item)
            condawrap_script_list.append(condawrap_script)
        cmdargs = new_cmdargs

    subp_result = run_cmdargs_base(
        cmdargs,
        logpath=logpath,
        max_run=max_run,
        use_ssh=use_ssh, 
        hostname=hostname,
        use_slurm=use_slurm,
        slurm_kwargs=slurm_kwargs,
    )
    if use_condaenv and remove_condawrap_script_dir:
        for x in condawrap_script_list:
            os.remove(x)

    if raise_with_failure and (not subp_result['all_success']):
        logpath_given = any((x is not None) for x in logpath)
        if logpath_given:
            logpath_string = '\n'.join(logpath)
            raise Exception(
                f'One or more of the jobs finished unsuccessfully. Refer to '
                f'the log files.\nLog file paths:\n{logpath_string}'
            )
        else:
            raise Exception(
                f'One or more of the jobs finished unsuccessfully. Log files '
                f'were not created because "logpath" argument was not set.'
            )

    return subp_result


####################
# running a script #
####################

@deco.get_deco_broadcast(
    ['script_path', 'logpath', 'hostname'],
    nargs_name='cmdargs',
    check_length=True,
)
def run_script(
    script_path, 
    logpath=None,
    max_run=None, # only works with non-slurm mode

    use_condaenv=False,
    conda_prefix=None, 
    condawrap_script_dir=DEFAULT_TMPFILE_DIR, 
    remove_condawrap_script_dir=True,

    use_ssh=False,
    hostname=None,

    use_slurm=False,
    slurm_kwargs=dict(),

    raise_with_failure=True,
):
    """
    ####################
    # For end-user use #
    ####################
    """
    if any((not os.access(x, os.R_OK | os.X_OK)) for x in script_path):
        raise Exception(f'read/execute permission is not set for the script file.')
    if any((not workflow_utils.check_has_shebang(script_path)) for x in script_path):
        raise Exception(f'The script file does not contain a shebang line.')

    subp_result = run_cmdargs(
        cmdargs=script_path, 
        logpath=logpath,
        max_run=max_run,

        use_condaenv=use_condaenv,
        conda_prefix=conda_prefix,
        condawrap_script_dir=condawrap_script_dir,
        remove_condawrap_script_dir=remove_condawrap_script_dir,

        use_ssh=use_ssh,
        hostname=hostname,

        use_slurm=use_slurm,
        slurm_kwargs=slurm_kwargs,

        raise_with_failure=raise_with_failure,
    )
    return subp_result
    

#############################
# running a python function #
#############################

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


def prepare_funcrun_write_pyscript(
    func, 
    python, 
    func_pklpath,
    pyscript_path,
    args_pklpath,
    kwargs_pklpath,
    result_pklpath,
):
    ###################################
    # make python script key contents #
    ###################################

    funcspec = tools.get_callable_spec(func)

    # sanitycheck
    if (not funcspec['isfunction']) and (funcspec['isclsmthd'] is False):  
        '''
        (not funcspec['isclsmthd']) is not appropriate because 
        funcspec['isclsmthd'] is None in case of builtin constructor methods
        '''
        raise Exception(f'If "func" argument is a method, it must be a classmethod.')

    '''
    If function cannot be loaded from a fresh python instance, 
    make the function into pickle
    '''
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
                    # If SystemExit is raised from this call, "result"
                    # is not created and below code is not executed.
                with open({repr(result_pklpath)}, 'wb') as f:
                    pickle.dump(result, f)
            ''')
        )

    os.chmod(pyscript_path, stat.S_IRUSR | stat.S_IXUSR)


#@deco.get_deco_broadcast(['func', 'args', 'kwargs'])
def prepare_funcrun(
    func, 
    *, 
    args=tuple(), 
    kwargs=dict(), 
    python=DEFAULT_PYTHON, 
    funcrun_tmpfile_dir=DEFAULT_TMPFILE_DIR,
):
    """Args:
        func: Can be a function object or a list of functions.
        args: If a tuple, interpreted as positional arguments for a single 
            instance. If a list, each element must be a tuple, representing 
            a single instance.
        kwargs: Similar to args but dict rather than tuple represents a single
            instance.
    """
    ###############
    # sanitycheck #
    ###############


    assert all(callable(x) for x in func)
    assert all(isinstance(x, tuple) for x in args)
    assert all(isinstance(x, dict) for x in kwargs)

    ##################################
    # make tmpdir and tmp file paths #
    ##################################

    tmpdir = tempfile.mkdtemp(
        prefix='python_func_running_tmpdir_', 
        dir=funcrun_tmpfile_dir,
    )
    tmpdir_tree = dict()
    tmpdir_tree['top'] = tmpdir
    tmpdir_tree['subdirs'] = dict()
    
    for zidx, (one_func, one_args, one_kwargs) in tools.zenumerate(
        zip(func, args, kwargs)
    ):
        subdir = os.path.join(tmpdir_tree['top'], zidx)
        os.mkdir(subdir)

        subdir_tree = dict()
        subdir_tree['top'] = subdir
        subdir_tree['func'] = os.path.join(subdir, 'function.pickle')
        subdir_tree['args'] = os.path.join(subdir, 'args.pickle')
        subdir_tree['kwargs'] = os.path.join(subdir, 'kwargs.pickle')
        subdir_tree['result'] = os.path.join(subdir, 'result.pickle')
        subdir_tree['pyscript'] = os.path.join(subdir, 'run.py')

        with open(subdir_tree['args'], 'wb') as f:
            pickle.dump(one_args, f)
        with open(subdir_tree['kwargs'], 'wb') as f:
            pickle.dump(one_kwargs, f)

        prepare_funcrun_write_pyscript(
            func=one_func, 
            python=python, 
            func_pklpath=subdir_tree['func'],
            pyscript_path=subdir_tree['pyscript'],
            args_pklpath=subdir_tree['args'],
            kwargs_pklpath=subdir_tree['kwargs'],
            result_pklpath=subdir_tree['result'],
        )

        tmpdir_tree['subdirs'][zidx] = subdir_tree

    return tmpdir_tree


@deco.get_deco_broadcast(
    ['func', 'args', 'kwargs', 'logpath', 'hostname'],
    nargs_name='func',
    check_length=True,
)
@deco.get_deco_num_set_differently(['use_ssh', 'use_slurm'], 1, how='le')
def run_func(
    func, 

    *, 

    args=tuple(), 
    kwargs=dict(), 
    logpath=None,
    max_run=None, # only works with non-slurm mode

    raise_with_failure=True,

    ###

    use_condaenv=False,
    conda_prefix=None, 
    condawrap_script_dir=DEFAULT_TMPFILE_DIR, 
    remove_condawrap_script_dir=True, 

    use_ssh=False, 
    hostname=None,

    use_slurm=False,
    slurm_kwargs=dict(),

    ### 

    python=None,
    funcrun_tmpfile_dir=DEFAULT_TMPFILE_DIR,
    remove_funcrun_tmpfile_dir=True,
):
    """
    ####################
    # For end-user use #
    ####################
    """
    # set python executable
    if python is None:
        if use_condaenv:
            conda_params = get_conda_params(conda_prefix=conda_prefix)
            python = conda_params['python_path']
        else:
            python = DEFAULT_PYTHON

    # make script
    tmpdir_tree = prepare_funcrun(
        func, 
        args=args, 
        kwargs=kwargs,
        python=python, 
        funcrun_tmpfile_dir=funcrun_tmpfile_dir, 
    )
    cmdargs = MultiArgsList(
        subdir_tree['pyscript']
        for subdir_tree in tmpdir_tree['subdirs'].values()
    )

    # run cmdargs
    subp_result = run_cmdargs(
        cmdargs=cmdargs,
        logpath=logpath,
        max_run=max_run,
        use_condaenv=use_condaenv,
        conda_prefix=conda_prefix,
        condawrap_script_dir=condawrap_script_dir,
        remove_condawrap_script_dir=remove_condawrap_script_dir,
        use_ssh=use_ssh,
        hostname=hostname,
        use_slurm=use_slurm,
        slurm_kwargs=slurm_kwargs,
        raise_with_failure=raise_with_failure,
    )

    # obtain result object  
    result = list()
    for subdir_tree in tmpdir_tree['subdirs'].values():
        if os.path.exists(subdir_tree['result']):
            with open(subdir_tree['result'], 'rb') as infile:
                result_item = pickle.load(infile)
        else:
            '''
            When the function "func" raises SystemExit, for example when
            argparse is given "-h" and help message is printed,
            tmppaths['result'] is not created.
            '''
            result_item = None
        result.append(result_item)

    # remove funcrun temporary directory
    if remove_funcrun_tmpfile_dir:
        shutil.rmtree(tmpdir_tree['top'])

    return result, subp_result



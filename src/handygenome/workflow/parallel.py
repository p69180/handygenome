import os
import subprocess

import handygenome.deco as deco
import handygenome.workflow.shell_wrapper as shell_wrapper
import handygenome.workflow.slurm as libslurm
from handygenome.workflow.slurm import JobList


############################
# running multiple scripts #
############################

def run_scripts_local(script_path_list, log_paths=None):
    if not all(os.access(x, os.R_OK | os.X_OK) for x in script_path_list):
        raise Exception(f'read/execute permission is not set for job script files')
    if log_paths is not None:
        assert len(log_paths) == len(script_path_list)

    plist = list()
    if log_paths is not None:
        log_files = [open(x, 'wt') for x in log_paths]
        for jobscript_path, logfile in zip(script_path_list, log_files):
            p = subprocess.Popen([jobscript_path], stdout=logfile, stderr=logfile, text=True)
            plist.append(p)
    else:
        for jobscript_path in script_path_list:
            p = subprocess.Popen([jobscript_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            plist.append(p)

    for p in plist:
        p.wait()

    if log_paths is not None:
        for x in log_files:
            x.close()

    exitcode_list = [p.returncode for p in plist]
    success = all(x == 0 for x in exitcode_list)

    return success, exitcode_list


@deco.get_deco_arg_choices({'sched': ['local', 'slurm']})
def run_scripts(
    jobscript_paths, 
    sched, 
    log_dir=None, 
    raise_on_failure=True,
    log_paths=None,
    **kwargs,

    #intv_check, 
    #intv_submit, 
    #max_submit, 
    #logger, 
    #job_status_logpath=None,
    #title=None,
):
    if log_paths is not None:
        assert all(x.endswith('.log') for x in log_paths)

    if sched == 'local':
        success, exitcode_list = run_scripts_local(
            jobscript_paths, log_paths=log_paths,
        )
    elif sched == 'slurm':
        joblist = JobList(
            jobscript_paths, 
            **kwargs,
        )
        joblist.run()

        success = joblist.success
        exitcode_list = joblist.get_exitcodes()

    # rename suffixes of log files
    if log_paths is not None:
        for code, filepath in zip(exitcode_list, log_paths):
            if code == 0:
                os.rename(filepath, re.sub(r'\.log$', '.success', filepath))
            else:
                os.rename(filepath, re.sub(r'\.log$', '.failure', filepath))

    if raise_on_failure:
        if not success:
            #raise SystemExit(
            #    (f'One or more jobs have finished unsuccessfully. Refer to '
            #     f'log files in {log_dir}'))
            #raise SystemExit(f'One or more jobs have finished unsuccessfully.')
            raise Exception(f'One or more jobs have finished unsuccessfully.')
        #else:
        #    return None
    #else:

    return success, exitcode_list


##########################################################
# running a function with multiple argument combinations #
##########################################################



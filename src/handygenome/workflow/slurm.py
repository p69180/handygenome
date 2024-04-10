import os
import re
import subprocess
import time
import textwrap
import itertools
import tempfile
import shutil

import numpy as np

import handygenome.deco as deco
import handygenome.logutils as logutils
import handygenome.tools as tools
import handygenome.workflow as workflow
import handygenome.utils.workflow_utils as workflow_utils


SLURMBIN = '/usr/local/slurm/bin'
SQUEUE = os.path.join(SLURMBIN, 'squeue')
SBATCH = os.path.join(SLURMBIN, 'sbatch')
SCONTROL = os.path.join(SLURMBIN, 'scontrol')
SCANCEL = os.path.join(SLURMBIN, 'scancel')
SINFO = os.path.join(SLURMBIN, 'sinfo')
assert all(
    os.path.exists(x) for x in [
        SQUEUE, SBATCH,SCONTROL, SCANCEL, SINFO,
    ]
)

#SQUEUE = tools.which('squeue')
#SBATCH = tools.which('sbatch')
#SCONTROL = tools.which('scontrol')
#SCANCEL = tools.which('scancel')
#SINFO = tools.which('sinfo')


class Job:
    """Designed for slurm

    Attributes:
        jobscript_path: Path of job script file to run.
        jobscript_string: Multi-line string to be written to stdin of sbatch.
        jobid: JobId as integer. None before job submission.
        submitted: True if job has been submitted, False otherwise.
        jobstate
        exitcode: type int
        success: 
            None: 
                1) before job has been finished
                2) job has been finished but its state cannot be checked 
                    because job id has disappeared from slurm
            True if the job has finished successfully
            False if the job has finished not successfully
        no_jobid: 
            None: before job has been submitted
            True: job id is no longer remaining on slurm
            False: job id exists on slurm
        scontrol_result: Dictionary which contains result of 
            "scontrol show job <jobid>" command.
        sbatch_err_info: An informative string created when sbatch returns 
            a nonzero exit code.
        stdout_path: job stdout path
        stderr_path: job stderr path
        status
    """
    jobstates_pending = [
        'CONFIGURING', 'PENDING',
    ]
    jobstates_running = [
        'COMPLETING', 'RUNNING', 'RESV_DEL_HOLD', 
        'REQUEUE_FED', 'REQUEUE_HOLD', 'RESIZING', 
        'SIGNALING', 'SPECIAL_EXIT', 'STAGE_OUT', 'STOPPED', 
        'SUSPENDED',
    ]
    jobstates_finished = [
        'CANCELLED', 'COMPLETED', 'BOOT_FAIL', 'DEADLINE', 
        'FAILED', 'NODE_FAIL', 'OUT_OF_MEMORY', 'PREEMPTED', 
        'TIMEOUT',
    ]
    jobstates_unknown = [
        'REVOKED',
    ]

    @deco.get_deco_num_set_differently(('jobscript_path', 'jobscript_string'), 1)
    def __init__(
        self, 
        jobscript_path=None, 
        jobscript_string=None, 
        verbose=True, 
    ):
        self.jobscript_path = jobscript_path
        self.jobscript_string = jobscript_string
        self.jobid = None
        self.submitted = False
        self.jobstate = None
        self.exitcode = None
        self.success = None
        self.no_jobid = None
        self.scontrol_result = None
        self.sbatch_err_info = None
        self.stdout_path = None
        self.stderr_path = None
        self.status = {
            'pending': False, 
            'running': False, 
            'finished': False,
        }

    def __repr__(self):
        infostr = ', '.join(
            f'{key}: {getattr(self, key)}'
            for key in ('jobid', 'jobstate', 'exitcode', 'jobscript_path')
        )
        return f'<Job ({infostr})>'

    def submit(self):
        if self.jobscript_path is None:
            self._submit_string()
        else:
            self._submit_path()

    def cancel(self):
        assert self.submitted, f'The job is not yet submitted.'
        p = subprocess.run([SCANCEL, str(self.jobid)],
                           capture_output=True, text=True)
        if p.returncode != 0:
            self.scancel_err_info = textwrap.dedent(f"""\
                Job cancellation using scancel failed.

                job id: {self.jobid}
                scancel stdout: {p.stdout}
                scancel stderr: {p.stderr}
                scancel exit code: {p.returncode}""")
            raise Exception(self.scancel_err_info)
        else:
            self.update()
            logutils.log(f'Cancelled a job: JobID - {self.jobid}', level='info')

    def set_status(self, key):
        self.status[key] = True
        other_keys = set(self.status.keys()).difference([key])
        for x in other_keys:
            self.status[x] = False

    def update(self):
        assert self.submitted

        scontrol_result, no_jobid = get_scontrol_job_result(self.jobid)[:2]
        self.no_jobid = no_jobid

        if self.no_jobid:
            # self.success and self.scontrol_result are not touched
            self.set_status('finished')
        else:
            self.scontrol_result = scontrol_result

            # update self.stderr/stdout path
            self.stdout_path = self.scontrol_result['StdOut']
            self.stderr_path = self.scontrol_result['StdErr']

            # update self.jobstate and self.status
            self.jobstate = self.scontrol_result['JobState']
            if self.jobstate in self.__class__.jobstates_pending:
                self.set_status('pending')
            elif self.jobstate in self.__class__.jobstates_running:
                self.set_status('running')
            elif self.jobstate in self.__class__.jobstates_finished:
                self.set_status('finished')
            else:
                raise Exception(
                    textwrap.dedent(f"""\
                        Unable to interpret JobState.
                        JobState: {self.jobstate}
                        scontrol result: {self.scontrol_result}"""
                    )
                )

            # update self.exitcode, self.success
            if self.status['finished']:
                self.success = (self.jobstate == 'COMPLETED')
                self.exitcode = int(scontrol_result['ExitCode'].split(':')[0])

    ##########################

    def _submit_string(self):
        p = subprocess.run([SBATCH], input=self.jobscript_string,
                           capture_output=True, text=True)

        if p.returncode != 0:
            self.sbatch_err_info = textwrap.dedent(f"""\
                Job submission using sbatch failed.

                string written to the stdin of sbatch: {self.jobscript_string}
                sbatch stdout: {p.stdout}
                sbatch stderr: {p.stderr}
                sbatch exit code: {p.returncode}"""
            )
            raise Exception(self.sbatch_err_info)
        else:
            self.jobid = int(p.stdout.split()[-1])
            self.submitted = True
            logutils.log(f'Submitted a job: JobID - {self.jobid}', level='info')

    def _submit_path(self):
        p = subprocess.run([SBATCH, self.jobscript_path],
                           capture_output=True, text=True)

        if p.returncode != 0:
            self.sbatch_err_info = textwrap.dedent(f"""\
                Job submission using sbatch failed.

                job script path: {self.jobscript_path}
                sbatch stdout: {p.stdout}
                sbatch stderr: {p.stderr}
                sbatch exit code: {p.returncode}"""
            )
            raise Exception(self.sbatch_err_info)
        else:
            self.jobid = int(p.stdout.split()[-1])
            self.submitted = True
            logutils.log(f'Submitted a job: JobID - {self.jobid}', level='info')


class JobList(list):
    """
    Attributes:
        success
    """
    def __init__(
        self, 
        jobscript_path_list, 
        intv_submit=workflow.DEFAULT_INTV_SUBMIT,
        intv_check=workflow.DEFAULT_INTV_CHECK,
        max_submit=workflow.DEFAULT_MAX_SUBMIT,
        verbose=True, 
        #logpath=None, 
        job_status_logpath=None,
        title=None,
    ):
        self.verbose = verbose

        # job status summary logger
        self.job_status_logpath = job_status_logpath

        # other attrs
        for jobscript_path in jobscript_path_list:
            self.append(
                Job(jobscript_path=jobscript_path, verbose=verbose)
            )

        self.title = title
        self.intv_submit = intv_submit
        self.intv_check = intv_check
        self.max_submit = max_submit
        self.success = None
        self.sublists = dict()
        self.update()

    def write_job_status_log(self, msg):
        if self.job_status_logpath is not None:
            #with open(self.job_status_logpath, 'wt') as f:
            with tools.openfile(self.job_status_logpath, 'a') as f:
                f.write(f'[{logutils.get_timestamp()}] {msg}')

    def get_num_pending_running(self):
        return len(self.sublists['pending']) + len(self.sublists['running'])

    def submit(self):
        if self.max_submit is None:
            num_to_submit = len(self.sublists['notsubmit'])
        else:
            num_pending_running = self.get_num_pending_running()
            num_to_submit = min(
                self.max_submit - num_pending_running, 
                len(self.sublists['notsubmit']),
            )

        if num_to_submit > 0:
            for job in self.sublists['notsubmit'][:num_to_submit]:
                job.submit()
                time.sleep(self.intv_submit)
            self.update()

    def cancel_all_running(self):
        for job in self.sublists['running']:
            job.cancel()


    def make_infostring(self, sublist_key):
        n_jobs = len(self.sublists[sublist_key])
        details = ', '.join(
            job.__repr__()  # 'jobid', 'jobstate', 'exitcode'
            for job in self.sublists[sublist_key]
        )
        return f'{n_jobs} ({details})'

    def log_status(self):
        n_notsubmit = len(self.sublists['notsubmit'])
        n_pending = len(self.sublists['pending'])
        n_running = len(self.sublists['running'])
        n_finished = len(self.sublists['finished'])

        info_pending = self.make_infostring('pending')
        info_running = self.make_infostring('running')
        info_finished = self.make_infostring('finished')

        # more verbose information
        self.write_job_status_log(
            textwrap.dedent(
                f"""\
                Current job status (title: {self.title}):
                    Not submitted yet: {n_notsubmit}

                    Pending: {info_pending}

                    Running: {info_running}

                    Finished: {info_finished}
                """
            )
        )

        # concise information (to be printed to stderr)
        logutils.log(
            textwrap.dedent(
                f"""\
                Current job status (title: {self.title}):
                    Not submitted yet: {n_notsubmit}
                    Pending: {n_pending}
                    Running: {n_running}
                    Finished: {n_finished}"""
            ),
            level='info',
        )

    def log_epilogue(self):
        n_success = len(self.sublists['success'])
        n_failure = len(self.sublists['failure'])
        n_unknown = len(self.sublists['unknown'])

        info_success = self.make_infostring('success')
        info_failure = self.make_infostring('failure')
        info_unknown = self.make_infostring('unknown')

        self.write_job_status_log(
            textwrap.dedent(
                f"""\
                All finished.
                    Successful jobs: {info_success}

                    Failed jobs: {info_failure}

                    Jobs with unknown exit statuses: {info_unknown}
                """
            )
        )
        logutils.log(
            textwrap.dedent(
                f"""\
                All finished.
                    Successful jobs: {n_success}
                    Failed jobs: {n_failure}
                    Jobs with unknown exit statuses: {n_unknown}"""
            ),
            level='info',
        )

    def submit_and_wait(self):
        # main
        try:
            while True:
                self.update()
                self.submit()
                self.log_status()
                if all(job.status['finished'] for job in self):
                    break
                else:
                    time.sleep(self.intv_check)
                    continue
        except KeyboardInterrupt:
            logutils.log(
                (
                    f'RECEIVED A KEYBOARD INTERRUPT; '
                    f'will cancel all pending and running jobs with scancel,'
                    f' then exit immediately.'
                ),
                level='info',
            )
            for job in itertools.chain(
                self.sublists['pending'], self.sublists['running']
            ):
                job.cancel()
            raise SystemExit(1)
        else:
            self.set_sublists()
            self.log_epilogue()
            self.success = all(job.success for job in self)

    def recover_cancelled_jobs(self):
        cancelled_jobs = [
            (idx, job) for (idx, job) in enumerate(self) if job.jobstate == 'CANCELLED'
        ]
        for idx, job in cancelled_jobs:
            del self[idx]

        for idx, job in cancelled_jobs:
            self.append(
                Job(
                    jobscript_path=job.jobscript_path, 
                    verbose=self.verbose, 
                )
            )

    def run(self):
        while True:
            self.submit_and_wait()
            if any((job.jobstate == 'CANCELLED') for job in self):
                self.recover_cancelled_jobs()
                continue
            else:
                break

    def get_exitcodes(self):
        return [job.exitcode for job in self]

    def update(self):
        self.update_jobs()
        self.set_sublists()

    def update_jobs(self):
        for job in self:
            if job.submitted:
                job.update()

    def set_sublists(self):
        for key in (
            'notsubmit', 
            'pending', 
            'running', 
            'finished', 
            'success', 
            'failure', 
            'unknown',
        ):
            self.sublists[key] = list()

        for job in self:
            if not job.submitted:
                self.sublists['notsubmit'].append(job)
            if job.status['pending']:
                self.sublists['pending'].append(job)
            if job.status['running']:
                self.sublists['running'].append(job)
            if job.status['finished']:
                self.sublists['finished'].append(job)

            if job.success is True:
                self.sublists['success'].append(job)
            if job.success is False:
                self.sublists['failure'].append(job)
            if job.success is None:
                self.sublists['unknown'].append(job)

#    def set_sublists_statuses(self):
#        self.sublists['notsubmit'] = [job for job in self 
#                                      if not job.submitted]
#        self.sublists['pending'] = [job for job in self 
#                                    if job.status['pending']]
#        self.sublists['running'] = [job for job in self 
#                                    if job.status['running']]
#        self.sublists['finished'] = [job for job in self 
#                                     if job.status['finished']]
#
#    def set_sublists_exitcodes(self):
#        self.sublists['success'] = [job for job in self 
#                                    if job.success is True]
#        self.sublists['failure'] = [job for job in self 
#                                    if job.success is False]
#        self.sublists['unknown'] = [job for job in self 
#                                    if job.success is None]


def get_scontrol_job_result(jobid):
    """
    Returns:
        A tuple (scontrol_result, no_jobid, returncode, stderr)

        scontrol_result: A dict with keys: 'JobId', 'JobName', 'UserId', 
            'GroupId', 'MCS_label', 'Priority', 'Nice', 'Account', 'QOS', 
            'JobState', 'Reason', 'Dependency', 'Requeue', 'Restarts', 
            'BatchFlag', 'Reboot', 'ExitCode', 'RunTime', 'TimeLimit', 
            'TimeMin', 'SubmitTime', 'EligibleTime', 'AccrueTime', 
            'StartTime', 'EndTime', 'Deadline', 'SuspendTime', 
            'SecsPreSuspend', 'LastSchedEval', 'Partition', 'AllocNode:Sid', 
            'ReqNodeList', 'ExcNodeList', 'NodeList', 'BatchHost', 'NumNodes', 
            'NumCPUs', 'NumTasks', 'CPUs/Task', 'ReqB:S:C:T', 'TRES', 
            'Socks/Node', 'NtasksPerN:B:S:C', 'CoreSpec', 'MinCPUsNode', 
            'MinMemoryNode', 'MinTmpDiskNode', 'Features', 'DelayBoot', 
            'OverSubscribe', 'Contiguous', 'Licenses', 'Network', 'Command', 
            'WorkDir', 'StdErr', 'StdIn', 'StdOut', 'Power', 'NtasksPerTRES:0'
        no_jobid: True or False 
        returncode
        stderr
    """

    cmdargs = [ SCONTROL, '-o', 'show', 'job', str(jobid) ]
    p = subprocess.run(args=cmdargs, text=True, capture_output=True)

    returncode = p.returncode
    stderr = p.stderr

    if p.returncode != 0:
        if p.stderr.strip() == 'slurm_load_jobs error: Invalid job id specified':
            scontrol_result = None
            no_jobid = True
        else:
            e_msg = textwrap.dedent(f'''\
                scontrol show job finished with nonzero exit code.
                commands: {cmdargs}
                stdout: {p.stdout}
                stderr: {p.stderr}
                exit code: {p.returncode}'''
            )
            raise Exception(e_msg)

    else:
        no_jobid = False

        scontrol_result = dict()
        for word in p.stdout.split():
            wordsp = word.split('=', maxsplit = 1)
            if len(wordsp) == 1:
                key = wordsp[0]
                val = None
            elif len(wordsp) == 2:
                key, val = wordsp
                """
            else:
                e_msg = f'''\
Field with more than 2 "=" character found from "scontrol show job" output.
Field: {word}
scontrol command: {cmdargs}'''
                raise Exception(e_msg)
                """
    
            scontrol_result[key] = val
    
    return scontrol_result, no_jobid, returncode, stderr


#########################
# JOB SCRIPT GENERATION #
#########################

def make_multiline_command(lines, leading_taps = 0):
    new_lines = list()
    new_lines.append( '\t'*leading_taps + lines[0] + ' \\' )
    new_lines.extend( [ '\t'*(leading_taps+1) + x + ' \\' for x in lines[1:-1] ] )
    new_lines.append( '\t'*(leading_taps+1) + lines[-1] )

    return '\n'.join(new_lines)


def make_jobscript_string(lines, shell = False, python = False, **kwargs):
    string_list = list()


    string_list.append('#!' + tools.which('bash'))

    if 'N' not in kwargs:
        kwargs['N'] = 1
    if 'n' not in kwargs:
        kwargs['n'] = 1
    if 'o' not in kwargs:
        kwargs['o'] = '/dev/null'
    if 'c' not in kwargs:
        kwargs['c'] = 1

    for k,v in kwargs.items():
        if v is not None:
            string_list.append(f'#SBATCH -{k} {v}')

    string_list.append('')

    for line in lines:
        string_list.append(line)

    return '\n'.join(string_list)


def make_jobscript(
    jobscript_path, 
    lines, 
    shebang=None, 
    nproc=1,
    stdout=None,
    stderr=None,
    jobname=None,
):
    workflow.check_outfile_validity(jobscript_path)
    if stdout is None:
        stdout = os.devnull
    if stderr is None:
        stderr = os.devnull

    contents = list()

    if shebang is None:
        shebang = '#!' + tools.which('bash')
    contents.append(shebang)
    contents.append(f'#SBATCH -N 1 -n 1')
    contents.append(f'#SBATCH -c {nproc}')
    contents.append(f'#SBATCH -o {stdout}')
    contents.append(f'#SBATCH -e {stderr}')
    if jobname is not None:
        contents.append(f'#SBATCH -J {jobname}')

    contents.append('')

    contents.extend(lines)

    with open(jobscript_path, 'wt') as outfile:
        outfile.write('\n'.join(contents))


########################################
# running command-line args with slurm #
########################################

@deco.get_deco_broadcast(
    ['cmdargs', 'logpath', 'nproc', 'jobname']
)
def run_cmdargs_with_slurm(
    cmdargs, 

    logpath=None,
    nproc=1,
    jobname=None,

    intv_submit=workflow.DEFAULT_INTV_SUBMIT,
    intv_check=workflow.DEFAULT_INTV_CHECK,
    max_submit=workflow.DEFAULT_MAX_SUBMIT,

    raise_on_failure=False,
    remove_jobscript=True, 
    tmpfile_dir=os.getcwd(), 
):
    """Args:
        cmdargs: Must be a list. Each element must be a list or str representing
            command-line arguments.
    """
    # broadcast
    n_item = len(cmdargs)
    logpath = np.broadcast_to(logpath, n_item)  # even when n_item is 1, a 1d array results
    nproc = np.broadcast_to(nproc, n_item)
    jobname = np.broadcast_to(jobname, n_item)

    # sanity check
    if any(
        ((x is not None) and (not x.endswith('.log')))
        for x in logpath
    ):
        raise Exception(f'All "logpath" must end with ".log".')

    # modify cmdargs
    cmdargs = [
        workflow_utils.handle_cmdargs_arg(x, to_list=False)
        for x in cmdargs
    ]

    # make script files
    jobscript_dir = tempfile.mkdtemp(prefix='slurm_script_tmpdir_', dir=tmpfile_dir)
    script_path_list = list()
    jobinfo_iter = zip(cmdargs, logpath, nproc, jobname)
    for zidx, iter_item in tools.zenumerate(jobinfo_iter):
        one_cmdargs, one_logpath, one_nproc, one_jobname = iter_item
        script_path = os.path.join(jobscript_dir, f'{zidx}.sbatch')
        script_path_list.append(script_path)
        make_jobscript(
            jobscript_path=script_path, 
            lines=[one_cmdargs], 
            stdout=one_logpath,
            stderr=one_logpath,
            nproc=one_nproc,
            jobname=one_jobname,
        )

    # run slurm job
    joblist = JobList(
        script_path_list, 
        intv_submit=intv_submit,
        intv_check=intv_check,
        max_submit=max_submit,
    )
    joblist.run()
    total_success = joblist.success
    exitcode_list = joblist.get_exitcodes()
    each_success_list = [(x == 0) for x in exitcode_list]

    # rename log files
    for each_logpath, each_success in zip(logpath, each_success_list):
        if each_logpath is not None:
            workflow_utils.rename_logfile(each_logpath, each_success)

    # raise if failed
    if (not total_success) and raise_on_failure:
        raise Exception(f'One or more jobs were finished unsuccessfully.')

    # finish
    if remove_jobscript:
        shutil.rmtree(jobscript_dir)
    return each_success_list, exitcode_list, script_path_list



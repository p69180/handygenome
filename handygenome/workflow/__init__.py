import sys
import os
import re
import argparse
import subprocess
import time
import logging
import uuid
import textwrap
import itertools
import tempfile

import pyranges as pr

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


SLURMBIN = '/usr/local/slurm/bin'
SQUEUE = os.path.join(SLURMBIN, 'squeue')
SBATCH = os.path.join(SLURMBIN, 'sbatch')
SCONTROL = os.path.join(SLURMBIN, 'scontrol')
SCANCEL = os.path.join(SLURMBIN, 'scancel')
SINFO = os.path.join(SLURMBIN, 'sinfo')

PBSBIN = '/usr/local/pbs/bin/'
PBSNODES = os.path.join(PBSBIN, 'pbsnodes')
QSTAT = os.path.join(PBSBIN, 'qstat')

RE_PATS = dict()
RE_PATS['scontrol_split_job'] = re.compile('^([^=]+)(=(.*))?$')
RE_PATS['scontrol_split_nodes'] = re.compile(r'(?<=\S)\s+(?=\S+=)')

DEFAULT_INTV_CHECK = 60  # seconds
DEFAULT_INTV_SUBMIT = 1  # seconds
DEFAULT_MAX_SUBMIT = None  # no limits

HELP_WIDTH = 50
DESCRIPTION_WIDTH = 80


# logging #

DEFAULT_DATEFMT = '%Z %Y-%m-%d %H:%M:%S' # KST 2022-03-23 22:12:34
DEFAULT_LOG_FORMATTERS = {
    'without_name': logging.Formatter(
        fmt='[%(asctime)s %(levelname)s] %(message)s', 
        datefmt=DEFAULT_DATEFMT,
    ),
    'with_name': logging.Formatter(
        fmt='[%(asctime)s %(levelname)s] %(name)s: %(message)s', 
        datefmt=DEFAULT_DATEFMT,
    ),
}


def get_logger(name=None, formatter=None, level='info', stderr=True, 
               filenames=None, append=False):
    if name is None:
        #name = str(uuid.uuid4())
        name = __name__.split('.')[-1]

    if formatter is None:
        formatter = DEFAULT_LOG_FORMATTERS['with_name']

    loglevel = getattr(logging, level.upper())

    logger = logging.getLogger(name)
    logger.setLevel(loglevel)
    logger.propagate = False

    if stderr:
        sh = logging.StreamHandler()
        sh.setLevel(loglevel)
        sh.setFormatter(formatter)
        logger.addHandler(sh)

    if filenames is not None:
        assert isinstance(filenames, (tuple, list))
        for fname in filenames:
            fh = logging.FileHandler(fname, mode=('a' if append else 'w'))
            fh.setLevel(loglevel)
            fh.setFormatter(formatter)
            logger.addHandler(fh)

    return logger


def get_debugging_logger(title, verbose):
    return get_logger(
        name=str(uuid.uuid4()),
        level=('debug' if verbose else 'info'),
        formatter=logging.Formatter(
            fmt=f'[%(asctime)s.%(msecs)03d {title}] line %(lineno)d: %(message)s', 
            datefmt='%Z %Y-%m-%d %H:%M:%S'
        ),
    )


# filesystem-related functions

def get_split_filenames(n_file, outdir, prefix, suffix):
    """
    Returns:
        A list of file names each of which looks like 
            <outdir>/<prefix>000<suffix>.
        File creation is not done.
    """

    result = list()
    padded_indices = common.get_padded_indices(n_file)
    result = [os.path.join(outdir, prefix + idx_pad + suffix)
              for idx_pad in padded_indices]
    
    return result


def check_infile_validity(infile):
    """
    Raise:
        If infile does not exist or user does not have read permission
    """

    if os.path.exists(infile):
        if os.path.isfile(infile):
            if not os.access(infile, os.R_OK):
                raise Exception(f'You do not have read permission on '
                                f'"{infile}".')
        else:
            raise Exception(f"'{infile}' is not a regular file.")
    else:
        raise Exception(f"'{infile}' does not exist.")


def check_indir_validity(indir):
    if os.path.exists(indir):
        if os.path.isdir(indir):
            if not os.access(indir, os.R_OK | os.X_OK):
                raise Exception(f'You do not have read and execution '
                                f'permission on "{indir}".')
        else:
            raise Exception(f"'{indir}' is not a directory.")
    else:
        raise Exception(f"'{indir}' does not exist.")


def check_outdir_validity(outdir):
    """
    Return: 
        True if outdir already exists, otherwise False.

    Raise:
        If 1) outdir is an existing non-directory file, or 
            2) outdir is an existing non-empty directory
    """

    if not os.path.exists(outdir):
        if os.access(os.path.dirname(outdir), os.W_OK | os.X_OK):
            return False
        else:
            raise Exception(f'You do not have w/x permission on '
                            f'"dirname" of "{outdir}".')

    elif os.path.isfile(outdir):
        raise Exception(f"'{outdir}' is an existing regular file.")

    elif os.path.isdir(outdir):
        if os.access(outdir, os.W_OK | os.X_OK):
            if len(os.listdir(outdir)) != 0:
                raise Exception(
                    f'"{outdir}" is an existing directory and is not empty.')
            else:
                return True
        else:
            raise Exception(
                f'"{outdir}" is an existing directory which you do not '
                f'have w/x permission on.')

    else:
        raise Exception(f'"{outdir}" is an existing file which is neither '
                        f'a regular file nor a directory.')


def check_outfile_validity(outfile, must_not_exist=False):
    """
    Return: 
        True if outfile already exists, otherwise False.
    Raise:
        permission problem
    """

    if not os.path.exists(outfile):
        if os.access(os.path.dirname(outfile), os.W_OK | os.X_OK):
            return False
        else:
            raise Exception(f'You do not have w/x permission on dirname '
                            f'of "{outfile}".')

    elif os.path.isfile(outfile):
        if must_not_exist:
            raise Exception(f'Specified file "{outfile}" must not exist '
                            f'in advance.')
        else:
            return True

    elif os.path.isdir(outfile):
        raise Exception(f"'{outfile}' is an existing directory.")

    else:
        raise Exception(f'"{outfile}" is an existing file which is neither '
                        f'a regular file nor a directory.')


def get_tmpfile_path(prefix=None, suffix=None, dir=None, where=None, delete=False, is_dir=False):
    """Args:
        where: alias for dir
    """
    # alias handling
    if where is not None:
        dir = where
    # main
    if dir is None:
        dir = os.getcwd()

    if delete:
        fd, path = tempfile.mkstemp(prefix=prefix, suffix=suffix, dir=dir)
        os.close(fd)
        os.remove(path)
    else:
        if is_dir:
            path = tempfile.mkdtemp(prefix=prefix, suffix=suffix, dir=dir)
        else:
            fd, path = tempfile.mkstemp(prefix=prefix, suffix=suffix, dir=dir)
            os.close(fd)

    return path


def get_tmpdir_paths(subdirs, prefix=None, suffix=None, where=os.getcwd(), 
                     top_path=None, root_path=None):
    """Args:
        subdirs: An iterable which contains subdirectory basenames
    """

    root_name = 'root'
    assert root_name not in subdirs, (
        f'"subdirs" argument must not include the name "{root_name}".')

    # alias handling
    if (top_path is not None) and (root_path is None):
        root_path = top_path

    if root_path is None:
        root_path = os.path.abspath(
            get_tmpfile_path(prefix=prefix, suffix=suffix, where=where, 
                             delete=False, is_dir=True))
    else:
        exists = check_outdir_validity(root_path)
        if not exists:
            os.mkdir(root_path)

    tmpdir_paths = dict()
    tmpdir_paths[root_name] = root_path
    for basename in subdirs:
        tmpdir_paths[basename] = os.path.join(tmpdir_paths[root_name], 
                                              basename)
        os.mkdir(tmpdir_paths[basename])

    return tmpdir_paths


# ARGPARSE SETUP FUNCTIONS

class CustomFormatter(
        argparse.ArgumentDefaultsHelpFormatter,
        argparse.RawTextHelpFormatter):
    pass


def init_parser(description=None):
    if description is None:
        parser = argparse.ArgumentParser(
            description=None,
            formatter_class=CustomFormatter, add_help=False)
    else:
        parser = argparse.ArgumentParser(
            description=textwrap.fill(description, width=DESCRIPTION_WIDTH),
            formatter_class=CustomFormatter, add_help=False)

    required = parser.add_argument_group(
        title='REQUIRED', 
        description=textwrap.fill(
            'Required ones which accept 1 or more arguments.',
            width=DESCRIPTION_WIDTH))

    optional = parser.add_argument_group(
        title='OPTIONAL', 
        description=textwrap.fill(
            'Optional ones which accept 1 or more arguments.',
            width=DESCRIPTION_WIDTH))

    flag = parser.add_argument_group(
        title='FLAG', 
        description=textwrap.fill(
            'Optional ones which accept 0 argument.',
            width=DESCRIPTION_WIDTH))

    add_help_arg(flag)

    parser_dict = {'main': parser, 'required': required, 
                   'optional': optional, 'flag': flag}

    return parser_dict


def get_basic_parser():
    parser_dict = init_parser()

    add_infile_arg(parser_dict['required'], required=True)
    add_outfile_arg(parser_dict['required'], required=True, 
                    must_not_exist=True)
    add_fasta_arg(parser_dict['required'], required=True)

    add_refver_arg(parser_dict['optional'], required=False, choices='all')
    add_outfmt_arg(parser_dict['optional'], required=False)
    add_scheduler_args(parser_dict['optional'], default_parallel=1,
                       default_sched='slurm',
                       default_check_interval=DEFAULT_INTV_CHECK,
                       default_submit_interval=DEFAULT_INTV_SUBMIT)

    return parser_dict['main']


def add_infile_arg(parser, required=True, help=f'Input vcf file path.'):
    parser.add_argument(
        '-i', '--infile', dest='infile_path', required=required, 
        type=arghandler_infile, metavar='<input file path>', 
        help=textwrap.fill(help, width=HELP_WIDTH))


def add_infilelist_arg(
        parser, required=True, 
        help=f'One or more input file paths, separated with whitespaces.'):
    parser.add_argument(
        '--infilelist', dest='infile_path_list', required=required, 
        nargs='+', type=arghandler_infile, metavar='<input file path>', 
        help=textwrap.fill(help, width=HELP_WIDTH))


def add_indir_arg(parser, required=True, help=f'Input directory path.'):
    parser.add_argument(
        '--dir', dest='indir_path', required=required, type=arghandler_indir,
        metavar='<input directory path>',
        help=textwrap.fill(help, width=HELP_WIDTH))


def add_outfile_arg(
    parser, required=True, must_not_exist='ask',
    help='Output vcf file path. Must not exist in advance.',
):
    handler = get_arghandler_outfile(must_not_exist)
    parser.add_argument(
        '-o', '--outfile', dest='outfile_path', required=required, 
        type=handler, metavar='<output file path>',
        help=textwrap.fill(help, width=HELP_WIDTH),
    )


def add_outdir_arg(
        parser, required=True,
        help=('Output directory path. It will be created if it does not '
              'exist; otherwise it must be an empty directory.')):
    parser.add_argument(
        '--outdir', dest='outdir_path', required=True,
        type=arghandler_outdir, metavar='<output directory>',
        help=textwrap.fill(help, width=HELP_WIDTH))


def add_fasta_arg(parser, required=True, help=None):
    available_versions = ", ".join(common.DEFAULT_FASTA_PATHS.keys())

    if help is None:
        help=(f'A fasta file path or, alternatively, a reference genome '
              f'version (one of {available_versions}), '
              f'in which case a preset fasta file for the reference version '
              f'is used.')

    parser.add_argument(
        '-f', '--fasta', dest='fasta_path', required=required,
        type=arghandler_fasta, metavar='<fasta path>', 
        help=textwrap.fill(help, width=HELP_WIDTH))


def add_refver_arg(parser, required=True, choices='all', help=None):
    if choices == 'all':
        #allowed_vals = tuple(common.RefverDict.aliases.keys()) + tuple(common.RefverDict.converter.keys())
        allowed_vals = common.RefverDict.known_refvers
    elif choices == 'mouse':
        allowed_vals = ('mm9', 'mm10', 'mm39')
    elif choices == 'human':
        allowed_vals = ('hg18', 'hg19', 'hg38')
    else:
        allowed_vals = choices

    if help is None:
        help = f'Reference genome version. Must be one of {allowed_vals}.'

    parser.add_argument(
        '--refver', dest='refver', required=required, default=None, 
        choices=allowed_vals, metavar='<reference genome version>', 
        help=textwrap.fill(help, width=HELP_WIDTH))


def add_outfmt_arg(
        parser, required=False, default='z',
        help=(f'Output vcf file format (bcftools style). Must be one of: '
              f'v (uncompressed VCF), z (compressed VCF), '
              f'u (uncompressed BCF), b (compressed BCF).')):
    assert default in ('v', 'z', 'u', 'b')
    parser.add_argument(
        '-O', '--output-format', dest='mode_pysam', required=required, 
        default=default, choices=('v', 'z', 'u', 'b'), 
        type=(lambda x: common.PYSAM_FORMAT_DICT[x]),
        metavar='<output format>',
        help=textwrap.fill(help, width=HELP_WIDTH))


def add_parallel_arg(parser, required=False, default=1,
                     help=f'Number of parallelization.'):
    parser.add_argument(
        '-p', '--parallel', dest='parallel', required=required,
        default=default, type=int,
        metavar='<number of parallelization>', 
        help=textwrap.fill(help, width=HELP_WIDTH))


def add_sched_arg(
        parser, required=False, default='slurm',
        help=f'Parallelization method. Must be "slurm" or "local".'):
    assert default in ('slurm', 'local')
    parser.add_argument(
        '--sched', dest='sched', required=required, default=default, 
        choices=('slurm', 'local'), metavar='<"slurm" or "local">', 
        help=textwrap.fill(help, width=HELP_WIDTH))


def add_checkint_arg(
        parser, required=False, default=DEFAULT_INTV_CHECK,
        help=f'Slurm job status check interval in seconds.'):
    parser.add_argument(
        '--slurm-check-interval', dest='intv_check', required=required,
        default=default, type=float, metavar='<slurm check interval>', 
        help=textwrap.fill(help, width=HELP_WIDTH))


def add_submitint_arg(
        parser, required=False, default=DEFAULT_INTV_SUBMIT,
        help=f'Slurm job submission interval in seconds.'):
    parser.add_argument(
        '--slurm-submit-interval', dest='intv_submit', required=required,
        default=default, type=float, metavar='<slurm submit interval>', 
        help=textwrap.fill(help, width=HELP_WIDTH))


def add_maxsubmit_arg(
    parser, 
    required=False, 
    default=DEFAULT_MAX_SUBMIT,
    help=f'Allowed maximum number of pending and running jobs at a time.',
):
    parser.add_argument(
        '--slurm-max-submit', dest='max_submit', required=required,
        default=default, type=int, metavar='<slurm maximum submit>', 
        help=textwrap.fill(help, width=HELP_WIDTH)
    )


def add_index_arg(parser_dict,
                  help='If set, output vcf file is not indexed.'):
    parser_dict['flag'].add_argument(
        '--donot-index', dest='donot_index', action='store_true',
        help=textwrap.fill(help, width=HELP_WIDTH))

#######

def add_incl_region_arg(
        parser, required=False,
        help=(f'Bed file (gzipped or not) representing regions to be included.'
              f'If not set, whole genomic regions are used.')):
    parser.add_argument('--included-region-bed', dest='incl_bed_path', 
                        required=required, help=help)


def add_excl_region_arg(
        parser, required=False,
        help=(f'Bed file (gzipped or not) representing regions to be excluded.'
              f'If not set, no region is excluded.')):
    parser.add_argument('--excluded-region-bed', dest='excl_bed_path', 
                        required=required, help=help)


def get_incl_gr_from_args(args):
    if args.incl_region_path is None:
        if args.fasta_path is not None:
            chromdict = common.ChromDict(fasta_path=fasta_path)
            gr_incl = chromdict.to_interval_list().to_gr()
        elif args.refver is not None:
            chromdict = common.ChromDict(refver=refver)
            gr_incl = chromdict.to_interval_list().to_gr()
        else:
            raise Exception(
                f'"incl_region_path", "fasta_path", "refver" are all not set.')
    else:
        gr_incl = pr.read_bed(args.incl_region_path)

    return gr_incl


def get_merged_intervals_from_args(args):
    gr_incl = get_incl_gr_from_args(args)
    if args.excl_region_path is None:
        gr_merged = gr_incl
    else:
        gr_excl = pr.read_bed(args.excl_region_path)
        gr_merged = gr_incl.subtract(gr_excl)

    return common.IntervalList.from_gr(gr_merged)
            

########

def add_help_arg(parser, help='show this help message and exit'):
    parser.add_argument('-h', '--help', action='help',
                        help=textwrap.fill(help, width=HELP_WIDTH))


# functions which receive 'parser_dict' argument

def add_logging_args(parser_dict):
    parser_dict['optional'].add_argument(
        '--log', required=False, 
        help=textwrap.fill(
            f'If used, progress messages will be written to this file. '
            f'Existing file will be truncated.',
            width=HELP_WIDTH))

    parser_dict['flag'].add_argument(
        '--append-log', dest='append_log', action='store_true',
        help=textwrap.fill('If set, existing log file is not trucated.',
                           width=HELP_WIDTH))

    parser_dict['flag'].add_argument(
        '--silent', action='store_true',
        help=textwrap.fill(
            'If set, progress messages will not be printed to the terminal.',
            width=HELP_WIDTH))


def add_scheduler_args(
    parser_dict, 
    default_parallel=1, 
    default_sched='slurm',
    default_check_interval=DEFAULT_INTV_CHECK,
    default_submit_interval=DEFAULT_INTV_SUBMIT,
    default_max_submit=DEFAULT_MAX_SUBMIT,
):
    add_parallel_arg(parser_dict['optional'], required=False, default=default_parallel)
    add_sched_arg(parser_dict['optional'], required=False, default=default_sched)
    add_checkint_arg(parser_dict['optional'], required=False, default=default_check_interval)
    add_submitint_arg(parser_dict['optional'], required=False, default=default_submit_interval)
    add_maxsubmit_arg(parser_dict['optional'], required=False, default=default_max_submit)


def add_rmtmp_arg(parser_dict):
    parser_dict['flag'].add_argument(
        '--donot-remove-tmp', dest='donot_rm_tmp', action='store_true',
        help=textwrap.fill(
            'If set, temporary files are not removed (which are removed by '
            'default).',
            width=HELP_WIDTH))


# arghandlers

def arghandler_infile(arg):
    arg = os.path.abspath(arg)

    try:
        check_infile_validity(arg)
    except Exception as e:
        raise argparse.ArgumentTypeError(str(e))
    else:
        return arg


def arghandler_infilelist(arg):
    """Args:
        arg: A list of input file paths 
    """

    return [arghandler_infile(x) for x in arg]


def arghandler_indir(arg):
    arg = os.path.abspath(arg)

    try:
        check_indir_validity(arg)
    except Exception as e:
        raise argparse.ArgumentTypeError(str(e))
    else:
        return arg


def arghandler_outfile_ask(arg):
    arg = os.path.abspath(arg)

    try:
        exists = check_outfile_validity(arg, must_not_exist=False)
    except Exception as e:
        raise argparse.ArgumentTypeError(str(e))
    else:
        if exists:
            msg = (f'Specified output file "{arg}" already exists. Would '
                   f'you like to proceed with overwriting? (y/n) ')
            ans = input(msg)
            while True:
                if ans == 'y':
                    return arg
                elif ans == 'n':
                    sys.exit(1)
                else:
                    ans = input('Please enter "y" or "n". ')
                    continue
        else:
            return arg


def arghandler_outfile_mayexist(arg):
    arg = os.path.abspath(arg)

    try:
        check_outfile_validity(arg, must_not_exist=False)
    except Exception as e:
        raise argparse.ArgumentTypeError(str(e))
    else:
        return arg


def arghandler_outfile_mustbeabsent(arg):
    arg = os.path.abspath(arg)

    try:
        check_outfile_validity(arg, must_not_exist=True)
    except Exception as e:
        raise argparse.ArgumentTypeError(str(e))
    else:
        return arg


def get_arghandler_outfile(must_not_exist):
    assert must_not_exist in ('yes', 'no', 'ask'), (
        '"must_not_exist" argument must be on of "yes", "no", or "ask".'
    )
    if must_not_exist == 'yes':
        return arghandler_outfile_mustbeabsent
    elif must_not_exist == 'no':
        return arghandler_outfile_mayexist
    elif must_not_exist == 'ask':
        return arghandler_outfile_ask


def arghandler_outdir(arg):
    arg = os.path.abspath(arg)

    try:
        exists = check_outdir_validity(arg)
    except Exception as e:
        raise argparse.ArgumentTypeError(str(e))
    else:
        if not exists:
            os.mkdir(arg)
        return arg


def arghandler_fasta(arg):
    if arg in common.DEFAULT_FASTA_PATHS.get_valid_keys():
        return common.DEFAULT_FASTA_PATHS[arg]
    else:
        return arghandler_infile(arg)


# decorators

def get_deco_arghandler_infile(argname):
    def decorator(func):
        sig = inspect.signature(func)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.arguments[argname] = arghandler_infile(ba.arguments[argname])

            return func(*ba.args, **ba.kwargs)

        return wrapper

    return decorator


def get_deco_arghandler_infilelist(argname):
    def decorator(func):
        sig = inspect.signature(func)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.arguments[argname] = arghandler_infilelist(
                ba.arguments[argname])

            return func(*ba.args, **ba.kwargs)

        return wrapper

    return decorator


##########################################################


# SLURM JOB HANDLING

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
        logger
    """

    jobstates_pending = ('CONFIGURING', 'PENDING')
    jobstates_running = ('COMPLETING', 'RUNNING', 'RESV_DEL_HOLD', 
                         'REQUEUE_FED', 'REQUEUE_HOLD', 'RESIZING', 
                         'SIGNALING', 'SPECIAL_EXIT', 'STAGE_OUT', 'STOPPED', 
                         'SUSPENDED')
    jobstates_finished = ('CANCELLED', 'COMPLETED', 'BOOT_FAIL', 'DEADLINE', 
                          'FAILED', 'NODE_FAIL', 'OUT_OF_MEMORY', 'PREEMPTED', 
                          'TIMEOUT')
    jobstates_unknown = ('REVOKED',)

    @common.get_deco_num_set_differently(
        ('jobscript_path', 'jobscript_string'), 1)
    def __init__(self, jobscript_path=None, jobscript_string=None, 
                 verbose=True, logger=None):
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

        if logger is None:
            self.logger = get_logger(
                formatter=DEFAULT_LOG_FORMATTERS['without_name'], 
                stderr=verbose)
        else:
            self.logger = logger

    def __repr__(self):
        infostr = ', '.join(f'{key}: {getattr(self, key)}' 
                            for key in ('jobid', 'jobstate', 'exitcode'))
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
            self.logger.info(f'Cancelled a job: JobID - {self.jobid}')

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
                raise Exception(textwrap.dedent(f"""\
                    Unable to interpret JobState.
                    JobState: {self.jobstate}
                    scontrol result: {self.scontrol_result}"""))

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
                sbatch exit code: {p.returncode}""")
            raise Exception(self.sbatch_err_info)
        else:
            self.jobid = int(p.stdout.split()[-1])
            self.submitted = True
            self.logger.info(f'Submitted a job: JobID - {self.jobid}')

    def _submit_path(self):
        p = subprocess.run([SBATCH, self.jobscript_path],
                           capture_output=True, text=True)

        if p.returncode != 0:
            self.sbatch_err_info = textwrap.dedent(f"""\
                Job submission using sbatch failed.

                job script path: {self.jobscript_path}
                sbatch stdout: {p.stdout}
                sbatch stderr: {p.stderr}
                sbatch exit code: {p.returncode}""")
            raise Exception(self.sbatch_err_info)
        else:
            self.jobid = int(p.stdout.split()[-1])
            self.submitted = True
            self.logger.info(f'Submitted a job: JobID - {self.jobid}')


class JobList(list):
    """
    Attributes:
        logger
        success
        
    """

    def __init__(
        self, 
        jobscript_path_list, 
        intv_submit=DEFAULT_INTV_SUBMIT,
        intv_check=DEFAULT_INTV_CHECK,
        max_submit=DEFAULT_MAX_SUBMIT,
        verbose=True, 
        logpath=None, 
        logger=None,
    ):
        for jobscript_path in jobscript_path_list:
            self.append(
                Job(jobscript_path=jobscript_path, verbose=verbose, logger=logger)
            )

        self.intv_submit = intv_submit
        self.intv_check = intv_check
        self.max_submit = max_submit

        if logger is None:
            self.logger = get_logger(
                formatter=DEFAULT_LOG_FORMATTERS['without_name'], 
                stderr=verbose, filenames=[logpath],
            )
        else:
            self.logger = logger

        self.success = None
        self.sublists = dict()
        self.update()

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

    def submit_and_wait(self):
        def make_infostring(sublist_key):
            n_jobs = len(self.sublists[sublist_key])
            details = ', '.join(job.__repr__()
                                for job in self.sublists[sublist_key])
            return f'{n_jobs} ({details})'

        def log_status():
            n_notsubmit = len(self.sublists['notsubmit'])
            info_pending = make_infostring('pending')
            info_running = make_infostring('running')
            info_finished = make_infostring('finished')

            msg = textwrap.dedent(f"""\
                Current job status:
                  Not submitted yet: {n_notsubmit}

                  Pending: {info_pending}

                  Running: {info_running}

                  Finished: {info_finished}
                """)
            self.logger.info(msg)

        def log_epilogue():
            info_success = make_infostring('success')
            info_failure = make_infostring('failure')
            info_unknown = make_infostring('unknown')

            msg = textwrap.dedent(f"""\
                All finished.

                Successful jobs: {info_success}

                Failed jobs: {info_failure}

                Jobs with unknown exit statuses: {info_unknown}
                """)
            self.logger.info(msg)

        # main
        try:
            while True:
                self.update()
                self.submit()
                log_status()
                if all(job.status['finished'] for job in self):
                    break
                else:
                    time.sleep(self.intv_check)
                    continue
        except KeyboardInterrupt:
            self.logger.info(
                f'RECEIVED A KEYBOARD INTERRUPT; '
                f'will cancel all pending and running jobs with scancel,'
                f' then exit immediately.')
            for job in itertools.chain(self.sublists['pending'], 
                                       self.sublists['running']):
                job.cancel()
            raise SystemExit(1)
        else:
            self.set_sublists_exitcodes()
            log_epilogue()
            self.success = all(job.success for job in self)

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

    def set_sublists_statuses(self):
        self.sublists['notsubmit'] = [job for job in self 
                                      if not job.submitted]
        self.sublists['pending'] = [job for job in self 
                                    if job.status['pending']]
        self.sublists['running'] = [job for job in self 
                                    if job.status['running']]
        self.sublists['finished'] = [job for job in self 
                                     if job.status['finished']]

    def set_sublists_exitcodes(self):
        self.sublists['success'] = [job for job in self 
                                    if job.success is True]
        self.sublists['failure'] = [job for job in self 
                                    if job.success is False]
        self.sublists['unknown'] = [job for job in self 
                                    if job.success is None]


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
            e_msg = f'''\
scontrol show job finished with nonzero exit code.
commands: {cmdargs}
stdout: {p.stdout}
stderr: {p.stderr}
exit code: {p.returncode}'''
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


###############################


def run_jobs(
    jobscript_paths, 
    sched, 
    intv_check, 
    intv_submit, 
    max_submit, 
    logger, 
    log_dir, 
    raise_on_failure=True,
):
    assert sched in ('local', 'slurm')

    if sched == 'local':
        success, exitcode_list = run_jobs_local(jobscript_paths)
    elif sched == 'slurm':
        joblist = JobList(
            jobscript_paths, 
            intv_check=intv_check,
            intv_submit=intv_submit,
            max_submit=max_submit,
            logger=logger,
        )
        joblist.submit_and_wait()

        success = joblist.success
        exitcode_list = joblist.get_exitcodes()

    if raise_on_failure:
        if not success:
            raise SystemExit(
                (f'One or more jobs have finished unsuccessfully. Refer to '
                 f'log files in {log_dir}'))
        else:
            return None
    else:
        return success, exitcode_list


def run_jobs_local(jobscript_path_list):
    plist = list()
    for jobscript_path in jobscript_path_list:
        p = subprocess.Popen([jobscript_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        plist.append(p)

    for p in plist:
        p.wait()

    exitcode_list = [p.returncode for p in plist]
    success = all(x == 0 for x in exitcode_list)

    return success, exitcode_list


"""
JOB SCRIPT GENERATION
"""

def make_multiline_command(lines, leading_taps = 0):
    new_lines = list()
    new_lines.append( '\t'*leading_taps + lines[0] + ' \\' )
    new_lines.extend( [ '\t'*(leading_taps+1) + x + ' \\' for x in lines[1:-1] ] )
    new_lines.append( '\t'*(leading_taps+1) + lines[-1] )

    return '\n'.join(new_lines)


def make_jobscript_string(lines, shell = False, python = False, **kwargs):
    string_list = list()


    string_list.append('#!' + common.BASH)

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


def make_jobscript(jobscript_path, lines, **kwargs):
    check_outfile_validity(jobscript_path)
    jobscript_string = make_jobscript_string(lines, **kwargs)
    with open(jobscript_path, 'w') as outfile:
        outfile.write(jobscript_string)


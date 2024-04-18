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
import shlex

import pyranges as pr

import handygenome.logutils as logutils
import handygenome.interval as libinterval
import handygenome.tools as tools
import handygenome.deco as deco
import handygenome.refgenome.refgenome as refgenome
import handygenome.vcfeditor.misc as vcfmisc
from handygenome.utils.workflow_utils import MultiArgsList



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


def get_debugging_logger(title=None, verbose=True):
    if title is None:
        title = '%(module)s.%(funcName)s'

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
    return [
        os.path.join(outdir, prefix + idx_pad + suffix)
        for idx_pad in tools.get_padded_indices(n_file)
    ]


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


def get_tmpdir_paths(
    subdirs, 
    prefix=None, 
    suffix=None, 
    where=os.getcwd(), 
    top_path=None, 
    root_path=None,
):
    """Args:
        subdirs: An iterable which contains subdirectory basenames
    """
    root_name = 'root'
    assert root_name not in subdirs, (
        f'"subdirs" argument must not include the name "{root_name}".'
    )

    # alias handling
    if (top_path is not None) and (root_path is None):
        root_path = top_path

    if root_path is None:
        root_path = os.path.abspath(
            get_tmpfile_path(
                prefix=prefix, suffix=suffix, where=where, 
                delete=False, is_dir=True,
            )
        )
    else:
        exists = check_outdir_validity(root_path)
        if not exists:
            os.mkdir(root_path)

    tmpdir_paths = dict()
    tmpdir_paths[root_name] = root_path
    for basename in subdirs:
        tmpdir_paths[basename] = os.path.join(tmpdir_paths[root_name], basename)
        os.mkdir(tmpdir_paths[basename])

    return tmpdir_paths


def get_tmpdir_paths_vcfeditor(
    *,

    prefix=None, 
    suffix=None, 
    where=os.getcwd(), 
    top_path=None, 
    root_path=None,

    nproc=1,
    make_outfiles=True,
    make_logs=True,
    other_subdirs=[], 
):
    """Args:
        subdirs: An iterable which contains subdirectory basenames
    """
    root_name = 'root'
    assert root_name not in other_subdirs, (
        f'"other_subdirs" argument must not include the name "{root_name}".'
    )

    # alias handling
    if (top_path is not None) and (root_path is None):
        root_path = top_path

    # make root path
    if root_path is None:
        root_path = os.path.abspath(
            get_tmpfile_path(
                prefix=prefix, suffix=suffix, where=where, 
                delete=False, is_dir=True,
            )
        )
    else:
        exists = check_outdir_validity(root_path)
        if not exists:
            os.mkdir(root_path)

    # main
    tmpdir_paths = dict()
    tmpdir_paths[root_name] = root_path

    if make_outfiles:
        subdic = dict()
        subdic['top'] = os.path.join(root_path, 'split_outfiles')
        os.mkdir(subdic['top'])
        subdic['files'] = MultiArgsList(
            os.path.join(subdic['top'], f'{zidx}.vcf.gz')
            for zidx in tools.zrange(nproc)
        )
        tmpdir_paths['split_outfiles'] = subdic

    if make_logs:
        subdic = dict()
        subdic['top'] = os.path.join(root_path, 'logs')
        os.mkdir(subdic['top'])
        subdic['files'] = MultiArgsList(
            os.path.join(subdic['top'], f'{zidx}.log')
            for zidx in tools.zrange(nproc)
        )
        tmpdir_paths['logs'] = subdic

    if set(other_subdirs).intersection(tmpdir_paths.keys()):
        raise Exception(
            f'Names of custom subdirs overlap with preset '
            f'subdir names. custom: {other_subdirs}, preset: '
            f'{list(tmpdir_paths.keys())}'
        )
    for key in other_subdirs:
        top = os.path.join(root_path, key)
        os.mkdir(top)
        tmpdir_paths[key] = top

    return tmpdir_paths


############################
# ARGPARSE SETUP FUNCTIONS #
############################

class CustomFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    #argparse.RawDescriptionHelpFormatter, 
    argparse.RawTextHelpFormatter
):
    pass


def init_parser(description=None, formatter_class=None):
    if formatter_class is None:
        formatter_class = CustomFormatter

    if description is None:
        parser = argparse.ArgumentParser(
            description=None,
            formatter_class=formatter_class, add_help=False)
    else:
        parser = argparse.ArgumentParser(
            #description=textwrap.fill(description, width=DESCRIPTION_WIDTH),
            description=description,
            formatter_class=formatter_class, 
            add_help=False,
        )

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
    parser, 
    required=True, 
    help=f'One or more input file paths, separated with whitespaces.',
):
    parser.add_argument(
        '--infilelist', dest='infile_path_list', required=required, 
        nargs='+', type=arghandler_infile, metavar='<input file path>', 
        help=textwrap.fill(help, width=HELP_WIDTH),
    )


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
    available_versions = ", ".join(refgenome.REFVERINFO.list_known_refvers())

    if help is None:
        help=(f'A fasta file path or, alternatively, a reference genome '
              f'version (one of {available_versions}), '
              f'in which case a preset fasta file for the reference version '
              f'is used.')

    parser.add_argument(
        '-f', '--fasta', dest='fasta_path', required=required,
        type=arghandler_fasta, metavar='<fasta path>', 
        help=textwrap.fill(help, width=HELP_WIDTH))


def add_refver_arg(parser, default='hg19', required=False, choices='all', help=None, standardize=False):
    if choices == 'all':
        allowed_vals = refgenome.REFVERINFO.list_known_refvers()
    elif choices == 'mouse':
        allowed_vals = ('mm9', 'mm10', 'mm39')
    elif choices == 'human':
        allowed_vals = ('hg18', 'hg19', 'hg38')
    else:
        allowed_vals = choices

    if help is None:
        help = (
            f'Reference genome version. Must be one of: '
            f'{refgenome.REFVERINFO.list_known_refvers()}'
        )

    if standardize:
        type_arg = refgenome.standardize
    else:
        type_arg = str
    parser.add_argument(
        '--refver', 
        dest='refver', 
        required=required, 
        default=default, 
        #choices=allowed_vals, 
        metavar='<reference genome version>', 
        help=textwrap.fill(help, width=HELP_WIDTH),
        type=type_arg,
    )


def add_outfmt_arg(
        parser, required=False, default='z',
        help=(f'Output vcf file format (bcftools style). Must be one of: '
              f'v (uncompressed VCF), z (compressed VCF), '
              f'u (uncompressed BCF), b (compressed BCF).')):
    assert default in ('v', 'z', 'u', 'b')
    parser.add_argument(
        '-O', '--output-format', dest='mode_pysam', required=required, 
        default=default, choices=('v', 'z', 'u', 'b'), 
        type=(lambda x: vcfmisc.PYSAM_MODE_DICT[x]),
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


def add_jobname_arg(
    parser, 
    required=False, 
    default=None,
    help=(
        f'Slurm job name prefix. There are preset default values '
        f'for each application.'
    ),
):
    parser.add_argument(
        '--jobname', dest='jobname', required=required, default=default, 
        metavar='<Slurm job name prefix>', 
        help=textwrap.fill(help, width=HELP_WIDTH),
    )


def add_index_arg(
    parser_dict,
    help='If set, output vcf file is not indexed.',
):
    parser_dict['flag'].add_argument(
        '--donot-index', dest='donot_index', action='store_true',
        help=textwrap.fill(help, width=HELP_WIDTH),
    )


def add_tmpdirloc_arg(
    parser_dict,
    help='Location where temporary directory will be created.',
):
    parser_dict['optional'].add_argument(
        '--tmpdir-loc', 
        default=os.getcwd(),
        dest='tmpdirloc', 
        help=textwrap.fill(help, width=HELP_WIDTH),
    )


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
            chromdict = refgenome.ChromDict.from_fasta_path(args.fasta_path)
        elif args.refver is not None:
            chromdict = refgenome.ChromDict.from_refver(refver)
        else:
            raise Exception(
                f'"incl_region_path", "fasta_path", "refver" are all not set.'
            )

        #gr_incl = chromdict.to_interval_list().to_gr()
        gr_incl = chromdict.to_gr(assembled_only=False, as_gr=True)
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

    return libinterval.IntervalList.from_gr(gr_merged)
            

########

def add_help_arg(parser, help='show this help message and exit'):
    parser.add_argument('-h', '--help', action='help',
                        help=textwrap.fill(help, width=HELP_WIDTH))


# functions which receive 'parser_dict' argument

def add_logging_args(parser_dict):
    parser_dict['optional'].add_argument(
        '--log', 
        required=False, 
        help=textwrap.fill(
            f'If used, progress messages will be written to this file. '
            f'Existing file will be truncated.',
            width=HELP_WIDTH,
        )
    )

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
    add_jobname_arg(parser_dict['optional'], required=False)


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
    if arg in refgenome.REFVERINFO.list_known_refvers():
    #if arg in common.DEFAULT_FASTA_PATHS.get_valid_keys():
        return refgenome.get_fasta_path(arg)
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


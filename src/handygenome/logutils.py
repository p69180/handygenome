import os
import sys
import re
import logging
import datetime
import inspect
import itertools

import handygenome.tools as tools


LOCATION_LOGLEVELS = [logging.DEBUG]
LOCATION_EXCL_PAT = re.compile(
    '|'.join(
        (
            r'^<.*>$',
            r'^/tmp/ipykernel_',
            r'/IPython/core/interactiveshell\.py$',
            r'/IPython/core/async_helpers\.py$',
            r'/ipykernel/zmqshell\.py$',
            r'/ipykernel/ipkernel\.py$',
            r'/ipykernel/kernelbase\.py$',
            r'/asyncio/events\.py$',
            r'/asyncio/base_events\.py$',
            r'/tornado/platform/asyncio\.py$',
            r'/ipykernel/kernelapp\.py$',
            r'/traitlets/config/application\.py$',
            r'/ipykernel_launcher\.py$',
        )
    )
)


def _make_logger():
    LOGGER = logging.getLogger('handygenome_logger')
    LOGGER.propagate = False
    LOGGER.setLevel(logging.DEBUG)
    
    STREAMHANDLER = logging.StreamHandler()
    LOGGER.addHandler(STREAMHANDLER)
    
    STREAMHANDLER.setLevel(logging.DEBUG)

    return LOGGER, STREAMHANDLER

LOGGER, STREAMHANDLER = _make_logger()


# timestamp

def get_datestring():
    """Returns a string like '2021-12-06 11:49:55'"""
    return str(datetime.datetime.now()).split('.')[0]


def get_timestamp():
    """Returns a string like 'KST 2021-12-06 11:51:36'"""
    dt = datetime.datetime.now().astimezone()
    return f'{str(dt.tzinfo)} {str(dt).split(".")[0]}'


def print_err(*args, stderr=True, files=None, **kwargs):
    """Args:
        stderr: (Bool) Whether to write to stderr
        files: A list of file paths to which message is written
    """
    if stderr:
        print(*args, file=sys.stderr, flush=True, **kwargs)

    if files is not None:
        for fname in files:
            with open(fname, 'a') as f:
                print(*args, file=f, flush=True, **kwargs)


def print_timestamp(*args, **kwargs):
    print_err(f'[{get_timestamp()}]', *args, **kwargs)


printlog = print_timestamp


# main logger

#def log(msg, level='info', add_locstring=None, verbose_locstring=None):
def log(
    msg, 
    level='info', 
    add_locstring=None, locstring_mode=None,
    filename=None, append=False,
):
    # postprocess params
    level = getattr(logging, level.upper())
    if add_locstring is None:
        add_locstring = True
    if locstring_mode is None:
        if level == getattr(logging, 'DEBUG'):
            locstring_mode = 'verbose'
        else:
            locstring_mode = 'default'

    # main
    formatter = logging.Formatter(
        fmt=make_logformat(level, add_locstring=add_locstring),
        datefmt='%Z %Y-%m-%d %H:%M:%S', # KST 2022-03-23 22:12:34
    )
    # streamhandler
    STREAMHANDLER.setFormatter(formatter)
    # filehandler
    if filename is not None:
        fh = make_filehandler(
            filename, 
            level=level, 
            formatter=formatter, 
            append=append,
        )
        LOGGER.addHandler(fh)
    
    if add_locstring:
        locstring = make_locstring(mode=locstring_mode)
        LOGGER.log(level, msg, extra={'locstring': locstring})
    else:
        LOGGER.log(level, msg)

    if filename is not None:
        LOGGER.removeHandler(fh)


def log_old(msg, level='info', add_locstring=None, verbose_locstring=None):
    """Deprecated in order not to use the 'root' logger; 
    using 'root' logger makes rpy2 pacakge to emit log messages.
    """

    # postprocess params
    level = getattr(logging, level.upper())
    if add_locstring is None:
        add_locstring = True
    if verbose_locstring is None:
        verbose_locstring = (level == getattr(logging, 'DEBUG'))

    # main
    logging.basicConfig(
        format=make_logformat(level, add_locstring=add_locstring),
        datefmt='%Z %Y-%m-%d %H:%M:%S', # KST 2022-03-23 22:12:34
        level=level,
        force=True,
    )
    if add_locstring:
        locstring = make_locstring(verbose=verbose_locstring)
        logging.log(level, msg, extra={'locstring': locstring})
    else:
        logging.log(level, msg)


def make_filehandler(filename, level, formatter, append=False):
    fh = logging.FileHandler(filename, mode=('a' if append else 'w'))
    fh.setLevel(level)
    fh.setFormatter(formatter)
    return fh


def make_locstring(mode='default'):
    assert mode in ['default', 'verbose', 'last']

    if mode in ['verbose', 'last']:
        verbose_loc_list = [
            f'{os.path.basename(finfo.filename)}: {finfo.function} (lineno {finfo.lineno})'
            for finfo in inspect.stack()
            if LOCATION_EXCL_PAT.search(finfo.filename) is None
        ]

    if mode == 'verbose':
        return (
            '\n'
            + '\n--> '.join(verbose_loc_list)
            + '\n\n'
        )
    elif mode == 'last':
        return verbose_loc_list[-1] + ' |'
    elif mode == 'default':
        finfo = get_calling_frameinfo()
        return f'{os.path.basename(finfo.filename)}: {finfo.function} (lineno {finfo.lineno}) |'


def make_logformat(level, add_locstring=False):
    """Helper of 'log'"""
    if level == logging.DEBUG:
        levelcol = tools.COLORS['cyan']
    elif level == logging.INFO:
        levelcol = tools.COLORS['green']
    elif level == logging.WARNING:
        levelcol = tools.COLORS['yellow']
    elif level == logging.ERROR:
        levelcol = tools.COLORS['orange']
    elif level == logging.CRITICAL:
        levelcol = tools.COLORS['red']

    levelname = levelcol + '%(levelname)s' + tools.COLORS['end']

    if add_locstring:
        return f'[%(asctime)s.%(msecs)03d {levelname}] %(locstring)s %(message)s'
    else:
        return f'[%(asctime)s.%(msecs)03d {levelname}] %(message)s'


def get_calling_frameinfo():
    """Helper of 'log'"""
    frameinfo_groups = tuple(
        tuple(subiter) for key, subiter in
        itertools.groupby(inspect.stack(), key=(lambda x: x.filename))
    )
    assert len(frameinfo_groups) >= 2
    frameinfo = frameinfo_groups[1][0]
    return frameinfo


# line number logging

def iter_lineno_logging(line_iterator, logger, logging_lineno, msgfunc=None):
    if msgfunc is None:
        def msgfunc(NR):
            return f'Processing {NR:,}th line'

    if logging_lineno is None:
        for line in line_iterator:
            yield line
    else:
        NR = 0
        for line in line_iterator:
            NR += 1
            if NR % logging_lineno == 0:
                #logger.info(msgfunc(NR))
                log(msgfunc(NR), level='info')
            yield line


# making custom logger instance


#DEFAULT_DATEFMT = '%Z %Y-%m-%d %H:%M:%S' # KST 2022-03-23 22:12:34
#DEFAULT_LOG_FORMATTERS = {
#    'without_name': logging.Formatter(
#        fmt='[%(asctime)s %(levelname)s] %(message)s', 
#        datefmt=DEFAULT_DATEFMT,
#    ),
#    'with_name': logging.Formatter(
#        fmt='[%(asctime)s %(levelname)s] %(name)s: %(message)s', 
#        datefmt=DEFAULT_DATEFMT,
#    ),
#}
#
#def get_logger(
#    name=None, 
#    formatter=None, 
#    level='info', 
#    stderr=True, 
#    filenames=None, 
#    append=False,
#):
#    if name is None:
#        #name = str(uuid.uuid4())
#        name = __name__.split('.')[-1]
#
#    if formatter is None:
#        formatter = DEFAULT_LOG_FORMATTERS['with_name']
#
#    loglevel = getattr(logging, level.upper())
#
#    logger = logging.getLogger(name)
#    logger.setLevel(loglevel)
#    logger.propagate = False
#
#    if stderr:
#        sh = logging.StreamHandler()
#        sh.setLevel(loglevel)
#        sh.setFormatter(formatter)
#        logger.addHandler(sh)
#
#    if filenames is not None:
#        assert isinstance(filenames, (tuple, list))
#        for fname in filenames:
#            fh = logging.FileHandler(fname, mode=('a' if append else 'w'))
#            fh.setLevel(loglevel)
#            fh.setFormatter(formatter)
#            logger.addHandler(fh)
#
#    return logger
#
#
## funclogger
#
#def make_funclogger(level, name):
#    formatter = logging.Formatter(
#        #fmt=(
#        #    '[%(asctime)s.%(msecs)03d %(levelname)s] '
#        #    '%(module)s.%(funcName) line %(lineno)d: %(message)s'
#        #), 
#        fmt=(
#            '[%(asctime)s.%(msecs)s %(levelname)s] '
#            '%(module)s.%(funcName) line %(lineno)s: %(message)s'
#        ), 
#        datefmt='%Z %Y-%m-%d %H:%M:%S',
#    )
#    return get_logger(
#        level=level, 
#        name=name, 
#        formatter=formatter,
#        stderr=True, 
#        filenames=None, 
#        append=False,
#    )
#
#
#FUNCLOGGER_DEBUG = make_funclogger(level='debug', name='FUNCLOGGER_DEBUG')
#FUNCLOGGER_INFO = make_funclogger(level='info', name='FUNCLOGGER_INFO')
#
#
#def get_funclogger(verbose):
#    if verbose:
#        return FUNCLOGGER_DEBUG
#    else:
#        return FUNCLOGGER_INFO
#
#
#def printlog_funclogger(msg):
#    FUNCLOGGER_INFO.info(msg)


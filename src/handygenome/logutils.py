import sys
import logging
import datetime


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
                logger.info(msgfunc(NR))
            yield line


# making custom logger instance


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

def get_logger(
    name=None, 
    formatter=None, 
    level='info', 
    stderr=True, 
    filenames=None, 
    append=False,
):
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


# funclogger

def make_funclogger(level, name):
    formatter = logging.Formatter(
        #fmt=(
        #    '[%(asctime)s.%(msecs)03d %(levelname)s] '
        #    '%(module)s.%(funcName) line %(lineno)d: %(message)s'
        #), 
        fmt=(
            '[%(asctime)s.%(msecs)s %(levelname)s] '
            '%(module)s.%(funcName) line %(lineno)s: %(message)s'
        ), 
        datefmt='%Z %Y-%m-%d %H:%M:%S',
    )
    return get_logger(
        level=level, 
        name=name, 
        formatter=formatter,
        stderr=True, 
        filenames=None, 
        append=False,
    )


FUNCLOGGER_DEBUG = make_funclogger(level='debug', name='FUNCLOGGER_DEBUG')
FUNCLOGGER_INFO = make_funclogger(level='info', name='FUNCLOGGER_INFO')


def get_funclogger(verbose):
    if verbose:
        return FUNCLOGGER_DEBUG
    else:
        return FUNCLOGGER_INFO


def printlog_funclogger(msg):
    FUNCLOGGER_INFO.info(msg)


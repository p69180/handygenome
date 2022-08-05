import os
import inspect

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))


def setup_logger(args, logger_name=None, tmpdir_root=None, 
                 with_genlog=False):
    if logger_name is None:
        caller_frame = inspect.stack()[1].frame
        caller_module = inspect.getmodule(caller_frame)
        logger_name = caller_module.__name__.split('.')[-1]

    filenames = list()
    if args.log is not None:
        filenames.append(args.log)
    if with_genlog:
        genlog_path = os.path.join(tmpdir_root, 'genlog.txt')
        filenames.append(genlog_path)

    logger = workflow.get_logger(name=logger_name, 
                                 stderr=(not args.silent),
                                 filenames=filenames, 
                                 append=args.append_log)

    return logger

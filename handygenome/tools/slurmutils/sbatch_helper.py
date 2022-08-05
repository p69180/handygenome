import logging

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
toolsetup = importlib.import_module('.'.join([top_package_name, 'workflow', 'toolsetup']))


DEFAULT_NUM_SUBMITTED = 400


def argument_parser(cmdargs):
    def sanity_check(args):
        pass

    parser_dict = workflow.init_parser()
    workflow.add_indir_arg(parser_dict['required'], 
                           required=True,
                           help=f'Directory where slurm job scripts exist.')
    workflow.add_checkint_arg(parser_dict['optional'], required=False)
    workflow.add_submitint_arg(parser_dict['optional'], required=False)
    workflow.add_logging_args(parser_dict)

    parser_dict['optional'].add_argument(
        '-n', required=False, default=DEFAULT_NUM_SUBMITTED, type=int, 
        help=(f'The number of jobs to submitted simultaneously. '
              f'Default: {DEFAULT_NUM_SUBMITTED}'))
    parser_dict['flag'].add_argument(
        '--nowait', action='store_true',
        help=(f'If set, this program will not wait for the submitted '
              f'jobs to finish.'))

    args = parser_dict['main'].parse_args(cmdargs)
    sanity_check(args)

    return args


def main(cmdargs):
    args = argument_parser(cmdargs)
    logger = toolsetup.setup_logger(args)

    jobscript_path_list = common.listdir(args.indir_path)
    joblist = workflow.JobList(jobscript_path_list, logger=logger)
    joblist.submit(intv_submit=args.intv_submit)
    if not args.nowait:
        joblist.wait(intv_check=args.intv_check, edit_log_suffix=True)

    logger.info('ALL FINISHED SUCCESSFULLY')

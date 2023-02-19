import os
import re
import inspect
import textwrap
import stat

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


def handle_region_args(chromdict, incl_bed_path, excl_bed_path, num_split, regionfiles_dir):
    incl_intvlist = get_included_intvlist(chromdict, incl_bed_path, excl_bed_path)
    incl_intvlist_split = incl_intvlist.split(num=num_split)
    # write
    for zidx, intvlist in common.zenumerate(incl_intvlist_split):
        fname = os.path.join(regionfiles_dir, f"{zidx}.bed")
        intvlist.write_bed(fname)
    # make padded versions and write
    for intvlist in incl_intvlist_split:
        intvlist.slop(chromdict, b=5000)
    for zidx, intvlist in common.zenumerate(incl_intvlist_split):
        fname = os.path.join(regionfiles_dir, f"{zidx}.padded.bed")
        intvlist.write_bed(fname)


def get_included_intvlist(chromdict, incl_bed_path, excl_bed_path):
    if incl_bed_path is None:
        incl_intvlist = common.IntervalList.from_chromdict(chromdict)
    else:
        incl_intvlist = common.IntervalList.from_bed(incl_bed_path)

    if excl_bed_path is not None:
        excl_intvlist = common.IntervalList.from_bed(excl_bed_path)
        incl_intvlist = incl_intvlist.subtract(excl_intvlist)

    incl_intvlist.sort_intervals(chromdict)

    return incl_intvlist


def write_region_files(regionfiles_dir, incl_intvlist_split, incl_intvlist_split_padded):
    for zidx, intvlist in common.zenumerate(incl_intvlist_split):
        fname = os.path.join(regionfiles_dir, f"{zidx}.bed")
        intvlist.write_bed(fname)

    for zidx, intvlist in common.zenumerate(incl_intvlist_split_padded):
        fname = os.path.join(regionfiles_dir, f"{zidx}.padded.bed")
        intvlist.write_bed(fname)


def write_jobscripts(
    script_path_list,
    log_path_list,
    module_name,
    unit_job_func_name,
    kwargs_single,
    kwargs_multi,
    jobname_prefix,
    nproc,
):
    # sanity check
    arglens = list()
    arglens.append(len(script_path_list))
    arglens.append(len(log_path_list))
    arglens.extend(len(val) for val in kwargs_multi.values())
    if len(set(arglens)) != 1:
        raise Exception(
            f'The lengths of "script_path_list", "log_path_list", '
            f'and components of "kwargs_multi" are not all the same.'
        )

    # main
    tab = ' ' * 4
    for zidx, (script_path, log_path) in common.zenumerate(zip(script_path_list, log_path_list)):
        idx = common.zidx_to_idx(zidx)
        success_log_path = re.sub("\.log$", ".success", log_path)
        failure_log_path = re.sub("\.log$", ".failure", log_path)

        unit_job_func_args = list()
        for key, val in kwargs_single.items():
            unit_job_func_args.append(f"{key}={repr(val)},")
        for key, val in kwargs_multi.items():
            unit_job_func_args.append(f"{key}={repr(val[idx])},")
        unit_job_func_args = ('\n' + (7 * tab)).join(unit_job_func_args)

        script_contents = list()
        script_contents.append(
            textwrap.dedent(f"""\
                #!{common.PYTHON}

                #SBATCH -N 1
                #SBATCH -n 1

                #SBATCH -c {nproc}
                #SBATCH -o {os.devnull}
                #SBATCH -e {os.devnull}
                #SBATCH -J {jobname_prefix}_{zidx}

                import os
                import contextlib
                import traceback
                import sys
                sys.path.append({repr(common.PACKAGE_LOCATION)})
                from {module_name} import {unit_job_func_name}

                log = open({repr(log_path)}, 'w')
                with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                    try:
                        {unit_job_func_name}(
                            {unit_job_func_args}
                        )
                    except:
                        print(traceback.format_exc())
                        success = False
                    else:
                        success = True
                log.close()

                if success:
                    os.rename({repr(log_path)}, {repr(success_log_path)})
                else:
                    os.rename({repr(log_path)}, {repr(failure_log_path)})
                    raise SystemExit(1)"""
            )
        )

        script_contents = '\n'.join(script_contents)
        with open(script_path, "w") as outfile:
            outfile.write(script_contents)

        os.chmod(script_path, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)


def get_outfile_script_log_paths(split_infile_path_list, outfile_dir, script_dir, log_dir):
    split_outfile_path_list = list()
    script_path_list = list()
    log_path_list = list()

    for zidx, split_infile_path in common.zenumerate(split_infile_path_list):
        basename = os.path.basename(split_infile_path)
        split_outfile_path_list.append(
            os.path.join(outfile_dir, basename)
        )
        script_path_list.append(
            os.path.join(script_dir, f"{zidx}.sbatch")
        )
        log_path_list.append(
            os.path.join(log_dir, f"{zidx}.sbatch.log")
        )

    return split_outfile_path_list, script_path_list, log_path_list


def get_script_log_paths(script_dir, log_dir, num_split):
    script_path_list = list()
    log_path_list = list()

    for zidx in common.zrange(num_split):
        script_path_list.append(os.path.join(script_dir, f"{zidx}.sbatch"))
        log_path_list.append(os.path.join(log_dir, f"{zidx}.sbatch.log"))

    return script_path_list, log_path_list



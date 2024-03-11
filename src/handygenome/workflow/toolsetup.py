import os
import re
import inspect
import textwrap
import stat
import string

import pysam
import pyranges as pr

import handygenome
import handygenome.tools as tools
import handygenome.interval as libinterval
import handygenome.workflow as workflow
import handygenome.vcfeditor.misc as vcfmisc
import handygenome.ucscdata as ucscdata
from handygenome.genomedf.genomedf_base import GenomeDataFrameBase
#import handygenome.blacklist as blacklist


SPLIT_INFILE_PAT = re.compile(r'([0-9]+)([A-Za-z]*)(\.vcf\.gz)')
LETTERS = sorted(string.ascii_letters)


def make_next_split_vcf_path(fname):
    """Examples:
        012.vcf.gz => 012A.vcf.gz => 012B.vcf.gz => ...
        0114Z.vcf.gz => 0114a.vcf.gz => 0114b.vcf.gz => ...
        0234z.vcf.gz => 0234zA.vcf.gz => 0234zB.vcf.gz => ...
    """
    bname = os.path.basename(fname)
    dname = os.path.dirname(fname)
    mat = SPLIT_INFILE_PAT.fullmatch(bname)
    assert mat is not None

    chars = mat.group(2)
    if chars == '':
        next_chars = LETTERS[0]
    else:
        lastchar = chars[-1]
        if lastchar == LETTERS[-1]:
            next_chars = chars + LETTERS[0]
        else:
            next_chars = chars[:-1] + LETTERS[LETTERS.index(lastchar) + 1]

    next_bname = mat.group(1) + next_chars + mat.group(3)
    return os.path.join(dname, next_bname)

make_next_infile = make_next_split_vcf_path


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
    for zidx, intvlist in tools.zenumerate(incl_intvlist_split):
        fname = os.path.join(regionfiles_dir, f"{zidx}.bed")
        intvlist.write_bed(fname)
    # make padded versions and write
    for intvlist in incl_intvlist_split:
        intvlist.slop(chromdict, b=5000)
    for zidx, intvlist in tools.zenumerate(incl_intvlist_split):
        fname = os.path.join(regionfiles_dir, f"{zidx}.padded.bed")
        intvlist.write_bed(fname)


def get_included_intvlist(chromdict, incl_bed_path, excl_bed_path):
    if incl_bed_path is None:
        incl_intvlist = libinterval.IntervalList.from_chromdict(chromdict)
    else:
        incl_intvlist = libinterval.IntervalList.from_bed(incl_bed_path)

    if excl_bed_path is not None:
        excl_intvlist = libinterval.IntervalList.from_bed(excl_bed_path)
        incl_intvlist = incl_intvlist.subtract(excl_intvlist)

    incl_intvlist.sort_intervals(chromdict)

    return incl_intvlist


def write_region_files(regionfiles_dir, incl_intvlist_split, incl_intvlist_split_padded):
    for zidx, intvlist in tools.zenumerate(incl_intvlist_split):
        fname = os.path.join(regionfiles_dir, f"{zidx}.bed")
        intvlist.write_bed(fname)

    for zidx, intvlist in tools.zenumerate(incl_intvlist_split_padded):
        fname = os.path.join(regionfiles_dir, f"{zidx}.padded.bed")
        intvlist.write_bed(fname)


def write_jobscripts(
    script_path_list,
    log_path_list,
    slurm_log_pf_list,
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
    for key, val in kwargs_single.items():
        if val != eval(repr(val)):
            raise Exception(f'Keyword argument {key}: {val} cannot be represented as string')
    for key, vallist in kwargs_multi.items():
        for val in vallist:
            if val != eval(repr(val)):
                raise Exception(f'Keyword argument {key}: {val} cannot be represented as string')

    # main
    tab = ' ' * 4
    for zidx, (script_path, log_path, slurm_log_pf) in tools.zenumerate(
        zip(script_path_list, log_path_list, slurm_log_pf_list)
    ):
        idx = tools.zidx_to_idx(zidx)
        success_log_path = re.sub("\.log$", ".success", log_path)
        failure_log_path = re.sub("\.log$", ".failure", log_path)
        slurmlog_out_path = slurm_log_pf + '.out'
        slurmlog_err_path = slurm_log_pf + '.err'

        unit_job_func_args = list()
        for key, val in kwargs_single.items():
            unit_job_func_args.append(f"{key}={repr(val)},")
        for key, val in kwargs_multi.items():
            unit_job_func_args.append(f"{key}={repr(val[idx])},")
        unit_job_func_args = ('\n' + (7 * tab)).join(unit_job_func_args)

        script_contents = list()
        script_contents.append(
            textwrap.dedent(f"""\
                #!{handygenome.PARAMS["python"]}

                #SBATCH -N 1
                #SBATCH -n 1

                #SBATCH -c {nproc}
                #SBATCH -o {slurmlog_out_path}
                #SBATCH -e {slurmlog_err_path}
                #SBATCH -J {jobname_prefix}_{zidx}

                import os
                import contextlib
                import traceback
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

    for zidx, split_infile_path in tools.zenumerate(split_infile_path_list):
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

    for zidx in tools.zrange(num_split):
        script_path_list.append(os.path.join(script_dir, f"{zidx}.sbatch"))
        log_path_list.append(os.path.join(log_dir, f"{zidx}.sbatch.log"))

    return script_path_list, log_path_list


def make_infile_copy(infile_path, tmpdir_root, logger):
    infile_link_path = os.path.join(tmpdir_root, 'infile_link.vcf.gz')
    is_vcf, comp, is_bgzf = vcfmisc.get_vcf_format(infile_path)
    if is_bgzf:
        os.symlink(infile_path, infile_link_path)
        indexfile_path = vcfmisc.get_indexfile_path(infile_path)
        if indexfile_path is None:
            logger.info(
                f'Making index of input VCF (which will be saved in the temporary directory)'
            )
            infile_link_index_path = vcfmisc.make_index(infile_link_path)
        else:
            infile_link_index_path = f'{infile_link_path}.csi'
            os.symlink(indexfile_path, infile_link_index_path)
    else:
        logger.info(f'Making bgzipped copy of input VCF')
        in_vcf = pysam.VariantFile(infile_path)
        out_vcf = pysam.VariantFile(infile_link_path, mode='wz', header=in_vcf.header.copy())
        for vr in in_vcf.fetch():
            out_vcf.write(vr)

        out_vcf.close()
        in_vcf.close()

        logger.info(f'Making index of the copied VCF')
        infile_link_index_path = vcfmisc.make_index(infile_link_path)

    return infile_link_path, infile_link_index_path
            
        
def get_blacklist_gr(refver):
    cytoband_gdf = ucscdata.get_cytoband(refver=refver)
    centromere_gdf = cytoband_gdf.loc[cytoband_gdf['Stain'] == 'acen', :]
    blacklist_gdf = GenomeDataFrameBase.concat(
        [
            centromere_gdf, 
            GenomeDataFrameBase.from_frame(blacklist.HIGHDEPTH_BLACKLIST_GR, refver=refver),
        ]
    )
    blacklist_gdf.sort()
    blacklist_gr = blacklist_gdf.merge().gr

    #centromere_gr = cytoband_gr[(cytoband_gr.Stain == 'acen').to_list()]
    #blacklist_gr = pr.concat([centromere_gr, blacklist.HIGHDEPTH_BLACKLIST_GR])
    #blacklist_gr = blacklist_gr.sort().merge()

    return blacklist_gr


import re
import os
import argparse
import subprocess
import tempfile
import itertools
import shutil

import pysam

import handygenome.common as common
import handygenome.workflow as workflow
import handygenome.workflow.toolsetup as toolsetup
import handygenome.vcfeditor.varscan_editing as varscan_editing


PROGNAME = __name__.split(".")[-1]
JAVA_PATH = '/home/users/pjh/conda_bin/java'
VARSCAN2_JAR_PATH = '/home/users/tools/varscan2.4.2/VarScan.v2.4.2.jar'
SAMTOOLS_PATH = '/home/users/pjh/conda_bin/samtools'


def unit_job_1(fasta_path, regionfile_path, in_bam_path, out_pileup_path):
    common.print_timestamp('BEGINNING')

    common.print_timestamp('regionfile_path', regionfile_path)
    common.print_timestamp('in_bam_path', in_bam_path)
    common.print_timestamp('out_pileup_path', out_pileup_path)

    with open(regionfile_path) as infile:
        for line in infile:
            linesp = line.replace('\n', '').split('\t')
            chrom = linesp[0]
            start0 = int(linesp[1])
            end0 = int(linesp[2])
            start1 = start0 + 1
            end1 = end0

            args = [
                SAMTOOLS_PATH,
                'mpileup',
                '--no-BAQ',
                '-q', '0',
                '-Q', '0',
                '-f', f'{fasta_path}',
                '-r', f'{chrom}:{start1}-{end1}',
                '-o', f'{out_pileup_path}',
                f'{in_bam_path}',
            ]
            p = subprocess.run(args, capture_output=True, text=True)
            if p.returncode != 0:
                common.print_timestamp('Subprocess finished unsuccessfully. Stderr is:\n{p.stderr}')
                p.check_returncode()

    common.print_timestamp('FINISHED')


def unit_job_2(
    normal_pileup_path, 
    tumor_pileup_path, 
    vcf_prefix, 
    tumor_purity_opt,
    somatic_p_value_opt,
    min_coverage_normal_opt,
    min_coverage_tumor_opt,
):
    common.print_timestamp('BEGINNING')

    args = [
        JAVA_PATH, '-jar', VARSCAN2_JAR_PATH, 'somatic',
        normal_pileup_path,
        tumor_pileup_path,
        vcf_prefix,
        '--tumor-purity', str(tumor_purity_opt),
        '--somatic-p-value', str(somatic_p_value_opt),
        '--min-coverage-normal', str(min_coverage_normal_opt),
        '--min-coverage-tumor', str(min_coverage_tumor_opt),
        '--output-vcf', '1',
    ]
    p = subprocess.run(args, capture_output=True, text=True)
    if p.returncode != 0:
        common.print_timestamp('Subprocess finished unsuccessfully. Stderr is:\n{p.stderr}')
        p.check_returncode()

    common.print_timestamp('FINISHED')


def argument_parser(cmdargs):
    parser_dict = workflow.init_parser(
        description=(f"Runs gatk HaplotypeCaller with N parallel jobs.")
    )
    # required
    parser_dict['required'].add_argument(
        '--tbam',
        required=True,
        metavar='<tumor bam path>',
        dest='tbam_path',
        type=workflow.arghandler_infile,
        help=f'Tumor bam file path',
    )
    parser_dict['required'].add_argument(
        '--nbam',
        required=True,
        metavar='<normal bam path>',
        dest='nbam_path',
        type=workflow.arghandler_infile,
        help=f'Normal bam file path',
    )
    parser_dict['required'].add_argument(
        '--outfile-somatic',
        required=True,
        metavar='<somatic vcf path>',
        dest='somatic_outfile_path',
        type=workflow.get_arghandler_outfile(must_not_exist='ask'),
        help=f'Output VCF file path with "SOMATIC=TRUE" records.',
    )
    workflow.add_fasta_arg(parser_dict["required"], required=True)

    # optional
    parser_dict['optional'].add_argument(
        '--outfile-germline',
        required=False,
        metavar='<germline vcf path>',
        dest='germline_outfile_path',
        type=workflow.get_arghandler_outfile(must_not_exist='ask'),
        help=f'Output VCF file path with "SOMATIC=FALSE" records. If not set, germline VCF file is not written.',
    )

    workflow.add_outfmt_arg(parser_dict["optional"], required=False)
    workflow.add_incl_region_arg(parser_dict["optional"], required=False)
    workflow.add_excl_region_arg(parser_dict["optional"], required=False)
    workflow.add_logging_args(parser_dict)
    workflow.add_scheduler_args(parser_dict, default_parallel=1, default_sched="slurm")

    parser_dict["optional"].add_argument(
        '--java',
        required=False,
        dest='java_path',
        default=JAVA_PATH,
        help=f'Java executable path.',
    )
    parser_dict["optional"].add_argument(
        '--varscan-jar',
        required=False,
        dest='varscan_jar_path',
        default=VARSCAN2_JAR_PATH,
        help=f'VarScan2 jar file path.',
    )
    parser_dict["optional"].add_argument(
        '--tumor-purity',
        required=False,
        dest='tumor_purity_opt',
        default=1,
        type=float,
        help=f'"--tumor-purity" option given to "VarScan somatic" command. 1 is VarScan default.',
    )
    parser_dict["optional"].add_argument(
        '--somatic-p-value',
        required=False,
        dest='somatic_p_value_opt',
        default=0.05,
        type=float,
        help=f'"--somatic-p-value" option given to "VarScan somatic" command. 0.05 is VarScan default.',
    )
    parser_dict["optional"].add_argument(
        '--min-coverage-normal',
        required=False,
        dest='min_coverage_normal_opt',
        default=8,
        type=int,
        help=f'"--min-coverage-normal" option given to "VarScan somatic" command. 8 is VarScan default.',
    )
    parser_dict["optional"].add_argument(
        '--min-coverage-tumor',
        required=False,
        dest='min_coverage_tumor_opt',
        default=6,
        type=int,
        help=f'"--min-coverage-tumor" option given to "VarScan somatic" command. 8 is VarScan default.',
    )
    # flag
    workflow.add_rmtmp_arg(parser_dict)
    workflow.add_index_arg(parser_dict)
    # main
    args = parser_dict["main"].parse_args(cmdargs)
    return args


def make_tmpdir(infile_path):
    tmpdir_paths = workflow.get_tmpdir_paths(
        [
            "regions",
            "job1_scripts", "job1_logs", 
            "job2_scripts", "job2_logs", 
            "tumor_pileups", "normal_pileups",
            "split_results", 
        ],
        prefix=(os.path.basename(infile_path) + "_" + PROGNAME),
        where=os.path.dirname(infile_path),
    )

    return tmpdir_paths


def get_varscan2_version(java_path, varscan_jar_path):
    """Returns:
        VarScan2 version string, e.g. v2.4.2
    """
    p = subprocess.run(
        [java_path, '-jar', varscan_jar_path],
        capture_output=True,
        text=True,
    )
    return p.stderr.split('\n')[0].split()[1]


def write_job1_scripts(args, tmpdir_paths, num_split):
    script_path_list, log_path_list = toolsetup.get_script_log_paths(
        script_dir=tmpdir_paths['job1_scripts'], 
        log_dir=tmpdir_paths['job1_logs'], 
        num_split=(num_split * 2),
    )

    regionfile_path_list = list()
    for x in sorted(os.listdir(tmpdir_paths['regions'])):
        if 'padded' not in x:
            regionfile_path_list.append(os.path.join(tmpdir_paths['regions'], x))
    regionfile_path_list *= 2

    in_bam_path_list = (
        ([args.tbam_path] * num_split) 
        + ([args.nbam_path] * num_split)
    )

    out_pileup_path_list = list()
    for zidx in common.zrange(num_split):
        out_pileup_path_list.append(os.path.join(tmpdir_paths['tumor_pileups'], f'{zidx}.pileup'))
    for zidx in common.zrange(num_split):
        out_pileup_path_list.append(os.path.join(tmpdir_paths['normal_pileups'], f'{zidx}.pileup'))

    toolsetup.write_jobscripts(
        script_path_list=script_path_list,
        log_path_list=log_path_list,
        module_name=__name__,
        unit_job_func_name='unit_job_1',
        kwargs_single={
            'fasta_path': args.fasta_path,
        },
        kwargs_multi={
            'regionfile_path': regionfile_path_list,
            'in_bam_path': in_bam_path_list,
            'out_pileup_path': out_pileup_path_list,
        },
        jobname_prefix='VarScan2_wrapper_step1_makePileup',
        ncore_per_job=1,
    )

    return script_path_list


def write_job2_scripts(args, tmpdir_paths, num_split):
    script_path_list, log_path_list = toolsetup.get_script_log_paths(
        script_dir=tmpdir_paths['job2_scripts'], 
        log_dir=tmpdir_paths['job2_logs'], 
        num_split=num_split,
    )

    normal_pileup_path_list = list()
    tumor_pileup_path_list = list()
    vcf_prefix_list = list()
    for basename in sorted(os.listdir(tmpdir_paths['normal_pileups'])):
        normal_pileup_path_list.append(os.path.join(tmpdir_paths['normal_pileups'], basename))
    for basename in sorted(os.listdir(tmpdir_paths['tumor_pileups'])):
        tumor_pileup_path_list.append(os.path.join(tmpdir_paths['tumor_pileups'], basename))
    for npup_path, tpup_path in zip(normal_pileup_path_list, tumor_pileup_path_list):
        assert os.path.basename(npup_path) == os.path.basename(tpup_path)
        zidx = os.path.basename(npup_path).split('.')[0]
        vcf_prefix_list.append(os.path.join(tmpdir_paths['split_results'], zidx))

    toolsetup.write_jobscripts(
        script_path_list=script_path_list,
        log_path_list=log_path_list,
        module_name=__name__,
        unit_job_func_name='unit_job_2',
        kwargs_single={
            'tumor_purity_opt': args.tumor_purity_opt,
            'somatic_p_value_opt': args.somatic_p_value_opt,
            'min_coverage_normal_opt': args.min_coverage_normal_opt,
            'min_coverage_tumor_opt': args.min_coverage_tumor_opt,
        },
        kwargs_multi={
            'normal_pileup_path': normal_pileup_path_list,
            'tumor_pileup_path': tumor_pileup_path_list,
            'vcf_prefix': vcf_prefix_list,
        },
        jobname_prefix='VarScan2_wrapper_step2_runVarScan',
        ncore_per_job=1,
    )

    return script_path_list


def postprocess_and_merge(tmpdir_paths, args, fasta):
    infile_path_list = common.listdir(tmpdir_paths['split_results'])
    varscan_editing.modify_varscan_vcf(
        infile_path_list, 
        somatic_outfile_path=args.somatic_outfile_path, 
        germline_outfile_path=args.germline_outfile_path,
        fasta=fasta,
        mode_pysam=args.mode_pysam, 
        pass_only=True,
        index=(not args.donot_index),
    )


def main(cmdargs):
    args = argument_parser(cmdargs)

    # setup parameters
    fasta = pysam.FastaFile(args.fasta_path)
    chromdict = common.ChromDict(fasta=fasta)

    # make temporary directory tree
    tmpdir_paths = make_tmpdir(args.tbam_path)

    # setup logger
    logger = toolsetup.setup_logger(
        args=args, tmpdir_root=tmpdir_paths["root"], with_genlog=True,
    )
    logger.info("Beginning")

    # get VarScan version
    version = get_varscan2_version(args.java_path, args.varscan_jar_path)

    # make split intervals and write bed files
    toolsetup.handle_region_args(
        chromdict, 
        incl_bed_path=args.incl_bed_path,
        excl_bed_path=args.excl_bed_path,
        num_split=args.parallel,
        regionfiles_dir=tmpdir_paths["regions"],
    )

    # job1
    logger.info("Beginning job1 - making pileup for normal and tumor bam")
    script_path_list_job1 = write_job1_scripts(args=args, tmpdir_paths=tmpdir_paths, num_split=args.parallel)
    workflow.run_jobs(
        script_path_list_job1,
        sched=args.sched,
        intv_check=args.intv_check,
        intv_submit=args.intv_submit,
        logger=logger,
        log_dir=tmpdir_paths['job1_logs'],
        raise_on_failure=True,
    )
    logger.info("Finished job1")

    # job2
    logger.info("Beginning job2 - running VarScan")
    script_path_list_job2 = write_job2_scripts(args=args, tmpdir_paths=tmpdir_paths, num_split=args.parallel)
    workflow.run_jobs(
        script_path_list_job2,
        sched=args.sched,
        intv_check=args.intv_check,
        intv_submit=args.intv_submit,
        logger=logger,
        log_dir=tmpdir_paths['job2_logs'],
        raise_on_failure=True,
    )
    logger.info("Finished job2")
    
    # merge & clearance
    logger.info("Postprocess and merge")
    postprocess_and_merge(tmpdir_paths, args, fasta)

    if (not args.donot_rm_tmp):
        shutil.rmtree(tmpdir_paths['root'])

    logger.info("All successfully finished.")



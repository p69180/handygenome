'''
- gatk MarkDuplicatesSpark cannot take named pipe as infile
- gatk SortSam requires ~7G memory
'''

import re
import os
import stat
import argparse
import subprocess
import multiprocessing
import shutil
import textwrap
import pprint
import time
import itertools
import shlex

import handygenome
import handygenome.logutils as logutils
import handygenome.tools as tools
import handygenome.refgenome.refgenome as refgenome
import handygenome.workflow as workflow
import handygenome.workflow.parallel as libparallel
import handygenome.workflow.toolsetup as toolsetup
import handygenome.align.bwa as libbwa


BAI_SUFFIX = '.bai'
MKDUP_METRIC_SUFFIX = '.mkdup_metric'


def argument_parser():
    parser_dict = workflow.init_parser(
        description=textwrap.dedent(
            f'''\
            Runs bwa mem -> gatk MarkDuplicatesSpark -> gatk SetNmMdAndUqTags
            Created files:
                1) output bam: same as --outfile value
                2) bam index: <outfile>{BAI_SUFFIX}
                3) MarkDuplicatesSpark metric file: <outfile>{MKDUP_METRIC_SUFFIX} 
                    (created only if --make-mkdup-metric option is used)
            '''
        ),
    )

    # required
    parser_dict['required'].add_argument(
        '--fq1', 
        type=workflow.arghandler_infile,
        metavar=f'<read 1 fastq file>',
        help=f'Fastq file for read 1. May be gzipped.',
        required=True,
    )
    parser_dict['required'].add_argument(
        '--fq2', 
        type=workflow.arghandler_infile,
        metavar=f'<read 2 fastq file>',
        help=f'Fastq file for read 2. May be gzipped.',
        required=True,
    )
    parser_dict['required'].add_argument(
        '--outfile', 
        type=workflow.get_arghandler_outfile('yes'),
        metavar=f'<output bam path>',
        help=f'Path to the output bam. Must end with ".bam". Recommended name: <ID>.*.bam',
        required=True,
    )
    parser_dict['required'].add_argument(
        '--ID', 
        metavar=f'<ID tag value>',
        help=f'This will be used as ID tag for output bam file.',
        required=True,
    )

    # optional
    workflow.add_refver_arg(parser_dict['optional'], required=False)
    parser_dict['optional'].add_argument(
        '--fasta', 
        metavar=f'<fasta path>',
        help=f'If used, this fasta will be used. Otherwise, preset fasta file for "--refver" will be used.',
        required=False,
    )
    parser_dict['optional'].add_argument(
        '--SM', 
        metavar=f'<SM tag value>',
        help=f'This will be used as SM tag for output bam file. If not given, --ID value will be used.',
        required=False,
    )
    parser_dict['optional'].add_argument(
        '--job-prefix', 
        metavar=f'<prefix>',
        help=f'This will be used as prefix to slurm job name. If not given, --ID value plus "_" will be used.',
        required=False,
    )
    parser_dict['optional'].add_argument(
        '--odpd', 
        metavar=f'<optical density distance>',
        help=f'Argument for --optical-duplicate-pixel-distance option of MarkDuplicatesSpark',
        default=2500,
        type=int,
        required=False,
    )
    parser_dict['optional'].add_argument(
        '--partition-size', 
        metavar=f'<int>',
        help=f'Argument for --bam-partition-size option of MarkDuplicatesSpark. Seems to regulate memory usage.',
        default=0,
        type=int,
        required=False,
    )
    workflow.add_logging_args(parser_dict)
    workflow.add_scheduler_args(parser_dict, default_parallel=6, default_sched='slurm')

    # flag
    workflow.add_rmtmp_arg(parser_dict)
    parser_dict['flag'].add_argument(
        '-M', 
        action='store_true',
        help=f'If set, bwa mem will run with -M option.',
    )
    parser_dict['flag'].add_argument(
        '--make-mkdup-metric', 
        action='store_true',
        help=f'If set, metrics file from MarkDuplicatesSpark will be created. This is said to slow down process.',
    )

    # main
    args = parser_dict['main'].parse_args()
    return args


def sanitycheck_and_postprocess(args):
    # sanitycheck
    if not args.outfile.endswith('.bam'):
        raise Exception(f'--outfile value must end with ".bam".')
    if (args.fasta is None) and (args.refver is None):
        raise Exception(f'At least one of --fasta or --refver must be used.')

    # postprocess
    if args.fasta is None:
        args.fasta = refgenome.get_fasta_path(args.refver)
    if args.SM is None:
        args.SM = args.ID
    if args.job_prefix is None:
        args.job_prefix = args.ID + '_'

    # sanitycheck after processing
    if not libbwa.check_index(args.fasta):
        raise Exception(f'Input fasta is not indexed with bwa.')


#########################
# job script generation #
#########################

def make_markdup_args(
):
    """Markdup cannot take named pipe as infile"""
    markdup_args = [
        handygenome.PARAMS['gatk'], 
        'MarkDuplicatesSpark',
        '-I', pipe1,
        '-O', output_path,
        #'--MAX_RECORDS_IN_RAM', str(500000),
        '--bam-partition-size', str(partition_size),
        '--duplicate-scoring-strategy', 'SUM_OF_BASE_QUALITIES',
        '--optical-duplicate-pixel-distance', str(odpd),
        '--tmp-dir', mkdup_tmpdir_path,
        '--duplicate-tagging-policy', 'All',
        '--conf', f'spark.executor.cores={args.parallel}',
    ]
    if args.make_mkdup_metric:
        markdup_args.extend(
            ['-M', mkdup_metric_path]
        )


def make_job1_script(
    args,
    script_path, 

    bwa_log_path, 
    gatk_log_path,
    script_log_path, 

    output_path, 
    tmp_output_dir,
    sortsam_tmpdir_path,
):
    # make pipes
    pipe1 = os.path.join(tmp_output_dir, 'job1_pipe1')
    os.mkfifo(pipe1)

    # bwa
    bwa_args = [
        handygenome.PARAMS['bwa'], 'mem',
        '-Y', 
        '-t', str(args.parallel), 
        '-R', f'@RG\\tID:{args.ID}\\tSM:{args.SM}',
    ]
    if args.M:
        bwa_args.append('-M')
    bwa_args.extend(
        [
            args.fasta, 
            args.fq1, 
            args.fq2,
        ]
    )

    # gatk SortSam
    gatk_sortsam_args = [
        handygenome.PARAMS['gatk'],
        'SortSam',
        '-I', pipe1,
        '-O', output_path,
        '-SO', 'queryname',
        '--TMP_DIR', sortsam_tmpdir_path,
    ]
        
    with open(script_path, 'wt') as outfile:
        outfile.write(
            textwrap.dedent(
                f'''\
                #!{handygenome.PARAMS["bash"]}
                #SBATCH -J {args.job_prefix}ALIGNMENT_job1_bwa_namesort
                #SBATCH -N 1 -n 1 -c {args.parallel}
                #SBATCH -o {script_log_path}

                set -euo pipefail

                {shlex.join(bwa_args)} 1> {pipe1} 2> {bwa_log_path} &
                {shlex.join(gatk_sortsam_args)} &> {gatk_log_path} &
                wait
                '''
            )
        )

    os.system(f'chmod 700 {script_path}')


def make_job2_script(
    args,
    script_path, 

    job1_output_path,

    gatk_log_path,
    script_log_path, 

    mkdup_tmpdir_path,
    mkdup_metric_path,

    output_path, 

    odpd,  # optical duplicate distance
    partition_size,
):
    """Markdup cannot take named pipe as infile"""
    markdup_args = [
        handygenome.PARAMS['gatk'], 
        'MarkDuplicatesSpark',
        '-I', job1_output_path,
        '-O', output_path,
        '--bam-partition-size', str(partition_size),
        '--duplicate-scoring-strategy', 'SUM_OF_BASE_QUALITIES',
        '--optical-duplicate-pixel-distance', str(odpd),
        '--tmp-dir', mkdup_tmpdir_path,
        '--duplicate-tagging-policy', 'All',
        '--conf', f'spark.executor.cores={args.parallel}',
        '--create-output-bam-index',
    ]
    if args.make_mkdup_metric:
        markdup_args.extend(
            ['-M', mkdup_metric_path]
        )

    with open(script_path, 'wt') as outfile:
        outfile.write(
            textwrap.dedent(
                f'''\
                #!{handygenome.PARAMS["bash"]}
                #SBATCH -J {args.job_prefix}ALIGNMENT_job2_mkdup_sort
                #SBATCH -N 1 -n 1 -c {args.parallel}
                #SBATCH -o {script_log_path}

                set -euo pipefail

                {shlex.join(markdup_args)} &> {gatk_log_path}
                '''
            )
        )

    os.system(f'chmod 700 {script_path}')


def make_job3_script(
    args,
    script_path, 

    job2_output_path,

    gatk_log_path,
    script_log_path, 

    output_path, 
):
    gatk_args = [
        handygenome.PARAMS['gatk'], 
        'SetNmMdAndUqTags',
        '-I', job2_output_path,
        '-R', args.fasta,
        '-O', output_path,
        '--CREATE_INDEX', 'true',
    ]

    with open(script_path, 'wt') as outfile:
        outfile.write(
            textwrap.dedent(
                f'''\
                #!{handygenome.PARAMS["bash"]}
                #SBATCH -J {args.job_prefix}ALIGNMENT_job3_NmMdUq
                #SBATCH -N 1 -n 1 -c 1
                #SBATCH -o {script_log_path}

                set -euo pipefail

                {shlex.join(gatk_args)} &> {gatk_log_path}
                '''
            )
        )

    os.system(f'chmod 700 {script_path}')


##################
# tmp file paths #
##################

def make_log_paths(log_topdir):
    return {
        'job1_bwa': os.path.join(log_topdir, 'job1_bwa.log'),
        'job1_samtools': os.path.join(log_topdir, 'job1_samtools.log'),
        'job1_gatk': os.path.join(log_topdir, 'job1_gatk.log'),
        'job1_script': os.path.join(log_topdir, 'job1_script.log'),

        'job2_gatk': os.path.join(log_topdir, 'job2_gatk.log'),
        'job2_samtools': os.path.join(log_topdir, 'job2_samtools.log'),
        'job2_script': os.path.join(log_topdir, 'job2_script.log'),

        'job3_gatk': os.path.join(log_topdir, 'job3_gatk.log'),
        'job3_script': os.path.join(log_topdir, 'job3_script.log'),
    }


def make_script_paths(script_topdir):
    return {
        'job1': os.path.join(script_topdir, 'job1.sh'),
        'job2': os.path.join(script_topdir, 'job2.sh'),
        'job3': os.path.join(script_topdir, 'job3.sh'),
    }


def make_output_paths(topdir):
    return {
        'job1_tmp': os.path.join(topdir, 'job1_output.bam'),
        'job2_tmp': os.path.join(topdir, 'job2_output.bam'),
        'job2_mkdup_metric': os.path.join(topdir, 'job2_mkdup_metric.bam'),
        'job3_tmp': os.path.join(topdir, 'job3_output.bam'),
    }

########
# main #
########

def main():
    args = argument_parser()
    sanitycheck_and_postprocess(args)

    tmpdir_paths = workflow.get_tmpdir_paths(
        [
            'scripts', 
            'logs', 
            'tmp_outputs', 
            'markdup', 
            'sortsam',
        ],
        prefix = (
            os.path.basename(args.outfile)
            + f'_alignment_tmpdir_'
        ),
        where=os.path.dirname(args.outfile),
    )
    log_paths = make_log_paths(tmpdir_paths['logs'])
    script_paths = make_script_paths(tmpdir_paths['scripts'])
    output_paths = make_output_paths(tmpdir_paths['tmp_outputs'])

    logger = toolsetup.setup_logger(
        args=args,
        tmpdir_root=tmpdir_paths['root'],
        with_genlog=True,
    )
    logutils.log('Beginning')

    # make job scripts
    make_job1_script(
        args=args,
        script_path=script_paths['job1'], 

        bwa_log_path=log_paths['job1_bwa'], 
        gatk_log_path=log_paths['job1_gatk'],
        script_log_path=log_paths['job1_script'], 

        output_path=output_paths['job1_tmp'], 
        tmp_output_dir=tmpdir_paths['tmp_outputs'],
        sortsam_tmpdir_path=tmpdir_paths['sortsam'],
    )

    make_job2_script(
        args=args,
        script_path=script_paths['job2'], 

        job1_output_path=output_paths['job1_tmp'],

        gatk_log_path=log_paths['job2_gatk'],
        script_log_path=log_paths['job2_script'], 

        mkdup_tmpdir_path=tmpdir_paths['markdup'],
        mkdup_metric_path=output_paths['job2_mkdup_metric'],

        output_path=output_paths['job2_tmp'], 

        odpd=args.odpd,  # optical duplicate distance
        partition_size=args.partition_size,
    )

    make_job3_script(
        args=args,
        script_path=script_paths['job3'], 

        job2_output_path=output_paths['job2_tmp'],

        gatk_log_path=log_paths['job3_gatk'],
        script_log_path=log_paths['job3_script'], 

        output_path=output_paths['job3_tmp'], 
    )

    # run jobs
    def run_jobs(script_path, script_log_path, title):
        libparallel.run_jobs(
            jobscript_paths=[script_path], 
            sched=args.sched, 
            intv_check=args.intv_check, 
            intv_submit=args.intv_submit, 
            max_submit=args.max_submit, 
            logger=logger, 
            log_dir=None,
            job_status_logpath=os.path.join(tmpdir_paths['root'], 'job_status_log.gz'),
            raise_on_failure=True,
            log_paths=[script_log_path],
            title=title,
        )

    job1_title = f'Alignment step 1 - bwa -> gatk SortSam'
    logutils.log(f'Beginning {job1_title}')
    run_jobs(script_paths['job1'], log_paths['job1_script'], job1_title)

    job2_title = f'Alignment step 2 - gatk MarkDuplicatesSpark'
    logutils.log(f'Beginning {job2_title}')
    run_jobs(script_paths['job2'], log_paths['job2_script'], job2_title)

    job3_title = f'Alignment step 3 - gatk SetNmMdAndUqTags'
    logutils.log(f'Beginning {job3_title}')
    run_jobs(script_paths['job3'], log_paths['job3_script'], job3_title)

    # final
    os.rename(output_paths['job3_tmp'], args.outfile)
    os.rename(
        re.sub(r'\.bam$', '.bai', output_paths['job3_tmp']), 
        args.outfile + BAI_SUFFIX,
    )
    if os.path.exists(output_paths['job2_mkdup_metric']):
        os.rename(
            output_paths['job2_mkdup_metric'],
            args.outfile + MKDUP_METRIC_SUFFIX,
        )

    # rmtmp
    if not args.donot_rm_tmp:
        shutil.rmtree(tmpdir_paths['root'])
        
    logutils.log(f'All finished.')



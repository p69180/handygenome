import re
import os
import stat
import argparse
import subprocess
import shutil
import textwrap

import pysam
import numpy as np

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
toolsetup = importlib.import_module('.'.join([top_package_name, 'workflow', 'toolsetup']))
libreadstats = importlib.import_module('.'.join([top_package_name, 'annotation', 'readstats']))
varianthandler = importlib.import_module('.'.join([top_package_name, 'variantplus', 'varianthandler']))
variantplus = importlib.import_module('.'.join([top_package_name, 'variantplus', 'variantplus']))
vpfilter = importlib.import_module('.'.join([top_package_name, 'variantplus', 'vpfilter']))
infoformat = importlib.import_module('.'.join([top_package_name, 'variantplus', 'infoformat']))
indexing = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'indexing']))
split = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'split']))
concat = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'concat']))


def unit_job(split_infile_path, split_outfile_path, bam_path_list, id_list, 
             refver, no_matesearch):
    # basic setup
    bam_list = [pysam.AlignmentFile(x) for x in bam_path_list]
    fasta = pysam.FastaFile(common.DEFAULT_FASTA_PATHS[refver])
    chromdict = common.ChromDict(fasta=fasta)

    # edit vcf header
    in_vcf = pysam.VariantFile(split_infile_path, 'r')
    new_header = in_vcf.header.copy()
    update_header(new_header, id_list)
    out_vcf = pysam.VariantFile(split_outfile_path, 'wz', 
                                header=new_header)

    # loop over variant records
    for vr in in_vcf.fetch():
        new_vr = varianthandler.reheader(vr, new_header)
        update_new_vr(new_vr, refver, fasta, chromdict, bam_list, id_list,
                      no_matesearch)
        out_vcf.write(new_vr)

    out_vcf.close()
    in_vcf.close()


def mean(array, round_digits=3):
    if len(array) == 0:
        return None
    else:
        return round(np.mean(array), round_digits)


def update_header(vcfheader, id_list):
    libreadstats.ReadStats.add_meta(vcfheader)
    vpfilter.add_format_filter_meta(vcfheader)
    for sampleid in id_list:
        if sampleid not in vcfheader.samples:
            vcfheader.samples.add(sampleid)


def update_new_vr(new_vr, refver, fasta, chromdict, bam_list, id_list,
                  no_matesearch):
    vcfspec = varianthandler.get_vcfspec(new_vr)
    print(vcfspec, flush=True)  # for logging

    # get irrelevant sample ids: these need to be set '.' afterward
    FORMAT_key = libreadstats.ReadStats.meta['ID']
    irrelevant_samples = set(x for x in new_vr.samples.keys()
                             if infoformat.check_NA_format(new_vr, x, FORMAT_key))
    irrelevant_samples.difference_update(id_list)

    for bam, sampleid in zip(bam_list, id_list):
        readstats = libreadstats.get_readstats(
            vcfspec, bam, fasta, chromdict,
            rpplist_kwargs={'no_matesearch': no_matesearch},
        )
        readstats.write(new_vr, sampleid)

    for sampleid in irrelevant_samples:
        infoformat.set_NA_format(new_vr, sampleid, FORMAT_key)


#########################################################


def argument_parser(cmdargs):
    def sanity_check(args):
        # bam file arguments usage
        if not (((args.bamlist is not None) and 
                 (args.idlist is not None) and
                 (args.bamlist_file_path is None)) or
                ((args.bamlist is None) and 
                 (args.idlist is None) and
                 (args.bamlist_file_path is not None))):
            raise Exception(
                f'Allowed argument usages: '
                f'1) use --bamlist and --idlist, without --bamlist-file ; '
                f'2) use --bamlist-file, without --bamlist and --idlist')

        # --bamlist and --idlist length comparison
        if args.bamlist is not None and args.idlist is not None:
            if len(args.bamlist) != len(args.idlist):
                raise Exception(f'The number of --bamlist and --idlist '
                                f'arguments must be the same.')

    def postprocess(args):
        # bamlist file loading
        if args.bamlist_file_path is not None:
            bamlist = list()
            idlist = list()
            with open(args.bamlist_file_path, 'r') as infile:
                for line in infile:
                    linesp = line.replace('\n', '').split('\t')
                    if len(linesp) != 2:
                        raise Exception(
                            f'The number of tab-seperated fields is not 2 '
                            f'in the bamlist file.')

                    sampleid = linesp[0]
                    bam_path = linesp[1]
                    if not bam_path.endswith('.bam'):
                        raise Exception(
                            f'Bam file path does not end with ".bam" in the '
                            f'bamlist file. Please check if column 1 is '
                            f'sample ID and column 2 is bam file path.')

                    idlist.append(sampleid)
                    bamlist.append(bam_path)

            args.bamlist = bamlist
            args.idlist = idlist

    parser_dict = workflow.init_parser(
        description=(
            f'Calculates allele-supporting read count information for each '
            f'input bam, then writes to FORMAT. To specify input bam '
            f'files and sample IDs, use one of the following options: '
            f'1) use --bamlist and --idlist, without --bamlist-file ; '
            f'2) use --bamlist-file, without --bamlist and --idlist'))

    # required
    workflow.add_infile_arg(parser_dict['required'], required=True)
    workflow.add_outfile_arg(parser_dict['required'], required=True, 
                             must_not_exist='ask')
    workflow.add_refver_arg(parser_dict['required'], required=True, 
                            choices=('hg19', 'hg38'))

    # optional
    parser_dict['optional'].add_argument(
        '--bamlist', dest='bamlist', required=False,
        nargs='+', type=workflow.arghandler_infile, 
        metavar='<Input bam file path>',
        help=f'One or more input bam file paths separated by whitespaces.')
    parser_dict['optional'].add_argument(
        '--idlist', dest='idlist', required=False,
        nargs='+', 
        metavar='<Sample ID>',
        help=(f'Sample IDs of the input bam files, in the same order, '
              f'separated by whitespaces.'))
    parser_dict['optional'].add_argument(
        '--bamlist-file', dest='bamlist_file_path', required=False,
        type=workflow.arghandler_infile,
        metavar='<bam list file path>',
        help=(f'A 2-column tab-separated file which contains sample IDs on '
              f'the first column and bam file paths on the second column.'))

    workflow.add_outfmt_arg(parser_dict['optional'], required=False)
    workflow.add_logging_args(parser_dict)
    workflow.add_scheduler_args(parser_dict, default_parallel=1, 
                                default_sched='slurm')

    # flag
    workflow.add_rmtmp_arg(parser_dict)
    workflow.add_index_arg(parser_dict)
    parser_dict['flag'].add_argument(
        '--no-matesearch', dest='no_matesearch', action='store_true',
        help=(f'If set, "no_matesearch" argument is set to True when '
              f'creating ReadPlusPairList objects. It can be much faster, '
              f'but mate reads far away from the variant site are missed.'))

    args = parser_dict['main'].parse_args(cmdargs)
    sanity_check(args)
    postprocess(args)

    return args


def write_jobscripts(tmpdir_paths, split_infile_path_list, 
                     bam_path_list, id_list, refver, no_matesearch,
                     jobname_prefix='readstatsannot', ncore_perjob=1):
    jobscript_path_list = list()
    split_outfile_path_list = list()

    for zidx, split_infile_path in common.zenumerate(split_infile_path_list):
        basename = os.path.basename(split_infile_path)

        split_outfile_path = os.path.join(tmpdir_paths['split_outfiles'], 
                                          basename)
        split_outfile_path_list.append(split_outfile_path)

        jobscript_path = os.path.join(tmpdir_paths['scripts'], 
                                      f'{zidx}.sbatch')
        jobscript_path_list.append(jobscript_path)
        logpath = os.path.join(tmpdir_paths['logs'], 
                               os.path.basename(jobscript_path) + '.log')
        success_logpath = re.sub('\.log$', '.success', logpath)
        failure_logpath = re.sub('\.log$', '.failure', logpath)

        script_contents = textwrap.dedent(f"""\
            #!{common.PYTHON}

            #SBATCH -N 1
            #SBATCH -n 1

            #SBATCH -c {ncore_perjob}
            #SBATCH -o {os.devnull}
            #SBATCH -e {os.devnull}
            #SBATCH -J {jobname_prefix}_{zidx}

            import os
            import contextlib
            import traceback
            import sys
            sys.path.append('{common.PACKAGE_LOCATION}')
            from {__name__} import unit_job

            log = open('{logpath}', 'w')
            with contextlib.redirect_stdout(log), \\
                    contextlib.redirect_stderr(log):
                try:
                    unit_job(
                        split_infile_path='{split_infile_path}',
                        split_outfile_path='{split_outfile_path}',
                        bam_path_list={str(bam_path_list)},
                        id_list={str(id_list)},
                        refver='{refver}',
                        no_matesearch={no_matesearch},
                        )
                except:
                    print(traceback.format_exc())
                    success = False
                else:
                    success = True
            log.close()

            if success:
                os.rename('{logpath}', '{success_logpath}')
            else:
                os.rename('{logpath}', '{failure_logpath}')
                raise SystemExit(1)
            """)

        with open(jobscript_path, 'w') as outfile:
            outfile.write(script_contents)

        os.chmod(jobscript_path, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)

    return jobscript_path_list, split_outfile_path_list


def main(cmdargs):
    args = argument_parser(cmdargs)

    # make tmpdir tree
    tmpdir_paths = workflow.get_tmpdir_paths(
        ['scripts', 'logs', 'split_infiles', 'split_outfiles'],
        prefix='_'.join([os.path.basename(args.infile_path),
                         __name__.split('.')[-1]]),
        where=os.path.dirname(args.infile_path))

    # setup logger
    logger = toolsetup.setup_logger(args=args,
                                    tmpdir_root=tmpdir_paths['root'],
                                    with_genlog=True)
    logger.info('Beginning')

    # split the input file
    logger.info('Splitting the input file')
    split_infile_path_list = split.main(vcf_path=args.infile_path, 
                                        outdir=tmpdir_paths['split_infiles'], 
                                        n_file=args.parallel, 
                                        mode_bcftools='z', prefix='', 
                                        suffix='.vcf.gz')

    # make job scripts and run
    jobscript_path_list, split_outfile_path_list = write_jobscripts(
        tmpdir_paths, split_infile_path_list, args.bamlist, args.idlist, 
        args.refver, args.no_matesearch)
    logger.info('Running annotation jobs for each split file')
    workflow.run_jobs(jobscript_path_list, sched=args.sched, 
                      intv_check=args.intv_check, 
                      intv_submit=args.intv_submit, 
                      logger=logger, log_dir=tmpdir_paths['logs'],
                      raise_on_failure=True)

    # concatenates split files
    logger.info('Merging split files')
    concat.main(infile_path_list=split_outfile_path_list,
                outfile_path=args.outfile_path, mode_pysam=args.mode_pysam,
                outfile_must_not_exist='no')

    # rmtmp, index
    if not args.donot_rm_tmp:
        shutil.rmtree(tmpdir_paths['root'])
    if not args.donot_index:
        indexing.index_vcf(args.outfile_path)

    logger.info('All successfully finished')

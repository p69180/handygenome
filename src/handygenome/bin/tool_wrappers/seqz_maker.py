import os
import re
import argparse
import gzip
import textwrap

import pysam
import pyranges as pr

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
sequenza = importlib.import_module('.'.join([top_package_name, 'cnv', 'sequenza']))


PROGNAME = 'seqz_maker'


def unit_job(tbam_path, nbam_path, regionfile_path, split_outfile_subdir,
             refver, fasta_path, gcwiggle_path, hom, het):
    intvlist = common.IntervalList.from_gr(pr.read_bed(regionfile_path))
    for zidx, intv in common.zenumerate(intvlist):
        split_outfile_path = os.path.join(split_outfile_subdir, 
                                          f'{zidx}.seqz.gz')
        sequenza.seqzutils_bam2seqz(
            tbam_path=tbam_path, nbam_path=nbam_path, 
            outfile_path=split_outfile_path, 
            refver=refver, fasta_path=fasta_path, 
            gcwiggle_path=gcwiggle_path, 
            force_outfile_name=False,
            chrom=intv.chrom, start0=intv.start0, end0=intv.end0,
            hom=hom, het=het)


def argument_parser(cmdargs):
    def sanity_check(args):
        sequenza.bam2seqz_sanity_check(
            args.force_outfile_name, args.outfile_path, args.refver, 
            args.fasta_path, args.gcwiggle_path)

    parser_dict = workflow.init_parser(
        description=(
            f'1) Whole genome is split into N intervals; '
            f'2) sequenza-utils bam2seqz is run for each interval; '
            f'3) Split results are concatenated; '
            f'4) sequenza-utils seqz_binning is done.\n\n'
            f'Preset fasta and gcwiggle files are used when --refver option '
            f'is used.\n'
            f'When --refver is set, --fasta and --gcwiggle must not be set.\n'
            f'When --refver is not set, --fasta and --gcwiggle must be set.\n'
            ))
    # required
    workflow.add_outfile_arg(parser_dict['required'], required=True, 
                             must_not_exist='ask')
    parser_dict['required'].add_argument(
        '-t', '--tbam', dest='tbam_path', required=True,
        help=f'Tumor bam file path.')
    parser_dict['required'].add_argument(
        '-n', '--nbam', dest='nbam_path', required=True,
        help=f'Normal bam file path.')

    # optional
    workflow.add_refver_arg(parser_dict['optional'], required=False, 
                            choices=('hg19', 'hg38'))
    workflow.add_fasta_arg(parser_dict['optional'], required=False)

    workflow.add_incl_region_arg(parser_dict['optional'], required=False)
    workflow.add_excl_region_arg(parser_dict['optional'], required=False)

    workflow.add_logging_args(parser_dict)
    workflow.add_scheduler_args(parser_dict, default_parallel=1, 
                                default_sched='slurm')

    parser_dict['optional'].add_argument(
        '--gcwiggle', dest='gcwiggle', required=False,
        help=f'GC wiggle file path.')
    parser_dict['optional'].add_argument(
        '--hom-threshold', dest='hom', required=False,
        type=float, default=sequenza.DEFAULT_HOM,
        help=(f'--hom option value to be used in sequenza-utils bam2seqz. '
              f'Represents VAF threshold to treat as homozygous SNP.'))
    parser_dict['optional'].add_argument(
        '--het-threshold', dest='het', required=False,
        type=float, default=sequenza.DEFAULT_HET,
        help=(f'--het option value to be used in sequenza-utils bam2seqz. '
              f'Represents VAF threshold to treat as heterozygous SNP.'))

    # flag
    workflow.add_rmtmp_arg(parser_dict)
    parser_dict['flag'].add_argument(
        '--force-outfile-name', dest='force_outfile_name', 
        action='store_true',
        help=f'When not set, output file name must end with ".seqz.gz".')

    args = parser_dict['main'].parse_args(cmdargs)
    return args


def make_tmpdir(infile_path):
    tmpdir_paths = workflow.get_tmpdir_paths(
        ['scripts', 'logs', 'regions', 'split_outfiles'],
        prefix=(os.path.basename(infile_path) + '_' + PROGNAME),
        where=os.path.dirname(infile_path))

    return tmpdir_paths


def make_regionfiles(args, regions_dir):
    regionfile_path_list = list()
    if args.refver is not None:
        chromdict = common.ChromDict(refver=refver)
    #elif args.fasta_path is not None:
    #    chromdict = common.ChromDict(fasta_path=fasta_path)

    merged_intvlist = workflow.get_merged_intervals_from_args(args)
    split_intvlists = merged_intvlist.get_split_interval_lists(
        chromdict, num=args.parallel)

    for zidx, intvlist in common.zenumerate(split_intvlists):
        outfile_path = os.path.join(regions_dir, f'{zidx}.bed')
        regionfile_path_list.append(outfile_path)
        intvlist.write_bed(outfile_path)

    return regionfile_path_list


def make_split_outfiles_subdirs(split_outfiles_dir, parallel):
    split_outfiles_subdirs = list()
    for zidx in common.zrange(parallel):
        subdir = os.path.join(split_outfiles_dir, zidx)
        os.mkdir(subdir)
        split_outfiles_subdirs.append(subdir)

    return split_outfiles_subdirs


def write_jobscripts(tmpdir_paths, split_outfiles_subdirs,

                        tbam_path
                        nbam_path
                        regionfile_path
                        split_outfile_subdir={repr(split_outfile_subdir)},
                        refver={repr(refver)},
                        fasta_path={repr(fasta_path)}, 
                        gcwiggle_path={repr(gcwiggle_path)}, 
                        hom={repr(hom)}, 
                        het={repr(het)},

                     jobname_prefix=PROGNAME, ncore_perjob=1):
    jobscript_path_list = list()

    for zidx, outfile_subdir in common.zenumerate(split_outfiles_subdirs):
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
            with \\
                    contextlib.redirect_stdout(log), \\
                    contextlib.redirect_stderr(log):
                try:
                    unit_job(
                        tbam_path={repr(tbam_path)}, 
                        nbam_path={repr(nbam_path)}, 
                        regionfile_path={repr(regionfile_path)},
                        split_outfile_subdir={repr(split_outfile_subdir)},
                        refver={repr(refver)},
                        fasta_path={repr(fasta_path)}, 
                        gcwiggle_path={repr(gcwiggle_path)}, 
                        hom={repr(hom)}, 
                        het={repr(het)},
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

    return jobscript_path_list


def main(cmdargs):
    args = argument_parser(cmdargs)

    # setups logger
    logger = workflow.get_logger(name=PROGNAME, 
                                 stderr=(not args.silent),
                                 filenames=[args.log], append=False)
    logger.info('Beginning')

    # make 
    tmpdir_paths = make_tmpdir(args.tbam_path)
    regionfile_path_list = make_regionfiles(args, tmpdir_paths['regions_dir'])
    split_outfiles_subdirs = make_split_outfiles_subdirs(split_outfiles_dir, 
                                                         parallel)

    jobscript_path_list = write_jobscripts(
        tmpdir_paths, split_infile_path_list, fasta_path, 
        species, assembly, distance, refver,
        do_features, do_cosmic, do_popfreq,
        j)


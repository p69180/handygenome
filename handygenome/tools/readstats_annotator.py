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
import gc
import itertools

import pysam
import numpy as np
import pyranges as pr
import pandas as pd

import handygenome.common as common
import handygenome.workflow as workflow
import handygenome.workflow.toolsetup as toolsetup
import handygenome.annotation.readstats as libreadstats
import handygenome.variant.varianthandler as varianthandler
import handygenome.variant.variantplus as variantplus
import handygenome.variant.filter as libfilter
import handygenome.variant.infoformat as infoformat
import handygenome.vcfeditor.indexing as indexing
import handygenome.vcfeditor.split as libsplit
import handygenome.vcfeditor.concat as libconcat
import handygenome.variant.vcfspec as libvcfspec
import handygenome.variant.ponbams as libponbams
import handygenome.bameditor as bameditor
import handygenome.ucscdata as ucscdata
import handygenome.blacklist as blacklist
import handygenome.vcfeditor.misc as vcfmisc


LOGGER_NAME = __name__.split('.')[-1]


def unit_job(
    split_outfile_path,
    fetchregion,
    infile_copy_path, 
    bam_path_list, 
    id_list, 
    pon_bam_path_list, 
    pon_id_list, 
    refver, 
    no_matesearch, 
    countonly, 
    memuse_limit_gb,
    depth_limits,
    mq_limits,
    include_blacklist,
    monitor_interval=10,
):
    # setup manager
    manager = multiprocessing.Manager()
    read_vrspecs = manager.list()
    shareddict = manager.dict()
    shareddict['finished'] = False
    shareddict['outfile_path'] = split_outfile_path
    shareddict['parent_memuse_gb'] = common.get_rss(mode='total', unit='g')
    common.print_timestamp(f"Memory usage: {shareddict['parent_memuse_gb']} GB")

    # setup paramaters
    blacklist_gr = toolsetup.get_blacklist_gr(refver)

    # make subprocess arguments
    subproc_kwargs = {
        'shareddict': shareddict,
        'read_vrspecs': read_vrspecs,

        'infile_copy_path': infile_copy_path, 
        'fetchregion': fetchregion,
        'bam_path_list': bam_path_list, 
        'id_list': id_list, 
        'pon_bam_path_list': pon_bam_path_list, 
        'pon_id_list': pon_id_list, 
        'refver': refver, 
        'no_matesearch': no_matesearch, 
        'countonly': countonly, 
        'memuse_limit_gb': memuse_limit_gb,
        'depth_limits': depth_limits,
        'mq_limits': mq_limits,
        'include_blacklist': include_blacklist,
        'blacklist_gr': blacklist_gr,
    }

    # main
    while True:
        p = multiprocessing.Process(target=unit_job_core, kwargs=subproc_kwargs)
        p.start()
        while True:
            time.sleep(monitor_interval)
            shareddict['parent_memuse_gb'] = common.get_rss(mode='total', unit='g')
            common.print_timestamp(f"Memory usage: {shareddict['parent_memuse_gb']} GB")
            if not p.is_alive():
                break

        if p.exitcode != 0:
            raise Exception(f'The subprocess finished with nonzero exit code.')
        p.close()

        if shareddict['finished']:
            break
        else:
            common.print_timestamp(
                f"Beginning next cycle with new output file: {shareddict['outfile_path']}"
            )
            continue


def unit_job_core(
    shareddict,
    read_vrspecs,
    infile_copy_path,
    fetchregion,
    bam_path_list,
    id_list,
    pon_bam_path_list,
    pon_id_list,
    refver,
    no_matesearch,
    countonly,
    memuse_limit_gb,
    depth_limits,
    mq_limits,
    include_blacklist,
    blacklist_gr,
):
    # basic setup
    bam_dict = {
        sampleid: pysam.AlignmentFile(bam_path)
        for sampleid, bam_path in zip(id_list, bam_path_list)
    }
    pon_bam_dict = {
        sampleid: pysam.AlignmentFile(bam_path)
        for sampleid, bam_path in zip(pon_id_list, pon_bam_path_list)
    }
    fasta = common.DEFAULT_FASTAS[refver]
    chromdict = common.DEFAULT_CHROMDICTS[refver]

    # edit vcf header
    in_vcf = pysam.VariantFile(infile_copy_path, 'r')
    new_header = in_vcf.header.copy()
    added_new_samples = update_header(new_header, id_list, pon_id_list)
    out_vcf = pysam.VariantFile(shareddict['outfile_path'], 'wz', header=new_header)

    # loop over variant records
    vr_iterator = itertools.chain.from_iterable(
        vcfmisc.get_vr_fetcher(
            in_vcf, refver, fetchargs[0], fetchargs[1], fetchargs[2],
            respect_refregion=False,
        )
        for fetchargs in fetchregion
    )
    shareddict['finished'] = True
    for vr in vr_iterator:
        vrspec = '\t'.join([vr.contig, str(vr.pos), vr.ref, ','.join(vr.alts)])
        if vrspec in read_vrspecs:
            continue

        read_vrspecs.append(vrspec)
        common.print_timestamp(f'Processing {vrspec}')  # for logging

        if added_new_samples:
            new_vr = varianthandler.reheader(vr, new_header)
        else:
            new_vr = vr

        update_new_vr(
            new_vr, fasta, chromdict, bam_dict, pon_bam_dict, no_matesearch, countonly,
            depth_limits, mq_limits, blacklist_gr,
        )
        out_vcf.write(new_vr)

        if shareddict['parent_memuse_gb'] > memuse_limit_gb:
            next_outfile_path = toolsetup.make_next_split_vcf_path(shareddict['outfile_path'])
            shareddict['outfile_path'] = next_outfile_path
            shareddict['finished'] = False
            break

    out_vcf.close()
    in_vcf.close()


##################################


def unit_job_old(
    split_infile_path, 
    split_outfiles_dir, 
    bam_path_list, 
    id_list, 
    pon_bam_path_list, 
    pon_id_list, 
    refver, 
    no_matesearch, 
    countonly, 
    memuse_limit_gb,
    depth_limits,
    mq_limits,
    include_blacklist,
    monitor_interval=10,
):
    # setup manager
    manager = multiprocessing.Manager()
    shareddict = manager.dict()
    shareddict['parent_memuse_gb'] = common.get_rss(mode='total', unit='g')
    common.print_timestamp(f"Memory usage: {shareddict['parent_memuse_gb']} GB")

    # setup paramaters
#    bam_dict = {
#        sampleid: pysam.AlignmentFile(bam_path)
#        for sampleid, bam_path in zip(id_list, bam_path_list)
#    }
#    pon_bam_dict = {
#        sampleid: pysam.AlignmentFile(bam_path)
#        for sampleid, bam_path in zip(pon_id_list, pon_bam_path_list)
#    }
#    fasta = common.DEFAULT_FASTAS[refver]
#    chromdict = common.DEFAULT_CHROMDICTS[refver]

    cytoband_gr = ucscdata.get_cytoband_gr(refver=refver, rename_hg19=True)
    centromere_gr = cytoband_gr[cytoband_gr.Stain == 'acen']
    blacklist_gr = pr.concat([centromere_gr, blacklist.CURATED_BLACKLIST_GRS[refver]])
    blacklist_gr = blacklist_gr.sort().merge()

    # make subprocess arguments
    subproc_kwargs = {
        'split_infile_path': split_infile_path, 
        'split_outfiles_dir': split_outfiles_dir, 
        'bam_path_list': bam_path_list, 
        'id_list': id_list, 
        'pon_bam_path_list': pon_bam_path_list, 
        'pon_id_list': pon_id_list, 
        'refver': refver, 
        'no_matesearch': no_matesearch, 
        'countonly': countonly, 
        'shareddict': shareddict,
        'memuse_limit_gb': memuse_limit_gb,
        'depth_limits': depth_limits,
        'mq_limits': mq_limits,
        'include_blacklist': include_blacklist,
        'blacklist_gr': blacklist_gr,
    }

    # main
    while True:
        p = multiprocessing.Process(target=unit_job_core, kwargs=subproc_kwargs)
        shareddict['next_infile_path'] = None
        p.start()
        while True:
            time.sleep(monitor_interval)
            shareddict['parent_memuse_gb'] = common.get_rss(mode='total', unit='g')
            common.print_timestamp(f"Memory usage: {shareddict['parent_memuse_gb']} GB")
            if not p.is_alive():
                break

        if p.exitcode != 0:
            raise Exception(f'The subprocess finished with nonzero exit code.')
        p.close()

        if shareddict['next_infile_path'] is None:
            break
        else:
            common.print_timestamp(f"Next infile created: {shareddict['next_infile_path']}")
            subproc_kwargs['split_infile_path'] = shareddict['next_infile_path']
            continue


def unit_job_core_old(
    split_infile_path, 
    split_outfiles_dir, 
    bam_path_list, 
    id_list, 
    pon_bam_path_list, 
    pon_id_list, 
    refver, 
    no_matesearch, 
    countonly, 
    shareddict,
    memuse_limit_gb,
    depth_limits,
    mq_limits,
    include_blacklist,
    blacklist_gr,
):
    # basic setup
    bam_dict = {
        sampleid: pysam.AlignmentFile(bam_path)
        for sampleid, bam_path in zip(id_list, bam_path_list)
    }
    pon_bam_dict = {
        sampleid: pysam.AlignmentFile(bam_path)
        for sampleid, bam_path in zip(pon_id_list, pon_bam_path_list)
    }
    fasta = common.DEFAULT_FASTAS[refver]
    chromdict = common.DEFAULT_CHROMDICTS[refver]

    split_outfile_path = os.path.join(
        split_outfiles_dir, 
        os.path.basename(split_infile_path),
    )

    # edit vcf header
    in_vcf = pysam.VariantFile(split_infile_path, 'r')
    new_header = in_vcf.header.copy()
    added_new_samples = update_header(new_header, id_list, pon_id_list)
    out_vcf = pysam.VariantFile(split_outfile_path, 'wz', header=new_header)

    # loop over variant records
    vr_iterator = in_vcf.fetch()
    for vr in vr_iterator:
        vrspec = '\t'.join([vr.contig, str(vr.pos), vr.ref, ','.join(vr.alts)])
        common.print_timestamp(f'Processing {vrspec}')  # for logging

        if added_new_samples:
            new_vr = varianthandler.reheader(vr, new_header)
        else:
            new_vr = vr

        update_new_vr(
            new_vr, fasta, chromdict, bam_dict, pon_bam_dict, no_matesearch, countonly,
            depth_limits, mq_limits, blacklist_gr,
        )
        out_vcf.write(new_vr)

        if shareddict['parent_memuse_gb'] > memuse_limit_gb:
            #next_split_infile_path = re.sub('\.vcf\.gz$', 'a.vcf.gz', split_infile_path)
            next_split_infile_path = toolsetup.make_next_infile(split_infile_path)
            next_out_vcf = pysam.VariantFile(next_split_infile_path, 'wz', header=in_vcf.header.copy())
            shareddict['next_infile_path'] = next_split_infile_path
            for vr in vr_iterator:
                next_out_vcf.write(vr)
            next_out_vcf.close()
            break

    out_vcf.close()
    in_vcf.close()


def update_header(vcfheader, id_list, pon_id_list):
    libreadstats.ReadStatsSampledict.add_meta(vcfheader)
    added_new_samples = False
    for sampleid in itertools.chain(id_list, pon_id_list):
        if sampleid not in vcfheader.samples:
            added_new_samples = True
            vcfheader.samples.add(sampleid)

    return added_new_samples


def update_new_vr(
    new_vr, fasta, chromdict, bam_dict, pon_bam_dict, no_matesearch, countonly,
    depth_limits, mq_limits, blacklist_gr,
):
    vcfspec = libvcfspec.Vcfspec.from_vr(new_vr)
    if vcfspec.to_gr().intersect(blacklist_gr).empty:
        init_invalid = False
    else:
        init_invalid = True

    readstats_dict = libreadstats.ReadStatsSampledict.from_bam_dict(
        bam_dict, vcfspec, fasta, chromdict,
        rpplist_kwargs={'no_matesearch': no_matesearch},
        countonly=countonly,
        depth_limits=depth_limits, 
        mq_limits=mq_limits,
        init_invalid=init_invalid,
    )
    readstats_dict.update_bams(
        bam_dict=pon_bam_dict,
        rpplist_kwargs={'no_matesearch': no_matesearch},
        countonly=True,
        depth_limits=depth_limits, 
        mq_limits=mq_limits,
        init_invalid=init_invalid,
    )
    readstats_dict.write(new_vr)

    del readstats_dict
    gc.collect()

    # get irrelevant sample ids: these need to be set '.' afterward
    #FORMAT_key = libreadstats.ReadStatsSampledict.annotkey
#    irrelevant_samples = set(
#        x for x in new_vr.samples.keys()
#        if infoformat.check_NA_format(new_vr, x, FORMAT_key)
#    )
#    irrelevant_samples.difference_update(id_list)


    #for sampleid in irrelevant_samples:
    #    infoformat.set_NA_format(new_vr, sampleid, FORMAT_key)


#########################################################


def argument_parser(cmdargs):
    def sanity_check(args):
        if all(
            [
                (args.bamlist is None),
                (args.bamlist_file_path is None),
                (args.pon_cohorts is None),
            ]
        ):
            raise Exception(f'You must use at least one of "--bamlist", "--bamlist-file", or "--pon-cohorts".')

        if (args.bamlist is None) != (args.idlist is None):
            raise Exception(f'--bamlist and --idlist options must be used together.')

        # --bamlist and --idlist length comparison
        if (args.bamlist is not None) and (args.idlist is not None):
            if len(args.bamlist) != len(args.idlist):
                raise Exception(f'The number of --bamlist and --idlist arguments must be the same.')

    def bam_path_processing(args):
        # all information specified by --bamlist/--idlist, --bamlist-file, --pon-cohorts are added
        # handling --bamlist-file
        if args.bamlist_file_path is not None:
            with open(args.bamlist_file_path, 'r') as infile:
                for line in infile:
                    linesp = line.replace('\n', '').split('\t')
                    if len(linesp) != 2:
                        raise Exception(f'The number of tab-seperated fields must be 2 in the bamlist file.')

                    sampleid = linesp[0]
                    bam_path = linesp[1]
                    if not bam_path.endswith('.bam'):
                        raise Exception(
                            f'Bam file path does not end with ".bam" in the '
                            f'bamlist file. Please check if column 1 is '
                            f'sample ID and column 2 is bam file path.'
                        )

                    args.idlist.append(sampleid)
                    args.bamlist.append(bam_path)
        # handling --pon-cohorts
        args.pon_idlist = list()
        args.pon_bamlist = list()
        if args.pon_cohorts is not None:
            for cohort in args.pon_cohorts:
                for sampleid, bam_path in libponbams.PON_BAM_PATHS[args.refver][cohort].items():
                    args.pon_idlist.append(sampleid)
                    args.pon_bamlist.append(bam_path)

    parser_dict = workflow.init_parser(
        description=(
            f'Calculates allele-supporting read count information for each '
            f'input bam, then writes to FORMAT. To specify input bam '
            f'files and sample IDs, you can use any combination of:\n'
            f'1) --bamlist and --idlist ; '
            f'2) --bamlist-file ; '
            f'3) --pon-cohorts\n'
            f'All these bam files will be used for annotation.\n'
            f'Bam files specified with "--pon-cohorts" will be used with '
            f'"countonly" mode, regardless of the usage of "--countonly" argument.'
        )
    )

    # required
    workflow.add_infile_arg(parser_dict['required'], required=True)
    workflow.add_outfile_arg(parser_dict['required'], required=True, must_not_exist='ask')
    workflow.add_refver_arg(parser_dict['required'], required=True)

    # optional
    parser_dict['optional'].add_argument(
        '--bamlist', 
        dest='bamlist', 
        required=False,
        nargs='+', 
        type=workflow.arghandler_infile, 
        default=list(),
        #metavar='<Input bam file path>',
        help=f'One or more input bam file paths separated by whitespaces.',
    )
    parser_dict['optional'].add_argument(
        '--idlist', 
        dest='idlist', 
        required=False,
        nargs='+', 
        default=list(),
        #metavar='<Sample ID>',
        help=(
            f'Sample IDs of the input bam files, in the same order, '
            f'separated by whitespaces.'
        ),
    )
    parser_dict['optional'].add_argument(
        '--bamlist-file', dest='bamlist_file_path', required=False,
        type=workflow.arghandler_infile,
        #metavar='<bam list file path>',
        help=(
            f'A 2-column tab-separated file which contains sample IDs on '
            f'the first column and bam file paths on the second column.'
        ),
    )

    allowed_pon_cohorts = dict()
    for refver, subdict in libponbams.PON_BAM_PATHS_WITHOUT_NAMES.items():
        allowed_pon_cohorts[refver] = list(subdict.keys())
    parser_dict['optional'].add_argument(
        '--pon-cohorts', dest='pon_cohorts', required=False,
        nargs='+',
        #metavar='<PON cohort name>',
        help=(
            f'One or more known Panel Of Normal cohort names. Allowed values:\n'
            f'{textwrap.indent(pprint.pformat(allowed_pon_cohorts), " " * 4)}'
        ),
    )

    parser_dict['optional'].add_argument(
        '--memlimit', dest='memuse_limit_gb', required=False,
        default=4.0,
        type=float,
        #metavar='<memory limit per job>',
        help=f'Maximum memory limit per job in gigabytes.',
    )

    parser_dict['optional'].add_argument(
        '--depth-limits', dest='depth_limits', required=False,
        default=None,
        #type=float,
        #metavar='<depth upper limit>',
        nargs='+',
        help=(
            f'Limits (both inclusive) of valid depth range. '
            f'If depth of certain site is out of this range, ReadStats '
            f'information will be all set as NaN. -1 given as the second value '
            f'is interpreted as positive infinity. '
            f'Default: lower is 0, upper is '
            f'(--depthlimit-cov-ratio value) * (average depth of the bam file)'
        ),
    )
    parser_dict['optional'].add_argument(
        '--depthlimit-cov-ratio', dest='depthlimit_cov_ratio', required=False,
        default=4,
        type=float,
        #metavar='<bam coverage multiplication factor>',
        help=(
            f'Upper depth limit is set as (this value) * (average bam depth).'
            f'If "--depth-limits" option is used, this option is ignored.'
        ),
    )
    parser_dict['optional'].add_argument(
        '--target-bed', dest='target_bed_path', required=False,
        default=None,
        #metavar='<target region bed file>',
        help=(
            f'If this option is used, bam file average depth calculation will '
            f'be based on the region length derived from this file. '
            f'If "--depth-limits" option is used, this option is ignored.'
        ),
    )
    parser_dict['optional'].add_argument(
        '--mq-limits', dest='mq_limits', required=False,
        default=[
            ','.join(str(x) for x in libreadstats.DEFAULT_MQ_LIMITS)
        ],
        #type=float,
        #metavar='<MQ limits>',
        nargs='+',
        help=(
            f'Limits (both inclusive) of valid mapping quality range. '
            f'If mapping quality of certain site is out of this range, ReadStats '
            f'information will be all set as NaN. -1 given as the second value '
            f'is interpreted as positive infinity.'
        ),
    )

    workflow.add_outfmt_arg(parser_dict['optional'], required=False)
    workflow.add_logging_args(parser_dict)
    workflow.add_scheduler_args(parser_dict, default_parallel=1, default_sched='slurm')

    # flag
    workflow.add_rmtmp_arg(parser_dict)
    workflow.add_index_arg(parser_dict)
    parser_dict['flag'].add_argument(
        '--do-matesearch', dest='do_matesearch', action='store_true',
        help=(f'If set, "no_matesearch" argument is set as False when '
              f'creating ReadPlusPairList objects. It can be much slower, '
              f'but mate reads far away from the variant site can be captured.')
    )
    parser_dict['flag'].add_argument(
        '--countonly', dest='countonly', action='store_true',
        help=f'If set, only read count informations is calculated. Suitable for PON data.'
    )
    parser_dict['flag'].add_argument(
        '--include-blacklist', dest='include_blacklist', action='store_true',
        help=(
            f'By default, variants within blacklist region '
            f'(centromere + curated blacklist regions) '
            f'will be annotated as NaN. '
            f'If this option is set, blacklist region variants will be fully annotated.'
        ),
    )

    # main
    args = parser_dict['main'].parse_args(cmdargs)
    sanity_check(args)
    bam_path_processing(args)

    return args


def depth_mq_limit_processing(args):
    def rawarg_handler(rawarg, nsample):
        if len(rawarg) not in (1, nsample):
            raise Exception(f'Invalid "--depth-limits" or "--mq-limits" arguments')

        if len(rawarg) == 1:
            rawarg = rawarg * nsample
            
        result = list()
        for x in rawarg:
            x = x.split(',')
            minval = float(x[0])
            maxval = float(x[1])
            result.append([minval, maxval])

        return result

    # calculate region length
    if args.depth_limits is None:
        if args.target_bed_path is None:
            sample_bam_path = next(itertools.chain(args.bamlist, args.pon_bamlist))
            with pysam.AlignmentFile(sample_bam_path) as in_bam:
                region_length = sum(in_bam.lengths)
        else:
            region_gr = pr.read_bed(args.target_bed_path)
            region_length = region_gr.merge().length

    # modify raw argument
    sampleid_list = list(itertools.chain(args.idlist, args.pon_idlist))
    bampath_list = list(itertools.chain(args.bamlist, args.pon_bamlist))
    nsample = len(sampleid_list)

    mq_limits = dict(zip(sampleid_list, rawarg_handler(args.mq_limits, nsample)))

    if args.depth_limits is None:
        depth_limits = dict()
        for sampleid, bam_path in zip(sampleid_list, bampath_list):
            with pysam.AlignmentFile(bam_path) as in_bam:
                upper = (
                    bameditor.get_average_depth(in_bam, aligned_region_length=region_length) 
                    * args.depthlimit_cov_ratio
                )
            depth_limits[sampleid] = [0, upper]
    else:
        depth_limits = dict(zip(sampleid_list, rawarg_handler(args.depth_limits, nsample)))

#    mq_limits = dict()
#    depth_limits = dict()
#    for sampleid, bam_path in zip(
#        itertools.chain(args.idlist, args.pon_idlist),
#        itertools.chain(args.bamlist, args.pon_bamlist),
#    ):
#        mq_limits[sampleid] = list(args.mq_limits)
#        if args.depth_limits is None:
#            with pysam.AlignmentFile(bam_path) as in_bam:
#                upper = (
#                    bameditor.get_average_depth(in_bam, aligned_region_length=region_length) 
#                    * args.depthlimit_cov_ratio
#                )
#                depth_limits[sampleid] = [0, upper]
#        else:
#            depth_limits[sampleid] = list(args.depth_limits)

    return depth_limits, mq_limits


def write_jobscripts(
    tmpdir_paths, 
    fetchregion_list,
    parallel, 
    infile_copy_path,
    bam_path_list, 
    id_list, 
    pon_bam_path_list, 
    pon_id_list, 
    refver, 
    no_matesearch,
    countonly,
    memuse_limit_gb,
    depth_limits,
    mq_limits,
    include_blacklist,
    jobname_prefix=__name__.split('.')[-1], 
    nproc=1,
):
    jobscript_path_list = list()
    log_path_list = list()
    slurm_log_pf_list = list()
    split_outfile_path_list = list()
    for zidx in common.zrange(parallel):
        jobscript_path = os.path.join(tmpdir_paths['scripts'], f'{zidx}.sbatch')
        jobscript_path_list.append(jobscript_path)

        bname = os.path.basename(jobscript_path)
        log_path_list.append(
            os.path.join(tmpdir_paths['logs'], bname + '.log')
        )
        slurm_log_pf_list.append(
            os.path.join(tmpdir_paths['slurm_logs'], bname)
        )

        split_outfile_path = os.path.join(tmpdir_paths['split_outfiles'], f'{zidx}.vcf.gz')
        split_outfile_path_list.append(split_outfile_path)

    kwargs_single = {
        'infile_copy_path': infile_copy_path,
        'bam_path_list': bam_path_list,
        'id_list': id_list,
        'pon_bam_path_list': pon_bam_path_list,
        'pon_id_list': pon_id_list,
        'refver': refver,
        'no_matesearch': no_matesearch,
        'countonly': countonly,
        'memuse_limit_gb': memuse_limit_gb,
        'depth_limits': depth_limits,
        'mq_limits': mq_limits,
        'include_blacklist': include_blacklist,
    }
    kwargs_multi = {
        'split_outfile_path': split_outfile_path_list,
        'fetchregion': fetchregion_list,
    }

    toolsetup.write_jobscripts(
        script_path_list=jobscript_path_list,
        log_path_list=log_path_list,
        slurm_log_pf_list=slurm_log_pf_list,
        module_name=__name__,
        unit_job_func_name='unit_job',
        kwargs_single=kwargs_single,
        kwargs_multi=kwargs_multi,
        jobname_prefix=jobname_prefix,
        nproc=nproc,
    )

    return jobscript_path_list


#def write_jobscripts_old(
#    tmpdir_paths, 
#    split_infile_path_list, 
#    bam_path_list, 
#    id_list, 
#    pon_bam_path_list, 
#    pon_id_list, 
#    refver, 
#    no_matesearch,
#    countonly,
#    memuse_limit_gb,
#    depth_limits,
#    mq_limits,
#    include_blacklist,
#    jobname_prefix=__name__.split('.')[-1], 
#    nproc=1,
#):
#    jobscript_path_list = list()
#    log_path_list = list()
#    for zidx in common.zrange(len(split_infile_path_list)):
#        jobscript_path = os.path.join(tmpdir_paths['scripts'], f'{zidx}.sbatch')
#        jobscript_path_list.append(jobscript_path)
#        log_path = os.path.join(tmpdir_paths['logs'], os.path.basename(jobscript_path) + '.log')
#        log_path_list.append(log_path)
#
#    kwargs_single = {
#        'split_outfiles_dir': tmpdir_paths['split_outfiles'],
#        'bam_path_list': bam_path_list,
#        'id_list': id_list,
#        'pon_bam_path_list': pon_bam_path_list,
#        'pon_id_list': pon_id_list,
#        'refver': refver,
#        'no_matesearch': no_matesearch,
#        'countonly': countonly,
#        'memuse_limit_gb': memuse_limit_gb,
#        'depth_limits': depth_limits,
#        'mq_limits': mq_limits,
#        'include_blacklist': include_blacklist,
#    }
#    kwargs_multi = {
#        'split_infile_path': split_infile_path_list,
#    }
#
#    toolsetup.write_jobscripts(
#        script_path_list=jobscript_path_list,
#        log_path_list=log_path_list,
#        module_name=__name__,
#        unit_job_func_name='unit_job',
#        kwargs_single=kwargs_single,
#        kwargs_multi=kwargs_multi,
#        jobname_prefix=jobname_prefix,
#        nproc=nproc,
#    )
#
#    return jobscript_path_list


def main(cmdargs):
    args = argument_parser(cmdargs)

    # make tmpdir tree
    tmpdir_paths = workflow.get_tmpdir_paths(
        ['scripts', 'logs', 'slurm_logs', 'split_outfiles'],
        prefix = (
            os.path.basename(args.infile_path)
            + '_'
            + __name__.split('.')[-1]
            + '_'
        ),
        where=os.path.dirname(args.infile_path),
    )

    # setup logger
    logger = toolsetup.setup_logger(
        args=args,
        tmpdir_root=tmpdir_paths['root'],
        with_genlog=True,
    )
    logger.info('Beginning')

    # postprocess depth and mq limits
    logger.info(f'Processing depth and MQ limits (may take a while calculating average depths of bam files)')
    depth_limits, mq_limits = depth_mq_limit_processing(args)
    logger.info(f'Depth limits: {depth_limits}')
    logger.info(f'MQ limits: {mq_limits}')

    # make infile copy
    infile_copy_path, infile_copy_index_path = toolsetup.make_infile_copy(
        args.infile_path, tmpdir_paths['root'], logger,
    )

    # make split fetchregions
    logger.info(f'Getting split fetch regions of input VCF')
    fetchregion_list = vcfmisc.get_vcf_fetchregions(
        infile_copy_path, args.parallel, args.refver, verbose=False,
    )
    #logger.info(f'Split fetch regions: {fetchregion_list}')

    # make job scripts and run
    logger.info(f'Creating split job scripts')
    jobscript_path_list = write_jobscripts(
        tmpdir_paths, 
        fetchregion_list,
        args.parallel,
        infile_copy_path,
        args.bamlist, 
        args.idlist, 
        args.pon_bamlist, 
        args.pon_idlist, 
        args.refver, 
        (not args.do_matesearch),
        args.countonly,
        args.memuse_limit_gb,
        depth_limits,
        mq_limits,
        args.include_blacklist,
    )
    logger.info('Running annotation jobs for each split file')
    workflow.run_jobs(
        jobscript_path_list, 
        sched=args.sched, 
        intv_check=args.intv_check, 
        intv_submit=args.intv_submit, 
        max_submit=args.max_submit, 
        logger=logger, 
        log_dir=tmpdir_paths['logs'],
        job_status_logpath=os.path.join(tmpdir_paths['root'], 'job_status_log.gz'),
        raise_on_failure=True,
    )

    # concatenates split files
    logger.info('Merging split files')
    outfile_path_list = common.listdir(tmpdir_paths['split_outfiles'])
    libconcat.main(
        infile_path_list=outfile_path_list,
        outfile_path=args.outfile_path, 
        mode_pysam=args.mode_pysam,
        outfile_must_not_exist='no',
    )

    # rmtmp, index
    if not args.donot_rm_tmp:
        shutil.rmtree(tmpdir_paths['root'])
    if not args.donot_index:
        indexing.index_vcf(args.outfile_path)

    logger.info('All successfully finished')


#def main_old(cmdargs):
#    args = argument_parser(cmdargs)
#
#    # postprocess depth and mq limits
#    common.print_timestamp(f'Processing depth and MQ limits (may take a while calculating average depths of bam files)')
#    depth_limits, mq_limits = depth_mq_limit_processing(args)
#    common.print_timestamp(f'Depth limits: {depth_limits}')
#    common.print_timestamp(f'MQ limits: {mq_limits}')
#
#    # make tmpdir tree
#    tmpdir_paths = workflow.get_tmpdir_paths(
#        ['scripts', 'logs', 'split_infiles', 'split_outfiles'],
#        prefix = (
#            os.path.basename(args.infile_path)
#            + '_'
#            + __name__.split('.')[-1]
#            + '_'
#        ),
#        where=os.path.dirname(args.infile_path),
#    )
#
#    # setup logger
#    logger = toolsetup.setup_logger(args=args,
#                                    tmpdir_root=tmpdir_paths['root'],
#                                    with_genlog=True)
#    logger.info('Beginning')
#
#    # split the input file
#    logger.info('Splitting the input file')
#    split_infile_path_list = libsplit.main(vcf_path=args.infile_path, 
#                                        outdir=tmpdir_paths['split_infiles'], 
#                                        n_file=args.parallel, 
#                                        mode_bcftools='z', prefix='', 
#                                        suffix='.vcf.gz')
#
#    # make job scripts and run
#    jobscript_path_list = write_jobscripts(
#        tmpdir_paths, 
#        split_infile_path_list, 
#        args.bamlist, 
#        args.idlist, 
#        args.pon_bamlist, 
#        args.pon_idlist, 
#        args.refver, 
#        (not args.do_matesearch),
#        args.countonly,
#        args.memuse_limit_gb,
#        depth_limits,
#        mq_limits,
#        args.include_blacklist,
#    )
#    logger.info('Running annotation jobs for each split file')
#    workflow.run_jobs(
#        jobscript_path_list, 
#        sched=args.sched, 
#        intv_check=args.intv_check, 
#        intv_submit=args.intv_submit, 
#        max_submit=args.max_submit, 
#        logger=logger, 
#        log_dir=tmpdir_paths['logs'],
#        job_status_logpath=os.path.join(tmpdir_paths['root'], 'job_status_log.gz'),
#        raise_on_failure=True,
#    )
#
#    # concatenates split files
#    logger.info('Merging split files')
#    outfile_path_list = common.listdir(tmpdir_paths['split_outfiles'])
#    print(outfile_path_list)
#    libconcat.main(
#        infile_path_list=outfile_path_list,
#        outfile_path=args.outfile_path, 
#        mode_pysam=args.mode_pysam,
#        outfile_must_not_exist='no',
#    )
#
#    # rmtmp, index
#    if not args.donot_rm_tmp:
#        shutil.rmtree(tmpdir_paths['root'])
#    if not args.donot_index:
#        indexing.index_vcf(args.outfile_path)
#
#    logger.info('All successfully finished')

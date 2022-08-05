import os
import argparse
import subprocess
import tempfile

import pysam

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
workflow = importlib.import_module(".".join([top_package_name, "workflow"]))
libgatk = importlib.import_module(".".join([top_package_name, "gatk"]))


PROGNAME = __name__.split(".")[-1]


def unit_job(split_outfile_path, fasta_path, region_path, padded_region_path, mbq):
    # run HaplotypeCaller
    tmpdir = tempfile.TemporaryDirectory(where=os.getcwd())
    split_outfile_path_prefilter = split_outfile_path + ".prefilter"
    libgatk.run_haplotypecaller(
        fasta_path=fasta_path,
        infile_path=infile_path,
        outfile_path=split_outfile_path_prefilter,
        tmpdir=tmpdir,
        incl_bed_path=padded_region_path,
        mbq=mbq,
    )

    # filter
    region_intvlist = common.IntervalList.from_bed(region_path)
    in_vcf = pysam.VariantFile(split_outfile_path_prefilter, "r")
    out_vcf = pysam.VariantFile(split_outfile_path, "wz", header=in_vcf.header.copy())

    for vr in in_vcf.fetch():
        if region_intvlist.includes_vr(vr):
            out_vcf.write(vr)


def argument_parser(cmdargs):
    parser_dict = workflow.init_parser(
        description=(f"Runs gatk HaplotypeCaller with N parallel jobs.")
    )
    # required
    workflow.add_infile_arg(
        parser_dict["required"], required=True, help=f"Input bam file path."
    )
    workflow.add_outfile_arg(
        parser_dict["required"], required=True, must_not_exist="ask"
    )
    workflow.add_fasta_arg(parser_dict["required"], required=True)
    # optional
    workflow.add_outfmt_arg(parser_dict["optional"], required=False)
    workflow.add_incl_region_arg(parser_dict["optional"], required=False)
    workflow.add_excl_region_arg(parser_dict["optional"], required=False)
    workflow.add_logging_args(parser_dict)
    workflow.add_scheduler_args(parser_dict, default_parallel=1, default_sched="slurm")
    # flag
    workflow.add_rmtmp_arg(parser_dict)
    workflow.add_index_arg(parser_dict)


def make_tmpdir(infile_path):
    tmpdir_paths = workflow.get_tmpdir_paths(
        ["scripts", "logs", "regions", "split_results"],
        prefix=(os.path.basename(infile_path) + "_" + PROGNAME),
        where=os.path.dirname(infile_path),
    )

    return tmpdir_paths


def get_included_intvlist(fasta_path, incl_bed_path, excl_bed_path, chromdict):
    if incl_bed_path is None:
        incl_intvlist = common.IntervalList.from_chromdict(
            common.ChromDict(fasta_path=fasta_path)
        )
    else:
        incl_intvlist = common.IntervalList.from_bed(incl_bed_path)

    if excl_bed_path is not None:
        excl_intvlist = common.IntervalList.from_bed(excl_bed_path)
        incl_intvlist = incl_intvlist.subtract(excl_intvlist)

    incl_intvlist.sort_intervals(chromdict)

    return incl_intvlist


def split_intvlist(incl_intvlist, parallel):
    incl_intvlist_split = incl_intvlist.split(num=parallel)
    incl_intvlist_split_padded = incl_intvlist_split.slop(chromdict, b=5000)

    return incl_intvlist_split, incl_intvlist_split_padded


def write_region_files(
    regionfiles_dir, incl_intvlist_split, incl_intvlist_split_padded
):
    for zidx, intvlist in common.zenumerate(incl_intvlist_split):
        fname = os.path.join(regionfiles_dir, f"{zidx}.bed")
        intvlist.write_bed(fname)

    for zidx, intvlist in common.zenumerate(incl_intvlist_split_padded):
        fname = os.path.join(regionfiles_dir, f"{zidx}.padded.bed")
        intvlist.write_bed(fname)


def main(cmdargs):
    args = argument_parser(cmdargs)

    # make temporary directory tree
    tmpdir_paths = make_tmpdir(args.infile_path)

    # setup logger
    logger = toolsetup.setup_logger(
        args=args, tmpdir_root=tmpdir_paths["root"], with_genlog=True
    )
    logger.info("Beginning")

    # setup other parameters
    chromdict = common.ChromDict(fasta_path=args.fasta_path)

    # make split intervals and write bed files
    incl_intvlist = get_included_intvlist(
        args.fasta_path, args.incl_bed_path, args.excl_bed_path, chromdict
    )
    incl_intvlist_split, incl_intvlist_split_padded = split_intvlist(
        incl_intvlist, args.parallel
    )
    write_region_files(
        tmpdir_paths["regions"], incl_intvlist_split, incl_intvlist_split_padded
    )

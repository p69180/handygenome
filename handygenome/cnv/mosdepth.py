import os
import tempfile
import shutil
import subprocess

import pyranges as pr
import pandas as pd

import handygenome.common as common
import handygenome.bameditor as bameditor
import handygenome.cnv.misc as cnvmisc


MOSDEPTH_PATH = '/home/users/pjh/tools/miniconda/221104/miniconda3/envs/mosdepth/bin/mosdepth'


def get_mosdepth_args(prefix, bam_path, t, use_median, no_perbase=True, bed_path=None, window_size=None):
    if (
        (bed_path is not None)
        and (window_size is not None)
    ):
        raise Exception(f'"window_size" and "bed_path" arguments must not be used at the same time.')

    args = [MOSDEPTH_PATH, '-x', '-t', str(t)]

    if use_median:
        args.append('-m')

    if no_perbase:
        args.append('-n')

    if bed_path is not None:
        args.extend(('-b', bed_path))
    elif window_size is not None:
        args.extend(('-b', str(window_size)))

    args.extend([prefix, bam_path])

    return args


def load_mosdepth_output(filename, depth_colname='mean_depth', as_gr=True):
    """Only for *.regions.bed.gz file"""
    df = pd.read_csv(
        filename, 
        sep='\t', 
        names=['Chromosome', 'Start', 'End', depth_colname], 
        dtype={'Chromosome': str, 'Start': int, 'End': int, depth_colname: float},
    )
    if as_gr:
        return pr.PyRanges(df)
    else:
        return df


def run_mosdepth(bam_path, t=8, use_median=False, region_bed_path=None, region_gr=None, window_size=None, donot_subset_bam=False, as_gr=True, load_perbase=False):
    # sanity check
    if (
        (window_size is not None)
        and (
            (region_bed_path is not None)
            or (region_gr is not None)
        )
    ):  
        raise Exception(f'"window_size" and "region*" arguments must not be used at the same time.')

    if (
        (region_bed_path is not None)
        and (region_gr is not None)
    ):  
        raise Exception(f'Only one of "region_bed_path" and "region_gr" argument must be set.')

    # handle region_gr argument
    region_df = cnvmisc.arg_into_df(region_gr)
    cnvmisc.genome_df_sanitycheck(region_df)

    # make tmp directory
    tmpdir = tempfile.mkdtemp(dir=os.getcwd(), prefix='mosdepth_tmpdir_')

    # write region bed file
    if (region_bed_path is None) and (region_df is None):
        mosdepth_input_bed_path = None
    else:
        if region_df is not None:
            unedited_input_bed_path = os.path.join(tmpdir, 'unedited_region.bed')
            region_df.to_csv(unedited_input_bed_path, sep='\t', header=False, index=False)
        elif region_bed_path is not None:
            unedited_input_bed_path = region_bed_path

        mosdepth_input_bed_path = os.path.join(tmpdir, 'region.bed')
        with common.openfile(unedited_input_bed_path, 'r') as infile:
            with open(mosdepth_input_bed_path, 'wt') as outfile:
                for line in infile:
                    linesp = line.strip().split('\t')
                    if not (linesp[1].isdigit() and linesp[2].isdigit()):
                        continue
                    outfile.write('\t'.join(linesp[:3]) + '\n')

    if mosdepth_input_bed_path is None:
        mosdepth_input_bed_gr = None
    else:
        mosdepth_input_bed_gr = pr.PyRanges(pr.read_bed(mosdepth_input_bed_path, as_df=True), int64=False)

    # subset input bam file
    if mosdepth_input_bed_path is not None:
        if donot_subset_bam:
            mosdepth_input_bam_path = bam_path
        else:
            mosdepth_input_bam_path = os.path.join(tmpdir, 'mosdepth_input.bam')
            bameditor.samtools_view(
                in_bam_path=bam_path, 
                out_bam_path=mosdepth_input_bam_path, 
                region_bed_path=mosdepth_input_bed_path, 
                index=True,
            )
    else:
        mosdepth_input_bam_path = bam_path

    # make mosdepth outdir
    outdir = os.path.join(tmpdir, 'outdir')
    os.mkdir(outdir)

    # run mosdepth
    prefix = 'PREFIX'
    mosdepth_args = get_mosdepth_args(
        prefix=os.path.join(outdir, prefix), 
        bam_path=mosdepth_input_bam_path, 
        t=t, 
        use_median=use_median, 
        no_perbase=(not load_perbase),
        bed_path=mosdepth_input_bed_path, 
        window_size=window_size,
    )
    p = subprocess.run(mosdepth_args, capture_output=True, text=True, check=False)
    if p.returncode != 0:
        print(f'stdout: {p.stdout}')
        print(f'stderr: {p.stderr}')
        p.check_returncode()

    # load mosdepth outputs
    outfile_path = os.path.join(outdir, f'{prefix}.regions.bed.gz')
    if use_median:
        depth_colname = 'median_depth'
    else:
        depth_colname = 'mean_depth'
    df = load_mosdepth_output(outfile_path, depth_colname=depth_colname, as_gr=False)

    if load_perbase:
        perbase_outfile_path = os.path.join(outdir, f'{prefix}.per-base.bed.gz')
        df_perbase = pd.read_csv(
            perbase_outfile_path, 
            sep='\t', 
            names=['Chromosome', 'Start', 'End', 'depth'], 
            dtype={'Chromosome': str, 'Start': int, 'End': int, 'depth': int},
        )

        if mosdepth_input_bed_gr is not None:
            df_perbase = mosdepth_input_bed_gr.window(1).join(
                pr.PyRanges(df_perbase, int64=False),
                apply_strand_suffix=False,
            )[['depth']].df
    else:
        df_perbase = None

    # prepare result
    if as_gr:
        result = pr.PyRanges(df, int64=False)
        if load_perbase:
            result_perbase = pr.PyRanges(df_perbase, int64=False)
        else:
            result_perbase = None
    else:
        result = df
        result_perbase = df_perbase

    # remove tmpdir
    shutil.rmtree(tmpdir)

    return result, result_perbase



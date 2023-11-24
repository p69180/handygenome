import os
import tempfile
import shutil
import subprocess
import multiprocessing

import pandas as pd

import handygenome
import handygenome.deco as deco
import handygenome.refgenome.refgenome as refgenome
import handygenome.logutils as logutils
from handygenome.genomedf.genomedf import GenomeDataFrame as GDF
from handygenome.cnv.depth import DepthRawDataFrame as DepthRawDF
import handygenome.tools as tools
import handygenome.bameditor as bameditor


def get_mosdepth_path():
    if handygenome.PARAMS['mosdepth'] is None:
        raise Exception(f'mosdepth executable cannot be found.')
    return handygenome.PARAMS['mosdepth']


def get_mosdepth_args(
    prefix, bam_path, t, use_median, no_perbase=True, bed_path=None, window_size=None,
):
    if (
        (bed_path is not None)
        and (window_size is not None)
    ):
        raise Exception(f'"window_size" and "bed_path" arguments must not be used at the same time.')

    args = [get_mosdepth_path(), '-x', '-t', str(t)]

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


@deco.get_deco_nproc_limit(3)
@deco.get_deco_num_set_differently(('region_bed_path', 'region_gdf'), 2, how='lt')
def run_mosdepth(
    bam_path, 
    t=1, 
    nproc=1,
    split_width=500,
    use_median=False, 

    region_bed_path=None, 
    region_gdf=None, 

    window=None,

    verbose=True,
    prefix='PREFIX'
):
    # set refver
    refver = refgenome.infer_refver_bampath(bam_path)

    # make bam index
    bam_idx_path = bameditor.get_index_path(bam_path)
    if bam_idx_path is None:
        _ = pysam.index(bam_path)

    region_gdf_list = run_mosdepth_prepare_region_gdfs(
        region_bed_path=region_bed_path,
        region_gdf=region_gdf,
        window=window,
        refver=refver,
        split_width=split_width,
        nosplit=(nproc == 1),
    )

    # run parallel jobs
    tmpdir = tempfile.mkdtemp(dir=os.getcwd(), prefix='mosdepth_tmpdir_')
    bam_basename = os.path.basename(bam_path)
    if verbose:
        logutils.log(f'Running mosdepth parallel jobs (bam: {bam_basename})', add_locstring=False)

    args = (
        (
            bam_path, 
            region_gdf,
            refver,
            t, 
            use_median, 
            prefix,
            False,  # "verbose" argument
            tmpdir,
        )
        for region_gdf in region_gdf_list
    )

    with multiprocessing.Pool(nproc) as pool:
        pool_result = pool.starmap(run_mosdepth_base, args)

    shutil.rmtree(tmpdir)

    if verbose:
        logutils.log(f'Finished mosdepth parallel jobs (bam: {bam_basename})', add_locstring=False)

    result = DepthRawDF.concat(pool_result)

    return result


def run_mosdepth_prepare_region_gdfs(
    region_bed_path,
    region_gdf,
    window,
    refver,
    split_width,
    nosplit=False,
):
    region_is_given = (region_bed_path is not None) or (region_gdf is not None)
    if region_is_given:
        if region_bed_path is not None:
            region_gdf = GDF.read_tsv(region_bed_path, refver)
    else:
        region_gdf = GDF.all_regions(refver, assembled_only=False)

    if window is not None:
        region_gdf = region_gdf.window(window)

    region_gdf.sort()
    if nosplit:
        region_gdf_list = [region_gdf]
    else:
        region_gdf_list = region_gdf.equal_nrow_split(width=split_width)
    return region_gdf_list


#@deco.get_deco_num_set_differently(('region_bed_path', 'region_gdf'), 1)
def run_mosdepth_base(
    bam_path, 
    region_gdf,
    refver,
    t=8, 
    use_median=False, 
    prefix='PREFIX',
    verbose=False,
    tmpdir_location=None,
):
    """- Input bam must be indexed.
    - region bed file: May be gzipped. May not have header row.
    """
    # make tmp directory
    if tmpdir_location is None:
        tmpdir_location = os.getcwd()
    tmpdir = tempfile.mkdtemp(dir=tmpdir_location, prefix='mosdepth_tmpdir_')

    # write region bed file
    if region_gdf is None:
        mosdepth_input_bed_path = None
    else:
        mosdepth_input_bed_path = os.path.join(tmpdir, 'region.bed')
        region_gdf.drop_annots().df.to_csv(
            mosdepth_input_bed_path, 
            sep='\t', 
            header=False, 
            index=False,
        )

    # prepare mosdepth input bam
    mosdepth_input_bam_path = os.path.join(tmpdir, 'mosdepth_input.bam')
    os.symlink(bam_path, mosdepth_input_bam_path)

    # make index of the symlink
    bam_idx_path = bameditor.get_index_path(bam_path)
    if bam_idx_path is None:
        _ = pysam.index(mosdepth_input_bam_path)
    else:
        os.symlink(bam_idx_path, mosdepth_input_bam_path + '.bai')

    # make mosdepth outdir
    outdir = os.path.join(tmpdir, 'outdir')
    os.mkdir(outdir)

    # run mosdepth
    bam_basename = os.path.basename(bam_path)
    if verbose:
        logutils.log(f'Running mosdepth (bam: {bam_basename})', add_locstring=False)
    mosdepth_args = get_mosdepth_args(
        prefix=os.path.join(outdir, prefix), 
        bam_path=mosdepth_input_bam_path, 
        t=t, 
        use_median=use_median, 
        no_perbase=True,
        bed_path=mosdepth_input_bed_path, 
    )
    p = subprocess.check_call(mosdepth_args)
    if verbose:
        logutils.log(f'Finished running mosdepth (bam: {bam_basename})', add_locstring=False)

    # load mosdepth outputs
    outfile_path = os.path.join(outdir, f'{prefix}.regions.bed.gz')
    result_gdf = DepthRawDF.load_mosdepth(outfile_path, refver, use_median=use_median)

    # remove tmpdir
    shutil.rmtree(tmpdir)

    # return
    return result_gdf


def run_mosdepth_old(
    bam_path, 
    refver=None,
    t=8, 
    use_median=False, 
    region_bed_path=None, 
    region_gdf=None, 
    window_size=None, 
    #donot_subset_bam=False, 
    subset_bam=None, 
    load_perbase=False,
    prefix='PREFIX'
):
    """- Input bam must be indexed.
    - region bed file: May be gzipped. May not have header row.
    """
    # sanity check
    if (
        (window_size is not None)
        and (
            (region_bed_path is not None)
            or (region_gdf is not None)
        )
    ):
        raise Exception(f'"window_size" and "region*" arguments must not be used at the same time.')

    if (
        (region_bed_path is not None)
        and (region_gdf is not None)
    ):  
        raise Exception(f'Only one of "region_bed_path" and "region_gdf" argument must be set.')

    # make tmp directory
    tmpdir = tempfile.mkdtemp(dir=os.getcwd(), prefix='mosdepth_tmpdir_')

    # set refver
    if refver is None:
        refver = refgenome.infer_refver_bampath(bam_path)

    # write region bed file
    region_is_given = (region_bed_path is not None) or (region_gdf is not None)
    if region_is_given:
        mosdepth_input_bed_path = os.path.join(tmpdir, 'region.bed')
        if region_bed_path is not None:
            region_gdf = GDF.read_tsv(region_bed_path, refver)
        region_gdf.df.iloc[:, :3].to_csv(
            mosdepth_input_bed_path, 
            sep='\t', 
            header=False, 
            index=False,
        )
    else:
        mosdepth_input_bed_path = None

    # prepare mosdepth input bam
    mosdepth_input_bam_path = os.path.join(tmpdir, 'mosdepth_input.bam')
    if subset_bam is None:
        subset_bam = load_perbase
    if (region_is_given and subset_bam): 
        # subset input bam file
        logutils.log(f'Subsetting input bam file', add_locstring=False)
        bameditor.samtools_view(
            in_bam_path=bam_path, 
            out_bam_path=mosdepth_input_bam_path, 
            region_bed_path=mosdepth_input_bed_path, 
            index=True,
        )
    else:  
        # do not subset bam
        os.symlink(bam_path, mosdepth_input_bam_path)

        # make index of the symlink
        bam_idx_path = bameditor.get_index_path(bam_path)
        if bam_idx_path is None:
            _ = pysam.index(mosdepth_input_bam_path)
        else:
            os.symlink(bam_idx_path, mosdepth_input_bam_path + '.bai')

    # make mosdepth outdir
    outdir = os.path.join(tmpdir, 'outdir')
    os.mkdir(outdir)

    # run mosdepth
    logutils.log(f'Running mosdepth', add_locstring=False)
    mosdepth_args = get_mosdepth_args(
        prefix=os.path.join(outdir, prefix), 
        bam_path=mosdepth_input_bam_path, 
        t=t, 
        use_median=use_median, 
        no_perbase=(not load_perbase),
        bed_path=mosdepth_input_bed_path, 
        window_size=window_size,
    )
    p = subprocess.check_call(mosdepth_args)

    # load mosdepth outputs
    region_outfile_path = os.path.join(outdir, f'{prefix}.regions.bed.gz')
    region_result_gdf = DepthRawDF.load_mosdepth(region_outfile_path, refver, use_median=use_median)

    # load mosdepth outputs - perbase
    perbase_outfile_path = os.path.join(outdir, f'{prefix}.per-base.bed.gz')
    if load_perbase:
        perbase_result_gdf = DepthRawDF.load_mosdepth_perbase(perbase_outfile_path, refver)
        if region_is_given:
            perbase_result_gdf = perbase_result_gdf.intersect(region_gdf)
    else:
        perbase_result_gdf = None

    # remove tmpdir
    shutil.rmtree(tmpdir)

    # return
    if load_perbase:
        return region_result_gdf, perbase_result_gdf
    else:
        return region_result_gdf



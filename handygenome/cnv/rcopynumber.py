import os
import tempfile
import subprocess
import shutil

import pyranges as pr
import numpy as np
import pandas as pd

import handygenome.common as common
import handygenome.ucscdata as ucscdata
import handygenome.pyranges_helper as pyranges_helper
import handygenome.cnv.misc as cnvmisc


RSCRIPT_COPYNUMBER = '/home/users/pjh/tools/miniconda/221104/miniconda3/envs/r-copynumber/lib/R/bin/Rscript'
COPYNUMBER_SCRIPT_PATH = os.path.join(common.R_DIR, 'copynumber.R')


def add_arm_info(gr, cytoband_gr, as_gr=True):
    if isinstance(gr, pd.DataFrame):
        gr = pr.PyRanges(gr, int64=False)

    new_cytoband_gr = pr.from_dict(
        {
            'Chromosome': cytoband_gr.Chromosome,
            'Start': cytoband_gr.Start,
            'End': cytoband_gr.End,
            'arm': [x[0] for x in cytoband_gr.Name],
        },
        int64=False,
    )

    joined_gr = pyranges_helper.join(gr, new_cytoband_gr, how='left', merge='first')
    if as_gr:
        return joined_gr
    else:
        return joined_gr.df


#def vaf_df_to_baf_df(vaf_df):
#    # make BAF series (all must be <= 0.5)
#    bafs = vaf_df['vaf'].copy()
#    bafs.where((bafs <= 0.5), (1 - bafs), inplace=True)
#    return pd.DataFrame.from_dict(
#        {
#            'Chromosome': vaf_df['CHROM'],
#            'Start': (vaf_df['POS'] - 1),
#            'End': vaf_df['POS'],
#            'baf': bafs,
#        }
#    )


def run_rcopynumber(
    depth_df, 
    baf_df=None, 
    refver=None, 
    cytoband_gr=None, 
    as_gr=True, 
    winsorize=True, 
    verbose=False, 
    **kwargs,
):
    def merge_depth_baf_respect_bafs(depth_df, baf_df):
        depth_gr = pr.PyRanges(depth_df, int64=False)
        baf_gr = pr.PyRanges(baf_df, int64=False)

        result = pyranges_helper.isec_union(depth_gr, baf_gr)
        result = pyranges_helper.join(result, depth_gr, how='left', merge='first')
            # Now "depth" column is added
        result = pyranges_helper.join(result, baf_gr, how='left', merge='first')
            # Now "baf" column is added. Missing values are np.nan

        return result

    def merge_depth_baf_average_bafs(depth_df, baf_df):
        depth_gr = pr.PyRanges(depth_df, int64=False)
        baf_gr = pr.PyRanges(baf_df, int64=False)
        result = pyranges_helper.join(depth_gr, baf_gr, how='left', merge='mean')

        return result

    def write_input_file(input_gr, cytoband_gr, infile_path):
        input_gr = add_arm_info(input_gr, cytoband_gr, as_gr=True)

        df_written_source = {
            'chrom': input_gr.Chromosome.array,
            'pos': (((input_gr.Start + 1) + input_gr.End) / 2).astype('int').array,
                # to make into 1-based coordinate
            'arm': input_gr.arm.array,
            'depth': input_gr.depth.array,
        }
        if 'baf' in input_gr.columns:
            df_written_source['baf'] = input_gr.baf.array

        df_written = pd.DataFrame.from_dict(df_written_source)

        # remove lines with missing 'depth' or 'arm'
        selector = np.logical_and(
            ~(df_written['depth'].isna().array), 
            ~(df_written['arm'].isna().array),
        )
        df_written = df_written.loc[selector, :]
        df_written.to_csv(infile_path, sep='\t', index=False)
        
        return df_written

    # sanity checks
    if (refver is None) and (cytoband_gr is None):
        raise Exception(f'At least one of "refver" or "cytoband_gr" argument must be set.')
    required_columns = {
        'depth': {'Chromosome', 'Start', 'End', 'depth'},
        'baf': {'Chromosome', 'Start', 'End', 'baf'},
    }
    if not required_columns['depth'].issubset(depth_df.columns):
        raise Exception(f'"depth_df" columns must include all of these: {required_columns["depth"]}')
    if baf_df is not None:
        if not required_columns['baf'].issubset(baf_df.columns):
            raise Exception(f'"baf_df" columns must include all of these: {required_columns["baf"]}')

    # make tmpdir and cytoband_gr
    tmpdir = tempfile.mkdtemp(dir=os.getcwd(), prefix=f'tmpdir_Rcopynumber_')
    if cytoband_gr is None:
        cytoband_gr = ucscdata.get_cytoband_gr(refver)

    # remove unnecessary columns
    depth_df = depth_df.loc[:, ['Chromosome', 'Start', 'End', 'depth']]
    if baf_df is not None:
        baf_df = baf_df.loc[:, ['Chromosome', 'Start', 'End', 'baf']]

    # make input data for r-copynumber
    if baf_df is None:
        input_gr = pr.PyRanges(depth_df, int64=False)
    else:
        input_gr = merge_depth_baf_average_bafs(depth_df, baf_df)

    infile_path = os.path.join(tmpdir, 'input.tsv.gz')
    df_written = write_input_file(input_gr, cytoband_gr, infile_path)

    # run R copynumber
    args = [RSCRIPT_COPYNUMBER, COPYNUMBER_SCRIPT_PATH, '--infile', infile_path]
    if winsorize:
        args.append('--winsorize')
    if verbose:
        args.append('--verbose')

    for key, val in kwargs.items():
        args.append(f'--{key}')
        args.append(str(val))

    p = subprocess.run(args, capture_output=False)  
        # "capture_output=False" lets stdout and stderr to be printed immediately
    p.check_returncode()
        
    # load copynumber output
    outfile_path = infile_path + '.out.gz'
    seg = pd.read_table(
        outfile_path, sep='\t', header=0, 
        dtype={'Chromosome': str, 'Start': int, 'End': int},
    )

    # remove tmpdir
    shutil.rmtree(tmpdir)

    # return
    if as_gr:
        seg = pr.PyRanges(seg, int64=False)
    return seg


def run_rcopynumber_compact(
    depth_df, 
    baf_df=None, 
    refver=None, 
    cytoband_gr=None, 
    as_gr=True, 
    winsorize=True, 
    verbose=False, 
    **kwargs,
):
    """Args:
        depth_df: Actually depth ratio df
            columns: Chromosome, Start, End, depth, (baf)
    """
    def merge_depth_baf_respect_bafs(depth_df, baf_df):
        depth_gr = pr.PyRanges(depth_df, int64=False)
        baf_gr = pr.PyRanges(baf_df, int64=False)

        result = pyranges_helper.isec_union(depth_gr, baf_gr)
        result = pyranges_helper.join(result, depth_gr, how='left', merge='first')
            # Now "depth" column is added
        result = pyranges_helper.join(result, baf_gr, how='left', merge='first')
            # Now "baf" column is added. Missing values are np.nan

        return result

    def merge_depth_baf_average_bafs(depth_df, baf_df):
        depth_gr = pr.PyRanges(depth_df, int64=False)
        baf_gr = pr.PyRanges(baf_df, int64=False)
        result = pyranges_helper.join(depth_gr, baf_gr, how='left', merge='mean')

        return result

    def make_merged_gr(depth_df, baf_df, cytoband_gr):
        if baf_df is None:
            merged_gr = pr.PyRanges(depth_df, int64=False)
        else:
            merged_gr = merge_depth_baf_average_bafs(depth_df, baf_df)
        merged_gr = add_arm_info(merged_gr, cytoband_gr, as_gr=True)
        return merged_gr

    def make_compact(merged_gr):
        merged_gr = merged_gr.sort()

        poslist = list()  # starts with 1
        coord_converter = dict()
        for idx, row in merged_gr.df.iterrows():
            key = (row['Start'], row['End'])
            pos = idx + 1
            poslist.append(pos)
            coord_converter[key] = pos

        rev_coord_converter = {val: key for key, val in coord_converter.items()}

        compacted_df = pd.DataFrame.from_dict({
            'chrom': merged_gr.Chromosome.array,
            'pos': poslist,
            'arm': merged_gr.arm.array,
            'depth': merged_gr.depth.array,
        })
        if 'baf' in merged_gr.columns:
            compacted_df['baf'] = merged_gr.baf.array

        return compacted_df, coord_converter, rev_coord_converter

    def write_input_file(compacted_df, infile_path):
        # remove lines with missing 'depth' or 'arm'
        selector = np.logical_and(
            ~(compacted_df['depth'].isna().array), 
            ~(compacted_df['arm'].isna().array),
        )
        df_written = compacted_df.loc[selector, :]
        df_written.to_csv(infile_path, sep='\t', index=False)
        
        return df_written

    def load_rcopynumber_output(outfile_path, rev_coord_converter):
        raw_seg = pd.read_table(outfile_path, sep='\t', header=0, dtype={'chrom': str})

        #chromlist = raw_seg['chrom'].array
        if 'logR.mean' in raw_seg.columns:
            depthlist = raw_seg['logR.mean'].array
            baflist = raw_seg['BAF.mean'].array
        else:
            depthlist = raw_seg['mean'].array
            baflist = None

        start0list = list()
        end0list = list()
        for idx, row in raw_seg.iterrows():
            start_coords = rev_coord_converter[row['start.pos']]
            end_coords = rev_coord_converter[row['end.pos']]
            start0list.append(start_coords[0])
            end0list.append(end_coords[1])

        source_dict = {
            'Chromosome': raw_seg['chrom'].array, 
            'Start': start0list,
            'End': end0list,
            'depth_mean': depthlist,
        }
        if baflist is not None:
            source_dict['baf_mean'] = baflist

        seg = pd.DataFrame.from_dict(source_dict)
        return seg

    # sanity checks
    if (refver is None) and (cytoband_gr is None):
        raise Exception(f'At least one of "refver" or "cytoband_gr" argument must be set.')
    required_columns = {
        'depth': {'Chromosome', 'Start', 'End', 'depth'},
        'baf': {'Chromosome', 'Start', 'End', 'baf'},
    }
    if not required_columns['depth'].issubset(depth_df.columns):
        raise Exception(f'"depth_df" columns must include all of these: {required_columns["depth"]}')
    if baf_df is not None:
        if not required_columns['baf'].issubset(baf_df.columns):
            raise Exception(f'"baf_df" columns must include all of these: {required_columns["baf"]}')

    # make tmpdir and cytoband_gr
    tmpdir = tempfile.mkdtemp(dir=os.getcwd(), prefix=f'tmpdir_Rcopynumber_')
    if cytoband_gr is None:
        cytoband_gr = ucscdata.get_cytoband_gr(refver)

    # remove unnecessary columns
    depth_df = depth_df.loc[:, ['Chromosome', 'Start', 'End', 'depth']]
    if baf_df is not None:
        baf_df = baf_df.loc[:, ['Chromosome', 'Start', 'End', 'baf']]

    # make input data for r-copynumber
    merged_gr = make_merged_gr(depth_df, baf_df, cytoband_gr)
    compacted_df, coord_converter, rev_coord_converter = make_compact(merged_gr)
    infile_path = os.path.join(tmpdir, 'input.tsv.gz')
    df_written = write_input_file(compacted_df, infile_path)

    # run R copynumber
    args = [RSCRIPT_COPYNUMBER, COPYNUMBER_SCRIPT_PATH, '--infile', infile_path, '--asis']
    if winsorize:
        args.append('--winsorize')
    if verbose:
        args.append('--verbose')

    for key, val in kwargs.items():
        args.append(f'--{key}')
        args.append(str(val))

    p = subprocess.run(args, capture_output=False)  
        # "capture_output=False" lets stdout and stderr to be printed immediately
    p.check_returncode()
        
    # load copynumber output
    outfile_path = infile_path + '.out.gz'
    seg = load_rcopynumber_output(outfile_path, rev_coord_converter)

    # remove tmpdir
    shutil.rmtree(tmpdir)

    # return
    if as_gr:
        seg = pr.PyRanges(seg, int64=False)
    return seg


def run_rcopynumber_unified(
    depth_df, 
    baf_df=None, 
    refver=None, 
    cytoband_gr=None, 
    as_gr=True, 
    winsorize=True, 
    compact=False,
    verbose=False, 
    **kwargs,
):
    """Args:
        depth_df: Actually depth ratio df
            columns: Chromosome, Start, End, depth, (baf)
    """
    def merge_depth_baf_respect_bafs(depth_df, baf_df):
        depth_gr = pr.PyRanges(depth_df, int64=False)
        baf_gr = pr.PyRanges(baf_df, int64=False)

        result = pyranges_helper.isec_union(depth_gr, baf_gr)
        result = pyranges_helper.join(result, depth_gr, how='left', merge='first')
            # Now "depth" column is added
        result = pyranges_helper.join(result, baf_gr, how='left', merge='first')
            # Now "baf" column is added. Missing values are np.nan

        return result

    def merge_depth_baf_average_bafs(depth_df, baf_df):
        depth_gr = pr.PyRanges(depth_df, int64=False)
        baf_gr = pr.PyRanges(baf_df, int64=False)
        result = pyranges_helper.join(depth_gr, baf_gr, how='left', merge='mean')

        return result

    def make_merged_gr(depth_df, baf_df, cytoband_gr):
        if baf_df is None:
            merged_gr = pr.PyRanges(depth_df, int64=False)
        else:
            merged_gr = merge_depth_baf_average_bafs(depth_df, baf_df)

        merged_gr = add_arm_info(merged_gr, cytoband_gr, as_gr=True)

        selector = np.logical_and(
            merged_gr.depth.notna().to_numpy(), 
            merged_gr.arm.notna().to_numpy(),
        )
        merged_gr = merged_gr[pd.Series(selector)] 

        return merged_gr

    def make_compact(merged_gr):
        merged_gr = merged_gr.sort()

        poslist = list(range(1, merged_gr.df.shape[0] + 1))  # starts with 1
        key_iter = zip(merged_gr.df['Start'], merged_gr.df['End'])
        coord_converter = dict(zip(key_iter, iter(poslist)))
        rev_coord_converter = {val: key for key, val in coord_converter.items()}

        compacted_df = pd.DataFrame.from_dict({
            'chrom': merged_gr.Chromosome.array,
            'pos': poslist,
            'arm': merged_gr.arm.array,
            'depth': merged_gr.depth.array,
        })
        if 'baf' in merged_gr.columns:
            compacted_df['baf'] = merged_gr.baf.array

        return compacted_df, coord_converter, rev_coord_converter

    def make_noncompact_input(merged_gr):
        df_written = pd.DataFrame.from_dict({
            'chrom': merged_gr.Chromosome.array,
            'pos': (((merged_gr.Start + 1) + merged_gr.End) / 2).astype('int').array,
                # to make into 1-based coordinate
            'arm': merged_gr.arm.array,
            'depth': merged_gr.depth.array,
        })
        if 'baf' in merged_gr.columns:
            df_written['baf'] = merged_gr.baf.array

        return df_written

    def load_rcopynumber_output(outfile_path, rev_coord_converter):
        raw_seg = pd.read_table(outfile_path, sep='\t', header=0, dtype={'chrom': str})

        #chromlist = raw_seg['chrom'].array
        if 'logR.mean' in raw_seg.columns:
            depthlist = raw_seg['logR.mean'].array
            baflist = raw_seg['BAF.mean'].array
        else:
            depthlist = raw_seg['mean'].array
            baflist = None

        tmp = raw_seg.apply(
            lambda row: (
                rev_coord_converter[row['start.pos']][0],
                rev_coord_converter[row['end.pos']][1],
            ),
            axis=1,
        )
        zipped = zip(*tmp)
        start0list = next(zipped)
        end0list = next(zipped)

        source_dict = {
            'Chromosome': raw_seg['chrom'].array, 
            'Start': start0list,
            'End': end0list,
            'depth_mean': depthlist,
        }
        if baflist is not None:
            source_dict['baf_mean'] = baflist

        seg = pd.DataFrame.from_dict(source_dict)
        return seg

    # sanity checks
    if (refver is None) and (cytoband_gr is None):
        raise Exception(f'At least one of "refver" or "cytoband_gr" argument must be set.')
    required_columns = {
        'depth': {'Chromosome', 'Start', 'End', 'depth'},
        'baf': {'Chromosome', 'Start', 'End', 'baf'},
    }
    if not required_columns['depth'].issubset(depth_df.columns):
        raise Exception(f'"depth_df" columns must include all of these: {required_columns["depth"]}')
    if baf_df is not None:
        if not required_columns['baf'].issubset(baf_df.columns):
            raise Exception(f'"baf_df" columns must include all of these: {required_columns["baf"]}')

    # make tmpdir and cytoband_gr
    tmpdir = tempfile.mkdtemp(dir=os.getcwd(), prefix=f'tmpdir_Rcopynumber_')
    if cytoband_gr is None:
        cytoband_gr = ucscdata.get_cytoband_gr(refver)

    # remove unnecessary columns
    depth_df = depth_df.loc[:, ['Chromosome', 'Start', 'End', 'depth']]
    if baf_df is not None:
        baf_df = baf_df.loc[:, ['Chromosome', 'Start', 'End', 'baf']]

    # make input data for r-copynumber
    merged_gr = make_merged_gr(depth_df, baf_df, cytoband_gr)
    if compact:
        df_written, coord_converter, rev_coord_converter = make_compact(merged_gr)
    else:
        df_written = make_noncompact_input(merged_gr)

    infile_path = os.path.join(tmpdir, 'input.tsv.gz')
    df_written.to_csv(infile_path, sep='\t', index=False)

    # run R copynumber
    args = [RSCRIPT_COPYNUMBER, COPYNUMBER_SCRIPT_PATH, '--infile', infile_path, '--asis']
    if winsorize:
        args.append('--winsorize')
    if verbose:
        args.append('--verbose')

    for key, val in kwargs.items():
        args.append(f'--{key}')
        args.append(str(val))

    p = subprocess.run(args, capture_output=False)  
        # "capture_output=False" lets stdout and stderr to be printed immediately
    p.check_returncode()
        
    # load copynumber output
    outfile_path = infile_path + '.out.gz'
    if compact:
        seg = load_rcopynumber_output(outfile_path, rev_coord_converter)
    else:
        seg = pd.read_table(
            outfile_path, sep='\t', header=0, 
            dtype={'Chromosome': str, 'Start': int, 'End': int},
        )

    # remove tmpdir
    shutil.rmtree(tmpdir)

    # return
    if as_gr:
        seg = pr.PyRanges(seg, int64=False)
    return seg


def add_CNn_to_targetseq_segment_gr(segment_gr, targetregion_gr, refver, is_female):
    """With targeted sequencing analysis, start and end coordinates of 
    "run_rcopynumber*" function output does not match actual target regions.
    To add CNn values to such a segment PyRanges, 1) CNn values are added
    to the target region PyRanges, 2) segment PyRanges is joined to the 
    annotated target region PyRanges to acquire CNn values
    Doing like this can handle an exceptional case where non-existent segment 
    subregion (a region included in the segment df but not part of the target 
    regions) spans PAR.
    """
    targetregion_gr = pyranges_helper.join(
        targetregion_gr, cnvmisc.get_CNn_gr(refver, is_female),
        how='left', merge='first', as_gr=True,
    )
    segment_gr = pyranges_helper.join(
        segment_gr, targetregion_gr, how='left', merge='longest', as_gr=True,
    )
    return segment_gr


def add_CNn_to_wgs_segment_gr(segment_gr, refver, is_female):
    return pyranges_helper.join(
        segment_gr, cnvmisc.get_CNn_gr(refver, is_female), 
        how='left', merge='longest', as_gr=True,
    )



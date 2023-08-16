import os
import tempfile
import subprocess
import shutil

import pyranges as pr
import numpy as np
import pandas as pd

import handygenome
import handygenome.workflow as workflow
import handygenome.ucscdata as ucscdata
import handygenome.pyranges_helper as pyranges_helper
import handygenome.cnv.misc as cnvmisc


LOGGER_INFO = workflow.get_debugging_logger(verbose=False)
LOGGER_DEBUG = workflow.get_debugging_logger(verbose=True)

RSCRIPT_COPYNUMBER = '/home/users/pjh/tools/miniconda/221104/miniconda3/envs/r-copynumber/lib/R/bin/Rscript'
COPYNUMBER_SCRIPT_PATH = os.path.join(handygenome.DIRS['R'], 'copynumber.R')


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


def add_arm_info_new(df, cytoband_df, refver, as_gr=False):
    new_cytoband_df = cytoband_df.loc[:, ['Chromosome', 'Start', 'End']]
    new_cytoband_df['arm'] = [x[0] for x in cytoband_df.Name]
    joined_df = pyranges_helper.join(
        df, new_cytoband_df, how='left', merge='first',
        sort=True, refver=refver,
    )

    if as_gr:
        return pr.PyRanges(joined_df)
    else:
        return joined_df


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
    result = pyranges_helper.join_new(depth_gr, baf_gr, how='left', merge='mean')

    return result


def make_merged_gr(depth_df, baf_df, cytoband_gr=None, refver=None):
    if cytoband_gr is None:
        cytoband_gr = ucscdata.get_cytoband_gr(refver)

    if baf_df is None:
        merged_gr = pr.PyRanges(depth_df, int64=False)
    else:
        merged_gr = merge_depth_baf_average_bafs(depth_df, baf_df)
        if getattr(merged_gr, 'baf_raw').isna().all():
            merged_gr = merged_gr[merged_gr.columns.drop(['baf_raw']).to_list()]

    merged_gr = add_arm_info(merged_gr, cytoband_gr, as_gr=True)

    selector = np.logical_and(
        getattr(merged_gr, 'depth_raw').notna().to_numpy(), 
        getattr(merged_gr, 'arm').notna().to_numpy(),
    )
    merged_gr_wona = merged_gr[pd.Series(selector)].sort()

    return merged_gr, merged_gr_wona


def make_merged_df(depth_df, baf_df, cytoband_df, refver, logger, nproc):
    if cytoband_df is None:
        cytoband_df = ucscdata.get_cytoband_gr(refver=refver, as_gr=False)

    if baf_df is None:
        merged_df = depth_df
    else:
        logger.debug(f'Joining depth and baf dataframes')
        merged_df = pyranges_helper.join(
            depth_df, baf_df, 
            how='left', 
            merge='mean', 
            find_nearest=False,
            as_gr=False,
            sort=True,
            refver=refver,
            verbose=True,
            nproc=nproc,
        )
        if merged_df['baf_raw'].isna().all():
            merged_df = merged_df.drop('baf_raw')

    logger.debug(f'Adding chromosome arm information')
    merged_df = add_arm_info_new(merged_df, cytoband_df, refver, as_gr=False)

    logger.debug(f'Removing rows with missing depth or arm information')
    selector = np.logical_and(
        merged_df['depth_raw'].notna().to_numpy(), 
        merged_df['arm'].notna().to_numpy(), 
    )
    merged_df_wona = merged_df.loc[selector, :]

    return merged_df, merged_df_wona


#def make_compact(merged_gr):
#    poslist = list(range(1, merged_gr.df.shape[0] + 1))  # starts with 1
#    key_iter = zip(merged_gr.Start, merged_gr.End)
#    coord_converter = dict(zip(key_iter, iter(poslist)))
#    rev_coord_converter = {val: key for key, val in coord_converter.items()}
#
#    compacted_df = pd.DataFrame.from_dict({
#        'chrom': merged_gr.Chromosome.array,
#        'pos': poslist,
#        'arm': merged_gr.arm.array,
#        'depth': getattr(merged_gr, 'depth_raw').array,
#    })
#    if 'baf_raw' in merged_gr.columns:
#        compacted_df['baf'] = getattr(merged_gr, 'baf_raw').array
#
#    return compacted_df, coord_converter, rev_coord_converter


def make_compact_df(merged_df, refver):
    assert isinstance(merged_df, pd.DataFrame)

    # sort input df
    merged_df = cnvmisc.sort_genome_df(merged_df, refver)

    # make mapping
    compact_pos1_list = list()

    original_to_compact = dict()
    for chrom, subdf in merged_df.groupby(merged_df['Chromosome'].to_numpy()):
        compact_pos1_sublist = np.arange(subdf.shape[0]) + 1
        key_iter = zip(subdf['Start'], subdf['End'])
        original_to_compact[chrom] = dict(zip(key_iter, compact_pos1_sublist))
        compact_pos1_list.extend(compact_pos1_sublist)

    compact_to_original = dict()
    for chrom, subdic in original_to_compact.items():
        compact_to_original[chrom] = {
            compact_pos1: original_tup for 
            original_tup, compact_pos1 in subdic.items()
        }

    # make compact df
    compact_df = pd.DataFrame.from_dict(
        {
            'chrom': merged_df.Chromosome.array,
            'pos': compact_pos1_list,
            'arm': merged_df.arm.array,
            'depth': merged_df.depth_raw.array,
        }
    )
    if 'baf_raw' in merged_df.columns:
        compact_df['baf'] = merged_df['baf_raw'].array

    return compact_df, compact_to_original


def make_noncompact_input(merged_gr):
    df_written = pd.DataFrame.from_dict({
        'chrom': merged_gr.Chromosome.array,
        'pos': (((merged_gr.Start + 1) + merged_gr.End) / 2).astype('int').array,
            # to make into 1-based coordinate
        'arm': merged_gr.arm.array,
        'depth': getattr(merged_gr, 'depth_raw').array,
    })
    if 'baf_raw' in merged_gr.columns:
        df_written['baf'] = getattr(merged_gr, 'baf_raw').array

    return df_written


def make_noncompact_input_df(merged_df):
    df_written = pd.DataFrame.from_dict({
        'chrom': merged_df.Chromosome.array,
        'pos': (((merged_df.Start + 1) + merged_df.End) / 2).astype('int').array,
            # to make into 1-based coordinate
        'arm': merged_df.arm.array,
        'depth': merged_df.depth_raw.array,
    })
    if 'baf_raw' in merged_df.columns:
        df_written['baf'] = merged_df.baf_raw.array

    return df_written


def load_compact_rcopynumber_output(outfile_path, compact_to_original):
    # load text file
    raw_seg = pd.read_table(
        outfile_path, 
        sep='\t', 
        header=0, 
        dtype={'chrom': str},
    )

    # trim segment limits
    for chrom in set(raw_seg['chrom']):
        compact_pos1_min = min(compact_to_original[chrom].keys())
        compact_pos1_max = max(compact_to_original[chrom].keys())
        selector = (raw_seg['chrom'] == chrom)
        raw_seg.loc[selector, 'start.pos'] = np.clip(
            raw_seg.loc[selector, 'start.pos'], 
            compact_pos1_min, 
            compact_pos1_max,
        )
        raw_seg.loc[selector, 'end.pos'] = np.clip(
            raw_seg.loc[selector, 'end.pos'], 
            compact_pos1_min, 
            compact_pos1_max,
        )

    raw_seg.drop_duplicates(
        ['chrom', 'start.pos', 'end.pos'], 
        keep='first', 
        inplace=True,
    )

    # set parameters
    if 'logR.mean' in raw_seg.columns:
        depthlist = raw_seg['logR.mean'].array
        baflist = raw_seg['BAF.mean'].array
    else:
        depthlist = raw_seg['mean'].array
        baflist = None
                
    # do conversion
    def start_converter(row):
        return compact_to_original[row['chrom']][row['start.pos']][0]

    def end_converter(row):
        return compact_to_original[row['chrom']][row['end.pos']][1]

    start0_list = raw_seg.loc[:, ['chrom', 'start.pos']].apply(start_converter, axis=1)
    end0_list = raw_seg.loc[:, ['chrom', 'end.pos']].apply(end_converter, axis=1)

    assert start0_list.notna().all()
    assert end0_list.notna().all()

    # make result
    source_dict = {
        'Chromosome': raw_seg['chrom'].array, 
        'Start': start0_list.array,
        'End': end0_list.array,
        'depth_segment_mean': depthlist,
    }
    if baflist is not None:
        source_dict['baf_segment_mean'] = baflist

    seg = pd.DataFrame.from_dict(source_dict)

    return seg


def sanitycheck(
    argtype1, 
    argtype2, 
    argtype3,
    required_columns,
    depth_df,
    baf_df,
    refver, 
    cytoband_gr, 
    merged_df, 
    merged_df_wona,
):
    assert refver is not None, f'"refver" must be given'

    if depth_df is not None:
        if not set(required_columns['depth']).issubset(depth_df.columns):
            raise Exception(
                f'"depth_df" columns must include all of these: {required_columns["depth"]}'
            )
    if baf_df is not None:
        if not set(required_columns['baf']).issubset(baf_df.columns):
            raise Exception(
                f'"baf_df" columns must include all of these: {required_columns["baf"]}'
            )

    if argtype1 or argtype3:
        if (refver is None) and (cytoband_gr is None):
            raise Exception(f'At least one of "refver" or "cytoband_gr" argument must be set.')
    elif argtype2:
        if (merged_df is None) or (merged_df_wona is None):
            raise Exception(f'Both "merged_df", and "merged_df_wona" must be given.')


def input_df_handling(input_df, remove_unassembled_contigs, req_cols):
    if input_df is None:
        return None
    else:
        input_df = input_df.copy()
        input_df = cnvmisc.arg_into_df(input_df)
        if remove_unassembled_contigs:
            input_df = cnvmisc.remove_unassembled_contigs(input_df)

        # remove unnecessary columns
        input_df = input_df.loc[:, req_cols]

        return input_df


def argtype1_prepare_input(
    depth_df, baf_df, remove_unassembled_contigs, required_columns,
    cytoband_gr, refver, logger, join_nproc,
):
    depth_df = input_df_handling(
        depth_df, remove_unassembled_contigs, required_columns['depth'],
    )
    baf_df = input_df_handling(
        baf_df, remove_unassembled_contigs, required_columns['baf'],
    )

    # make input data for r-copynumber
    merged_df, merged_df_wona = make_merged_df(
        depth_df, baf_df, cytoband_gr, refver, logger, join_nproc,
    )

    return merged_df, merged_df_wona


def argtype3_prepare_input(
    baf_df, remove_unassembled_contigs, required_columns,
    cytoband_gr, refver, logger, join_nproc,
):
    merged_df = input_df_handling(
        baf_df, remove_unassembled_contigs, required_columns['baf'],
    )
    merged_df['depth_raw'] = 1

    if cytoband_gr is None:
        cytoband_gr = ucscdata.get_cytoband_gr(refver=refver, as_gr=False)
    merged_df = add_arm_info_new(merged_df, cytoband_gr, refver, as_gr=False)

    selector = merged_df['arm'].notna().to_numpy()
    merged_df_wona = merged_df.loc[selector, :]

    return merged_df, merged_df_wona


###########################################################################


def run_rcopynumber(
    *,
    depth_df=None, 
    baf_df=None, 
    refver=None, 
    cytoband_gr=None, 

    merged_df=None, 
    merged_df_wona=None,

    as_gr=True, 
    winsorize=False, 
    compact=False,
    verbose=False, 
    remove_unassembled_contigs=True,
    join_nproc=1,
    **kwargs,

    #df_written=None,
):
    """Args:
        depth_df:
            columns: Chromosome, Start, End, depth_raw
        baf_df:
            columns: Chromosome, Start, End, baf_raw
        refver: A required argument

    Returns:
        seg: Columns ("Chromosome", "Start", "End", "depth_segment_mean", +/- "baf_segment_mean")
        merged_df: Columns ("Chromosome", "Start", "End", "depth_raw", +/- "baf_raw")
    """
    # argtype
    argtype1 = (depth_df is not None)
    argtype2 = (merged_df is not None) or (merged_df_wona is not None)
    argtype3 = (depth_df is None) and (baf_df is not None)
    if (argtype1 + argtype2 + argtype3) != 1:
        raise Exception(f'Invalid argument usage')

    # sanity checks
    required_columns = {
        'depth': ['Chromosome', 'Start', 'End', 'depth_raw'],
        'baf': ['Chromosome', 'Start', 'End', 'baf_raw'],
    }
    sanitycheck(
        argtype1, 
        argtype2, 
        argtype3,
        required_columns,
        depth_df,
        baf_df,
        refver, 
        cytoband_gr, 
        merged_df, 
        merged_df_wona,
    )

    # set logger
    logger = (LOGGER_DEBUG if verbose else LOGGER_INFO)

    # make tmpdir
    tmpdir = tempfile.mkdtemp(dir=os.getcwd(), prefix=f'tmpdir_Rcopynumber_')

    if argtype1:
        merged_df, merged_df_wona = argtype1_prepare_input(
            depth_df, baf_df, remove_unassembled_contigs, required_columns,
            cytoband_gr, refver, logger, join_nproc,
        )
    elif argtype3:
        merged_df, merged_df_wona = argtype3_prepare_input(
            baf_df, remove_unassembled_contigs, required_columns,
            cytoband_gr, refver, logger, join_nproc,
        )
    elif argtype2:
        pass

    # write input file for R copynumber
    if compact:
        logger.debug(f'Making compact')
        df_written, compact_to_original = make_compact_df(merged_df_wona, refver)
    else:
        df_written = make_noncompact_input_df(merged_df_wona)

    logger.debug(f'Writing dataframe into file for R-copynumber input')
    infile_path = os.path.join(tmpdir, 'input.tsv.gz')
    df_written.to_csv(infile_path, sep='\t', index=False)

    # run R copynumber
    logger.debug(f'Running R-copynumber')
    args = [RSCRIPT_COPYNUMBER, COPYNUMBER_SCRIPT_PATH, '--infile', infile_path]
    if compact:
        args.append('--asis')

    if winsorize:
        args.append('--winsorize')
    if verbose:
        args.append('--verbose')

    for key, val in kwargs.items():
        if val is not None:
            args.append(f'--{key}')
            args.append(str(val))

    p = subprocess.run(args, capture_output=False)  
        # "capture_output=False" lets stdout and stderr to be printed immediately
    p.check_returncode()
        
    # load copynumber output
    logger.debug(f'Loading R-copynumber output')
    outfile_path = infile_path + '.out.gz'
    if compact:
        seg = load_compact_rcopynumber_output(outfile_path, compact_to_original)
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
        merged_df = pr.PyRanges(merged_df, int64=False)
    #else:
        #merged_gr = merged_gr.df

    return seg, merged_df

run_rcopynumber_unified = run_rcopynumber


def add_CNn_to_targetseq_segment_gr(segment_gr, targetregion_gr, refver, is_female, as_gr=True):
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
        how='left', merge='first', as_gr=False,
    )
    segment_gr = pyranges_helper.join(
        segment_gr, targetregion_gr, how='left', merge='longest', as_gr=as_gr,
        sort=True, refver=refver,
    )
    return segment_gr

add_CNn_to_targetseq_segment = add_CNn_to_targetseq_segment_gr


def add_CNn_to_wgs_segment_gr(segment_gr, refver, is_female, as_gr=True):
    return pyranges_helper.join(
        segment_gr, cnvmisc.get_CNn_gr(refver, is_female), 
        how='left', merge='longest', as_gr=as_gr,
        sort=True, refver=refver,
    )

add_CNn_to_wgs_segment = add_CNn_to_wgs_segment_gr



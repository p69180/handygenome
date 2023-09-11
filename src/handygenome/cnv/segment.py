import pandas as pd
import numpy as np
import pyranges as pr

from handygenome.genomedf import GenomeDataFrame as GDF









##########################


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

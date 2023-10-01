import pandas as pd
import numpy as np

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri


'''rpy2 package seems to utilize 'root' logger. 
In logutils module, when 'root' logger is used, ("log_old" function), rpy2 emits some log messages.
When logutils module use a logger with another name, (current "log" function), the rpy2 log message is not seen anymore.
'''


DNACOPY = importr('DNAcopy')


#def run_segmentation(
#    gdf, 
#    colname='value', 
#    out_colname=None, 
#    only_coords=True,
#    #join_raw=False,
#    #N_colname='N',
#
#    smoothing=False, 
#    verbose=True, 
#    smooth_kwargs=dict(), 
#    segment_kwargs=dict(),
#):
#    """kwargs for 'DNAcopy::segment' function (https://bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf)
#        - weights
#        - alpha
#        - nperm
#        - p.method
#        - min.width
#        - kmax
#        - nmin
#        - eta
#        - sbdry
#        - trim
#        - undo.splits
#        - undo.prune
#        - undo.SD
#        - verbose
#        
#    - Input position to R dnacopy is calculated as floor((start0 + end0 - 1) / 2)
#    """ 
#    # sanitycheck & argument handling
#    assert colname in gdf.columns
#    if out_colname is None:
#        out_colname = colname
#
#    # run R function
#    seg_df = run_dnacopy_segment(
#        chromosomes=gdf['Chromosome'], 
#        positions=gdf.get_midpoints(), 
#        values=gdf[colname],
#
#        smoothing=smoothing, 
#        verbose=verbose, 
#        smooth_kwargs=smooth_kwargs,
#        segment_kwargs=segment_kwargs,
#    )
#
#    # postprocess df
#    seg_df.drop('ID', axis=1, inplace=True)
#    seg_df.rename(
#        columns={
#            'chrom': 'Chromosome', 
#            'loc.start': 'Start',  # 0-based closed interval
#            'loc.end': 'End',  # 0-based closed interval
#            'num.mark': 'num_data', 
#            'seg.mean': out_colname,
#        }, 
#        inplace=True
#    )
#    seg_df['End'] = seg_df['End'] + 1  # to make into 0-based half-open system
#
#    if only_coords:
#        seg_df.drop(['num_data', out_colname], axis=1, inplace=True)
#
#    return seg_df
#
#    #result = GDF.from_frame(seg_df, gdf.refver)
#
##    if join_raw:
##        result = result.join(
##            gdf, 
##            other_cols=[colname], 
##            how='left', 
##            merge=['mean', 'std'], 
##            ddof=0,
##            overlapping_length=False,
##            N_colname=N_colname,
##        )
#
#    #return result


def rdataframe_to_df(rdataframe):
    with (ro.default_converter + pandas2ri.converter).context():
        return ro.conversion.get_conversion().rpy2py(rdataframe)


def run_dnacopy_segment(
    chromosomes, positions, values,
    *,
    N_colname='N',
    value_colname='value',
    smoothing=False, 
    verbose=True, 
    smooth_kwargs=dict(), 
    segment_kwargs=dict(),
):
    """kwargs for 'DNAcopy::segment' function (https://bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf)
        - weights
        - alpha
        - nperm
        - p.method
        - min.width
        - kmax
        - nmin
        - eta
        - sbdry
        - trim
        - undo.splits
        - undo.prune
        - undo.SD
        - verbose
        
    - Works well when 'positions' argument contains 0 or negative integers.
    """ 
    # prepare R function arguments
    default_segment_kwargs = {
        'alpha': 0.01,
        'nperm': 10000,
        'p_method': 'hybrid',
        'min_width': 2,
        'kmax': 25,
        'nmin': 200,
        'eta': 0.05,
        'trim': 0.025,
        'undo_splits': 'none',
        'undo_prune': 0.05,
        'undo_SD': 3,
    }
    assert set(segment_kwargs.keys()).issubset(default_segment_kwargs.keys())
    segment_kwargs = default_segment_kwargs | segment_kwargs
    segment_kwargs['verbose'] = int(verbose)

    default_smooth_kwargs = {
        'smooth_region': 10, 
        'outlier_SD_scale': 4, 
        'smooth_SD_scale': 2, 
        'trim': 0.025, 
    }
    assert set(smooth_kwargs.keys()).issubset(default_smooth_kwargs.keys())
    smooth_kwargs = default_smooth_kwargs | smooth_kwargs

    # main
    arg_chrom = ro.StrVector(chromosomes)
    arg_pos = ro.IntVector(positions)
    arg_value = ro.FloatVector(values)

    cnaobj = DNACOPY.CNA(arg_value, arg_chrom, arg_pos, data_type='logratio')
    if smoothing:
        cnaobj = DNACOPY.smooth_CNA(cnaobj, **smooth_kwargs)
    segresult = DNACOPY.segment(cnaobj, **segment_kwargs)

    # turn into pandas dataframe
    segdf = rdataframe_to_df(segresult.rx2['output'])
    # segdf is an R dataframe which has following columns: 
        # ID, chrom, loc.start, loc.end, num.mark, seg.mean

    segdf.drop('ID', axis=1, inplace=True)
    segdf.rename(
        columns={
            'chrom': 'Chromosome', 
            'loc.start': 'Start',  # 0-based closed interval
            'loc.end': 'End',  # 0-based closed interval
            'num.mark': N_colname, 
            'seg.mean': value_colname,
        }, 
        inplace=True
    )
    segdf['End'] = segdf['End'] + 1  # to make into 0-based half-open system

    return segdf

    #return {
    #    'Chromosome': np.array(segdf.rx2['chrom']),
    #    'Start': np.array(segdf.rx2['loc.start']),
    #    'End': (np.array(segdf.rx2['loc.end']) + 1),
    #    'N': np.array(segdf.rx2['num.mark']),
    #    'value': np.array(segdf.rx2['seg.mean']),
    #}



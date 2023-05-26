import os
import collections
import random
import itertools
import functools
import multiprocessing
import operator
import warnings
import inspect
import pickle

import Bio.SeqUtils
import numpy as np
import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
#from scipy.stats.mstats import winsorize

import handygenome.common as common
import handygenome.workflow as workflow
import handygenome.pyranges_helper as pyranges_helper
import handygenome.cnv.mosdepth as libmosdepth
import handygenome.cnv.gcfraction as gcfraction
import handygenome.assemblyspec as libassemblyspec
import handygenome.deco as deco


LOGGER_DEBUG = workflow.get_debugging_logger(verbose=True)
LOGGER_INFO = workflow.get_debugging_logger(verbose=False)

MALE_HAPLOID_CHROMS = ('X', 'Y', 'chrX', 'chrY')
DEFAULT_CNT_WEIGHT = 5


class CCFInfo(
    collections.namedtuple('CCFInfo', ('CNm', 'ccf', 'mutated_allele'))
):
    pass


class CPPair(
    collections.namedtuple('CPPair', ('cellularity', 'ploidy'))
):
    pass


####################
# argument handler #
####################

def check_having_coord_cols(arg):
    if not (
        {'Chromosome', 'Start', 'End'}.issubset(arg.columns)
        or arg.index.names == ['Chromosome', 'Start', 'End']
    ):
        raise Exception(f'Input dataframe must have columns "Chromosome", "Start", and "End".')


def arg_into_df(arg):
    if arg is None:
        return None
    elif isinstance(arg, pd.DataFrame):
        check_having_coord_cols(arg)
        return arg
    elif isinstance(arg, pr.PyRanges):
        return arg.df
    else:
        raise Exception(f'Argument must be either None, pd.DataFrame, or pr.PyRanges')


def arg_into_gr(arg):
    if arg is None:
        return None
    elif isinstance(arg, pd.DataFrame):
        check_having_coord_cols(arg)
        return pr.PyRanges(arg)
    elif isinstance(arg, pr.PyRanges):
        return arg
    else:
        raise Exception(f'Argument must be either None, pd.DataFrame, or pr.PyRanges')


def genome_df_sanitycheck(df):
    if not isinstance(df, pr.PyRanges):
        if not {'Chromosome', 'Start', 'End'}.issubset(df.columns):
            raise Exception(f'Genome DataFrame must contain these columns: "Chromosome", "Start", "End".')


def get_genome_df_annotcols(df):
    return df.columns.drop(['Chromosome', 'Start', 'End'])


def genome_df_groupkey(df, refver):
    chromdict = common.DEFAULT_CHROMDICTS[refver]
    chrom_indexes = df['Chromosome'].apply(chromdict.contigs.index)
    coord_arr = np.stack(
        [chrom_indexes.to_numpy(), df['Start'].to_numpy(), df['End'].to_numpy()], 
        axis=1,
    )
    return common.array_grouper(coord_arr, omit_values=True)[2]


def sort_genome_df(df, refver):
    chromdict = common.DEFAULT_CHROMDICTS[refver]
    chrom_indexes = df['Chromosome'].apply(lambda x: chromdict.contigs.index(x))
    sortkey = np.lexsort([df['End'], df['Start'], chrom_indexes])
    result = df.iloc[sortkey, :]
    result.reset_index(drop=True, inplace=True)
    return result


def remove_unassembled_contigs(df):
    selector = df['Chromosome'].apply(
        lambda x: common.RE_PATS['assembled_chromosome'].fullmatch(x) is not None
    )
    return df.loc[selector, :]


def group_df_bychrom(df, as_gr=True):
    gr = arg_into_gr(df).sort()
    chroms = set(gr.Chromosome)
    if as_gr:
        result = {chrom: gr[chrom] for chrom in chroms}
    else:
        result = {chrom: gr[chrom].df for chrom in chroms}
    return result


################
# peak finding #
################

def get_1d_peaks(row, invert=False):
    #assert not hasattr(row, '__next__'), f'"row" must not be an iterator'

    row = np.array(row)
    if invert:
        row = -1 * row

    peaks_result, _ = scipy.signal.find_peaks(row)
    if len(peaks_result) == 0:  # uniformly increasing
        peak_indexes = list()
        if len(row) > 1:
            if row[0] > row[1]:
                peak_indexes.append(0)
            if row[-1] > row[-2]:
                peak_indexes.append(len(row) - 1)
        elif len(row) == 1:
            peak_indexes.append(0)
            
#        peak_indexes = sorted(
#            (
#                x[0] for x in
#                common.multi_max(enumerate(row), key=operator.itemgetter(1))
#            )
#        )
    else:
        diff = np.diff(row)

        peak_indexes = list()
        for within_row_idx in peaks_result:
            peak_indexes.append(within_row_idx)
            # leftward
            current_idx = within_row_idx
            while True:
                current_idx -= 1
                if diff[current_idx] != 0:
                    break
                else:
                    peak_indexes.append(current_idx)

                if current_idx == 0:
                    break
            # rightward
            current_idx = within_row_idx - 1
            while True:
                current_idx += 1
                if diff[current_idx] != 0:
                    break
                else:
                    peak_indexes.append(current_idx + 1)

                if current_idx == len(diff) - 1:
                    break

        peak_indexes.sort()

    return np.array(peak_indexes)


def get_hist_peaks(data, weights=None, bins=10, threshold=None):
    hist, histbins = np.histogram(data, bins=bins, weights=weights, density=True)
    bin_midpoints = 0.5 * (histbins[:-1] + histbins[1:])

    peak_indexes = get_1d_peaks(hist, invert=False)
    if threshold is not None:
        peak_indexes = peak_indexes[hist[peak_indexes] > threshold]

    peak_values = bin_midpoints[peak_indexes]

    result = {
        'peak_indexes': peak_indexes,
        'histbins': histbins,
        'hist': hist,
        'bin_midpoints': bin_midpoints,
        'peak_values': peak_values,
    }
    return result


def get_density_peaks(data, weights=None, xs=None, threshold=None, bw_method=None, invert=False):
    density = scipy.stats.gaussian_kde(data, weights=weights, bw_method=bw_method)
    if xs is None:
        xs = np.linspace(data.min(), data.max(), 100)
    ys = density(xs)

    peak_indexes = get_1d_peaks(ys, invert=invert)
    if threshold is not None:
        peak_indexes = peak_indexes[ys[peak_indexes] > threshold]

    peak_values = xs[peak_indexes]
    peak_densities = ys[peak_indexes]

    return peak_values, peak_densities, density


def get_df_diagonals(df, upslope=True):
    """Sorting of diagonal rows: top-down
    Sorting of values in a row: left-right
    """
    def get_start_coord_upslope(n, nrow):
        if n < nrow:
            return (n, 0)
        else:
            return (nrow - 1, n - (nrow - 1))
        
    def get_diags_upslope(df, N, nrow, ncol):
        result = list()
        for n in range(N):
            pairs = list() 
            coord = get_start_coord_upslope(n, nrow)
            while True:
                if (coord[0] < 0) or (coord[1] >= ncol):
                    break

                val = df.iloc[coord[0], coord[1]]
                pairs.append((coord, val))
                coord = (coord[0] - 1, coord[1] + 1)
            result.append(pairs)
        return result

    def get_start_coord_downslope(n, ncol):
        if n < ncol:
            return (0, ncol - 1 - n)
        else:
            return (n - (ncol - 1), 0)
        
    def get_diags_downslope(df, N, nrow, ncol):
        result = list()
        for n in range(N):
            pairs = list() 
            coord = get_start_coord_downslope(n, ncol)
            while True:
                if (coord[0] >= nrow) or (coord[1] >= ncol):
                    break

                val = df.iloc[coord[0], coord[1]]
                pairs.append((coord, val))
                coord = (coord[0] + 1, coord[1] + 1)
            result.append(pairs)
        return result

    N = sum(df.shape) - 1
    nrow = df.shape[0]
    ncol = df.shape[1]

    if upslope:
        return get_diags_upslope(df, N, nrow, ncol)
    else:
        return get_diags_downslope(df, N, nrow, ncol)


def get_peak_coords_from_diagonals(diagonals, invert=False):
    peak_coords = list()
    for diag in diagonals:
        peak_indexes = get_1d_peaks([x[1] for x in diag], invert=invert)
        peak_coords.extend(diag[idx][0] for idx in peak_indexes)
    return peak_coords


def get_df_diagonal_peaks(df, invert=False):
    upslope_diags = get_df_diagonals(df, upslope=True)
    upslope_peak_coords = get_peak_coords_from_diagonals(upslope_diags, invert=invert)

    downslope_diags = get_df_diagonals(df, upslope=False)
    downslope_peak_coords = get_peak_coords_from_diagonals(downslope_diags, invert=invert)
        
    return upslope_peak_coords, downslope_peak_coords


def find_df_peaks(df, invert=False):
    row_peaks = df.apply(lambda x: get_1d_peaks(x, invert=invert), axis=1)
    row_peak_indexes = list()
    for row_idx, peaks in enumerate(row_peaks):
        row_peak_indexes.extend((row_idx, col_idx) for col_idx in peaks)

    col_peaks = df.apply(lambda x: get_1d_peaks(x, invert=invert), axis=0)
    col_peak_indexes = list()
    for col_idx, peaks in enumerate(col_peaks):
        col_peak_indexes.extend((row_idx, col_idx) for row_idx in peaks)

    return list(set(row_peak_indexes).intersection(set(col_peak_indexes)))


def find_df_peaks_4directions(df, invert=False):
    row_peaks = df.apply(lambda x: get_1d_peaks(x, invert=invert), axis=1)
    row_peak_indexes = list()
    for row_idx, peaks in enumerate(row_peaks):
        row_peak_indexes.extend((row_idx, col_idx) for col_idx in peaks)

    col_peaks = df.apply(lambda x: get_1d_peaks(x, invert=invert), axis=0)
    col_peak_indexes = list()
    for col_idx, peaks in enumerate(col_peaks):
        col_peak_indexes.extend((row_idx, col_idx) for row_idx in peaks)

    upslope_peak_coords, downslope_peak_coords = get_df_diagonal_peaks(df, invert=invert)

    return set.intersection(
        set(row_peak_indexes), set(col_peak_indexes),
        set(upslope_peak_coords), set(downslope_peak_coords),
    )


def find_df_peak_cpvalues(df, invert=False):
    assert df.index.name == 'cellularity'
    assert df.columns.name == 'ploidy'

    #peaks = find_df_peaks(df, invert=invert)
    peaks = find_df_peaks_4directions(df, invert=invert)
    cpvalues = [CPPair(df.index[x], df.columns[y]) for (x, y) in peaks]
    return sorted(cpvalues)


def get_peak_info(cp_score_dict, key='segfit', invert=True):
    dfs = make_cpscore_dfs(cp_score_dict)
    peaks_cpvalues = find_df_peak_cpvalues(dfs[key], invert=invert)
    peak_values = list()
    for c, p in peaks_cpvalues:
        data = {
            k: v for k, v in cp_score_dict[(c, p)].items()
            if k != 'CNt_list'
        }
        data['cellularity'] = c
        data['ploidy'] = p
        peak_values.append(data)
    #peak_values = sorted(peak_values, key=(lambda x: x['segfit_score']))

    return peak_values, dfs


###############
# CNn related #
###############

@functools.cache
def get_CNn_gr(refver, is_female):
    """Only accepts PAR-available reference versions"""
    chromdict = common.DEFAULT_CHROMDICTS[refver]
    autosomal_chroms = [
        x for x in chromdict.contigs 
        if (
            (common.RE_PATS['assembled_chromosome'].fullmatch(x) is not None)
            and (x not in MALE_HAPLOID_CHROMS)
        )
    ]  # 1, 2, ...,
    X_chrom, Y_chrom = chromdict.XY_names
    par_gr = libassemblyspec.get_par_gr(refver)

    # autosomal
    all_chrom_gr = chromdict.to_gr()
    autosomal_gr = all_chrom_gr[all_chrom_gr.Chromosome.isin(autosomal_chroms)]
    autosomal_gr.CNn = 2
    # sex
    if is_female:
        sex_gr = all_chrom_gr[all_chrom_gr.Chromosome.isin([X_chrom, Y_chrom])]
        sex_gr.CNn = 2
    else:
        all_sex_gr = all_chrom_gr[all_chrom_gr.Chromosome.isin([X_chrom, Y_chrom])]
        nonpar_sex_gr = all_sex_gr.subtract(par_gr)
        nonpar_sex_gr.CNn = 1
        par_sex_gr = par_gr.copy()
        par_sex_gr.CNn = [(2 if x == X_chrom else 0) for x in par_sex_gr.Chromosome]
        sex_gr = pr.concat([nonpar_sex_gr, par_sex_gr]).sort()
    # result
    return pr.concat([autosomal_gr, sex_gr]).sort()


def add_CNn_to_gr(gr, refver, is_female):
    return pyranges_helper.join(
        gr, get_CNn_gr(refver, is_female), how='left', merge='first',
        as_gr=True,
    )


def get_normal_mean_ploidy(refver, is_female, target_region_gr=None):
    """Only supports hg19 or hg38"""
    if is_female:
        return 2
    else:
        # sanity check
        if target_region_gr is not None:
            assert 'CNn' not in target_region_gr.columns

        chromdict = common.DEFAULT_CHROMDICTS[refver]
        if target_region_gr is None:
            N_region_gr = libassemblyspec.get_N_region_gr(refver)
            target_region_gr = chromdict.to_gr(assembled_only=True, as_gr=True).subtract(N_region_gr)

        CNn_gr = get_CNn_gr(refver, is_female)
        target_region_gr = pyranges_helper.join(
            target_region_gr, CNn_gr, how='left', merge='longest', as_gr=True,
        )

        indexer = target_region_gr.CNn.notna()
        CNn_list = target_region_gr.CNn.array[indexer]
        weights = target_region_gr.lengths().array[indexer]
        return np.average(CNn_list, weights=weights)


################
# miscellanous #
################

def check_haploid(is_female, chrom):
    return (not is_female) and (chrom in MALE_HAPLOID_CHROMS)


def get_quantile_df(df):
    return df.applymap(lambda x: (df < x).to_numpy().mean())


def get_bafs(vaf_list):
    bafs = np.array(vaf_list)
    gt_half_indexes = bafs > 0.5
    bafs[gt_half_indexes] = 1 - bafs[gt_half_indexes]
    return bafs


def flip_gt_half(array):
    gt_half_indexes = array > 0.5
    array[gt_half_indexes] = 1 - array[gt_half_indexes]
    return array


def get_nearest_integers(x):
    x = np.atleast_1d(x)

    lower = np.floor(x)
    upper = np.ceil(x)
    same_ind = (upper == lower)
    upper[upper == lower] += 1

    #lower = np.squeeze(lower)
    #upper = np.squeeze(upper)

    return upper, lower


##################################################
# theoretical value calculation helper functions # 
##################################################

def baf_function_wrapper(func):
    sig = inspect.signature(func)
    params = sig.parameters
    assert 'CNn' in params, f'There must be an argument named "CNn"'

    @functools.wraps(func)
    def wrapper(*args, **kwargs): 
        # setup
        ba = sig.bind(*args, **kwargs)
        ba.apply_defaults()
        #keys = list(ba.arguments.keys())
        #vals = [np.atleast_1d(x) for x in ba.arguments.values()]
        #argdict = dict(zip(keys, vals))

        ## get mask for invalid CNn values
        #mask = (argdict['CNn'] < 2)

        # replace invalid CNn values with nan
        argdict = ba.arguments
        argdict['CNn'] = np.where(argdict['CNn'] < 2, np.nan, argdict['CNn'])

        # main
        result = func(**argdict)

        return result

    return wrapper


def baf_function_wrapper_old(func):
    sig = inspect.signature(func)
    params = sig.parameters
    assert 'CNn' in params, f'There must be an argument named "CNn"'

    @functools.wraps(func)
    def wrapper(*args, **kwargs): 
        ba = sig.bind(*args, **kwargs)
        ba.apply_defaults()
        keys = list(ba.arguments.keys())
        vals = common.broadcast_args(list(ba.arguments.values()))
        argdict = dict(zip(keys, vals))

        CNn = argdict['CNn']
        result = np.full(CNn.shape, np.nan)
        valid_CNn_indexes = CNn[(0,) * (CNn.ndim - 1) + (slice(None, None),)] >= 2
        valid_argdict = {
            key: val[..., valid_CNn_indexes]
            for key, val in argdict.items()
        }

        # main
        result[..., valid_CNn_indexes] = func(**valid_argdict)

        # reduce dimension
        if isinstance(result, (list, tuple)):
            result = tuple(
                x[0] if x.shape == (1,) else x
                for x in result
            )
        else:
            if result.shape == (1,):
                result = result[0]

        return result

    return wrapper


def get_normal_norm_depth(normal_avg_depth, CNn, normal_ploidy):
    return normal_avg_depth * (CNn / normal_ploidy)


def depthratio_tumor_numerator(cellularity, CNt, CNn):
    '''
    <Derivation>
    tumor_upper_term = (CNt * cellularity) + (CNn * (1 - cellularity))
    tumor_upper_term = cellularity * (CNt - CNn) + CNn
    '''
    return cellularity * (CNt - CNn) + CNn


def depthratio_tumor_denominator(cellularity, tumor_ploidy, normal_ploidy):
    '''
    <Derivation>
    tumor_lower_term = (tumor_ploidy * cellularity) + (normal_ploidy * (1 - cellularity))
    tumor_lower_term = cellularity * (tumor_ploidy - normal_ploidy) + normal_ploidy
    '''
    return cellularity * (tumor_ploidy - normal_ploidy) + normal_ploidy


def depthratio_tumor_numerator_subclone(cellularity, ccf, clonal_CNt, subclonal_CNt, CNn):
    '''
    <Derivation>
    tumor_upper_term = (
        (clonal_CNt * cellularity * (1 - ccf)) 
        + (subclonal_CNt * cellularity * ccf) 
        + (CNn * (1 - cellularity))
    ) 
    tumor_upper_term = (
        ccf * (subclonal_CNt * cellularity - clonal_CNt * cellularity)
        + clonal_CNt * cellularity
        + (CNn * (1 - cellularity))
    )
    tumor_upper_term = (
        ccf * cellularity * (subclonal_CNt - clonal_CNt)
        + clonal_CNt * cellularity
        + (CNn - cellularity * CNn)
    )
    tumor_upper_term = (
        cellularity * (
            ccf * (subclonal_CNt - clonal_CNt)
            + clonal_CNt
            - CNn
        )
        + CNn
    )
    '''
    return (
        cellularity * (
            ccf * (subclonal_CNt - clonal_CNt)
            + (clonal_CNt - CNn)
        )
        + CNn
    )


def baf_numerator_subclone(cellularity, ccf, clonal_B, subclonal_B, Bn=1):
    '''
    <Derivation>
    numerator = (
        (clonal_B * cellularity * (1 - ccf))
        + (subclonal_B * cellularity * ccf)
        + Bn * (1 - cellularity)
    )
    numerator = (
        ccf * (subclonal_B * cellularity - clonal_B * cellularity)
        + clonal_B * cellularity
        + Bn * (1 - cellularity)
    )
    numerator = (
        ccf * cellularity * (subclonal_B - clonal_B)
        + cellularity * (clonal_B - Bn)
        + Bn
    )
    numerator = (
        cellularity * (
            ccf * (subclonal_B - clonal_B)
            + (clonal_B - Bn)
        )
        + Bn
    )
    '''
    return (
        cellularity * (
            ccf * (subclonal_B - clonal_B)
            + (clonal_B - Bn)
        )
        + Bn
    )


#################################
# theoretical value calculation # 
#################################

@deco.get_deco_asarray(['CNt', 'CNn'])
def theoretical_depth_ratio(
    CNt, cellularity, tumor_ploidy, CNn, normal_ploidy, 
    tumor_avg_depth=1, 
    normal_avg_depth=1,
):
    """Args:
        CNn, normal_ploidy must be nonzero
    Can be vectorized
    """
    tumor_norm_depth = tumor_avg_depth * (
        depthratio_tumor_numerator(cellularity, CNt, CNn)
        / depthratio_tumor_denominator(cellularity, tumor_ploidy, normal_ploidy)
    )
    normal_norm_depth = normal_avg_depth * (CNn / normal_ploidy)

    return tumor_norm_depth / normal_norm_depth


@deco.get_deco_broadcast(['depth_ratio', 'CNn'])
def inverse_theoretical_depth_ratio(
    depth_ratio, cellularity, tumor_ploidy, CNn, normal_ploidy, 
    tumor_avg_depth=1, 
    normal_avg_depth=1,
):
    """Args:
        cellularity must be nonzero
    Returns:
        CNt estimate
    """
    '''
    <Derivation>
    tumor_upper_term = (CNt * cellularity) + (CNn * (1 - cellularity))
    tumor_lower_term = (tumor_ploidy * cellularity) + (normal_ploidy * (1 - cellularity))
    tumor_norm_depth = tumor_avg_depth * (tumor_upper_term / tumor_lower_term)
    normal_norm_depth = normal_avg_depth * (CNn / normal_ploidy)
    depth_ratio = tumor_norm_depth / normal_norm_depth
    ---------------------------------
    tumor_norm_depth = depth_ratio * normal_norm_depth
    tumor_avg_depth * (tumor_upper_term / tumor_lower_term) = depth_ratio * normal_norm_depth
    tumor_upper_term = (depth_ratio * normal_norm_depth * tumor_lower_term) / tumor_avg_depth
    (CNt * cellularity) + (CNn * (1 - cellularity)) = (depth_ratio * normal_norm_depth * tumor_lower_term) / tumor_avg_depth
    CNt * cellularity = ((depth_ratio * normal_norm_depth * tumor_lower_term) / tumor_avg_depth) - (CNn * (1 - cellularity))
    CNt = (
        ((depth_ratio * normal_norm_depth * tumor_lower_term) / tumor_avg_depth) 
        - (CNn * (1 - cellularity))
    ) / cellularity
    CNt = (
        (depth_ratio * ((normal_norm_depth * tumor_lower_term) / tumor_avg_depth)) 
        - (CNn * (1 - cellularity))
    ) / cellularity
    CNt = (
        depth_ratio * (
            (normal_norm_depth * tumor_lower_term) 
            / (tumor_avg_depth * cellularity)
        )
        - (CNn * (1 / cellularity - 1))
    )
    CNt = (
        depth_ratio * (
            (normal_avg_depth * (CNn / normal_ploidy) * tumor_lower_term) 
            / (tumor_avg_depth * cellularity)
        )
        - (CNn * (1 / cellularity - 1))
    )
    CNt = (
        depth_ratio * (
            (normal_avg_depth * CNn * tumor_lower_term) 
            / (tumor_avg_depth * cellularity * normal_ploidy)
        )
        - (CNn * (1 / cellularity - 1))
    )
    CNt = CNn * (
        depth_ratio * (
            (normal_avg_depth * tumor_lower_term) 
            / (tumor_avg_depth * cellularity * normal_ploidy)
        )
        - (1 / cellularity - 1)
    )
    CNt = CNn * (
        depth_ratio * (
            (
                normal_avg_depth * (
                    (tumor_ploidy * cellularity) + (normal_ploidy * (1 - cellularity))
                )
            )
            / (tumor_avg_depth * cellularity * normal_ploidy)
        )
        - (1 / cellularity - 1)
    )
    CNt = CNn * (
        depth_ratio * (
            (
                normal_avg_depth * (
                    tumor_ploidy + (normal_ploidy * (1 / cellularity - 1))
                )
            )
            / (tumor_avg_depth * normal_ploidy)
        )
        - (1 / cellularity - 1)
    )
    CNt = depth_ratio * (
        CNn * (
            (
                normal_avg_depth * (
                    tumor_ploidy + (normal_ploidy * (1 / cellularity - 1))
                )
            )
            / (tumor_avg_depth * normal_ploidy)
        )
    ) - CNn * (1 / cellularity - 1)
    CNt = depth_ratio * (
        (
            CNn 
            * normal_avg_depth 
            * (tumor_ploidy + (normal_ploidy * (1 / cellularity - 1)))
        )
        / (
            tumor_avg_depth * normal_ploidy
        )
    ) - CNn * (1 / cellularity - 1)
    CNt = CNn * (
        depth_ratio * (
            (
                normal_avg_depth 
                * (tumor_ploidy + (normal_ploidy * (1 / cellularity - 1)))
            )
            / (
                tumor_avg_depth * normal_ploidy
            )
        ) - (1 / cellularity - 1)
    )
    '''
    K = (1 / cellularity - 1)
    CNt = CNn * (
        depth_ratio * (
            (
                normal_avg_depth 
                * (tumor_ploidy + (normal_ploidy * K))
            )
            / (
                tumor_avg_depth * normal_ploidy
            )
        ) - K
    )

    return CNt


def delta_depth_ratio(
    cellularity, tumor_ploidy, CNn, normal_ploidy, 
    tumor_avg_depth=1, normal_avg_depth=1,
):
    """Returns:
        Difference of depth ratio caused by 1 CNt change
    """
    '''
    <Derivation>
    tumor_upper_term = (CNt * cellularity) + (CNn * (1 - cellularity))
    tumor_lower_term = (tumor_ploidy * cellularity) + (normal_ploidy * (1 - cellularity))
    tumor_norm_depth = tumor_avg_depth * (tumor_upper_term / tumor_lower_term)
    normal_norm_depth = normal_avg_depth * (CNn / normal_ploidy)
    depth_ratio = tumor_norm_depth / normal_norm_depth
    ------------------------------------------
    depth_ratio = (tumor_avg_depth * (tumor_upper_term / tumor_lower_term)) / normal_norm_depth
    depth_ratio = tumor_upper_term * (tumor_avg_depth / (tumor_lower_term * normal_norm_depth))
    depth_ratio = (
        ((CNt * cellularity) + (CNn * (1 - cellularity))) 
        * 
        (tumor_avg_depth / (tumor_lower_term * normal_norm_depth))
    )
    delta_depth_ratio = (
        ((CNt_1 * cellularity) - (CNt_2 * cellularity)) 
        * 
        (tumor_avg_depth / (tumor_lower_term * normal_norm_depth))
    )
    delta_depth_ratio = (
        cellularity
        * 
        (tumor_avg_depth / (tumor_lower_term * normal_norm_depth))
    )
    delta_depth_ratio = (cellularity * tumor_avg_depth) / (tumor_lower_term * normal_norm_depth)
    '''
    tumor_lower_term = depthratio_tumor_denominator(cellularity, tumor_ploidy, normal_ploidy)
    normal_norm_depth = normal_avg_depth * (CNn / normal_ploidy)
    delta_depth_ratio = (
        (cellularity * tumor_avg_depth) 
        / (tumor_lower_term * normal_norm_depth)
    )

    return delta_depth_ratio


def get_cp_from_twodata(depthratio1, CNt1, depthratio2, CNt2, normal_ploidy, CNn, tumor_avg_depth=1, normal_avg_depth=1):
    '''
    <Derivation>
    tumor_upper_term = (CNt * cellularity) + (CNn * (1 - cellularity))
    tumor_lower_term = (tumor_ploidy * cellularity) + (normal_ploidy * (1 - cellularity))
    tumor_norm_depth = tumor_avg_depth * (tumor_upper_term / tumor_lower_term)
    normal_norm_depth = normal_avg_depth * (CNn / normal_ploidy)
    depthratio = tumor_norm_depth / normal_norm_depth
    ------------------------------------------
    cellularity = term1 / (tumor_ploidy + term2)
    '''
    term1_data1 = get_cp_from_twodata_term1(depthratio1, normal_ploidy, tumor_avg_depth, normal_avg_depth)
    term1_data2 = get_cp_from_twodata_term1(depthratio2, normal_ploidy, tumor_avg_depth, normal_avg_depth)

    term2_data1 = get_cp_from_twodata_term2(depthratio1, CNt1, CNn, normal_ploidy, tumor_avg_depth, normal_avg_depth)
    term2_data2 = get_cp_from_twodata_term2(depthratio2, CNt2, CNn, normal_ploidy, tumor_avg_depth, normal_avg_depth)

    tumor_ploidy = (term2_data2 * term1_data1 - term2_data1 * term1_data2) / (term1_data2 - term1_data1)
    cellularity = term1_data1 / (tumor_ploidy + term2_data1)

    if not (
        (cellularity > 0)
        and (cellularity < 1)
        and (tumor_ploidy > 0)
    ):
        valid_CNt_pairs = get_valid_CNts_from_depthratios(
            depthratio1, 
            depthratio2, 
            normal_ploidy, 
            CNn, 
            tumor_avg_depth=tumor_avg_depth, 
            normal_avg_depth=normal_avg_depth, 
            CNt1_range=range(0, 6), 
            CNt2_maxdiff=6,
        )
        raise Exception(
            f'Input depthratio and CNt values resulted in invalid cellularity/ploidy values:\n'
            f'cellularity={cellularity}, ploidy={tumor_ploidy}\n'
            f'Valid (CNt1, CNt2) pairs for given depthratio values are:\n{valid_CNt_pairs}'
        )

    return cellularity, tumor_ploidy


def get_valid_CNts_from_depthratios(depthratio1, depthratio2, normal_ploidy, CNn, tumor_avg_depth=1, normal_avg_depth=1, CNt1_range=range(0, 6), CNt2_maxdiff=6):
    '''Args:
        depthratio1: lesser one
        depthratio2: greater one

    <Sanity condition>
    For cellularity to be between 0 and 1, following conditions are required:
        (depthratio1 designates the smaller one of two depthratio values)
        1) term2_data1 > term2_data2
        2) term2_data1 - term2_data2 > term1_data1 - term1_data2

        #1) (depthratio2 / depthratio1) < (CNt2 / CNt1)
        #2) (depthratio2 / depthratio1) < ((CNt2 - CNn) / (CNt1 - CNn))
    For ploidy to be greater than 0, following condition is required:
        term2_data1 * term1_data2 > term2_data2 * term1_data1
    '''
    if depthratio2 <= depthratio1:
        raise Exception(f'"depthratio2" must be greater than "depthratio1"')

    drr = depthratio2 / depthratio1
    result = list()
    for CNt1 in CNt1_range:
        for CNt2 in range(CNt1 + 1, CNt1 + CNt2_maxdiff + 1):

            term1_data1 = get_cp_from_twodata_term1(depthratio1, normal_ploidy, tumor_avg_depth, normal_avg_depth)
            term1_data2 = get_cp_from_twodata_term1(depthratio2, normal_ploidy, tumor_avg_depth, normal_avg_depth)

            term2_data1 = get_cp_from_twodata_term2(depthratio1, CNt1, CNn, normal_ploidy, tumor_avg_depth, normal_avg_depth)
            term2_data2 = get_cp_from_twodata_term2(depthratio2, CNt2, CNn, normal_ploidy, tumor_avg_depth, normal_avg_depth)

            if (
                (term2_data1 > term2_data2)
                and (term2_data1 - term2_data2 > term1_data1 - term1_data2)
                and (term2_data1 * term1_data2 > term2_data2 * term1_data1)
            ):
                result.append((CNt1, CNt2))

    return result


def get_cp_from_twodata_term1(depthratio, normal_ploidy, tumor_avg_depth, normal_avg_depth):
    return normal_ploidy * (tumor_avg_depth / (depthratio * normal_avg_depth) - 1)


def get_cp_from_twodata_term2(depthratio, CNt, CNn, normal_ploidy, tumor_avg_depth, normal_avg_depth):
    return (
        -1 
        * normal_ploidy 
        * (
            1 
            + (
                (tumor_avg_depth * (CNt - CNn)) 
                / (CNn * normal_avg_depth * depthratio)
            )
        )
    )


def get_c_from_ddr_p(delta_depth_ratio, tumor_ploidy, normal_ploidy, CNn, tumor_avg_depth=1, normal_avg_depth=1):
    """c: cellularity; ddr: delta depthratio; p: ploidy
    """
    '''
    <Derivation>
    tumor_lower_term = (tumor_ploidy * cellularity) + (normal_ploidy * (1 - cellularity))
    normal_norm_depth = normal_avg_depth * (CNn / normal_ploidy)
    delta_depth_ratio = (cellularity * tumor_avg_depth) / (tumor_lower_term * normal_norm_depth)

    delta_depth_ratio * tumor_lower_term = (cellularity * tumor_avg_depth) / normal_norm_depth

    left term
        = delta_depth_ratio * ( (tumor_ploidy * cellularity) + (normal_ploidy * (1 - cellularity)) )
        = cellularity * (delta_depth_ratio * tumor_ploidy) + (1 - cellularity) * (delta_depth_ratio * normal_ploidy)
        = cellularity * (delta_depth_ratio * tumor_ploidy - delta_depth_ratio * normal_ploidy) + (delta_depth_ratio * normal_ploidy)
        = cellularity * delta_depth_ratio * (tumor_ploidy - normal_ploidy) + (delta_depth_ratio * normal_ploidy)
    right term
        = (cellularity * tumor_avg_depth) / ( normal_avg_depth * (CNn / normal_ploidy) )
        = cellularity * (tumor_avg_depth / ( normal_avg_depth * (CNn / normal_ploidy) ) )
        = cellularity * ((tumor_avg_depth * normal_ploidy) / (normal_avg_depth * CNn))

    cellularity * (
        ((tumor_avg_depth * normal_ploidy) / (normal_avg_depth * CNn))
        - delta_depth_ratio * (tumor_ploidy - normal_ploidy)
    ) = (delta_depth_ratio * normal_ploidy)
    ...
    numer = 1
    denom = (
        tumor_avg_depth / (delta_depth_ratio * normal_avg_depth * CNn)
        - (tumor_ploidy / normal_ploidy)
        + 1
    )
    '''
    denom = (
        tumor_avg_depth / (delta_depth_ratio * normal_avg_depth * CNn)
        - (tumor_ploidy / normal_ploidy)
        + 1
    )
    return 1 / denom


def theoretical_depth_ratio_sequenza(CNt, cellularity, tumor_ploidy, CNn=2, normal_ploidy=2, avg_depth_ratio=1):
    """Args:
        CNn, normal_ploidy must be nonzero
    """
    cellularity_term = ((CNt / CNn) * cellularity) + (1 - cellularity)
    ploidy_term = ((tumor_ploidy / normal_ploidy) * cellularity) + (1 - cellularity)
    return avg_depth_ratio * (cellularity_term / ploidy_term)


def inverse_theoretical_depth_ratio_sequenza(depth_ratio, cellularity, tumor_ploidy, CNn=2, normal_ploidy=2, avg_depth_ratio=1):
    """Args:
        cellularity must be nonzero
    Returns:
        CNt estimate
    """
    ploidy_term = ((tumor_ploidy / normal_ploidy) * cellularity) + (1 - cellularity)
    # depth_ratio * ploidy_term == avg_depth_ratio * cellularity_term
    return (((depth_ratio * ploidy_term) / avg_depth_ratio) - (1 - cellularity)) * (CNn / cellularity)


def delta_depth_ratio_sequenza(cellularity, tumor_ploidy, CNn=2, normal_ploidy=2, avg_depth_ratio=1):
    ploidy_term = ((tumor_ploidy / normal_ploidy) * cellularity) + (1 - cellularity)
    return (avg_depth_ratio / ploidy_term) * (cellularity / CNn)


@baf_function_wrapper
def theoretical_baf(CNt, B, cellularity, CNn, Bn):
    '''
    <Derivation>
    numerator = (B * cellularity) + Bn * (1 - cellularity)
    numerator = cellularity * (B - Bn) + Bn

    denominator = (CNt * cellularity) + (CNn * (1 - cellularity))
    denominator = cellularity * (CNt - CNn) + CNn
    '''
    return (
        (cellularity * (B - Bn) + Bn)
        / (cellularity * (CNt - CNn) + CNn)
    )


@baf_function_wrapper
def inverse_theoretical_baf(baf, cellularity, CNt, CNn, Bn):
    """Args:
        cellularity must be nonzero
    Returns:
        B estimate
    """
    '''
    numerator = (B * cellularity) + Bn * (1 - cellularity)
    denominator = (CNt * cellularity) + (CNn * (1 - cellularity))
    baf = numerator / denominator

    ###
    baf = ((B * cellularity) + Bn * (1 - cellularity)) / denominator
    (B * cellularity) + Bn * (1 - cellularity) = baf * denominator
    B * cellularity = (baf * denominator) - Bn * (1 - cellularity)
    B = ((baf * denominator) - Bn * (1 - cellularity)) / cellularity
    B = (
        (baf * denominator) - Bn * (1 - cellularity)
    ) / cellularity
    B = (
        (baf * ((CNt * cellularity) + (CNn * (1 - cellularity)))) - Bn * (1 - cellularity)
    ) / cellularity
    B = (
        (baf * CNt * cellularity + baf * CNn * (1 - cellularity)) - Bn * (1 - cellularity)
    ) / cellularity
    B = (
        baf * CNt * cellularity + baf * CNn * (1 - cellularity) - Bn * (1 - cellularity)
    ) / cellularity
    B = (
        baf * CNt * cellularity + (1 - cellularity) * (baf * CNn - Bn)
    ) / cellularity
    B = baf * CNt + (1 / cellularity - 1) * (baf * CNn - Bn)
    '''
    return baf * CNt + (1 / cellularity - 1) * (baf * CNn - Bn)


#############################################
# find an optimal copy number from cp value # 
#############################################

def get_B(CNt, cellularity, baf, CNn):
    max_B = int(CNt / 2)
    B_estimate = inverse_theoretical_baf(baf, cellularity, CNt, CNn)
    if B_estimate <= 0:
        B = 0
        theo_baf = theoretical_baf(CNt, B, cellularity, CNn)
        if np.isnan(theo_baf):
            B = np.nan
            diff = np.nan
        else:
            diff = abs(baf - theo_baf)
    else:
        B_candidate_upper = int(np.ceil(B_estimate))
        B_candidate_lower = int(np.floor(B_estimate))

        theo_baf_upper = (
            theoretical_baf(CNt, B_candidate_upper, cellularity, CNn)
            if B_candidate_upper <= max_B else
            np.nan
        )
        theo_baf_lower = (
            theoretical_baf(CNt, B_candidate_lower, cellularity, CNn)
            if B_candidate_lower <= max_B else
            np.nan
        )

        if np.isnan(theo_baf_upper) and np.isnan(theo_baf_lower):
            B = np.nan
            diff = np.nan
        elif (not np.isnan(theo_baf_upper)) and np.isnan(theo_baf_lower):
            B = B_candidate_upper
            diff = abs(baf - theo_baf_upper)
        elif np.isnan(theo_baf_upper) and (not np.isnan(theo_baf_lower)):
            B = B_candidate_lower
            diff = abs(baf - theo_baf_lower)
        else:
            diff_upper = abs(theo_baf_upper - baf)
            diff_lower = abs(theo_baf_lower - baf)
            if diff_upper <= diff_lower:
                B = B_candidate_upper
                diff = diff_upper
            else:
                B = B_candidate_lower
                diff = diff_lower

    return B, diff


def get_CN_from_cp(cellularity, tumor_ploidy, depth_ratio, baf, CNt_weight, CNn, normal_ploidy):
    def save_cache(CNt_candidate, B_cache, segfit_score_cache, baf_diff_cache):
        B, baf_diff = get_B(CNt_candidate, cellularity, baf, CNn)
        if np.isnan(B):
            # drop this CNt candidate if B cannot be calculated
            return

        ratio_diff = abs(
            theoretical_depth_ratio(CNt_candidate, cellularity, tumor_ploidy, CNn, normal_ploidy) 
            - depth_ratio
        )
        B_cache[CNt_candidate] = B
        baf_diff_cache[CNt_candidate] = baf_diff
        segfit_score_cache[CNt_candidate] = ratio_diff * CNt_weight + baf_diff

    # get initial CNt candidates
    CNt_estimate = inverse_theoretical_depth_ratio(
        depth_ratio, cellularity, tumor_ploidy, CNn, normal_ploidy,
    )
    if CNt_estimate <= 0:
        initial_candidate_upper = 0
        initial_candidate_lower = None
    else:
        initial_candidate_upper = int(np.ceil(CNt_estimate))
        initial_candidate_lower = int(np.floor(CNt_estimate))

    # set caches
    B_cache = dict()
    segfit_score_cache = dict()
    baf_diff_cache = dict()

    delta_ratio = delta_depth_ratio(cellularity, tumor_ploidy, CNn, normal_ploidy) * CNt_weight

    # upper
    save_cache(initial_candidate_upper, B_cache, segfit_score_cache, baf_diff_cache)
    if initial_candidate_upper in baf_diff_cache:
        n_further_candidates = int(baf_diff_cache[initial_candidate_upper] / delta_ratio)
        for offset in range(n_further_candidates):
            CNt_candidate = initial_candidate_upper + (offset + 1)
            save_cache(CNt_candidate, B_cache, segfit_score_cache, baf_diff_cache)
    # lower
    if initial_candidate_lower is not None:
        save_cache(initial_candidate_lower, B_cache, segfit_score_cache, baf_diff_cache)
        if initial_candidate_lower in baf_diff_cache:
            try:
                n_further_candidates = int(baf_diff_cache[initial_candidate_lower] / delta_ratio)
            except:
                print(cellularity, tumor_ploidy, depth_ratio, baf)
                print(baf, baf_diff_cache[initial_candidate_lower], initial_candidate_lower, delta_ratio)
                raise

            for offset in range(n_further_candidates):
                CNt_candidate = initial_candidate_lower - (offset + 1)
                if CNt_candidate < 0:
                    break
                save_cache(CNt_candidate, B_cache, segfit_score_cache, baf_diff_cache)

    # result
    if len(segfit_score_cache) == 0:
        return np.nan, np.nan, np.nan
    else:
        CNt, segfit_score = min(segfit_score_cache.items(), key=operator.itemgetter(1))
        B = B_cache[CNt]
        #selected_score = segfit_score_cache[CNt]

        return CNt, B, segfit_score


def get_CN_from_cp_wo_baf(cellularity, tumor_ploidy, depth_ratio, CNt_weight, CNn, normal_ploidy):
    CNt_estimate = inverse_theoretical_depth_ratio(
        depth_ratio, cellularity, tumor_ploidy, CNn, normal_ploidy,
    )
    if CNt_estimate <= 0:
        CNt = 0
        diff = abs(
            theoretical_depth_ratio(CNt, cellularity, tumor_ploidy, CNn, normal_ploidy) 
            - depth_ratio
        )
        segfit_score = diff * CNt_weight
    else:
        upper_candidate = int(np.ceil(CNt_estimate))
        lower_candidate = int(np.floor(CNt_estimate))

        upper_diff = abs(
            theoretical_depth_ratio(upper_candidate, cellularity, tumor_ploidy, CNn, normal_ploidy) 
            - depth_ratio
        )
        lower_diff = abs(
            theoretical_depth_ratio(lower_candidate, cellularity, tumor_ploidy, CNn, normal_ploidy) 
            - depth_ratio
        )
        if upper_diff <= lower_diff:
            CNt = upper_candidate
            segfit_score = upper_diff * CNt_weight
        else:
            CNt = lower_candidate
            segfit_score = lower_diff * CNt_weight

    return CNt, np.nan, segfit_score 


#def calc_clonal_CNt(
#    depth_ratio, 
#    CNn, 
#    cellularity, 
#    tumor_ploidy, 
#    normal_ploidy,
#):
#    estimated_CNts = inverse_theoretical_depth_ratio(
#        depth_ratio, cellularity, tumor_ploidy, CNn, normal_ploidy, 
#        tumor_avg_depth=1, normal_avg_depth=1,
#    )
#    return estimated_CNts
#    #return np.rint(estimated_CNts).astype(int)


@deco.get_deco_broadcast(['depth_ratio', 'baf', 'CNn', 'Bn'])
def find_clonal_solution(
    depth_ratio,
    baf,
    CNn,
    cellularity,
    tumor_ploidy,
    normal_ploidy,
    Bn=1,
    average_CNt=None,
):
    # get CNt
    if average_CNt is None:
        average_CNt = inverse_theoretical_depth_ratio(
            depth_ratio=depth_ratio, 
            cellularity=cellularity, 
            tumor_ploidy=tumor_ploidy, 
            CNn=CNn, 
            normal_ploidy=normal_ploidy, 
            tumor_avg_depth=1, 
            normal_avg_depth=1,
        )
    else:
        average_CNt = np.asarray(average_CNt)
        assert average_CNt.shape == depth_ratio.shape

    CNts = np.clip(np.rint(average_CNt), 0, None)

    ## mask invalid CNt positions
    #CNt_mask = (CNts < 0)
    #CNts[CNt_mask] = np.nan

    # get B
    estimated_Bs = inverse_theoretical_baf(
        baf=baf,
        cellularity=cellularity, 
        CNt=CNts, 
        CNn=CNn, 
        Bn=Bn,
    )
    Bs = np.clip(
        np.rint(estimated_Bs), 
        0, 
        np.floor(0.5 * CNts),
    )

#    half_CNts = 0.5 * CNts
#    assert (estimated_Bs <= half_CNts).all()
#    Bs = np.rint(estimated_Bs)
#    Bs = np.minimum(Bs, np.floor(half_CNts))  # forces B to be less than A
#
#    B_mask = (Bs < 0)  # CNn < 2 positions are not included
#    Bs[B_mask] = np.nan

    # get theoretical values
    expected_depthratio = theoretical_depth_ratio(
        CNt=CNts, 
        cellularity=cellularity, 
        tumor_ploidy=tumor_ploidy, 
        CNn=CNn, 
        normal_ploidy=normal_ploidy, 
        tumor_avg_depth=1, 
        normal_avg_depth=1,
    )
    expected_baf = theoretical_baf(
        CNt=CNts, 
        B=Bs, 
        cellularity=cellularity, 
        CNn=CNn, 
        Bn=Bn,
    )

    # result
    result = dict()
    result['CNt'] = CNts
    result['B'] = Bs
    #result['invalid_CNt_flag'] = CNt_mask
        # where resulting CNt is less than 0
    #result['invalid_B_flag'] = B_mask
        # where resulting B is less than 0
        # invalid CNt and (CNn < 2) are not included
    result['theoretical_depth_ratio'] = expected_depthratio
    result['theoretical_baf'] = expected_baf
    result['depth_ratio_diff'] = expected_depthratio - depth_ratio
    result['baf_diff'] = expected_baf - baf

    return result


#####################################
# subclonal copy number calculation #
#####################################

def theoretical_depth_ratio_subclone(
    clonal_CNt, 
    subclonal_CNt,
    ccf,
    cellularity, 
    tumor_ploidy, 
    CNn, 
    normal_ploidy, 
    tumor_avg_depth=1, 
    normal_avg_depth=1,
    rescue_nan_ccf=False,
):
    """Can be vectorized
    """
    # argument modification
    if rescue_nan_ccf:
        ccf_isnan = np.isnan(ccf)
        ccf = np.where(ccf_isnan, 0, ccf)
        subclonal_CNt = np.where(ccf_isnan, 0, subclonal_CNt)
            # without this, subclonal_CNt values remain as nan

    # main
    tumor_upper_term = depthratio_tumor_numerator_subclone(
        cellularity, ccf, clonal_CNt, subclonal_CNt, CNn,
    )
    tumor_lower_term = depthratio_tumor_denominator(cellularity, tumor_ploidy, normal_ploidy)
    tumor_norm_depth = tumor_avg_depth * (tumor_upper_term / tumor_lower_term)

    normal_norm_depth = get_normal_norm_depth(normal_avg_depth, CNn, normal_ploidy)

    return tumor_norm_depth / normal_norm_depth


def get_ccf_from_depthratio(
    depth_ratio,
    cellularity,
    tumor_ploidy,
    normal_ploidy,
    clonal_CNt,
    subclonal_CNt,
    CNn,
    tumor_avg_depth=1,
    normal_avg_depth=1,
):
    """
    calculate ccf from copy number, cellularity, ploidy, etc.
    """

    '''
    <Derivation>
    depth_ratio = tumor_norm_depth / normal_norm_depth
    depth_ratio = (
        tumor_avg_depth * (tumor_upper_term / tumor_lower_term) 
        / normal_norm_depth
    )
    depth_ratio = (tumor_avg_depth * tumor_upper_term) / (tumor_lower_term * normal_norm_depth)
    tumor_upper_term = (depth_ratio * tumor_lower_term * normal_norm_depth) / tumor_avg_depth

    ###

    tumor_upper_term = (
        ccf * cellularity * (subclonal_CNt - clonal_CNt)
        + clonal_CNt * cellularity
        + (CNn - cellularity * CNn)
    )
    tumor_upper_term = (
        ccf * cellularity * (subclonal_CNt - clonal_CNt)
        + cellularity * (clonal_CNt - CNn)
        + CNn
    )

    ###

    (
        ccf * cellularity * (subclonal_CNt - clonal_CNt)
        + cellularity * (clonal_CNt - CNn)
        + CNn
    ) = (depth_ratio * tumor_lower_term * normal_norm_depth) / tumor_avg_depth
    ccf * cellularity * (subclonal_CNt - clonal_CNt) = (
        (depth_ratio * tumor_lower_term * normal_norm_depth) / tumor_avg_depth
        - cellularity * (clonal_CNt - CNn)
        - CNn
    )
    ccf = (
        (depth_ratio * tumor_lower_term * normal_norm_depth) / tumor_avg_depth
        - cellularity * (clonal_CNt - CNn)
        - CNn
    ) / (
        cellularity * (subclonal_CNt - clonal_CNt)
    )
    ccf = (
        (depth_ratio * tumor_lower_term * normal_norm_depth) / tumor_avg_depth
        - cellularity * (clonal_CNt - CNn)
        - CNn
    ) / (
        cellularity * delta_CNt
    )  (delta_CNt == subclonal_CNt - clonal_CNt)
    '''
    tumor_lower_term = depthratio_tumor_denominator(cellularity, tumor_ploidy, normal_ploidy)
    normal_norm_depth = get_normal_norm_depth(normal_avg_depth, CNn, normal_ploidy)
    delta_CNt = subclonal_CNt - clonal_CNt
    ccf = (
        (depth_ratio * tumor_lower_term * normal_norm_depth) / tumor_avg_depth
        - cellularity * (clonal_CNt - CNn)
        - CNn
    ) / (
        cellularity * delta_CNt
    )
    return ccf


def get_ccf_from_baf(
    baf,
    cellularity,
    clonal_CNt,
    subclonal_CNt,
    clonal_B,
    subclonal_B,
    CNn,
    Bn,
):
    '''
    <Derivation>
    numerator = (
        (clonal_B * cellularity * (1 - ccf))
        + (subclonal_B * cellularity * ccf)
        + Bn * (1 - cellularity)
    )
    denominator = (
        (clonal_CNt * cellularity * (1 - ccf)) 
        + (subclonal_CNt * cellularity * ccf) 
        + (CNn * (1 - cellularity))
    )
    baf = numerator / denominator

    ###

    numerator 
        = ccf * (
            subclonal_B * cellularity
            - clonal_B * cellularity
        ) + (
            clonal_B * cellularity
            + Bn * (1 - cellularity)
        )
        = ccf * cellularity * (subclonal_B - clonal_B)
        + (
            cellularity * (clonal_B - Bn) + Bn
        )

    ###

    denominator
        = ccf * (
            subclonal_CNt * cellularity
            - clonal_CNt * cellularity
        ) + (
            clonal_CNt * cellularity
            + (CNn * (1 - cellularity))
        )
        = ccf * cellularity * (subclonal_CNt - clonal_CNt)
        + (
            cellularity * (clonal_CNt - CNn) + CNn
        )

    ###

    baf * denominator = numerator

    (
        baf * ccf * cellularity * (subclonal_CNt - clonal_CNt)
        + baf * (
            cellularity * (clonal_CNt - CNn) + CNn
        )
    ) = (
        ccf * cellularity * (subclonal_B - clonal_B)
        + (
            cellularity * (clonal_B - Bn) + Bn
        )
    )

    (
        baf * ccf * cellularity * (subclonal_CNt - clonal_CNt)
        - ccf * cellularity * (subclonal_B - clonal_B)
    ) = (
        (cellularity * (clonal_B - Bn) + Bn)
        - baf * (cellularity * (clonal_CNt - CNn) + CNn)
    ) 

    ccf * cellularity * (
        baf * (subclonal_CNt - clonal_CNt)
        - (subclonal_B - clonal_B)
    ) = (
        (cellularity * (clonal_B - Bn) + Bn)
        - baf * (cellularity * (clonal_CNt - CNn) + CNn)
    )

    ccf = (
        (cellularity * (clonal_B - Bn) + Bn)
        - baf * (cellularity * (clonal_CNt - CNn) + CNn)
    ) / (
        cellularity * (
            baf * (subclonal_CNt - clonal_CNt)
            - (subclonal_B - clonal_B)
        )
    )

    ccf = (
        ((clonal_B - Bn) + Bn / cellularity)
        - baf * ((clonal_CNt - CNn) + CNn / cellularity)
    ) / (
        baf * (subclonal_CNt - clonal_CNt)
        - (subclonal_B - clonal_B)
    )
    
    '''
    return (
        ((clonal_B - Bn) + Bn / cellularity)
        - baf * ((clonal_CNt - CNn) + CNn / cellularity)
    ) / (
        baf * (subclonal_CNt - clonal_CNt)
        - (subclonal_B - clonal_B)
    )


@baf_function_wrapper
def theoretical_baf_subclone(
    clonal_CNt, 
    clonal_B, 
    subclonal_CNt, 
    subclonal_B, 
    ccf,
    cellularity, 
    CNn,
    Bn,
    rescue_nan_ccf=False,
):
    # arg modification
    if rescue_nan_ccf:
        ccf_isnan = np.isnan(ccf)
        ccf = np.where(ccf_isnan, 0, ccf)
        subclonal_CNt = np.where(ccf_isnan, 0, subclonal_CNt)
        subclonal_B = np.where(ccf_isnan, 0, subclonal_B)
            # without this, subclonal_CNt values remain as nan

    # main
    numerator = baf_numerator_subclone(cellularity, ccf, clonal_B, subclonal_B, Bn=Bn)
    denominator = depthratio_tumor_numerator_subclone(
        cellularity, ccf, clonal_CNt, subclonal_CNt, CNn,
    )
    result = numerator / denominator

    return result


@deco.args_into_array
def calc_subclonal_CNt(
    depth_ratio, 
    clonal_CNt,
    CNn, 

    cellularity, 
    tumor_ploidy, 
    normal_ploidy,
    tumor_avg_depth=1,
    normal_avg_depth=1,
    max_CNt_change=5,
    only_max_ccf=True,
):
    # get delta-CNt (subclonal CNt - clonal CNt) values
    clonal_onecp_ddr = delta_depth_ratio(cellularity, tumor_ploidy, CNn, normal_ploidy)
    baseline_depthratio = theoretical_depth_ratio(
        clonal_CNt, cellularity, tumor_ploidy, CNn, normal_ploidy, 
    )
    ddr_ratios = np.atleast_1d(
        (depth_ratio - baseline_depthratio) / clonal_onecp_ddr
    )
    initial_deltaCNts = np.ceil(np.abs(ddr_ratios)).astype(int)
    initial_deltaCNts[initial_deltaCNts == 0] = 1

    if only_max_ccf:
        delta_CNt = initial_deltaCNts
    else:
        delta_CNt = initial_deltaCNts[np.newaxis, :] + np.arange(max_CNt_change)[:, np.newaxis]

    delta_CNt = np.copysign(delta_CNt, ddr_ratios)

    # get ccf from delta_CNt
    ccf = get_ccf_from_depthratio(
        depth_ratio=depth_ratio,
        cellularity=cellularity,
        tumor_ploidy=tumor_ploidy,
        normal_ploidy=normal_ploidy,
        clonal_CNt=clonal_CNt,
        delta_CNt=delta_CNt,
        CNn=CNn,
        tumor_avg_depth=tumor_avg_depth,
        normal_avg_depth=normal_avg_depth,
    )
    subclonal_CNt = clonal_CNt + delta_CNt

    return subclonal_CNt, ccf


def calc_subclonal_B(
    cellularity, 
    ccf, 
    clonal_CNt, 
    subclonal_CNt, 
    CNn,
    baf,
    clonal_B,
    Bn,
):
    '''
    <Derivation>
    baf = baf_numerator_subclone / depthratio_tumor_numerator_subclone
    baf_numerator_subclone = baf * depthratio_tumor_numerator_subclone

    baf_numerator_subclone = (
        (clonal_B * cellularity * (1 - ccf))
        + (subclonal_B * cellularity * ccf)
        + Bn * (1 - cellularity)
    )
    (subclonal_B * cellularity * ccf) = (
        (baf * depthratio_tumor_numerator_subclone)
        - (clonal_B * cellularity * (1 - ccf))
        - Bn * (1 - cellularity)
    )
    (subclonal_B * cellularity * ccf) = (
        (baf * depthratio_tumor_numerator_subclone)
        - (cellularity * (clonal_B * (1 - ccf) - Bn))
        - Bn
    )
    subclonal_B = (
        (baf * depthratio_tumor_numerator_subclone)
        - (cellularity * (clonal_B * (1 - ccf) - Bn))
        - Bn
    ) / (cellularity * ccf)
    '''
    dtns = depthratio_tumor_numerator_subclone(
        cellularity, ccf, clonal_CNt, subclonal_CNt, CNn,
    )
    return (
        (baf * dtns)
        - (cellularity * (clonal_B * (1 - ccf) - Bn))
        - Bn
    ) / (cellularity * ccf)


#def get_subclonal_B(
#    cellularity, 
#    ccf, 
#    clonal_CNt, 
#    subclonal_CNt, 
#    CNn,
#    baf,
#    clonal_B,
#):
#    baf_gthalf = 1 - baf
#    subclonal_B_lthalf = calc_subclonal_B(
#        cellularity=cellularity, 
#        ccf=ccf, 
#        clonal_CNt=clonal_CNt, 
#        subclonal_CNt=subclonal_CNt, 
#        CNn=CNn,
#        baf=baf,
#        clonal_B=clonal_B,
#    )
#    subclonal_B_gthalf = calc_subclonal_B(
#        cellularity=cellularity, 
#        ccf=ccf, 
#        clonal_CNt=clonal_CNt, 
#        subclonal_CNt=subclonal_CNt, 
#        CNn=CNn,
#        baf=baf_gthalf,
#        clonal_B=clonal_B,
#    )
#
#    return subclonal_B_lthalf, subclonal_B_gthalf


def get_average_B(baf, average_CNt, CNn, Bn, cellularity):
    '''
    baf numerator 
        = (
            (clonal_B * cellularity * (1 - ccf))
            + (subclonal_B * cellularity * ccf)
            + Bn * (1 - cellularity)
        )
        = (
            cellularity * average_B
            + Bn * (1 - cellularity)
        )
    baf denominator 
        = (
            (clonal_CNt * cellularity * (1 - ccf)) 
            + (subclonal_CNt * cellularity * ccf) 
            + (CNn * (1 - cellularity))
        ) 
        = (
            cellularity * average_CNt
            + (CNn * (1 - cellularity))
        )
        = cellularity * (average_CNt - CNn) + CNn

    baf numerator = baf * baf denominator
    (
        cellularity * average_B 
        + Bn * (1 - cellularity)
    ) = baf * (cellularity * (average_CNt - CNn) + CNn)
    cellularity * average_B = (
        baf * (cellularity * (average_CNt - CNn) + CNn)
        - Bn * (1 - cellularity)
    )
    average_B = (
        baf * (cellularity * (average_CNt - CNn) + CNn)
        - Bn * (1 - cellularity)
    ) / cellularity
    '''
    return (
        baf * (cellularity * (average_CNt - CNn) + CNn)
        - Bn * (1 - cellularity)
    ) / cellularity


def get_CNt_candidates(average_CNt, min_N_CNt_candidates, N_CNt_candidates_fraction):
    average_CNt = np.squeeze(average_CNt)
    assert average_CNt.ndim == 1

    # make lower_CNts, upper_CNts
    upperCNt_border, lowerCNt_border = get_nearest_integers(average_CNt)
    num_CNt_candidate = np.maximum(
        min_N_CNt_candidates, 
        np.rint(average_CNt * N_CNt_candidates_fraction),
    )
    max_dimsize = np.nanmax(num_CNt_candidate, axis=None)
    CNt_shifts = np.arange(max_dimsize)
    lower_CNts = lowerCNt_border[:, np.newaxis] - CNt_shifts
    upper_CNts = upperCNt_border[:, np.newaxis] + CNt_shifts

    # masking irrelevant values
    upper_mask = np.stack(
        [
            (
                np.repeat(True, max_dimsize)
                if np.isnan(x) else
                np.repeat([False, True], (x, max_dimsize - x)) 
            )
            for x in num_CNt_candidate
        ], 
        axis=0,
    )
    lower_mask = np.logical_or(upper_mask, lower_CNts < 0)

    lower_CNts = np.where(lower_mask, np.nan, lower_CNts)
    upper_CNts = np.where(upper_mask, np.nan, upper_CNts)

    return lower_CNts, upper_CNts


def get_min_diff_args(diffs, global_mask=None):
    """raveled argmin is 0 for global_masked positions"""
    if global_mask is None:
        global_mask = np.full(diffs.shape[0], False)

    diffs_nomask = diffs[~global_mask]
    minargs = np.zeros(diffs.shape[0], dtype=int)
    minargs[~global_mask] = np.nanargmin(
        diffs_nomask.reshape(diffs_nomask.shape[0], -1), 
        axis=1,
    )
    minargs_unraveled = (
        (np.arange(diffs.shape[0]),)
        + np.unravel_index(minargs, shape=diffs.shape[1:])
    )
    return minargs_unraveled


def get_min_diff_args_debug(diffs, global_mask):
    """raveled argmin is 0 for global_masked positions"""
    all_nan_indexes = np.isnan(diffs).reshape(diffs.shape[0], -1).all(axis=1)
    nonglobal_all_nan_indexes = np.logical_and(all_nan_indexes, ~global_mask)

    diffs_nomask = diffs[~all_nan_indexes]
    minargs = np.zeros(diffs.shape[0], dtype=int)
    minargs[~all_nan_indexes] = np.nanargmin(
        diffs_nomask.reshape(diffs_nomask.shape[0], -1), 
        axis=1,
    )
    minargs_unraveled = (
        (np.arange(diffs.shape[0]),)
        + np.unravel_index(minargs, shape=diffs.shape[1:])
    )
    return minargs_unraveled, nonglobal_all_nan_indexes


def average_CNt_arghandler(
    average_CNt, depth_ratio, CNn, cellularity, tumor_ploidy, normal_ploidy,
):
    if average_CNt is None:
        average_CNt = inverse_theoretical_depth_ratio(
            depth_ratio=depth_ratio, 
            CNn=CNn, 
            cellularity=cellularity, 
            tumor_ploidy=tumor_ploidy, 
            normal_ploidy=normal_ploidy,
        )
    else:
        average_CNt = np.asarray(average_CNt)
        assert average_CNt.shape == depth_ratio.shape

    return average_CNt


@deco.get_deco_broadcast(['depth_ratio', 'CNn'])
def find_subclonal_solution_without_B(
    depth_ratio,
    CNn,
    cellularity,
    tumor_ploidy,
    normal_ploidy,
    ccf_candidates,
    average_CNt=None,
    min_N_CNt_candidates=5,
    N_CNt_candidates_fraction=0.5,
):
    """For positions where B-allele frequency doesn't make sense because CNn < 2
    """
    # modify args
    ccf_candidates = np.asarray(ccf_candidates)

    # get average CNt
    average_CNt = average_CNt_arghandler(
        average_CNt, depth_ratio, CNn, cellularity, tumor_ploidy, normal_ploidy,
    )

    # find (average_CNt < 0) regions and treat differently
    lt0_CNts = (average_CNt < 0)  # lt0: less than 0
    valid_CNts = np.logical_not(lt0_CNts)

    # make valid position subsets of input data
    average_CNt_valid = average_CNt[valid_CNts]
    depth_ratio_valid = depth_ratio[valid_CNts]
    CNn_valid = CNn[valid_CNts]

    # make lower_CNts, upper_CNts
    lower_CNt_cand, upper_CNt_cand = get_CNt_candidates(
        average_CNt=average_CNt_valid, 
        min_N_CNt_candidates=min_N_CNt_candidates, 
        N_CNt_candidates_fraction=N_CNt_candidates_fraction,
    )
        # ndim == 2

    # make ccfs
    calculated_ccf = get_ccf_from_depthratio(
        depth_ratio=np.expand_dims(depth_ratio_valid, axis=(1, 2)),
        cellularity=cellularity,
        tumor_ploidy=tumor_ploidy,
        normal_ploidy=normal_ploidy,
        clonal_CNt=np.expand_dims(lower_CNt_cand, axis=2),
        subclonal_CNt=np.expand_dims(upper_CNt_cand, axis=1),
        CNn=np.expand_dims(CNn_valid, axis=(1, 2)),
    )
        # ndim == 3 

    # pick a minimal difference ccf candidate for each position
    ccfs_invertpair = np.stack((calculated_ccf, 1 - calculated_ccf), axis=3)  # ndim == 4
    ccf_diffs = np.abs(np.expand_dims(ccfs_invertpair, axis=4) - ccf_candidates)  # ndim == 5
    minargs = get_min_diff_args(ccf_diffs)
        # axis 0: positions
        # axis 1: lower CNt
        # axis 2: upper CNt
        # axis 3: original ccf, inverted ccf
        # axis 4: ccf candidates

    # check whether original or inverted ccf was chosen
    ccf_inverted = (minargs[3] == 1)

    # make result CNts
    chosen_lower_CNt = lower_CNt_cand[(minargs[0], minargs[1])]
    chosen_upper_CNt = upper_CNt_cand[(minargs[0], minargs[2])]

    result_clonal_CNt = np.zeros(ccf_diffs.shape[0], dtype=float)
    result_clonal_CNt[ccf_inverted] = chosen_upper_CNt[ccf_inverted]
    result_clonal_CNt[~ccf_inverted] = chosen_lower_CNt[~ccf_inverted]

    result_subclonal_CNt = np.zeros(ccf_diffs.shape[0], dtype=float)
    result_subclonal_CNt[ccf_inverted] = chosen_lower_CNt[ccf_inverted]
    result_subclonal_CNt[~ccf_inverted] = chosen_upper_CNt[~ccf_inverted]

    # final result
    result = dict()
    result['average_CNt'] = average_CNt

    result['clonal_CNt'] = result_clonal_CNt
    result['subclonal_CNt'] = result_subclonal_CNt
    result['ccf_chosen'] = ccf_candidates[minargs[4:]]
    result['ccf_calculated'] = ccfs_invertpair[minargs[:4]]
    result['ccf_diffs'] = ccf_diffs
    result['ccf_inverted'] = ccf_inverted

    return result


@deco.get_deco_broadcast(['depth_ratio', 'baf', 'CNn', 'Bn'])
def find_subclonal_solution_B_only(
    depth_ratio,
    baf,
    CNn,
    cellularity,
    tumor_ploidy,
    normal_ploidy,
    ccf_candidates,
    Bn=1,
    average_CNt=None,
):
    # modify args
    ccf_candidates = np.asarray(ccf_candidates)

    # average CNt
    average_CNt = average_CNt_arghandler(
        average_CNt, depth_ratio, CNn, cellularity, tumor_ploidy, normal_ploidy,
    )

    # determine CNt solutions (do not test over clonal/subclonal CNt combinations)
    average_CNt_rint = np.clip(np.rint(average_CNt), 0, None)
    clonal_CNt = average_CNt_rint.copy()  # ndim == 1
    subclonal_CNt = average_CNt_rint.copy()  # ndim == 1

    # get B candidates
    max_CNt = clonal_CNt.max(axis=None)
    clonal_B_cand = np.expand_dims(clonal_CNt, axis=1) - np.arange(max_CNt + 1)
    clonal_B_cand[clonal_B_cand < 0] = np.nan  # ndim == 2
    subclonal_B_cand = clonal_B_cand.copy()  # ndim == 2

    # calculated ccf
    with warnings.catch_warnings(): 
        warnings.simplefilter('ignore', category=RuntimeWarning)
            # divide by zero occurs when clonal_B == subclonal_B
        calculated_ccf = get_ccf_from_baf(
            baf=np.expand_dims(baf, axis=(1, 2)),
            cellularity=cellularity,
            clonal_CNt=np.expand_dims(clonal_CNt, axis=(1, 2)),
            subclonal_CNt=np.expand_dims(subclonal_CNt, axis=(1, 2)),
            clonal_B=np.expand_dims(clonal_B_cand, axis=2),
            subclonal_B=np.expand_dims(subclonal_B_cand, axis=1),
            CNn=np.expand_dims(CNn, axis=(1, 2)),
            Bn=np.expand_dims(Bn, axis=(1, 2)),
        )
        # ndim == 3 
        # axis 0: positions & CNts, axis 1: clonal B, axis 2: subclonal B
    #calculated_ccf = np.clip(calculated_ccf, 0, 1)

    # pick a minimal difference ccf candidate for each position
    ccf_diffs = np.abs(np.expand_dims(calculated_ccf, axis=3) - ccf_candidates) 
        # ndim == 4
        # axis 0: positions, axis 1: clonal B, axis 2: subclonal B, axis 3: ccf candidates
    minargs = get_min_diff_args(ccf_diffs)

    # result
    result = dict()
    result['clonal_CNt'] = clonal_CNt
    result['subclonal_CNt'] = subclonal_CNt
    result['clonal_B'] = clonal_B_cand[(minargs[0], minargs[1])]
    result['subclonal_B'] = subclonal_B_cand[(minargs[0], minargs[2])]
    result['ccf_chosen'] = ccf_candidates[minargs[3:]]

    return result


@deco.get_deco_broadcast(['depth_ratio', 'baf', 'CNn', 'Bn'])
def find_subclonal_solution_freeccf(
    depth_ratio,
    baf,
    CNn,
    cellularity,
    tumor_ploidy,
    normal_ploidy,
    Bn=1,
    average_CNt=None,
    min_N_CNt_candidates=5,
    N_CNt_candidates_fraction=0.5,
):
    # average_CNt = clonal_CNt * (1 - ccf) + subclonal_CNt * ccf
    average_CNt = average_CNt_arghandler(
        average_CNt, depth_ratio, CNn, cellularity, tumor_ploidy, normal_ploidy,
    )

    # find (average_CNt < 0) regions and treat differently
    lt0_CNts = (average_CNt < 0)  # lt0: less than 0
    valid_CNts = np.logical_not(lt0_CNts)

    # make valid position subsets of input data
    average_CNt_valid = average_CNt[valid_CNts]
    depth_ratio_valid = depth_ratio[valid_CNts]
    baf_valid = baf[valid_CNts]
    CNn_valid = CNn[valid_CNts]
    Bn_valid = Bn[valid_CNts]

    # make lower_CNts, upper_CNts
    lower_CNt_cand, upper_CNt_cand = get_CNt_candidates(
        average_CNt=average_CNt_valid, 
        min_N_CNt_candidates=min_N_CNt_candidates, 
        N_CNt_candidates_fraction=N_CNt_candidates_fraction,
    )
        # ndim == 2 for lower_CNt and upper_CNt
        # *_cand: candidate

    # make ccfs
    ccf_cand = get_ccf_from_depthratio(
        depth_ratio=np.expand_dims(depth_ratio_valid, axis=(1, 2)),
        cellularity=cellularity,
        tumor_ploidy=tumor_ploidy,
        normal_ploidy=normal_ploidy,
        clonal_CNt=np.expand_dims(lower_CNt_cand, axis=2),
        subclonal_CNt=np.expand_dims(upper_CNt_cand, axis=1),
        CNn=np.expand_dims(CNn_valid, axis=(1, 2)),
    )
        # ndim == 3 

    # make B
    lower_B_cand = (
        np.expand_dims(lower_CNt_cand, axis=2) 
        - np.arange(np.nanmax(lower_CNt_cand, axis=None) + 1)
    )
        # ndim == 3 (axis 0: positions, axis 1: lower CNt, axis 2: lower B)
    lower_B_cand[lower_B_cand < 0] = np.nan

    upper_B_cand = calc_subclonal_B(
        cellularity=cellularity,
        ccf=np.expand_dims(ccf_cand, axis=3),
        clonal_CNt=np.expand_dims(lower_CNt_cand, axis=(2, 3)),
        subclonal_CNt=np.expand_dims(upper_CNt_cand, axis=(1, 3)),
        CNn=np.expand_dims(CNn_valid, axis=(1, 2, 3)),
        baf=np.expand_dims(baf_valid, axis=(1, 2, 3)),
        clonal_B=np.expand_dims(lower_B_cand, axis=2),
        Bn=np.expand_dims(Bn_valid, axis=(1, 2, 3)),
    )
        # ndim == 4 
        # (axis 0: positions, axis 1: lower CNt, axis 2: upper CNt, axis 3: upper & lower B)
        # dtype float
    upper_B_cand = np.clip(
        np.rint(upper_B_cand), 
        0, 
        np.expand_dims(upper_CNt_cand, axis=(1, 3)),
    )

    # make upper B mask
#    upper_B_mask = np.logical_or(
#        upper_B_rint > np.expand_dims(upper_CNt, axis=(1, 3)),
#        upper_B_rint < 0,
#    )
#    upper_B_rint[upper_B_mask] = np.nan

    # pick optimal solutions
    #B_diffs = np.abs(upper_B - upper_B_rint)
    expected_bafs = theoretical_baf_subclone(
        clonal_CNt=np.expand_dims(lower_CNt_cand, axis=(2, 3)), 
        clonal_B=np.expand_dims(lower_B_cand, axis=2), 
        subclonal_CNt=np.expand_dims(upper_CNt_cand, axis=(1, 3)), 
        subclonal_B=upper_B_cand, 
        ccf=np.expand_dims(ccf_cand, axis=3),
        cellularity=cellularity, 
        CNn=np.expand_dims(CNn_valid, axis=(1, 2, 3)),
        Bn=np.expand_dims(Bn_valid, axis=(1, 2, 3)),
        rescue_nan_ccf=False,
    )
    baf_diffs = np.abs(expected_bafs - np.expand_dims(baf_valid, axis=(1, 2, 3)))
    minargs = get_min_diff_args(baf_diffs)
        # axis 0: positions
        # axis 1: lower CNt
        # axis 2: upper CNt
        # axis 3: lower B (& upper B)

    # result
        # setup
    result = dict()
    result['average_CNt'] = average_CNt
    result['lt0_CNt_flag'] = lt0_CNts
    result['valid_CNt_flag'] = valid_CNts

    result['lower_CNt'] = np.empty(len(average_CNt), dtype=float)
    result['upper_CNt'] = np.empty(len(average_CNt), dtype=float)
    result['lower_B'] = np.empty(len(average_CNt), dtype=float)
    result['upper_B'] = np.empty(len(average_CNt), dtype=float)
    result['ccf'] = np.empty(len(average_CNt), dtype=float)

        # set values for lt0_CNts position
    result['lower_CNt'][lt0_CNts] = 0
    result['upper_CNt'][lt0_CNts] = 0
    result['ccf'][lt0_CNts] = 0
    result['lower_B'][lt0_CNts] = 0
    result['upper_B'][lt0_CNts] = 0

        # set values for valid CNts position
    result['lower_CNt'][valid_CNts] = lower_CNt_cand[(minargs[0], minargs[1])]
    result['upper_CNt'][valid_CNts] = upper_CNt_cand[(minargs[0], minargs[2])]
    result['ccf'][valid_CNts] = ccf_cand[(minargs[0], minargs[1], minargs[2])]
    result['lower_B'][valid_CNts] = lower_B_cand[(minargs[0], minargs[1], minargs[3])]
    result['upper_B'][valid_CNts] = upper_B_cand[minargs]

    return result


def upper_from_lower_CNt(average_CNt, lower_CNt, ccf):
    '''
    <Derivation>
    lower_CNt * (1 - ccf) + upper_CNt * ccf = average_CNt
    upper_CNt * ccf = average_CNt - lower_CNt * (1 - ccf)
    upper_CNt = (average_CNt - lower_CNt * (1 - ccf)) / ccf
    upper_CNt = (average_CNt - lower_CNt + lower_CNt * ccf) / ccf
    upper_CNt = (average_CNt - lower_CNt) / ccf + lower_CNt
    '''
    return (average_CNt - lower_CNt) / ccf + lower_CNt


@deco.get_deco_broadcast(['depth_ratio', 'baf', 'CNn', 'Bn'])
def find_subclonal_solution_fixedccf(
    depth_ratio,
    baf,
    CNn,
    cellularity,
    tumor_ploidy,
    normal_ploidy,
    fixed_ccfs,
    depthratio_diff_factor,
    baf_diff_factor,

    Bn=1,
    average_CNt=None,
    min_N_CNt_candidates=5,
    N_CNt_candidates_fraction=0.5,
):
    # modify args
    fixed_ccfs = np.asarray(fixed_ccfs)

    # average_CNt = clonal_CNt * (1 - ccf) + subclonal_CNt * ccf
    average_CNt = average_CNt_arghandler(
        average_CNt, depth_ratio, CNn, cellularity, tumor_ploidy, normal_ploidy,
    )

    # find (average_CNt < 0) regions and treat differently
    lt0_CNts = (average_CNt < 0)  # lt0: less than 0
    valid_CNts = np.logical_not(lt0_CNts)

    # make valid position subsets of input data
    average_CNt_valid = average_CNt[valid_CNts]
    depth_ratio_valid = depth_ratio[valid_CNts]
    baf_valid = baf[valid_CNts]
    CNn_valid = CNn[valid_CNts]
    Bn_valid = Bn[valid_CNts]

    # expand dims
        # axis 0: positions, axis 1: lower & upper CNt, axis 2: lower & upper B, axis 3: ccfs
    average_CNt_valid = np.expand_dims(average_CNt_valid, axis=(1, 2, 3))
    depth_ratio_valid = np.expand_dims(depth_ratio_valid, axis=(1, 2, 3))
    baf_valid = np.expand_dims(baf_valid, axis=(1, 2, 3))
    CNn_valid = np.expand_dims(CNn_valid, axis=(1, 2, 3))
    Bn_valid = np.expand_dims(Bn_valid, axis=(1, 2, 3))
    fixed_ccfs = np.expand_dims(fixed_ccfs, axis=(0, 1, 2))

    # make clonal CNts
    lower_CNt_cand, upper_CNt_cand = get_CNt_candidates(
        average_CNt=average_CNt_valid, 
        min_N_CNt_candidates=min_N_CNt_candidates, 
        N_CNt_candidates_fraction=N_CNt_candidates_fraction,
    )
    clonal_CNt_cand = np.concatenate([lower_CNt_cand, upper_CNt_cand], axis=1)
    clonal_CNt_cand = np.expand_dims(clonal_CNt_cand, axis=(2, 3))

    # calculate upper CNts using fixed ccfs
    subclonal_CNt_cand = upper_from_lower_CNt(
        average_CNt=average_CNt_valid,
        lower_CNt=clonal_CNt_cand,
        ccf=fixed_ccfs,
    )
    subclonal_CNt_cand = np.rint(subclonal_CNt_cand)

    # make B
    clonal_B_cand = (
        clonal_CNt_cand 
        - np.expand_dims(
            np.arange(np.nanmax(clonal_CNt_cand, axis=None) + 1),
            axis=(0, 1, 3),
        )
    )
    clonal_B_cand[clonal_B_cand < 0] = np.nan

    subclonal_B_cand = calc_subclonal_B(
        cellularity=cellularity,
        ccf=fixed_ccfs,
        clonal_CNt=clonal_CNt_cand,
        subclonal_CNt=subclonal_CNt_cand,
        CNn=CNn_valid,
        baf=baf_valid,
        clonal_B=clonal_B_cand,
        Bn=Bn_valid,
    )
    subclonal_B_cand = np.clip(np.rint(subclonal_B_cand), 0, subclonal_CNt_cand)

    # pick optimal solutions
    expected_depthratios = theoretical_depth_ratio_subclone(
        clonal_CNt=clonal_CNt_cand, 
        subclonal_CNt=subclonal_CNt_cand,
        ccf=fixed_ccfs,
        cellularity=cellularity, 
        tumor_ploidy=tumor_ploidy, 
        CNn=CNn_valid, 
        normal_ploidy=normal_ploidy, 
        tumor_avg_depth=1, 
        normal_avg_depth=1,
        rescue_nan_ccf=False,
    )
    depthratio_diffs = np.abs(expected_depthratios - depth_ratio_valid)

    expected_bafs = theoretical_baf_subclone(
        clonal_CNt=clonal_CNt_cand,
        clonal_B=clonal_B_cand,
        subclonal_CNt=subclonal_CNt_cand,
        subclonal_B=subclonal_B_cand, 
        ccf=fixed_ccfs,
        cellularity=cellularity, 
        CNn=CNn_valid,
        Bn=Bn_valid,
        rescue_nan_ccf=False,
    )
    baf_diffs = np.abs(expected_bafs - baf_valid)
    #minargs = get_min_diff_args(baf_diffs)

    normalized_depthratio_diff_factor = depthratio_diff_factor / (CNn_valid / 2)
    diffsum = (
        depthratio_diffs / normalized_depthratio_diff_factor
        + baf_diffs / baf_diff_factor 
    )
    minargs = get_min_diff_args(diffsum)

    # result
        # setup
    result = dict()
    result['average_CNt'] = average_CNt
    result['clonal_CNt'] = np.empty_like(average_CNt, dtype=float)
    result['subclonal_CNt'] = np.empty_like(average_CNt, dtype=float)
    result['clonal_B'] = np.empty_like(average_CNt, dtype=float)
    result['subclonal_B'] = np.empty_like(average_CNt, dtype=float)
    result['ccf'] = np.empty_like(average_CNt, dtype=float)

        # set values for lt0_CNts position
    result['clonal_CNt'][lt0_CNts] = 0
    result['subclonal_CNt'][lt0_CNts] = 0
    result['ccf'][lt0_CNts] = 0
    result['clonal_B'][lt0_CNts] = 0
    result['subclonal_B'][lt0_CNts] = 0

        # set values for valid CNts position
    result['clonal_CNt'][valid_CNts] = clonal_CNt_cand[minargs[0], minargs[1], 0, 0]
    result['subclonal_CNt'][valid_CNts] = subclonal_CNt_cand[minargs[0], minargs[1], 0, 0]
    result['ccf'][valid_CNts] = fixed_ccfs[0, 0, 0, minargs[3]]
    result['clonal_B'][valid_CNts] = clonal_B_cand[minargs[0], minargs[1], minargs[2], 0]
    result['subclonal_B'][valid_CNts] = subclonal_B_cand[minargs]

    return result


#@deco.get_deco_broadcast(['depth_ratio', 'baf', 'CNn', 'Bn'])
#def find_solution_freeccf(
#    depth_ratio,
#    baf,
#    CNn,
#    cellularity,
#    tumor_ploidy,
#    normal_ploidy,
#    Bn=1,
#    depth_ratio_diff=0.05,
#    baf_diff=0.01,
#    min_N_CNt_candidates=5,
#    N_CNt_candidates_fraction=0.5,
#):
#    # average CNt
#    average_CNt = inverse_theoretical_depth_ratio(
#        depth_ratio=depth_ratio, 
#        cellularity=cellularity, 
#        tumor_ploidy=tumor_ploidy, 
#        CNn=CNn, 
#        normal_ploidy=normal_ploidy, 
#        tumor_avg_depth=1, 
#        normal_avg_depth=1,
#    )
#
#    # clonal solution
#    clonal_solution = find_clonal_solution(
#        depth_ratio=depth_ratio,
#        baf=baf,
#        CNn=CNn,
#        cellularity=cellularity,
#        tumor_ploidy=tumor_ploidy,
#        normal_ploidy=normal_ploidy,
#        Bn=Bn,
#        average_CNt=average_CNt,
#    )
#
#    # get positions to find subclonal solution
#    monoploid_flag = (CNn < 2)
#    depthratio_fit_flag = (np.abs(clonal_solution['depth_ratio_diff']) < depth_ratio_diff)
#    baf_fit_flag = (np.abs(clonal_solution['baf_diff']) < baf_diff)
#    depth_baf_fit_flag = np.logical_and(depthratio_fit_flag, baf_fit_flag)
#
#    monoploid_fit_flag = np.logical_and(monoploid_flag, depthratio_fit_flag)
#    monoploid_unfit_flag = np.logical_and(monoploid_flag, ~depthratio_fit_flag)
#    polyploid_fit_flag = np.logical_and(~monoploid_flag, depth_baf_fit_flag)
#    #polyploid_unfit_flag = np.logical_and(~monoploid_flag, ~depth_baf_fit_flag)
#
#    polyploid_unfit_flag = np.logical_and(~monoploid_flag, ~depthratio_fit_flag)
#    polyploid_unfit_bafonly_flag = functools.reduce(
#        np.logical_and,
#        (
#            ~monoploid_flag,
#            depthratio_fit_flag,
#            ~baf_fit_flag,
#        ),
#    )
#    
#    fit_flag = np.logical_or(monoploid_fit_flag, polyploid_fit_flag)
#
#    # get subclonal solution with B
#    subclonal_solution_withB = find_subclonal_solution_freeccf(
#        depth_ratio=depth_ratio[polyploid_unfit_flag],
#        baf=baf[polyploid_unfit_flag],
#        CNn=CNn[polyploid_unfit_flag],
#        cellularity=cellularity,
#        tumor_ploidy=tumor_ploidy,
#        normal_ploidy=normal_ploidy,
#        Bn=Bn[polyploid_unfit_flag],
#        average_CNt=average_CNt[polyploid_unfit_flag],
#        min_N_CNt_candidates=min_N_CNt_candidates,
#        N_CNt_candidates_fraction=N_CNt_candidates_fraction,
#    )
#
#    # make clonal/subclonal results
#    clonal_CNt = np.full(len(depth_ratio), np.nan)
#    clonal_CNt[fit_flag] = clonal_solution['CNt'][fit_flag]
#    #clonal_CNt[monoploid_unfit_flag] = subclonal_solution_woB['clonal_CNt']
#    #clonal_CNt[polyploid_unfit_bafonly_flag] = subclonal_solution_Bonly['clonal_CNt']
#    clonal_CNt[polyploid_unfit_flag] = subclonal_solution_withB['lower_CNt']
#
#    subclonal_CNt = np.full(len(depth_ratio), np.nan)
#    subclonal_CNt[fit_flag] = np.nan
#    #subclonal_CNt[monoploid_unfit_flag] = subclonal_solution_woB['subclonal_CNt']
#    #subclonal_CNt[polyploid_unfit_bafonly_flag] = subclonal_solution_Bonly['subclonal_CNt']
#    subclonal_CNt[polyploid_unfit_flag] = subclonal_solution_withB['upper_CNt']
#
#    clonal_B = np.full(len(depth_ratio), np.nan)
#    clonal_B[monoploid_flag] = np.nan
#    clonal_B[polyploid_fit_flag] = clonal_solution['B'][polyploid_fit_flag]
#    #clonal_B[polyploid_unfit_bafonly_flag] = subclonal_solution_Bonly['clonal_B']
#    clonal_B[polyploid_unfit_flag] = subclonal_solution_withB['lower_B']
#
#    subclonal_B = np.full(len(depth_ratio), np.nan)
#    subclonal_B[np.logical_or(monoploid_flag, polyploid_fit_flag)] = np.nan
#    #subclonal_B[polyploid_unfit_bafonly_flag] = subclonal_solution_Bonly['subclonal_B']
#    subclonal_B[polyploid_unfit_flag] = subclonal_solution_withB['upper_B']
#
#    ccf = np.full(len(depth_ratio), np.nan)
#    ccf[fit_flag] = np.nan
#    #ccf[monoploid_unfit_flag] = subclonal_solution_woB['ccf_chosen']
#    #ccf[polyploid_unfit_bafonly_flag] = subclonal_solution_Bonly['ccf_chosen']
#    ccf[polyploid_unfit_flag] = subclonal_solution_withB['ccf']
#
#    # make theoretical values and diffs
#    calculated_depth_ratio = theoretical_depth_ratio_subclone(
#        clonal_CNt=clonal_CNt, 
#        subclonal_CNt=subclonal_CNt,
#        ccf=ccf,
#        cellularity=cellularity, 
#        tumor_ploidy=tumor_ploidy, 
#        CNn=CNn, 
#        normal_ploidy=normal_ploidy, 
#        tumor_avg_depth=1, 
#        normal_avg_depth=1,
#        rescue_nan_ccf=True,
#    )
#    calculated_baf = theoretical_baf_subclone(
#        clonal_CNt=clonal_CNt, 
#        clonal_B=clonal_B, 
#        subclonal_CNt=subclonal_CNt, 
#        subclonal_B=subclonal_B, 
#        ccf=ccf,
#        cellularity=cellularity, 
#        CNn=CNn,
#        Bn=Bn,
#        rescue_nan_ccf=True,
#    )
#
#    # result
#    result = dict()
#    result['average_CNt'] = average_CNt
#    result['clonal_CNt'] = clonal_CNt
#    result['subclonal_CNt'] = subclonal_CNt
#    result['clonal_B'] = clonal_B
#    result['subclonal_B'] = subclonal_B
#    result['ccf'] = ccf
#
#    result['clonal_theoretical_depth_ratio'] = clonal_solution['theoretical_depth_ratio']
#    result['clonal_theoretical_baf'] = clonal_solution['theoretical_baf']
#    result['clonal_depth_ratio_diff'] = clonal_solution['depth_ratio_diff']
#    result['clonal_baf_diff'] = clonal_solution['baf_diff']
#
#    result['theoretical_depth_ratio'] = calculated_depth_ratio
#    result['theoretical_baf'] = calculated_baf
#    result['flag'] = {
#        'monoploid': monoploid_flag,
#        'depthratio_fit': depthratio_fit_flag,
#        'baf_fit': baf_fit_flag,
#        'depth_baf_fit': depth_baf_fit_flag,
#
#        'monoploid_fit': monoploid_fit_flag,
#        'monoploid_unfit': monoploid_unfit_flag,
#        'polyploid_fit': polyploid_fit_flag,
#        'polyploid_unfit': polyploid_unfit_flag,
#
#        'fit': fit_flag,
#    }
#
#    return result


def choose_ccfs(peak_ccfs, peak_densities):
    lehalf_index = peak_ccfs <= 0.5  # lehalf: less than or equal to half
    peak_ccfs_lthalf = peak_ccfs[lehalf_index]
    peak_densities_lthalf = peak_densities[lehalf_index]

    # sort by density
    peak_ccfs_lthalf = peak_ccfs_lthalf[np.argsort(peak_densities_lthalf)[::-1]]

    # choose as many values with which cumsum becomes less than 1
    chosen_ccfs = peak_ccfs_lthalf[np.cumsum(peak_ccfs_lthalf) < 1]

    return chosen_ccfs


def select_fixed_ccfs(ccfs, lengths, bandwidth=0.1):
    dup_ccfs = np.concatenate([ccfs, 1 - ccfs])
    dup_lengths = np.tile(lengths, 2)

    peak_values, peak_densities, density = get_density_peaks(
        dup_ccfs, 
        weights=dup_lengths, 
        xs=np.arange(0, 1, 0.01), 
        threshold=None,
        bw_method=bandwidth,
    )

    fixed_ccfs = choose_ccfs(peak_values, peak_densities)
    ccf_plotdata = {
        'ccfs': dup_ccfs,
        'lengths': dup_lengths,
        'density': density,
        'peak_values': peak_values,
    }

    return fixed_ccfs, ccf_plotdata


def get_default_depth_ratio_diff(cellularity, tumor_ploidy, normal_ploidy, CNn):
    assert CNn.ndim == 1

    ddr = delta_depth_ratio(
        cellularity=cellularity, 
        tumor_ploidy=tumor_ploidy, 
        CNn=CNn, 
        normal_ploidy=normal_ploidy, 
        tumor_avg_depth=1, 
        normal_avg_depth=1,
    )
    return 0.1 * ddr


@deco.get_deco_broadcast(['depth_ratio', 'baf', 'CNn', 'Bn', 'lengths'])
def find_solution_before_ccfs(
    depth_ratio,
    baf,
    CNn,
    lengths,
    cellularity,
    tumor_ploidy,
    normal_ploidy,
    depth_ratio_diff,
    baf_diff,
    average_CNt=None,
    Bn=1,
    min_N_CNt_candidates=5,
    N_CNt_candidates_fraction=0.5,
    ccf_bw=0.1,
):
    # average CNt
    average_CNt = average_CNt_arghandler(
        average_CNt, depth_ratio, CNn, cellularity, tumor_ploidy, normal_ploidy,
    )

    # clonal solution
    clonal_solution = find_clonal_solution(
        depth_ratio=depth_ratio,
        baf=baf,
        CNn=CNn,
        cellularity=cellularity,
        tumor_ploidy=tumor_ploidy,
        normal_ploidy=normal_ploidy,
        Bn=Bn,
        average_CNt=average_CNt,
    )

    # get positions to find subclonal solution
    monoploid_flag = (CNn < 2)
    depthratio_fit_flag = (np.abs(clonal_solution['depth_ratio_diff']) < depth_ratio_diff)
    baf_fit_flag = (np.abs(clonal_solution['baf_diff']) < baf_diff)
    depth_baf_fit_flag = np.logical_and(depthratio_fit_flag, baf_fit_flag)

    monoploid_fit_flag = np.logical_and(monoploid_flag, depthratio_fit_flag)
    monoploid_unfit_flag = np.logical_and(monoploid_flag, ~depthratio_fit_flag)
    polyploid_fit_flag = np.logical_and(~monoploid_flag, depth_baf_fit_flag)
    polyploid_unfit_flag = np.logical_and(~monoploid_flag, ~depthratio_fit_flag)
    polyploid_unfit_bafonly_flag = functools.reduce(
        np.logical_and,
        (
            ~monoploid_flag,
            depthratio_fit_flag,
            ~baf_fit_flag,
        ),
    )
    
    fit_flag = np.logical_or(monoploid_fit_flag, polyploid_fit_flag)

    # determine ccf peaks from polyploid depth-baf-unfit regions
    #subclonal_solution_withB = find_subclonal_solution_freeccf(
    freeccf_solution = find_subclonal_solution_freeccf(
        depth_ratio=depth_ratio[polyploid_unfit_flag],
        baf=baf[polyploid_unfit_flag],
        CNn=CNn[polyploid_unfit_flag],
        cellularity=cellularity,
        tumor_ploidy=tumor_ploidy,
        normal_ploidy=normal_ploidy,
        Bn=Bn[polyploid_unfit_flag],
        average_CNt=average_CNt[polyploid_unfit_flag],
        min_N_CNt_candidates=min_N_CNt_candidates,
        N_CNt_candidates_fraction=N_CNt_candidates_fraction,
    )
    valid_CNt_indexes = freeccf_solution['valid_CNt_flag'].nonzero()
    fixed_ccfs, ccf_plotdata = select_fixed_ccfs(
        ccfs=freeccf_solution['ccf'][valid_CNt_indexes], 
        lengths=lengths[polyploid_unfit_flag][valid_CNt_indexes],
        bandwidth=ccf_bw,
    )

    flags = {
        'monoploid': monoploid_flag,
        'depthratio_fit': depthratio_fit_flag,
        'baf_fit': baf_fit_flag,
        'depth_baf_fit': depth_baf_fit_flag,
        'monoploid_fit': monoploid_fit_flag,
        'monoploid_unfit': monoploid_unfit_flag,
        'polyploid_fit': polyploid_fit_flag,
        'polyploid_unfit': polyploid_unfit_flag,
        'polyploid_unfit_bafonly': polyploid_unfit_bafonly_flag,
        'fit': fit_flag,
    }

    # make theoretical values and diffs
    calculated_depth_ratio = theoretical_depth_ratio_subclone(
        clonal_CNt=freeccf_solution['lower_CNt'], 
        subclonal_CNt=freeccf_solution['upper_CNt'],
        ccf=freeccf_solution['ccf'],
        cellularity=cellularity, 
        tumor_ploidy=tumor_ploidy, 
        CNn=CNn[polyploid_unfit_flag], 
        normal_ploidy=normal_ploidy, 
        tumor_avg_depth=1, 
        normal_avg_depth=1,
        rescue_nan_ccf=True,
    )
    calculated_baf = theoretical_baf_subclone(
        clonal_CNt=freeccf_solution['lower_CNt'], 
        clonal_B=freeccf_solution['lower_B'], 
        subclonal_CNt=freeccf_solution['upper_CNt'], 
        subclonal_B=freeccf_solution['upper_B'], 
        ccf=freeccf_solution['ccf'],
        cellularity=cellularity, 
        CNn=CNn[polyploid_unfit_flag],
        Bn=Bn[polyploid_unfit_flag],
        rescue_nan_ccf=True,
    )

    return fixed_ccfs, ccf_plotdata, clonal_solution, flags, freeccf_solution, calculated_depth_ratio, calculated_baf


@deco.get_deco_broadcast(['depth_ratio', 'baf', 'CNn', 'Bn', 'lengths', 'average_CNt'])
def find_solution_after_ccfs(
    depth_ratio,
    baf,
    CNn,
    lengths,
    cellularity,
    tumor_ploidy,
    normal_ploidy,
    average_CNt,

    fixed_ccfs, clonal_solution, flags,

    Bn=1,
    min_N_CNt_candidates=5,
    N_CNt_candidates_fraction=0.5,
):
    # make diff factors
    baf_selector = ~flags['monoploid']
    depth_selector = (CNn == 2)
    depthratio_diff_factor = common.get_diffmean(
        depth_ratio[depth_selector], 
        lengths[depth_selector],
    )  # based on CNn == 2 regions
    baf_diff_factor = common.get_diffmean(
        baf[baf_selector], 
        lengths[baf_selector],
    )

    # main
    if flags['polyploid_unfit'].any():
        subclonal_solution_both = find_subclonal_solution_fixedccf(
            depth_ratio=depth_ratio[flags['polyploid_unfit']],
            baf=baf[flags['polyploid_unfit']],
            CNn=CNn[flags['polyploid_unfit']],
            cellularity=cellularity,
            tumor_ploidy=tumor_ploidy,
            normal_ploidy=normal_ploidy,
            fixed_ccfs=fixed_ccfs,
            depthratio_diff_factor=depthratio_diff_factor,
            baf_diff_factor=baf_diff_factor,
            Bn=Bn[flags['polyploid_unfit']],
            average_CNt=average_CNt[flags['polyploid_unfit']],
            min_N_CNt_candidates=min_N_CNt_candidates,
            N_CNt_candidates_fraction=N_CNt_candidates_fraction,
        )

    # get subclonal solution without B
    if flags['monoploid_unfit'].any():
        subclonal_solution_woB = find_subclonal_solution_without_B(
            depth_ratio=depth_ratio[flags['monoploid_unfit']],
            CNn=CNn[flags['monoploid_unfit']],
            cellularity=cellularity,
            tumor_ploidy=tumor_ploidy,
            normal_ploidy=normal_ploidy,
            ccf_candidates=fixed_ccfs,
            average_CNt=average_CNt[flags['monoploid_unfit']],
            min_N_CNt_candidates=min_N_CNt_candidates,
            N_CNt_candidates_fraction=N_CNt_candidates_fraction,
        )

    # get subclonal solution B-only
    if flags['polyploid_unfit_bafonly'].any():
        subclonal_solution_Bonly = find_subclonal_solution_B_only(
            depth_ratio=depth_ratio[flags['polyploid_unfit_bafonly']],
            baf=baf[flags['polyploid_unfit_bafonly']],
            CNn=CNn[flags['polyploid_unfit_bafonly']],
            cellularity=cellularity,
            tumor_ploidy=tumor_ploidy,
            normal_ploidy=normal_ploidy,
            ccf_candidates=fixed_ccfs,
            Bn=Bn[flags['polyploid_unfit_bafonly']],
            average_CNt=average_CNt[flags['polyploid_unfit_bafonly']],
        )

    # make clonal/subclonal results
    N_position = depth_ratio.shape[0]

    clonal_CNt = np.full(N_position, np.nan)
    subclonal_CNt = np.full(N_position, np.nan)
    clonal_B = np.full(N_position, np.nan)
    subclonal_B = np.full(N_position, np.nan)
    ccf = np.full(N_position, np.nan)

    if flags['fit'].any():
        clonal_CNt[flags['fit']] = clonal_solution['CNt'][flags['fit']]
        subclonal_CNt[flags['fit']] = np.nan
        ccf[flags['fit']] = np.nan

    if flags['monoploid_unfit'].any():
        clonal_CNt[flags['monoploid_unfit']] = subclonal_solution_woB['clonal_CNt']
        subclonal_CNt[flags['monoploid_unfit']] = subclonal_solution_woB['subclonal_CNt']
        ccf[flags['monoploid_unfit']] = subclonal_solution_woB['ccf_chosen']

    if flags['polyploid_unfit_bafonly'].any():
        clonal_CNt[flags['polyploid_unfit_bafonly']] = subclonal_solution_Bonly['clonal_CNt']
        subclonal_CNt[flags['polyploid_unfit_bafonly']] = subclonal_solution_Bonly['subclonal_CNt']
        clonal_B[flags['polyploid_unfit_bafonly']] = subclonal_solution_Bonly['clonal_B']
        subclonal_B[flags['polyploid_unfit_bafonly']] = subclonal_solution_Bonly['subclonal_B']
        ccf[flags['polyploid_unfit_bafonly']] = subclonal_solution_Bonly['ccf_chosen']

    if flags['polyploid_unfit'].any():
        clonal_CNt[flags['polyploid_unfit']] = subclonal_solution_both['clonal_CNt']
        subclonal_CNt[flags['polyploid_unfit']] = subclonal_solution_both['subclonal_CNt']
        clonal_B[flags['polyploid_unfit']] = subclonal_solution_both['clonal_B']
        subclonal_B[flags['polyploid_unfit']] = subclonal_solution_both['subclonal_B']
        ccf[flags['polyploid_unfit']] = subclonal_solution_both['ccf']

    if flags['monoploid'].any():
        clonal_B[flags['monoploid']] = np.nan
        subclonal_B[np.logical_or(flags['monoploid'], flags['polyploid_fit'])] = np.nan

    if flags['polyploid_fit'].any():
        clonal_B[flags['polyploid_fit']] = clonal_solution['B'][flags['polyploid_fit']]
        subclonal_B[np.logical_or(flags['monoploid'], flags['polyploid_fit'])] = np.nan

    subclonal_B[np.logical_or(flags['monoploid'], flags['polyploid_fit'])] = np.nan

    #clonal_CNt = np.full(N_position, np.nan)
    #clonal_CNt[flags['fit']] = clonal_solution['CNt'][flags['fit']]
    #clonal_CNt[flags['monoploid_unfit']] = subclonal_solution_woB['clonal_CNt']
    #clonal_CNt[flags['polyploid_unfit_bafonly']] = subclonal_solution_Bonly['clonal_CNt']
    #clonal_CNt[flags['polyploid_unfit']] = subclonal_solution_both['lower_CNt']

    #subclonal_CNt = np.full(N_position, np.nan)
    #subclonal_CNt[flags['fit']] = np.nan
    #subclonal_CNt[flags['monoploid_unfit']] = subclonal_solution_woB['subclonal_CNt']
    #subclonal_CNt[flags['polyploid_unfit_bafonly']] = subclonal_solution_Bonly['subclonal_CNt']
    #subclonal_CNt[flags['polyploid_unfit']] = subclonal_solution_both['upper_CNt']

    #clonal_B = np.full(N_position, np.nan)
    #clonal_B[flags['monoploid']] = np.nan
    #clonal_B[flags['polyploid_fit']] = clonal_solution['B'][flags['polyploid_fit']]
    #clonal_B[flags['polyploid_unfit_bafonly']] = subclonal_solution_Bonly['clonal_B']
    #clonal_B[flags['polyploid_unfit']] = subclonal_solution_both['lower_B']

    #subclonal_B = np.full(N_position, np.nan)
    #subclonal_B[np.logical_or(flags['monoploid'], flags['polyploid_fit'])] = np.nan
    #subclonal_B[flags['polyploid_unfit_bafonly']] = subclonal_solution_Bonly['subclonal_B']
    #subclonal_B[flags['polyploid_unfit']] = subclonal_solution_both['upper_B']

    #ccf = np.full(N_position, np.nan)
    #ccf[flags['fit']] = np.nan
    #ccf[flags['monoploid_unfit']] = subclonal_solution_woB['ccf_chosen']
    #ccf[flags['polyploid_unfit_bafonly']] = subclonal_solution_Bonly['ccf_chosen']
    #ccf[flags['polyploid_unfit']] = subclonal_solution_both['ccf']

    # make theoretical values and diffs
    calculated_depth_ratio = theoretical_depth_ratio_subclone(
        clonal_CNt=clonal_CNt, 
        subclonal_CNt=subclonal_CNt,
        ccf=ccf,
        cellularity=cellularity, 
        tumor_ploidy=tumor_ploidy, 
        CNn=CNn, 
        normal_ploidy=normal_ploidy, 
        tumor_avg_depth=1, 
        normal_avg_depth=1,
        rescue_nan_ccf=True,
    )
    calculated_baf = theoretical_baf_subclone(
        clonal_CNt=clonal_CNt, 
        clonal_B=clonal_B, 
        subclonal_CNt=subclonal_CNt, 
        subclonal_B=subclonal_B, 
        ccf=ccf,
        cellularity=cellularity, 
        CNn=CNn,
        Bn=Bn,
        rescue_nan_ccf=True,
    )

    # result
    result = dict()
    result['average_CNt'] = average_CNt
    result['clonal_CNt'] = clonal_CNt
    result['subclonal_CNt'] = subclonal_CNt
    result['clonal_B'] = clonal_B
    result['subclonal_B'] = subclonal_B
    result['ccf'] = ccf

    result['clonal_theoretical_depth_ratio'] = clonal_solution['theoretical_depth_ratio']
    result['clonal_theoretical_baf'] = clonal_solution['theoretical_baf']
    result['clonal_depth_ratio_diff'] = clonal_solution['depth_ratio_diff']
    result['clonal_baf_diff'] = clonal_solution['baf_diff']

    result['theoretical_depth_ratio'] = calculated_depth_ratio
    result['theoretical_baf'] = calculated_baf
    result['flag'] = flags

    result['fixed_ccfs'] = fixed_ccfs
    result['clonal_solution'] = clonal_solution

    return result


@deco.get_deco_broadcast(['depth_ratio', 'baf', 'CNn', 'Bn', 'lengths'])
def find_solution(
    depth_ratio,
    baf,
    CNn,
    lengths,
    cellularity,
    tumor_ploidy,
    normal_ploidy,
    depth_ratio_diff=None,
    baf_diff=0.05,
    Bn=1,
    min_N_CNt_candidates=5,
    N_CNt_candidates_fraction=0.5,
):
    # arg handling
    if depth_ratio_diff is None:
        depth_ratio_diff = get_default_depth_ratio_diff(
            cellularity, tumor_ploidy, normal_ploidy, CNn,
        )

    # average CNt
    average_CNt = inverse_theoretical_depth_ratio(
        depth_ratio=depth_ratio, 
        CNn=CNn, 
        cellularity=cellularity, 
        tumor_ploidy=tumor_ploidy, 
        normal_ploidy=normal_ploidy,
    )

    # determine ccf peaks from polyploid depth-baf-unfit regions
    (
        fixed_ccfs, 
        ccf_plotdata, 
        clonal_solution, 
        flags, 
        freeccf_solution,
        calculated_depth_ratio, 
        calculated_baf,
    ) = find_solution_before_ccfs(
        depth_ratio=depth_ratio,
        baf=baf,
        CNn=CNn,
        lengths=lengths,
        cellularity=cellularity,
        tumor_ploidy=tumor_ploidy,
        normal_ploidy=normal_ploidy,
        depth_ratio_diff=depth_ratio_diff,
        baf_diff=baf_diff,
        average_CNt=average_CNt,
        Bn=Bn,
        min_N_CNt_candidates=min_N_CNt_candidates,
        N_CNt_candidates_fraction=N_CNt_candidates_fraction,
    )

    # make final solution with fixed ccfs
    solution = find_solution_after_ccfs(
        depth_ratio=depth_ratio,
        baf=baf,
        CNn=CNn,
        lengths=lengths,
        cellularity=cellularity,
        tumor_ploidy=tumor_ploidy,
        normal_ploidy=normal_ploidy,
        average_CNt=average_CNt,

        fixed_ccfs=fixed_ccfs, 
        clonal_solution=clonal_solution, 
        flags=flags,

        Bn=Bn,
        min_N_CNt_candidates=min_N_CNt_candidates,
        N_CNt_candidates_fraction=N_CNt_candidates_fraction,
    )

    return solution, ccf_plotdata


#########################
# find optimal cp value #
#########################

#def get_depthratio_peaks(depthratio_segment_means, lengths, limits=(0, 2), step=0.01, peak_cutoff=1e8):
#    bins = np.arange(limits[0], limits[1], step)
#    peak_indexes, histbins, hist, bin_midpoints = get_hist_peaks(depthratio_segment_means, weights=lengths, bins=bins)
#    peak_indexes = peak_indexes[hist[peak_indexes] > peak_cutoff]
#    peak_xs = bin_midpoints[peak_indexes]
#
#    return peak_xs


def calc_cp_score(
    segment_df, cellularity, tumor_ploidy, is_female, CNt_weight, normal_ploidy,
):
    def applied_func(row):
        #if check_haploid(is_female, row.iloc[0]) or np.isnan(row.iloc[4]):
            # Male X chromosome PAR should be treated as non-haploid
        if (row.iloc[5] == 1) or np.isnan(row.iloc[4]):
            CNt, B, segfit_score = get_CN_from_cp_wo_baf(
                cellularity, 
                tumor_ploidy, 
                row.iloc[3], 
                CNt_weight, 
                row.iloc[5], 
                normal_ploidy,
            )
        else:
            CNt, B, segfit_score = get_CN_from_cp(
                cellularity, 
                tumor_ploidy, 
                row.iloc[3], 
                row.iloc[4], 
                CNt_weight, 
                row.iloc[5], 
                normal_ploidy,
            )

        return CNt, B, segfit_score

    #assert {'depthratio_segment_mean', 'baf_segment_mean', 'CNn'}.issubset(segment_df.columns)
    assert {'depthratio_segment_mean', 'corrected_baf_segment_mean', 'CNn'}.issubset(segment_df.columns)
    assert segment_df['depthratio_segment_mean'].notna().all()

    # set column order
    segment_df = segment_df.loc[
        :, 
        ['Chromosome', 'Start', 'End', 'depthratio_segment_mean', 'corrected_baf_segment_mean', 'CNn']
    ]

    apply_result = segment_df.apply(applied_func, axis=1)
    segfit_score_sum = sum(x[2] for x in apply_result if not np.isnan(x[2]))
    CNt_list = np.fromiter((x[0] for x in apply_result), dtype=float)
    B_list = np.fromiter((x[1] for x in apply_result), dtype=float)

    return {'segfit_score': segfit_score_sum, 'CNt_list': CNt_list, 'B_list': B_list}


def get_cp_score_dict(
    segment_df, 
    refver, 
    is_female, 
    target_region_gr=None, 
    CNt_weight=DEFAULT_CNT_WEIGHT, 
    normal_ploidy=None,
    nproc=None,
):
    def get_calculated_tumor_ploidies_wgs(segment_df, scorelist):
        weights = segment_df['End'] - segment_df['Start']
        CNt_values = np.array([x['CNt_list'] for x in scorelist])
        calculated_tumor_ploidies = np.average(CNt_values, axis=1, weights=weights)
        return calculated_tumor_ploidies

    def get_calculated_tumor_ploidies_new(segment_df, scorelist, target_region_gr):
        stripped_segment_df = segment_df.loc[:, ['Chromosome', 'Start', 'End']].copy()
        all_indexes = list(range(stripped_segment_df.shape[0]))  # means index of segments
        stripped_segment_df.insert(3, 'index', all_indexes)
        index_annotated_targetregion_gr = annotate_region_with_segment(  # each region is annotated with corresponding segment index
            target_region_gr[['Chromosome', 'Start', 'End']], 
            pr.PyRanges(stripped_segment_df),
        )
        # make weights of CNt values
        index_annotated_targetregion_gr.length = index_annotated_targetregion_gr.lengths()
        weights_dict = index_annotated_targetregion_gr.df.loc[:, ['length', 'index']].groupby('index').sum().to_dict()['length']
        weights = [
            (weights_dict[x] if x in weights_dict else 0)
            for x in all_indexes   
        ]

        CNt_values = np.array([x['CNt_list'] for x in scorelist])
        calculated_tumor_ploidies = np.average(CNt_values, axis=1, weights=weights)

        return calculated_tumor_ploidies

    def get_calculated_tumor_ploidies(segment_df, scorelist, target_region_gr):
        stripped_segment_df = segment_df.loc[:, ['Chromosome', 'Start', 'End']]
        CNt_col_list = [
            pd.Series(x['CNt_list'], name=f'CNt_{idx}') 
            for idx, x in enumerate(scorelist)
        ]
        CNt_annotated_segment_df = pd.concat(
            ([stripped_segment_df] + CNt_col_list), axis=1
        )
        CNt_annotated_segment_gr = pr.PyRanges(CNt_annotated_segment_df)
        CNt_annotated_targetregion_gr = annotate_region_with_segment(
            target_region_gr[['Chromosome', 'Start', 'End']], CNt_annotated_segment_gr
        )
        calculated_tumor_ploidies = np.average(
            CNt_annotated_targetregion_gr.df.iloc[:, 3:],
            axis=0,
            weights=CNt_annotated_targetregion_gr.lengths(),
        )
        return calculated_tumor_ploidies
   
    # df arg handling
    target_region_gr = arg_into_gr(target_region_gr)
    segment_df = arg_into_df(segment_df)

    # get calculated CNt and segment fit cores
    if normal_ploidy is None:
        normal_ploidy = get_normal_mean_ploidy(refver, is_female, target_region_gr)

    c_candidates = np.round(np.arange(0.01, 1.00, 0.01), 2)  # 0 and 1 are excluded
    p_candidates = np.round(np.arange(0.1, 7.1, 0.1), 1)  # 0.1 to 7.0
    cp_pairs = tuple(CPPair(x, y) for (x, y) in itertools.product(c_candidates, p_candidates))

    with multiprocessing.Pool(nproc) as pool:
        scorelist = pool.starmap(
            calc_cp_score, 
            (   
                (segment_df, x.cellularity, x.ploidy, is_female, CNt_weight, normal_ploidy) 
                for x in cp_pairs
            ),
        )

    # get calculated tumor mean ploidy for each cp pair
    if target_region_gr is None:
        calculated_tumor_ploidies = get_calculated_tumor_ploidies_wgs(segment_df, scorelist)
    else:
        calculated_tumor_ploidies = get_calculated_tumor_ploidies_new(
            segment_df, scorelist, target_region_gr
        )

    # calculate ploidy fitting scores
    for cppair, calc_ploidy, dic in zip(cp_pairs, calculated_tumor_ploidies, scorelist):
        dic['calculated_tumor_ploidy'] = calc_ploidy
        dic['ploidy_diff'] = calc_ploidy - cppair.ploidy
        dic['ploidy_diff_ratio'] = dic['ploidy_diff'] / cppair.ploidy
        dic['delta_depth_ratio'] = delta_depth_ratio(
            cellularity=cppair.cellularity, 
            tumor_ploidy=cppair.ploidy, 
            CNn=2, 
            normal_ploidy=normal_ploidy, 
            tumor_avg_depth=1, 
            normal_avg_depth=1,
        )

    cp_scores = dict(zip(cp_pairs, scorelist))
    return cp_scores


def make_cpscore_dfs(cpscore_dict):
    clist = sorted(set(x.cellularity for x in cpscore_dict.keys()))
    plist = sorted(set(x.ploidy for x in cpscore_dict.keys()))
    data = {
        'segfit': list(), 
        'ploidy_diff': list(), 
        'ploidy_diff_abs': list(),
        'ploidy_diff_ratio': list(),
        'ploidy_diff_ratio_abs': list(),
    }
    for c in clist:
        values = [
            val for key, val in 
            sorted(
                ((k, v) for (k, v) in cpscore_dict.items() if k.cellularity == c),
                key=(lambda y: y[0].ploidy),
            )
        ]
        data['segfit'].append([x['segfit_score'] for x in values])
        data['ploidy_diff'].append([x['ploidy_diff'] for x in values])
        data['ploidy_diff_abs'].append([abs(x['ploidy_diff']) for x in values])
        data['ploidy_diff_ratio'].append([x['ploidy_diff_ratio'] for x in values])
        data['ploidy_diff_ratio_abs'].append([abs(x['ploidy_diff_ratio']) for x in values])

    dfs = dict()
    for key, val in data.items():
        df = pd.DataFrame.from_records(val, index=clist, columns=plist)
        df.index.name = 'cellularity'
        df.columns.name = 'ploidy'
        dfs[key] = df

    return dfs


def show_heatmap_peaks(score_df, invert=True, quantile=True, figsize=None):
    #peaks = find_df_peaks(score_df, invert=invert)
    peaks = find_df_peaks_4directions(score_df, invert=invert)
    ys, xs = zip(*peaks)

    fig, ax = plt.subplots(figsize=figsize)
    if quantile:
        ax = sns.heatmap(get_quantile_df(score_df), ax=ax)
    else:
        ax = sns.heatmap(score_df, ax=ax)
    ax.plot(xs, ys, linestyle='', marker='o', markersize=1.5)

    return fig, ax


def show_heatmap_peaks_new(dfs, invert=True, quantile=True, figsize=None, fourway=True):
    if fourway:
        peaks = find_df_peaks_4directions(dfs['segfit'], invert=invert)
    else:
        peaks = find_df_peaks(dfs['segfit'], invert=invert)

    ys, xs = zip(*peaks)
    ploidy_ratios = [dfs['ploidy_diff_ratio'].iloc[y, x] for y, x in zip(ys, xs)]
    segscores = [dfs['segfit'].iloc[y, x] for y, x in zip(ys, xs)]

    fig, ax = plt.subplots(figsize=figsize)
    if quantile:
        ax = sns.heatmap(get_quantile_df(dfs['segfit']), ax=ax)
    else:
        ax = sns.heatmap(dfs['segfit'], ax=ax)
    ax.plot(xs, ys, linestyle='', marker='o', markersize=3)
    for x, y, pratio, segscore in zip(xs, ys, ploidy_ratios, segscores):
        cellularity = dfs['segfit'].index[y]
        ploidy = dfs['segfit'].columns[x]
        text = '\n'.join([
            f'cellularity: {cellularity}',
            f'ploidy: {ploidy}',
            f'segfit score: {round(segscore, 6)}',
            f'ploidy diff ratio: {round(pratio, 6)}',
        ])
        ax.text(
            x + 0.5, y - 0.5, text, 
            fontsize=6,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
        )

    return fig, ax


##################################
# segment dataframe modification #
##################################


def annotate_region_with_segment(target_region_gr, segment_gr, as_gr=True):
    """Add values of segment_gr, such as CNt or B, to target_region_gr"""
    return pyranges_helper.join(
        target_region_gr, segment_gr, how='left', merge='first', 
        find_nearest=True, as_gr=as_gr,
    )


def add_CNt_to_segment(
    segment_df, cellularity, tumor_ploidy, normal_ploidy, is_female, 
    CNt_weight=DEFAULT_CNT_WEIGHT,
):
    assert isinstance(segment_df, pd.DataFrame)

    data = calc_cp_score(
        segment_df, cellularity, tumor_ploidy, is_female, CNt_weight, normal_ploidy,
    )
    result = segment_df.copy()
    result['CNt'] = data['CNt_list']
    result['B'] = data['B_list']
    return result


def add_theoreticals_to_segment(
    segment_df, cellularity, tumor_ploidy, normal_ploidy, is_female,
):
    assert isinstance(segment_df, pd.DataFrame)
    assert {'CNn', 'CNt', 'B'}.issubset(segment_df.columns)

    #def depthr_getter(row):
    #    return theoretical_depth_ratio(
    #        row['CNt'], cellularity, tumor_ploidy, row['CNn'], normal_ploidy, 
    #        tumor_avg_depth=1, normal_avg_depth=1,
    #    )

    #def baf_getter(row):
    #    if np.isnan(row['B']):
    #        return np.nan
    #    else:
    #        return theoretical_baf(row['CNt'], row['B'], cellularity, row['CNn'])

#    theo_depthr_list = segment_df.apply(depthr_getter, axis=1)
#    theo_baf_list = segment_df.apply(baf_getter, axis=1)
#
#    return segment_df.assign(
#        **{'depthratio_predicted': theo_depthr_list, 'baf_predicted': theo_baf_list}
#    )

    #segment_df['depthratio_predicted'] = segment_df.apply(depthr_getter, axis=1)
    segment_df['depthratio_predicted'] = theoretical_depth_ratio(
        segment_df['CNt'], 
        cellularity, 
        tumor_ploidy, 
        segment_df['CNn'], 
        normal_ploidy, 
        tumor_avg_depth=1, 
        normal_avg_depth=1,
    )

    #segment_df['baf_predicted'] = segment_df.apply(baf_getter, axis=1)
    segment_df['baf_predicted'] = theoretical_baf(
        segment_df['CNt'], 
        segment_df['B'], 
        cellularity, 
        segment_df['CNn'],
    )

    return segment_df
    

###################################


def _get_mean_ploidy_from_values(start0_list, end0_list, CNt_list):
    lengths = pd.Series(end0_list) - pd.Series(start0_list)
    return np.average(CNt_list, weights=lengths)


def get_mean_ploidy(seg_df):
    seg_lengths = seg_df.End - seg_df.Start
    numerator = (seg_df.CNt * seg_lengths).sum()
    denominator = seg_lengths.sum()
    return numerator / denominator


def get_tumor_mean_ploidy(CNt_annotated_gr):
    return np.average(CNt_annotated_gr.CNt, weights=CNt_annotated_gr.lengths())


def theoretical_somatic_vaf(CNt, CNm, cellularity, ccf, CNn=2):
    """Assumes no overlapping subclonal CNV.
    Args:
        CNm: Mutated copy number
        CNn: Germline copy number (should be 1 in haploid area such as male sex chromosome)
        ccf: Cancer cell fraction
    """
    numerator = CNm * cellularity * ccf
    denominator = CNt * cellularity + CNn * (1 - cellularity)
    return numerator / denominator


def ccf_from_params(vaf, CNt, CNm, cellularity, CNn):
    """Reverse of 'theoretical_somatic_vaf'"""
    numerator = vaf * (CNt * cellularity + CNn * (1 - cellularity))
    denominator = CNm * cellularity
    return numerator / denominator


def CNm_likelihood_difference(observed_vaf, CNt, CNm, cellularity, CNn):
    expected_vaf = theoretical_somatic_vaf(
        CNt=CNt, CNm=CNm, cellularity=cellularity, ccf=1, CNn=CNn,
    )
    return -abs(expected_vaf - observed_vaf)


def CNm_likelihood_binom(total_read, var_read, CNt, CNm, cellularity, CNn):
    expected_vaf = theoretical_somatic_vaf(
        CNt=CNt, CNm=CNm, cellularity=cellularity, ccf=1, CNn=CNn,
    )
    return scipy.stats.binom.pmf(var_read, total_read, expected_vaf)


def get_ccf_CNm(vaf, cellularity, CNt, CNB=None, likelihood='diff', total_read=None, var_read=None, CNn=2):
    """Notes:
        - Assumes no overlapping subclonal CNV
        - If subclonal, mutation may reside in only one copy. CNm must be 1.
        - If clonal, CNm can be any value between 1 and max(CNA, CNB).
    Args:
        CNn: Should be 1 in haploid area such as male sex chromosome
        CNB: None when haploid
    """
    assert likelihood in ('diff', 'binom')

    # NA cases
    if (
        CNt == 0 or
        np.isnan(vaf)
    ):
        return CCFInfo(None, None, None)

    # get CNA
    if CNB is None:
        CNA = None
    else:
        CNA = CNt - CNB
    # max_subclonal_vaf results when CNm == 1 and ccf == 1
    max_subclonal_vaf = theoretical_somatic_vaf(
        CNt=CNt, CNm=1, cellularity=cellularity, ccf=1, CNn=CNn,
    )

    if vaf >= max_subclonal_vaf:
        # Assumes clonal mutation
        ccf = 1
        max_CNm = CNt if (CNB is None) else max(CNA, CNB)
        CNm_candidates = list(range(1, max_CNm + 1))
        if likelihood == 'diff':
            CNm = max(
                CNm_candidates, 
                key=(lambda x: CNm_likelihood_difference(
                    observed_vaf=vaf, CNt=CNt, CNm=x, cellularity=cellularity, CNn=CNn,
                )),
            )
        elif likelihood == 'binom':
            CNm = max(
                CNm_candidates, 
                key=(lambda x: CNm_likelihood_binom(
                    total_read=total_read, 
                    var_read=var_read, 
                    CNt=CNt, 
                    CNm=x, 
                    cellularity=cellularity, 
                    CNn=CNn,
                )),
            )
    else:
        # Assumes subclonal mutation
        CNm = 1
        ccf = ccf_from_params(vaf=vaf, CNt=CNt, CNm=CNm, cellularity=cellularity, CNn=CNn)

    # get which allele is mutated
    if CNB is None:
        mutated_allele = None
    else:
        if CNA == CNB:
            mutated_allele = None
        else:
            if CNm > min(CNA, CNB):
                if CNA > CNB:
                    mutated_allele = 'A'
                else:
                    mutated_allele = 'B'
            else:
                mutated_allele = None

    return CCFInfo(CNm, ccf, mutated_allele)


#################################
# depth df processing functions #
#################################

def get_gc_breaks(n_bin):
    breaks = np.linspace(0, 1, (n_bin + 1), endpoint=True)
    breaks[0] = -0.001
    intervals = pd.IntervalIndex.from_breaks(breaks, closed='right')

    return breaks, intervals


def handle_outliers_with_preset_region(depth_df, outlier_region):
    outlier_region = arg_into_gr(outlier_region)
    depth_gr = arg_into_gr(depth_df)

    output_depth_df = depth_gr.count_overlaps(outlier_region, overlap_col='NumberOverlaps').df
    output_depth_df['excluded'] = output_depth_df['NumberOverlaps'] > 0
    output_depth_df.drop('NumberOverlaps', axis=1, inplace=True)

    return output_depth_df


# 28.75 sec
def handle_outliers_with_preset_included_region_old(depth_df, included_region):
    included_region = arg_into_gr(included_region)
    depth_gr = arg_into_gr(depth_df)

    output_depth_df = depth_gr.count_overlaps(included_region, overlap_col='NumberOverlaps').df
    output_depth_df['excluded'] = (output_depth_df['NumberOverlaps'] == 0)
    output_depth_df.drop('NumberOverlaps', axis=1, inplace=True)

    return output_depth_df


# 18.85 sec
def handle_outliers_with_preset_included_region(depth_df, included_region):
    assert isinstance(depth_df, pd.DataFrame)

    depth_df['_index'] = range(depth_df.shape[0])

    included_region = arg_into_gr(included_region)
    depth_gr = pr.PyRanges(depth_df)

    overlap = depth_gr.overlap(included_region)

    depth_df['excluded'] = True
    depth_df.iloc[overlap._index, depth_df.columns.get_loc('excluded')] = False

    return depth_df


def handle_outliers(depth_df, trim_limits=None, lower_cutoff=None, upper_cutoff=None, exclude_y=False):
    """Args:
        trim_limits: Tuple with length 2 (lower proportion, upper proportion).
            Used as an argument to "scipy.stats.mstats.trim".
            Example: (0.05, 0.05)
    Returns:
        Tuple (depth_df, selector)
            depth_df: The input dataframe, with a new column "excluded" added
            selector: boolean array indicating non-outlier rows
        Outlier handling processes:
            1) Removes/masks oringinal nan values
            2) Removes/masks based on relative proportions
            3) Removes/masks based on absolute cutoff values
    """
    assert 'mean_depth' in depth_df.columns

    # remove nan
    selector = depth_df['mean_depth'].notna()

    # remove y
    if exclude_y:
        selector = np.logical_and(
            selector, 
            ~depth_df['Chromosome'].isin(['Y', 'chrY']),
        )

    # trim by proportion
    if trim_limits is not None:
        trimmed_depth = scipy.stats.mstats.trim(
            depth_df['mean_depth'], limits=trim_limits, relative=True,
        )
        selector = np.logical_and(selector, ~trimmed_depth.mask)

    # trim by absolute cutoff
    if lower_cutoff is not None:
        selector = np.logical_and(selector, depth_df['mean_depth'] >= lower_cutoff)
    if upper_cutoff is not None:
        selector = np.logical_and(selector, depth_df['mean_depth'] <= upper_cutoff)

    # result
    selector = selector.to_numpy()
    #depth_df.loc[~selector, 'mean_depth'] = np.nan
    depth_df['excluded'] = ~selector

    return depth_df, selector


def get_depth_df_from_file(chroms, start0s, end0s, depth_file_path):
    preset_depth_df = pd.read_csv(
        depth_file_path, 
        sep='\t', 
        names=('Chromosome', 'Start', 'End', 'depth'), 
        dtype={'Chromosome': str, 'Start': int, 'End': int, 'depth': float},
    )
    preset_depth_gr = pr.PyRanges(preset_depth_df)
    result_gr = pr.PyRanges(chromosomes=chroms, starts=start0s, ends=end0s)
    result_gr = pyranges_helper.join(result_gr, preset_depth_gr, how='left', merge='weighted_mean', as_gr=True)
    return result_gr


# NOT USED
def get_gc_depth_data(depth_df, n_bin):
    assert 'GC' in depth_df.columns

    breaks = get_gc_breaks(n_bin)
    cutresult = pd.cut(depth_df['GC'], bins=breaks)

    gcbin_average_depths = dict()
    for intv, subdf in depth_df.groupby(cutresult, axis=0):
        # since "cutresult" is a categorical series, iteration over groupby object yields all bins including empty ones
        if subdf.shape[0] == 0:
            gcbin_average_depths[intv] = np.nan
        else:
            gcbin_average_depths[intv] = np.average(subdf['mean_depth'], weights=(subdf['End'] - subdf['Start']))
                # weighted average depth for the given GC bin

    gcbin_average_depths = pd.Series(gcbin_average_depths)
    gcbin_norm_average_depths = gcbin_average_depths / gcbin_average_depths.mean(skipna=True)

    return gcbin_average_depths, gcbin_norm_average_depths, cutresult


def get_gcbin_data(depth_df, gc_breaks, make_raw_data=False):
    """This does not apply weights when calculating averages. Suitable for even-sized bins"""
    # cut by gc bins
    cutresult = pd.cut(depth_df['GC'], bins=gc_breaks, include_lowest=False)
    all_intervals = list(cutresult.dtype.categories)

    cutresult_filtered = cutresult.loc[~depth_df['excluded']]
    depth_df_filtered = depth_df.loc[~depth_df['excluded'], :]

    # make results
    if make_raw_data:
        gcbin_depth_data = {intv: list() for intv in all_intervals}
        for intv, depth in zip(cutresult_filtered, depth_df_filtered['mean_depth']):
            gcbin_depth_data[intv].append(depth)

        gcbin_average_depths = {
            intv: (np.nan if len(data) == 0 else np.mean(data))
            for intv, data in gcbin_depth_data.items()
        }
    else:
        gcbin_depth_data = None
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)

            gcbin_average_depths = dict()
            for intv in all_intervals:
                gcbin_average_depths[intv] = np.mean(
                    depth_df_filtered['mean_depth'].loc[(cutresult_filtered == intv).array]
                )

    gcbin_average_depths = pd.Series(gcbin_average_depths)
    gcbin_norm_average_depths = gcbin_average_depths / gcbin_average_depths.mean(skipna=True)

    return gcbin_depth_data, gcbin_average_depths, gcbin_norm_average_depths, cutresult


#def get_normal_wgs_depth_cutoffs(depth_df):
#    mean = depth_df['mean_depth'].mean()

def depth_df_gc_addition(depth_df, gc_df=None, refver=None, fasta=None, window=None):
    if 'GC' in depth_df:
        depth_df = depth_df.drop('GC', axis=1)

    if gc_df is None:
        result = gcfraction.add_gc_calculating(
            depth_df, refver=refver, fasta=fasta, window=window,
        )
    else:
        gc_df = arg_into_df(gc_df)
        assert 'GC' in gc_df.columns
        assert gc_df.index.names == ['Chromosome', 'Start', 'End']
            # gc_df must have a MultiIndex ['Chromosome', 'Start', 'End']

        result = depth_df.join(gc_df, on=['Chromosome', 'Start', 'End'], how='left')

    return result


@deco.get_deco_num_set_differently(('fasta', 'refver', 'gc_df'), 1)
@deco.get_deco_num_set_differently(('outlier_region', 'included_region'), 1)
def postprocess_depth_df(
    depth_df, 
    *,
    refver=None,
    fasta=None, 
    gc_df=None,

    gc_window=None,

    outlier_region=None,
    included_region=None,
    preset_cutoffs=None,
    lower_cutoff=None,
    upper_cutoff=None,
    trim_limits=None,
    exclude_y=False,

    add_norm_depth=False,
    add_gccorr_depth=False,

    #nan_lower_cutoff=None,
    #nan_upper_cutoff=None,

    #gcdata_trim_limits=(0.05, 0.05),
    #gcdata_lower_cutoff=None,
    #gcdata_upper_cutoff=None,

    n_gcbin=100, 
    as_gr=True, 

    verbose=False,
):
    """Args:
        fasta, refver: Only used for gc fraction calculation. Only one of
            "fasta" or "refver" must be set. Using "refver" is recommended.

    Returns:
        postprocessed df & gcbin average depth dict
        Adds 4 columns:
            - GC
            - norm_mean_depth
            - gc_corrected_mean_depth
            - sequenza_style_norm_mean_depth
    """
    if verbose:
        def printlog(msg):
            funcname = inspect.stack()[0].function
            common.print_timestamp(f'{funcname}: {msg}')
    else:
        def printlog(msg):
            pass

    # sanity check
    depth_df = arg_into_df(depth_df)
    assert 'mean_depth' in depth_df.columns, f'"depth_df" must include a column named "mean_depth"'
    assert preset_cutoffs in ('wgs', 'normal_wgs', 'panel', None)

    # calculate average depth
    avg_depth = depth_df['mean_depth'].mean()

    # set outlier cutoffs
    if preset_cutoffs == 'wgs':
        printlog(f'Running "postprocess_depth_df" function with preset cutoff mode "wgs"')
        lower_cutoff = 0
        upper_cutoff = np.inf
    elif preset_cutoffs == 'panel':
        printlog(f'Running "postprocess_depth_df" function with preset cutoff mode "panel"')
        lower_cutoff = 50
        upper_cutoff = np.inf
    elif preset_cutoffs == 'normal_wgs':
        printlog(f'Running "postprocess_depth_df" function with preset cutoff mode "normal_wgs"')
        lower_cutoff = avg_depth * 0.7
        upper_cutoff = avg_depth * 1.3

    # add GC fractions
    printlog(f'Adding GC fraction values')
    depth_df = depth_df_gc_addition(
        depth_df, gc_df=gc_df, refver=refver, fasta=fasta, window=gc_window,
    )

    ## replace outliers with nan; separately for output and gc data 
    # mark outlier positions with 'excluded' column
    printlog(f'Marking outlier positions')
    if outlier_region is not None:
        output_depth_df = handle_outliers_with_preset_region(depth_df, outlier_region)
    elif included_region is not None:
        output_depth_df = handle_outliers_with_preset_included_region(depth_df, included_region)
    else:
        output_depth_df, selector = handle_outliers(
            depth_df, 
            trim_limits=trim_limits, 
            lower_cutoff=lower_cutoff, 
            upper_cutoff=upper_cutoff,
            exclude_y=exclude_y,
        )
    #gcdata_depth_df, selector2 = handle_outliers(
    #    depth_df, trim_limits=None, lower_cutoff=0, upper_cutoff=(avg_depth * 10),
    #)

    # get gc depth data
    printlog(f'Getting depth statistics stratified by gc ranges')
    gc_breaks, gc_intervals = get_gc_breaks(n_gcbin)
    (
        gcbin_depth_data, 
        gcbin_average_depths, 
        gcbin_norm_average_depths, 
        cutresult 
    #) = get_gcbin_data(gcdata_depth_df, gc_breaks, make_raw_data=False)
    ) = get_gcbin_data(
        output_depth_df, 
        gc_breaks, 
        make_raw_data=False,
    )
    cutresult_idx = common.get_indexes_of_array(cutresult, gcbin_norm_average_depths.index)

    # norm_mean_depth
    if add_norm_depth:
        printlog(f'Getting normalized depths')
        filterd_depth_df = output_depth_df.loc[~output_depth_df['excluded'], :]
        avg_depth_wo_outliers = common.nanaverage(
            filterd_depth_df['mean_depth'].to_numpy(), 
            weights=(filterd_depth_df['End'] - filterd_depth_df['Start']).to_numpy(),
        )
        output_depth_df['norm_mean_depth'] = (output_depth_df['mean_depth'] / avg_depth_wo_outliers).array

    # gc_corrected_mean_depth
    if add_gccorr_depth:
        printlog(f'Getting GC-corrected depths')
        output_depth_df['gc_corrected_mean_depth'] = (
            output_depth_df['mean_depth'].array 
            / gcbin_norm_average_depths.iloc[cutresult_idx].array
        )

    # sequenza_style_norm_mean_depth
    printlog(f'Getting GC-corrected normalized depths (sequenza style)')
    output_depth_df['sequenza_style_norm_mean_depth'] = (
        output_depth_df['mean_depth'].array 
        / gcbin_average_depths.iloc[cutresult_idx].array
    )

    # result
    if as_gr:
        output_depth_df = pr.PyRanges(output_depth_df, int64=False)

    printlog(f'Finished')

    return output_depth_df, gcbin_average_depths


def get_processed_depth_df(
    bam_path, fasta, region_gr, 
    gc_window=None, 
    preset_cutoffs='wgs',
    lower_cutoff=None,
    upper_cutoff=None,
    trim_limits=None,
    donot_subset_bam=False, 
    as_gr=False,
):
    # run mosdepth
    depth_df, _ = libmosdepth.run_mosdepth(
        bam_path,
        t=8, 
        region_gr=region_gr, 
        donot_subset_bam=donot_subset_bam, 
        as_gr=False, 
        load_perbase=False,
    )

    # postprocess
    depth_gr, gcbin_average_depths = postprocess_depth_df(
        depth_df, 

        fasta=fasta, 
        gc_window=gc_window,

        preset_cutoffs=preset_cutoffs,
        lower_cutoff=lower_cutoff,
        upper_cutoff=upper_cutoff,
        trim_limits=trim_limits,

        n_gcbin=100,
        as_gr=as_gr, 
    )

    return depth_gr, gcbin_average_depths


def make_depth_ratio(
    tumor_depth_df, normal_depth_df, 
    make_depthratio_mystyle=False,
    as_gr=False,
    #verbose=False,
):
    def make_ratio(numerator, denominator):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)
            arr = numerator / denominator
            arr[np.isinf(arr)] = np.nan
            return arr

    def copy_annotations(result, input_df, prefix):
        annot_cols = get_genome_df_annotcols(input_df).drop(['GC', 'excluded'])
        for colname in annot_cols:
            result[f'{prefix}_{colname}'] = input_df.loc[:, colname].array

    # sanity check
    tumor_depth_df = arg_into_df(tumor_depth_df)
    normal_depth_df = arg_into_df(normal_depth_df)

    genome_df_sanitycheck(tumor_depth_df)
    genome_df_sanitycheck(normal_depth_df)

    required_cols = {
        'mean_depth', 
        #'norm_mean_depth', 
        #'gc_corrected_mean_depth', 
        'GC', 
        'sequenza_style_norm_mean_depth',
    }
    assert all(
        required_cols.issubset(df.columns)
        for df in (tumor_depth_df, normal_depth_df)
    )

    # set logger
    #logger = (LOGGER_DEBUG if verbose else LOGGER_INFO)

    # tumor/normal dfs shape check
    compared_columns = ['Chromosome', 'Start', 'End', 'GC', 'excluded']
    tumor_common_cols = tumor_depth_df.loc[:, compared_columns]
    normal_common_cols = normal_depth_df.loc[:, compared_columns]
    if not (normal_common_cols == tumor_common_cols).all(axis=None):
        raise Exception(
            f'The columns {compared_columns} of the two input dataframes must be identical.'
        )

    # main
    result = tumor_common_cols

    # make depth ratios
    result['depth_ratio_sequenzastyle'] = make_ratio(
        tumor_depth_df['sequenza_style_norm_mean_depth'].array, 
        normal_depth_df['sequenza_style_norm_mean_depth'].array,
    )
    if make_depthratio_mystyle:
        result['depth_ratio_mystyle'] = make_ratio(
            tumor_depth_df['gc_corrected_mean_depth'].array, 
            normal_depth_df['gc_corrected_mean_depth'].array,
        )

    # copy annotation columns from input dfs
    copy_annotations(result, tumor_depth_df, 'tumor')
    copy_annotations(result, normal_depth_df, 'normal')

    if as_gr:
        return pr.PyRanges(result)
    else:
        return result


def merge_normal_tumor_bafs(tumor_baf_gr, normal_baf_gr, region_gr, as_gr=False):
    # convert df types
    #tumor_baf_gr = arg_into_gr(tumor_baf_gr)
    #normal_baf_gr = arg_into_gr(normal_baf_gr)
    #region_gr = arg_into_gr(region_gr)

    # sanity check
    assert 'baf' in tumor_baf_gr.columns
    assert 'baf' in normal_baf_gr.columns

    # remove annotation columns
    result = region_gr[[]]  

    # tumor
    result = pyranges_helper.join_newest(
        result, tumor_baf_gr,
        how='left', find_nearest=False, merge='mean',
        as_gr=False,
    )
    result.rename(columns={'baf': 'tumor_baf'}, inplace=True)
    result = pr.PyRanges(result)

    # normal
    result = pyranges_helper.join_newest(
        result, normal_baf_gr,
        how='left', find_nearest=False, merge='mean',
        as_gr=False,
    )
    result.rename(columns={'baf': 'normal_baf'}, inplace=True)

    # return
    if as_gr:
        return pr.PyRanges(result)
    else:
        return result


def handle_merged_baf_df(merged_baf_df):
    """Find regions with normal sample baf deviating from 0.5"""
    pass


@deco.get_deco_num_set_differently(('refver', 'chromdict'), 1)
def upsize_depth_df_bin(depth_df, size, refver=None, chromdict=None, annot_cols=None):
    """Input must be sorted"""
    def get_chrom_df_groupers(chrom_df, chromlen, new_binsize, old_binsize):
        nrow = chrom_df.shape[0]
        assert nrow == np.ceil(chromlen / old_binsize).astype(int)

        n_newbin = np.ceil(chromlen / new_binsize).astype(int)
        olddf_idx_ends = np.fromiter(
            (
                min(int((new_binsize * x) / old_binsize), nrow)
                for x in range(1, n_newbin + 1)
            ),
            dtype=int,
        )
        group_lengths = np.diff(olddf_idx_ends, prepend=0)
        groupers = np.repeat(range(len(group_lengths)), group_lengths)
        return groupers, n_newbin

    def upsize_chrom_df(chrom_df, new_binsize, groupers, n_newbin, chromlen, annot_cols):
        starts = np.arange(n_newbin) * new_binsize
        ends = np.concatenate((starts[1:], [chromlen]))
        nonannot_subdf = pd.DataFrame({
            'Chromosome': chrom_df['Chromosome'][0],
            'Start': starts,
            'End': ends,
        })

        annot_subdf = chrom_df.groupby(by=groupers, axis=0, sort=False)[annot_cols].mean()
        annot_subdf.reset_index(inplace=True, drop=True)

        result = pd.concat([nonannot_subdf, annot_subdf], axis=1)

        if 'excluded' in chrom_df:
            merged_exclude = pd.Series(chrom_df['excluded']).groupby(
                by=groupers, sort=False,
            ).any()
            result['excluded'] = merged_exclude.to_numpy()

        return result

    # arg handling
    if annot_cols is None:
        #annot_cols = get_genome_df_annotcols(depth_df)
        annot_cols = [
            k for (k, v) in depth_df.dtypes.to_dict().items()
            if (
                (str(v) != 'bool')
                and (k not in ['Chromosome', 'Start', 'End'])
            )
        ]
    if chromdict is None:
        chromdict = common.DEFAULT_CHROMDICTS[refver]

    dfs_bychrom = pyranges_helper.group_df_bychrom(depth_df)
    new_binsize = size

    # get input df bin size
    first_df = next(iter(dfs_bychrom.values()))
    assert first_df.shape[0] >= 2
    old_binsize = first_df['End'][0] - first_df['Start'][0]
    if size <= old_binsize:
        raise Exception(f'Upsized bin size ({size}) must be greater than old bin size ({old_binsize})')

    # process each subdf by chrom
    newbin_dfs_bychrom = dict()
    for chrom, chrom_df in dfs_bychrom.items():
        chromlen = chromdict[chrom]
        groupers, n_newbin = get_chrom_df_groupers(chrom_df, chromlen, new_binsize, old_binsize)
        newbin_dfs_bychrom[chrom] = upsize_chrom_df(
            chrom_df, new_binsize, groupers, n_newbin, chromlen, annot_cols
        )

    # sort and concat
    sorted_chroms = sorted(newbin_dfs_bychrom.keys(), key=(lambda x: chromdict.contigs.index(x)))
    result = pd.concat((newbin_dfs_bychrom[x] for x in sorted_chroms), axis=0)

    return result


#########################################
# plotting GCfrac - depth relationships #
#########################################

def arghandler_ylims(ydata, ylims):
    if ylims is None:
        return None
    else:
        return (
            min(ylims[0], min(ydata)), 
            max(ylims[1], max(ydata)),
        )


def plot_gc_depth(
    fasta, bam_path, 
    gc_window=50, 
    xwidth=1000, 
    chrom=None,
    start0=None,
    end0=None,
    return_gr=True, 
    depth_ylims=None, 
    gc_ylims=None,
    alpha=0.6,
    linewidth=1,
):
    """Intended to be used with short regions (~10kb)
    GC fraction value is calculated for every single base, with a window "gc_window" wide.
    """
    pad_left = int(gc_window / 2)
    if gc_window % 2 == 0:
        pad_right = pad_left
    else:
        pad_right = pad_left + 1    

    if any(x is None for x in [chrom, start0, end0]):
        chrom = random.choice(fasta.references[:24])
        chromlen = fasta.lengths[fasta.references.index(chrom)]
        start0 = random.choice(range(pad_left, (chromlen - xwidth - pad_right) + 1))
        end0 = start0 + xwidth

    fetch_start0 = start0 - pad_left
    fetch_end0 = end0 + pad_right
    
    ### GC ###
    seq = fasta.fetch(chrom, fetch_start0, fetch_end0)
    x = list()
    y_gc = list()
    for idx in range(0, end0 - start0):
        subseq = seq[idx:(idx + gc_window)]
        gcfrac = Bio.SeqUtils.gc_fraction(subseq)
        x.append(start0 + idx)
        y_gc.append(gcfrac)

    # modify ylims
    gc_ylims = arghandler_ylims(y_gc, gc_ylims)
    
    # plot
    fig, ax_gc = plt.subplots(figsize=(20, 5))
    if gc_ylims is not None:
        ax_gc.set_ylim(*gc_ylims)
    ax_gc.set_title(f'{chrom}:{start0 + 1:,}-{end0:,}')
    ax_gc.plot(
        x, y_gc, 
        color='green', marker=None, linestyle='-', markersize=1, linewidth=linewidth, alpha=alpha
    )
    
    ### depth ###
    region_gr = pr.from_dict({'Chromosome': [chrom], 'Start': [start0], 'End': [end0]})
    gr, gr_perbase = libmosdepth.run_mosdepth(
        bam_path, 
        use_median=False,
        region_gr=region_gr, 
        as_gr=True, 
        load_perbase=True,
    )
    y_depth = list(gr_perbase.depth)

    # modify ylims
    depth_ylims = arghandler_ylims(y_depth, depth_ylims)

    # plot
    ax_depth = ax_gc.twinx()
    if depth_ylims is not None:
        ax_depth.set_ylim(*depth_ylims)
    ax_depth.plot(
        x, y_depth, 
        color='black', marker=None, linestyle='-', markersize=1, linewidth=linewidth, alpha=alpha
    )

    # return
    if return_gr:
        gr_perbase.GC = y_gc
        return gr_perbase


def plot_gc_distribution(
    depth_df, 

    n_gcbin=100,

    trim_limits=(0.05, 0.05),
    lower_cutoff=None,
    upper_cutoff=None,

    histogram_mode='2d',
    scatter_ylims=None,
    heatmap_ylims=None, 
    heatmap_vmax=None,
):
    def make_average_scatter(gcbin_average_depths, ylims=None):
        x = list()
        y = list()
        for key, val in gcbin_average_depths.items():
            x.append(key.mid)
            y.append(np.mean(val))

        ylims = arghandler_ylims(y, ylims)

        fig, ax = plt.subplots(figsize=(10, 5))
        if ylims is not None:
            ax.set_ylim(*ylims)
        ax.set_xlabel('GC fraction')
        ax.set_ylabel('mean depth')
        ax.plot(x, y, 'ko', markersize=2)
        return fig

    def make_heatmap(gcbin_depth_data, ylims=None, vmax=None):
        ydata = list(itertools.chain.from_iterable(gcbin_depth_data.values()))
        if ylims is None:
            ylims = (min(ydata), max(ydata)) 
        else:
            ylims = arghandler_ylims(ydata, ylims)

        # set breaks
        breaks = np.linspace(ylims[0], ylims[1], 100, endpoint=True)

        # get histogram
        data = dict()
        for gc_interval, val in gcbin_depth_data.items():
            data[gc_interval] = np.histogram(val, bins=breaks, density=True)[0]
        df_plot = pd.DataFrame.from_dict(data)
        df_plot.index = pd.IntervalIndex.from_breaks(np.round(breaks, 1))

        # draw
        fig, ax = plt.subplots(figsize=(15, 10))
        sns.heatmap(df_plot, ax=ax, vmax=vmax)
        ax.set_xlabel('GC fraction')
        ax.set_ylabel('mean depth')
        return fig

    def make_heatmap_hist2d(depth_df, n_gcbin, ylims=None, vmax=None):
        # set breaks
        if ylims is None:
            ylims = (min(depth_df['mean_depth']), max(depth_df['mean_depth'])) 
        else:
            ylims = arghandler_ylims(depth_df['mean_depth'], ylims)
        depth_breaks = np.linspace(ylims[0], ylims[1], 100, endpoint=True)
        #gc_breaks = np.linspace(0, 1, 101, endpoint=True)
        gc_breaks, gc_intervals = get_gc_breaks(n_gcbin)

        # get histogram
        H, xedges, yedges = np.histogram2d(
            x=depth_df['mean_depth'], y=depth_df['GC'], bins=(depth_breaks, gc_breaks), density=True,
        )
        plot_df = pd.DataFrame(
            H, 
            index=pd.IntervalIndex.from_breaks(np.round(xedges, 1)), 
            columns=pd.IntervalIndex.from_breaks(np.round(yedges, 2)),
        )
      
        # draw
        fig, ax = plt.subplots(figsize=(15, 10))
        sns.heatmap(plot_df, ax=ax, vmax=vmax)
        ax.set_xlabel('GC fraction')
        ax.set_ylabel('mean depth')
        return fig

    # main
    assert histogram_mode in ('1d', '2d')

    depth_df, selector = handle_outliers(depth_df, trim_limits=trim_limits, lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff)
    trimmed_depth_df = depth_df.iloc[selector, :]
    gc_breaks, gc_intervals = get_gc_breaks(n_gcbin)
    gcbin_depth_data, gcbin_average_depths, gcbin_norm_average_depths, cutresult = get_gcbin_data(trimmed_depth_df, gc_breaks, make_raw_data=True)
    fig_avg_scatter = make_average_scatter(gcbin_average_depths, ylims=scatter_ylims)
    if histogram_mode == '2d':
        fig_heatmap = make_heatmap_hist2d(trimmed_depth_df, n_gcbin, ylims=heatmap_ylims, vmax=heatmap_vmax)
    elif histogram_mode == '1d':
        fig_heatmap = make_heatmap(gcbin_depth_data, ylims=heatmap_ylims, vmax=heatmap_vmax)

    return fig_avg_scatter, fig_heatmap



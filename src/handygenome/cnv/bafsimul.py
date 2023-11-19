import pickle
import functools
import itertools
import multiprocessing

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize
import scipy.interpolate

import handygenome
import handygenome.tools as tools
import handygenome.stats as stats
import handygenome.cnv.baf as libbaf
import handygenome.plot.misc as plotmisc


EMPIRICAL_ALPHASUM = 193.8113381099127
EMPIRICAL_OMEGA = -0.5547888600614372
INTERP_SAVEPATH = handygenome.USERDATA_DIR / 'baf_interpolator.pickle'


#########
# basic #
#########

def make_simulated_total_depths(mean_depth, N, omega=None, rng=None, mode='gpois'):
    """Only non-zero values"""
    assert mode in ('pois', 'gpois', 'repeat')
    
    if mode == 'repeat':
        return np.repeat(np.rint(mean_depth).astype(int), N)

    ########

    if rng is None:
        rng = np.random.default_rng()

    result = list()
    while True:
        size = N - sum(map(len, result))

        if mode == 'pois':
            arr = rng.poisson(lam=mean_depth, size=size)
        elif mode == 'gpois':
            arr = stats.rgpois(size, lam=mean_depth, omega=omega)

        arr = arr[arr > 0]
        result.append(arr)
        if sum(map(len, result)) == N:
            break

    return np.concatenate(result)


def make_simulated_allele_depths(total_depths, allele_portions, alpha_sum=None, rng=None):
    assert np.isclose(allele_portions.sum(), 1)
    if alpha_sum is not None:  # dirichlet multinomial
        dirichlet_alpha = allele_portions * alpha_sum
        allele_depths = tools.dirichlet_multinomial_rvs(total_depths, dirichlet_alpha, rng=rng)
    else:
        allele_depths = rng.multinomial(n=total_depths, pvals=allele_portions)

    return allele_depths


def simulate_variant_data(
    mean_depth, allele_portions, size, 
    omega=EMPIRICAL_OMEGA,
    alpha_sum=EMPIRICAL_ALPHASUM, 
    dirichlet=True, 
    rng=None, 
    mode='gpois',
    ALT_biases=None,
): 
    """Args:
        allele_portions: Must be ordered. First item should indicate REF, next one ALT1, then ALT2, ...
        REF_bias: length must be len(allele_portions) - 1. First item indicates ALT1/REF bias, next one ALT2/REF, ...

    Help:
        - "allele_portions" is first normalized such that it sums to 1
        - The order of "allele_portions" is not important
        - Larger alpha_sum leads to less variation
        - Total depths are generated from poisson distribution
        - If "dirichlet" is False, plain multinomial distribution is used.
    """
    # allele_portions sanitycheck
    allele_portions = np.asarray(allele_portions)
    assert allele_portions.ndim == 1

    # apply ALT_biases to allele_portions
    if ALT_biases is not None:
        allele_portions = allele_portions * np.insert(ALT_biases, 0, 1)

    # normalize allele_portions such that sum becomes 1
    allele_portions = allele_portions / allele_portions.sum()

    # make total depths
    if rng is None:
        rng = np.random.default_rng()

    total_depths = make_simulated_total_depths(mean_depth=mean_depth, N=size, rng=rng, omega=omega, mode=mode)

    if dirichlet:
        assert alpha_sum is not None
        allele_depths = make_simulated_allele_depths(total_depths, allele_portions, alpha_sum=alpha_sum, rng=rng)
    else:
        allele_depths = make_simulated_allele_depths(total_depths, allele_portions, alpha_sum=None)

    # make vafs and bafs
    vafs = allele_depths / total_depths[:, np.newaxis]
    bafs = libbaf.get_baf_from_vaf(vafs)

    simul_vardata = {
        'total_depth': total_depths,
        'allele_depth': allele_depths,
        'vaf': vafs,
        'baf': bafs,
    }
    return simul_vardata


def simulate_corrected_variant_data(
    mean_depth, 
    corrector_mean_depth,
    allele_portions, 
    size, 
    omega=EMPIRICAL_OMEGA,
    alpha_sum=EMPIRICAL_ALPHASUM, 
    #dirichlet=True, 
    rng=None, 
    mode='gpois',
    ALT_biases=None,
): 
    """Args:
        allele_portions: Must be ordered. First item should indicate REF, next one ALT1, then ALT2, ...
        REF_bias: length must be len(allele_portions) - 1. First item indicates ALT1/REF bias, next one ALT2/REF, ...

    Help:
        - "allele_portions" is first normalized such that it sums to 1
        - The order of "allele_portions" is not important
        - Larger alpha_sum leads to less variation
        - Total depths are generated from poisson distribution
        - If "dirichlet" is False, plain multinomial distribution is used.
    """
    # allele_portions preprocess
    allele_portions = np.asarray(allele_portions)
    assert allele_portions.ndim == 1

    # make corrector allele portions
    corrector_allele_portions = np.repeat(1, len(allele_portions))

    # apply ALT_biases to allele_portions
    if ALT_biases is not None:
        ALT_biases = np.asarray(ALT_biases)
        assert ALT_biases.shape == (len(allele_portions) - 1,)
        ALT_biases = np.insert(ALT_biases, 0, 1)

        allele_portions = allele_portions * ALT_biases
        corrector_allele_portions = corrector_allele_portions * ALT_biases

    # normalize allele_portions
    allele_portions = allele_portions / allele_portions.sum()
    corrector_allele_portions = corrector_allele_portions / corrector_allele_portions.sum()

    # make depths
    if rng is None:
        rng = np.random.default_rng()

    total_depths = make_simulated_total_depths(mean_depth=mean_depth, N=size, rng=rng, omega=omega, mode=mode)
    allele_depths = make_simulated_allele_depths(total_depths, allele_portions, alpha_sum=alpha_sum, rng=rng)

    corrector_total_depths = make_simulated_total_depths(mean_depth=corrector_mean_depth, N=size, rng=rng, omega=omega, mode=mode)
    corrector_allele_depths = make_simulated_allele_depths(corrector_total_depths, corrector_allele_portions, alpha_sum=alpha_sum, rng=rng)

    # make vafs and bafs
    vafs = allele_depths / total_depths[:, np.newaxis]
    corrector_vafs = corrector_allele_depths / corrector_total_depths[:, np.newaxis]

    selector = (corrector_vafs != 0).all(axis=1)
    assert selector.mean() > 0.99
    vafs = vafs[selector, :]
    corrector_vafs = corrector_vafs[selector, :]

    corrected_vafs = vafs / corrector_vafs
    corrected_vafs = corrected_vafs / corrected_vafs.sum(axis=1)[:, np.newaxis]
    bafs = libbaf.get_baf_from_vaf(vafs)
    corrected_bafs = libbaf.get_baf_from_vaf(corrected_vafs)

    simul_vardata = {
        'total_depth': total_depths,
        'allele_depth': allele_depths,
        'corrector_total_depth': corrector_total_depths,
        'corrector_allele_depth': corrector_allele_depths,
        'vaf': vafs,
        'corrector_vaf': corrector_vafs,
        'corrected_vaf': corrected_vafs,
        'baf': bafs,
        'corrected_baf': corrected_bafs,
    }
    return simul_vardata


########################
# error vaf simulation #
########################

def simulate_error_variants(
    error_rate=0.005,
    alpha_sum=10,
    mean_depth=40,
    size=10000,
):
    return simulate_variant_data(
        mean_depth=mean_depth, 
        allele_portions=[error_rate, 1 - error_rate], 
        size=size, 
        alpha_sum=alpha_sum,
    )


def show_simulated_error_baf(
    baf_rawdata=None,
    bins=np.arange(0, 0.505, 0.01),
    **kwargs,
):
    if baf_rawdata is None:
        baf_rawdata = np.squeez(simulate_error_variants(**kwargs)['baf'])

    fig, ax = plt.subplots()
    ax.hist(
        baf_rawdata,
        bins=bins,
        density=True,
    )
    return fig, ax, baf_rawdata


###################################
# empirical alpha_sum calculation #
###################################

DEFAULT_DIGITIZE_BINS = np.arange(0.3, 0.5 + 0.025, 0.025)
DEFAULT_DIGITIZE_BINS_VAF = np.arange(0.3, 0.7 + 0.025, 0.025)
TRIM_CUTOFF = 0.3

def trim_baf_rawdata(bafdata):
    """
    with cutoff == 0.2, solution cannot be found because kurtosis or skewness
    obtained from simulated data cannot reach that of real data
    """
    return bafdata[bafdata > TRIM_CUTOFF]


def trim_vaf_rawdata(vafdata):
    """
    with cutoff == 0.2, solution cannot be found because kurtosis or skewness
    obtained from simulated data cannot reach that of real data
    """
    return vafdata[(vafdata > TRIM_CUTOFF) & (vafdata < (1 - TRIM_CUTOFF))]


def get_bafdata_stats(bafdata, trim=True, digitize=True):
    if trim:
        bafdata = trim_baf_rawdata(bafdata)
    if digitize:
        bafdata = tools.digitize(bafdata, DEFAULT_DIGITIZE_BINS)
    decal = libbaf.decal_baf(bafdata)

    bafdata_stats = {
        'skew': scipy.stats.skew(bafdata),
        'kurt': scipy.stats.kurtosis(decal),
        'mean': np.mean(bafdata),
        'var': np.var(decal, ddof=0),
    }
    return bafdata_stats


def get_vafdata_stats(vafdata, trim=True, digitize=False):
    assert vafdata.ndim == 1

    if trim:
        vafdata = trim_vaf_rawdata(vafdata)
    if digitize:
        vafdata = tools.digitize(vafdata, DEFAULT_DIGITIZE_BINS_VAF)
    bafdata = libbaf.get_baf_from_vaf(vafdata)[:, 0]

    vafdata_stats = {
        'skew': scipy.stats.skew(bafdata),
        'kurt': scipy.stats.kurtosis(vafdata),
        'mean': np.mean(bafdata),
        'var': np.var(vafdata, ddof=0),
    }
    return vafdata_stats


def make_simulated_halfbaf_stats_subjob(
    mean_depth,
    size,
    alpha_sum,
    mode,
    omega,
    trim,
    digitize,
):
    simul_bafs = simulate_variant_data(
        mean_depth=mean_depth, 
        allele_portions=[1, 1], 
        size=size, 
        alpha_sum=alpha_sum, 
        dirichlet=True,
        mode=mode,
        omega=omega,
    )['baf'][:, 0]
    bafdata_stats = get_bafdata_stats(simul_bafs, trim=trim, digitize=digitize)

    return bafdata_stats


def make_simulated_halfbaf_stats_subjob_usevaf(
    mean_depth,
    size,
    alpha_sum,
    mode,
    omega,
    trim,
    digitize,
):
    simul_vafs = simulate_variant_data(
        mean_depth=mean_depth, 
        allele_portions=[1, 1], 
        size=size, 
        alpha_sum=alpha_sum, 
        dirichlet=True,
        mode=mode,
        omega=omega,
    )['vaf'][:, 0]  # REF vafs
    vafdata_stats = get_vafdata_stats(simul_vafs, trim=trim, digitize=digitize)

    return vafdata_stats


def make_simulated_halfbaf_stats(
    mean_depth, 
    size=100000,
    rep=100, 
    alphas=np.arange(100, 1000, 50),
    #bins=DEFAULT_DIGITIZE_BINS,
    use_baf=False,
    mode='gpois',
    omega=None,
    trim=True,
    digitize=False,
    verbose=True,
    nproc=1,
):
    """Simulate baf rawdata with true_baf=0.5 (germline het, ploidy=2) 
    with various alpha_sum values. From the simulated baf rawdata,
    skewness and kurtosis are calculated.
    """
    stat_keys = ('skew', 'kurt', 'mean', 'var')

    simul_halfbaf_stats = (
        {'alpha': np.asarray(alphas)}
        | {
            key: {'mean': list(), 'std': list()} 
            for key in stat_keys
        }
    )
    with multiprocessing.Pool(nproc) as pool:
        for alpha_sum in alphas:
            if verbose:
                print(f'alpha value {alpha_sum}', flush=True)

            subdata = {x: list() for x in stat_keys}
            args = (
                (mean_depth, size, alpha_sum, mode, omega, trim, digitize)
                for _ in range(rep)
            )
            if use_baf:
                mp_result = pool.starmap(make_simulated_halfbaf_stats_subjob, args)
            else:
                mp_result = pool.starmap(make_simulated_halfbaf_stats_subjob_usevaf, args)

            for stat_values in mp_result:
                for key, val in stat_values.items():
                    subdata[key].append(val)

            for key, val in subdata.items():
                simul_halfbaf_stats[key]['mean'].append(np.mean(val))
                simul_halfbaf_stats[key]['std'].append(np.std(val, ddof=0))

    for key in stat_keys:
        for subkey in tuple(simul_halfbaf_stats[key].keys()):
            simul_halfbaf_stats[key][subkey] = np.asarray(simul_halfbaf_stats[key][subkey])
        
    return simul_halfbaf_stats


def show_simulated_halfbaf_stats(
    simul_halfbaf_stats, 
    real_germline_bafs=None, 
    real_germline_vafs=None, 
    figsize=(30, 8),
    trim=True,
    digitize=False,
):
    fig, axd = plt.subplot_mosaic(
        [
            ['skew', 'kurt'],
            ['mean', 'var'],
        ],
        figsize=figsize,
    )
    for key, ax in axd.items():
        ax.plot(
            simul_halfbaf_stats['alpha'], 
            simul_halfbaf_stats[key]['mean'], 
            marker='o', 
            linewidth=0.5,
            label=key,
        )
        ax.vlines(
            simul_halfbaf_stats['alpha'], 
            simul_halfbaf_stats[key]['mean'] - simul_halfbaf_stats[key]['std'],
            simul_halfbaf_stats[key]['mean'] + simul_halfbaf_stats[key]['std'],
            color='black',
        )
        ax.set_xlabel('alpha sum')
        ax.set_ylabel(key)

    if real_germline_vafs is not None:
        stat_values = get_vafdata_stats(real_germline_vafs, trim=trim, digitize=digitize)
        for key, val in stat_values.items():
            axd[key].axhline(val)
    elif real_germline_bafs is not None:
        stat_values = get_bafdata_stats(real_germline_bafs, trim=trim, digitize=digitize)
        for key, val in stat_values.items():
            axd[key].axhline(val)

    for ax in axd.values():
        ax.legend()

    return fig, axd


def get_alpha_estimates(
    simul_halfbaf_stats, 
    real_germline_data, 
    use_baf=False,
    trim=True,
    digitize=False,
):
    if use_baf:
        realdata_stats = get_bafdata_stats(real_germline_data, trim=trim, digitize=digitize)
    else:
        realdata_stats = get_vafdata_stats(real_germline_data, trim=trim, digitize=digitize)

    def helper(simul_stats, realdata_stat):
        spl = scipy.interpolate.CubicSpline(
            simul_halfbaf_stats['alpha'],
            simul_stats - realdata_stat,
        )
        estimates = spl.roots(extrapolate=False)
        if len(estimates) == 0:
            return np.nan
        elif len(estimates) > 1:
            raise Exception(f'Unique root could not be found')
        else:
            return estimates[0]

    alpha_estimates = dict()
    for key, val in realdata_stats.items():
        alpha_estimates[key] = helper(
            simul_halfbaf_stats[key]['mean'], 
            realdata_stats[key],
        )

    return alpha_estimates


def make_empirical_alpha_sum(depth_rawdata, baf_rawdata, rep=100):
    """- Rawdata given as arguments must be from a germline wgs sample (e.g. peripheral blood)
    - ploidy must be 2
    """

    simul_halfbaf_stats = make_simulated_halfbaf_stats(
        mean_depth=np.nanmean(depth_rawdata), 
        size=len(trim_baf_rawdata(baf_rawdata)),
        rep=rep, 
        alphas=np.arange(100, 1000, 50),
    )
    alpha_estimates = get_alpha_estimates(simul_halfbaf_stats, baf_rawdata)
    #return np.mean([alpha_estimates['skew'], alpha_estimates['kurt']])
    return alpha_estimates['var']


def validate_empirical_alphasum(empirical_alphasum, realdata_bafs, realdata_depths):
    # prepare data
    realdata_bafs_decal = np.concatenate([realdata_bafs, 1 - realdata_bafs])

    simul_bafs = simulate_variant_data(
        mean_depth=np.nanmean(realdata_depths), 
        allele_portions=[1, 1], 
        size=len(trim_baf_rawdata(realdata_bafs)), 
        alpha_sum=empirical_alphasum, 
        dirichlet=True,
    )['baf']
    simul_bafs = np.squeeze(simul_bafs)
    simul_bafs_decal = np.concatenate([simul_bafs, 1 - simul_bafs])

    # make axd
    fig, axd = plt.subplot_mosaic(
        [
            ['hist_real', 'qq'],
            ['hist_simul', 'qq'],
        ],
        figsize=(30, 10),
        gridspec_kw=dict(
            #wspace=0.4,
        ),
    )

    # draw hist
    hist_bins = np.arange(0, 1.01, 0.01)
    _ = axd['hist_real'].hist(realdata_bafs_decal, bins=hist_bins)
    axd['hist_real'].set_title('real')

    hist_ys, bins, patches = axd['hist_simul'].hist(simul_bafs_decal, bins=hist_bins)
    axd['hist_simul'].set_title('simulated')

    ymax = np.max(hist_ys) * 1.2
    for key in ['hist_real', 'hist_simul']:
        axd[key].set_ylim(0, ymax)

    # draw qqplot
    plotmisc.qqplot(
        trim_baf_rawdata(realdata_bafs),
        trim_baf_rawdata(simul_bafs),
        ax=axd['qq'],
    )
    axd['qq'].set_xlabel('real')
    axd['qq'].set_ylabel('simulated')
    axd['qq'].set_title('QQ plot')

    return fig, axd


def validate_empirical_alphasum_new(
    empirical_alphasum, 
    mean_depth, 
    realdata_vafs,
    omega=EMPIRICAL_OMEGA,
):
    # prepare data
    trimmed_real_vafs = trim_vaf_rawdata(realdata_vafs)

    simul_vafs = simulate_variant_data(
        mean_depth=mean_depth,
        allele_portions=[1, 1], 
        size=int(len(trimmed_real_vafs) * 1.1), 

        omega=omega,
        alpha_sum=empirical_alphasum, 
        dirichlet=True,
        mode='gpois',
    )['vaf'][:, 0]
    trimmed_simul_vafs = trim_vaf_rawdata(simul_vafs)

    # make axd
    fig, axd = plt.subplot_mosaic(
        [
            ['hist_real', 'qq'],
            ['hist_simul', 'qq'],
        ],
        figsize=(30, 10),
        gridspec_kw=dict(
            #wspace=0.4,
        ),
    )

    # draw hist
    hist_bins = np.arange(0, 1.01, 0.01)
    _ = axd['hist_real'].hist(
        trimmed_real_vafs, 
        bins=hist_bins, 
        density=True,
    )
    axd['hist_real'].set_title('REAL')

    hist_ys, bins, patches = axd['hist_simul'].hist(
        trimmed_simul_vafs, 
        bins=hist_bins, 
        density=True,
    )
    axd['hist_simul'].set_title('SIMULATED')

    ymax = np.max(hist_ys) * 1.2
    for key in ['hist_real', 'hist_simul']:
        axd[key].set_ylim(0, ymax)

    # draw qqplot
    plotmisc.qqplot(
        trimmed_real_vafs,
        trimmed_simul_vafs,
        ax=axd['qq'],
    )
    axd['qq'].set_xlabel('REAL')
    axd['qq'].set_ylabel('SIMULATED')
    axd['qq'].set_title('QQ plot')

    return fig, axd


##########################################################
# making source data for constructing true baf estimator #
##########################################################

def make_stat_values_unitjob(true_baf, mean_depth, reps, size, alpha_sum, omega):
    stat_values = {key: list() for key in ('mean',)}
    for _ in range(reps):
        simul_vardata = simulate_variant_data(
            mean_depth=mean_depth, 
            allele_portions=[true_baf, 1 - true_baf], 
            size=size, 
            alpha_sum=alpha_sum, 
            omega=omega,
            dirichlet=True,
            mode='gpois',
            ALT_biases=None,
        )
        vafdata_stats = get_vafdata_stats(simul_vardata['vaf'][:, 0], trim=False, digitize=False)
        stat_values['mean'].append(vafdata_stats['mean'])

    return stat_values


def make_truebaf_estimation_data(
    alpha_sum=EMPIRICAL_ALPHASUM, 
    omega=EMPIRICAL_OMEGA,
    reps=100, 

    true_baf_cands=np.round(np.arange(0.01, 0.51, 0.01), 2),
    mean_depth_cands=np.arange(10, 110, 10),
    size_cands=np.array([10, 20, 50, 100, 200, 500, 1000, 2000]), 

    verbose=True,
    nproc=1,
):
    true_baf_grid, mean_depth_grid, size_grid = np.meshgrid(
        true_baf_cands, mean_depth_cands, size_cands,
        indexing='ij',
    )
    stat_data = {
        key: {
            'values': np.empty(
                (len(true_baf_cands), len(mean_depth_cands), len(size_cands), reps), 
                dtype=float,
            ),
            'mean': None,
            'std': None,
        }
        for key in ('mean',)
    }
    with multiprocessing.Pool(nproc) as pool:
        args = (
            (baf, depth, reps, size, alpha_sum, omega)
            for (baf, depth, size) in zip(
                true_baf_grid.flat, mean_depth_grid.flat, size_grid.flat
            )
        )
        mp_result = pool.starmap(make_stat_values_unitjob, args)

    for idx, stat_values in zip(np.ndindex(*true_baf_grid.shape), mp_result):
        for key in stat_values.keys():
            stat_data[key]['values'][idx + (...,)] = stat_values[key]

    for subdic in stat_data.values():
        subdic['mean'] = subdic['values'].mean(axis=-1)
        subdic['std'] = subdic['values'].std(axis=-1, ddof=0)

    estim_data = {
        'stat_data': stat_data, 

        'true_bafs': true_baf_grid, 
        'mean_depths': mean_depth_grid,
        'sizes': size_grid,

        'reps': reps,
    }

    return estim_data


def draw_truebaf_estimation(
    estim_data, 
    mean_depths=None, 
    sizes=None, 
    statkey='mean', 
    figsize=(30, 20),
    gridspec_kw=dict(),
):
    all_true_bafs = list(estim_data['true_bafs'][:, 0, 0])
    all_mean_depths = list(estim_data['mean_depths'][0, :, 0])
    all_sizes = list(estim_data['sizes'][0, 0, :])

    if mean_depths is None:
        mean_depths = all_mean_depths
    if sizes is None:
        sizes = all_sizes

    mean_depth_indexes = [all_mean_depths.index(x) for x in mean_depths]
    size_indexes = [all_sizes.index(x) for x in sizes]

    fig, axs = plt.subplots(
        len(mean_depths), len(sizes), 
        figsize=figsize,
        gridspec_kw=gridspec_kw,
    )
    #fig.supxlabel(f'{statkey} of baf rawdata', y=0)
    #fig.supylabel('true baf', x=0)
    for (idx_depth, idx_size), ax in zip(
        itertools.product(mean_depth_indexes, size_indexes),
        axs.flat,
    ):
        depth_value = all_mean_depths[idx_depth]
        size_value = all_sizes[idx_size]
        ax.set_title(
            f'mean_depth={depth_value}, size={size_value}', 
            size=15,
        )

        dataset = np.transpose(
            estim_data['stat_data'][statkey]['values'][:, idx_depth, idx_size, :]
        )
        ax.violinplot(
            dataset,
            positions=all_true_bafs, 
            widths=0.005, 
            showmeans=False, 
            quantiles=np.tile(
                np.array([0.25, 0.5, 0.75])[:, np.newaxis],
                [1, dataset.shape[1]],
            ),
        )
    
    return fig, ax


######################
# true baf estimator #
######################

def make_interp(estim_data, size=1000):
    size_idx = list(estim_data['sizes'][0, 0, :]).index(size)
    true_bafs = estim_data['true_bafs'][:, :, size_idx]
    mean_depths = estim_data['mean_depths'][:, :, size_idx]
    baf_means = estim_data['stat_data']['mean']['mean'][:, :, size_idx]
    
    interp = scipy.interpolate.LinearNDInterpolator(
        np.stack([baf_means.ravel(), mean_depths.ravel()], axis=1),
        true_bafs.ravel(),
    )
    return interp


def save_interp(interp):
    with open(INTERP_SAVEPATH, 'wb') as outfile:
        pickle.dump(interp, outfile)


@functools.cache
def load_interp():
    with open(INTERP_SAVEPATH, 'rb') as infile:
        return pickle.load(infile)


@functools.cache
def make_meanbaf_limits(interp):
    # sort by depth, ascending
    points_sorted = interp.points[
        np.argsort(interp.points[:, 1]),
        :,
    ]

    # split by depth values
    points_bydepth = np.split(
        points_sorted,
        (
            np.insert(
                np.diff(points_sorted[:, 1]), 0, 0
            )
            .nonzero()[0]
        ),
    )

    bafmean_minmax_bydepth = list()
    for x in points_bydepth:
        baf_means = x[:, 0]
        baf_means_min = baf_means.min()
        baf_means_max = baf_means.max()
        depth = x[0, 1]
        bafmean_minmax_bydepth.append(
            (depth, baf_means_min, baf_means_max)
        )
    bafmean_minmax_bydepth = np.stack(bafmean_minmax_bydepth, axis=0)
        # axis 1: depth, baf_mean min, baf_mean max

    return bafmean_minmax_bydepth


def replace_nan_args(nan_data, bafmean_minmax_bydepth):
    mean_bafs = nan_data[:, 0]
    depths = nan_data[:, 1]
    
    #depth_bins = bafmean_minmax_bydepth[:, 0]
    #depth_bin_indexes = np.digitize(depths, depth_bins)
    #depth_bin_indexes[depth_bin_indexes == 0] = 1
    #depth_bin_indexes[:] = depth_bin_indexes - 1

    depth_bins = bafmean_minmax_bydepth[:, 0]
    depth_diffs = np.abs(depths[:, np.newaxis] - depth_bins)
    depth_bin_indexes = depth_diffs.argmin(axis=1)

    new_depths = bafmean_minmax_bydepth[depth_bin_indexes, 0]
    # mean_baf_mins = bafmean_minmax_bydepth[depth_bin_indexes, 1]
    mean_baf_mins = bafmean_minmax_bydepth[:, 1].max()
    mean_baf_maxs = bafmean_minmax_bydepth[depth_bin_indexes, 2]
    new_mean_bafs = np.clip(mean_bafs, mean_baf_mins, mean_baf_maxs)

    return np.stack([new_mean_bafs, new_depths], axis=1)


def predict_true_baf(data, interp=None):
    """Args:
        data: array with shape (ndata, 2); rows indicate data index; column 0 is raw baf; column 1 is depth
    """
    assert data.ndim == 2

    if interp is None:
        interp = load_interp()
    bafmean_minmax_bydepth = make_meanbaf_limits(interp)

    input_nan_selector = np.isnan(data).any(axis=1)

    interp_result = interp(data)
    nan_selector = (
        np.isnan(interp_result)
        & (~input_nan_selector)
    )
    new_data = replace_nan_args(
        data[nan_selector, :], bafmean_minmax_bydepth,
    )
    new_results = interp(new_data)
    assert not np.isnan(new_results).any()

    interp_result[nan_selector] = new_results
    return interp_result



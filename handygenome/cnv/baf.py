import os
import pickle
import itertools

import numpy as np
import scipy.stats
import scipy.optimize
import scipy.interpolate
import sklearn.mixture

import handygenome
import handygenome.cnv.misc as cnvmisc


BAFCORRECT_PATH = os.path.join(handygenome.DIRS['data'], f'baf_correction_func.pickle')

def save_bafcorrect_func(func):
    with open(BAFCORRECT_PATH, 'wb') as outfile:
        pickle.dump(func, outfile)


def load_bafcorrect_func(x_cutoff=None):
    with open(BAFCORRECT_PATH, 'rb') as infile:
        func = pickle.load(infile)

    return winsorize_pointfive(func, x_cutoff=x_cutoff)

#

def save_bafcorrect_func_xcutoff(func, x_cutoff):
    with open(BAFCORRECT_PATH, 'wb') as outfile:
        pickle.dump((func, x_cutoff), outfile)

    
def load_bafcorrect_func_xcutoff():
    with open(BAFCORRECT_PATH, 'rb') as infile:
        func, x_cutoff = pickle.load(infile)

    return winsorize_pointfive(func, x_cutoff)


def get_bafs(vaf_list):
    bafs = np.array(vaf_list)
    gt_half_indexes = bafs > 0.5
    bafs[gt_half_indexes] = 1 - bafs[gt_half_indexes]
    return bafs


# baf simulation

def make_vafs_bafs(baf, mean_depth, size):
    halfsize = int(size / 2)
    depths = scipy.stats.poisson.rvs(mu=mean_depth, size=size)
    alts = np.concatenate(
        (
            scipy.stats.binom.rvs(n=depths[:halfsize], p=baf), 
            scipy.stats.binom.rvs(n=depths[halfsize:], p=(1 - baf)), 
        )
    )

    vafs = alts / depths
    np.random.shuffle(vafs)

    bafs = vafs.copy()
    bafs[bafs > 0.5] = 1 - bafs[bafs > 0.5]

    return vafs, bafs


def make_simulated_germline_vafs(mean_depth, size, mode, a=None, b=None):
    depths = scipy.stats.poisson.rvs(mu=mean_depth, size=size)

    if mode == 'betabinom':
        alts = scipy.stats.betabinom.rvs(n=depths, a=a, b=b)
    elif mode == 'binom':
        alts = scipy.stats.binom.rvs(n=depths, p=0.5)

    vafs = alts / depths

    return vafs


def make_simulated_germline_vafs_normal(size, std=0.1):
    result = scipy.stats.norm.rvs(loc=0.5, scale=std, size=size)
    result[result < 0] = 0
    result[result > 1] = 1
    return result


def infer_baf_minimize(baf_values, density=None):
    if density is None:
        density = scipy.stats.gaussian_kde(baf_values)

    argmin = scipy.optimize.minimize_scalar(lambda x: -density(x)[0], bounds=(0, 0.5))
    return argmin.x


def infer_baf_mixture(vaf_values):
    vals = vaf_values.reshape(-1, 1)
    gm = sklearn.mixture.GaussianMixture(n_components=2).fit(vals)
    means = sorted(gm.means_.ravel())
    return (means[0] + (1 - means[1])) / 2


def make_simulation_data(mean_depth=34, size=50, reps=100, true_bafs=None):
    if true_bafs is None:
        true_bafs = np.round(np.arange(0.01, 0.51, 0.01), 2)

    skew_data = dict()
    ebafmean_data = dict()
    ebafmin_data = dict()
    ebafmix_data = dict()
    baf_data = dict()
    vaf_data = dict()

    for baf in true_bafs:
        print(baf)

        skew_values = list()
        ebafmean_values = list()
        ebafmin_values = list()
        ebafmix_values = list()
        vafs_list = list()
        bafs_list = list()

        for _ in range(reps):
            vafs, bafs = make_vafs_bafs(baf=baf, mean_depth=mean_depth, size=size)
            vafs_list.append(vafs)
            bafs_list.append(bafs)
            skew_values.append(scipy.stats.skew(bafs))
            ebafmean_values.append(bafs.mean())
            ebafmin_values.append(infer_baf_minimize(bafs))
            ebafmix_values.append(infer_baf_mixture(vafs))

        skew_data[baf] = np.array(skew_values)
        ebafmean_data[baf] = np.array(ebafmean_values)
        ebafmin_data[baf] = np.array(ebafmin_values)
        ebafmix_data[baf] = np.array(ebafmix_values)
        baf_data[baf] = np.array(bafs_list)
        vaf_data[baf] = np.array(vafs_list)

    return {
        'skew': skew_data,
        'mean': ebafmean_data,
        'minimize': ebafmin_data,
        'mixture': ebafmix_data,
        'baf': baf_data,
        'vaf': vaf_data,
    }


def winsorize_pointfive(func, x_cutoff=None):
    def newfunc(x):
        result = func(x)
        result[result > 0.5] = 0.5
        if x_cutoff is not None:
            result[x > x_cutoff] = 0.5
        return result

    return newfunc


def get_xs_ratio(data, realdata, size=10000, rep=10000):
    realdata_means = np.array([np.random.choice(realdata, size).mean() for x in range(rep)])
    xs_ratio = realdata_means.mean() / data['mean'][0.5].mean()
    return xs_ratio


def get_x_cutoff(realdata, size=100, rep=10000, q=0.05):
    realdata_means = np.array([np.random.choice(realdata, size).mean() for x in range(rep)])
    return np.quantile(realdata_means, q)


def xs_ratio_correction(xs_ratio, xs, ys):
    if xs_ratio is None:
        corrected_xs = None
        interp = scipy.interpolate.interp1d(xs, ys, fill_value='extrapolate')
    else:
        corrected_xs = xs * xs_ratio
        interp = scipy.interpolate.interp1d(corrected_xs, ys, fill_value='extrapolate')

    return interp, corrected_xs


def make_bafcorrector_type1(data, xs_ratio=None):
    xs = list()
    ys = list()
    for true_baf_val, ebaf_vals in data['mean'].items():
        xs.append(ebaf_vals.mean())
        ys.append(true_baf_val)

    xs = np.array(xs)
    ys = np.array(ys)

    interp, corrected_xs = xs_ratio_correction(xs_ratio, xs, ys)
    return interp, xs, ys, corrected_xs


def make_bafcorrector_type2(data, xs_ratio=None):
    xs = list()
    ys = list()
    for true_baf_val, ebaf_vals in data['mean'].items():
        xs.extend(ebaf_vals)
        ys.extend(itertools.repeat(true_baf_val, len(ebaf_vals)))

    xs = np.array(xs)
    ys = np.array(ys)

    interp, corrected_xs = xs_ratio_correction(xs_ratio, xs, ys)

    return interp, xs, ys, corrected_xs


def make_bafcorrector_from_realdata(data, realdata, mode='type1'):
    xs_ratio = get_xs_ratio(data, realdata, size=10000, rep=10000)
    x_cutoff = get_x_cutoff(realdata, size=100, rep=10000, q=0.05)
    if mode == 'type1':
        interp, xs, ys, corrected_xs = make_bafcorrector_type1(data, xs_ratio=xs_ratio)
    elif mode == 'type2':
        interp, xs, ys, corrected_xs = make_bafcorrector_type2(data, xs_ratio=xs_ratio)

    bafcorrector = winsorize_pointfive(interp, x_cutoff)

    return bafcorrector, interp, x_cutoff

        
# baf simulation v2 (230525)

def make_simulated_total_depths(mean_depth, N):
    result = list()
    n = N
    while True:
        arr = scipy.stats.poisson.rvs(mean_depth, size=n)
        arr = arr[arr > 0]
        result.extend(arr)
        if len(result) == N:
            break
        else:
            n = N - len(result)
            continue
    return np.asarray(result)


def make_simulated_het_bafs(vaf, mean_depth, N):
    td = make_simulated_total_depths(mean_depth, N)
    ad = scipy.stats.binom.rvs(n=td, p=vaf)
    vafs = ad / td
    bafs = np.where(vafs > 0.5, 1 - vafs, vafs)
    return bafs


def make_simulated_hom_bafs(p_error, mean_depth, N):
    td = make_simulated_total_depths(mean_depth, N)
    ad = scipy.stats.binom.rvs(n=td, p=p_error)
    vafs = ad / td
    bafs = np.where(vafs > 0.5, 1 - vafs, vafs)
    return bafs


def make_simulated_germline_call_bafs(mean_depth, vaf, N, p_error, hom_portion):
    N_hom = int(N * hom_portion)
    N_het = N - N_hom
    return np.concatenate([
        make_simulated_het_bafs(vaf, mean_depth, N_het),
        make_simulated_hom_bafs(p_error, mean_depth, N_hom),
    ])


def infer_baf_density(bafs, bw, rmzero=True):
    if rmzero:
        bafs = bafs[bafs > 0]
    assert len(bafs) > 0
    if len(bafs) == 1:
        return bafs[0]
    else:
        peak_values, peak_densities, density = cnvmisc.get_density_peaks(bafs, bw_method=bw)
        #return np.average(peak_values, weights=peak_densities)
        if peak_values is None:
            return np.nan
        else:
            return max(peak_values)


def infer_baf_mean(bafs):
    bafs = bafs[bafs > 0]
    return bafs.mean()


def make_simulated_baf_dataset(mean_depth=30, N=10000, p_error=0.001, reps=10, hom_portion=0.6, verbose=False):
    true_bafs = np.arange(0, 0.51, 0.01)
    simbaf_data = list()
    for bafval in true_bafs:
        if verbose:
            print(f'baf value {bafval}')
        subdata = list()
        for _ in range(reps):
            #print(f'repetition {_}')
            simulated_bafs = make_simulated_germline_call_bafs(
                mean_depth, bafval, N, p_error, hom_portion,
            )
            subdata.append(simulated_bafs)
        simbaf_data.append(np.array(subdata))

    simbaf_data = np.stack(simbaf_data, axis=0)

    return simbaf_data


def make_inferred_baf_dataset(simbaf_data, bw=1, verbose=False, infer_method='density'):
    inferred_baf_data = list()
    for truebaf_group in simbaf_data:
        subdata = list()
        for bafvals in truebaf_group:
            if infer_method == 'density':
                inferred_baf = infer_baf_density(bafvals, bw=bw)
            elif infer_method == 'mean':
                inferred_baf = infer_baf_mean(bafvals)
            subdata.append(inferred_baf)
        inferred_baf_data.append(np.array(subdata))
    return np.stack(inferred_baf_data, axis=0)


def find_best_bandwidth(mean_depth):
    true_bafs = np.arange(0, 0.51, 0.01)
    simbaf_data = make_simulated_baf_dataset(
        mean_depth=mean_depth, N=1000, p_error=0.001, reps=10, hom_portion=0.6, verbose=False
    )
    def target(bw):
        inferred_baf_data = make_inferred_baf_dataset(simbaf_data, bw=bw, verbose=False, infer_method='density')
        stdsum = np.std(inferred_baf_data, axis=1).sum()
        means = np.mean(inferred_baf_data, axis=1)
        devsum = np.sqrt(np.sum((true_bafs - means) ** 2))
        return stdsum + devsum

    bwlist = np.arange(0.5, 1, 0.05)
    targetvals = [target(x) for x in bwlist]
    return bwlist[np.argmin(targetvals)]



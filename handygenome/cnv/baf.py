import os
import pickle
import itertools

import numpy as np
import scipy.stats
import scipy.optimize
import scipy.interpolate
import sklearn.mixture

import handygenome.common as common


BAFCORRECT_PATH = os.path.join(common.DATA_DIR, f'baf_correction_func.pickle')

def save_bafcorrect_func(func):
    with open(BAFCORRECT_PATH, 'wb') as outfile:
        pickle.dump(func, outfile)

    
def load_bafcorrect_func():
    with open(BAFCORRECT_PATH, 'rb') as infile:
        result = pickle.load(infile)

    return winsorize_pointfive(result)


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


def winsorize_pointfive(func):
    def newfunc(x):
        result = func(x)
        result[result > 0.5] = 0.5
        return result

    return newfunc


def make_bafcorrector_type1(data):
    xs = list()
    ys = list()
    for true_baf_val, ebaf_vals in data['mean'].items():
        xs.append(ebaf_vals.mean())
        ys.append(true_baf_val)

    interp = scipy.interpolate.interp1d(xs, ys, fill_value='extrapolate')
    return interp, xs, ys
    #return winsorize_pointfive(interp), xs, ys


def make_bafcorrector_type2(data):
    xs = list()
    ys = list()
    for true_baf_val, ebaf_vals in data['mean'].items():
        xs.extend(ebaf_vals)
        ys.extend(itertools.repeat(true_baf_val, len(ebaf_vals)))

    interp = scipy.interpolate.interp1d(xs, ys, fill_value='extrapolate')
    return interp, xs, ys
    #return winsorize_pointfive(interp), xs, ys


import pprint
import os
import re
import gzip
import shutil
import itertools

import scipy.stats
import numpy as np

import handygenome.deco as deco


PAT_INT = re.compile('-?[0-9]+')
PAT_FLOAT = re.compile('(-?[0-9]+\.[0-9]+)|(-?[0-9]+(\.[0-9]+)?e-?[0-9]+)')


# custom class base

def repr_base(obj, keylist, comma_sep_int=False):
    string_list = list()
    for key in keylist:
        val = getattr(obj, key)
        if isinstance(val, int) and comma_sep_int:
            string_list.append(f'{key}={val:,}')
        else:
            string_list.append(f'{key}={repr(val)}')

    return ', '.join(string_list)


# colors

COLORS = {
    'red':     '\033[38;5;196m',
    'magenta': '\033[38;5;201m',
    'pink':    '\033[38;5;213m',
    'orange':  '\033[38;5;9m',
    'yellow':  '\033[38;5;214m',
    'gold':    '\033[38;5;11m',
    'green':   '\033[38;5;40m',
    'blue':    '\033[38;5;33m',
    'cyan':    '\033[38;5;14m',
    'purple':  '\033[38;5;93m',
    'gray':    '\033[38;5;8m',
    'white':   '\033[38;5;15m',
    'end':     '\033[0m',
    }


class ColorsQQ:
    mine="\033[48;5;6m"
    busy="\033[48;5;244m"
    free="\033[48;5;238m"
    end="\033[0m"
    nor="\033[48;5;160m"
    nor2="\033[48;5;52m"
    danger1="\033[38;5;208m"
    danger2="\033[38;5;196m"
    darkgray="\033[38;5;240m"


def visualize_colors():
    for i in range(256):
        print(f'{i:<3d} \033[38;5;{i}m\\033[38;5;{i}m\033[0m')


def cpformat(obj, **kwargs):
    result = pprint.pformat(obj, **kwargs)
    result = re.sub('(True)', COLORS['green'] + '\\1' + COLORS['end'], result)
    result = re.sub('(False)', COLORS['red'] + '\\1' + COLORS['end'], result)
    result = re.sub('(None)', COLORS['purple'] + '\\1' + COLORS['end'], result)
    return result


def cpprint(obj):
    print(cpformat(obj))


# string manipulators

def str_to_nonstr(val):
    #assert isinstance(val, str)

    if val.lower() in ('none', 'null'):
        return None
    elif val.lower() == 'nan':
        return np.nan
    elif val.lower() == 'true':
        return True
    elif val.lower() == 'false':
        return False
    elif PAT_INT.fullmatch(val) is not None:
        return int(val)
    elif PAT_FLOAT.fullmatch(val) is not None:
        return float(val)
    else:
        return val


def nonstr_to_str(val):
    assert not isinstance(val, str)

    if val is None:
        return None
    else:
        return str(val)


def zrange(n):
    assert n > 0
    width = len(str(n - 1))
    for idx in range(n):
        zidx = str(idx).zfill(width)
        yield zidx


def zenumerate(iterable):
    iterable_tup = tuple(iterable)
    length = len(iterable_tup)
    for zidx, item in zip(zrange(length), iterable_tup):
        yield zidx, item


def zidx_to_idx(zidx):
    if re.match('^0+$', zidx):
        return 0
    else:
        return int(re.sub(r'^0*', '', zidx))


def get_padded_indices(n):
    """Begins with 0"""
    width = len(str(n-1))
    for idx in range(n):
        yield str(idx).zfill(width)
    #result = [str(idx).zfill(width) for idx in range(n)]
    #return result


def shorten_int(numlist, n_after_dot=3):
    mapping = {0: '', 1: 'K', 2: 'M', 3: 'G'}

    log = np.log10(numlist)
    assert (log < 12).all(), f'Numbers greater than or equal to 10^12 are not allowed.'

    #qs = [int(y) for y in (log / 3)]
    qs = np.floor(log / 3)

    suffixes = [mapping[x] for x in qs]
    #new_numlist = numlist / (10 ** (np.array(qs) * 3))
    new_numlist = numlist / (10 ** (qs * 3))
    formatter = f'.{n_after_dot}f'
    result = [
        (
            f'{int(x)}' 
            if y == '' else 
            f'{{x:{formatter}}} {{y}}'.format(x=x, y=y)
        )
        for x, y in zip(new_numlist, suffixes)
    ]
    return result


"""
http://code.activestate.com/recipes/578019
https://psutil.readthedocs.io/en/latest/#bytes-conversion
"""
_bytes2human_symbols = ('K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
_bytes2human_symbols_rev = tuple(reversed(_bytes2human_symbols))
_bytes2human_prefix = {
    s: 1 << (i + 1) * 10
    for i, s in enumerate(_bytes2human_symbols)
}

def bytes2human(n):
    for s in _bytes2human_symbols_rev:
        if abs(n) >= _bytes2human_prefix[s]:
            value = float(n) / _bytes2human_prefix[s]
            return '%.1f%s' % (value, s)
    return "%sB" % n


# listdir

def listdir(path):
    return sorted(os.path.join(path, x) for x in os.listdir(path))


# file handling

def unzip(src, dest, rm_src=False):
    with gzip.open(src, 'rb') as infile:
        with open(dest, 'wb') as outfile:
            shutil.copyfileobj(infile, outfile)      
    if rm_src:
        os.remove(src)


@deco.get_deco_arg_choices({'mode': ('r', 'w', 'a')})
def openfile(fname, mode='r'):
    mode = mode + 't'

    if fname.endswith('.gz'):
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)


def fileiter(path, sep='\t', remove_leading_hashes=True, skip_double_hashes=True):
    infile = openfile(path, 'r')

    while True:
        line = next(infile)
        if skip_double_hashes:
            if line.startswith('##'):
                continue
            else:
                break
        else:
            break
        
    headerline = line
    if remove_leading_hashes:
        headerline = re.sub('^#*', '', headerline)
    header = get_linesp(headerline, sep=sep)

    for line in infile:
        linesp = get_linesp(line, sep=sep)
        if len(linesp) != len(header):
            raise Exception(
                f'Field numbers of the header line and the current '
                f'line are different:\n'
                f'header: {header}\n'
                f'line: {linesp}')
        linedict = dict(zip(header, linesp))
        yield linedict

    infile.close()


def rm_newline(line):
    return re.sub('(\r)?\n$', '', line)


def rm_newline_byte(byteline):
    return re.sub(b'(\r)?\n$', b'', byteline)


def get_linesp(line, sep='\t'):
    return rm_newline(line).split(sep)


def get_linesp_byte(byteline, sep=b'\t'):
    return rm_newline_byte(byteline).split(sep)


# printwidth

def printwidth_get_width_list(df):
    '''
    df: [
    [line1_field1, line1_field2, ... ],
    [line2_field1, line2_field2, ... ],
    ...,
    ]
    '''
    width_list = list()
    for i in range(len(df[0])):
        width_list.append(list())

    for line in df:
        for idx, field in enumerate(line):
            width_list[idx].append(len(field))

    for idx, e in enumerate(width_list):
        width_list[idx] = max(e)

    return width_list


def printwidth_print_line(line, width_list, margin, target):
    printresult = ''
    for idx, e in enumerate(line):
        printresult += f'{e:>{width_list[idx] + margin}s}'
    if target == 'out':
        print(printresult, flush = True)
    elif target == 'err':
        print(printresult, flush = True, file = sys.stderr)


def printwidth(df, margin = 2, target = 'out'):
    for line in df:
        for idx, e in enumerate(line):
            line[idx] = str(e)

    width_list = printwidth_get_width_list(df)
    for line in df:
        printwidth_print_line(line, width_list, margin, target)


# memory usage check
def get_rss(mode='total', unit='b'):
    assert mode in ('self', 'children', 'total')
    assert unit in ('b', 'k', 'm', 'g')

    proc = psutil.Process()
    children = proc.children(recursive=True)
    rss_self = proc.memory_info().rss  # in bytes
    rss_children = sum(p.memory_info().rss for p in children)  # in bytes

    if mode == 'self':
        rss = rss_self
    elif mode == 'children':
        rss = rss_children
    elif mode == 'total':
        rss = rss_self + rss_children

    if unit == 'b':
        rss = rss
    elif unit == 'k':
        rss = rss / 1024
    elif unit == 'm':
        rss = rss / 1024**2
    elif unit == 'g':
        rss = rss / 1024**3

    return rss


# multi-min, multi-max

def _multi_selector_base(sequence, comparing_func, key=None, with_target_val=False):
    sequence = list(sequence)
    if key is None:
        values = sequence
    else:
        values = [key(x) for x in sequence]
    target_val = comparing_func(values)
    is_target_val = [x == target_val for x in values]

    if with_target_val:
        return list(itertools.compress(sequence, is_target_val)), target_val
    else:
        return list(itertools.compress(sequence, is_target_val))


def multi_min(sequence, key=None, with_target_val=False):
    return _multi_selector_base(sequence, min, key=key, with_target_val=with_target_val)


def multi_max(sequence, key=None, with_target_val=False):
    return _multi_selector_base(sequence, max, key=key, with_target_val=with_target_val)



# itertools plus utilities

def pairwise(iterable):
    iterable = iter(iterable)
    try:
        x0 = next(iterable)
        x1 = next(iterable)
    except StopIteration:
        return

    while True:
        yield (x0, x1)
        try:
            x2 = next(iterable)
        except StopIteration:
            return
        x0 = x1
        x1 = x2
        continue


def chunk_iter(iterable, n):
    args = [iter(iterable)] * n
    for subiter in itertools.zip_longest(*args, fillvalue=None):
        yield tuple(x for x in subiter if x is not None)


# from Itertools Recipes; similar as chunk_iter
def grouper_Itertools_Recipes(iterable, n, *, incomplete='fill', fillvalue=None):
    "Collect data into non-overlapping fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, fillvalue='x') --> ABC DEF Gxx
    # grouper('ABCDEFG', 3, incomplete='strict') --> ABC DEF ValueError
    # grouper('ABCDEFG', 3, incomplete='ignore') --> ABC DEF
    args = [iter(iterable)] * n
    if incomplete == 'fill':
        return itertools.zip_longest(*args, fillvalue=fillvalue)
    if incomplete == 'strict':
        return zip(*args, strict=True)
    if incomplete == 'ignore':
        return zip(*args)
    else:
        raise ValueError('Expected fill, strict, or ignore')



# array handlers & mathematical ones

def round_sig(num, n):
    """round 'num' with 'n' significant digits"""

    assert n > 0
    from math import log10, floor
    return round(num, -floor(log10(abs(num))) + (n-1))


def get_split_nums(total_length, num):
    arrsp = np.array_split(np.zeros(total_length), num)
    return np.fromiter(
        (x.shape[0] for x in arrsp if x.shape[0] != 0), 
        dtype=int,
    )


def array_grouper(arr, omit_values=False):
    assert arr.ndim in (1, 2)

    diff = np.empty(arr.shape[0], dtype=bool)
    diff[0] = True
    if arr.ndim == 1:
        diff[1:] = np.diff(arr)
    elif arr.ndim == 2:
        diff[1:] = np.diff(arr, axis=0).any(axis=1)

    indexes = np.nonzero(diff)[0]

    if omit_values:
        values = None
    else:
        values = arr[indexes]

    counts = np.empty(indexes.shape, dtype=int)
    counts[:-1] = np.diff(indexes)
    counts[-1] = arr.shape[0] - indexes[-1]

    groupkey = np.repeat(np.arange(len(counts)), counts)
        
    return values, counts, groupkey
                

def get_indexes_of_array(values, ordered_keys):
    """Used for handling pandas.cut function output, in cnv.misc module"""
    result = np.zeros(len(values))
    for idx, key in enumerate(ordered_keys):
        result[np.where(values == key)[0]] = idx

    return result


def nanaverage(values, weights):
    #assert isinstance(values, np.ndarray)
    #assert isinstance(weights, np.ndarray)

    values = np.array(values)
    weights = np.array(weights)

    selector = ~np.isnan(values)
    new_values = values[selector]
    new_weights = weights[selector]
    return np.average(new_values, weights=new_weights)


def get_ranks(arr):
    ranks = scipy.stats.rankdata(data, method='max')
    return ranks / len(ranks)


def get_mode(data, xs=np.arange(0, 0.51, 0.01)):
    kernel = scipy.stats.gaussian_kde(data)
    ys = kernel(xs)
    return xs[np.argmax(ys)]


def bernoulli_iterator(iterable, p, block_size=int(1e5)):
    """Args:
        p: probability of selection
    """
    for x, rv in zip(iter(iterable), bernoulli_rv_generator(p, block_size=block_size)):
        if rv:  # True if 1
            yield x


def bernoulli_rv_generator(p, block_size=int(1e5)):
    while True:
        rvs = scipy.stats.bernoulli.rvs(p=p, loc=0, size=block_size)
        for x in rvs:
            yield x


def mean_mad(values):
    values = np.array(values)
    return np.mean(np.abs(values - np.mean(values)))


def median_mad(values):
    values = np.array(values)
    return np.median(np.abs(values - np.median(values)))


def get_diffmean(values, weights=None):
    assert values.ndim == 1

    if weights is None:
        weights = np.ones_like(values)
    else:
        assert weights.shape == values.shape

    indexes = np.triu_indices(values.shape[0], k=1)
    diffs = values[indexes[0]] - values[indexes[1]]
    diff_weights = weights[indexes[0]] + weights[indexes[1]]

    return np.average(np.abs(diffs), weights=diff_weights)


##############################
# Genomic coordinate sorting #
##############################

def coord_sortkey(chrom, pos, chromdict): 
    """
    Args:
        pos: 1-based
        chromdict: ChromDict class instance
    """

    return (chromdict.contigs.index(chrom), pos)


#def get_vcfspec_sortkey(chromdict):
#    def sortkey(vcfspec):
#        return (chromdict.contigs.index(vcfspec.chrom), vcfspec.pos, vcfspec.ref) + vcfspec.alts
#        #return coord_sortkey(vcfspec.chrom, vcfspec.pos, chromdict)
#
#    return sortkey


def compare_coords(chrom1, pos1, chrom2, pos2, chromdict):
    """
    Returns:
        0: equal; -1: chrom1/pos1 comes first; 1: chrom2/pos2 comes first
    """

    if chrom1 == chrom2 and pos1 == pos2:
        return 0
    else:
        if (coord_sortkey(chrom1, pos1, chromdict) 
            < coord_sortkey(chrom2, pos2, chromdict)):
            return -1
        else:
            return 1



##########################################
# Check if file is plain text or gzipped #
##########################################

def check_file_is_plain(fname):
    with open(fname, 'r') as f:
        try:
            _ = f.read(1)
        except UnicodeDecodeError:
            return False
        else:
            return True


def check_file_is_gzipped(fname):
    with gzip.open(fname, 'rb') as f:
        try:
            _ = f.read(1)
        except gzip.BadGzipFile as exc:
            if str(exc).startswith('Not a gzipped file'):
                return False
            else:
                raise
        else:
            return True
    

##################################
# temporary file path generation #
##################################

def get_tmpfile_path(prefix=None, suffix=None, dir=None, where=None, delete=False, is_dir=False):
    """Args:
        where: alias for dir
    """
    # alias handling
    assert not ((dir is not None) and (where is not None))
    if where is not None:
        dir = where
    # main
    if dir is None:
        dir = os.getcwd()

    if delete:
        fd, path = tempfile.mkstemp(prefix=prefix, suffix=suffix, dir=dir)
        os.close(fd)
        os.remove(path)
    else:
        if is_dir:
            path = tempfile.mkdtemp(prefix=prefix, suffix=suffix, dir=dir)
        else:
            fd, path = tempfile.mkstemp(prefix=prefix, suffix=suffix, dir=dir)
            os.close(fd)

    return path

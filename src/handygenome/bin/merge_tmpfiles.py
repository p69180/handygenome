#!/opt/pydev/mamba/231124/mambaforge/envs/hdgn-edit/bin/python

import pprint
import os
import argparse
import gzip
import operator
import datetime


def argument_parsing():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'resultdir', 
        help=f'Output directory of a single storage watchdog run. It must contain a subdirectory "tmpfiles".',
    )
    args = parser.parse_args()

    return args


def parse_tmpfile(filepath, resultdict):
    with open(filepath, 'rt') as infile:
        try:
            colname = tuple(next(infile).rstrip().split('\t'))
        except StopIteration:
            return None
        except UnicodeDecodeError:
            return None
        else:
            try:
                for line in infile:
                    linedict = dict(zip(colname, line.rstrip().split('\t')))
                    if len(linedict) != len(colname):
                        continue
                    linedict['filesize'] = int(linedict['filesize'])
                    resultdict.setdefault(linedict['username'], list())
                    resultdict[linedict['username']].append(linedict)
            except UnicodeDecodeError:
                pass
            finally:
                return colname


def make_resultdict(tmpfile_dir):
    resultdict = dict()
    colname_candidates = set()
    for fname in os.listdir(tmpfile_dir):
        print(f'Reading from temporary file: {fname}', flush=True)
        filepath = os.path.join(tmpfile_dir, fname)
        colname = parse_tmpfile(filepath, resultdict)
        if colname is not None:
            if not colname_candidates:
                colname_candidates.add(colname)
            else:
                assert colname in colname_candidates, (
                    f'Different colname encountered: new={colname}, old={colname_candidates}'
                )

    return resultdict, colname_candidates.pop()


def write_result(resultdict, outdir, colname):
    for user, sublist in resultdict.items():
        outfile_path = os.path.join(outdir, f'{user}.tsv.gz')
        with gzip.open(outfile_path, 'wt') as outfile:
            outfile.write('\t'.join(colname) + '\n')
            for linedict in sorted(sublist, key=(lambda x: x['filesize']), reverse=True):
                outfile.write(
                    '\t'.join(str(linedict[key]) for key in colname) + '\n'
                )


def write_sizebyuser(resultdict, outdir):
    outfile_path = os.path.join(outdir, f'_SIZE_BY_USER.tsv.gz')
    assert not os.path.exists(outfile_path)

    sizebyuser = {
        key: sum(x['filesize'] for x in val) 
        for key, val in resultdict.items()
    }
    with gzip.open(outfile_path, 'wt') as outfile:
        for user, sizesum in sorted(sizebyuser.items(), key=operator.itemgetter(1), reverse=True):
            sizesum_gb = round(sizesum / 1024**3, 3)
            outfile.write(f'{user}\t{sizesum_gb} GB\n')


def main():
    args = argument_parsing()

    # parse tmpfiles
    tmpfile_dir = os.path.join(args.resultdir, 'tmpfiles')
    assert os.path.exists(tmpfile_dir)
    resultdict, colname = make_resultdict(tmpfile_dir)

    # write result
    interim_topdir = os.path.join(args.resultdir, 'interim_results')
    os.makedirs(interim_topdir, exist_ok=True)
    outdir = os.path.join(interim_topdir, datetime.datetime.now().isoformat())
    os.makedirs(outdir, exist_ok=True)
    print(f'Writing output to {outdir}')
    write_result(resultdict, outdir, colname)
    write_sizebyuser(resultdict, outdir)


if __name__ == '__main__':
    main()



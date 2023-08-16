import re
import itertools
import collections
import functools

import numpy as np
import pandas as pd
import pyranges as pr

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
readhandler = importlib.import_module(
    ".".join([top_package_name, "read", "readhandler"])
)
alignhandler = importlib.import_module(".".join([top_package_name, "align", "alignhandler"]))
librealign = importlib.import_module(".".join([top_package_name, "align", "realign"]))
fetchcache = importlib.import_module(".".join([top_package_name, "read", "fetchcache"]))


DEL_VALUE = "*"
EMPTY_VALUE = ""
DEFAULT_EXTEND_FETCHEDREADS_BY = 300
DEFAULT_EXTEND_PILEUP_BY = 30
#MAX_READ_LENGTH = 151
#
#
#def fetch_arghandler(query_start0, query_end0, db_start0, db_end0):
#    if query_start0 is None:
#        result_start0 = db_start0
#    else:
#        result_start0 = max(query_start0, db_start0)
#
#    if query_end0 is None:
#        result_end0 = db_end0
#    else:
#        result_end0 = min(query_end0, db_end0)
#
#    return result_start0, result_end0
#
#
#class FetchedReadsCluster:
#    def __init__(self, bam, chrom, readfilter=None):
#        self.bam = bam
#        self.chrom = chrom
#
#        if readfilter is None:
#            self.readfilter = readhandler.readfilter_pileup
#        else:
#            self.readfilter = readfilter
#
#        self.fetch_start0 = None
#        self.fetch_end0 = None
#        self.dict = collections.OrderedDict()
#
#    def __repr__(self):
#        return f'<FetchedReadsCluster object (chrom={self.chrom}, fetch_start0={self.fetch_start0:,}, fetch_end0={self.fetch_end0:,})>'
#
#    @classmethod
#    def from_fetch(cls, bam, chrom, start0, end0, readfilter=None):
#        result = cls(bam, chrom, readfilter)
#        result.fetch_start0 = start0
#        result.fetch_end0 = end0
#
#        for read in readhandler.get_fetch(
#            result.bam, result.chrom, start0, end0, readfilter=result.readfilter
#        ):
#            read_uid = readhandler.get_uid(read)
#            result.dict[read_uid] = read
#
#        return result
#
#    ###########
#    # yielder #
#    ###########
#    @staticmethod
#    def yielder_with_uid(uid, read):
#        return uid, read
#
#    @staticmethod
#    def yielder_without_uid(uid, read):
#        return read
#
#    def get_yielder(self, with_uid):
#        if with_uid:
#            return self.yielder_with_uid
#        else:
#            return self.yielder_without_uid
#    ###########
#
#    @property
#    def fetch_range0(self):
#        return range(self.fetch_start0, self.fetch_end0)
#
#    @property
#    def readspan_start0(self):
#        if len(self.dict) == 0:
#            return None
#        else:
#            return min(read.reference_start for read in self.dict.values())
#
#    @property
#    def readspan_end0(self):
#        if len(self.dict) == 0:
#            return None
#        else:
#            return max(read.reference_end for read in self.dict.values())
#
##    def _check_uid_uniqueness(self):
##        if len(set(self.uids)) != len(self.uids):
##            reads_string = "\n".join(x.to_string() for x in self.list)
##            uids_string = "\n".join(x for x in self.uids)
##            raise Exception(
##                f"Duplicate read uids.\n"
##                f"Read list:\n{reads_string}.\n"
##                f"uid list:\n{uids_string}."
##            )
#
#    def get_id(self):
#        return (self.chrom, self.fetch_start0, self.fetch_end0)
#
#    def get_read(self, read_uid):
#        return self.dict[read_uid]
#
#    def copy(self):
#        result = self.__class__(self.bam, self.chrom)
#        result.fetch_start0 = self.fetch_start0
#        result.fetch_end0 = self.fetch_end0
#        result.dict.update(self.dict)
#        return result
#
#    # merge, extend #
#    def merge(self, other):
#        if (other.fetch_start0 > self.fetch_end0) or (self.fetch_start0 > other.fetch_end0):
#            raise Exception(f'Fetch ranges of two objects do not overlap.')
#
#        if self.fetch_start0 <= other.fetch_start0 and self.fetch_end0 >= other.fetch_end0:
#            return self.copy()
#        elif other.fetch_start0 <= self.fetch_start0 and other.fetch_end0 >= self.fetch_end0:
#            return other.copy()
#        else:
#            # now self.fetch_start0 != other.fetch_start0
#            result = self.__class__(self.bam, self.chrom)
#            if self.fetch_start0 < other.fetch_start0:
#                left = self
#                right = other
#                result.fetch_start0 = self.fetch_start0
#                result.fetch_end0 = other.fetch_end0
#            else:
#                left = other
#                right = self
#                result.fetch_start0 = other.fetch_start0
#                result.fetch_end0 = self.fetch_end0
#
#            left_only_width = right.fetch_start0 - left.fetch_start0
#            right_only_width = right.fetch_end0 - left.fetch_end0
#            if left_only_width >= right_only_width:
#                result.dict.update(left.dict)
#
#                buffer = list(
#                    itertools.takewhile(
#                        (lambda x: x[0] not in left.dict),
#                        reversed(right.dict.items())
#                    )
#                )
#                result.dict.update(reversed(buffer))
#            else:
#                result.dict.update(
#                    itertools.takewhile(
#                        (lambda x: x[0] not in right.dict),
#                        left.dict.items()
#                    )
#                )
#                result.dict.update(right.dict)
#
#            return result
#
#    def extend_rightward(self, width):
#        for read, uid in itertools.dropwhile(
#            (lambda x: x[1] in self.dict),
#            readhandler.get_fetch(self.bam, self.chrom, self.fetch_end0, self.fetch_end0 + width, readfilter=self.readfilter, with_uid=True),
#        ):
#            self.dict[uid] = read
#
#        self.fetch_end0 += width
#
#    def extend_leftward(self, width):
#        buffer = list(
#            itertools.takewhile(
#                (lambda x: x[1] not in self.dict),
#                readhandler.get_fetch(self.bam, self.chrom, self.fetch_start0 - width, self.fetch_start0, readfilter=self.readfilter, with_uid=True),
#            )
#        )
#        for read, uid in reversed(buffer):
#            self.dict[uid] = read
#            self.dict.move_to_end(uid, last=False)
#
#        self.fetch_start0 -= width
#            
#    # fetch #
#    def fetch(self, start0=None, end0=None, with_uid=False):
#        if len(self.dict) == 0:
#            return iter(())
#        else:
#            start0, end0 = fetch_arghandler(start0, end0, self.fetch_start0, self.fetch_end0)
#            if start0 >= end0:
#                return iter(())
#            else:
#                yielder = self.get_yielder(with_uid)
#
#                for uid, read in self.dict.items():
#                    if read.reference_start >= end0:
#                        break
#                    if read.reference_end > start0:
#                        yield yielder(uid, read)
#
#    def _fetch_forward(self, fetch_start0, fetch_end0, with_uid=False):
#        yielder = self.get_yielder(with_uid)
#            
#        is_area1 = False
#        is_area2 = False
#        for uid, read in self.dict.items():
#            # pre-area1
#            if not is_area1:
#                # check if within area1
#                if read.reference_start > fetch_start0 - MAX_READ_LENGTH:
#                    is_area1 = True
#                else:
#                    # job pre-area1 (null)
#                    # continue
#                    continue
#            # area1
#            if not is_area2:
#                # check if within area2
#                if read.reference_start >= fetch_start0:
#                    is_area2 = True
#                else:
#                    # job area1
#                    if read.reference_end > fetch_start0:
#                        yield yielder(uid, read)
#                    # continue
#                    continue
#            # check if beyond area2
#            if read.reference_start >= fetch_end0:
#                break
#            # job area2
#            yield yielder(uid, read)
#                
#    def _fetch_reverse(self, fetch_start0, fetch_end0, with_uid=False):
#        yielder = self.get_yielder(with_uid)
#            
#        is_area1 = False
#        is_area2 = False
#        for uid, read in reversed(self.dict.items()):
#        #for uid, read in zip(reversed(self.uids), reversed(self.list)):
#            # post-area2
#            if not is_area2:
#                # check if within area2
#                if read.reference_start < fetch_end0:
#                    is_area2 = True
#                else:
#                    # job post-area2 (null)
#                    # continue
#                    continue
#            # area2
#            if not is_area1:
#                # check if within area1
#                if read.reference_start < fetch_start0:
#                    is_area1 = True
#                else:
#                    # job area2
#                    yield yielder(uid, read)
#                    # continue
#                    continue
#            # check if within pre-area1
#            if read.reference_start <= fetch_start0 - MAX_READ_LENGTH:
#                break
#            # job area1
#            if read.reference_end > fetch_start0:
#                yield yielder(uid, read)
#                                    
#    def fetch_choosing_side(self, start0=None, end0=None, with_uid=False):
#        if len(self.dict) == 0:
#            return iter(())
#        else:
#            start0, end0 = fetch_arghandler(start0, end0, self.fetch_start0, self.fetch_end0)
#            if start0 >= end0:
#                return iter(())
#            else:
#                # return self._fetch_forward(start0, end0, with_uid=with_uid)                
#                dist_from_start = start0 - self.fetch_start0
#                dist_from_end = self.fetch_end0 - end0
#
#                if dist_from_start <= dist_from_end:
#                    return self._fetch_forward(
#                        start0, end0, with_uid=with_uid
#                    )
#                else:
#                    return reversed(
#                        tuple(
#                            self._fetch_reverse(
#                                start0, end0, with_uid=with_uid
#                            )
#                        )
#                    )
#
#
#class FetchedReadsSinglechrom:
#    def __init__(self, chrom, bam, readfilter=None):
#        self.chrom = chrom
#        self.bam = bam
#
#        if readfilter is None:
#            self.readfilter = readhandler.readfilter_pileup
#        else:
#            self.readfilter = readfilter
#
#        self.clusters = list()
#
#    def __repr__(self):
#        return f'<FetchedReadsSinglechrom object (chrom={self.chrom}, clusters={self.clusters})>'
#
#    @classmethod
#    def from_fetch(cls, bam, chrom, start0, end0, readfilter=None):
#        result = cls(chrom, bam, readfilter)
#        result.clusters.append(FetchedReadsCluster.from_fetch(result.bam, result.chrom, start0, end0, result.readfilter))
#
#        return result
#
#    @property
#    def gr(self):
#        starts = list()
#        ends = list()
#        ids = list()
#        for x in self.clusters:
#            starts.append(x.fetch_start0)
#            ends.append(x.fetch_end0)
#            ids.append(x.get_id())
#        chroms = [self.chrom] * len(self.clusters)
#        return pr.from_dict(
#            {'Chromosome': chroms, 'Start': starts, 'End': ends, 'Id': ids}
#        )
#
##        return pr.PyRanges(
##            chromosomes=([self.chrom] * len(self.clusters)),
##            starts=[x.fetch_start0 for x in self.clusters],
##            ends=[x.fetch_end0 for x in self.clusters],
##        )
#
#    def get_clusters_by_ids(self, ids):
#        result = dict()
#        for x in self.clusters:
#            id = x.get_id()
#            if id in ids:
#                result[id] = x
#        return result
#
#    def get_read(self, read_uid):
#        found_read = False
#        for frcluster in self.clusters:
#            try:
#                read = frcluster.dict[read_uid]
#            except KeyError:
#                continue
#            else:
#                found_read = True
#                break
#
#        if not found_read:
#            raise Exception(f'Given ReadUID is not present.')
#
#        return read
#
#    def iter_reads(self):
#        return itertools.chain.from_iterable(x.dict.items() for x in self.clusters)
#
#    def sort_clusters(self):
#        self.clusters.sort(key=(lambda x: x.fetch_start0))
#
#    def merge(self):
#        new_clusters = list()
#        self.sort_clusters()
#
#        def unit_func(cluster1, cluster2):
#            if cluster1.fetch_end0 >= cluster2.fetch_start0:
#                return cluster1.merge(cluster2)
#            else:
#                new_clusters.append(cluster1)
#                return cluster2
#
#        final = functools.reduce(unit_func, self.clusters)
#        new_clusters.append(final)
#        self.clusters = new_clusters
#
#    def add_reads(self, start0, end0):
#        if len(self.clusters) == 0:
#            self.clusters.append(FetchedReadsCluster.from_fetch(self.bam, self.chrom, start0, end0, self.readfilter))
#        else:
#            self_gr = self.gr
#            input_gr = pr.PyRanges(chromosomes=[self.chrom], starts=[start0], ends=[end0])
#            subtracted_gr = input_gr.subtract(self_gr)
#                # When there is nothing left after subtraction, "subtracted_gr" becomes an empty PyRanges. 
#            nearest_gr = subtracted_gr.nearest(self_gr)
#                # "nearest_gr" becomes also an empty PyRanges. 
#                # Iteration with .df.iterrows() results in an empty iteration.
#            # get adjacent clusters
#            adj_frcluster_ids = list()
#            for idx, row in nearest_gr.df.iterrows():
#                if abs(row['Distance']) == 1:
#                    adj_frcluster_ids.append(row['Id'])
#            adj_frclusters = self.get_clusters_by_ids(adj_frcluster_ids)
#            # do read addition
#            for idx, row in nearest_gr.df.iterrows():
#                if abs(row['Distance']) == 1:
#                    frcluster = adj_frclusters[row['Id']]
#                    added_region_width = row['End'] - row['Start']
#                    if row['Start'] == row['End_b']:
#                        frcluster.extend_rightward(row['End'] - row['Start'])
#                    elif row['End'] == row['Start_b']:
#                        frcluster.extend_leftward(row['End'] - row['Start'])
#                else:
#                    self.clusters.append(FetchedReadsCluster.from_fetch(self.bam, self.chrom, row['Start'], row['End'], self.readfilter))
#
#            if len(adj_frclusters) > 0:
#                self.merge()
#
#    def fetch(self, start0, end0, with_uid=False):
#        input_gr = pr.PyRanges(chromosomes=[self.chrom], starts=[start0], ends=[end0])
#        joined_gr = input_gr.join(self.gr, how=None)
#
#        valid_frcluster_ids = list()
#        for idx, row in joined_gr.df.iterrows():
#            if row['Start_b'] <= row['Start'] and row['End_b'] >= row['End']:
#                valid_frcluster_ids.append(row['Id'])
#
#        if len(valid_frcluster_ids) == 0:
#            raise Exception(f'Input fetch range is not fully contained within database.')
#        if len(valid_frcluster_ids) > 1:
#            raise Exception(f'There are more than one FetchedReadsCluster objects containing input fetch range.')
#
#        frcluster_id = valid_frcluster_ids[0]
#        frcluster = self.get_clusters_by_ids([frcluster_id])[frcluster_id]
#        return frcluster.fetch(start0, end0, with_uid=with_uid)
#
#    def check_fully_includes(self, start0, end0):
#        input_gr = pr.PyRanges(chromosomes=[self.chrom], starts=[start0], ends=[end0])
#        intersect_gr = input_gr.intersect(self.gr, how='containment')
#        return not intersect_gr.empty
#
#    def choose_cluster(self, start0, end0):
#        self_gr = self.gr
#        input_gr = pr.PyRanges(chromosomes=[self.chrom], starts=[start0], ends=[end0])
#        joined_gr = input_gr.join(self_gr, how=None)
#        if joined_gr.df.shape[0] == 0:
#            return None
#        else:
#            overlapping_cluster_ids = {row['Id'] for idx, row in joined_gr.df.iterrows()}
#            overlapping_clusters = list()
#            for idx, row in self_gr.df.iterrows():
#                if row['Id'] in overlapping_cluster_ids:
#                    overlapping_clusters.append(self.clusters[idx])
#
#            if len(overlapping_clusters) == 1:
#                return overlapping_clusters[0]
#            else:
#                raise Exception(f'More than one clusters overlap with input range. Overlapping clusters: {overlapping_clusters}')
#        
#
#class FetchedReads:
#    def __init__(self, bam, readfilter=None):
#        self.bam = bam
#
#        if readfilter is None:
#            self.readfilter = readhandler.readfilter_pileup
#        else:
#            self.readfilter = readfilter
#
#        self.subdata = dict()
#
#    def __repr__(self):
#        buffer = list()
#        buffer.append(f'<FetchedReads object, which consists of:')
#        for chrom, frsinglechrom in self.subdata.items():
#            buffer.append(f'\t{frsinglechrom}')
#        buffer.append(f'>')
#        return '\n'.join(buffer)
#
#    @classmethod
#    def from_fetch(cls, bam, chrom, start0, end0, readfilter=None):
#        result = cls(bam, readfilter)
#        result.subdata[chrom] = FetchedReadsSinglechrom.from_fetch(result.bam, chrom, start0, end0, result.readfilter)
#
#        return result
#
#    def get_read(self, read_uid):
#        found_read = False
#        for frcluster in itertools.chain.from_iterable(
#            (cluster for cluster in frsinglechr.clusters)
#            for frsinglechr in self.subdata.values()
#        ):
#            try:
#                read = frcluster.dict[read_uid]
#            except KeyError:
#                continue
#            else:
#                found_read = True
#                break
#
#        if not found_read:
#            raise Exception(f'Given ReadUID is not present.')
#
#        return read
#
#    def iter_reads(self, chrom=None):
#        if chrom is None:
#            return itertools.chain.from_iterable(x.iter_reads() for x in self.subdata.values())
#        else:
#            return self.subdata[chrom].iter_reads()
#
#    def add_reads(self, chrom, start0, end0):
#        if chrom not in self.subdata:
#            self.subdata[chrom] = FetchedReadsSinglechrom(chrom, self.bam, self.readfilter)
#        self.subdata[chrom].add_reads(start0, end0)
#
#    def fetch(self, chrom, start0, end0, with_uid=False):
#        if chrom not in self.subdata:
#            raise Exception(f'Input chrom is not included in the database.')
#        return self.subdata[chrom].fetch(start0, end0, with_uid=with_uid)
#
#    def check_fully_includes(self, chrom, start0, end0):
#        if chrom not in self.subdata:
#            return False
#        else:
#            return self.subdata[chrom].check_fully_includes(start0, end0)
#
#    def choose_cluster(self, chrom, start0, end0):
#        if chrom not in self.subdata:
#            return False
#        else:
#            return self.subdata[chrom].choose_cluster(start0, end0)


########################


class PileupBase:
    def _coord_arg_sanitycheck(self, start0, end0):
        #if start0 is not None:
        if start0 < self.start0:
            raise Exception('Input "start0" argument is out of pileup range.')
        #if end0 is not None:
        if end0 > self.end0:
            raise Exception('Input "end0" argument is out of pileup range.')

    def _coord_arg_sanitycheck_pos0(self, pos0):
        if pos0 not in self.df.columns:
            raise Exception('Input "pos0" argument is out of pileup range.')

    @property
    def start0(self):
        return self.df.columns[0]

    @property
    def end0(self):
        return self.df.columns[-1] + 1

    @property
    def range0(self):
        return range(self.start0, self.end0)

    def get_ref_seq(self, start0=None, end0=None):
        if start0 is None:
            start0 = self.start0
        if end0 is None:
            end0 = self.end0
        self._coord_arg_sanitycheck(start0, end0)

        key = (start0, end0)
        if key not in self._ref_seq_cache:
            self._ref_seq_cache[key] = self.fasta.fetch(self.chrom, key[0], key[1])
        return self._ref_seq_cache[key]

#    @staticmethod
#    def get_max_vaf_from_pileupcol(col, with_counts=False, empty_value=EMPTY_VALUE):
#        counts = collections.Counter(x for x in col if x != empty_value)
#        total = sum(counts.values())
#        if total == 0:
#            raise Exception(f"Pileup column allele count sum is 0.")
#
#        max_allele, max_count = max(counts.items(), key=(lambda x: x[1]))
#        max_vaf = max_count / total
#
#        if with_counts:
#            return max_allele, max_vaf, counts
#        else:
#            return max_allele, max_vaf
#
#    def get_max_vaf_colindex(self, col_idx, with_counts=False):
#        col = self.df.iloc[:, col_idx]
#        return self.get_max_vaf_from_pileupcol(col, with_counts=with_counts)
#
#    def get_max_vaf_colpos0(self, col_pos0, with_counts=False):
#        col = self.df.loc[:, col_pos0]
#        return self.get_max_vaf_from_pileupcol(col, with_counts=with_counts)

    def get_allele_counter(self, pos0):
        #self._coord_arg_sanitycheck_pos0(pos0)
        col = self.df.loc[:, pos0]
        return collections.Counter(x for x in col if x != EMPTY_VALUE)

    def get_allele_portions(self, pos0):
        counter = self.get_allele_counter(pos0)
        counter_sum = sum(counter.values())
        portions = dict()
        for key, val in counter.items():
            portions[key] = val / counter_sum
        return portions

    # active info related ones #
    def iter_active_info_old_deprecated(self, start0, end0, reverse=False):
        self._coord_arg_sanitycheck(start0, end0)

        ref_seq = self.get_ref_seq(start0, end0)
        pos0_range = range(start0, end0)
        if reverse:
            ref_seq = ref_seq[::-1]
            pos0_range = pos0_range[::-1]

        for ref_base, pos0 in zip(ref_seq, pos0_range):
            max_allele, max_vaf = self.get_max_vaf_colpos0(pos0, with_counts=False)
            if max_vaf >= self.active_threshold:
                is_active = not (max_allele == ref_base)
            else:
                is_active = True
            yield is_active

    def iter_active_info(self, start0, end0, reverse=False):
        self._coord_arg_sanitycheck(start0, end0)

        ref_seq = self.get_ref_seq(start0, end0)
        pos0_range = range(start0, end0)
        if reverse:
            ref_seq = ref_seq[::-1]
            pos0_range = pos0_range[::-1]

        for ref_base, pos0 in zip(ref_seq, pos0_range):
            counter = self.get_allele_counter(pos0)
            non_ref_portion = 1 - (counter[ref_base] / sum(counter.values()))
            is_active = (non_ref_portion >= self.active_threshold)
            yield is_active 

    def _set_active_info(self, start0=None, end0=None):
        if start0 is None:
            start0 = self.start0
        if end0 is None:
            end0 = self.end0

        self._coord_arg_sanitycheck(start0, end0)
        if not hasattr(self, '_active_info'):
            setattr(self, '_active_info', list())

        active_info_idx_start = start0 - self.start0
        active_info_idx_end = end0 - self.start0

        self._active_info[
            active_info_idx_start:active_info_idx_end
        ] = self.iter_active_info(start0, end0, reverse=False)

    def get_active_info(self, start0=None, end0=None):
        if start0 is None:
            start0 = self.start0
        if end0 is None:
            end0 = self.end0

        self._coord_arg_sanitycheck(start0, end0)

        active_info_idx_start = start0 - self.start0
        active_info_idx_end = end0 - self.start0

        return pd.Series(
            self._active_info[active_info_idx_start:active_info_idx_end], 
            index=self.df.columns[active_info_idx_start:active_info_idx_end],
        )

    def get_active_positions(self):
        return tuple(
            itertools.compress(self.df.columns, self._active_info)
        )

    def check_inactive_margin_right(self, length):
        return not any(self._active_info[-length:])

    def check_inactive_margin_left(self, length):
        return not any(self._active_info[:length])

    # row specs related ones
    @staticmethod
    def row_spec_matcher(query, target):
        if query['left_filled']:
            if query['right_filled']:
                if target['left_filled'] and target['right_filled']:
                    return query['seq'] in target['seq']
                else:
                    return False
            else:
                if target['left_filled']:
                    return query['seq'] == target['seq'][:len(query['seq'])]
                else:
                    return False
        else:
            if query['right_filled']:
                if target['right_filled']:
                    return query['seq'] == target['seq'][-len(query['seq']):]
                else:
                    return False
            else:
                return query['seq'] in target['seq']

    def set_row_spec_groups(self):
        """row_spec's are grouped by whether one is a substring of another."""
        self.row_spec_groups = collections.OrderedDict()
        for query_rowid, query in sorted(
            self.row_specs.items(), 
            key=(lambda x: x[1]['match_length']), 
            reverse=True,
        ):
            superseq_candidates = list()
            for target_rowid in self.row_spec_groups.keys():
                target = self.row_specs[target_rowid]
                if self.row_spec_matcher(query, target):
                    superseq_candidates.append(target_rowid)
                if len(superseq_candidates) >= 2:
                    break
                    
            if len(superseq_candidates) == 0:
                self.row_spec_groups[query_rowid] = {
                    'subseq_rowids': set(),
                    'subseq_hits': 0,
                }
                self.row_spec_groups[query_rowid]['subseq_rowids'].add(query_rowid)
                self.row_spec_groups[query_rowid]['subseq_hits'] += 1
            elif len(superseq_candidates) == 1:
                superseq_rowid = superseq_candidates[0]
                self.row_spec_groups[superseq_rowid]['subseq_rowids'].add(query_rowid)
                self.row_spec_groups[superseq_rowid]['subseq_hits'] += 1
            else:
                for superseq_rowid in superseq_candidates:
                    self.row_spec_groups[superseq_rowid]['subseq_rowids'].add(query_rowid)

    # for debugging
    def show_row_spec_alignments(self):
        subseq_hits_sum = sum(x['subseq_hits'] for x in self.row_spec_groups.values())
        for superseq_id, groupinfo in sorted(
            self.row_spec_groups.items(),
            key=(lambda x: x[1]['subseq_hits']),
            reverse=True,
        ):
            superseq_row_spec = self.row_specs[superseq_id]
            lstrip_query_gaps = not superseq_row_spec["left_filled"]
            rstrip_query_gaps = not superseq_row_spec["right_filled"]

            superseq_aln = self._alignment_cache[superseq_id]
            superseq_vcfspec_list = alignhandler.alignment_to_vcfspec(
                superseq_aln, 
                self.start0, 
                self.chrom, 
                self.fasta,
                lstrip_query_gaps=lstrip_query_gaps,
                rstrip_query_gaps=rstrip_query_gaps,
            )

            print('subseq hits:', groupinfo['subseq_hits'])
            print('subseq hits portion:', groupinfo['subseq_hits'] / subseq_hits_sum)
            print(superseq_row_spec)
            print(superseq_vcfspec_list)
            print(superseq_aln)
            print()


class Pileup(PileupBase):
    pat_insclip = re.compile("(\(.+\))?([^()]+)(\(.+\))?")
    pat_parenthesis = re.compile("[()]")
    row_spec_exluded_vals = {DEL_VALUE, EMPTY_VALUE}

    def __init__(self, df, chrom, fasta, fetchedreads, active_threshold=None):
        self.df = df
        self.chrom = chrom
        self.fasta = fasta
        self.fetchedreads = fetchedreads
        self.bam = self.fetchedreads.bam

        if active_threshold is None:
            self.active_threshold = librealign.DEFAULT_ACTIVE_THRESHOLD
        else:
            self.active_threshold = active_threshold

        self._ref_seq_cache = dict()
        self._alignment_cache = dict()

        self._set_active_info()

    def pick_frcluster(self):
        return self.fetchedreads.choose_cluster(self.chrom, self.start0, self.end0)

    def subset(self, start0, end0, inplace=False):
        self._coord_arg_sanitycheck(start0, end0)
        # subset dataframe
        new_df = self.df.loc[:, start0:(end0 - 1)]
        # remove out-of-range rows
        row_within_range = list()
        for row_id, row in new_df.iterrows():
            read = self.get_read(row_id)
            if read.reference_end <= start0 or read.reference_start >= end0:
                row_within_range.append(False)
            else:
                row_within_range.append(True)
        new_df = new_df.loc[row_within_range, :]
        # active_info
        active_info_idx_start = start0 - self.start0
        active_info_idx_end = end0 - self.start0
        new_active_info = self._active_info[active_info_idx_start:active_info_idx_end]

        # result
        if inplace:
            self.df = new_df
            self._active_info = new_active_info
        else:
            result = self.__class__(
                df=new_df, chrom=self.chrom, fasta=self.fasta,
                fetchedreads=self.fetchedreads, active_threshold=self.active_threshold,
            )
            result._active_info = new_active_info
            result._ref_seq_cache = self._ref_seq_cache
            return result

    def get_read(self, row_id):
        return self.fetchedreads.get_read(row_id)

    ### extend ###  
    def extend_rightward(self, width):
        # keep original starts and ends
        border_col_idx_left = self.df.shape[1] - 1
        # make new pileup
        new_pileup = self._extend_helper_make_new_pileup(
            chrom=self.chrom,
            start0=self.end0,
            end0=(self.end0 + width),
            fetchedreads=self.fetchedreads,
        )
        # join
        self.df = self.df.join(new_pileup.df, how="outer")
        self.df[self.df.isnull()] = EMPTY_VALUE  # turn NaN into EMPTY_VALUE
        self._active_info.extend(new_pileup._active_info)
        # handle insclips facing each other
        self._extend_helper_handle_facing_insclip(border_col_idx_left)

    def extend_leftward(self, width):
        # keep original starts and ends
        border_col_idx_left = width - 1
        # make new pileup
        new_pileup = self._extend_helper_make_new_pileup(
            chrom=self.chrom,
            start0=(self.start0 - width),
            end0=self.start0,
            fetchedreads=self.fetchedreads,
        )
        # join
        self.df = new_pileup.df.join(self.df, how="outer")
        self.df[self.df.isnull()] = EMPTY_VALUE  # turn NaN into EMPTY_VALUE
        self._active_info = new_pileup._active_info + self._active_info
        # handle insclips facing each other
        self._extend_helper_handle_facing_insclip(border_col_idx_left)

    def _extend_helper_handle_facing_insclip(self, border_col_idx_left):
        border_col_idx_right = border_col_idx_left + 1
        border_col_left = self.df.iloc[:, border_col_idx_left]
        border_col_right = self.df.iloc[:, border_col_idx_right]
        for row_idx, (val_left, val_right) in enumerate(
            zip(border_col_left, border_col_right)
        ):
            if val_right.startswith("(") and val_left.endswith(")"):
                mat_left = self.__class__.pat_insclip.fullmatch(val_left)
                mat_right = self.__class__.pat_insclip.fullmatch(val_right)
                if mat_left.group(3) != mat_right.group(1):
                    raise Exception(
                        f"Insclip seqs of adjacent entries are different.\n{self.df}"
                    )
                self.df.iloc[
                    row_idx, border_col_idx_right
                ] = self.__class__.pat_insclip.sub("\\2\\3", "", val_right)

    def _extend_helper_make_new_pileup(self, chrom, start0, end0, fetchedreads):
        return get_pileup(
            chrom, 
            start0, 
            end0,
            fetchedreads=fetchedreads,
            fasta=self.fasta,
            active_threshold=self.active_threshold,
            truncate=True,
            as_array=False,
            return_range=False,
        )

    ################################
    # realignment-related features #
    ################################

    # secure_inactive_padding
    def secure_inactive_padding_rightward(self, inactive_padding, extend_pileup_by=DEFAULT_EXTEND_PILEUP_BY):
        librealign.secure_inactive_padding_rightward(
            pileup=self,
            inactive_padding=inactive_padding, 
            extend_pileup_by=extend_pileup_by, 
            inplace=True,
        )

    def secure_inactive_padding_leftward(self, inactive_padding, extend_pileup_by=DEFAULT_EXTEND_PILEUP_BY):
        librealign.secure_inactive_padding_leftward(
            pileup=self, 
            inactive_padding=inactive_padding, 
            extend_pileup_by=extend_pileup_by, 
            inplace=True,
        )

    def secure_inactive_padding(self, inactive_padding, extend_pileup_by=DEFAULT_EXTEND_PILEUP_BY):
        self.secure_inactive_padding_rightward(inactive_padding, extend_pileup_by=extend_pileup_by)
        self.secure_inactive_padding_leftward(inactive_padding, extend_pileup_by=extend_pileup_by)

    ### row specs ###
    def get_row_spec(self, row_id):
        row = self.df.loc[row_id, :]
        return {
            "seq": "".join(
                self.__class__.pat_parenthesis.sub("", x) 
                for x in row if x not in self.__class__.row_spec_exluded_vals
            ),
            "match_length": sum(x != EMPTY_VALUE for x in row),
            "left_filled": (row.iloc[0] != EMPTY_VALUE),
            "right_filled": (row.iloc[-1] != EMPTY_VALUE),
            "id": row.name,
        }

    def set_row_specs(self):
        self.row_specs = dict()
        for row_id, row in self.df.iterrows():
            self.row_specs[row_id] = self.get_row_spec(row_id)

    # others
    def save_alignments(self, superseq_rowid_list, aligner=None):
        if aligner is None:
            aligner = librealign.ALIGNER_MAIN
        ref_seq = self.get_ref_seq()
        ref_seq_reversed = ref_seq[::-1]

        for row_id in superseq_rowid_list:
            self._alignment_cache[row_id] = librealign.get_alignment_from_rowid(row_id, self, ref_seq, ref_seq_reversed, aligner)

    def augment_range(self):
        active_region_pileup.set_row_specs()
        active_region_pileup.set_row_spec_groups()
        vcfspecs = active_region_pileup.get_vcfspecs()

    def get_realigned_reads(self, aligner=None):
        """Returns:
            dict (keys ReadUID, values pysam.AlignedSegment)
        """
        if aligner is None:
            aligner = librealign.ALIGNER_MAIN
        ref_seq = self.get_ref_seq()
        ref_seq_reversed = ref_seq[::-1]
        return librealign.get_realigned_reads(self, self.range0, ref_seq, ref_seq_reversed, aligner)

    def get_vcfspecs(self, allele_portion_threshold=None):
        if allele_portion_threshold is None:
            allele_portion_threshold = librealign.DEFAULT_VCFSPEC_EXTRACTION_THRESHOLD
        return librealign.get_vcfspecs_from_pileup(self, allele_portion_threshold=allele_portion_threshold)


class MultisamplePileup(PileupBase):
    def __init__(self, pileup_dict):
        """Args:
            pileup_dict: keys - sampleid; values - Pileup object
        """
        self.pileup_dict = pileup_dict

        #self.active_threshold = self._get_active_threshold()
        self._alignment_cache = dict()

        first_pileup = next(iter(self.pileup_dict.values()))
        self.chrom = first_pileup.chrom
        self.fasta = first_pileup.fasta
        self._ref_seq_cache = first_pileup._ref_seq_cache

    def _get_active_threshold(self):
        depth_sum = sum(x.df.shape[0] for x in self.pileup_dict.values())
        corrected_active_thresholds = (
            x.active_threshold * (x.df.shape[0] / depth_sum)
            for x in self.pileup_dict.values()
        )
        return min(corrected_active_thresholds)

    def _get_active_threshold_simple(self):
        return sum(x.active_threshold for x in self.pileup_dict.values()) / (len(self.pileup_dict) ** 2)

    def set_df(self):
        assert len(set(x.range0 for x in self.pileup_dict.values())) == 1, f'Genomic ranges of pileup objects are different.'
        self.df = pd.concat(
            {key: val.df for key, val in self.pileup_dict.items()},
            names=['SampleID', 'ReadUID'],
        )
        self._set_active_info()

    def get_read(self, row_id):
        sampleid, read_uid = row_id
        return self.pileup_dict[sampleid].get_read(read_uid)

    ### extend ###
    def extend_rightward(self, width, extend_fetchedreads_by=DEFAULT_EXTEND_FETCHEDREADS_BY):
        for pileup in self.pileup_dict.values():
            pileup.extend_rightward(width, extend_fetchedreads_by=extend_fetchedreads_by)

    def extend_leftward(self, width, extend_fetchedreads_by=DEFAULT_EXTEND_FETCHEDREADS_BY):
        for pileup in self.pileup_dict.values():
            pileup.extend_leftward(width, extend_fetchedreads_by=extend_fetchedreads_by)

    # equalize
    def equalize_margins_leftward(self):
        max_range_start0 = min(x.start0 for x in self.pileup_dict.values())
        for pileup in self.pileup_dict.values():
            if pileup.start0 > max_range_start0:
                pileup.extend_leftward(pileup.start0 - max_range_start0)

    def equalize_margins_rightward(self):
        max_range_end0 = max(x.end0 for x in self.pileup_dict.values())
        for pileup in self.pileup_dict.values():
            if pileup.end0 < max_range_end0:
                pileup.extend_rightward(max_range_end0 - pileup.end0)

    def equalize_margins(self):
        self.equalize_margins_leftward()
        self.equalize_margins_rightward()

    # secure inactive padding
    def secure_inactive_padding_leftward(self, inactive_padding, extend_pileup_by=DEFAULT_EXTEND_PILEUP_BY, extend_fetchedreads_by=DEFAULT_EXTEND_FETCHEDREADS_BY):
        while True:
            self.equalize_margins_leftward()
            # check inactive margins
            margin_insufficient = list()
            for sampleid, pileup in self.pileup_dict.items():
                if not pileup.check_inactive_margin_left(inactive_padding):
                    margin_insufficient.append(sampleid)
            # check whether to break
            if len(margin_insufficient) == 0:
                break
            # secure inactive margins
            for sampleid in margin_insufficient:
                self.pileup_dict[sampleid].secure_inactive_padding_leftward(inactive_padding, extend_pileup_by=extend_pileup_by, extend_fetchedreads_by=extend_fetchedreads_by)

    def secure_inactive_padding_rightward(self, inactive_padding, extend_pileup_by=DEFAULT_EXTEND_PILEUP_BY, extend_fetchedreads_by=DEFAULT_EXTEND_FETCHEDREADS_BY):
        while True:
            self.equalize_margins_rightward()
            # check inactive margins
            margin_insufficient = list()
            for sampleid, pileup in self.pileup_dict.items():
                if not pileup.check_inactive_margin_right(inactive_padding):
                    margin_insufficient.append(sampleid)
            # check whether to break
            if len(margin_insufficient) == 0:
                break
            # secure inactive margins
            for sampleid in margin_insufficient:
                self.pileup_dict[sampleid].secure_inactive_padding_rightward(inactive_padding, extend_pileup_by=extend_pileup_by, extend_fetchedreads_by=extend_fetchedreads_by)

    ### row specs ###
    def set_row_specs(self):
        for pileup in self.pileup_dict.values():
            pileup.set_row_specs()

        self.row_specs = dict()
        for idx, row in self.df.iterrows():
            sampleid, read_uid = idx
            self.row_specs[(sampleid, read_uid)] = self.pileup_dict[sampleid].row_specs[read_uid]

    def divide_row_spec_groups(self):
        # initialize
        row_spec_groups_by_sample = dict()
        superseq_rowids_sorted = [
            superseq_rowid for superseq_rowid, groupinfo in sorted(
                self.row_spec_groups.items(),
                key=(lambda x: x[1]['subseq_hits']),
                reverse=True,
            )
        ]
        for sampleid, pileup in self.pileup_dict.items():
            row_spec_groups_by_sample[sampleid] = collections.OrderedDict()
            for superseq_rowid in superseq_rowids_sorted:
                row_spec_groups_by_sample[sampleid][superseq_rowid] = {
                    #'subseq_rowids': set(),
                    'subseq_hits': 0,
                }
        # fill in values
        for superseq_rowid, groupinfo in self.row_spec_groups.items():
            for subseq_rowid in groupinfo['subseq_rowids']:
                subseq_sampleid, subseq_readuid = subseq_rowid
                #row_spec_groups_by_sample[subseq_sampleid][superseq_rowid]['subseq_rowids'].add(subseq_readuid)
                row_spec_groups_by_sample[subseq_sampleid][superseq_rowid]['subseq_hits'] += 1
        # assign to each single sample pileup
        for sampleid, groups in row_spec_groups_by_sample.items():
            self.pileup_dict[sampleid].row_spec_groups = groups

    # get realigned reads #
    def get_realigned_reads(self, aligner=None):
        """Returns:
            dict (keys sampleid, values dict (keys ReadUID, values pysam.AlignedSegment))
        """
        if aligner is None:
            aligner = librealign.ALIGNER_MAIN
        ref_seq = self.get_ref_seq()
        ref_seq_reversed = ref_seq[::-1]

        realigned_reads_allsamples = librealign.get_realigned_reads(self, self.range0, ref_seq, ref_seq_reversed, aligner)
        realigned_reads_multisample = dict()
        for row_id, read in realigned_reads_allsamples.items():
            sampleid, readuid = row_id
            realigned_reads_multisample.setdefault(sampleid, dict())
            realigned_reads_multisample[sampleid][readuid] = read

        return realigned_reads_multisample

    def get_vcfspecs(self, allele_portion_threshold=None):
        """Must be run after 'divide_row_spec_groups' and 'get_realigned_reads' 
            methods have been run.
            'divide_row_spec_groups' is needed for setting 'row_spec_groups' for
                each singlesample Pileup.
            'get_realigned_reads' is needed for initiating '_alignment_cache' attribute.
        """
        if allele_portion_threshold is None:
            allele_portion_threshold = librealign.DEFAULT_VCFSPEC_EXTRACTION_THRESHOLD
        vcfspecs_allsamples = librealign.get_vcfspecs_multisample_pileup(self, allele_portion_threshold)
        return vcfspecs_allsamples


@common.get_deco_num_set(('bam', 'fetchedreads'), 1)
def get_pileup(
    chrom,
    start0,
    end0,
    bam=None,
    fetchedreads=None,
    fasta=None,
    active_threshold=None,
    truncate=True,
    as_array=False,
    return_range=False,
    del_value=DEL_VALUE,
    empty_value=EMPTY_VALUE,
):
    def sanity_check(chrom, start0, end0, bam, fetchedreads, fasta, as_array):
        if (not as_array) and (fasta is None):
            raise Exception(f'If "as_array" is False, "fasta" must be set.')

#    def prepare_fetchedreads(bam, fetchedreads):
#        if fetchedreads is None:
#            fetchedreads = fetchcache.FetchedReads(bam, readfilter=readhandler.readfilter_pileup)
#        return fetchedreads

    def make_readlist(fetchedreads, chrom, start0, end0):
        readlist = list()
        start_list = list()
        end_list = list()
        uid_list = list()
        for uid, read in fetchedreads.fetch(chrom, start0, end0, with_uid=True):
            readlist.append(read)
            start_list.append(read.reference_start)
            end_list.append(read.reference_end)
            uid_list.append(uid)
        tmp_pileup_range = range(min(start_list), max(end_list))

        return readlist, tmp_pileup_range, uid_list

    def raise_cigarpattern_error(read):
        raise Exception(f"Unexpected cigar pattern:\n{read.to_string()}")

    def cigar_sanitycheck(read):
        """Assumes M.* or [IS]M.*"""
        error = False
        first_cigarop = read.cigartuples[0][0]
        if first_cigarop != 0:
            if first_cigarop not in (1, 4):
                error = True
            else:
                if read.cigartuples[1][0] != 0:
                    error = True
        if error:
            raise_cigarpattern_error(read)

    def write_read_to_array_row(read, arr, arr_row_idx, initial_arr_col_idx, del_value):
        def handle_leading_insclip(read):
            leading_insclip_len = 0
            cigartup_idx = 0
            for cigarop, cigarlen in alignhandler.iter_leading_queryonly(read.cigartuples):
                leading_insclip_len += cigarlen
                cigartup_idx += 1

            if leading_insclip_len == 0:
                leading_insseq = None
            else:
                leading_insseq = read.query_sequence[:leading_insclip_len]
            query_idx = leading_insclip_len

            return leading_insseq, cigartup_idx, query_idx

        def handle_queryonly_only_case(read, arr, arr_col_idx, arr_row_idx):
            arr[arr_row_idx, arr_col_idx] = f"({read.query_sequence})"

        def handle_targetonly_after_match(
            read,
            arr,
            query_idx,
            arr_col_idx,
            arr_row_idx,
            trailing_nonM_cigarlen_sum,
            del_value,
        ):
            arr[arr_row_idx, arr_col_idx] = read.query_sequence[
                query_idx
            ]  # query_idx indicates the last base of M
            arr[
                arr_row_idx,
                (arr_col_idx + 1) : (arr_col_idx + 1 + trailing_nonM_cigarlen_sum),
            ] = del_value

        def handle_queryonly_after_match(
            read,
            arr,
            query_idx,
            arr_col_idx,
            arr_row_idx,
            trailing_nonM_cigarlen_sum,
        ):
            last_match_base = read.query_sequence[query_idx]
            insseq = read.query_sequence[
                (query_idx + 1) : (query_idx + 1 + trailing_nonM_cigarlen_sum)
            ]
            arr[arr_row_idx, arr_col_idx] = f"{last_match_base}({insseq})"

        def handle_nonlast_match(
            read,
            arr,
            cigarlen,
            query_idx,
            arr_col_idx,
            arr_row_idx,
            cigartup_idx,
        ):
            # adds all match seqs except the last one
            # query_idx and arr_col_idx updated
            if cigarlen > 1:
                sl_query = slice(query_idx, query_idx + cigarlen - 1)
                sl_arr_col = slice(arr_col_idx, arr_col_idx + cigarlen - 1)
                arr[arr_row_idx, sl_arr_col] = tuple(read.query_sequence[sl_query])
                query_idx += cigarlen - 1
                arr_col_idx += cigarlen - 1
            # get all subsequent non-M cigar units
            # cigartup_idx updated
            trailing_nonMs = list()
            while True:
                cigartup_idx += 1
                if cigartup_idx >= len(read.cigartuples):
                    break
                else:
                    next_cigartup = read.cigartuples[cigartup_idx]
                    if next_cigartup[0] != 0:
                        trailing_nonMs.append(next_cigartup)
                    else:
                        break
            # add trailing non-M cigar units
            # query_idx and arr_col_idx updated
            trailing_nonM_cigarops = set(x[0] for x in trailing_nonMs)
            trailing_nonM_cigarlen_sum = sum(x[1] for x in trailing_nonMs)
            if trailing_nonM_cigarops.issubset(alignhandler.CIGAROPS_TARGETONLY):  # D, N
                handle_targetonly_after_match(
                    read,
                    arr,
                    query_idx,
                    arr_col_idx,
                    arr_row_idx,
                    trailing_nonM_cigarlen_sum,
                    del_value,
                )
                query_idx += 1
                arr_col_idx += trailing_nonM_cigarlen_sum + 1
            elif trailing_nonM_cigarops.issubset(alignhandler.CIGAROPS_QUERYONLY):  # I, S
                handle_queryonly_after_match(
                    read,
                    arr,
                    query_idx,
                    arr_col_idx,
                    arr_row_idx,
                    trailing_nonM_cigarlen_sum,
                )
                query_idx += trailing_nonM_cigarlen_sum + 1
                arr_col_idx += 1
            else:
                raise Exception(
                    f"Consecutive Non-M cigar units are composed of both target-only and query-only ones:\n{read.to_string()}"
                )

            return cigartup_idx, query_idx, arr_col_idx

        def handle_last_match(
            read,
            arr,
            cigarlen,
            query_idx,
            arr_col_idx,
            arr_row_idx,
        ):
            # adds all match seqs to array
            sl_query = slice(query_idx, query_idx + cigarlen)
            sl_arr_col = slice(arr_col_idx, arr_col_idx + cigarlen)
            arr[arr_row_idx, sl_arr_col] = tuple(read.query_sequence[sl_query])

        # main
        """1) First, leading query-only cigar units are collected.
        2) Then, cigartuple is iterated on. Initial cigartup_idx is at the 
            first cigar unit which is not a query-only one. This cigar unit
            is assumed to be M. (If not, raises an exception)
        3) What is done in a loop cycle:
            The cigar unit which cigartup_idx indicates (which should be M)
            and all subsequent non-M cigarop are treated at single loop cycle.
            The subsequent non-M cigar units are assumed to be composed of only
            query-only ones or target-only ones. (If not, raises an exception)
        """

        # set mutable arr_col_idx. initial one is kept separately.
        arr_col_idx = initial_arr_col_idx

        # handles leading query-only cigar units. cigartup_idx and query_idx are initialized
        leading_insseq, cigartup_idx, query_idx = handle_leading_insclip(read)

        if cigartup_idx == len(read.cigartuples):
            # cigartup is only composed of query-only operations (no M)
            arr[arr_row_idx, arr_col_idx] = "(" + read.query_sequence + ")"            
        else:

            
            if read.cigartuples[cigartup_idx][0] != 0:
                raise Exception(
                    f"The first cigar operation after stripping leading I and S is not M:\n{read.to_string()}"
                )

            # begins loop
            while True:
                if cigartup_idx > len(read.cigartuples) - 1:
                    raise Exception(
                        f'"cigartup_idx" became greater than "len(read.cigartuples) - 1" while looping.'
                    )

                cigarop, cigarlen = read.cigartuples[cigartup_idx]
                assert (
                    cigarop == 0
                ), f'The cigar unit indicated by "cigarop_idx" at the beginning of a loop cycle is not M.'

                if (
                    cigartup_idx == len(read.cigartuples) - 1
                ):  # current cigartup is the last one
                    handle_last_match(
                        read,
                        arr,
                        cigarlen,
                        query_idx,
                        arr_col_idx,
                        arr_row_idx,
                    )
                    
                    break
                else:
                    cigartup_idx, query_idx, arr_col_idx = handle_nonlast_match(
                        read,
                        arr,
                        cigarlen,
                        query_idx,
                        arr_col_idx,
                        arr_row_idx,
                        cigartup_idx,
                    )
                    if cigartup_idx == len(read.cigartuples):
                        break

            # add leading query-only seq
            if leading_insseq is not None:
                arr[arr_row_idx, initial_arr_col_idx] = (
                    f"({leading_insseq})" + arr[arr_row_idx, initial_arr_col_idx]
                )

    def truncate_array(arr, tmp_pileup_range, pileup_range):
        sl_start = tmp_pileup_range.index(pileup_range.start)
        if pileup_range.stop == tmp_pileup_range.stop:
            sl_end = None
        else:
            sl_end = tmp_pileup_range.index(pileup_range.stop)
        return arr[:, slice(sl_start, sl_end)]

    def make_pileup_from_array(
        chrom,
        arr,
        fetchedreads,
        uid_list,
        pileup_range,
        tmp_pileup_range,
        truncate,
        active_threshold,
        fasta,
    ):
        df_columns = pileup_range if truncate else tmp_pileup_range
        df = pd.DataFrame(arr, columns=list(df_columns), index=uid_list)
        pileup = Pileup(
            df=df,
            chrom=chrom,
            fasta=fasta,
            fetchedreads=fetchedreads,
            active_threshold=active_threshold,
        )
        return pileup

    # main
    # parameter handling
    sanity_check(chrom, start0, end0, bam, fetchedreads, fasta, as_array)
    if active_threshold is None:
        active_threshold = librealign.DEFAULT_ACTIVE_THRESHOLD
    if fetchedreads is None:
        fetchedreads = fetchcache.FetchedReads(bam, readfilter=readhandler.readfilter_pileup)
    # prepare pileup_range
    pileup_range = range(start0, end0)
    #pileup_range, fetchedreads = set_params(chrom, start0, end0, bam, fetchedreads)
    readlist, tmp_pileup_range, uid_list = make_readlist(fetchedreads, chrom, start0, end0)
    # create array
    arr = np.full(
        shape=(len(readlist), len(tmp_pileup_range)),
        fill_value=empty_value,
        dtype=object,
    )
    for arr_row_idx, read in enumerate(readlist):
        initial_arr_col_idx = read.reference_start - tmp_pileup_range.start
        try:
            write_read_to_array_row(read, arr, arr_row_idx, initial_arr_col_idx, del_value)
        except Exception as exc:
            try:
                cigar_sanitycheck(read)
            except Exception as exc_cigarpattern:
                raise exc_cigarpattern from exc
            else:
                raise exc
    # truncate array
    if truncate:
        arr = truncate_array(arr, tmp_pileup_range, pileup_range)
    # prepare results
    if as_array:
        if return_range:
            return (arr, pileup_range)
        else:
            return arr
    else:
        pileup = make_pileup_from_array(
            chrom,
            arr,
            fetchedreads,
            uid_list,
            pileup_range,
            tmp_pileup_range,
            truncate,
            active_threshold,
            fasta,
        )
        if return_range:
            return (pileup, pileup_range)
        else:
            return pileup


def get_pileup_multisample(
    chrom,
    start0,
    end0,
    bam_dict,
    fasta,
    active_threshold_onesample=None,
    del_value=DEL_VALUE,
    empty_value=EMPTY_VALUE,
):
    # parameter handling
    if active_threshold_onesample is None:
        active_threshold_onesample = librealign.DEFAULT_ACTIVE_THRESHOLD
    # initialize pileup_dict
    pileup_dict = dict()
    for sampleid, bam in bam_dict.items():
        pileup_dict[sampleid] = get_pileup(
            chrom,
            start0,
            end0,
            bam=bam,
            fasta=fasta,
            active_threshold=active_threshold_onesample,
            truncate=True,
            as_array=False,
            return_range=False,
            del_value=del_value,
            empty_value=empty_value,
        )
    # create MultisamplePileup object and postprocess
    mspileup = MultisamplePileup(pileup_dict)
    mspileup.set_df()

    return mspileup



# def get_pileup(bam, chrom, start0, end0, as_array=False, truncate=True, with_reads=False):
#    """Returns:
#        A numpy.ndarray object if "as_array" argument is True
#        A Pileup object if "as_array" argument is False
#    """
#
#    read_dict = {
#        readhandler.get_uid(x): x
#        for x in readhandler.get_fetch(bam, chrom, start0, end0, readfilter=readhandler.readfilter_pileup)
#    }
#    readlist = tuple(iter(read_dict.values()))
#    pileup, pileup_range = get_pileup_from_reads(readlist, as_array=as_array, return_range=True)
#
#    if truncate:
#        if as_array:
#            start0_idx = pileup_range.index(start0)
#            end0_idx = pileup_range.index(end0)
#            pileup = pileup[:, start0_idx:end0_idx]
#        else:
#            pileup.df = pileup.df.loc[:, list(range(start0, end0))]
#
#    if with_reads:
#        return pileup, read_dict
#    else:
#        return pileup


# def get_pileup_from_reads(readlist, as_array=False, pileup_range=None, return_range=False, del_char='*'):
#    def raise_err(read):
#        raise Exception(f"Unexpected cigar pattern:\n{read.to_string()}")
#
#    def cigar_sanitycheck(read):
#        """Assumes M.* or [IS]M.*"""
#        error = False
#        first_cigarop = read.cigartuples[0][0]
#        if first_cigarop != 0:
#            if first_cigarop not in (1, 4):
#                error = True
#            else:
#                if read.cigartuples[1][0] != 0:
#                    error = True
#        if error:
#            raise_err(read)
#
#    def get_pileup_range(readlist):
#        start0s = list()
#        end0s = list()
#        for read in readlist:
#            start0s.append(read.reference_start)
#            end0s.append(read.reference_end)
#        start0 = min(start0s)
#        end0 = max(end0s)
#
#        return range(start0, end0)
#
#    def check_overlaps_pileup_range(read, pileup_range):
#        if read.cigartuples[-1][0] in readhandler.CIGAROPS_QUERYONLY:
#            overlaps_upstream = (read.reference_end >= pileup_range.start)
#        else:
#            overlaps_upstream = (read.reference_end > pileup_range.start)
#
#        overlaps_downstream = (read.reference_start < pileup_range.stop)
#
#        return overlaps_upstream and overlaps_downstream
#
#    def handle_leading_insclip(read):
#        leading_insclip_len = 0
#        for cigartup_idx, (cigarop, cigarlen) in enumerate(read.cigartuples):
#            if cigarop in (1, 4):
#                leading_insclip_len += cigarlen
#            else:
#                break
#
#        if leading_insclip_len == 0:
#            leading_insseq = None
#        else:
#            leading_insseq = read.query_sequence[:leading_insclip_len]
#
#        query_idx = leading_insclip_len
#
#        return leading_insseq, cigartup_idx, query_idx
#
#    def handle_last_match(
#        read, arr, cigarlen, query_idx, arr_col_idx, arr_row_idx, cigartup_idx
#    ):
#        sl_query = slice(query_idx, query_idx + cigarlen)
#        seq_buffer = tuple(read.query_sequence[sl_query])
#        sl_arr_col = slice(arr_col_idx, arr_col_idx + cigarlen)
#        arr[arr_row_idx, sl_arr_col] = seq_buffer
#        cigartup_idx += 1
#
#        return cigartup_idx
#
#    def add_matches_before_last(
#        read, arr, cigarlen, query_idx, arr_col_idx, arr_row_idx
#    ):
#        sl_query = slice(query_idx, query_idx + cigarlen - 1)
#        seq = tuple(read.query_sequence[sl_query])
#        query_idx += cigarlen - 1
#
#        sl_arr_col = slice(arr_col_idx, arr_col_idx + cigarlen - 1)
#        arr[arr_row_idx, sl_arr_col] = seq
#        arr_col_idx += cigarlen - 1
#
#        return query_idx, arr_col_idx
#
#    def handle_del_after_match(
#        read, arr, query_idx, arr_col_idx, arr_row_idx, next_cigarlen
#    ):
#        seq = read.query_sequence[query_idx]
#        query_idx += 1
#        arr[arr_row_idx, arr_col_idx] = seq
#        arr_col_idx += 1
#
#        sl_arr_col = slice(arr_col_idx, arr_col_idx + next_cigarlen)
#        arr[arr_row_idx, sl_arr_col] = del_char
#        arr_col_idx += next_cigarlen
#
#        return query_idx, arr_col_idx
#
#    def handle_insclip_after_match(
#        read, arr, query_idx, arr_col_idx, arr_row_idx, next_cigarlen
#    ):
#        match_seq = read.query_sequence[query_idx]
#        query_idx += 1
#        insseq = read.query_sequence[query_idx : (query_idx + next_cigarlen)]
#        query_idx += next_cigarlen
#        seq = match_seq + f"({insseq})"
#
#        arr[arr_row_idx, arr_col_idx] = seq
#        arr_col_idx += 1
#
#        return query_idx, arr_col_idx
#
#    def add_seqs(read, arr, cigartup_idx, query_idx, arr_row_idx, arr_col_idx):
#        while True:
#            if cigartup_idx == len(read.cigartuples):
#                break
#
#            cigarop, cigarlen = read.cigartuples[cigartup_idx]
#            if cigarop != 0:
#                raise_err(read)
#
#            if cigartup_idx == len(read.cigartuples) - 1:
#                cigartup_idx = handle_last_match(
#                    read,
#                    arr,
#                    cigarlen,
#                    query_idx,
#                    arr_col_idx,
#                    arr_row_idx,
#                    cigartup_idx,
#                )
#                continue
#            else:
#                if cigarlen > 1:
#                    query_idx, arr_col_idx = add_matches_before_last(
#                        read, arr, cigarlen, query_idx, arr_col_idx, arr_row_idx
#                    )
#
#                next_cigarop, next_cigarlen = read.cigartuples[cigartup_idx + 1]
#                if next_cigarop == 2:  # deletion
#                    query_idx, arr_col_idx = handle_del_after_match(
#                        read, arr, query_idx, arr_col_idx, arr_row_idx, next_cigarlen
#                    )
#                elif next_cigarop in (1, 4):
#                    query_idx, arr_col_idx = handle_insclip_after_match(
#                        read, arr, query_idx, arr_col_idx, arr_row_idx, next_cigarlen
#                    )
#
#                cigartup_idx += 2
#                continue
#
#    def add_leading_insseq(read, leading_insseq, arr, arr_row_idx, arr_col_idx_init):
#        if leading_insseq is not None:
#            old_seq = arr[arr_row_idx, arr_col_idx_init]
#            new_seq = f"(l+{leading_insseq})" + old_seq
#            arr[arr_row_idx, arr_col_idx_init] = f"({leading_insseq})" + old_seq
#
#    # main
#    chrom = readlist[0].reference_name
#    if pileup_range is None:
#        pileup_range = get_pileup_range(readlist)
#
#    ncol = len(pileup_range)
#    nrow = len(readlist)
#    #arr = np.empty((nrow, ncol), dtype=object)
#    arr = np.full((nrow, ncol), None, dtype=object)
#
#    for arr_row_idx, read in enumerate(
#        read for read in readlist
#        if check_overlaps_pileup_range(read, pileup_range)
#    ):
#        try:
#            arr_col_idx_init = read.reference_start - pileup_range.start
#            arr_col_idx = arr_col_idx_init
#            leading_insseq, cigartup_idx, query_idx = handle_leading_insclip(read)
#            add_seqs(read, arr, cigartup_idx, query_idx, arr_row_idx, arr_col_idx)
#            add_leading_insseq(read, leading_insseq, arr, arr_row_idx, arr_col_idx_init)
#        except Exception as exc:
#            try:
#                cigar_sanitycheck(read)
#            except Exception as exc_cigarpattern:
#                raise exc_cigarpattern from exc
#            else:
#                raise exc
#
#    if as_array:
#        if return_range:
#            return (arr, pileup_range)
#        else:
#            return arr
#    else:
#        uids = [readhandler.get_uid(read) for read in readlist]
#        if len(uids) != len(set(uids)):
#            reads_as_string = '\n'.join(x.to_string() for x in readlist)
#            raise Exception(f'Duplicate read uid. Input reads are as following:\n{reads_as_string}')
#
#        df = pd.DataFrame(arr, columns=list(pileup_range), index=uids)
#        pileup = Pileup(df=df, uids=uids, chrom=chrom)
#        if return_range:
#            return (pileup, pileup_range)
#        else:
#            return pileup

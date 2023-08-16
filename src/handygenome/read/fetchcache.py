import itertools
import collections
import functools

import pyranges as pr

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
readhandler = importlib.import_module(
    ".".join([top_package_name, "read", "readhandler"])
)


MAX_READ_LENGTH = 151
DEFAULT_FETCHEDREADS_ADDITION_PADDING = 200
DEFAULT_FETCHEDREFSEQ_ADDITION_PADDING = 200


def cluster_fetch_arghandler(query_start0, query_end0, db_start0, db_end0):
    if (
        ((query_end0 is not None) and (query_end0 <= db_start0)) or
        ((query_start0 is not None) and (query_start0 >= db_end0))
    ):
        raise Exception(f'Input fetch range does not overlap with database range.')

    if query_start0 is None:
        result_start0 = db_start0
    else:
        result_start0 = max(query_start0, db_start0)

    if query_end0 is None:
        result_end0 = db_end0
    else:
        result_end0 = min(query_end0, db_end0)

    return result_start0, result_end0


class ClusterBase:
    def get_id(self):
        return (self.chrom, self.start0, self.end0)


class SinglechromBase:
    @property
    def gr(self):
        starts = list()
        ends = list()
        ids = list()
        for x in self.clusters:
            starts.append(x.start0)
            ends.append(x.end0)
            ids.append(x.get_id())
        chroms = [self.chrom] * len(self.clusters)
        return pr.from_dict(
            {'Chromosome': chroms, 'Start': starts, 'End': ends, 'Id': ids}
        )

    def get_clusters_by_ids(self, ids):
        result = dict()
        for x in self.clusters:
            id = x.get_id()
            if id in ids:
                result[id] = x
        return result

    def sort_clusters(self):
        self.clusters.sort(key=(lambda x: x.fetch_start0))

    def merge(self):
        new_clusters = list()
        self.sort_clusters()

        def unit_func(cluster1, cluster2):
            if cluster1.end0 >= cluster2.start0:
                return cluster1.merge(cluster2)
            else:
                new_clusters.append(cluster1)
                return cluster2

        final = functools.reduce(unit_func, self.clusters)
        new_clusters.append(final)
        self.clusters = new_clusters

    def add(self, start0, end0):
        if len(self.clusters) == 0:
            self.clusters.append(self.spawn_cluster_from_fetch(start0, end0))
        else:
            self_gr = self.gr
            input_gr = pr.PyRanges(chromosomes=[self.chrom], starts=[start0], ends=[end0])
            subtracted_gr = input_gr.subtract(self_gr)
                # When there is nothing left after subtraction, "subtracted_gr" becomes an empty PyRanges. 
            nearest_gr = subtracted_gr.nearest(self_gr)
                # "nearest_gr" becomes also an empty PyRanges. 
                # Iteration with .df.iterrows() results in an empty iteration.
            # get adjacent clusters
            adj_cluster_ids = list()
            for idx, row in nearest_gr.df.iterrows():
                if abs(row['Distance']) == 1:
                    adj_cluster_ids.append(row['Id'])
            adj_clusters = self.get_clusters_by_ids(adj_cluster_ids)
            # do read addition
            for idx, row in nearest_gr.df.iterrows():
                if abs(row['Distance']) == 1:
                    cluster = adj_clusters[row['Id']]
                    added_region_width = row['End'] - row['Start']
                    if row['Start'] == row['End_b']:
                        cluster.extend_rightward(row['End'] - row['Start'])
                    elif row['End'] == row['Start_b']:
                        cluster.extend_leftward(row['End'] - row['Start'])
                else:
                    self.clusters.append(self.spawn_cluster_from_fetch(row['Start'], row['End']))

            if len(adj_clusters) > 0:
                self.merge()

    def get_valid_cluster(self, start0, end0):
        """When query range is not contained within database, 'add' is automatically done."""
        input_gr = pr.PyRanges(chromosomes=[self.chrom], starts=[start0], ends=[end0])
        joined_gr = input_gr.join(self.gr, how=None, apply_strand_suffix=False)

        valid_cluster_ids = list()
        for idx, row in joined_gr.df.iterrows():
            if row['Start_b'] <= row['Start'] and row['End_b'] >= row['End']:
                valid_cluster_ids.append(row['Id'])

        if len(valid_cluster_ids) == 0:
            self.add(start0 - self.addition_padding, end0 + self.addition_padding)
            return self.get_valid_cluster(start0, end0)
        elif len(valid_cluster_ids) == 1:
            cluster_id = valid_cluster_ids[0]
            cluster = self.get_clusters_by_ids([cluster_id])[cluster_id]
            return cluster
        else:
            raise Exception(f'There are more than one clusters containing input fetch range.')


##############################################


class FetchedReadsCluster(ClusterBase):
    def __init__(self, bam, chrom, readfilter=None):
        self.bam = bam
        self.chrom = chrom

        if readfilter is None:
            self.readfilter = readhandler.readfilter_pileup
        else:
            self.readfilter = readfilter

        self.fetch_start0 = None
        self.fetch_end0 = None
        self.dict = collections.OrderedDict()

    def __repr__(self):
        return f'<FetchedReadsCluster object (chrom={self.chrom}, fetch_start0={self.fetch_start0:,}, fetch_end0={self.fetch_end0:,})>'

    @classmethod
    def from_fetch(cls, bam, chrom, start0, end0, readfilter=None):
        result = cls(bam, chrom, readfilter)
        result.fetch_start0 = start0
        result.fetch_end0 = end0

        for read in readhandler.get_fetch(
            result.bam, result.chrom, start0, end0, readfilter=result.readfilter
        ):
            read_uid = readhandler.get_uid(read)
            result.dict[read_uid] = read

        return result

    ###########
    # yielder #
    ###########
    @staticmethod
    def yielder_with_uid(uid, read):
        return uid, read

    @staticmethod
    def yielder_without_uid(uid, read):
        return read

    def get_yielder(self, with_uid):
        if with_uid:
            return self.yielder_with_uid
        else:
            return self.yielder_without_uid
    ###########

    @property
    def start0(self):
        return self.fetch_start0

    @property
    def end0(self):
        return self.fetch_end0

    @property
    def fetch_range0(self):
        return range(self.fetch_start0, self.fetch_end0)

    @property
    def readspan_start0(self):
        if len(self.dict) == 0:
            return None
        else:
            return min(read.reference_start for read in self.dict.values())

    @property
    def readspan_end0(self):
        if len(self.dict) == 0:
            return None
        else:
            return max(read.reference_end for read in self.dict.values())

#    def _check_uid_uniqueness(self):
#        if len(set(self.uids)) != len(self.uids):
#            reads_string = "\n".join(x.to_string() for x in self.list)
#            uids_string = "\n".join(x for x in self.uids)
#            raise Exception(
#                f"Duplicate read uids.\n"
#                f"Read list:\n{reads_string}.\n"
#                f"uid list:\n{uids_string}."
#            )

    def get_id(self):
        return (self.chrom, self.fetch_start0, self.fetch_end0)

    def get_read(self, read_uid):
        return self.dict[read_uid]

    def copy(self):
        result = self.__class__(self.bam, self.chrom)
        result.fetch_start0 = self.fetch_start0
        result.fetch_end0 = self.fetch_end0
        result.dict.update(self.dict)
        return result

    # merge, extend #
    def merge(self, other):
        if (other.fetch_start0 > self.fetch_end0) or (self.fetch_start0 > other.fetch_end0):
            raise Exception(f'Fetch ranges of two objects do not overlap.')

        if self.fetch_start0 <= other.fetch_start0 and self.fetch_end0 >= other.fetch_end0:
            return self.copy()
        elif other.fetch_start0 <= self.fetch_start0 and other.fetch_end0 >= self.fetch_end0:
            return other.copy()
        else:
            # now self.fetch_start0 != other.fetch_start0
            if self.fetch_start0 < other.fetch_start0:
                left = self
                right = other
            else:
                left = other
                right = self

            result = self.__class__(self.bam, self.chrom)
            result.fetch_start0 = left.fetch_start0
            result.fetch_end0 = right.fetch_end0

            left_only_width = right.fetch_start0 - left.fetch_start0
            right_only_width = right.fetch_end0 - left.fetch_end0
            if left_only_width >= right_only_width:
                result.dict.update(left.dict)

                buffer = list(
                    itertools.takewhile(
                        (lambda x: x[0] not in left.dict),
                        reversed(right.dict.items())
                    )
                )
                result.dict.update(reversed(buffer))
            else:
                result.dict.update(
                    itertools.takewhile(
                        (lambda x: x[0] not in right.dict),
                        left.dict.items()
                    )
                )
                result.dict.update(right.dict)

            return result

    def extend_rightward(self, width):
        for read, uid in itertools.dropwhile(
            (lambda x: x[1] in self.dict),
            readhandler.get_fetch(self.bam, self.chrom, self.fetch_end0, self.fetch_end0 + width, readfilter=self.readfilter, with_uid=True),
        ):
            self.dict[uid] = read

        self.fetch_end0 += width

    def extend_leftward(self, width):
        buffer = list(
            itertools.takewhile(
                (lambda x: x[1] not in self.dict),
                readhandler.get_fetch(self.bam, self.chrom, self.fetch_start0 - width, self.fetch_start0, readfilter=self.readfilter, with_uid=True),
            )
        )
        for read, uid in reversed(buffer):
            self.dict[uid] = read
            self.dict.move_to_end(uid, last=False)

        self.fetch_start0 -= width
            
    # fetch #
    def fetch(self, start0=None, end0=None, with_uid=False):
        if len(self.dict) == 0:
            return iter(())
        else:
            start0, end0 = cluster_fetch_arghandler(start0, end0, self.fetch_start0, self.fetch_end0)
            if start0 >= end0:
                return iter(())
            else:
                yielder = self.get_yielder(with_uid)

                for uid, read in self.dict.items():
                    if read.reference_start >= end0:
                        break
                    if read.reference_end > start0:
                        yield yielder(uid, read)

    def _fetch_forward(self, fetch_start0, fetch_end0, with_uid=False):
        yielder = self.get_yielder(with_uid)
            
        is_area1 = False
        is_area2 = False
        for uid, read in self.dict.items():
            # pre-area1
            if not is_area1:
                # check if within area1
                if read.reference_start > fetch_start0 - MAX_READ_LENGTH:
                    is_area1 = True
                else:
                    # job pre-area1 (null)
                    # continue
                    continue
            # area1
            if not is_area2:
                # check if within area2
                if read.reference_start >= fetch_start0:
                    is_area2 = True
                else:
                    # job area1
                    if read.reference_end > fetch_start0:
                        yield yielder(uid, read)
                    # continue
                    continue
            # check if beyond area2
            if read.reference_start >= fetch_end0:
                break
            # job area2
            yield yielder(uid, read)
                
    def _fetch_reverse(self, fetch_start0, fetch_end0, with_uid=False):
        yielder = self.get_yielder(with_uid)
            
        is_area1 = False
        is_area2 = False
        for uid, read in reversed(self.dict.items()):
        #for uid, read in zip(reversed(self.uids), reversed(self.list)):
            # post-area2
            if not is_area2:
                # check if within area2
                if read.reference_start < fetch_end0:
                    is_area2 = True
                else:
                    # job post-area2 (null)
                    # continue
                    continue
            # area2
            if not is_area1:
                # check if within area1
                if read.reference_start < fetch_start0:
                    is_area1 = True
                else:
                    # job area2
                    yield yielder(uid, read)
                    # continue
                    continue
            # check if within pre-area1
            if read.reference_start <= fetch_start0 - MAX_READ_LENGTH:
                break
            # job area1
            if read.reference_end > fetch_start0:
                yield yielder(uid, read)
                                    
    def fetch_choosing_side(self, start0=None, end0=None, with_uid=False):
        if len(self.dict) == 0:
            return iter(())
        else:
            start0, end0 = cluster_fetch_arghandler(start0, end0, self.fetch_start0, self.fetch_end0)
            if start0 >= end0:
                return iter(())
            else:
                # return self._fetch_forward(start0, end0, with_uid=with_uid)                
                dist_from_start = start0 - self.fetch_start0
                dist_from_end = self.fetch_end0 - end0

                if dist_from_start <= dist_from_end:
                    return self._fetch_forward(
                        start0, end0, with_uid=with_uid
                    )
                else:
                    return reversed(
                        tuple(
                            self._fetch_reverse(
                                start0, end0, with_uid=with_uid
                            )
                        )
                    )


class FetchedReadsSinglechrom(SinglechromBase):
    def __init__(self, chrom, bam, readfilter=None, addition_padding=DEFAULT_FETCHEDREADS_ADDITION_PADDING):
        self.chrom = chrom
        self.bam = bam
        self.addition_padding = addition_padding

        if readfilter is None:
            self.readfilter = readhandler.readfilter_pileup
        else:
            self.readfilter = readfilter

        self.clusters = list()

    def __repr__(self):
        return f'<FetchedReadsSinglechrom object (chrom={self.chrom}, clusters={self.clusters})>'

    @classmethod
    def from_fetch(cls, bam, chrom, start0, end0, readfilter=None, addition_padding=DEFAULT_FETCHEDREADS_ADDITION_PADDING):
        result = cls(chrom, bam, readfilter, addition_padding)
        result.clusters.append(FetchedReadsCluster.from_fetch(result.bam, result.chrom, start0, end0, result.readfilter))
        return result

#        return pr.PyRanges(
#            chromosomes=([self.chrom] * len(self.clusters)),
#            starts=[x.fetch_start0 for x in self.clusters],
#            ends=[x.fetch_end0 for x in self.clusters],
#        )

#    def get_clusters_by_ids(self, ids):
#        result = dict()
#        for x in self.clusters:
#            id = x.get_id()
#            if id in ids:
#                result[id] = x
#        return result

    def spawn_cluster_from_fetch(self, start0, end0):
        return FetchedReadsCluster.from_fetch(self.bam, self.chrom, start0, end0, self.readfilter)

    def get_read(self, read_uid):
        found_read = False
        for frcluster in self.clusters:
            try:
                read = frcluster.dict[read_uid]
            except KeyError:
                continue
            else:
                found_read = True
                break

        if not found_read:
            raise Exception(f'Given ReadUID is not present.')

        return read

    def iter_reads(self):
        return itertools.chain.from_iterable(x.dict.items() for x in self.clusters)

#    def sort_clusters(self):
#        self.clusters.sort(key=(lambda x: x.fetch_start0))

    def add_reads(self, start0, end0):
        self.add(start0, end0)

    def fetch(self, start0, end0, with_uid=False):
        valid_cluster = self.get_valid_cluster(start0, end0)
        return valid_cluster.fetch(start0, end0, with_uid=with_uid)

    def naive_fetch(self, start0, end0, with_uid=False):
        input_gr = pr.PyRanges(chromosomes=[self.chrom], starts=[start0], ends=[end0])
        joined_gr = input_gr.join(self.gr, how=None, apply_strand_suffix=False)

        valid_frcluster_ids = list()
        for idx, row in joined_gr.df.iterrows():
            if row['Start_b'] <= row['Start'] and row['End_b'] >= row['End']:
                valid_frcluster_ids.append(row['Id'])

        if len(valid_frcluster_ids) == 0:
            raise Exception(f'Input fetch range is not fully contained within database.')
        if len(valid_frcluster_ids) > 1:
            raise Exception(f'There are more than one FetchedReadsCluster objects containing input fetch range.')

        frcluster_id = valid_frcluster_ids[0]
        frcluster = self.get_clusters_by_ids([frcluster_id])[frcluster_id]
        return frcluster.fetch(start0, end0, with_uid=with_uid)

    def check_fully_includes(self, start0, end0):
        input_gr = pr.PyRanges(chromosomes=[self.chrom], starts=[start0], ends=[end0])
        intersect_gr = input_gr.intersect(self.gr, how='containment')
        return not intersect_gr.empty

    def choose_cluster(self, start0, end0):
        self_gr = self.gr
        input_gr = pr.PyRanges(chromosomes=[self.chrom], starts=[start0], ends=[end0])
        joined_gr = input_gr.join(self_gr, how=None, apply_strand_suffix=False)
        if joined_gr.df.shape[0] == 0:
            return None
        else:
            overlapping_cluster_ids = {row['Id'] for idx, row in joined_gr.df.iterrows()}
            overlapping_clusters = list()
            for idx, row in self_gr.df.iterrows():
                if row['Id'] in overlapping_cluster_ids:
                    overlapping_clusters.append(self.clusters[idx])

            if len(overlapping_clusters) == 1:
                return overlapping_clusters[0]
            else:
                raise Exception(f'More than one clusters overlap with input range. Overlapping clusters: {overlapping_clusters}')
        

class FetchedReads:
    def __init__(self, bam, readfilter=None, addition_padding=DEFAULT_FETCHEDREADS_ADDITION_PADDING):
        self.bam = bam
        self.addition_padding = addition_padding

        if readfilter is None:
            self.readfilter = readhandler.readfilter_pileup
        else:
            self.readfilter = readfilter

        self.subdata = dict()

    def __repr__(self):
        buffer = list()
        buffer.append(f'<FetchedReads object, which consists of:')
        for chrom, frsinglechrom in self.subdata.items():
            buffer.append(f'\t{frsinglechrom}')
        buffer.append(f'>')
        return '\n'.join(buffer)

    @classmethod
    def from_fetch(cls, bam, chrom, start0, end0, readfilter=None, addition_padding=DEFAULT_FETCHEDREADS_ADDITION_PADDING):
        result = cls(bam, readfilter, addition_padding)
        result.subdata[chrom] = FetchedReadsSinglechrom.from_fetch(result.bam, chrom, start0, end0, result.readfilter, result.addition_padding)

        return result

    def get_read(self, read_uid):
        found_read = False
        for frcluster in itertools.chain.from_iterable(
            (cluster for cluster in frsinglechr.clusters)
            for frsinglechr in self.subdata.values()
        ):
            try:
                read = frcluster.dict[read_uid]
            except KeyError:
                continue
            else:
                found_read = True
                break

        if not found_read:
            raise Exception(f'Given ReadUID is not present.')

        return read

    def iter_reads(self, chrom=None):
        if chrom is None:
            return itertools.chain.from_iterable(x.iter_reads() for x in self.subdata.values())
        else:
            return self.subdata[chrom].iter_reads()

    def add_reads(self, chrom, start0, end0):
        if chrom not in self.subdata:
            self.subdata[chrom] = FetchedReadsSinglechrom(chrom, self.bam, self.readfilter, self.addition_padding)
        self.subdata[chrom].add_reads(start0, end0)

    def fetch(self, chrom, start0, end0, with_uid=False):
        if chrom not in self.subdata:
            self.subdata[chrom] = FetchedReadsSinglechrom(chrom, self.bam, self.readfilter, self.addition_padding)
        return self.subdata[chrom].fetch(start0, end0, with_uid=with_uid)

    def naive_fetch(self, chrom, start0, end0, with_uid=False):
        if chrom not in self.subdata:
            raise Exception(f'Input chrom is not included in the database.')
        return self.subdata[chrom].naive_fetch(start0, end0, with_uid=with_uid)

    def check_fully_includes(self, chrom, start0, end0):
        if chrom not in self.subdata:
            return False
        else:
            return self.subdata[chrom].check_fully_includes(start0, end0)

    def choose_cluster(self, chrom, start0, end0):
        if chrom not in self.subdata:
            return False
        else:
            return self.subdata[chrom].choose_cluster(start0, end0)


####################################


class FetchedRefseqCluster(ClusterBase):
    def __init__(self, fasta, chrom):
        self.fasta = fasta
        self.chrom = chrom
        self.start0 = None
        self.end0 = None
        self.seq = None

    def __repr__(self):
        return f'<FetchedRefseqCluster object (chrom:{self.chrom}, start0:{self.start0}, end0:{self.end0})>'

    @classmethod
    def from_fetch(cls, fasta, chrom, start0, end0):
        result = cls(fasta, chrom)
        result.start0 = start0
        result.end0 = end0
        result.seq = list(fasta.fetch(chrom, start0, end0))
        print('fasta fetch')

        return result

    @property
    def range0(self):
        return range(self.start0, self.end0)

    def get_id(self):
        return (self.chrom, self.start0, self.end0)

    def copy(self):
        result = self.__class__(self.fasta, self.chrom)
        result.start0 = self.start0
        result.end0 = self.end0
        result.seq = self.seq.copy()
        return result

    def fetch(self, start0, end0):
        start0, end0 = cluster_fetch_arghandler(start0, end0, self.start0, self.end0)
        start_idx = start0 - self.start0
        end_idx = end0 - self.start0
        return ''.join(self.seq[start_idx:end_idx])

    def extend_rightward(self, width):
        self.seq.extend(self.fasta.fetch(self.chrom, self.end0, self.end0 + width))
        print('fasta fetch')
        self.end0 += width

    def extend_leftward(self, width):
        for base in reversed(self.fasta.fetch(self.chrom, self.start0 - width, self.start0)):
            self.seq.insert(0, base)
        print('fasta fetch')
        self.start0 -= width

    def merge(self, other):
        if (other.start0 > self.end0) or (self.start0 > other.end0):
            raise Exception(f'Two objects do not overlap.')

        if self.start0 <= other.start0 and self.end0 >= other.end0:
            return self.copy()
        elif other.start0 <= self.start0 and other.end0 >= self.end0:
            return other.copy()
        else:
            # now self.start0 != other.start0
            if self.start0 < other.start0:
                left = self
                right = other
            else:
                left = other
                right = self

            result = self.__class__(self.fasta, self.chrom)
            result.start0 = left.start0
            result.end0 = right.end0
            result.seq = left.seq[:(right.start0 - left.start0)] + right.seq
            return result


class FetchedRefseqSinglechrom(SinglechromBase):
    def __init__(self, chrom, fasta, addition_padding=DEFAULT_FETCHEDREFSEQ_ADDITION_PADDING):
        self.chrom = chrom
        self.fasta = fasta
        self.addition_padding = addition_padding
        self.clusters = list()

    def __repr__(self):
        return f'<FetchedRefseqSinglechrom object (chrom={self.chrom}, clusters={self.clusters})>'

    @classmethod
    def from_fetch(cls, fasta, chrom, start0, end0, addition_padding=DEFAULT_FETCHEDREFSEQ_ADDITION_PADDING):
        result = cls(chrom, fasta, addition_padding)
        result.clusters.append(FetchedRefseqCluster.from_fetch(result.fasta, result.chrom, start0, end0))
        return result

    def spawn_cluster_from_fetch(self, start0, end0):
        return FetchedRefseqCluster.from_fetch(self.fasta, self.chrom, start0, end0)

    def add_seq(self, start0, end0):
        self.add(start0, end0)

    def fetch(self, start0, end0):
        valid_cluster = self.get_valid_cluster(start0, end0)
        return valid_cluster.fetch(start0, end0)


class FetchedRefseq:
    def __init__(self, fasta, addition_padding=DEFAULT_FETCHEDREFSEQ_ADDITION_PADDING):
        self.fasta = fasta
        self.addition_padding = addition_padding
        self.subdata = dict()

    def __repr__(self):
        buffer = list()
        buffer.append(f'<FetchedRefseq object, which consists of:')
        for chrom, frsinglechrom in self.subdata.items():
            buffer.append(f'\t{frsinglechrom}')
        buffer.append(f'>')
        return '\n'.join(buffer)

    def fetch(self, chrom, start0, end0):
        if chrom not in self.subdata:
            self.subdata[chrom] = FetchedRefseqSinglechrom(chrom, self.fasta, self.addition_padding)
        return self.subdata[chrom].fetch(start0, end0)


###########################



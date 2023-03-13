import os
import itertools
import pprint
import collections
import functools

import numpy as np
import pysam

import handygenome.common as common
import handygenome.workflow as workflow
import handygenome.read.readhandler as readhandler
import handygenome.read.alleleinfo as liballeleinfo
import handygenome.read.alleleinfo_sv as liballeleinfo_sv
import handygenome.align.alignhandler as alignhandler
import handygenome.bameditor as bameditor


RPPLIST_MAX_SHOW_LEN = 30
THRESHOLD_TEMPLATE_LENGTH = 1000

FETCH_PADDING_COMMON = 50  # In order to fetch clip-overlapping reads
FETCH_PADDING_SV = 850
FETCH_PADDING_VIEW = 1000
NEW_FETCH_PADDING = 0
LONG_INSERT_THRESHOLD = 1000

LOGGER_RPPLIST = workflow.get_logger(name='ReadPlusPairList', level='warning')

ALLELECLASS_TAG_RPP = 'a2'
ALLELECLASS_TAG_RP = 'a1'


'''
Meaning of range
    - Assumes range.step is either 1 or -1
    - "current position" travels, beginning from "range.start" toward "range.stop".
    - "cursor" position
        - When range.step == 1: immediately left to "current position"
        - When range.step == -1: immediately right to "current position"
    - Area swept by "cursor" is the area represented by range object.
    - When len(range) == 0:
        - When range.step == 1: 0-length interval between (range.start - 1) and (range.start)
        - When range.step == -1: 0-length interval between (range.start) and (range.start + 1)
'''


Clipspec = collections.namedtuple(
    "Clipspec", 
    #("pos0", "is_forward", "seq", "qual", "qname"),
    ("pos0", "is_5prime", "seq", "qual"),
)


class ReadPlus:
    def __init__(self, read, fasta=None, minimal=False, skip_attrs_setting=False, recalc_NMMD=False):
        self.read = read
        # fasta
        if minimal:
            self.fasta = fasta
        else:
            if fasta is None:
                self.fasta = common.DEFAULT_FASTAS[
                    common.infer_refver_bamheader(self.read.header)
                ]
            else:
                self.fasta = fasta
        # pairs_dict
        self.pairs_dict = readhandler.get_pairs_dict(
            self.read, fasta=self.fasta, skip_refseq=minimal, set_cigarop=False,
        )
        # NMMD
        if (
            recalc_NMMD or 
            (
                (not minimal) and 
                (
                    (not self.read.has_tag('NM')) or 
                    (not self.read.has_tag('MD'))
                )
            )
        ):
            self._set_NMMD()
        # others
        self.alleleclass = dict()
        self.monoalt_support = dict()

    def __del__(self):
        del self.read
        del self.pairs_dict
        del self.alleleclass
        del self.monoalt_support

    def __repr__(self):
        qname = self.read.query_name
        chrom = self.read.reference_name
        start1 = self.read.reference_start + 1
        end1 = self.read.reference_end
        region = f'{chrom}:{start1:,}-{end1:,}'
        infostring = f'{qname}; {region}; alleleclass: {self.alleleclass}'
        return (f'<ReadPlus object ({infostring})>')

    def _set_NMMD(self):
        readhandler.set_NMMD(self.read, pairs_dict=self.pairs_dict)

    def _set_pairs_dict(self, skip_refseq=False, set_cigarop=False):
        self.pairs_dict = readhandler.get_pairs_dict(
            self.read, fasta=self.fasta, skip_refseq=skip_refseq, set_cigarop=set_cigarop,
        )
            # keys: 'querypos0', 'refpos0','refseq'

    ##############
    # properties #
    ##############

    @property
    def fiveprime_end(self):
        return readhandler.get_fiveprime_end(self.read)

    @property
    def threeprime_end(self):
        return readhandler.get_threeprime_end(self.read)

    @functools.cached_property
    def cigarstats(self):
        return self.read.get_cigar_stats()[0]

    def get_range0(self):
        return range(self.read.reference_start, self.read.reference_end)

    @property
    def query_name(self):
        return self.read.query_name

    @property
    def range0(self):
        return self.get_range0()

    @property
    def ref_range0(self):  # alias
        return self.range0

    @functools.cache
    def softclip_range0(self):
        return readhandler.get_softclip_ends_range0(self.read)

    #############################################
    # ReadStats non-rppcount attributes-related #
    #############################################

    def get_BQlist(self, vcfspec):
        """For insertion: BQ of the base right before the insertion site
        For deletion: BQ of the base rigth before the deleted seq
        """
        BQlist = list()
        REF_range0 = vcfspec.REF_range0
        if self.check_spans(REF_range0):
            pairs_indexes = self.get_pairs_indexes(
                REF_range0,
                flanking_queryonly_default_mode=False,
                include_leading_queryonly=True, 
                include_trailing_queryonly=True,
            )
            for idx in pairs_indexes:
                querypos0 = self.pairs_dict['querypos0'][idx]
                if querypos0 is not None:
                    BQlist.append(self.read.query_qualities[querypos0])

        return BQlist

    def get_mNM_clipspec_data(self, vcfspec):
        """mNM
            - Softclip is not included
            - Insertion on the right or left border is not included
            - Consecutive mismatches/deletions are split into single-base positions
            - Multi-base insertion is treated as one
        clipspec
            - Leading/trailing insertions are included
            - Those on the right border are assigned to the position on the RIGHT

        Result:
            mNM_data, clipspec_data
            mNM_data:
                list of tuples
                tuple format:
                    - Mismatch: (reference pos0, 0, read base)
                    - Insertion: (reference pos0, 1, inserted bases)
                    - Deletion: (reference pos0, 2)
            clipspec_data:
                list of Clipspec instances
        """
        # set parameters
        mNM_data = list()
        clipspec_data = list()
        vcfspec_ref_range0 = vcfspec.get_range0()
        current_refpos0 = self.read.reference_start

        # main
        for key, subiter in itertools.groupby(
            zip(
                self.pairs_dict['refpos0'],
                self.pairs_dict['querypos0'],
                self.pairs_dict['refseq'],
            ),
            key=(lambda x: ((x[0] is None), (x[1] is None))),
        ):
            if (not key[0]) and key[1]:  # D
                for refpos0, querypos0, refseq in subiter:
                    current_refpos0 += 1
                    if refpos0 in vcfspec_ref_range0:
                        continue

                    mNM_data.append((refpos0, 2))

            elif (not key[0]) and (not key[1]):  # match
                for refpos0, querypos0, refseq in subiter:
                    current_refpos0 += 1
                    if refpos0 in vcfspec_ref_range0:
                        continue

                    if refseq.islower():
                        current_read_base = self.read.query_sequence[querypos0]
                        mNM_data.append((refpos0, 0, current_read_base))

            elif key[0] and (not key[1]):  # I, S
                if current_refpos0 in vcfspec_ref_range0:
                    continue
                else:
                    subiter_list = list(subiter)
                    read_slice = slice(subiter_list[0][1], subiter_list[-1][1] + 1)
                    inserted_seq = self.read.query_sequence[read_slice]

                    if current_refpos0 == self.read.reference_start:
                        clipspec = Clipspec(
                            pos0=current_refpos0, 
                            is_5prime=True,
                            seq=inserted_seq,
                            qual=tuple(self.read.query_qualities)[read_slice],
                        )
                        clipspec_data.append(clipspec)
                    elif current_refpos0 == self.read.reference_end:
                        clipspec = Clipspec(
                            pos0=current_refpos0, 
                            is_5prime=False,
                            seq=inserted_seq,
                            qual=tuple(self.read.query_qualities)[read_slice],
                        )
                        clipspec_data.append(clipspec)
                    else:
                        mNM_data.append((current_refpos0, 1, inserted_seq))

        return mNM_data, clipspec_data


    ######################################################
    # index within pairs (result of get_aligned_pairs()) #
    ######################################################

    def get_pairs_indexes_old(self, range0):
        """Args:
            range0: A directional range
        Returns:
            A directional range composed of pairs_dict indexes.
                (0-length range if the input range0 is 0-length)
        """

        if not self.check_spans(range0):
            raise Exception(f'The read does not span over the input range.')

        start = self.pairs_dict['refpos0'].index(range0.start)
        step = range0.step

        if range0.stop in self.ref_range0:
            stop = self.pairs_dict['refpos0'].index(range0.stop)
        else:
            if range0.step == 1:
                stop = len(self.pairs_dict['refpos0'])
            elif range0.step == -1:
                stop = -1
            
        return range(start, stop, step)

    def get_pairs_indexes(
        self, range0, 
        flanking_queryonly_default_mode=True,
        include_leading_queryonly=False, 
        include_trailing_queryonly=False,
    ):
        """Args:
            range0: A range object. Its step must be +1.
        Returns:
            If read range and target range("range0" argument) overlaps,
                returns a range object representing pairs_dict indexes.
                Partial overlap is okay.
            Otherwise, returns None.
        """
        assert range0.step == 1
        assert len(range0) > 0

        # check if read range contains target range
        if not readhandler.check_overlaps_forward_nonzero(self.range0, range0):
            return None

        # handle all-queryonly cases (added on 230127)
        if self.check_all_queryonly():
            # now overlap is guaranteed
            return range(0, len(self.pairs_dict['refpos0']))

        # handle flanking_queryonly_default_mode
        if flanking_queryonly_default_mode:
            include_leading_queryonly = True
            include_trailing_queryonly = True

        min_ref_pos0 = max(self.read.reference_start, range0.start)
        max_ref_pos0 = min(self.read.reference_end - 1, range0.stop - 1)
        first = self.pairs_dict['refpos0'].index(min_ref_pos0)
        last = self.pairs_dict['refpos0'].index(max_ref_pos0)

        # search for leading query-only cigarops
        if include_leading_queryonly:
            for offset in itertools.count(1):
                pairs_idx = first - offset
                if pairs_idx < 0:
                    break
                if not (
                    (self.pairs_dict['refpos0'][pairs_idx] is None) and
                    (self.pairs_dict['querypos0'][pairs_idx] is not None)
                ):
                    # current position is not query-only
                    break

            first = pairs_idx + 1

        # search for trailing query-only cigarops
        if include_trailing_queryonly:
            for offset in itertools.count(1):
                pairs_idx = last + offset
                if pairs_idx == len(self.pairs_dict['refpos0']):
                    break
                if not (
                    (self.pairs_dict['refpos0'][pairs_idx] is None) and
                    (self.pairs_dict['querypos0'][pairs_idx] is not None)
                ):
                    # current position is not query-only
                    break

            if flanking_queryonly_default_mode:
                if pairs_idx == len(self.pairs_dict['refpos0']):
                    last = pairs_idx - 1
            else:
                last = pairs_idx - 1

        return range(first, last + 1)

    ########################
    # span, match, overlap #
    ########################

    def check_cigarN_includes_range(self, range0):
        return readhandler.check_cigarN_includes_range(
            self.read, range0.start, range0.stop,
        )

    def check_spans(self, range0):
        num_cigarM = self.get_num_cigarM()
        if num_cigarM == 0:
            return False
        else:
            return readhandler.check_spans(self.range0, range0)

    def check_matches(
        self, range0, 
        flanking_queryonly_default_mode=True,
        include_leading_queryonly=False, 
        include_trailing_queryonly=False,
    ):
        pairs_indexes = self.get_pairs_indexes(
            range0, 
            flanking_queryonly_default_mode=flanking_queryonly_default_mode,
            include_leading_queryonly=include_leading_queryonly,
            include_trailing_queryonly=include_trailing_queryonly,
        )
        return self.check_matches_with_pairs_indexes(pairs_indexes)

    def check_matches_with_pairs_indexes(self, pairs_indexes):
        if pairs_indexes is None:
            return False

        for pairs_idx in pairs_indexes:
            querypos0 = self.pairs_dict['querypos0'][pairs_idx]
            refseq = self.pairs_dict['refseq'][pairs_idx]
            matches = ((querypos0 is not None) and 
                       (refseq is not None) and
                       refseq.isupper())
            if not matches:
                return False
        return True

    def check_spans_and_matches(
        self, range0,
        flanking_queryonly_default_mode=True,
        include_leading_queryonly=False, 
        include_trailing_queryonly=False,
    ):
        if self.check_spans(range0):
            pairs_indexes = self.get_pairs_indexes(
                range0,
                flanking_queryonly_default_mode=flanking_queryonly_default_mode,
                include_leading_queryonly=include_leading_queryonly,
                include_trailing_queryonly=include_trailing_queryonly,
            )
            return self.check_matches_with_pairs_indexes(pairs_indexes)
        else:
            return False

    def check_spans_and_matches_vcfspec_flanks(self, vcfspec, flanklen=liballeleinfo.DEFAULT_FLANKLEN):
        preflank_range0, postflank_range0 = vcfspec.get_flank_range0s_equivalents(flanklen=flanklen)
        spans = (self.check_spans(preflank_range0) and self.check_spans(postflank_range0))
        if spans:
            matches = (
                self.check_matches(
                    preflank_range0,
                    flanking_queryonly_default_mode=False,
                    include_leading_queryonly=True, 
                    include_trailing_queryonly=False,
                ) and 
                self.check_matches(
                    postflank_range0,
                    flanking_queryonly_default_mode=False,
                    include_leading_queryonly=False, 
                    include_trailing_queryonly=True,
                )
            )
        else:
            matches = False

        return spans, matches

    def check_overlaps(self, range0):
        return readhandler.check_overlaps(self.range0, range0)

    def check_overlaps_vcfspec(self, vcfspec):
        return (
            self.read.reference_name == vcfspec.chrom and
            self.check_overlaps(vcfspec.REF_range0)
        )

    @functools.cache
    def get_softclip_range0(self):
        return readhandler.get_softclip_ends_range0(self.read)
        
    def check_softclip_overlaps(self, range0):
        softclip_range0 = self.get_softclip_range0()
        return readhandler.check_overlaps(softclip_range0, range0)

    def check_softclip_overlaps_vcfspec(self, vcfspec):
        return (
            self.read.reference_name == vcfspec.chrom and
            self.check_softclip_overlaps(vcfspec.REF_range0)
        )

    def check_softclip_spans_vcfspec_flanks(self, vcfspec, flanklen=liballeleinfo.DEFAULT_FLANKLEN):
        preflank_range0, postflank_range0 = vcfspec.get_flank_range0s_equivalents(flanklen=flanklen)
        softclip_range0 = self.get_softclip_range0()
        return (
            readhandler.check_spans(softclip_range0, preflank_range0)
            and readhandler.check_spans(softclip_range0, postflank_range0)
        )

    #################################
    # get sequence from coordinates #
    #################################

    def get_seq_from_coord(self, start0, length, forward):
        """If there is not enough query bases compared to "length" argument,
            return string length may be shorter than "length" argument.
        """

        assert start0 in self.ref_range0, f'The read does not span "start0".'

        valid_idx_range = range(0, len(self.pairs_dict['refpos0']))
        pairs_idx = self.pairs_dict['refpos0'].index(start0)
        step = 1 if forward else -1
        querypos0_list = list()
        while True:
            current = self.pairs_dict['querypos0'][pairs_idx]
            if current is not None:
                querypos0_list.append(current)
            if len(querypos0_list) >= length:
                break

            pairs_idx += step
            if pairs_idx not in valid_idx_range:
                break

        return self.get_seq_from_querypos0_list(querypos0_list)

    def get_seq_from_range0(
        self, range0,
        flanking_queryonly_default_mode=True,
        include_leading_queryonly=False, 
        include_trailing_queryonly=False,
    ):
        pairs_indexes = self.get_pairs_indexes(
            range0,
            flanking_queryonly_default_mode=flanking_queryonly_default_mode,
            include_leading_queryonly=include_leading_queryonly,
            include_trailing_queryonly=include_trailing_queryonly,
        )
        return self.get_seq_from_pairs_indexes(pairs_indexes)

    def get_seq_from_pairs_indexes(self, pairs_indexes):
        if pairs_indexes is None:
            return None

        query_idx_list = list()
        for pairs_idx in pairs_indexes:
            query_idx = self.pairs_dict['querypos0'][pairs_idx]
            if query_idx is not None:
                query_idx_list.append(query_idx)

        if len(query_idx_list) == 0:
            return ''
        else:
            return self.read.query_sequence[query_idx_list[0]:(query_idx_list[-1] + 1)]
        
    def get_seq_from_querypos0_list(self, querypos0_list):
        if len(querypos0_list) == 0:
            return ''
        else:
            start = min(querypos0_list)
            stop = max(querypos0_list) + 1
            return self.read.query_sequence[slice(start, stop, 1)]

    #########################################################################
    # get read-oriented coordinate(querypos0) of a certain genomic location #
    #########################################################################

    def get_querypos0_of_range0(self, range0, mode='left', fraction=False):
        querypos0_fromleft = self._get_querypos0_of_range0_fromleft(range0)
        if querypos0_fromleft is None:
            return None
        else:
            querypos0 = self._convert_querypos0_of_range0(
                querypos0_fromleft, mode=mode, fraction=fraction)
            return querypos0

    def get_querypos0_of_range0_allmodes(self, range0):
        querypos0_fromleft = self._get_querypos0_of_range0_fromleft(range0)
        if querypos0_fromleft is None:
            result = {'left': None, 'right': None, 
                      '5prime': None, '3prime': None,
                      'left_fraction': None, 'right_fraction': None,
                      '5prime_fraction': None, '3prime_fraction': None}
        else:
            querypos0_fromright = self._get_querypos0_of_range0_fromright(
                querypos0_fromleft)
            if self.read.is_forward:
                querypos0_from5 = querypos0_fromleft
                querypos0_from3 = querypos0_fromright
            else:
                querypos0_from5 = querypos0_fromright
                querypos0_from3 = querypos0_fromleft

            readlen = len(self.read.query_sequence)
            result = {
                'left': querypos0_fromleft, 
                'right': querypos0_fromright, 
                '5prime': querypos0_from5,
                '3prime': querypos0_from3,
                'left_fraction': querypos0_fromleft / readlen,
                'right_fraction': querypos0_fromright / readlen,
                '5prime_fraction': querypos0_from5 / readlen,
                '3prime_fraction': querypos0_from3 / readlen,
            }

        return result

    def _get_querypos0_of_range0_fromleft(self, range0):
        if not self.check_spans(range0):
            return None
        else:
            start_idx = self.pairs_dict['refpos0'].index(range0.start)
            querypos0_fromleft = self.pairs_dict['querypos0'][start_idx]
            return querypos0_fromleft

    def _get_querypos0_of_range0_fromright(self, querypos0_fromleft):
        querypos0_fromright = (len(self.read.query_sequence)
                                     - querypos0_fromleft - 1)
        return querypos0_fromright

    def _convert_querypos0_of_range0(self, querypos0_fromleft,
                                     mode='left', fraction=False):
        assert mode in ('left', 'right', '5prime', '3prime'), (
            f'"mode" argument must be one of '
            f'"left", "right", "5prime", "3prime".')

        if mode == 'left':
            result = querypos0_fromleft
        else:
            querypos0_fromright = self._get_querypos0_of_range0_fromright(
                querypos0_fromleft)
            if mode == 'right':
                result = querypos0_fromright
            elif mode == '5prime':
                if self.read.is_forward:
                    result = querypos0_fromleft
                else:
                    result = querypos0_fromright
            elif mode == '3prime':
                if self.read.is_reverse:
                    result = querypos0_fromleft
                else:
                    result = querypos0_fromright

        if fraction:
            result = result / len(self.read.query_sequence)

        return result

    #########
    # cigar #
    #########

    def walk_cigar(self):
        return alignhandler.walk_cigar(self.read.cigartuples, self.read.reference_start)

    def split_cigar(self, split_range0):
        return readhandler.split_cigar(self.read.cigartuples, self.read.reference_start, split_range0)

    def get_num_cigarM(self):
        return self.read.get_cigar_stats()[0][0]

    def check_all_queryonly(self):
        return readhandler.check_all_queryonly(self.read)

    ###############
    # alleleclass #
    ###############

    def update_alleleclass(self, vcfspec, **kwargs):
        #self.alleleclass[vcfspec] = liballeleinfo.get_alleleclass_asis_readplus(
        #    vcfspec=vcfspec, rp=self, flanklen=flanklen,
        #)
        self.alleleclass[vcfspec] = liballeleinfo.get_alleleclass_asis_readplus_new(
            vcfspec=vcfspec, rp=self, **kwargs, 
        )

    def set_alleleclass_tag(self, vcfspec):
        value = str(self.alleleclass[vcfspec])
        self.read.set_tag(ALLELECLASS_TAG_RP, value, 'Z', replace=True)

    def update_monoalt_support(self, vcfspec, flanklen=liballeleinfo.DEFAULT_FLANKLEN):
        self.monoalt_support[vcfspec] = liballeleinfo.get_normalized_monoalt_supports_readplus(
            vcfspec, self, flanklen=flanklen,
        )

    def update_alleleinfo_sv(
            self, bnds, 
            flanklen_parside=liballeleinfo_sv.DEFAULT_FLANKLEN_PARSIDE,
            flanklen_bndside=liballeleinfo_sv.DEFAULT_FLANKLEN_BNDSIDE):
        aiitem_bnd1 = liballeleinfo_sv.make_alleleinfoitem_readplus(
            bnds=bnds, is_bnd1=True, rp=self, 
            flanklen_parside=flanklen_parside,
            flanklen_bndside=flanklen_bndside)
        aiitem_bnd2 = liballeleinfo_sv.make_alleleinfoitem_readplus(
            bnds=bnds, is_bnd1=False, rp=self, 
            flanklen_parside=flanklen_parside,
            flanklen_bndside=flanklen_bndside)

        # sanity checks
        if aiitem_bnd1['alt_support'] and aiitem_bnd2['alt_support']:
            raise Exception(
                f'UNEXPECTED SV SUPPORT PATTERN: '
                f'This read supports both bnd1 and bnd2!:\n'
                f'read: {self.read.to_string()}\n'
                f'Breakends: {bnds}')

        # update alleleinfo dict
        self.alleleinfo[bnds] = {'bnd1': aiitem_bnd1, 'bnd2': aiitem_bnd2}

    def set_alleleinfo_tag_sv(self, bnds):
        alleleclass_list_bnd1 = [
            k 
            for (k, v) in self.alleleinfo[bnds]['bnd1'].items()
            if v
        ]
        alleleclass_list_bnd2 = [
            k 
            for (k, v) in self.alleleinfo[bnds]['bnd2'].items()
            if v
        ]

        if len(alleleclass_list_bnd1) == 0 or len(alleleclass_list_bnd2) == 0:
            raise Exception(
                f'Invalid alleleinfoitem.\n'
                f'ReadPlusPair: {self}'
            )
        else:
            alleleclass_bnd1 = '&'.join(alleleclass_list_bnd1)
            alleleclass_bnd2 = '&'.join(alleleclass_list_bnd2)

        value = '_'.join([f'bnd1={alleleclass_bnd1}', f'bnd2={alleleclass_bnd2}'])

        self.read.set_tag(ALLELECLASS_TAG_RP, value, 'Z', replace=True)

    #################
    # miscellaneous #
    #################

    def get_distance(self, range0):
        if self.check_overlaps(range0):
            distance = 0
        else:
            if self.read.reference_end <= min(range0):
                distance = min(range0) - self.read.reference_end + 1
            elif self.read.reference_start > max(range0):
                distance = self.read.reference_start - max(range0)

        return distance

    def _set_SAlist(self):
        """cigartuples pattern check is done:
            asserts that one end is softclip and the other end is match.
        """

        def SA_cigartuples_sanitycheck(cigartuples, read):
            if not (
                    (cigartuples[0][0] == 4 and cigartuples[-1][0] == 0) or
                    (cigartuples[-1][0] == 4 and cigartuples[0][0] == 0)):
                raise Exception(
                    f'Unexpected SA cigarstring pattern:\n{read.to_string()}')

        if self.read.has_tag('SA'):
            self.SAlist = list()
            SAtag_split = self.read.get_tag('SA').strip(';').split(';')
            for field in SAtag_split:
                field_sp = field.split(',')

                # 
                MQ = int(field_sp[4])
                #if MQ == 0:
                #    continue

                chrom = field_sp[0]
                pos = int(field_sp[1])

                if field_sp[2] == '+':
                    is_forward = True
                elif field_sp[2] == '-':
                    is_forward = False
                else:
                    raise Exception(
                        f'Unexpected strand string from SA tag. read:\n'
                        f'{self.read.to_string()}')

                cigarstring = field_sp[3]
                cigartuples = readhandler.get_cigartuples(cigarstring)
                #SA_cigartuples_sanitycheck(cigartuples, self.read)

                NM = int(field_sp[5])

                SAitem = {'chrom': chrom, 'pos': pos, 
                          'is_forward': is_forward,
                          'MQ': MQ, 'NM': NM,
                          'cigarstring': cigarstring, 
                          'cigartuples': cigartuples}
                self.SAlist.append(SAitem)

            if len(self.SAlist) == 0:
                self.SAlist = None
        else:
            self.SAlist = None


class ReadPlusPair:
    """
    Requires two readplus objects which are mates of each other and all of 
        primary alignment.
    Non-primary reads with the same queryname is stored in 
        'rplist_nonprimary' attribute.
    """

    def __init__(
        self, 
        rplist_primary, 
        rplist_nonprimary, 
        chromdict,
        threshold_tlen=THRESHOLD_TEMPLATE_LENGTH,
    ):
        self._set_rp1_rp2(rplist_primary, chromdict)
        self.rplist_nonprimary = rplist_nonprimary
        self.alleleclass = dict()
        self.monoalt_support = dict()

        #self._set_is_proper_pair()
        #self._set_sv_supporting()
        #self.irrelevant = (self.rp1.irrelevant or self.rp2.irrelevant)

    def __del__(self):
        del self.rp1
        del self.rp2
        for read in self.rplist_nonprimary:
            del read
        del self.rplist_nonprimary
        del self.alleleclass
        del self.monoalt_support

    def __repr__(self):
        qname = self.query_name

        rp1_chrom = self.rp1.read.reference_name
        rp1_start1 = self.rp1.read.reference_start + 1
        rp1_end1 = self.rp1.read.reference_end
        rp1_region = f'{rp1_chrom}:{rp1_start1:,}-{rp1_end1:,}'

        if self.rp2 is None:
            rp2_region = 'None'
        else:
            rp2_chrom = self.rp2.read.reference_name
            rp2_start1 = self.rp2.read.reference_start + 1
            rp2_end1 = self.rp2.read.reference_end
            rp2_region = f'{rp2_chrom}:{rp2_start1:,}-{rp2_end1:,}'

        infostring = (f'{qname}; rp1: {rp1_region}; rp2: {rp2_region}; '
                      f'alleleclass: {self.alleleclass}')

        return (f'<ReadPlusPair object ({infostring})>')

    def _set_rp1_rp2(self, rplist_primary, chromdict):
        if len(rplist_primary) == 1:
            self.rp1 = rplist_primary[0]
            self.rp2 = None
        elif len(rplist_primary) == 2:
            order = common.compare_coords(
                rplist_primary[0].read.reference_name, 
                rplist_primary[0].fiveprime_end,
                rplist_primary[1].read.reference_name, 
                rplist_primary[1].fiveprime_end,
                chromdict)
            if order <= 0:
                self.rp1 = rplist_primary[0]
                self.rp2 = rplist_primary[1]
            else:
                self.rp1 = rplist_primary[1]
                self.rp2 = rplist_primary[0]

    ##############
    # properties #
    ##############

    @property
    def query_name(self):
        return self.rp1.read.query_name

    @property
    def mate_chroms_differ(self):
        if self.rp2 is None:
            return None
        else:
            return self.rp1.read.reference_name != self.rp2.read.reference_name

    is_TRA = mate_chroms_differ  # alias

    @property
    def tlen(self):
        if self.rp2 is None:
            return None
        else:
            if self.mate_chroms_differ:
                return None
            else:
                return abs(self.rp1.read.template_length)

    template_length = tlen  # alias

    @property
    def pairorient(self):
        return readhandler.get_pairorient(self.rp1.read)
            # may be None when: mate unmapped, TLEN == 0

    ##########################
    # non-alleleinfo methods #
    ##########################

    def get_range0(self, chrom):
        relevant_rps = list()
        if self.rp1.read.reference_name == chrom:
            relevant_rps.append(self.rp1)
        if self.rp2 is not None:
            if self.rp2.read.reference_name == chrom:
                relevant_rps.append(self.rp2)

        if len(relevant_rps) == 0:
            return None
        else:
            return range(
                min(x.read.reference_start for x in relevant_rps),
                max(x.read.reference_end for x in relevant_rps),
            )

    def check_softclip_overlaps_vcfspec(self, vcfspec):
        if self.rp2 is None:
            return self.rp1.check_softclip_overlaps_vcfspec(vcfspec)
        else:
            if self.rp1.check_overlaps_vcfspec(vcfspec):
                return False
            elif self.rp2.check_overlaps_vcfspec(vcfspec):
                return False
            else:
                return (self.rp1.check_softclip_overlaps_vcfspec(vcfspec) or
                        self.rp2.check_softclip_overlaps_vcfspec(vcfspec))

    ##############
    # alleleinfo #
    ##############

    def update_alleleclass(self, vcfspec, **kwargs):
        # rp1
        self.rp1.update_alleleclass(vcfspec, **kwargs)
        # rp2
        if self.rp2 is None:
            self.alleleclass[vcfspec] = self.rp1.alleleclass[vcfspec]
        else:
            self.rp2.update_alleleclass(vcfspec, **kwargs)
            self.alleleclass[vcfspec] = liballeleinfo.merge_alleleclasses(
                self.rp1.alleleclass[vcfspec],
                self.rp2.alleleclass[vcfspec],
            )

    def set_alleleclass_tag(self, vcfspec):
        value = str(self.alleleclass[vcfspec])
        self.rp1.read.set_tag(ALLELECLASS_TAG_RPP, value, 'Z', replace=True)
        if self.rp2 is not None:
            self.rp2.read.set_tag(ALLELECLASS_TAG_RPP, value, 'Z', replace=True)

    def update_monoalt_support(self, vcfspec, flanklen=liballeleinfo.DEFAULT_FLANKLEN):
        self.rp1.update_monoalt_support(vcfspec, flanklen=flanklen)
        if self.rp2 is None:
            self.monoalt_support[vcfspec] = self.rp1.monoalt_support[vcfspec]
        else:
            self.rp2.update_monoalt_support(vcfspec, flanklen=flanklen)
            self.monoalt_support[vcfspec] = liballeleinfo.merge_monoalt_supports(
                self.rp1.monoalt_support[vcfspec], 
                self.rp2.monoalt_support[vcfspec],
            )

    def update_alleleinfo_sv(
            self, bnds, 
            flanklen_parside=liballeleinfo_sv.DEFAULT_FLANKLEN_PARSIDE,
            flanklen_bndside=liballeleinfo_sv.DEFAULT_FLANKLEN_BNDSIDE):
        self.rp1.update_alleleinfo_sv(bnds,
                                      flanklen_parside=flanklen_parside,
                                      flanklen_bndside=flanklen_bndside)
        self.rp2.update_alleleinfo_sv(bnds,
                                      flanklen_parside=flanklen_parside,
                                      flanklen_bndside=flanklen_bndside)
        aiitem = liballeleinfo_sv.make_alleleinfoitem_readpluspair(
            aiitem_bnd1_rp1=self.rp1.alleleinfo[bnds]['bnd1'], 
            aiitem_bnd1_rp2=self.rp2.alleleinfo[bnds]['bnd1'],
            aiitem_bnd2_rp1=self.rp1.alleleinfo[bnds]['bnd2'], 
            aiitem_bnd2_rp2=self.rp2.alleleinfo[bnds]['bnd2'],
            bnds=bnds)

        self.alleleinfo[bnds] = aiitem

    def set_alleleinfo_tag_sv(self, bnds):
        alleleclass_list = [k for (k, v) 
                            in self.alleleinfo[bnds].items()
                            if v]
        if len(alleleclass_list) == 0:
            raise Exception(
                f'Invalid alleleinfoitem: {self.alleleinfo[bnds]}\n'
                f'ReadPlusPair: {self}')
        else:
            alleleclass = '&'.join(alleleclass_list)

        #value = f'{bnds.get_id()}_{alleleclass}'
        value = str(alleleclass)

        self.rp1.read.set_tag(ALLELECLASS_TAG_RPP, value, 'Z', replace=True)
        self.rp2.read.set_tag(ALLELECLASS_TAG_RPP, value, 'Z', replace=True)


    #############################################
    # ReadStats non-rppcount attributes-related #
    #############################################

    def get_MQ(self):
        if self.rp2 is None:
            return self.rp1.read.mapping_quality
        else:
            return (self.rp1.read.mapping_quality + self.rp2.read.mapping_quality) / 2

    def get_cliplen(self):
        if self.rp2 is None:
            return self.rp1.cigarstats[4]
        else:
            return self.rp1.cigarstats[4] + self.rp2.cigarstats[4]

    def get_mNM_clipspec_data(self, vcfspec):
        if self.rp2 is None:
            return self.rp1.get_mNM_clipspec_data(vcfspec)
        else:
            mNM_data_rp1, clipspec_data_rp1 = self.rp1.get_mNM_clipspec_data(vcfspec)
            mNM_data_rp2, clipspec_data_rp2 = self.rp2.get_mNM_clipspec_data(vcfspec)

            mNM_data = sorted(
                set(itertools.chain(mNM_data_rp1, mNM_data_rp2)),
                key=(lambda x: x[0])
            )
            clipspec_data = sorted(
                set(itertools.chain(clipspec_data_rp1, clipspec_data_rp2)),
                key=(lambda x: x[0])
            )
            return mNM_data, clipspec_data



#    def _set_is_proper_pair(self):
#        if self.pairorient is None:
#            self.is_proper_pair = False
#        else:
#            self.is_proper_pair = (self.pairorient[0] == 'F' and 
#                                   self.pairorient[2] == 'R')

#    def _set_sv_supporting(self):
#        if self.mate_chroms_differ:
#            self.SV_supporting = True
#        else:
#            if self.tlen == 0:
#                self.SV_supporting = False
#            else:
#                if not self.is_proper_pair:
#                    self.SV_supporting = True
#                else:
#                    self.SV_supporting = (
#                        self.tlen > THRESHOLD_TEMPLATE_LENGTH)



##################################################
# ReadPlusPairList-related classes and functions #
##################################################

class ReadPlusPairList(list):

    def __init__(self, chromdict):
        self.chromdict = chromdict

    def __del__(self):
        for rpp in self:
            del rpp

    def select_by_qname(self, qname):
        for rpp in self:
            if rpp.query_name == qname:
                return rpp
        raise Exception(f'QNAME {qname} is not present')

    def sortby_rp1(self):
        read_sortkey = common.get_read_sortkey(self.chromdict)
        sortkey = lambda rpp: read_sortkey(rpp.rp1.read)
        self.sort(key=sortkey)

    def get_ref_range0(self, chrom):
        start_list = list()
        end_list = list()
        for rpp in self:
            if rpp.rp1.read.reference_name == chrom:
                start_list.append(rpp.rp1.read.reference_start)
                end_list.append(rpp.rp1.read.reference_end)
            if rpp.rp2 is not None:
                if rpp.rp2.read.reference_name == chrom:
                    start_list.append(rpp.rp2.read.reference_start)
                    end_list.append(rpp.rp2.read.reference_end)

        ref_range0 = range(min(start_list), max(end_list))
        
        return ref_range0

    get_range0 = get_ref_range0

    def get_ref_range0_sv(self, bnds):
        start_list_bnd1 = list()
        end_list_bnd1 = list()
        start_list_bnd2 = list()
        end_list_bnd2 = list()

        for rpp in self:
            for rp in (rpp.rp1, rpp.rp2):
                bnd1_relevant = (
                    rp.read.reference_name == bnds.chrom_bnd1 and
                    (rp.get_distance(bnds.get_pos_range0_bnd1())
                     < liballeleinfo_sv.BND_DISTANCE_THRESHOLD))
                bnd2_relevant = (
                    rp.read.reference_name == bnds.chrom_bnd2 and
                    (rp.get_distance(bnds.get_pos_range0_bnd2())
                     < liballeleinfo_sv.BND_DISTANCE_THRESHOLD))

                if bnd1_relevant:
                    start_list_bnd1.append(rp.read.reference_start)
                    end_list_bnd1.append(rp.read.reference_end)
                if bnd2_relevant:
                    start_list_bnd2.append(rp.read.reference_start)
                    end_list_bnd2.append(rp.read.reference_end)

        ref_range0_bnd1 = range(min(start_list_bnd1), max(end_list_bnd1))
        ref_range0_bnd2 = range(min(start_list_bnd2), max(end_list_bnd2))

        return ref_range0_bnd1, ref_range0_bnd2

    def update_alleleclass(self, vcfspec, **kwargs):
        for rpp in self:
            rpp.update_alleleclass(vcfspec, **kwargs)

    def update_monoalt_support(self, vcfspec, flanklen=liballeleinfo.DEFAULT_FLANKLEN):
        for rpp in self:
            rpp.update_monoalt_support(vcfspec, flanklen=flanklen)

    def summarize_monoalt_support(self):
        result = dict()
        for rpp in self:
            for vcfspec, monoalt_support in rpp.monoalt_support.items():
                result.setdefault(vcfspec, dict())
                for alt_index, support_val in monoalt_support.items():
                    result[vcfspec].setdefault(alt_index, 0)
                    if support_val:
                        result[vcfspec][alt_index] += 1

        return result

    def update_alleleinfo_sv(
            self, bnds, 
            flanklen_parside=liballeleinfo_sv.DEFAULT_FLANKLEN_PARSIDE,
            flanklen_bndside=liballeleinfo_sv.DEFAULT_FLANKLEN_BNDSIDE):
        for rpp in self:
            rpp.update_alleleinfo_sv(bnds, flanklen_parside=flanklen_parside,
                                     flanklen_bndside=flanklen_bndside)

    def set_alleleclass_tag_rpp(self, vcfspec):
        for rpp in self:
            rpp.set_alleleclass_tag(vcfspec)

    def set_alleleinfo_tag_rpp_sv(self, bnds):
        for rpp in self:
            rpp.set_alleleinfo_tag_sv(bnds)

    def set_alleleclass_tag_rp(self, vcfspec):
        for rpp in self:
            rpp.rp1.set_alleleclass_tag(vcfspec)
            if rpp.rp2 is not None:
                rpp.rp2.set_alleleclass_tag(vcfspec)

    def set_alleleinfo_tag_rp_sv(self, bnds):
        for rpp in self:
            rpp.rp1.set_alleleinfo_tag_sv(bnds)
            rpp.rp2.set_alleleinfo_tag_sv(bnds)

    def write_bam(self, outfile_path=None, outfile_dir=None):
        # set outfile_path 
        if outfile_path is None:
            if outfile_dir is None:
                outfile_path=workflow.get_tmpfile_path(suffix='.bam')
            else:
                outfile_path=workflow.get_tmpfile_path(suffix='.bam', where=outfile_dir)

        # get a list of reads
        readlist = list()
        for rpp in self:
            readlist.append(rpp.rp1.read)
            if rpp.rp2 is not None:
                readlist.append(rpp.rp2.read)

        # sort & write
        if not bameditor.check_header_compatibility(
            [rpp.rp1.read.header for rpp in self]
        ):
            raise Exception(f'Contig specs are different between read headers.')

        bamheader = self[0].rp1.read.header.copy()
        sortkey = common.get_read_sortkey(common.ChromDict(bamheader=bamheader))
        readlist.sort(key=sortkey)

        with pysam.AlignmentFile(outfile_path, mode='wb', header=bamheader) as out_bam:
            for read in readlist:
                out_bam.write(read)
        # index
        pysam.index(outfile_path)

        #self.bam_path = outfile_path


def get_rpplist_nonsv(
    bam, chrom, start0, end0, 
    fasta=None, 
    chromdict=None, 
    view=False, 
    no_matesearch=True,
    fetch_padding_common=FETCH_PADDING_COMMON,
    fetch_padding_view=FETCH_PADDING_VIEW,
    new_fetch_padding=NEW_FETCH_PADDING,
    long_insert_threshold=LONG_INSERT_THRESHOLD,
    recalc_NMMD=False,
    include_irrelevant_reads=False,
):
    if fasta is None or chromdict is None:
        refver = common.infer_refver_bamheader(bam.header)
    if fasta is None:
        fasta = common.DEFAULT_FASTAS[refver]
    if chromdict is None:
        chromdict = common.ChromDict(refver=refver)

    LOGGER_RPPLIST.info('Beginning initial fetch')
    chromlen = chromdict[chrom]
    (relevant_qname_set, new_fetch_range) = initial_fetch_nonsv(
        bam=bam, 
        chrom=chrom, 
        start0=start0, 
        end0=end0, 
        chromlen=chromlen,
        view=view,
        fetch_padding_common=fetch_padding_common, 
        fetch_padding_view=fetch_padding_view,
        new_fetch_padding=new_fetch_padding, 
        long_insert_threshold=long_insert_threshold,
    )

    rpplist = ReadPlusPairList(chromdict=chromdict)
    if (
        include_irrelevant_reads
        or ((not include_irrelevant_reads) and (len(relevant_qname_set) > 0))
    ):
    #if len(relevant_qname_set) > 0:
        LOGGER_RPPLIST.info('Beginning refined fetch')
        fetchresult_dict = refined_fetch_nonsv(
            bam, chrom, new_fetch_range, relevant_qname_set, include_irrelevant_reads,
        )

        LOGGER_RPPLIST.info('Beginning assembly into readpluspair')
        for readlist in fetchresult_dict.values():
            rpp = get_rpp_from_refinedfetch(
                readlist, bam, fasta, chromdict, no_matesearch, recalc_NMMD=recalc_NMMD,
            )
            del readlist
            if rpp is not None:
                rpplist.append(rpp)

        # sort
        rpplist.sortby_rp1()

    return rpplist


def get_rpplist_sv(bam, fasta, chromdict, bnds, view=False,
                   no_matesearch=False,
                   fetch_padding_common=FETCH_PADDING_COMMON,
                   fetch_padding_sv=FETCH_PADDING_SV,
                   fetch_padding_view=FETCH_PADDING_VIEW,
                   new_fetch_padding=NEW_FETCH_PADDING,
                   long_insert_threshold=LONG_INSERT_THRESHOLD):
    """Args:
        no_matesearch: Only for compatibility with VariantPlus.make_rpplist 
            method. Its input value is not effective.
    """
    LOGGER_RPPLIST.info('Beginning initial fetch')
    (relevant_qname_set_bnd1, 
     new_fetch_range_bnd1,
     relevant_qname_set_bnd2, 
     new_fetch_range_bnd2,
     relevant_qname_set_union) = initial_fetch_sv(
        bam, bnds, view, 
        fetch_padding_common, fetch_padding_sv, fetch_padding_view,
        new_fetch_padding, long_insert_threshold)

    rpplist = ReadPlusPairList(chromdict=chromdict)
    if len(relevant_qname_set_union) > 0:
        LOGGER_RPPLIST.info('Beginning refined fetch')
        fetchresult_dict = refined_fetch_sv(bam, bnds, new_fetch_range_bnd1,
                                            new_fetch_range_bnd2,
                                            relevant_qname_set_union)

        LOGGER_RPPLIST.info('Beginning assembly into readpluspair')
        for readlist in fetchresult_dict.values():
            rpp = get_rpp_from_refinedfetch(readlist, bam, fasta, chromdict,
                                            no_matesearch=False)
            if rpp is not None:
                rpplist.append(rpp)

        # sort
        rpplist.sortby_rp1()

    return rpplist


# deprecated
def get_rpplist_sv_separately(bam, fasta, chromdict, bnds, view=False,
                   fetch_padding_common=FETCH_PADDING_COMMON,
                   fetch_padding_sv=FETCH_PADDING_SV,
                   fetch_padding_view=FETCH_PADDING_VIEW,
                   new_fetch_padding=NEW_FETCH_PADDING,
                   long_insert_threshold=LONG_INSERT_THRESHOLD):
    LOGGER_RPPLIST.info('Beginning initial fetch')
    (relevant_qname_set_bnd1, 
     new_fetch_range_bnd1,
     relevant_qname_set_bnd2, 
     new_fetch_range_bnd2,
     relevant_qname_set_union) = initial_fetch_sv(bam, bnds, view, 
                                                  fetch_padding_common,
                                                  fetch_padding_sv, 
                                                  fetch_padding_view,
                                                  new_fetch_padding, 
                                                  long_insert_threshold)
    if len(relevant_qname_set_union) == 0:
        rpplist_bnd1 = ReadPlusPairList(chromdict=chromdict)
        rpplist_bnd2 = ReadPlusPairList(chromdict=chromdict)
    else:
        LOGGER_RPPLIST.info('Beginning refined fetch')
        fetchresult_dict = refined_fetch_sv(bam, bnds, new_fetch_range_bnd1, 
                                            new_fetch_range_bnd2,
                                            relevant_qname_set_union)

        LOGGER_RPPLIST.info('Beginning assembly into readpluspair')
        rpplist_bnd1 = ReadPlusPairList(chromdict=chromdict)
        rpplist_bnd2 = ReadPlusPairList(chromdict=chromdict)
        for readlist in fetchresult_dict.values():
            rpp = get_rpp_from_refinedfetch(readlist, bam, fasta, chromdict)
            if rpp is not None:
                if relevant_qname_list_bnd1 is not None:
                    if rpp.query_name in relevant_qname_list_bnd1:
                        rpplist_bnd1.append(rpp)
                if relevant_qname_list_bnd2 is not None:
                    if rpp.query_name in relevant_qname_list_bnd2:
                        rpplist_bnd2.append(rpp)

        # sort
        rpplist_sortkey = lambda rpp: min(rpp.rp1.read.reference_start,
                                          rpp.rp2.read.reference_start)
        rpplist_bnd1.sort(key=rpplist_sortkey)
        rpplist_bnd2.sort(key=rpplist_sortkey)

    return rpplist_bnd1, rpplist_bnd2


#### initial fetch ####

def readfilter_matesearch(read, long_insert_threshold):
    """For use in "initial_fetch" function.
    The read is selected if True.
    For decision of whether mate start position will be included in the 
        "new_fetch_range".
    """

    return (read.reference_name == read.next_reference_name and
            abs(read.template_length) < long_insert_threshold)


def initial_fetch_nonsv(
    bam, chrom, start0, end0, chromlen, view,
    fetch_padding_common, fetch_padding_view,
    new_fetch_padding, long_insert_threshold,
):
    """Read filtering done by readhandler.readfilter_bad_read"""

    def readfilter_qname(read, start0, end0):
        """The read is selected if True.
        For decision of whether the read query_name will be included in the
            "relevant_qname_list", for non-SV.
        """

        softclip_range0 = readhandler.get_softclip_ends_range0(read)
        return (
            softclip_range0.stop > start0 
            and softclip_range0.start < end0
        )
        #return True

    relevant_qname_set = set()
    start0_set = set()
    end0_set = set()

    initial_fetch_padding = (fetch_padding_view if view else fetch_padding_common)
    fetcher = readhandler.get_fetch(
        bam, chrom, 
        start=(start0 - initial_fetch_padding),
        end=(end0 + initial_fetch_padding),
        readfilter=readhandler.readfilter_bad_read,
    )

    for read in fetcher:
        read_is_relevant = (True if view else readfilter_qname(read, start0, end0))
        if read_is_relevant:
            relevant_qname_set.add(read.query_name)
            start0_set.add(read.reference_start)
            end0_set.add(read.reference_end)
            if readfilter_matesearch(read, long_insert_threshold):
                start0_set.add(read.next_reference_start)
                if read.template_length > 0:
                    end0_set.add(read.reference_start + read.template_length)

        del read

    if len(start0_set) == 0:
        new_fetch_range = None
    else:
        new_fetch_range = range(
            (min(start0_set) - new_fetch_padding),
            #(max(start0_set) + 1 + new_fetch_padding),
            (max(end0_set) + new_fetch_padding),
        )

    return relevant_qname_set, new_fetch_range


def initial_fetch_sv(bam, bnds, view, fetch_padding_common,
                     fetch_padding_sv, fetch_padding_view,
                     new_fetch_padding, long_insert_threshold):
    """Read filtering done by readhandler.readfilter_bad_read"""

    def get_initial_fetcher(bam, chrom_bnd, pos_range0, endis5, 
                               fetch_padding_common, fetch_padding_sv):
        if endis5:
            fetcher = readhandler.get_fetch(
                bam, chrom,
                start=(
                    min(pos_range0)
                    - fetch_padding_common
                ),
                end=(
                    max(pos_range0) + 1
                    + fetch_padding_common
                    + fetch_padding_sv
                ),
                readfilter=readhandler.readfilter_bad_read)
        else:
            fetcher = readhandler.get_fetch(
                bam, chrom, 
                start=(
                    min(pos_range0)
                    - fetch_padding_common
                    - fetch_padding_sv
                ),
                end=(
                    max(pos_range0) + 1
                    + fetch_padding_common
                ),
                readfilter=readhandler.readfilter_bad_read)

        return fetcher

    def get_initial_fetcher_view(bam, chrom_bnd, pos_range0, 
                                 fetch_padding_view):
        fetcher = readhandler.get_fetch(
            bam, chrom_bnd,
            start=(min(pos_range0) - fetch_padding_view),
            end=(max(pos_range0) + 1 + fetch_padding_view),
            readfilter=readhandler.readfilter_bad_read,
        )

        return fetcher

    def readfilter_qname(read, bndborder_start0, bndborder_end0, endis5):
        """The read is selected if True.
        For decision of whether the read query_name will be included in the
            "relevant_qname_list", for SV.
        """

        cond_overlapping = (read.reference_start < bndborder_end0 and 
                            read.reference_end > bndborder_start0)
        if endis5:
            cond_distant = (read.is_reverse and 
                            read.reference_start >= bndborder_end0)
        else:
            cond_distant = (read.is_forward and 
                            read.reference_end <= bndborder_start0)

        result = (cond_overlapping or cond_distant)

        return result

    def handle_initial_fetcher(fetcher, view, pos_range0, endis5, 
                               long_insert_threshold):
        relevant_qname_set = set()
        start0_set = set()

        for read in fetcher:
            if view:
                read_is_relevant = True
            else:
                read_is_relevant = readfilter_qname(
                    read, min(pos_range0), max(pos_range0) + 1, endis5)

            if read_is_relevant:
                relevant_qname_set.add(read.query_name)
                start0_set.add(read.reference_start)
                if readfilter_matesearch(read, long_insert_threshold):
                    start0_set.add(read.next_reference_start)

        return relevant_qname_set, start0_set

    def get_new_fetch_range(start0_set, new_fetch_padding):
        if len(start0_set) == 0:
            new_fetch_range = None
        else:
            new_fetch_range = range(
                min(start0_set) - new_fetch_padding, 
                max(start0_set) + 1 + new_fetch_padding)

        return new_fetch_range

    # MAIN
    # get fetchers
    if view:
        fetcher_bnd1 = get_initial_fetcher_view(bam, bnds.chrom_bnd1, 
                                                bnds.get_pos_range0_bnd1(),
                                                fetch_padding_view)
        fetcher_bnd2 = get_initial_fetcher_view(bam, bnds.chrom_bnd2, 
                                                bnds.get_pos_range0_bnd2(),
                                                fetch_padding_view)
    else:
        fetcher_bnd1 = get_initial_fetcher(
            bam, bnds.chrom_bnd1, bnds.get_pos_range0_bnd1(), bnds.endis5_bnd1,
            fetch_padding_common, fetch_padding_sv)
        fetcher_bnd2 = get_initial_fetcher(
            bam, bnds.chrom_bnd2, bnds.get_pos_range0_bnd2(), bnds.endis5_bnd2,
            fetch_padding_common, fetch_padding_sv)

    # extract information from initial fetchers
    relevant_qname_set_bnd1, start0_set_bnd1 = handle_initial_fetcher(
        fetcher_bnd1, view, bnds.get_pos_range0_bnd1(), 
        bnds.endis5_bnd1, long_insert_threshold)
    relevant_qname_set_bnd2, start0_set_bnd2 = handle_initial_fetcher(
        fetcher_bnd2, view, bnds.get_pos_range0_bnd2(), 
        bnds.endis5_bnd2, long_insert_threshold)

    # get new fetch ranges
    new_fetch_range_bnd1 = get_new_fetch_range(start0_set_bnd1, 
                                               new_fetch_padding)
    new_fetch_range_bnd2 = get_new_fetch_range(start0_set_bnd2, 
                                               new_fetch_padding)

    # qname list union
    relevant_qname_set_union = set.union(relevant_qname_set_bnd1,
                                         relevant_qname_set_bnd2)

    return (relevant_qname_set_bnd1, new_fetch_range_bnd1,
            relevant_qname_set_bnd2, new_fetch_range_bnd2,
            relevant_qname_set_union)


#### refined fetch ####

def refined_fetch_nonsv(
    bam, 
    chrom, 
    new_fetch_range, 
    relevant_qname_set, 
    include_irrelevant_reads,
):
    """Read filtering is done by readfilter_bad_read.
    #Store the fetched read only if its qname is included in the
    #    "relevant_qname_set".
    """
    fetchresult_dict = dict()
    for read in readhandler.get_fetch(
        bam, chrom, start=new_fetch_range.start, 
        end=new_fetch_range.stop, 
        readfilter=readhandler.readfilter_bad_read,
    ):
        if read.query_name in relevant_qname_set:
            fetchresult_dict.setdefault(read.query_name, list())
            fetchresult_dict[read.query_name].append(read)
        else:
            if include_irrelevant_reads:
                fetchresult_dict.setdefault(read.query_name, list())
                fetchresult_dict[read.query_name].append(read)
            else:
                del read

        #if read.query_name in fetchresult_dict:
        #    fetchresult_dict[read.query_name].append(read)
        #else:
        #    del read

    return fetchresult_dict


def refined_fetch_sv(bam, bnds, new_fetch_range_bnd1, new_fetch_range_bnd2,
                     relevant_qname_set_union):
    """Read filtering is done by readfilter_bad_read.
    Store the fetched read only if its qname is included in the
        "relevant_qname_list".
    """

    def populate_fetchresult_dict(fetchresult_dict, new_fetch_range, bam,
                                  chrom_bnd):
        if new_fetch_range is not None:
            fetch = readhandler.get_fetch(
                bam, chrom_bnd, 
                start=new_fetch_range.start, 
                end=new_fetch_range.stop,
                readfilter=readhandler.readfilter_bad_read,
            )
            for read in fetch:
                if read.query_name in fetchresult_dict:
                    if all(
                        (read.compare(x) != 0)
                        for x in fetchresult_dict[read.query_name]
                    ):
                        fetchresult_dict[read.query_name].append(read)

    # init fetchresult_dict
    fetchresult_dict = dict()
    for qname in relevant_qname_set_union:
        fetchresult_dict[qname] = list()

    # populate fetchresult_dict
    populate_fetchresult_dict(fetchresult_dict, new_fetch_range_bnd1, bam,
                              bnds.chrom_bnd1)
    populate_fetchresult_dict(fetchresult_dict, new_fetch_range_bnd2, bam,
                              bnds.chrom_bnd2)

    return fetchresult_dict


def get_rpp_from_refinedfetch(readlist, bam, fasta, chromdict, no_matesearch, recalc_NMMD=False):
    # classify reads into primary and non-primary ones
    readlist_primary = list()
    readlist_nonprimary = list()
    for read in readlist:
        if read.reference_name not in chromdict.contigs:
            del read
            continue

        if readhandler.check_primary_alignment(read):
            readlist_primary.append(read)
        else:
            readlist_nonprimary.append(read)

    if len(readlist_primary) > 2:
        e_msg_list = list()
        e_msg_list.append('More than two primary reads found:')
        for read in readlist_primary:
            e_msg_list.append(read.to_string())
        raise Exception('\n'.join(e_msg_list))
    elif len(readlist_primary) == 0: 
        # found reads are all non-primary
        # supplementary reads overlapping pos0 can be missed
        rpp = None

        for read in readlist_nonprimary:
            del read
        del readlist_nonprimary
    else:
        rplist_nonprimary = [
            ReadPlus(x, fasta, recalc_NMMD=recalc_NMMD)
            for x in readlist_nonprimary
        ]
        if len(readlist_primary) == 1: # mate read is somewhere far away
            if no_matesearch:
                rplist_primary = [ReadPlus(readlist_primary[0], fasta, recalc_NMMD=recalc_NMMD)]
            else:
                mate = readhandler.get_primary_mate(readlist_primary[0], bam)
                if mate is None:
                    raise Exception(
                        f'Mate not found for this read:\n'
                        f'{readlist_primary[0].to_string()}')

                if mate.reference_name in chromdict.contigs:
                    rplist_primary = [
                        ReadPlus(readlist_primary[0], fasta, recalc_NMMD=recalc_NMMD), 
                        ReadPlus(mate, fasta, recalc_NMMD=recalc_NMMD)
                    ]
                else:
                    rplist_primary = [ReadPlus(readlist_primary[0], fasta, recalc_NMMD=recalc_NMMD)]

        else: # len(readlist_primary) == 2
            rplist_primary = [
                ReadPlus(x, fasta, recalc_NMMD=recalc_NMMD) 
                for x in readlist_primary
            ]

        rpp = ReadPlusPair(rplist_primary, rplist_nonprimary, chromdict)

    return rpp




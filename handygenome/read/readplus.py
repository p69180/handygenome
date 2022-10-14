import os
import itertools
import pprint
import collections
import functools

import numpy as np
import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
readhandler = importlib.import_module('.'.join([top_package_name, 'read', 'readhandler']))
alleleinfosetup = importlib.import_module('.'.join([top_package_name, 'read', 'alleleinfosetup']))
alleleinfosetup_sv = importlib.import_module('.'.join([top_package_name, 'read', 'alleleinfosetup_sv']))


RPPLIST_MAX_SHOW_LEN = 30
THRESHOLD_TEMPLATE_LENGTH = 1000

FETCH_PADDING_COMMON = 50  # In order to fetch clip-overlapping reads
FETCH_PADDING_SV = 850
FETCH_PADDING_VIEW = 1000
NEW_FETCH_PADDING = 0
LONG_INSERT_THRESHOLD = 1000

LOGGER_RPPLIST = workflow.get_logger(name='ReadPlusPairList', 
                                     level='warning')

ALLELEINFO_TAG_RPP = 'a2'
ALLELEINFO_TAG_RP = 'a1'


class ReadPlus:

    def __init__(self, read, fasta=None, skip_attrs_setting=False, recalc_NMMD=False):
        self.read = read
        if not skip_attrs_setting:
            self._init_attrs(fasta, recalc_NMMD)

    def _init_attrs(self, fasta, recalc_NMMD):
        self.fasta = (
            pysam.FastaFile(
                common.DEFAULT_FASTA_PATHS[
                    common.infer_refver(bamheader=self.read.header)
                ]
            )
            if fasta is None else 
            fasta
        )

        if recalc_NMMD:
            self._set_NMMD()
        else:
            if (not self.read.has_tag('NM')) or (not self.read.has_tag('MD')):
                self._set_NMMD()

        #self._set_SAlist()
        self.pairs_dict = readhandler.get_pairs_dict(self.read, self.fasta)
            # keys: 'querypos0', 'refpos0','refseq'

        self.fiveprime_end = readhandler.get_fiveprime_end(self.read)
        self.threeprime_end = readhandler.get_threeprime_end(self.read)
        self.cigarstats = self.read.get_cigar_stats()[0] 
            # length 11 array which contains value for M, I, D, N, S, ...
        self.ref_range0 = range(self.read.reference_start, 
                                self.read.reference_end)
        self.range0 = self.ref_range0  # alias
        self.softclip_range0 = None
        self.alleleinfo = dict()

    def __repr__(self):
        qname = self.read.query_name
        chrom = self.read.reference_name
        start1 = self.read.reference_start + 1
        end1 = self.read.reference_end
        region = f'{chrom}:{start1:,}-{end1:,}'
        infostring = f'{qname}; {region}; alleleinfo: {self.alleleinfo}'
        return (f'<ReadPlus object ({infostring})>')

    # publics
    @functools.cache
    def get_softclip_range0(self):
        return readhandler.get_softclip_ends_range0(self.read)
#        if self.softclip_range0 is None:
#            self.set_softclip_range0()
#        return self.softclip_range0

#    def set_softclip_range0(self):
#        self.softclip_range0 = readhandler.get_softclip_ends_range0(self.read)

    def get_BQlist(self, vcfspec):
        """For insertion: BQ of the base right before the insertion site
        For deletion: BQ of the base rigth before the deleted seq
        """
        BQlist = list()
        REF_range0 = vcfspec.REF_range0
        if self.check_spans(REF_range0):
            pairs_indexes = self.get_pairs_indexes(REF_range0)
            for idx in pairs_indexes:
                querypos0 = self.pairs_dict['querypos0'][idx]
                if querypos0 is not None:
                    BQlist.append(self.read.query_qualities[querypos0])

        return BQlist

    def check_overlaps(self, range0):
        return (self.read.reference_start <= max(range0) and
                self.read.reference_end > min(range0))

    def check_overlaps_vcfspec(self, vcfspec):
        return (self.read.reference_name == vcfspec.chrom and
                self.check_overlaps(vcfspec.REF_range0))
        
    def check_softclip_overlaps(self, range0):
        softclip_range0 = self.get_softclip_range0()
        return (softclip_range0.start <= max(range0) and
                softclip_range0.stop > min(range0))

    def check_softclip_overlaps_vcfspec(self, vcfspec):
        return (self.read.reference_name == vcfspec.chrom and
                self.check_softclip_overlaps(vcfspec.REF_range0))

    def get_num_cigarM(self):
        return self.read.get_cigar_stats()[0][0]
        
    def check_spans(self, range0):
        num_cigarM = self.get_num_cigarM()
        if num_cigarM == 0:
            return False
        else:
            if len(range0) == 0:
                return (range0.stop >= self.read.reference_start and
                        range0.stop < self.read.reference_end)
            else:
                return (self.read.reference_start <= min(range0) and
                        self.read.reference_end > max(range0))

    def get_distance(self, range0):
        if self.check_overlaps(range0):
            distance = 0
        else:
            if self.read.reference_end <= min(range0):
                distance = min(range0) - self.read.reference_end + 1
            elif self.read.reference_start > max(range0):
                distance = self.read.reference_start - max(range0)

        return distance

    def get_pairs_indexes(self, range0):
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

    def check_matches_with_pairs_indexes(self, pairs_indexes):
        for pairs_idx in pairs_indexes:
            querypos0 = self.pairs_dict['querypos0'][pairs_idx]
            refseq = self.pairs_dict['refseq'][pairs_idx]
            matches = ((querypos0 is not None) and 
                       (refseq is not None) and
                       refseq.isupper())
            if not matches:
                return False
        return True

    def check_matches(self, range0):
        pairs_indexes = self.get_pairs_indexes(range0)
        return self.check_matches_with_pairs_indexes(pairs_indexes)

    def check_spans_and_matches(self, range0):
        if self.check_spans(range0):
            pairs_indexes = self.get_pairs_indexes(range0)
            return self.check_matches_with_pairs_indexes(pairs_indexes)
        else:
            return False

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

    def get_seq_from_range0(self, range0):
        pairs_indexes = self.get_pairs_indexes(range0)
        return self.get_seq_from_pairs_indexes(pairs_indexes)

    def get_seq_from_pairs_indexes(self, pairs_indexes):
        pairs_slice = slice(pairs_indexes.start, pairs_indexes.stop, 
                            pairs_indexes.step)
        querypos0_list = [x for x in self.pairs_dict['querypos0'][pairs_slice]
                          if x is not None]

        return self.get_seq_from_querypos0_list(querypos0_list)
        
    def get_seq_from_querypos0_list(self, querypos0_list):
        if len(querypos0_list) == 0:
            return ''
        else:
            start = min(querypos0_list)
            stop = max(querypos0_list) + 1
            return self.read.query_sequence[slice(start, stop, 1)]

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
            result = {'left': querypos0_fromleft, 
                      'right': querypos0_fromright, 
                      '5prime': querypos0_from5,
                      '3prime': querypos0_from3,
                      'left_fraction': querypos0_fromleft / readlen,
                      'right_fraction': querypos0_fromright / readlen,
                      '5prime_fraction': querypos0_from5 / readlen,
                      '3prime_fraction': querypos0_from3 / readlen}

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

    # cigar
    def walk_cigar(self):
        return readhandler.walk_cigar(self.read.cigartuples, self.read.reference_start)

    def split_cigar(self, split_range0):
        return readhandler.split_cigar(self.read.cigartuples, self.read.reference_start, split_range0)

    # alleleinfo
    def update_alleleinfo(self, vcfspec, 
                          flanklen=alleleinfosetup.DEFAULT_FLANKLEN):
        # aiitem: AlleleInfoItem
        aiitem = alleleinfosetup.make_alleleinfoitem_readplus(
            vcfspec=vcfspec, rp=self, flanklen=flanklen)
        self.alleleinfo[vcfspec] = aiitem

    def set_alleleinfo_tag(self, vcfspec):
        alleleclass = self.alleleinfo[vcfspec]['alleleclass']
        #value = f'{vcfspec.get_id()}_{alleleclass}'
        value = str(alleleclass)

        self.read.set_tag(ALLELEINFO_TAG_RP, value, 'Z', replace=True)

    def update_alleleinfo_sv(
            self, bnds, 
            flanklen_parside=alleleinfosetup_sv.DEFAULT_FLANKLEN_PARSIDE,
            flanklen_bndside=alleleinfosetup_sv.DEFAULT_FLANKLEN_BNDSIDE):
        aiitem_bnd1 = alleleinfosetup_sv.make_alleleinfoitem_readplus(
            bnds=bnds, is_bnd1=True, rp=self, 
            flanklen_parside=flanklen_parside,
            flanklen_bndside=flanklen_bndside)
        aiitem_bnd2 = alleleinfosetup_sv.make_alleleinfoitem_readplus(
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
            k for (k, v) in self.alleleinfo[bnds]['bnd1'].items()
            if v]
        alleleclass_list_bnd2 = [
            k for (k, v) in self.alleleinfo[bnds]['bnd2'].items()
            if v]

        if len(alleleclass_list_bnd1) == 0 or len(alleleclass_list_bnd2) == 0:
            raise Exception(
                f'Invalid alleleinfoitem.\n'
                f'ReadPlusPair: {self}')
        else:
            alleleclass_bnd1 = '&'.join(alleleclass_list_bnd1)
            alleleclass_bnd2 = '&'.join(alleleclass_list_bnd2)

        #value = '_'.join([bnds.get_id(),
        #                  f'bnd1={alleleclass_bnd1}',
        #                  f'bnd2={alleleclass_bnd2}'])
        value = '_'.join([f'bnd1={alleleclass_bnd1}',
                          f'bnd2={alleleclass_bnd2}'])

        self.read.set_tag(ALLELEINFO_TAG_RP, value, 'Z', replace=True)

    # others
    def _set_NMMD(self):
        (ref_seq_padded, read_seq_padded) = readhandler.get_padded_seqs(self.read, self.fasta)
        self.read.set_tag(
            tag='NM', value_type='i',
            value=readhandler.get_NM(
                self.read, ref_seq_padded=ref_seq_padded, 
                read_seq_padded=read_seq_padded))
        self.read.set_tag(
            tag='MD', value_type='Z',
            value=readhandler.get_MD(
                self.read, ref_seq_padded=ref_seq_padded, 
                read_seq_padded=read_seq_padded))

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

    def __init__(self, rplist_primary, rplist_nonprimary, chromdict,
                 threshold_tlen=THRESHOLD_TEMPLATE_LENGTH):
        self._set_rp1_rp2(rplist_primary, chromdict)
        self.rplist_nonprimary = rplist_nonprimary

        # query_name
        self.query_name = self.rp1.read.query_name

        # mate_chroms_differ
        if self.rp2 is None:
            self.mate_chroms_differ = None
        else:
            self.mate_chroms_differ = (
                self.rp1.read.reference_name != self.rp2.read.reference_name)
        self.is_TRA = self.mate_chroms_differ # alias

        self.pairorient = readhandler.get_pairorient(self.rp1.read)
            # may be None when: mate unmapped, TLEN == 0
        self.tlen = (None 
                     if self.mate_chroms_differ else 
                     abs(self.rp1.read.template_length))
        self.template_length = self.tlen  # alias

        self.alleleinfo = dict()

        #self._set_is_proper_pair()
        #self._set_sv_supporting()
        #self.irrelevant = (self.rp1.irrelevant or self.rp2.irrelevant)

    # non-alleleinfo methods
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

    # alleleinfo-related methods
    def update_alleleinfo(self, vcfspec, 
                          flanklen=alleleinfosetup.DEFAULT_FLANKLEN):
        # rp1
        #if vcfspec not in self.rp1.alleleinfo:
        self.rp1.update_alleleinfo(vcfspec, flanklen=flanklen)
        aiitem_rp1 = self.rp1.alleleinfo[vcfspec]

        # rp2
        #if vcfspec not in self.rp2.alleleinfo:
        if self.rp2 is None:
            self.alleleinfo[vcfspec] = aiitem_rp1
        else:
            self.rp2.update_alleleinfo(vcfspec, flanklen=flanklen)
            aiitem_rp2 = self.rp2.alleleinfo[vcfspec]
            aiitem = alleleinfosetup.make_alleleinfoitem_readpluspair(
                vcfspec, aiitem_rp1, aiitem_rp2)
            self.alleleinfo[vcfspec] = aiitem

    def set_alleleinfo_tag(self, vcfspec):
        alleleclass = self.alleleinfo[vcfspec]['alleleclass']
        #value = f'{vcfspec.get_id()}_{alleleclass}'
        value = str(alleleclass)

        self.rp1.read.set_tag(ALLELEINFO_TAG_RPP, value, 'Z', replace=True)
        if self.rp2 is not None:
            self.rp2.read.set_tag(ALLELEINFO_TAG_RPP, value, 'Z', replace=True)

    def update_alleleinfo_sv(
            self, bnds, 
            flanklen_parside=alleleinfosetup_sv.DEFAULT_FLANKLEN_PARSIDE,
            flanklen_bndside=alleleinfosetup_sv.DEFAULT_FLANKLEN_BNDSIDE):
        self.rp1.update_alleleinfo_sv(bnds,
                                      flanklen_parside=flanklen_parside,
                                      flanklen_bndside=flanklen_bndside)
        self.rp2.update_alleleinfo_sv(bnds,
                                      flanklen_parside=flanklen_parside,
                                      flanklen_bndside=flanklen_bndside)
        aiitem = alleleinfosetup_sv.make_alleleinfoitem_readpluspair(
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

        self.rp1.read.set_tag(ALLELEINFO_TAG_RPP, value, 'Z', replace=True)
        self.rp2.read.set_tag(ALLELEINFO_TAG_RPP, value, 'Z', replace=True)

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
                      f'alleleinfo: {self.alleleinfo}')

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


# rpplist classes and functions

class ReadPlusPairList(list):

    def __init__(self, chromdict):
        self.chromdict = chromdict
        self.bam_path = None

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
                     < alleleinfosetup_sv.BND_DISTANCE_THRESHOLD))
                bnd2_relevant = (
                    rp.read.reference_name == bnds.chrom_bnd2 and
                    (rp.get_distance(bnds.get_pos_range0_bnd2())
                     < alleleinfosetup_sv.BND_DISTANCE_THRESHOLD))

                if bnd1_relevant:
                    start_list_bnd1.append(rp.read.reference_start)
                    end_list_bnd1.append(rp.read.reference_end)
                if bnd2_relevant:
                    start_list_bnd2.append(rp.read.reference_start)
                    end_list_bnd2.append(rp.read.reference_end)

        ref_range0_bnd1 = range(min(start_list_bnd1), max(end_list_bnd1))
        ref_range0_bnd2 = range(min(start_list_bnd2), max(end_list_bnd2))

        return ref_range0_bnd1, ref_range0_bnd2

    def update_alleleinfo(self, vcfspec, 
                          flanklen=alleleinfosetup.DEFAULT_FLANKLEN):
        for rpp in self:
            rpp.update_alleleinfo(vcfspec, flanklen=flanklen)

    def update_alleleinfo_sv(
            self, bnds, 
            flanklen_parside=alleleinfosetup_sv.DEFAULT_FLANKLEN_PARSIDE,
            flanklen_bndside=alleleinfosetup_sv.DEFAULT_FLANKLEN_BNDSIDE):
        for rpp in self:
            rpp.update_alleleinfo_sv(bnds, flanklen_parside=flanklen_parside,
                                     flanklen_bndside=flanklen_bndside)

    def set_alleleinfo_tag_rpp(self, vcfspec):
        for rpp in self:
            rpp.set_alleleinfo_tag(vcfspec)

    def set_alleleinfo_tag_rpp_sv(self, bnds):
        for rpp in self:
            rpp.set_alleleinfo_tag_sv(bnds)

    def set_alleleinfo_tag_rp(self, vcfspec):
        for rpp in self:
            rpp.rp1.set_alleleinfo_tag(vcfspec)
            if rpp.rp2 is not None:
                rpp.rp2.set_alleleinfo_tag(vcfspec)

    def set_alleleinfo_tag_rp_sv(self, bnds):
        for rpp in self:
            rpp.rp1.set_alleleinfo_tag_sv(bnds)
            rpp.rp2.set_alleleinfo_tag_sv(bnds)

    def get_readcounts(self, vcfspec):
        """Designed only for non-sv"""

        counter = collections.Counter(rpp.alleleinfo[vcfspec]['alleleclass']
                                      for rpp in self)
        readcounts = {key: counter[key]
                      for key in (None, -1, 0, 1)}

        return readcounts

    def write_bam(self, outfile_path=None, outfile_dir=None):
        # sanity check
        if len(self) == 0:
            raise Exception(f'The ReadPlusPairList is empty.')
        # set outfile_path 
        if outfile_path is None:
            if outfile_dir is None:
                outfile_path=workflow.get_tmpfile_path(suffix='.bam')
            else:
                outfile_path=workflow.get_tmpfile_path(suffix='.bam',
                                                       where=outfile_dir)
        # get a list of reads
        readlist = list()
        for rpp in self:
            readlist.append(rpp.rp1.read)
            if rpp.rp2 is not None:
                readlist.append(rpp.rp2.read)
        # sort
        sortkey = common.get_read_sortkey(self.chromdict)
        readlist.sort(key=sortkey)
        # write bam
        bamheader = self[0].rp1.read.header
        with pysam.AlignmentFile(outfile_path, mode='wb', 
                                 header=bamheader) as out_bam:
            for read in readlist:
                out_bam.write(read)
        # index
        pysam.index(outfile_path)

        self.bam_path = outfile_path

    def rm_bam(self):
        if self.bam_path is not None:
            os.remove(self.bam_path)
            bai_path = self.bam_path + '.bai'
            if os.path.exists(bai_path):
                os.remove(bai_path)
            self.bam_path = None


def get_rpplist_nonsv(bam, fasta, chromdict, chrom, start0, end0, 
                      view=False, no_matesearch=True,
                      fetch_padding_common=FETCH_PADDING_COMMON,
                      fetch_padding_view=FETCH_PADDING_VIEW,
                      new_fetch_padding=NEW_FETCH_PADDING,
                      long_insert_threshold=LONG_INSERT_THRESHOLD):
    LOGGER_RPPLIST.info('Beginning initial fetch')
    (relevant_qname_set, new_fetch_range) = initial_fetch_nonsv(
        bam, chrom, start0, end0, view,
        fetch_padding_common, fetch_padding_view,
        new_fetch_padding, long_insert_threshold)

    rpplist = ReadPlusPairList(chromdict=chromdict)
    if len(relevant_qname_set) > 0:
        LOGGER_RPPLIST.info('Beginning refined fetch')
        fetchresult_dict = refined_fetch_nonsv(bam, chrom, new_fetch_range, 
                                               relevant_qname_set)

        LOGGER_RPPLIST.info('Beginning assembly into readpluspair')
        for readlist in fetchresult_dict.values():
            rpp = get_rpp_from_refinedfetch(readlist, bam, fasta, chromdict,
                                            no_matesearch)
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


def initial_fetch_nonsv(bam, chrom, start0, end0, view,
                        fetch_padding_common, fetch_padding_view,
                        new_fetch_padding, long_insert_threshold):
    """Read filtering done by readhandler.readfilter_bad_read"""

    def readfilter_qname(read, start0, end0):
        """The read is selected if True.
        For decision of whether the read query_name will be included in the
            "relevant_qname_list", for non-SV.
        """

        softclip_range0 = readhandler.get_softclip_ends_range0(read)
        return (softclip_range0.stop > start0 and
                softclip_range0.start < end0)
        #return True

    relevant_qname_set = set()
    start0_set = set()

    initial_fetch_padding = (fetch_padding_view 
                             if view else 
                             fetch_padding_common)
    fetcher = readhandler.get_fetch(bam, chrom, 
                                    start=(start0 - initial_fetch_padding),
                                    end=(end0 + initial_fetch_padding),
                                    readfilter=readhandler.readfilter_bad_read)

    for read in fetcher:
        read_is_relevant = (True
                            if view else
                            readfilter_qname(read, start0, end0))
        if read_is_relevant:
            relevant_qname_set.add(read.query_name)
            start0_set.add(read.reference_start)
            if readfilter_matesearch(read, long_insert_threshold):
                start0_set.add(read.next_reference_start)

    if len(start0_set) == 0:
        new_fetch_range = None
    else:
        new_fetch_range = range(min(start0_set) - new_fetch_padding, 
                                max(start0_set) + 1 + new_fetch_padding)

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
                start=(min(pos_range0)
                       - fetch_padding_common),
                end=(max(pos_range0) + 1
                     + fetch_padding_common
                     + fetch_padding_sv),
                readfilter=readhandler.readfilter_bad_read)
        else:
            fetcher = readhandler.get_fetch(
                bam, chrom, 
                start=(min(pos_range0)
                       - fetch_padding_common
                       - fetch_padding_sv),
                end=(max(pos_range0) + 1
                     + fetch_padding_common),
                readfilter=readhandler.readfilter_bad_read)

        return fetcher

    def get_initial_fetcher_view(bam, chrom_bnd, pos_range0, 
                                 fetch_padding_view):
        fetcher = readhandler.get_fetch(
            bam, chrom_bnd,
            start=(min(pos_range0) - fetch_padding_view),
            end=(max(pos_range0) + 1 + fetch_padding_view),
            readfilter=readhandler.readfilter_bad_read)

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

def refined_fetch_nonsv(bam, chrom, new_fetch_range, relevant_qname_set):
    """Read filtering is done by readfilter_bad_read.
    Store the fetched read only if its qname is included in the
        "relevant_qname_set".
    """

    fetchresult_dict = dict()

    for qname in relevant_qname_set:
        fetchresult_dict[qname] = list()

    for read in readhandler.get_fetch(
            bam, chrom, start=new_fetch_range.start, 
            end=new_fetch_range.stop, 
            readfilter=readhandler.readfilter_bad_read):
        if read.query_name in fetchresult_dict:
            fetchresult_dict[read.query_name].append(read)

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
                start=new_fetch_range.start, end=new_fetch_range.stop,
                readfilter=readhandler.readfilter_bad_read)
            for read in fetch:
                if read.query_name in fetchresult_dict:
                    if all((read.compare(x) != 0)
                           for x in fetchresult_dict[read.query_name]):
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


def get_rpp_from_refinedfetch(readlist, bam, fasta, chromdict, no_matesearch):
    # classify reads into primary and non-primary ones
    readlist_primary = list()
    readlist_nonprimary = list()
    for read in readlist:
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
    else:
        rplist_nonprimary = [ReadPlus(x, fasta)
                             for x in readlist_nonprimary]
        if len(readlist_primary) == 1: # mate read is somewhere far away
            if no_matesearch:
                rplist_primary = [ReadPlus(readlist_primary[0], fasta)]
            else:
                mate = readhandler.get_primary_mate(readlist_primary[0], bam)
                if mate is None:
                    raise Exception(
                        f'Mate not found for this read:\n'
                        f'{readlist_primary[0].to_string()}')
                rplist_primary = [
                    ReadPlus(readlist_primary[0], fasta), 
                    ReadPlus(mate, fasta)]

        else: # len(readlist_primary) == 2
            rplist_primary = [ReadPlus(x, fasta) 
                              for x in readlist_primary]

        rpp = ReadPlusPair(rplist_primary, rplist_nonprimary, chromdict)

    return rpp



#class ReadPlusPairList_old(list):
#
#    @common.get_deco_arg_choices({'mode': ('sv', 'nonsv', 'view')})
#    def __init__(self, bam, chrom, start0, end0, fasta, chromdict, 
#                 mode, endis5=None,
#                 fetch_padding_common=FETCH_PADDING_COMMON,
#                 fetch_padding_sv=FETCH_PADDING_SV,
#                 long_insert_threshold=LONG_INSERT_THRESHOLD,
#                 new_fetch_padding=NEW_FETCH_PADDING,
#                 verbose=False):
#        """Read filtering is done by "readhandler.readfilter_bad_read".
#        Mate-unmapped reads are filtered out.
#        Args:
#            mode: must be one of "sv", "nonsv", "view"
#                sv: for SV breakends
#                nonsv: for non-SV
#                view: for creating a new bam with tags added
#        """
#
#        # sanitycheck
#        if mode == 'sv':
#            if endis5 is None:
#                raise Exception(f'"endis5" must be set if "mode" is "sv".')
#
#        LOGGER_RPPLIST.info('Beginning initial fetch')
#        (relevant_qname_list, new_fetch_range) = self._initial_fetch(
#            bam, chrom, start0, end0, mode, endis5, 
#            fetch_padding_common, fetch_padding_sv,
#            new_fetch_padding, long_insert_threshold)
#
#        if relevant_qname_list is None:
#            self.qname_list = None
#            self.fetch_range = None
#            self.rpplist = None
#            self.start = None
#            self.end = None
#        else:
#            self.qname_list = relevant_qname_list
#            self.fetch_range = new_fetch_range
#    
#            LOGGER_RPPLIST.info('Beginning refined fetch')
#            fetchresult_dict = self._refined_fetch(bam, chrom, 
#                                                   new_fetch_range, 
#                                                   relevant_qname_list)
#    
#            LOGGER_RPPLIST.info('Beginning assembly into readpluspair')
#            for readlist in fetchresult_dict.values():
#                rpp = self._get_readpluspair(readlist, bam, fasta, chromdict)
#                if rpp is not None:
#                    self.append(rpp)
#
#            # sort
#            self.sort(key=(lambda rpp: min(rpp.rp1.read.reference_start,
#                                           rpp.rp2.read.reference_start)))
#
#    def update_alleleinfo(self, vcfspec):
#        for rpp in self:
#            rpp.update_alleleinfo(vcfspec)
#
#    def _readfilter_qname_nonsv(self, read, start0, end0):
#        """For use in "_initial_fetch" function.
#        The read is selected if True.
#        For decision of whether the read query_name will be included in the
#            "relevant_qname_list", for non-SV.
#        """
#
#        #softclip_range0 = readhandler.get_softclip_ends_range0(read)
#        #return (softclip_range0.stop > start0 and
#        #        softclip_range0.start < end0)
#
#        return True
#
#    def _readfilter_qname_sv(self, read, start0, end0, endis5):
#        """For use in "_initial_fetch" function.
#        The read is selected if True.
#        For decision of whether the read query_name will be included in the
#            "relevant_qname_list", for SV.
#        """
#
#        cond_overlapping = (read.reference_start < end0 and 
#                            read.reference_end > start0)
#        if endis5:
#            cond_distant = (read.is_reverse and read.reference_start >= end0)
#        else:
#            cond_distant = (read.is_forward and read.reference_end <= start0)
#
#        result = (cond_overlapping or cond_distant)
#
#        return result
#
#    def _readfilter_matesearch(self, read, long_insert_threshold):
#        """For use in "_initial_fetch" function.
#        The read is selected if True.
#        For decision of whether mate start position will be included in the 
#            "new_fetch_range".
#        """
#
#        return (read.reference_name == read.next_reference_name and
#                abs(read.template_length) < long_insert_threshold)
#
#    def _get_initial_fetcher(self, bam, chrom, start0, end0, mode, endis5,
#                             fetch_padding_common, fetch_padding_sv):
#        if mode == 'sv':
#            if endis5:
#                fetcher = readhandler.get_fetch(
#                    bam, chrom, 
#                    start=(start0 - fetch_padding_common),
#                    end=(end0 
#                         + fetch_padding_common
#                         + fetch_padding_sv),
#                    readfilter=readhandler.readfilter_bad_read)
#            else:
#                fetcher = readhandler.get_fetch(
#                    bam, chrom, 
#                    start=(start0 
#                           - fetch_padding_common
#                           - fetch_padding_sv),
#                    end=(end0 + fetch_padding_common),
#                    readfilter=readhandler.readfilter_bad_read)
#        elif mode in ('nonsv', 'view'):
#            fetcher = readhandler.get_fetch(
#                bam, chrom, 
#                start=(start0 - fetch_padding_common),
#                end=(end0 + fetch_padding_common),
#                readfilter=readhandler.readfilter_bad_read)
#
#        return fetcher
#
#    def _initial_fetch(self, bam, chrom, start0, end0, mode, endis5,
#                       fetch_padding_common, fetch_padding_sv,
#                       new_fetch_padding, long_insert_threshold):
#        """Read filtering done by readhandler.readfilter_bad_read"""
#
#        relevant_qname_list = list()
#        start0_list = list()
#
#        # fetch start/end is padded with "initial_fetch_padding" to include 
#            # reads which overlap range(start0, end0) with soft-clipped 
#            # bases on IGV.
#
#        fetcher = self._get_initial_fetcher(
#            bam, chrom, start0, end0, mode, endis5,
#            fetch_padding_common, fetch_padding_sv)
#
#        for read in fetcher:
#            if mode == 'sv':
#                read_is_selected = self._readfilter_qname_sv(
#                    read, start0, end0, endis5)
#            elif mode == 'nonsv':
#                read_is_selected = self._readfilter_qname_nonsv(
#                    read, start0, end0)
#            elif mode == 'view':
#                read_is_selected = True
#
#            if read_is_selected:
#                relevant_qname_list.append(read.query_name)
#                start0_list.append(read.reference_start)
#                if self._readfilter_matesearch(read, long_insert_threshold):
#                    start0_list.append(read.next_reference_start)
#
#        if len(relevant_qname_list) == 0:
#            relevant_qname_list = None
#            new_fetch_range = None
#        else:
#            relevant_qname_list = list(set(relevant_qname_list))
#            new_fetch_range = range(min(start0_list) - new_fetch_padding, 
#                                    max(start0_list) + 1 + new_fetch_padding)
#
#        return relevant_qname_list, new_fetch_range
#
#    def _refined_fetch(self, bam, chrom, new_fetch_range, 
#                       relevant_qname_list):
#        """Read filtering is done by readfilter_bad_read.
#        Store the fetched read only if its qname is included in the
#            "relevant_qname_list".
#        """
#
#        fetchresult_dict = dict()
#
#        for qname in relevant_qname_list:
#            fetchresult_dict[qname] = list()
#
#        for read in readhandler.get_fetch(
#                bam, chrom, start=new_fetch_range.start, 
#                end=new_fetch_range.stop, 
#                readfilter=readhandler.readfilter_bad_read):
#            if read.query_name in fetchresult_dict:
#                fetchresult_dict[read.query_name].append(read)
#
#        return fetchresult_dict
#
#    def _get_readpluspair(self, readlist, bam, fasta, chromdict):
#        # classify reads into primary and non-primary ones
#        readlist_primary = list()
#        readlist_nonprimary = list()
#        for read in readlist:
#            if readhandler.check_primary_alignment(read):
#                readlist_primary.append(read)
#            else:
#                readlist_nonprimary.append(read)
#
#        if len(readlist_primary) > 2:
#            e_msg_list = list()
#            e_msg_list.append('More than two primary reads found:')
#            for read in readlist_primary:
#                e_msg_list.append(read.to_string())
#            raise Exception('\n'.join(e_msg_list))
#        elif len(readlist_primary) == 0: 
#            # found reads are all non-primary
#            # supplementary reads overlapping pos0 can be missed
#            rpp = None
#        else:
#            rplist_nonprimary = [ReadPlus(x, fasta)
#                                 for x in readlist_nonprimary]
#            if len(readlist_primary) == 1: # mate read is somewhere far away
#                mate = readhandler.get_primary_mate(readlist_primary[0], bam)
#                if mate is None:
#                    raise Exception(
#                        f'Mate not found for this read:\n'
#                        f'{readlist_primary[0].to_string()}')
#                rplist_primary = [
#                    ReadPlus(readlist_primary[0], fasta), 
#                    ReadPlus(mate, fasta)]
#
#            else: # len(readlist_primary) == 2
#                rplist_primary = [ReadPlus(x, fasta) 
#                                  for x in readlist_primary]
#
#            rpp = ReadPlusPair(rplist_primary, rplist_nonprimary, chromdict)
#
#        return rpp



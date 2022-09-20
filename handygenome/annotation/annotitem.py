import re
import pprint
import sys
import collections
import copy
import json
import functools

import pysam

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))

infoformat = importlib.import_module(
    ".".join([top_package_name, "variantplus", "infoformat"])
)
#varianthandler = importlib.import_module(
#    ".".join([top_package_name, "variantplus", "varianthandler"])
#)

# annotation = importlib.import_module('.'.join([top_package_name, 'annotation']))
#ensembl_parser = importlib.import_module(
#    ".".join([top_package_name, "annotation", "ensembl_parser"])
#)
#ensembl_rest = importlib.import_module(
#    ".".join([top_package_name, "annotation", "ensembl_rest"])
#)
#veplib = importlib.import_module(".".join([top_package_name, "annotation", "veplib"]))
#customfile = importlib.import_module(
#    ".".join([top_package_name, "annotation", "customfile"])
#)
#rnalib = importlib.import_module(".".join([top_package_name, "annotation", "rnalib"]))

#structvars = importlib.import_module(
#    ".".join([top_package_name, "svlib", "structvars"])
#)


# annotitem codec functions #

ENCODING_SEPS = ("\x1c", "\x1d")
DEFAULT_SEPS = (",", ":")
PERCENT_TABLE = {
    "%": "%25",
    ":": "%3A",
    ";": "%3B",
    "=": "%3D",
    ",": "%2C",
}
TR_TABLE_DECODER = str.maketrans(dict(zip(ENCODING_SEPS, DEFAULT_SEPS)))


def _replace(string, tup):
    return string.replace(*tup)


def _encode_replacer_percent(string):
    return functools.reduce(_replace, PERCENT_TABLE.items(), string)


def _decode_replacer_percent(string):
    return functools.reduce(
        _replace, (reversed(x) for x in PERCENT_TABLE.items()), string
    )


def _decode_replacer_seps(string):
    return string.translate(TR_TABLE_DECODER)


def _convert_key_to_nonstr(dic):
    result = dict()
    for k, v in dic.items():
        new_k = common.str_to_nonstr(k)
        if isinstance(v, dict):
            result[new_k] = _convert_key_to_nonstr(v)
        else:
            result[new_k] = v
    return result


def _encode(data, handle_percent=True):
    if handle_percent:
        return _encode_replacer_percent(json.dumps(data, separators=ENCODING_SEPS))
    else:
        return json.dumps(data, separators=ENCODING_SEPS)


def _decode(string, handle_percent=True, convert_key=True):
    if handle_percent:
        json_loaded = json.loads(_decode_replacer_seps(_decode_replacer_percent(string)))
    else:
        json_loaded = json.loads(_decode_replacer_seps(string))

    if convert_key:
        return _convert_key_to_nonstr(json_loaded)
    else:
        return json_loaded

#############################


# CONSTANTS
# SEPS = {
#    "nonprint": {"top": "\x1c", "keyval": "\x1d", "seq": "\x1e", "subkeyval": "\x1f"},
#    "initial": {"top": "$$$", "keyval": ":::", "seq": "&&&", "subkeyval": "@@@"},
#    "percent": {"top": "%3A", "keyval": "%3B", "seq": "%3D", "subkeyval": "%25"},
# }
#DEFAULT_SEPTYPE = "nonprint"

# SEP_TOP = '$$$'
# SEP_KEYVAL = ':::'
# SEP_SEQ = '&&&'
# SEP_SUBKEYVAL = '@@@'

# SEP_TOP = '%3A'
# SEP_KEYVAL = '%3B'
# SEP_SEQ = '%3D'
# SEP_SUBKEYVAL = '%25'

# SEP_TOP = '\x1c'
# SEP_KEYVAL = '\x1d'
# SEP_SEQ = '\x1e'
# SEP_SUBKEYVAL = '\x1f'

#ATOMIC_TYPES = (type(None), str, int, float, bool)


#class AnnotDB:
#    """
#    Attributes:
#        featureDB
#        transcript
#        transcript_canonical
#        regulatory
#        motif
#        repeat
#        popfreq
#        cosmic
#    """
#
#    def __init__(
#        self, annotdb_type, refver, fasta, chromdict, vr=None, septype=DEFAULT_SEPTYPE
#    ):
#        assert (
#            annotdb_type in ANNOTDB_TYPES
#        ), f"annotdb_type must be one of: {ANNOTDB_TYPES}"
#
#        self.annotdb_type = annotdb_type
#        self.refver = refver
#        # assert self.refver in ensembl_rest.REFVERS_ALLOWED
#        self.fasta = fasta
#        self.chromdict = chromdict
#        self.septype = septype
#
#        if vr is None:
#            self.init_empty()
#        else:
#            self.load_vr(vr)
#
#    def __repr__(self):
#        return common.cpformat(self.get_showdict())
#
#    def get_showdict(self):
#        showdict = dict()
#        showdict["popfreq"] = dict(self.popfreq)
#        showdict["cosmic"] = dict(self.cosmic)
#
#        showdict["transcript"] = self.transcript.get_showdict()
#        showdict["regulatory"] = self.regulatory.get_showdict()
#        showdict["motif"] = self.motif.get_showdict()
#        showdict["repeat"] = self.repeat.get_showdict()
#
#        return showdict
#
#    def init_empty(self):
#        self.transcript = AnnotItemDict(septype=self.septype)
#        self.regulatory = AnnotItemDict(septype=self.septype)
#        self.motif = AnnotItemDict(septype=self.septype)
#        self.repeat = AnnotItemDict(septype=self.septype)
#
#        self.set_transcript_subsets()
#
#        self.popfreq = Popfreq(septype=self.septype)
#        self.cosmic = Cosmic(septype=self.septype)
#
#    def load_vr(self, vr):
#        self.transcript = AnnotItemDict(
#            vr,
#            infokey=INFOMETAS[self.annotdb_type]["transcript"]["ID"],
#            septype=self.septype,
#        )
#        self.regulatory = AnnotItemDict(
#            vr,
#            infokey=INFOMETAS[self.annotdb_type]["regulatory"]["ID"],
#            septype=self.septype,
#        )
#        self.motif = AnnotItemDict(
#            vr,
#            infokey=INFOMETAS[self.annotdb_type]["motif"]["ID"],
#            septype=self.septype,
#        )
#        self.repeat = AnnotItemDict(
#            vr,
#            infokey=INFOMETAS[self.annotdb_type]["repeat"]["ID"],
#            septype=self.septype,
#        )
#
#        self.set_transcript_subsets()
#
#        self.popfreq = Popfreq(
#            vr,
#            infokey=INFOMETAS[self.annotdb_type]["popfreq"]["ID"],
#            septype=self.septype,
#        )
#        self.cosmic = Cosmic(
#            vr,
#            infokey=INFOMETAS[self.annotdb_type]["cosmic"]["ID"],
#            septype=self.septype,
#        )
#
#    def write(self, vr, add_meta=False):
#        if add_meta:
#            add_infometas(vr.header)
#
#        self.transcript.write_info(
#            vr, key=INFOMETAS[self.annotdb_type]["transcript"]["ID"], add_meta=False
#        )
#        self.regulatory.write_info(
#            vr, key=INFOMETAS[self.annotdb_type]["regulatory"]["ID"], add_meta=False
#        )
#        self.motif.write_info(
#            vr, key=INFOMETAS[self.annotdb_type]["motif"]["ID"], add_meta=False
#        )
#        self.repeat.write_info(
#            vr, key=INFOMETAS[self.annotdb_type]["repeat"]["ID"], add_meta=False
#        )
#
#        self.popfreq.write_info(
#            vr, key=INFOMETAS[self.annotdb_type]["popfreq"]["ID"], add_meta=False
#        )
#        self.cosmic.write_info(
#            vr, key=INFOMETAS[self.annotdb_type]["cosmic"]["ID"], add_meta=False
#        )
#
#    # POPFREQ ###################################################
#
#    def update_popfreq(self, vcfspec, dbsnp_vcf, search_equivs=True, overwrite=False):
#        assert (
#            self.annotdb_type == "plain"
#        ), f'This method is available only for annotdb_type "plain".'
#
#        try:
#            dbsnp_vr = fetch_dbsnp_vr(vcfspec, dbsnp_vcf, self.fasta, search_equivs)
#        except Exception as e:
#            raise Exception(f"dbsnp VCF fetch failed: " f"vcfspec: {vcfspec}") from e
#
#        self.popfreq.load_other(parse_dbsnp_vr(dbsnp_vr), overwrite=overwrite)
#
#    # COSMIC #####################################################
#
#    def update_cosmic(
#        self,
#        vcfspec,
#        cosmic_coding_vcf,
#        cosmic_noncoding_vcf,
#        search_equivs=True,
#        overwrite=False,
#    ):
#        assert (
#            self.annotdb_type == "plain"
#        ), f'This method is available only for annotdb_type "plain".'
#
#        cosmic_vr, cosmic_vr_noncoding = fetch_cosmic_vr(
#            vcfspec,
#            cosmic_coding_vcf,
#            cosmic_noncoding_vcf,
#            self.fasta,
#            search_equivs=search_equivs,
#        )
#        self.cosmic.load_other(
#            parse_cosmic_vr(cosmic_vr, cosmic_vr_noncoding), overwrite=overwrite
#        )
#
#    # FEATURES ###################################################
#
#    # post-cmdline-vep modifiers
#
#    def update_features_postvep_plain(
#        self,
#        vcfspec,
#        tabixfile_geneset,
#        tabixfile_regulatory,
#        tabixfile_repeats,
#        distance=veplib.DEFAULT_DISTANCE,
#    ):
#        """To be run with non-SV variants"""
#        assert self.annotdb_type == "plain"
#
#        mttype = vcfspec.get_mttype_firstalt()
#        assert mttype != "sv", f"Input variant must not be an SV."
#
#        # update using customfiles
#        fetch_interval = common.Interval(
#            chrom=vcfspec.chrom,
#            start1=max(1, vcfspec.pos - distance),
#            end1=min(self.chromdict[vcfspec.chrom], vcfspec.pos + distance),
#        )
#        self.update_customfile_transcript(
#            fetch_interval, tabixfile_geneset, fill_missing_is_canonical=True
#        )
#        self.update_customfile_regulatory(fetch_interval, tabixfile_regulatory)
#        self.update_customfile_repeat(fetch_interval, tabixfile_repeats)
#
#        # modifier
#        self.modifier_for_plain(vcfspec, mttype)
#        self.set_transcript_subsets()
#
#    def update_features_postvep_interval(
#        self,
#        interval,
#        tabixfile_geneset,
#        tabixfile_regulatory,
#        tabixfile_repeats,
#        distance=veplib.DEFAULT_DISTANCE,
#    ):
#        assert self.annotdb_type == "plain"
#
#        # update using customfiles
#        fetch_interval = common.Interval(
#            chrom=interval.chrom,
#            start1=max(1, interval.start1 - distance),
#            end1=min(self.chromdict[interval.chrom], interval.end1 + distance),
#        )
#        self.update_customfile_transcript(
#            fetch_interval, tabixfile_geneset, fill_missing_is_canonical=True
#        )
#        self.update_customfile_regulatory(fetch_interval, tabixfile_regulatory)
#        self.update_customfile_repeat(fetch_interval, tabixfile_repeats)
#
#        # modifier
#        self.modifier_for_interval(interval)
#        self.set_transcript_subsets()
#
#    def update_features_postvep_bnd(
#        self,
#        chrom,
#        pos,
#        endis5,
#        tabixfile_geneset,
#        tabixfile_regulatory,
#        tabixfile_repeats,
#        distance=veplib.DEFAULT_DISTANCE,
#    ):
#        assert self.annotdb_type in ("bnd1", "bnd2")
#
#        # update using customfiles
#        if endis5:
#            fetch_interval = common.Interval(
#                chrom=chrom, start1=pos, end1=min(self.chromdict[chrom], pos + distance)
#            )
#        else:
#            fetch_interval = common.Interval(
#                chrom=chrom, start1=max(1, pos - distance), end1=pos
#            )
#
#        self.update_customfile_transcript(
#            fetch_interval, tabixfile_geneset, fill_missing_is_canonical=True
#        )
#        self.update_customfile_regulatory(fetch_interval, tabixfile_regulatory)
#        self.update_customfile_repeat(fetch_interval, tabixfile_repeats)
#
#        # modifier
#        self.modifier_for_bnd(pos, endis5)
#        self.set_transcript_subsets()
#
#    # ensembl rest wrappers
#
#    def update_features_ensembl_wrapper_plain(
#        self,
#        vcfspec,
#        distance=veplib.DEFAULT_DISTANCE,
#        vep=False,
#        overlap=False,
#        transcript=False,
#        regulatory=False,
#    ):
#        """To be run with non-SV variants"""
#
#        assert self.annotdb_type == "plain"
#        mttype = vcfspec.get_mttype_firstalt()
#        assert mttype != "sv", f"Input variant must not be an SV."
#
#        if vep:
#            self.update_ensembl_rest_vep_plain(vcfspec, distance)
#        if overlap:
#            self.update_ensembl_rest_overlap_plain(vcfspec, distance)
#        if transcript:
#            self.update_ensembl_rest_lookup_transcript()
#        if regulatory:
#            self.update_ensembl_rest_regulatory()
#
#        self.modifier_for_plain(vcfspec, mttype)
#        self.set_transcript_subsets()
#
#    def update_features_ensembl_wrapper_interval(
#        self,
#        interval,
#        distance=veplib.DEFAULT_DISTANCE,
#        vep=False,
#        overlap=False,
#        transcript=False,
#        regulatory=False,
#    ):
#        assert self.annotdb_type == "plain"
#
#        if vep:
#            self.update_ensembl_rest_vep_interval(interval, distance)
#        if overlap:
#            self.update_ensembl_rest_overlap_interval(interval, distance)
#        if transcript:
#            self.update_ensembl_rest_lookup_transcript()
#        if regulatory:
#            self.update_ensembl_rest_regulatory()
#
#        self.modifier_for_interval(interval)
#        self.set_transcript_subsets()
#
#    def update_features_bnd(
#        self,
#        chrom,
#        pos,
#        endis5,
#        distance=veplib.DEFAULT_DISTANCE,
#        vep=False,
#        overlap=False,
#        transcript=False,
#        regulatory=False,
#    ):
#        assert self.annotdb_type in ("bnd1", "bnd2")
#
#        if vep:
#            self.update_ensembl_rest_vep_bnd(chrom, pos, distance)
#        if overlap:
#            self.update_ensembl_rest_overlap_bnd(chrom, pos, endis5, distance)
#        if transcript:
#            self.update_ensembl_rest_lookup_transcript()
#        if regulatory:
#            self.update_ensembl_rest_regulatory()
#
#        self.modifier_for_bnd(pos, endis5)
#        self.set_transcript_subsets()
#
#    # FEATURE UPDATE MODIFIERS ##########################################
#
#    # for setting missing distances
#
#    def set_missing_distances(self, distance_setter):
#        """Sets distance for annotitems not derived from VEP (e.g. repeat
#        elements). Simulates the distance given by VEP.
#        """
#
#        for attrname in ANNOTDB_ATTRNAMES_FEATURES:
#            annotitemdict = getattr(self, attrname)
#            for key, annotitem in annotitemdict.items():
#                if "distance" not in annotitem:
#                    # annotitems derived from "ensembl rest overlap"
#                    #   do not have "distance" value
#                    assert "start0" in annotitem and "end0" in annotitem, (
#                        f'"start0" and "end0" values are not set for this '
#                        f"AnnotItem object:\n{annotitem}"
#                    )
#                    feature_range0 = range(annotitem["start0"], annotitem["end0"])
#                    annotitem["distance"] = distance_setter(feature_range0)
#
#    def get_distance_setter_ins(self, vcfspec):
#        pos0 = vcfspec.pos - 1
#
#        def distance_setter(feature_range0):
#            if pos0 <= feature_range0.start - 1:
#                distance = feature_range0.start - pos0 - 1
#            elif pos0 >= feature_range0.stop - 1:
#                distance = feature_range0.stop - 1 - pos0
#            else:
#                distance = None
#
#            return distance
#
#        return distance_setter
#
#    def get_distance_setter_snv(self, vcfspec):
#        pos0 = vcfspec.pos - 1
#
#        def distance_setter(feature_range0):
#            if pos0 < feature_range0.start:
#                distance = feature_range0.start - pos0
#            elif pos0 >= feature_range0.stop:
#                distance = feature_range0.stop - 1 - pos0
#            else:
#                distance = None
#
#            return distance
#
#        return distance_setter
#
#    def get_distance_setter_interval(self, var_range0):
#        def distance_setter(feature_range0):
#            if var_range0.stop <= feature_range0.start:
#                distance = feature_range0.start - var_range0.stop + 1
#            elif var_range0.start >= feature_range0.stop:
#                distance = feature_range0.stop - 1 - var_range0.start
#            else:
#                distance = None
#
#            return distance
#
#        return distance_setter
#
#    def get_distance_setter_del(self, vcfspec):
#        var_range0 = range(vcfspec.pos, vcfspec.pos + len(vcfspec.ref) - 1)
#        return self.get_distance_setter_interval(var_range0)
#
#    def get_distance_setter_bnd(self, pos):
#        vcfspec = common.Vcfspec(None, pos, None, None)
#        return self.get_distance_setter_snv(vcfspec)
#
#    # modifier for plain
#    def modifier_for_plain(self, vcfspec, mttype):
#        if mttype == "snv":
#            distance_setter = self.get_distance_setter_snv(vcfspec)
#        elif mttype == "ins":
#            distance_setter = self.get_distance_setter_ins(vcfspec)
#        elif mttype in ("del", "delins", "mnv", "cpgmet"):
#            distance_setter = self.get_distance_setter_del(vcfspec)
#        else:
#            raise Exception(f"Unexpected mttype: {mttype}")
#
#        self.set_missing_distances(distance_setter)
#
#    # modifier for interval
#    def modifier_for_interval(self, interval):
#        distance_setter = self.get_distance_setter_interval(interval.range0)
#        self.set_missing_distances(distance_setter)
#        self.set_is_enclosed(interval.range0)
#
#    def set_is_enclosed(self, var_range0):
#        """
#        Sets "is_enclosed" attribute. Must be run after "start0" and "end0"
#        values have been set for all annotitems.
#        """
#
#        for attrname in ANNOTDB_ATTRNAMES_FEATURES:
#            annotitemdict = getattr(self, attrname)
#            for key, annotitem in annotitemdict.items():
#                feature_range0 = range(annotitem["start0"], annotitem["end0"])
#                annotitem["is_enclosed"] = (
#                    feature_range0.start >= var_range0.start
#                    and feature_range0.stop <= var_range0.stop
#                )
#
#    # modifier for breakend
#    def modifier_for_bnd(self, pos, endis5):
#        distance_setter = self.get_distance_setter_bnd(pos)
#        self.set_missing_distances(distance_setter)
#
#        self.settle_border_issues(pos, endis5)
#
#    def settle_border_issues(self, pos, endis5):
#        """
#        Must be run after "set_missing_distances" and "lookup_transcript" or
#        "lookup_regulatory" has been run.
#
#        1) Removes features entirely on the other side
#        2) For features whose border coincide with the breakend border,
#        distance is changed from None to 0.
#        3) Sets "is_broken" attribute
#        """
#
#        (on_the_other_side, on_the_border) = self.get_feature_classifiers(endis5)
#
#        for attrname in ANNOTDB_ATTRNAMES_FEATURES:
#            annotitemdict = getattr(self, attrname)
#            # step 1) and 2)
#            keys_to_delete = list()
#            for key, annotitem in annotitemdict.items():
#                if on_the_other_side(annotitem):
#                    keys_to_delete.append(key)
#                else:
#                    if annotitem["distance"] is None:
#                        if on_the_border(annotitem, pos):
#                            annotitem["distance"] = 0
#            for key in keys_to_delete:
#                del annotitemdict[key]
#            # step 3)
#            for key, annotitem in annotitemdict.items():
#                annotitem["is_broken"] = annotitem["distance"] is None
#
#    def get_feature_classifiers(self, endis5):
#        """
#                (endis5 == True)
#                                   THE OTHER SIDE    THE BREAKEND SIDE
#                                   - - - - - - - -  -----------------
#        FEATURE ON THE OTHER SIDE  = = = = = =     | = = = = = = = =  : FEATURE ON THE BORDER
#                                   - - - - - - - -  -----------------
#                                   A C G T A C G T   A C G T A C G T
#                                                     *(pos)
#                ("= = =" denotes a feature)
#        """
#
#        if endis5:
#
#            def on_the_other_side(annotitem):
#                return (annotitem["distance"] is not None) and (
#                    annotitem["distance"] < 0
#                )
#
#            def on_the_border(annotitem, pos):
#                return annotitem["start1"] == pos
#
#        else:
#
#            def on_the_other_side(annotitem):
#                return (annotitem["distance"] is not None) and (
#                    annotitem["distance"] > 0
#                )
#
#            def on_the_border(annotitem, pos):
#                return annotitem["end1"] == pos
#
#        return on_the_other_side, on_the_border
#
#    # CUSTOM FILE FETCH FUNCTIONS #####################################
#
#    def update_customfile_repeat(self, fetch_interval, tabixfile_repeats):
#        repeat_dict = customfile.fetch_repeat(
#            fetch_interval.chrom,
#            fetch_interval.start0,
#            fetch_interval.end0,
#            tabixfile_repeats,
#        )
#        self.repeat.load_other(repeat_dict, overwrite=False, create_new=True)
#
#    def update_customfile_transcript(
#        self, fetch_interval, tabixfile_geneset, fill_missing_is_canonical=True
#    ):
#        transcript_dict = customfile.fetch_transcript(
#            fetch_interval.chrom,
#            fetch_interval.start0,
#            fetch_interval.end0,
#            tabixfile_geneset,
#        )
#        self.transcript.load_other(transcript_dict, overwrite=False, create_new=True)
#        if fill_missing_is_canonical:
#            for annotitem_id, annotitem in self.transcript.items():
#                if "is_canonical" not in annotitem:
#                    raw_result = ensembl_rest.lookup_id(
#                        ID=annotitem["id"], refver=self.refver, expand=False
#                    )
#                    if raw_result["is_canonical"] == 1:
#                        annotitem["is_canonical"] = True
#                    else:
#                        annotitem["is_canonical"] = False
#
#    def update_customfile_regulatory(self, fetch_interval, tabixfile_regulatory):
#        regulatory_dict = customfile.fetch_regulatory(
#            fetch_interval.chrom,
#            fetch_interval.start0,
#            fetch_interval.end0,
#            tabixfile_regulatory,
#        )
#        self.regulatory.load_other(regulatory_dict, overwrite=False, create_new=True)
#
#    # ENSEMBL UNIT FUNCTIONS ##########################################
#
#    # vep
#    def update_ensembl_rest_vep_plain(self, vcfspec, distance=veplib.DEFAULT_DISTANCE):
#        self._update_ensembl_rest_vep_helper(
#            vcfspec=vcfspec, hgvsg=None, distance=distance
#        )
#
#    def update_ensembl_rest_vep_interval(
#        self, interval, distance=veplib.DEFAULT_DISTANCE
#    ):
#        sv_del = structvars.Deletion(
#            interval.chrom, interval.start1, interval.end1, fasta=self.fasta
#        )
#        hgvsg = sv_del.get_hgvsg()
#        self._update_ensembl_rest_vep_helper(
#            vcfspec=None, hgvsg=hgvsg, distance=distance
#        )
#
#    def update_ensembl_rest_vep_bnd(self, chrom, pos, distance=veplib.DEFAULT_DISTANCE):
#        ref = self.fasta.fetch(chrom, pos - 1, pos)
#        vcfspec = common.Vcfspec(chrom, pos, ref, ["N"])
#        self._update_ensembl_rest_vep_helper(
#            vcfspec=vcfspec, hgvsg=None, distance=distance
#        )
#
#    # overlap
#    def update_ensembl_rest_overlap_plain(
#        self, vcfspec, distance=veplib.DEFAULT_DISTANCE
#    ):
#        chrom = vcfspec.chrom
#        start1 = vcfspec.pos - distance
#        end1 = vcfspec.pos + distance
#        interval = common.Interval(chrom, start1, end1)
#        self._update_ensembl_rest_overlap_helper(interval)
#
#    def update_ensembl_rest_overlap_interval(
#        self, interval, distance=veplib.DEFAULT_DISTANCE
#    ):
#        new_interval = common.Interval(
#            chrom=interval.chrom,
#            start1=max(1, interval.start1 - distance),
#            end1=min(self.chromdict[interval.chrom], interval.end1 + distance),
#        )
#        self._update_ensembl_rest_overlap_helper(new_interval)
#
#    def update_ensembl_rest_overlap_bnd(
#        self, chrom, pos, endis5, distance=veplib.DEFAULT_DISTANCE
#    ):
#        if endis5:
#            start1 = pos
#            end1 = min(self.chromdict[chrom], pos + distance)
#        else:
#            start1 = max(1, pos - distance)
#            end1 = pos
#        interval = common.Interval(chrom, start1, end1)
#        self._update_ensembl_rest_overlap_helper(interval)
#
#    # lookup_transcript
#    def update_ensembl_rest_lookup_transcript(self):
#        parsed = ensembl_parser.parse_rest_lookup_transcript_post(
#            ensembl_rest.lookup_id_post(
#                tuple(self.transcript.keys()), refver=self.refver, expand=True
#            ),
#            refver=self.refver,
#            set_gene_name=True,
#        )
#        self._load_ensembl_parsed(parsed, overwrite=False, create_new=False)
#
#    # regulatory
#    def update_ensembl_rest_regulatory(self):
#        for ID in self.regulatory.keys():
#            self._update_ensembl_rest_regulatory_helper(ID)
#
#    # cmdline VEP
#    def update_cmdline_vep(self, vr, overwrite=True, create_new=True):
#        try:
#            parsed = ensembl_parser.parse_cmdline_vep(vr)
#        except Exception as e:
#            raise Exception(
#                f"Failure of VEP output variant record parsing:" f"\n{vr}"
#            ) from e
#
#        self._load_ensembl_parsed(parsed, overwrite=overwrite, create_new=create_new)
#
#    # helpers
#    def _load_ensembl_parsed(self, parsed, overwrite=False, create_new=False):
#        for annotdb_attrname, other_annotitemdict in parsed.items():
#            annotitemdict = getattr(self, annotdb_attrname)
#            annotitemdict.load_other(
#                other_annotitemdict, overwrite=overwrite, create_new=create_new
#            )
#
#    def _update_ensembl_rest_vep_helper(self, vcfspec, hgvsg, distance):
#        parsed = ensembl_parser.parse_rest_vep(
#            ensembl_rest.vep(
#                refver=self.refver,
#                vcfspec=vcfspec,
#                hgvsg=hgvsg,
#                distance=distance,
#                with_CADD=True,
#                with_Phenotypes=False,
#                with_canonical=True,
#                with_mane=True,
#                with_miRNA=False,
#                with_numbers=True,
#                with_protein=True,
#                with_ccds=True,
#                with_hgvs=True,
#            )
#        )
#        self._load_ensembl_parsed(parsed, overwrite=False, create_new=True)
#
#    def _update_ensembl_rest_overlap_helper(self, interval):
#        parsed = ensembl_parser.parse_rest_overlap(
#            ensembl_rest.overlap(
#                chrom=interval.chrom,
#                start1=interval.start1,
#                end1=interval.end1,
#                refver=self.refver,
#                transcript=False,
#                regulatory=False,
#                motif=False,
#                repeat=True,
#            ),
#            refver=self.refver,
#            include_motif_without_evidence=False,
#        )
#        self._load_ensembl_parsed(parsed, overwrite=False, create_new=True)
#
#    def _update_ensembl_rest_lookup_transcript_helper(self, ID):
#        parsed = ensembl_parser.parse_rest_lookup_transcript(
#            ensembl_rest.lookup_id(ID, refver=self.refver, expand=True),
#            refver=self.refver,
#            set_gene_name=True,
#        )
#        self._load_ensembl_parsed(parsed, overwrite=False, create_new=False)
#
#    def _update_ensembl_rest_regulatory_helper(self, ID):
#        parsed = ensembl_parser.parse_rest_regulatory(
#            ensembl_rest.regulatory(ID, refver=self.refver)
#        )
#        self._load_ensembl_parsed(parsed, overwrite=False, create_new=False)
#
#    # FEATURE ADDITIONAL ATTRIBUTES GENERATION #####################
#
#    def set_transcript_subsets(self):
#        # canonical, overlap, both
#        self.transcript_canon = AnnotItemDict()
#        self.transcript_ovlp = AnnotItemDict()
#        self.transcript_canon_ovlp = AnnotItemDict()
#
#        for ID, annotitem in self.transcript.items():
#            if annotitem["is_canonical"]:
#                self.transcript_canon[ID] = annotitem
#            if annotitem["distance"] is None:
#                self.transcript_ovlp[ID] = annotitem
#            if annotitem["is_canonical"] and (annotitem["distance"] is None):
#                self.transcript_canon_ovlp[ID] = annotitem
#
#    def set_exon_intron_ranges(self):
#        def make_exon_ranges(annotitem):
#            exon_ranges = dict()
#            for key, val in annotitem["exons"].items():
#                keysp = key.split("_")
#                valsp = [int(x) for x in val.split("_")]
#                exon_number = int(keysp[0])
#                exon_ranges[exon_number] = range(valsp[0], valsp[1])
#            return exon_ranges
#
#        def make_intron_ranges(exon_ranges, annotitem):
#            if len(exon_ranges) == 1:
#                intron_ranges = None
#            else:
#                intron_ranges = dict()
#                ascending = sorted(exon_ranges.items(), key=lambda x: x[1][0])
#                for idx in range(len(ascending) - 1):
#                    intron_number = (
#                        ascending[idx][0]
#                        if annotitem["is_forward"]
#                        else ascending[idx + 1][0]
#                    )
#                    range_item = range(
#                        ascending[idx][1][-1] + 1, ascending[idx + 1][1][0]
#                    )
#                    intron_ranges[intron_number] = range_item
#                intron_ranges = dict(sorted(intron_ranges.items(), key=lambda x: x[0]))
#
#            return intron_ranges
#
#        for ID, annotitem in self.transcript.items():
#            if "exons" in annotitem.keys():
#                annotitem.exon_ranges = make_exon_ranges(annotitem)
#                annotitem.intron_ranges = make_intron_ranges(
#                    annotitem.exon_ranges, annotitem
#                )


class AnnotItemBaseAll:
    def __repr__(self):
        return common.cpformat(self.get_self_show(), sort_dicts=False)

#    @property
#    def annotkey(self):
#        return self.__class__.meta["ID"]

    @staticmethod
    def decode(infostring):
        return _decode(infostring)

    @classmethod
    def get_annotkey(cls):
        return cls.meta["ID"]

#    @classmethod
#    def add_meta_info(cls, vcfheader):
#        vcfheader.add_meta(key="INFO", items=cls.meta.items())
#
#    @classmethod
#    def add_meta_format(cls, vcfheader):
#        vcfheader.add_meta(key="FORMAT", items=cls.meta.items())

#    def write_header(self, vcfheader):
#        vcfheader.add_meta(key=self.annotkey, value=self.encode())

#    def write_info(self, vr, key=None, add_meta=False):
#        if key is None:
#            key = self.annotkey
#        if add_meta:
#            self.__class__.add_meta_info(vr.header)
#
#        infoformat.set_info(vr, key, self.encode())
#
#    def write_format(self, vr, sampleid, key=None, add_meta=False):
#        if key is None:
#            key = self.annotkey
#        if add_meta:
#            self.__class__.add_meta_format(vr.header)
#
#        infoformat.set_format(vr, sampleid, key, self.encode())


class AnnotItemBaseHeader:
    def write(self, vcfheader):
        vcfheader.add_meta(key=self.get_annotkey(), value=self.encode())


class AnnotItemBaseInfo:
    @classmethod
    def add_meta(cls, vcfheader):
        vcfheader.add_meta(key="INFO", items=cls.meta.items())

    def write(self, vr, key=None, add_meta=False):
        if key is None:
            key = self.get_annotkey()
        if add_meta:
            self.__class__.add_meta(vr.header)

        infoformat.set_info(vr, key, self.encode())


class AnnotItemBaseFormat:
    @classmethod
    def add_meta(cls, vcfheader):
        vcfheader.add_meta(key="FORMAT", items=cls.meta.items())

    def write(self, vr, sampleid, key=None, add_meta=False):
        if key is None:
            key = self.get_annotkey()
        if add_meta:
            self.__class__.add_meta(vr.header)

        infoformat.set_format(vr, sampleid, key, self.encode())


class AnnotItemSingle(AnnotItemBaseAll, dict):
    ##### constructors #####
    @classmethod
    def from_annotstring(cls, annotstring):
        result = cls()
        result.update(cls.decode(annotstring))
        return result

    from_infostring = from_annotstring

    @classmethod
    def from_dict(cls, data):
        return cls(data)
    ########################

    @classmethod
    def annotstring_parser(cls, annotstring):
        return cls.from_annotstring(annotstring)

    def get_self_show(self):
        return dict(self)

    def encode(self):
        return _encode(dict(self))

    def update_annotstring(self, annotstring, overwrite=True):
        self.update_dict(self.__class__.decode(annotstring), overwrite=overwrite)

    update_infostring = update_annotstring
    load_infostring = update_annotstring

    def update_dict(self, other, overwrite=False):
        if overwrite:
            for key, val in other.items():
                self[key] = val
        else:
            for key, new_val in other.items():
                try:
                    old_val = self[key]
                except KeyError:
                    self[key] = new_val
                else:
                    if old_val is None:
                        self[key] = new_val
                    else:
                        if new_val != old_val:
                            self._update_conflict_handler(key, old_val, new_val)

    load_other = update_dict

    def init_missing(cls):
        return None

    def _update_conflict_handler(self, key, old_val, new_val):
        values = {'3prime_overlapping_ncRNA', '3prime_overlapping_ncrna'}
        if (
            key == 'biotype' and 
            old_val in values and 
            new_val in values
        ):
            self[key] = '3prime_overlapping_ncRNA'
        else:
            raise Exception(f"Old value and new value are different; key: {key}; old_val: {old_val}; new_val: {new_val}")


class AnnotItemALTlist(AnnotItemBaseAll, list):
    #meta = None
    def get_self_show(self):
        return [(None if (x is None) else x.get_self_show()) for x in self]

    def encode(self):
        return tuple((None if (x is None) else x.encode()) for x in self)

#########

class AnnotItemHeader(AnnotItemSingle, AnnotItemBaseHeader):
    @classmethod
    def from_vcfheader(cls, vcfheader):
        """Returns:
            None if the vcfheader does not contain annotkey
        """
        result = None
        annotkey = cls.get_annotkey()
        for rec in vcfheader.records:
            if rec.key == annotkey:
                result = cls.from_annotstring(rec.value)

        return result


class AnnotItemVariantInfoSingle(AnnotItemSingle, AnnotItemBaseInfo):
    @classmethod
    def from_vr_base(cls, vr, annotstring_parser_kwargs=dict()):
        if infoformat.check_NA_info(vr, cls.get_annotkey()):
            return None
        else:
            annotstring = vr.info[cls.get_annotkey()]
            return cls.annotstring_parser(annotstring, **annotstring_parser_kwargs)


class AnnotItemVariantFormatSingle(AnnotItemSingle, AnnotItemBaseFormat):
    @classmethod
    def from_vr_base(cls, vr, sampleid, annotstring_parser_kwargs=dict()):
        if infoformat.check_NA_format(vr, sampleid, cls.get_annotkey()):
            return None
        else:
            annotstring = vr.samples[sampleid][cls.get_annotkey()]
            return cls.annotstring_parser(annotstring, **annotstring_parser_kwargs)


class AnnotItemVariantInfoALTlist(AnnotItemALTlist, AnnotItemBaseInfo):
    @classmethod
    def init_missing(cls, vr):
        annotkey = cls.get_annotkey()
        if annotkey in vr.header.info.keys():
            result = cls()
            transl_num = infoformat.get_translated_number_info(vr, annotkey)
            for _ in range(transl_num):
                result.append(None)
            return result
        else:
            return None

    @classmethod
    def from_vr_base(cls, vr, annotstring_parser_kwargs=dict()):
        if infoformat.check_NA_info(vr, cls.get_annotkey()):
            result = cls.init_missing(vr)
        else:
            result = cls()
            for annotstring in vr.info[cls.get_annotkey()]:
                if annotstring in infoformat.NA_VALUES:
                    result.append(None)
                else:
                    result.append(cls.annotstring_parser(annotstring, **annotstring_parser_kwargs))
        return result


class AnnotItemVariantFormatALTList(AnnotItemALTlist, AnnotItemBaseFormat):
    @classmethod
    def init_missing(cls, vr):
        annotkey = cls.get_annotkey()
        if annotkey in vr.header.formats.keys():
            result = cls()
            transl_num = infoformat.get_translated_number_format(vr, cls.get_annotkey())
            for _ in range(transl_num):
                result.append(None)
            return result
        else:
            return None

    @classmethod
    def from_vr_base(cls, vr, sampleid, annotstring_parser_kwargs=dict()):
        if infoformat.check_NA_format(vr, sampleid, cls.get_annotkey()):
            result = cls.init_missing(vr)
        else:
            result = cls()
            for annotstring in vr.samples[sampleid][cls.get_annotkey()]:
                if annotstring in infoformat.NA_VALUES:
                    result.append(None)
                else:
                    result.append(cls.annotstring_parser(annotstring, **annotstring_parser_kwargs))
        return result


##########################


#class AnnotItem(AnnotItemBase, dict):
#    meta = None
#    # to be used as "items" argument for pysam.VariantHeader.add_meta
#
#    ##### constructors #####
#    @classmethod
#    def from_infostring(cls, infostring):
#        result = cls()
#        result.update(cls.decode(infostring))
#        return result
#
#    @classmethod
#    def from_vcfheader(cls, vcfheader):
#        result = None
#        annotkey = cls.get_annotkey()
#        for rec in vcfheader.records:
#            if rec.key == annotkey:
#                result = cls.from_infostring(rec.value)
#        if result is None:
#            raise Exception(
#                f"The vcfheader does not contain annotation key {annotkey}."
#            )
#
#        return result
#
#    @classmethod
#    def from_dict(cls, data):
#        return cls(data)
#    ########################
#
#    def encode(self):
#        return _encode(dict(self))
#
#    def load_infostring(self, infostring, overwrite=True):
#        self.update_dict(self.__class__.decode(infostring), overwrite=overwrite)
#    
#    #    def load_vr_info(self, vr, infokey=None):
#    #        if infokey is None:
#    #            infokey = self.annotkey
#    #
#    #        assert vr.header.info[infokey].number == 1
#    #        assert vr.header.info[infokey].type == "String"
#    #
#    #        if not infoformat.check_NA_info(vr, infokey):
#    #            infostring = infoformat.get_info(vr, infokey)
#    #            self.load_infostring(infostring)
#    #
#    #    load_vr = load_vr_info  # alias
#    #
#    #    def load_vr_format(self, vr, sampleid, key=None):
#    #        if key is None:
#    #            key = self.annotkey
#    #
#    #        assert vr.header.formats[key].number == 1
#    #        assert vr.header.formats[key].type == "String"
#    #
#    #        if not infoformat.check_NA_format(vr, sampleid, key):
#    #            infostring = infoformat.get_format(vr, sampleid, key)
#    #            self.load_infostring(infostring)
#
#    def update_dict(self, other, overwrite=False):
#        if overwrite:
#            for key, val in other.items():
#                self[key] = val
#        else:
#            for key, new_val in other.items():
#                try:
#                    old_val = self[key]
#                except KeyError:
#                    self[key] = new_val
#                else:
#                    if old_val is None:
#                        self[key] = new_val
#                    else:
#                        if new_val != old_val:
#                            self._update_conflict_handler(key, old_val, new_val)
#
#    def _update_conflict_handler(self, key, old_val, new_val):
#        values = {'3prime_overlapping_ncRNA', '3prime_overlapping_ncrna'}
#        if (
#            key == 'biotype' and 
#            old_val in values and 
#            new_val in values
#        ):
#            self[key] = '3prime_overlapping_ncRNA'
#        else:
#            raise Exception(f"Old value and new value are different; key: {key}; old_val: {old_val}; new_val: {new_val}")
#
#
#    load_other = update_dict
#
#
#    def write_info(self, vr, key=None, add_meta=False):
#        if key is None:
#            key = self.annotkey
#        if add_meta:
#            self.__class__.add_meta_info(vr.header)
#
#        infoformat.set_info(vr, key, self.encode())
#
#    # write = write_info  # alias
#
#    def write_format(self, vr, sampleid, key=None, add_meta=False):
#        if key is None:
#            key = self.annotkey
#        if add_meta:
#            self.__class__.add_meta_format(vr.header)
#
#        infoformat.set_format(vr, sampleid, key, self.encode())




#class AnnotItemDict(dict):
#    def __init__(self, vr=None, infokey=None, septype=DEFAULT_SEPTYPE):
#        # super().__init__()
#        self.septype = septype
#        if (vr is not None) and (infokey is not None):
#            self.load_vr(vr, infokey)
#
#    def __setitem__(self, key, val):
#        assert isinstance(val, AnnotItem)
#        super().__setitem__(key, val)
#
#    def __repr__(self):
#        return common.cpformat(self.get_showdict())
#
#    def get_showdict(self):
#        showdict = {key: dict(val) for (key, val) in self.items()}
#        return showdict
#
#    def get_annotitem_id(self, annotitem):
#        pass
#
#    def load_vr(self, vr, infokey):
#        assert vr.header.info[infokey].number == "."
#        assert vr.header.info[infokey].type == "String"
#
#        if not infoformat.check_NA_info(vr, infokey):
#            infoval = infoformat.get_info(vr, infokey, collapse_tuple=False)
#            for infostring in infoval:
#                annotitem = AnnotItem(septype=self.septype)
#                annotitem.load_infostring(infostring)
#                ID = self.get_annotitem_id(annotitem)
#                self[ID] = annotitem
#
#    def load_other(self, other, overwrite=True, create_new=False):
#        for ID, annotitem in other.items():
#            if ID in self.keys():
#                self[ID].load_other(annotitem, overwrite=overwrite)
#            else:
#                if create_new:
#                    self[ID] = annotitem
#
#    def write_info(self, vr, key, add_meta=False):
#        if add_meta:
#            add_infometas(vr.header)
#
#        infoval = list()
#        for annotitem in self.values():
#            infoval.append(annotitem.encode())
#
#        vr.info[key] = infoval

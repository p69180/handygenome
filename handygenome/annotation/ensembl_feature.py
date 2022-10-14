import operator
import functools
import itertools
import time

import pysam

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))

infoformat = importlib.import_module(
    ".".join([top_package_name, "variant", "infoformat"])
)
varianthandler = importlib.import_module(
    ".".join([top_package_name, "variant", "varianthandler"])
)

annotitem = importlib.import_module(
    ".".join([top_package_name, "annotation", "annotitem"])
)
ensembl_parser = importlib.import_module(
    ".".join([top_package_name, "annotation", "ensembl_parser"])
)
ensembl_rest = importlib.import_module(
    ".".join([top_package_name, "annotation", "ensembl_rest"])
)
veplib = importlib.import_module(".".join([top_package_name, "annotation", "veplib"]))
customfile = importlib.import_module(
    ".".join([top_package_name, "annotation", "customfile"])
)
rnalib = importlib.import_module(".".join([top_package_name, "annotation", "rnalib"]))
libvcfspec = importlib.import_module('.'.join([top_package_name, 'variant', 'vcfspec']))


META_ID_PREFIX = {
    "transcript": "transcript",
    "regulatory": "regulatory",
    "motif": "motif",
    "repeat": "repeat",
}
META_ID_SUFFIX = {"plain": "", "bnd1": "_bnd1", "bnd2": "_bnd2"}
META_NUMS = {"transcript": ".", "regulatory": ".", "motif": ".", "repeat": "."}
META_DESCRIPTIONS = {
    "transcript": "Transcript feature annotations encoded as a string",
    "regulatory": "Regulatory element feature annotations encoded as a string",
    "motif": "Transcription factor binding motif feature annotations encoded as a string",
    "repeat": "Repeat feature annotations encoded as a string",
}
META_DESCRIPTIONS_SUFFIXES = {
    "plain": "for non-breakend variants (e.g. snv, indel, cnv)",
    "bnd1": "for breakend 1",
    "bnd2": "for breakend 2",
}


###############################################


#REGULATORY_FEATURE_TYPES = {
#    "regulatory": {
#        "Promoter": "promoter",
#        "Promoter Flanking Region": "promoter_flank",
#        "CTCF Binding Site": "CTCFBS",
#        "TF binding site": "TFBS",
#        "Enhancer": "enhancer",
#        "Open chromatin": "open_chromatin",
#    },
#    "overlap": {
#        "Predicted promoter": "promoter",
#        "Predicted promoter flanking region": "promoter_flank",
#        "CTCF binding site": "CTCFBS",
#        "Transcription factor binding site": "TFBS",
#        "Predicted enhancer region": "enhancer",
#        "Open chromatin region": "open_chromatin",
#    },
#    "vep": {
#        "promoter": "promoter",
#        "promoter_flanking_region": "promoter_flank",
#        "CTCF_binding_site": "CTCFBS",
#        "TF_binding_site": "TFBS",
#        "enhancer": "enhancer",
#        "open_chromatin_region": "open_chromatin",
#    },
#}

BIOTYPE_REGULATORY_CLASSES = {
    "promoter": (
        "Promoter",
        "Predicted promoter",
        "promoter",
    ),
    "promoter_flank": (
        "Promoter Flanking Region",
        "Predicted promoter flanking region",
        "promoter_flanking_region",
    ),
    "CTCFBS": (
        "CTCF Binding Site",
        "CTCF binding site",
        "CTCF_binding_site",
    ),
    "TFBS": (
        "TF binding site",
        "Transcription factor binding site",
        "TF_binding_site",
    ),
    "enhancer": (
        "Enhancer",
        "Predicted enhancer region",
        "enhancer",
    ),
    "open_chromatin": (
        "Open chromatin",
        "Open chromatin region",
        "open_chromatin_region",
    ),
}
BIOTYPE_REGULATORY_VALUES = set(
    itertools.chain.from_iterable(BIOTYPE_REGULATORY_CLASSES.values())
)


# reference: https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html
CONSEQUENCE_TRANSCRIPT_CLASSES = {
    "splice_region_involved": (
        "splice_acceptor_variant",
        "splice_donor_variant",
        "splice_donor_region_variant",  # added on 220422
        "splice_region_variant",
        "splice_polypyrimidine_tract_variant",  # added on 220329
        "splice_donor_5th_base_variant",  # added on 220401
    ),
    "splice_acceptor_involved": ("splice_acceptor_variant",),
    "splice_donor_involved": (
        "splice_donor_variant",
        "splice_donor_region_variant",
        "splice_donor_5th_base_variant",  # added on 220401
    ),
    "5pUTR_involved": ("5_prime_UTR_variant",),
    "3pUTR_involved": ("3_prime_UTR_variant",),
    "protein_altering": (
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "start_lost",
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "protein_altering_variant",
    ),
    "not_protein_altering": (
        "splice_acceptor_variant",
        "splice_donor_variant",
        "splice_region_variant",
        "start_retained_variant",
        "stop_retained_variant",
        "synonymous_variant",
        "mature_miRNA_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "non_coding_transcript_exon_variant",
        "intron_variant",
        "NMD_transcript_variant",
        "non_coding_transcript_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "TF_binding_site_variant",
        "regulatory_region_variant",
        "intergenic_variant",
    ),
    "synonymous": ("synonymous_variant",),
    "missense": ("missense_variant",),
    "frameshift": ("frameshift_variant",),
    "inframe": ("inframe_deletion",),
    "stopgain": ("stop_gained",),
    "stoplost": ("stop_lost",),
    "startlost": ("start_lost",),
    "unclassified": (
        "incomplete_terminal_codon_variant",
        "coding_sequence_variant",
    ),
    "SV_consequence": (
        "transcript_ablation",
        "transcript_amplification",
        "TFBS_ablation",
        "TFBS_amplification",
        "regulatory_region_ablation",
        "regulatory_region_amplification",
        "feature_elongation",
        "feature_truncation",
    ),
}


CONSEQUENCE_TRANSCRIPT_VALUES = set(
    itertools.chain.from_iterable(CONSEQUENCE_TRANSCRIPT_CLASSES.values())
)


BIOTYPE_TRANSCRIPT_CLASSES = {
    # coding
    "coding": (
        "IG_C_gene",
        "IG_D_gene",
        "IG_J_gene",
        "IG_V_gene",
        "ccds_gene",
        "mRNA",
        "protein_coding",
    ),
    "IG": (
        "IG_C_gene",
        "IG_D_gene",
        "IG_J_gene",
        "IG_V_gene",
    ),
    "TR": (
        "TR_C_gene",
        "TR_D_gene",
        "TR_J_gene",
        "TR_V_gene",
    ),
    # noncoding
    "noncoding": (
        "3prime_overlapping_ncRNA",  # updated on 220421; ENST00000606379
        #"3prime_overlapping_ncrna",
        "Mt_rRNA",
        "Mt_tRNA",
        "RNase_MRP_RNA",
        "RNase_P_RNA",
        "Y_RNA",
        "antisense_RNA",
        "antisense",  # updated on 220316; ENST00000442411
        "lncRNA",
        "lincRNA",
        "miRNA",
        "misc_RNA",
        "non_stop_decay",
        "nonsense_mediated_decay",
        "processed_transcript",
        "rRNA",
        "retained_intron",
        "ribozyme",
        "sRNA",
        "scRNA",
        "scaRNA",
        "snRNA",
        "snoRNA",
        "tRNA",
        "telomerase_RNA",
        "vault_RNA",
        "sense_overlapping",
        "sense_intronic",  # updated on 220421; ENST00000428191
    ),
    "rRNA": (
        "Mt_rRNA",
        "rRNA",
    ),
    "tRNA": (
        "Mt_tRNA",
        "tRNA",
    ),
    "miRNA": ("miRNA",),
    "lncRNA": (
        "lncRNA",
        "lincRNA",
    ),
    "NMD": ("nonsense_mediated_decay",),
    "NSD": ("non_stop_decay",),
    "antisense": (
        "antisense_RNA",
        "antisense",
    ),
    # pseudogenes
    "pseudogene": (
        "IG_C_pseudogene",
        "IG_J_pseudogene",
        "IG_V_pseudogene",
        "IG_pseudogene",
        "TR_J_pseudogene",
        "TR_V_pseudogene",
        "ncRNA_pseudogene",
        "polymorphic_pseudogene",
        "processed_pseudogene",
        "pseudogene",
        "rRNA_pseudogene",
        "transcribed_processed_pseudogene",
        "transcribed_pseudogene",
        "transcribed_unitary_pseudogene",
        "transcribed_unprocessed_pseudogene",
        "translated_processed_pseudogene",
        "translated_unprocessed_pseudogene",
        "unitary_pseudogene",
        "unprocessed_pseudogene",
    ),
    "IG_pseudogene": (
        "IG_C_pseudogene",
        "IG_J_pseudogene",
        "IG_V_pseudogene",
        "IG_pseudogene",
    ),
    "TR_pseudogene": (
        "TR_J_pseudogene",
        "TR_V_pseudogene",
    ),
    "processed_pseudogene": (
        "processed_pseudogene",
        "transcribed_processed_pseudogene",
        "transcribed_unprocessed_pseudogene",
        "translated_processed_pseudogene",
        "translated_unprocessed_pseudogene",
        "unprocessed_pseudogene",
    ),
    "unprocessed_pseudogene": (
        "transcribed_unprocessed_pseudogene",
        "translated_unprocessed_pseudogene",
        "unprocessed_pseudogene",
    ),
    "unitary_pseudogene": (
        "transcribed_unitary_pseudogene",
        "unitary_pseudogene",
    ),
    "translated_pseudogene": (
        "translated_processed_pseudogene",
        "translated_unprocessed_pseudogene",
    ),
    "transcribed_pseudogene": (
        "transcribed_processed_pseudogene",
        "transcribed_pseudogene",
        "transcribed_unitary_pseudogene",
        "transcribed_unprocessed_pseudogene",
    ),
    # unknown types
    "unknown_type": (
        "LRG_gene",
        "TEC",
        "aligned_transcript",
        "cdna_update",
        "guide_RNA",
        "other",
    ),
}

BIOTYPE_TRANSCRIPT_VALUES = set(
    itertools.chain.from_iterable(BIOTYPE_TRANSCRIPT_CLASSES.values())
)


###############################################


class EnsemblFeature(annotitem.AnnotItemVariantInfoSingle):
    @property
    def is_overlapping(self):
        return self["distance"] is None


class EnsemblFeatureSetBase(annotitem.AnnotItemVariantInfoSingle):
    @classmethod
    def from_annotstring_base(cls, annotstring, feature_class):
        result = cls()
        tmp = cls.decode(annotstring)
        for key, val in tmp.items():
            result[key] = feature_class.from_dict(val)
        return result

    @property
    def overlapping(self):
        result = self.__class__()
        for key, val in self.items():
            if val.is_overlapping:
                result[key] = val
        return result

    def get_self_show(self):
        return {key: val.get_self_show() for key, val in self.items()}

    def add_feature(self, feature):
        self[feature["id"]] = feature

    def update_other(self, other, overwrite=True, create_new=True):
        """Args:
            other: another EnsemblFeatureSetBase object
        """
        for ID, feature in other.items():
            if ID in self.keys():
                self[ID].update_dict(feature, overwrite=overwrite)
            else:
                if create_new:
                    self[ID] = feature

    def get_distance_setter_snv(self, vcfspec):
        pos0 = vcfspec.pos0
        def distance_setter(feature_range0):
            if pos0 < feature_range0.start:
                distance = feature_range0.start - pos0
            elif pos0 >= feature_range0.stop:
                distance = feature_range0.stop - 1 - pos0
            else:
                distance = None
            return distance
        return distance_setter

    def get_distance_setter_ins(self, vcfspec):
        pos0 = vcfspec.pos0
        def distance_setter(feature_range0):
            if pos0 <= feature_range0.start - 1:
                distance = feature_range0.start - pos0 - 1
            elif pos0 >= feature_range0.stop - 1:
                distance = feature_range0.stop - 1 - pos0
            else:
                distance = None
            return distance
        return distance_setter

    def get_distance_setter_interval(self, var_range0):
        def distance_setter(feature_range0):
            if var_range0.stop <= feature_range0.start:
                distance = feature_range0.start - var_range0.stop + 1
            elif var_range0.start >= feature_range0.stop:
                distance = feature_range0.stop - 1 - var_range0.start
            else:
                distance = None
            return distance
        return distance_setter

    def get_distance_setter_del(self, vcfspec):
        var_range0 = range(vcfspec.pos, vcfspec.pos + len(vcfspec.ref) - 1)
        return self.get_distance_setter_interval(var_range0)

    def get_distance_setter_bnd(self, pos):
        vcfspec = libvcfspec.Vcfspec(chrom=None, pos=pos, ref=None, alts=None)
        return self.get_distance_setter_snv(vcfspec)

    def set_missing_distances_base(self, distance_setter):
        """Sets distance for features not derived from VEP (e.g. repeat
        elements). Simulates the distance given by VEP.
        """
        for feature in self.values():
            if "distance" not in feature:
                # features derived from "ensembl rest overlap"
                #   do not have "distance" value
                assert "start0" in feature and "end0" in feature, (
                    f'"start0" or "end0" missing for this feature:\n{feature}'
                )
                feature_range0 = range(feature["start0"], feature["end0"])
                feature["distance"] = distance_setter(feature_range0)


class EnsemblFeatureSetPlain(EnsemblFeatureSetBase):
    def get_fetch_interval(self, vcfspec, distance, chromdict):
        return common.Interval(
            chrom=vcfspec.chrom,
            start1=max(1, vcfspec.pos - distance),
            end1=min(chromdict[vcfspec.chrom], vcfspec.pos + distance),
        )

    def set_missing_distances(self, vcfspec, alt_idx=0):
        mttype = vcfspec.get_mutation_type(alt_idx)
        if mttype == "snv":
            distance_setter = self.get_distance_setter_snv(vcfspec)
        elif mttype == "ins":
            distance_setter = self.get_distance_setter_ins(vcfspec)
        elif mttype in ("del", "delins", "mnv", "cpgmet"):
            distance_setter = self.get_distance_setter_del(vcfspec)
        else:
            raise Exception(f"Unexpected mutation type: {mttype}")

        self.set_missing_distances_base(distance_setter)


class EnsemblFeatureSetBreakend(EnsemblFeatureSetBase):
    pass


########


class Transcript(EnsemblFeature):
    def __getitem__(self, key):
        if key == 'consequence_flags':
            return self.consequence_flags
        elif key == 'subtype_flags':
            return self.subtype_flags
        elif key == 'exon_ranges':
            return self.exon_ranges
        elif key == 'intron_ranges':
            return self.intron_ranges
        else:
            return super().__getitem__(key)

    def get_self_show(self):
        self_show = dict(self)
        for key in (
            'consequence_flags', 
            'subtype_flags',
            'exon_ranges',
            'intron_ranges',
        ):
            try:
                self_show[key] = self[key]
            except KeyError:
                pass

        return self_show

    @functools.cached_property
    def consequence_flags(self):
        consequences = set(self["consequences"])
        if not consequences.issubset(CONSEQUENCE_TRANSCRIPT_VALUES):
            unexpected_consequences = consequences.difference(
                CONSEQUENCE_TRANSCRIPT_VALUES
            )
            raise Exception(
                f'Unexpected consequence value; consequence: {unexpected_consequences}; feature_id: {self["id"]}'
            )
        result = dict()
        for key, val in CONSEQUENCE_TRANSCRIPT_CLASSES.items():
            result[key] = bool(consequences.intersection(val))
        return result

    @functools.cached_property
    def subtype_flags(self):
        if self["biotype"] not in BIOTYPE_TRANSCRIPT_VALUES:
            raise Exception(
                f'Unexpected transcript biotype value; biotype: {self["biotype"]}; feature_id: {self["id"]}'
            )
        result = dict()
        for key, val in BIOTYPE_TRANSCRIPT_CLASSES.items():
            result[key] = self["biotype"] in val
        return result

    @functools.cached_property
    def exon_ranges(self):
        if 'exon_borders' in self.keys():
            if self['exon_borders'] is None:
                return None
            else:
                return [range(*x) for x in self['exon_borders']]
        else:
            return None

    @functools.cached_property
    def intron_ranges(self):
        exon_ranges = self.exon_ranges
        if exon_ranges is None:
            result = None
        else:
            if len(exon_ranges) == 1:
                result = None
            else:
                assert self['is_forward'] is not None
                result = list()
                if self['is_forward']:
                    for idx in range(len(exon_ranges) - 1):
                        result.append(
                            range(exon_ranges[idx].stop, exon_ranges[idx + 1].start)
                        )
                else:
                    for idx in range(len(exon_ranges) - 1):
                        result.append(
                            range(exon_ranges[idx + 1].stop, exon_ranges[idx].start)
                        )
        return result


class TranscriptSet(EnsemblFeatureSetPlain):
    # constructor
    @classmethod
    def from_annotstring(cls, annotstring):
        return cls.from_annotstring_base(annotstring, Transcript)

    @property
    def canonical(self):
        result = self.__class__()
        for key, val in self.items():
            if val['is_canonical']:
                result[key] = val
        return result

    @property
    def canon_ovlp(self):
        result = self.__class__()
        for key, val in self.items():
            if (val['is_canonical'] and val.is_overlapping):
                result[key] = val
        return result

    def get_gene_names(self, canonical=True, overlap=True):
        if canonical:
            if overlap:
                features = self.canon_ovlp.values()
            else:
                features = self.canonical.values()
        else:
            if overlap:
                features = self.overlapping.values()
            else:
                features = self.values()

        return set(x['gene_name'] for x in features)

    def update_ensembl_gff(self, vcfspec, chromdict, distance, tabixfile_geneset, refver=None, fill_missing_canonical=True):
        assert not ((refver is None) and fill_missing_canonical), (
            f'If "fill_missing_canonical" is True, "refver" must be given.'
        )

        fetch_interval = self.get_fetch_interval(vcfspec, distance, chromdict)
        transcript_set = customfile.fetch_transcript(
            fetch_interval.chrom,
            fetch_interval.start0,
            fetch_interval.end0,
            tabixfile_geneset,
        )

        if fill_missing_canonical:
            for ID, transcript in transcript_set.items():
                if "is_canonical" not in transcript.keys():
                    raw_result = ensembl_rest.lookup_id(
                        ID=transcript["id"], refver=refver, expand=False
                    )
                    time.sleep(0.1)
                        # You have exceeded the limit of 15 requests per second; please reduce your concurrent connections
                    transcript["is_canonical"] = (raw_result["is_canonical"] == 1)

        self.update_other(transcript_set, overwrite=False, create_new=False)

    def update_postvep(self, vcfspec, chromdict, distance, tabixfile_geneset):
        self.update_ensembl_gff(vcfspec, chromdict, distance, tabixfile_geneset, refver=None, fill_missing_canonical=False)


class TranscriptSetALTlist(annotitem.AnnotItemVariantInfoALTlist):
    meta = {
        "ID": "transcripts",
        "Number": "A",
        "Type": "String",
        "Description": "Transcript annotations encoded as a string, one for each ALT allele",
    }

    @classmethod
    def annotstring_parser(cls, annotstring):
        return TranscriptSet.from_annotstring(annotstring)

    @classmethod
    def from_vr(cls, vr):
        return cls.from_vr_base(vr)


class Regulatory(EnsemblFeature):
    def __getitem__(self, key):
        if key == 'subtype_flags':
            return self.subtype_flags
        else:
            return super().__getitem__(key)

    def get_self_show(self):
        self_show = dict(self)
        for key in (
            'subtype_flags',
        ):
            try:
                self_show[key] = self[key]
            except KeyError:
                pass

        return self_show

    @functools.cached_property
    def subtype_flags(self):
        if self["biotype"] not in BIOTYPE_REGULATORY_VALUES:
            raise Exception(
                f'Unexpected regulatory biotype value;\n{dict(self)}'
            )
        result = dict()
        for key, val in BIOTYPE_REGULATORY_CLASSES.items():
            result[key] = self["biotype"] in val
        return result


class RegulatorySet(EnsemblFeatureSetPlain):
    meta = {
        "ID": "regulatory_elements",
        "Number": "1",
        "Type": "String",
        "Description": "Regulatory element annotations encoded as a string.",
    }

    @classmethod
    def annotstring_parser(cls, annotstring):
        return cls.from_annotstring_base(annotstring, Regulatory)

    # constructor
    @classmethod
    def from_vr(cls, vr):
        return cls.from_vr_base(vr)

    def update_ensembl_gff(self, vcfspec, chromdict, distance, tabixfile_regulatory):
        fetch_interval = self.get_fetch_interval(vcfspec, distance, chromdict)
        regulatory_set = customfile.fetch_regulatory(
            fetch_interval.chrom,
            fetch_interval.start0,
            fetch_interval.end0,
            tabixfile_regulatory,
        )
        self.update_other(regulatory_set, overwrite=False, create_new=True)


class Motif(EnsemblFeature):
    def get_self_show(self):
        return dict(self)


class MotifSet(EnsemblFeatureSetPlain):
    meta = {
        "ID": "TF_binding_motifs",
        "Number": "1",
        "Type": "String",
        "Description": "Transcription factor binding motif annotations encoded as a string.",
    }

    @classmethod
    def annotstring_parser(cls, annotstring):
        return cls.from_annotstring_base(annotstring, Motif)

    # constructor
    @classmethod
    def from_vr(cls, vr):
        return cls.from_vr_base(vr)


class Repeat(EnsemblFeature):
    def __getitem__(self, key):
        if key == 'id':
            return self.id
        else:
            return super().__getitem__(key)

    def get_self_show(self):
        self_show = dict(self)
        for key in (
            'id',
        ):
            try:
                self_show[key] = self[key]
            except KeyError:
                pass

        return self_show

    @functools.cached_property
    def id(self):
        return '_'.join(
            [
                self['name'],
                self['chrom'], 
                str(self['start1']), 
                str(self['end1']), 
            ]
        )

class RepeatSet(EnsemblFeatureSetPlain):
    meta = {
        "ID": "repeat_elements",
        "Number": "1",
        "Type": "String",
        "Description": "Repeat sequence annotations encoded as a string.",
    }

    @classmethod
    def annotstring_parser(cls, annotstring):
        return cls.from_annotstring_base(annotstring, Repeat)

    # constructor
    @classmethod
    def from_vr(cls, vr):
        return cls.from_vr_base(vr)

    def update_repeatmasker_bed(self, vcfspec, chromdict, distance, tabixfile_repeats):
        fetch_interval = self.get_fetch_interval(vcfspec, distance, chromdict)
        repeat_set = customfile.fetch_repeat(
            fetch_interval.chrom,
            fetch_interval.start0,
            fetch_interval.end0,
            tabixfile_repeats,
        )
        self.update_other(repeat_set, overwrite=False, create_new=True)


###################


#class EnsemblFeatureDict(annotitem.AnnotItemDict):
#    @common.get_deco_arg_choices({"annotation_type": ANNOTATION_TYPES})
#    def __init__(self, annotation_type="plain", **kwargs):
#        super().__init__(**kwargs)
#        self.annotation_type = annotation_type
#
#    def __setitem__(self, key, val):
#        # assert isinstance(val, EnsemblFeature)
#        super().__setitem__(key, val)
#
#    def get_showdict(self):
#        showdict = dict()
#        for key, val in sorted(self.items(), key=(lambda x: x[1]["start0"])):
#            showdict[key] = dict(val)
#
#        return showdict
#
#    def get_list(self):
#        return sorted(self.values(), key=operator.attrgetter("start0"))
#
#    def get_annotitem_id(self, annotitem):
#        return annotitem["id"]
#
#    def add_meta(self, vcfheader):
#        vcfheader.add_meta(key="INFO", items=self.get_meta().items())
#
#
#class TranscriptDict(EnsemblFeatureDict):
#    def get_meta(self):
#        return INFOMETAS[self.annotation_type]["transcript"]
#
#    def get_canon(self):
#        result = self.__class__()
#        for ID, transcript in self.items():
#            if transcript.check_canonical():
#                result[ID] = transcript
#
#        return result
#
#    def get_ovlp(self):
#        result = self.__class__()
#        for ID, transcript in self.items():
#            if transcript.check_overlapping():
#                result[ID] = transcript
#
#        return result
#
#    def get_canon_ovlp(self):
#        result = self.__class__()
#        for ID, transcript in self.items():
#            if transcript.check_canonical() and transcript.check_overlapping():
#                result[ID] = transcript
#
#        return result
#
#
#class RegulatoryDict(EnsemblFeatureDict):
#    def get_meta(self):
#        return INFOMETAS[self.annotation_type]["regulatory"]
#
#
#class MotifDict(EnsemblFeatureDict):
#    def get_meta(self):
#        return INFOMETAS[self.annotation_type]["motif"]
#
#
#class RepeatDict(EnsemblFeatureDict):
#    def get_meta(self):
#        return INFOMETAS[self.annotation_type]["repeat"]

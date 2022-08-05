import re
import pprint
import sys
import collections
from operator import attrgetter
import copy
import functools

import pysam

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))

infoformat = importlib.import_module(
    ".".join([top_package_name, "variantplus", "infoformat"])
)
varianthandler = importlib.import_module(
    ".".join([top_package_name, "variantplus", "varianthandler"])
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


REGULATORY_FEATURE_TYPES = {
    "regulatory": {
        "Promoter": "promoter",
        "Promoter Flanking Region": "promoter_flank",
        "CTCF Binding Site": "CTCFBS",
        "TF binding site": "TFBS",
        "Enhancer": "enhancer",
        "Open chromatin": "open_chromatin",
    },
    "overlap": {
        "Predicted promoter": "promoter",
        "Predicted promoter flanking region": "promoter_flank",
        "CTCF binding site": "CTCFBS",
        "Transcription factor binding site": "TFBS",
        "Predicted enhancer region": "enhancer",
        "Open chromatin region": "open_chromatin",
    },
    "vep": {
        "promoter": "promoter",
        "promoter_flanking_region": "promoter_flank",
        "CTCF_binding_site": "CTCFBS",
        "TF_binding_site": "TFBS",
        "enhancer": "enhancer",
        "open_chromatin_region": "open_chromatin",
    },
}

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
        "3prime_overlapping_ncrna",
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


class EnsemblFeature(annotitem.AnnotItem):
    feature_type = None

    @property
    def is_overlapping(self):
        return self["distance"] is None


class Transcript(EnsemblFeature):
    feature_type = 'transcript'

    @property
    def is_canonical(self):
        return self["is_canonical"]

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

    @functools.cached_property
    def exon_ranges(self):
        pass

    @functools.cached_property
    def intron_ranges(self):
        exon_ranges = self.get_exon_ranges()
        if exon_ranges is not None:
            pass

        if len(exon_ranges) == 1:
            intron_ranges = None
        else:
            intron_ranges = dict()
            ascending = sorted(exon_ranges.items(), key=lambda x: x[1][0])
            for idx in range(len(ascending) - 1):
                intron_number = (
                    ascending[idx][0]
                    if annotitem["is_forward"]
                    else ascending[idx + 1][0]
                )
                range_item = range(ascending[idx][1][-1] + 1, ascending[idx + 1][1][0])
                intron_ranges[intron_number] = range_item
            intron_ranges = dict(sorted(intron_ranges.items(), key=lambda x: x[0]))

        return intron_ranges


class TranscriptInfo(annotitem.AnnotItem):
    pass


class TranscriptInfoList(annotitem.AnnotItemList):
    pass


class Regulatory(EnsemblFeature):
    feature_type = 'regulatory'

    @functools.cached_property
    def subtype_flags(self):
        if self["biotype"] not in BIOTYPE_REGULATORY_VALUES:
            raise Exception(
                f'Unexpected regulatory biotype value; biotype: {self["biotype"]}; feature_id: {self["id"]}'
            )
        result = dict()
        for key, val in BIOTYPE_REGULATORY_CLASSES.items():
            result[key] = self["biotype"] in val


class RegulatoryInfo(annotitem.AnnotItem):
    pass


class Motif(EnsemblFeature):
    feature_type = 'motif'


class MotifInfo(annotitem.AnnotItem):
    pass


class Repeat(EnsemblFeature):
    feature_type = 'repeat'


class RepeatInfo(annotitem.AnnotItem):
    pass


###################


class EnsemblFeatureDict(annotitem.AnnotItemDict):
    @common.get_deco_arg_choices({"annotation_type": ANNOTATION_TYPES})
    def __init__(self, annotation_type="plain", **kwargs):
        super().__init__(**kwargs)
        self.annotation_type = annotation_type

    def __setitem__(self, key, val):
        # assert isinstance(val, EnsemblFeature)
        super().__setitem__(key, val)

    def get_showdict(self):
        showdict = dict()
        for key, val in sorted(self.items(), key=(lambda x: x[1]["start0"])):
            showdict[key] = dict(val)

        return showdict

    def get_list(self):
        return sorted(self.values(), key=attrgetter("start0"))

    def get_annotitem_id(self, annotitem):
        return annotitem["id"]

    def add_meta(self, vcfheader):
        vcfheader.add_meta(key="INFO", items=self.get_meta().items())


class TranscriptDict(EnsemblFeatureDict):
    def get_meta(self):
        return INFOMETAS[self.annotation_type]["transcript"]

    def get_canon(self):
        result = self.__class__()
        for ID, transcript in self.items():
            if transcript.check_canonical():
                result[ID] = transcript

        return result

    def get_ovlp(self):
        result = self.__class__()
        for ID, transcript in self.items():
            if transcript.check_overlapping():
                result[ID] = transcript

        return result

    def get_canon_ovlp(self):
        result = self.__class__()
        for ID, transcript in self.items():
            if transcript.check_canonical() and transcript.check_overlapping():
                result[ID] = transcript

        return result


class RegulatoryDict(EnsemblFeatureDict):
    def get_meta(self):
        return INFOMETAS[self.annotation_type]["regulatory"]


class MotifDict(EnsemblFeatureDict):
    def get_meta(self):
        return INFOMETAS[self.annotation_type]["motif"]


class RepeatDict(EnsemblFeatureDict):
    def get_meta(self):
        return INFOMETAS[self.annotation_type]["repeat"]

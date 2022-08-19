import re

import pysam

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
infoformat = importlib.import_module(
    ".".join([top_package_name, "variantplus", "infoformat"])
)
annotitem = importlib.import_module(
    ".".join([top_package_name, "annotation", "annotitem"])
)
customfile = importlib.import_module(
    ".".join([top_package_name, "annotation", "customfile"])
)


class PopfreqInfo(annotitem.AnnotItemVariantInfoSingle):
    # constructors
    def __init__(self, metadata=None, **kwargs):
        super().__init__(**kwargs)
        self.metadata = metadata

    @classmethod
    def init_blank(cls):
        result = cls()
        result['id'] = None
        result['dbSNPBuild'] = None
        result['common'] = None
        result['freqs'] = dict()
        return result
    ##############

    def get_self_show(self):
        result = dict()
        for key in (
            'id',
            'dbSNPBuild',
            'common',
        ):
            if key in self:
                result[key] = self[key]
            else:
                result[key] = None
        result['freqs'] = self.get_freqs_show()

        return result

    def get_freq(self, popname):
        try:
            return self['freqs'][popname]
        except KeyError:
            if popname in self.metadata['popnames']:
                return 0
            else:
                raise Exception(f'Unsupported population name')

    def get_freqs_show(self):
        return {popname: self.get_freq(popname)
                for popname in self.metadata['popnames']}


class PopfreqInfoALTlist(annotitem.AnnotItemVariantInfoALTlist):
    meta = {
        "ID": "popfreq",
        "Number": "A",
        "Type": "String",
        "Description": "Population frequencies encoded as a string, one for each ALT allele",
    }

    @classmethod
    def annotstring_parser(cls, annotstring, metadata):
        result = PopfreqInfo.from_annotstring(annotstring)
        result.metadata = metadata
        return result

    @classmethod
    def from_vr(cls, vr, metadata=None):
        if metadata is None:
            metadata = PopfreqMetadata.from_vcfheader(vr.header)

        annotstring_parser_kwargs = {'metadata': metadata}
        result = cls.from_vr_base(vr, annotstring_parser_kwargs)

        return result

    @classmethod
    def from_vcfspec(cls, vcfspec, dbsnp_vcf, fasta, metadata=None, donot_init_metadata=False):
        if metadata is None:
            if not donot_init_metadata:
                metadata = PopfreqMetadata.from_vcfheader(dbsnp_vcf.header)

        dbsnp_vr_list = customfile.fetch_relevant_vr_multialt(
            vcfspec, dbsnp_vcf, fasta=fasta, search_equivs=True, allow_multiple=False
        )
        result = cls()
        for vr in dbsnp_vr_list:
            if vr is None:
                result.append(PopfreqInfo.init_blank())
            else:
                result.extend(cls.from_vr(vr, metadata=metadata))

        return result


class PopfreqMetadata(annotitem.AnnotItemHeader):
    meta = {"ID": "popfreq_metadata"}


def update_vcfheader(vcfheader, dbsnp_vcf):
    PopfreqInfoALTlist.add_meta(vcfheader)
    popfreqmeta = PopfreqMetadata.from_vcfheader(dbsnp_vcf.header)
    popfreqmeta.write(vcfheader)


def extract_population_names(dbsnp_vcf_header):
    pop_names = list()
    for key in dbsnp_vcf_header.info:
        if key.startswith("AF_"):
            pop_names.append(re.sub("^AF_", "", key))

    return pop_names


def fetch_dbsnp_vr(vcfspec, dbsnp_vcf, fasta, search_equivs=True):
    """Result may be None"""
    dbsnp_vr = customfile.fetch_relevant_vr(
        vcfspec,
        dbsnp_vcf,
        fasta=fasta,
        search_equivs=search_equivs,
        allow_multiple=False,
    )

    return dbsnp_vr

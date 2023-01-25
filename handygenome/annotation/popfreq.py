import re

import pysam
import numpy as np

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
infoformat = importlib.import_module(
    ".".join([top_package_name, "variant", "infoformat"])
)
annotitem = importlib.import_module(
    ".".join([top_package_name, "annotation", "annotitem"])
)
customfile = importlib.import_module(
    ".".join([top_package_name, "annotation", "customfile"])
)
annotation_data = importlib.import_module(
    ".".join([top_package_name, "annotation", "data"])
)


#from handygenome.annotation.annotitem import AnnotItemInfoALTlist


DBSNP_VCFS = annotation_data.VCFS_DBSNP


class PopfreqInfo(annotitem.AnnotItemInfoSingle):
    # constructors
    #def __init__(self, metadata=None, **kwargs):
    #    super().__init__(**kwargs)
    #    self.metadata = metadata

    @classmethod
    def init_blank(cls, metadata=None):
        result = cls.init_nonmissing()
        result['id'] = None
        result['dbSNPBuild'] = None
        result['common'] = None
        result['freqs'] = dict()
        result.metadata = metadata
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
                #return 0
                return np.nan
            else:
                raise Exception(f'Unsupported population name')

    def get_freqs_show(self):
        return {
            popname: self.get_freq(popname)
            for popname in self.metadata['popnames']
        }


class PopfreqMetadata(annotitem.AnnotItemHeader):
    meta = {"ID": "popfreq_metadata"}

    @classmethod
    def from_vcfheader(cls, vcfheader):
        return cls.from_vcfheader_base(vcfheader)

    def write(self, vcfheader):
        self.write_base(vcfheader)


class PopfreqInfoALTlist(annotitem.AnnotItemInfoALTlist):
#class PopfreqInfoALTlist(AnnotItemInfoALTlist):
    meta = {
        "ID": "popfreq",
        "Number": "A",
        "Type": "String",
        "Description": "Population frequencies encoded as a string, one for each ALT allele",
    }
    unit_class = PopfreqInfo
    metadata_class = PopfreqMetadata

    @classmethod
    def from_vr(cls, vr, metadata=None):
        if metadata is None:
            metadata = cls.metadata_class.from_vcfheader(vr.header)

        result = cls.from_vr_base(vr)
        for unit in result:
            unit.metadata = metadata

        return result

    @classmethod
    def from_vcfspec(cls, vcfspec, dbsnp_vcf, metadata=None, donot_init_metadata=False):
        if metadata is None:
            if not donot_init_metadata:
                metadata = cls.metadata_class.from_vcfheader(dbsnp_vcf.header)

        dbsnp_vr_list = customfile.fetch_relevant_vr_multialt(
            vcfspec, dbsnp_vcf, search_equivs=True, raise_with_multihit=True,
        )
        result = cls()
        for vr in dbsnp_vr_list:
            if vr is None:
                result.append(cls.unit_class.init_blank(metadata=metadata))
                #result.append(cls.unit_class.init_missing())
            else:
                result.extend(cls.from_vr(vr, metadata=metadata))

        return result

    def write(self, vr):
        self.write_base(vr)


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


def fetch_dbsnp_vr(vcfspec, dbsnp_vcf, search_equivs=True):
    """Result may be None"""
    return customfile.fetch_relevant_vr(
        vcfspec,
        dbsnp_vcf,
        search_equivs=search_equivs,
        raise_with_multihit=True,
    )


import re

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
infoformat = importlib.import_module('.'.join([top_package_name, 'variantplus', 'infoformat']))
annotitem = importlib.import_module('.'.join([top_package_name, 'annotation', 'annotitem']))
customfile = importlib.import_module('.'.join([top_package_name, 'annotation', 'customfile']))


class PopfreqInfo(annotitem.AnnotItem):
    def __init__(self, metadata=None, **kwargs):
        super().__init__(**kwargs)
        self.metadata = metadata

    def get_self_show(self):
        return dict(sorted(self.items(), key=(lambda x: x[1]), reverse=True))


class PopfreqInfoList(annotitem.AnnotItemList):
    meta = {
            'ID': 'popfreq', 'Number': 'A', 'Type': 'String', 
            'Description': 'Population frequencies encoded as a string, one for each ALT allele',
    }

    @classmethod
    def from_vr(cls, vr, metadata=None):
        if metadata is None:
            metadata = PopfreqMetadata.from_vcfheader(vr.header)

        popinfolist = cls()
        if infoformat.check_NA_info(vr, cls.get_annotkey()):
            transl_num = infoformat.get_translated_number_info(vr, cls.get_annotkey())
            popinfolist.extend([None] * transl_num)
        else:
            for infostring in vr.info[cls.get_annotkey()]:
                if infostring is None:
                    popinfolist.append(None)
                else:
                    popinfo = PopfreqInfo(metadata=metadata)
                    popinfo.load_infostring(infostring)
                    popinfolist.append(popinfo)

        return popinfolist

    @classmethod
    def from_vcfspec(cls, vcfspec, dbsnp_vcf, fasta, metadata=None):
        if metadata is None:
            metadata = PopfreqMetadata.from_vcfheader(vr.header)

        dbsnp_vr_list = customfile.fetch_relevant_vr_multialt(vcfspec, dbsnp_vcf, fasta=fasta, search_equivs=True, allow_multiple=False)
        popinfolist = cls()
        for vr in dbsnp_vr_list:
            if vr is None:
                popinfolist.append(None)
            else:
                popinfolist.extend(PopfreqInfoList.from_vr(vr, metadata=metadata))

        return popinfolist


class PopfreqMetadata(annotitem.AnnotItem):
    meta = {'ID': 'popfreq_metadata'}

    def get_self_show(self):
        return dict(self)


def extract_population_names(dbsnp_vcf_header):
    pop_names = list()
    for key in dbsnp_vcf_header.info:
        if key.startswith('AF_'):
            pop_names.append(re.sub('^AF_', '', key))

    return pop_names


def fetch_dbsnp_vr(vcfspec, dbsnp_vcf, fasta, search_equivs=True):
    """Result may be None"""
    dbsnp_vr = customfile.fetch_relevant_vr(vcfspec, dbsnp_vcf, fasta=fasta,
                                            search_equivs=search_equivs,
                                            allow_multiple=False)
    
    return dbsnp_vr



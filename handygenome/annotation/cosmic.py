import copy

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
infoformat = importlib.import_module('.'.join([top_package_name, 'variant', 'infoformat']))
annotitem = importlib.import_module('.'.join([top_package_name, 'annotation', 'annotitem']))
customfile = importlib.import_module('.'.join([top_package_name, 'annotation', 'customfile']))
annotation_misc = importlib.import_module(".".join([top_package_name, "annotation", "misc"]))

COSMIC_VCFS = annotation_misc.VCFS_COSMIC


class CosmicInfo(annotitem.AnnotItemVariantInfoSingle):
    # constructors
    def __init__(self, metadata=None, **kwargs):
        super().__init__(**kwargs)
        self.metadata = metadata

    @classmethod
    def init_blank(cls):
        result = cls()
        result['id'] = None
        result['occurrence'] = None
        result['occurrence_somatic'] = None
        result['coding_score'] = None
        result['noncoding_score'] = None
        return result
    ##############

    def get_self_show(self):
        self_show = dict()

        for key in (
            'id', 
            'occurrence', 'portion_by_site', 'total_occurrence', 'total_portion', 
            'occurrence_somatic', 'portion_by_site_somatic', 'total_occurrence_somatic', 'total_portion_somatic', 
        ):
            try:
                self_show[key] = self[key]
            except KeyError:
                self_show[key] = None

        for key in ('occurrence', 'occurrence_somatic',
                    'portion_by_site', 'portion_by_site_somatic'):
            try:
                if self_show[key] is not None:
                    self_show[key] = dict(sorted(self_show[key].items(),
                                                 key=(lambda x: x[1]),
                                                 reverse=True))
            except KeyError:
                self_show[key] = None

        return self_show

    def __getitem__(self, key):
        """Real dictionary keys are: "id", "occurrence", "occurrence_somatic", 
        "coding_score", "noncoding_score"
        """
        if key == 'portion_by_site':
            return self.get_portion_by_site(somatic=False)
        elif key == 'portion_by_site_somatic':
            return self.get_portion_by_site(somatic=True)
        elif key == 'total_occurrence':
            return self.get_total_occurrence(somatic=False)
        elif key == 'total_occurrence_somatic':
            return self.get_total_occurrence(somatic=True)
        elif key == 'total_portion':
            return self.get_total_portion(somatic=False)
        elif key == 'total_portion_somatic':
            return self.get_total_portion(somatic=True)
        else:
            return super().__getitem__(key)

    def get_portion_by_site(self, somatic=False):
        key = 'occurrence_somatic' if somatic else 'occurrence'
        if self[key] is None:
            return None
        else:
            portion_by_site = dict()
            for site, count in self[key].items():
                num_sample_site = self.metadata['num_sample_by_site'][site]
                frac = count / num_sample_site
                portion_by_site[site] = frac
            return portion_by_site

    def get_total_occurrence(self, somatic=False):
        key = 'occurrence_somatic' if somatic else 'occurrence'
        if self[key] is None:
            return None
        else:
            return sum(self[key].values())

    def get_total_portion(self, somatic=False):
        total_occurrence = self.get_total_occurrence(somatic=somatic)
        if total_occurrence is None:
            return None
        else:
            num_sample_allsites = sum(self.metadata['num_sample_by_site'].values())
            return total_occurrence / num_sample_allsites


class CosmicInfoALTlist(annotitem.AnnotItemVariantInfoALTlist):
    meta = {'ID': 'cosmic_info', 'Number': 'A', 'Type': 'String', 'Description': 'COSMIC information encoded as a string, one for each ALT allele'}

    @classmethod
    def annotstring_parser(cls, annotstring, metadata):
        result = CosmicInfo.from_annotstring(annotstring)
        result.metadata = metadata
        return result

    @classmethod
    def from_vr(cls, vr, metadata=None):
        if metadata is None:
            metadata = CosmicMetadata.from_vcfheader(vr.header)

        annotstring_parser_kwargs = {'metadata': metadata}
        result = cls.from_vr_base(vr, annotstring_parser_kwargs)

        return result

    @classmethod
    def from_vcfspec(cls, vcfspec, cosmic_vcf, fasta, metadata=None, donot_init_metadata=False):
        if metadata is None:
            if not donot_init_metadata:
                metadata = CosmicMetadata.from_vcfheader(cosmic_vcf.header)

        cosmic_vr_list = customfile.fetch_relevant_vr_multialt(vcfspec, cosmic_vcf, fasta=fasta, search_equivs=True, allow_multiple=False)
        result = cls()
        for vr in cosmic_vr_list:
            if vr is None:
                result.append(CosmicInfo.init_blank())
            else:
                result.extend(CosmicInfoALTlist.from_vr(vr, metadata=metadata))

        return result


class CosmicMetadata(annotitem.AnnotItemHeader):
    meta = {'ID': 'cosmic_metadata'}


def update_vcfheader(vcfheader, cosmic_vcf):
    CosmicInfoALTlist.add_meta(vcfheader)
    cosmicmeta = CosmicMetadata.from_vcfheader(cosmic_vcf.header)
    cosmicmeta.write(vcfheader)


#def annotate_vp(vp, cosmic_vcf, fasta, metadata=None):
#    vp.cosmic = CosmicInfoALTlist.from_vcfspec(vp.vcfspec, cosmic_vcf, fasta, metadata=metadata)



import copy

import pysam

import handygenome.annotation.annotitem as libannotitem
import handygenome.annotation.customfile as customfile


class CosmicInfo(libannotitem.AnnotItemVariantInfoSingle):
    @classmethod
    def init_blank(cls, metadata=None):
        result = cls.init_nonmissing()
        result['id'] = None
        result['occurrence'] = None
        result['occurrence_somatic'] = None
        result['coding_score'] = None
        result['noncoding_score'] = None
        result.metadata = metadata
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


class CosmicMetadata(libannotitem.AnnotItemHeader):
    meta = {'ID': 'cosmic_metadata'}

    @classmethod
    def from_vcfheader(cls, vcfheader):
        return cls.from_vcfheader_base(vcfheader)

    def write(self, vcfheader):
        self.write_base(vcfheader)


class CosmicInfoALTlist(libannotitem.AnnotItemVariantInfoALTlist):
    meta = {
        'ID': 'cosmic_info', 
        'Number': 'A', 
        'Type': 'String', 
        'Description': 'COSMIC information encoded as a string, one for each ALT allele',
    }

    unit_class = CosmicInfo
    metadata_class = CosmicMetadata

    @classmethod
    def from_vr(cls, vr, metadata=None):
        if metadata is None:
            metadata = cls.metadata_class.from_vcfheader(vr.header)

        result = cls.from_vr_base(vr)
        for unit in result:
            unit.metadata = metadata

        return result

    @classmethod
    def from_vcfspec(cls, vcfspec, cosmic_vcf, metadata=None, donot_init_metadata=False):
        if metadata is None:
            if not donot_init_metadata:
                metadata = cls.metadata_class.from_vcfheader(cosmic_vcf.header)

        cosmic_vr_list = customfile.fetch_relevant_vr_multialt(
            vcfspec, cosmic_vcf, search_equivs=True, raise_with_multihit=True,
        )
        result = cls()
        for vr in cosmic_vr_list:
            if vr is None:
                result.append(cls.unit_class.init_blank(metadata=metadata))
                #result.append(cls.unit_class.init_missing())
            else:
                result.extend(cls.from_vr(vr, metadata=metadata))

        return result

    def write(self, vr):
        self.write_base(vr)


def update_vcfheader(vcfheader, cosmic_vcf):
    CosmicInfoALTlist.add_meta(vcfheader)
    cosmicmeta = CosmicMetadata.from_vcfheader(cosmic_vcf.header)
    cosmicmeta.write(vcfheader)


#def annotate_vp(vp, cosmic_vcf, fasta, metadata=None):
#    vp.cosmic = CosmicInfoALTlist.from_vcfspec(vp.vcfspec, cosmic_vcf, fasta, metadata=metadata)



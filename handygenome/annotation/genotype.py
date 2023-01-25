import itertools

import Bio.Align

import handygenome.common as common
import handygenome.annotation.annotitem as annotitem


class GenotypeInfo(annotitem.AnnotItemInfoSingle):
    @classmethod
    def from_vp(
        cls, vp, sampleinfo, 
        germline_CNt=2,
        cutoff_somatic_vaf={
            'blood': 0.05, 'tumor': 0.05, 'clonal': 0.2,
        },
    ):
        """Args:
            sampleinfo: 
                format: dict(
                    sampleid=dict(
                        "type"="clonal"/"blood"/"tumor", 
                        "CNt"=<int>,
                    )
                )
                Assumes all samples are derived from the same individual
        """
        pass
        # set params
        possible_allele_indexes = set(range(len(vp.vcfspec.alts) + 1))
        sampletypes = {'clonal': list(), 'blood': list(), 'tumor': list()}
        for sampleid, subdic in sampleinfo.items():
            sampletypes[subdic['type']].append(sampleid)

        if (
            len(sampletypes['clonal']) == 0 and 
            len(sampletypes['blood']) == 0
        ):
            raise Exception(f'No blood or clonal samples; cannot determine germline genotype')

        # determine germline
        allele_index_polls = dict()
        for sampleid in itertools.chain(sampletypes['clonal'], sampletypes['blood']):
            rppcount_dict = vp.readstats_dict[sampleid]['rppcounts']
            allele_rppcounts_sorted = sorted(
                (
                    (allele_index, rppcount_dict[allele_index])
                    for allele_index in possible_allele_indexes
                ),
                key=(lambda x: x[1]),
                reverse=True,
            )
            for allele_index, rppcount in allele_rppcounts_sorted[:germline_CNt]:
                allele_index_polls.setdefault(allele_index, 0)
                allele_index_polls[allele_index] += 1

        germline_allele_indexes = sorted(
            allele_index_polls.items(),
            key=(lambda x: x[1]),
            reverse=True,
        )
        germline_allele_indexes = [
            x[0] for x in germline_allele_indexes[:germline_CNt]
        ]

        # determine somatic
        possible_somatic_indexes = possible_allele_indexes.difference(germline_allele_indexes)
        somatic_allele_indexes = list()

        for sample_type, sampleid_list in sampletypes.items():
            for sampleid in sampleid_list:
                readstats = vp.readstats_dict[sampleid]
                for allele_index in possible_somatic_indexes:
                    if allele_index in somatic_allele_indexes:
                        continue
                    vaf = readstats.get_vaf(alleleclass=allele_index)
                    if vaf >= cutoff_somatic_vaf[sample_type]:
                        somatic_allele_indexes.append(allele_index)
            
        # prepare result
        result = cls(is_missing=False)
        result['germline_allele_indexes'] = germline_allele_indexes
        result['somatic_allele_indexes'] = somatic_allele_indexes

        return result

    def get_possible_allele_indexes(self):
        return set(
            itertools.chain(
                self['germline_allele_indexes'],
                self['somatic_allele_indexes'],
            )
        )

    def get_clonal_genotype(self, vp, sampleid, CNt):
        """Assumes the sample is clonal"""
        readstats = vp.readstats_dict[sampleid]

        possible_allele_indexes = self.get_possible_allele_indexes()
        rppcount_dict = {
            allele_index: readstats['rppcounts'][allele_index]
            for allele_index in possible_allele_indexes
        }
        rppcount_sum = sum(rppcount_dict.values())

        allele_vafs_sorted = sorted(
            (
                (allele_index, rppcount / rppcount_sum)
                for allele_index, rppcount in rppcount_dict.items()
            ),
            key=(lambda x: x[1]),
            reverse=True,
        )

        vaf_cutoff = (1 / CNt) * 0.4
        allele_vafs_sorted = [
            x for x in allele_vafs_sorted
            if x[1] >= vaf_cutoff
        ]

        return [x[0] for x in allele_vafs_sorted[:CNt]]
        


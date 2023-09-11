import itertools

import pandas as pd

import handygenome.tools as tools
import handygenome.network as network
import handygenome.annotation.annotitem as annotitem
import handygenome.variant.vcfspec as libvcfspec


URL_allCuratedGenes = 'https://www.oncokb.org/api/v1/utils/allCuratedGenes'
URL_cancerGeneList = 'https://www.oncokb.org/api/v1/utils/cancerGeneList'
URL_HGVSG = 'https://www.oncokb.org/api/v1/annotate/mutations/byHGVSg'


class OncoKBInfo(annotitem.AnnotItemInfoSingle):
    pass


class OncoKBInfoALTlist(annotitem.AnnotItemInfoALTlist):
    meta = {
        "ID": "oncokb",
        "Number": "A",
        "Type": "String",
        "Description": "OncoKB data annotation, acquired throught web API(https://www.oncokb.org/swagger-ui/index.html).",
    }
    unit_class = OncoKBInfo

    @classmethod
    def from_vr(cls, vr):
        return cls.from_vr_base(vr)

    def write(self, vr):
        self.write_base(vr)


def get_allCuratedGenes(token):
    return pd.DataFrame.from_dict(
        network.http_get(
            url=URL_allCuratedGenes,
            headers={
                'Authorization': f'Bearer {token}',
                'accept': 'application/json',
            },
            params={
                'includeEvidence': True,
            },
        )
    )


def get_cancerGeneList(token):
    return pd.DataFrame.from_dict(
        network.http_get(
            url=URL_cancerGeneList,
            headers={
                'Authorization': f'Bearer {token}',
                'accept': 'application/json',
            },
        )
    )


def query_hgvsg(hgvsg, token, tumor_type=None, evidence_types=None):
    """Args:
        hgvsg: e.g. '7:g.140453136A>T'
    """
    params={
        'hgvsg': hgvsg,
    }
    if tumor_type is not None:
        params['tumorType'] = tumor_type
    if evidence_types is not None:
        params['evidenceType'] = ','.join(evidence_types)

    result = OncoKBInfo.init_nonmissing()
    result.update_dict(
        network.http_get(
            url=URL_HGVSG,
            params=params,
            headers={
                'Authorization': f'Bearer {token}',
                'accept': 'application/json',
            },
        )
    )
    return result


def query_hgvsg_post(hgvsg_list, token, tumor_type=None, evidence_types=None, chunk_size=50):
    result = list()
    #NR = 0
    for hgvsg_chunk in tools.chunk_iter(hgvsg_list, chunk_size):
        #NR += 1
        #print(f'{NR * chunk_size} entries being processed')
        data = list()
        for hgvsg in hgvsg_chunk:
            dic = {'hgvsg': hgvsg}
            if tumor_type is not None:
                dic['tumorType'] = tumor_type
            if evidence_types is not None:
                dic['evidenceTypes'] = evidence_types
            data.append(dic)

        for item in network.http_post(
            url=URL_HGVSG,
            data=data,
            headers={
                'Authorization': f'Bearer {token}',
                'accept': 'application/json',
                'Content-Type': 'application/json',
            },
        ):
            unit_oncokb = OncoKBInfo.init_nonmissing()
            unit_oncokb.update_dict(item)
            result.append(unit_oncokb)

    return result


def add_oncokb_info(vp_list, token, **kwargs):
    vp_list = list(vp_list)
    hgvsg_list = list()
    vp_indexes = list()
    for vp_idx, vp in enumerate(vp_list):
        for sub_vcfspec in vp.vcfspec.iter_monoalts():
            hgvsg_list.append(sub_vcfspec.to_hgvsg())
            vp_indexes.append(vp_idx)

    query_result = query_hgvsg_post(hgvsg_list=hgvsg_list, token=token, **kwargs)
    for key, subiter in itertools.groupby(
        zip(vp_indexes, query_result),
        key=(lambda x: x[0]),
    ):
        oncokb_ALTlist = OncoKBInfoALTlist()
        for vp_idx, oncokb_item in subiter:
            oncokb_ALTlist.append(oncokb_item)
        vp_list[vp_idx].oncokb = oncokb_ALTlist


def make_OncoKBInfoALTlist_list(vr_iterator, token, **kwargs):
    hgvsg_list = list()
    vr_indexes = list()
    for vr_idx, vr in enumerate(vr_iterator):
        vcfspec = libvcfspec.Vcfspec.from_vr(vr)
        for sub_vcfspec in vcfspec.iter_annotation_forms():
            hgvsg_list.append(sub_vcfspec.to_hgvsg())
            vr_indexes.append(vr_idx)

    result = list()
    query_result = query_hgvsg_post(hgvsg_list=hgvsg_list, token=token, **kwargs)
    for key, subiter in itertools.groupby(
        zip(vr_indexes, query_result),
        key=(lambda x: x[0]),
    ):
        oncokb_ALTlist = OncoKBInfoALTlist()
        for vr_idx, oncokb_item in subiter:
            oncokb_ALTlist.append(oncokb_item)
        result.append(oncokb_ALTlist)

    return result

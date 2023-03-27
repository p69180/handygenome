import itertools

import pandas as pd
import pyranges as pr

import handygenome.common as common


CURATED_BLACKLIST = common.RefverDict({
    'GRCh37': {
        'peri_centromere': [
            ('10', 38_000_000, 42_523_247),
            ('10', 42_527_152, 42_546_687),
            # 42_546_687 - 42_596_687 : N
            ('10', 42_596_687, 42_818_502),

            ('7', 61_967_157, 62_010_018),
            ('7', 62_015_012, 62_054_190),

            ('11', 51_583_335, 51_594_205),

            ('Y', 13_400_000, 13_573_805),
            ('Y', 13_634_339, 13_748_578),

            ('8', 43_092_656, 43_100_000),

            ('2', 89_830_436, 89_880_102),
        ],
        'non_centromere': [
            ('1', 144_810_724, 145_382_727),
            ('2', 33_091_515, 33_092_197),
            ('2', 33_093_197, 33_095_590),
            ('2', 33_141_000, 33_141_692),
        ],
    },
})


def make_blacklist_gr(data, refver):
    df_data = list(itertools.chain.from_iterable(data.values()))
    df = pd.DataFrame.from_records(df_data, columns=['Chromosome', 'Start', 'End'])
    assert set(df['Chromosome']).issubset(common.DEFAULT_CHROMDICTS[refver].contigs)
    gr = pr.PyRanges(df).merge().sort()

    return gr


CURATED_BLACKLIST_GRS = common.RefverDict(
    {refver: make_blacklist_gr(data, refver) for refver, data in CURATED_BLACKLIST.items()}
)

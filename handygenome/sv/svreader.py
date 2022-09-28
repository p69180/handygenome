import pandas as pd

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
libcircos = importlib.import_module('.'.join([top_package_name, 'plot', 'circos']))


def previous_sv_reader(sv_path):
    df = pd.read_table(sv_path, sep='\t', header=0, na_values=['.'])
    df = df.rename(
        columns={
            '#CHR1': 'Chromosome_site1',
            'POS1': 'Pos0_site1',
            'CHR2': 'Chromosome_site2',
            'POS2': 'Pos0_site2',
        }
    )
    df['Pos0_site1'] = df['Pos0_site1'] - 1
    df['Pos0_site2'] = df['Pos0_site2'] - 1
    df['Color'] = [libcircos.SV_COLORMAP[x] for x in df['svtype']]

    return df

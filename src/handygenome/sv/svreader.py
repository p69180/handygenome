import numpy as np
import pandas as pd

import handygenome.tools as tools
import handygenome.plot.circos as libcircos


def previous_sv_reader(sv_path):
    data = list()
    for linedict in tools.fileiter(
        sv_path, sep='\t', remove_leading_hashes=True, skip_double_hashes=True,
    ):
        keys = list(linedict.keys())
        for key in keys:
            if linedict[key] == ".":
                linedict[key] = np.nan
        data.append(linedict)

    df = pd.DataFrame.from_records(data)
    #df = pd.read_table(sv_path, sep='\t', header=0, na_values=['.'])
    df.rename(
        columns={
            'CHR1': 'CHROM_site1',
            'POS1': 'POS0_site1',
            'CHR2': 'CHROM_site2',
            'POS2': 'POS0_site2',
        },
        inplace=True,
    )

    df['POS0_site1'] = df['POS0_site1'].astype('float') - 1
    df['POS0_site2'] = df['POS0_site2'].astype('float') - 1

    svtype_keys = ('svtype', 'SVtype')
    svtype_col = None
    for svtype_key in svtype_keys:
        try:
            svtype_col = df[svtype_key]
        except KeyError:
            continue
        else:
            break
    if svtype_col is None:
        raise Exception(f'Input file does not have a svtype column name.')

    df['Color'] = [libcircos.SV_COLORMAP[x] for x in svtype_col]

    return df


def add_color_to_sv_df(df):
    if 'svtype' in df.columns:
        colors = list()
        for x in df['svtype']:
            try:
                col = libcircos.SV_COLORMAP[x]
            except KeyError:
                col = libcircos.SV_COLORMAP[np.nan]
            colors.append(col)
    elif 'CT' in df.columns:
        colors = list()
        for idx, row in df.iterrows():
            if row['CHROM_site1'] != row['CHROM_site2']:
                colors.append(libcircos.SV_COLORMAP['TRA'])
            else:
                if row['CT'] == '3to5':
                    colors.append(libcircos.SV_COLORMAP['DEL'])
                elif row['CT'] == '5to3':
                    colors.append(libcircos.SV_COLORMAP['DUP'])
                elif row['CT'] in ('3to3', '5to5'):
                    colors.append(libcircos.SV_COLORMAP['INV'])
                else:
                    raise Exception(f'Unexpected CT value:\n{row}')
    else:
        raise Exception(f'Input dataframe columns does not include "svtype" nor "CT".')

    df.loc[:, 'Color'] = colors


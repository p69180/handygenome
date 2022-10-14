import re
import functools
import itertools

import pysam
import pyranges as pr
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import Bio.Seq

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
libsigmisc = importlib.import_module('.'.join([top_package_name, 'signature', 'misc']))
#sequenza_loader = importlib.import_module('.'.join([top_package_name, 'cnv', 'sequenza_loader']))


"""Issues with pyranges.PyRanges.intersect:
    - "Start" and "End" columns of pyranges.PyRanges object are by 
        default "int32".
    - "intersect" and "join" methods results in error if "Start" and "End"
        columns of "self" gr and "other" gr are of different types.
    - When loading breakends data with missing values, 
        - pd.DataFrame is created with POS1/POS2 columns with nan.
        - After removing nan, the columns are coerced into integer type.
        - At that time, dtype must be "int32".
"""

SV_COLORMAP = {
    'DEL': 'red',
    'DUP': 'green',
    'INV': 'aquamarine',
    'TRA': 'gold',
    np.nan: 'black',
}

CYTOBAND_COLORMAP = {
    "gneg": "#FFFFFF00", 
    "gpos25": "#EEEEEE", 
    "gpos50": "#BBBBBB", 
    "gpos75": "#777777", 
    "gpos100": "#000000", 
    "gvar": "#FFFFFF00", 
    "stalk": "#C01E27", 
    "acen": "#D82322",
}

CYTOBAND_PATHS = common.RefverDict({
    'GRCh37': '/home/users/pjh/References/cytoband/cytoBand_hg19.bed',
    'GRCh38': '/home/users/pjh/References/cytoband/cytoBand_hg38.bed',
})


def read_cytoband(cytoband_path, refver):
    result = pd.read_table(
        cytoband_path, 
        header=None,
        names=['Chromosome', 'Start', 'End', 'Name', 'Stain']
    )
    result['Color'] = result['Stain'].apply(lambda x: CYTOBAND_COLORMAP[x])

    refver = common.RefverDict.converter[refver]
    if refver == 'GRCh37':
        result['Chromosome'] = result['Chromosome'].apply(
            lambda x: re.sub('^chr', '', x)
        )

    return result


def check_column_existence(gr, names):
    for key in names:
        if key not in gr.columns:
            raise Exception(f'"{key}" must be included in columns')


def linspace_rounded(start, stop, num):
    width = (stop - start) / num
    decimals = int(-np.log10(width))
    return np.linspace(start, stop, num).round(decimals)


class Sector:
    def __init__(self, chrom, start0, end0):
        self.chrom = chrom
        self.start0 = start0
        self.end0 = end0

    def __repr__(self):
        return f'<Sector object(chrom={self.chrom}, start0={self.start0}, end0={self.end0}>'

    @property
    def id(self):
        return (self.chrom, self.start0, self.end0)

    @property
    def length(self):
        return self.end0 - self.start0

    def get_label(self, short_for_full_chromosome=True, chromdict=None):
        if short_for_full_chromosome and (chromdict is None):
            raise Exception(f'If "short_for_full_chromosome" argument is True, "chromdict" argument must be set.')

        if short_for_full_chromosome:
            if self.start0 == 0 and self.end0 == chromdict[self.chrom]:
                return self.chrom
            else:
                return f'{self.chrom}_{self.start0}_{self.end0}'
        else:
            return f'{self.chrom}_{self.start0}_{self.end0}'


class SectorList(list):
    def __init__(self, chroms=list(), start0s=list(), end0s=list()):
        for chrom, start0, end0 in zip(chroms, start0s, end0s):
            self.append(Sector(chrom, start0, end0))

    def add_sector(self, chrom, start0, end0):
        self.append(Sector(chrom, start0, end0))

    @classmethod
    def from_chromdict(cls, chromdict, assembled_only=True):
        chroms = list()
        start0s = list()
        end0s = list()
        for chrom, length in chromdict.items():
            if assembled_only:
                if not common.RE_PATS['assembled_chromosome'].fullmatch(chrom):
                    continue

            chroms.append(chrom)
            start0s.append(0)
            end0s.append(length)

        return cls(chroms, start0s, end0s)

    @classmethod
    def from_refver(cls, refver, assembled_only=True):
        chromdict = common.ChromDict(refver=refver)
        return cls.from_chromdict(chromdict, assembled_only=assembled_only)

    @classmethod
    def from_gr(cls, gr):
        chroms = list()
        start0s = list()
        end0s = list()
        for idx, row in gr.df.iterrows():
            chroms.append(row['Chromosome'])
            start0s.append(row['Start'])
            end0s.append(row['End'])
        return cls(chroms, start0s, end0s)

    @functools.cached_property
    def dict(self):
        return {x.id: x for x in self}

    @functools.cached_property
    def gr(self):
        chroms = list()
        starts = list()
        ends = list()
        ids = list()
        
        for x in self:
            chroms.append(x.chrom)
            starts.append(x.start0)
            ends.append(x.end0)
            ids.append(x.id)

        result = pr.from_dict(
            {
                'Chromosome': chroms,
                'Start': starts,
                'End': ends,
                'ID': ids,
            },
            int64=False,
        )

        return result

    def get_ids(self):
        return [x.id for x in self]


class Circos:
    @common.get_deco_num_set(('refver', 'chromdict'), 1)
    def __init__(self, refver=None, chromdict=None, sectorlist=None, ylim=(0, 1), sector_interspace_deg=2, figsize=(10, 10), angle_range_deg=(0, 360), theta_offset_deg=90, theta_direction=-1):
        """Args:
            sector_interspace_deg: A scalar or a sequence type with the same 
                length as sectorlist. Each element represents the space BEFORE
                to the sector of the same index.
        """
        # chromdict
        if refver is not None:
            self.chromdict = common.ChromDict(refver=refver)
        elif chromdict is not None:
            self.chromdict = chromdict
        # sectorlist
        if sectorlist is None:
            self.sectorlist = SectorList.from_chromdict(self.chromdict, assembled_only=True)
        else:
            self.sectorlist = sectorlist
        # sector interspace
        self.sector_interspace_deg = self.arghandler_sector_interspace_deg(sector_interspace_deg)
        self.sector_interspace_rad = np.deg2rad(sector_interspace_deg)
        # others 
        self.ylim = ylim
        self.figsize = figsize
        self.angle_range_deg = self.get_angle_range(angle_range_deg)
        self.angle_length_deg = self.angle_range_deg[1] - self.angle_range_deg[0]
        self.theta_offset_deg = theta_offset_deg
        self.theta_direction = theta_direction
        # postprocess
        #self.init_figure()
        self.set_sector_angles()

    @staticmethod
    def get_angle_range(raw_angle_range_deg):
        if raw_angle_range_deg[0] == raw_angle_range_deg[1]:
            raise Exception(f'Start and end of "angle_range" must be different.')

        start_angle = np.mod(raw_angle_range_deg[0], 360)
        end_angle = np.mod(raw_angle_range_deg[1], 360)
        if start_angle == end_angle:
            return (start_angle, end_angle + 360)
        elif end_angle > start_angle:
            return (start_angle, end_angle)
        elif end_angle < start_angle:
            return (start_angle, end_angle + 360)

    def set_sector_angles(self):
        sector_interspace_deg_sum = self.sector_interspace_deg.sum()
        if sector_interspace_deg_sum > self.angle_length_deg:
            raise Exception(f'Sum of interspaces({sector_interspace_deg_sum}) is greater than the total circos angle({self.angle_length_deg}). Choose a smaller "sector_interspace_deg".')

        sector_length_angle_deg_sum = self.angle_length_deg - sector_interspace_deg_sum
        sector_length_sum = sum(x.length for x in self.sectorlist)

        cursor_deg = self.angle_range_deg[0]
        sector_angles_rad = dict()
        #sector_angles_deg = dict()
        for sector, interspace_deg in zip(self.sectorlist, self.sector_interspace_deg):
            # now cursor_deg indicates "sector start position" - "interspace"
            sector_start_angle_deg = cursor_deg + interspace_deg
            sector_end_angle_deg = sector_start_angle_deg + sector_length_angle_deg_sum * (sector.length / sector_length_sum)

            sector_angle_pair_deg = (sector_start_angle_deg, sector_end_angle_deg)
            sector_angle_pair_rad = tuple(np.deg2rad(sector_angle_pair_deg))
            #sector_angles_deg[sector.id] = sector_angle_pair_deg
            sector_angles_rad[sector.id] = sector_angle_pair_rad

            cursor_deg = sector_end_angle_deg

        self.sector_angles_rad = sector_angles_rad
        #self.sector_angles_deg = sector_angles_deg

    # drawing methods
#    def init_figure(self):  # deprecated
#        self.figure = plt.figure(figsize=self.figsize)
#        self.ax = self.figure.add_axes([0, 0, 1, 1], polar=True)
#
#        self.ax.set_theta_offset(np.deg2rad(self.theta_offset_deg))
#        self.ax.set_theta_direction(self.theta_direction)
#        self.ax.set_xlim(*np.deg2rad(self.angle_range_deg))
#        self.ax.set_ylim(self.ylim[0], self.ylim[1])
#        self.ax.spines['polar'].set_visible(False)
#        self.ax.xaxis.set_ticks([])
#        self.ax.xaxis.set_ticklabels([])
#        self.ax.yaxis.set_ticks([])
#        self.ax.yaxis.set_ticklabels([])
#
#        plt.close()

    def init_plot(self, ax):
        ax.set_theta_offset(np.deg2rad(self.theta_offset_deg))
        ax.set_theta_direction(self.theta_direction)
        ax.set_xlim(*np.deg2rad(self.angle_range_deg))
        ax.set_ylim(self.ylim[0], self.ylim[1])
        ax.spines['polar'].set_visible(False)
        ax.xaxis.set_ticks([])
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticks([])
        ax.yaxis.set_ticklabels([])

        plt.close()

    def draw_sector_titles(self, ax, short_for_full_chromosome=True, radial_pos=1.15, fontsize=20, perpendicular=False):
        for sector in self.sectorlist: 
            label = sector.get_label(short_for_full_chromosome=short_for_full_chromosome, chromdict=self.chromdict)
            sector_angle_pair_rad = self.sector_angles_rad[sector.id]
            position_rad = (sector_angle_pair_rad[0] + sector_angle_pair_rad[1]) / 2
            position_deg = self.to_canonical_degrees(position_rad)
            if perpendicular:
                rot = self.get_rotation_perpendicular(position_deg)
            else:
                rot = self.get_rotation_upright(position_deg)
            # draw text
            ax.text(position_rad, radial_pos, label, rotation=rot, ha='center', va='center', fontsize=fontsize)

    def draw_poslabels_within_sector(self, ax, interspace_deg=30, radial_pos=1.10, fontsize=10):
        tick_angles_deg = np.arange(start=self.angle_range_deg[0], stop=self.angle_range_deg[1], step=interspace_deg)
        tick_angles_rad = np.deg2rad(tick_angles_deg)
        for angle_rad in tick_angles_rad:
            data_position = self.get_data_position(angle_rad)
            if data_position is None:
                continue
            label = f'{data_position:,}'
            canonical_degree = self.to_canonical_degrees(angle_rad)
            rot = self.get_rotation_perpendicular(canonical_degree)
            ax.text(angle_rad, radial_pos, label, rotation=rot, ha='center', va='center', fontsize=fontsize)

    def draw_cytoband(self, ax, refver, radial_range=(0.94, 1)):
        cytoband_df = read_cytoband(CYTOBAND_PATHS[refver], refver)
        self.draw_frame(ax, radial_range)
        self.draw_bar(ax, cytoband_df, radial_range)

    def draw_point_mutation(self, ax, data_df, radial_range=(0.8, 0.94), size=5, draw_gridlabel=True):
        def make_colors(df):
            colormap = libsigmisc.COLORS_SBS6
            colorlist = list()
            for idx, row in df.iterrows():
                ref = row['REF']
                alt = row['ALT']
                if len(ref) == len(alt) and len(ref) == 1:
                    if ref in 'CT':
                        color = colormap[ref + '>' + alt]
                    elif ref in 'AG':
                        ref = Bio.Seq.complement(ref)
                        alt = Bio.Seq.complement(alt)
                        color = colormap[ref + '>' + alt]
                    else:
                        color = colormap['other']
                else:
                    color = colormap['other']

                colorlist.append(color)

            return colorlist

        # main
        where_to_draw_label = self.arghandler_draw_gridlabel(draw_gridlabel)
        check_column_existence(data_df, ('CHROM', 'POS', 'REF', 'ALT', 'VAF'))

        data_df['Chromosome'] = data_df['CHROM']
        data_df['Pos0'] = data_df['POS'] - 1
        data_df['Y'] = data_df['VAF']
        data_df['Color'] = make_colors(data_df)

        ymin = 0
        ymax = 1
        gridline_ys = [0.25, 0.5, 0.75]
        gridlabel_ys = [0, 0.5, 1] if where_to_draw_label.any() else None
        self.draw_frame(ax, radial_range, gridline_ys=gridline_ys, gridlabel_ys=gridlabel_ys, ymin=ymin, ymax=ymax, where_to_draw_label=where_to_draw_label)
        self.draw_scatter(ax, data_df, radial_range, ymin=ymin, ymax=ymax, size=size)

    def draw_copynumber(self, ax, data_df, radial_range, color_A='red', color_B='blue', ymin=None, ymax=None, ymax_quantile=0.95, y_offset=0.1, linewidth=0.5, draw_gridlabel=True):
        where_to_draw_label = self.arghandler_draw_gridlabel(draw_gridlabel)
        check_column_existence(data_df, ('Chromosome', 'Start', 'End', 'A', 'B'))
        # filter out NaN rows
        data_df = data_df.loc[np.logical_not(data_df['A'].isnull()), :]
        data_df = data_df.loc[np.logical_not(data_df['B'].isnull()), :]
        # A
        data_df_A = data_df.copy()
        data_df_A['Y'] = data_df_A['A']
        zero_rows = data_df_A['Y'] == 0
        data_df_A.loc[zero_rows, 'Y'] += 2 * y_offset
        data_df_A.loc[np.logical_not(zero_rows), 'Y'] += y_offset
        data_df_A['Color'] = color_A
        # B
        data_df_B = data_df.copy()
        data_df_B['Y'] = data_df_B['B']
        zero_rows = data_df_A['Y'] == 0
        data_df_B.loc[zero_rows, 'Y'] += y_offset
        data_df_B.loc[np.logical_not(zero_rows), 'Y'] -= y_offset
        data_df_B['Color'] = color_B
        # set limits
        if ymin is None:
            ymin = 0
        if ymax is None:
            data_df_modified = self.modify_data(data_df, make_radius=False)
            all_ys = pd.Series(itertools.chain(data_df_modified['A'], data_df_modified['B']))
            all_ys = all_ys.loc[np.logical_not(all_ys.isnull())]
            ymax = np.quantile(all_ys, ymax_quantile)
            ymax = max(ymax, 6)
        # draw frame
        gridline_ys = np.arange(start=np.rint(ymin), stop=np.rint(ymax), step=1, dtype='int')
        if where_to_draw_label.any():
            gridlabel_ys = np.arange(start=np.rint(ymin), stop=np.rint(ymax), step=2, dtype='int')
        else:
            gridlabel_ys = None

        self.draw_frame(ax, radial_range, gridline_ys=gridline_ys, gridlabel_ys=gridlabel_ys, ymin=ymin, ymax=ymax, where_to_draw_label=where_to_draw_label)
        # draw arcs
        self.draw_arc(ax, data_df_A, radial_range, linewidth=linewidth, ymin=ymin, ymax=ymax)
        self.draw_arc(ax, data_df_B, radial_range, linewidth=linewidth, ymin=ymin, ymax=ymax)

    def draw_depth_scatter(self, ax, data_df, radial_range, color='black', ymin=None, ymax=None, ymax_quantile=0.99, size=1, gridline_num=8, draw_gridlabel=True):
        where_to_draw_label = self.arghandler_draw_gridlabel(draw_gridlabel)
        check_column_existence(data_df, ('Chromosome', 'Start', 'End', 'Depth'))

        data_df['Pos0'] = ((data_df['Start'] + data_df['End']) * 0.5).apply(int)
        data_df['Y'] = data_df['Depth']
        data_df['Color'] = color
        # set limits
        if ymin is None:
            ymin = 0
        if ymax is None:
            data_df_modified = self.modify_data(data_df, make_radius=False)
            ymax = np.quantile(data_df_modified['Depth'], ymax_quantile)
            ymax = max(ymax, 60)
        # draw frame
        gridline_start = np.floor(ymin)
        gridline_stop = np.ceil(ymax)
        gridline_ys = np.rint(np.linspace(start=gridline_start, stop=gridline_stop, num=gridline_num, endpoint=True)).astype('int')
        if where_to_draw_label.any():
            gridlabel_ys = gridline_ys[::2]
        else:
            gridlabel_ys = None

        self.draw_frame(ax, radial_range, gridline_ys=gridline_ys, gridlabel_ys=gridlabel_ys, ymin=ymin, ymax=ymax, where_to_draw_label=where_to_draw_label)
        # draw dots
        self.draw_scatter(ax, data_df, radial_range, ymin=ymin, ymax=ymax, size=size)

    draw_depth = draw_depth_scatter

    def draw_breakends(self, ax, data_df, radial_pos, linewidth=0.5):
        check_column_existence(data_df, ('CHROM_site1', 'POS0_site1', 'CHROM_site2', 'POS0_site2', 'Color'))
        self.draw_chord(ax, data_df, radial_pos, linewidth=linewidth)

    #################

    @staticmethod
    def get_radius_list(y, radial_range, ymin=None, ymax=None):
        y = np.array(y)

        if ymin is None:
            ymin = min(y)
        if ymax is None:
            ymax = max(y)

        y_range = ymax - ymin
        radial_length = radial_range[1] - radial_range[0]
        bottom = radial_range[0]
        dataspace_portion = (y - ymin) / y_range
        return bottom + radial_length * dataspace_portion

    @staticmethod
    def get_angle(data_position, sector, sector_start_angle_rad, sector_length_angle_rad):
        portion_within_sector = (data_position - sector.start0) / sector.length
        angle_within_sector_rad = sector_length_angle_rad * portion_within_sector
        return sector_start_angle_rad + angle_within_sector_rad
    
    def get_data_position(self, angle_rad):
        # find a corresponding sector id
        found_a_sector = False
        for sector_id, angle_pair_rad in self.sector_angles_rad.items():
            if angle_rad >= angle_pair_rad[0] and angle_rad < angle_pair_rad[1]:
                found_a_sector = True
                break
        if not found_a_sector:
            return None

        # calculate data position
        angle_within_sector_rad = angle_rad - angle_pair_rad[0]
        sector_length_angle_rad = angle_pair_rad[1] - angle_pair_rad[0]
        portion_within_sector = angle_within_sector_rad / sector_length_angle_rad
        sector = self.sectorlist.dict[sector_id]
        data_position = int(sector.start0 + portion_within_sector * sector.length)
        return data_position

    def to_canonical_degrees(self, angle_rad):
        return self.theta_offset_deg + np.rad2deg(angle_rad) * self.theta_direction

    @staticmethod
    def get_rotation_upright(canonical_degree):
        rot = canonical_degree - 90
        canonical_degree_mod = np.mod(canonical_degree, 360)
        if canonical_degree_mod > 180 and canonical_degree_mod < 360:
            rot += 180
        return rot

    @staticmethod
    def get_rotation_perpendicular(canonical_degree):
        rot = canonical_degree
        canonical_degree_mod = np.mod(canonical_degree, 360)
        if canonical_degree_mod > 90 and canonical_degree_mod < 270:
            rot += 180
        return rot

    def modify_data(self, data_df, make_radius=False, radial_range=None, ymin=None, ymax=None):
        def sanity_check(data_df, make_radius, radial_range):
            check_column_existence(data_df, ('Chromosome', 'Start', 'End'))
            if make_radius:
                if radial_range is None:
                    raise Exception(f'If "make_radius" is True, "radial_range" must be given.')
                if 'Y' not in data_df.columns:
                    raise Exception(f'If "make_radius" is True, "Y" column must be included in "data_df".')

        # sanity check
        sanity_check(data_df, make_radius, radial_range)
        # save original data_df columns
        columns = data_df.columns
        # perform join
        data_gr = pr.PyRanges(df=data_df, int64=False)
        gr_joined = data_gr.join(self.sectorlist.gr).new_position('intersection')
        if gr_joined.empty:
            result_columns = list(columns)
            result_columns.extend(['Start_rad', 'End_rad'])
            if make_radius:
                result_columns.append('Radius')
            result = pd.DataFrame(columns=result_columns)
        else:
            data_start_angles_rad = list()
            data_end_angles_rad = list()
            for idx, row in gr_joined.df.iterrows():
                sector_id = row['ID']
                sector = self.sectorlist.dict[sector_id]
                sector_angle_pair_rad = self.sector_angles_rad[sector_id]
                sector_start_angle_rad = sector_angle_pair_rad[0]
                sector_length_angle_rad = sector_angle_pair_rad[1] - sector_angle_pair_rad[0]

                data_start_angles_rad.append(
                    self.get_angle(row['Start'], sector, sector_start_angle_rad, sector_length_angle_rad)
                )
                data_end_angles_rad.append(
                    self.get_angle(row['End'], sector, sector_start_angle_rad, sector_length_angle_rad)
                )

            result = gr_joined.df
            result.insert(result.shape[1], 'Start_rad', data_start_angles_rad)
            result.insert(result.shape[1], 'End_rad', data_end_angles_rad)

            if make_radius:
                radius_list = self.get_radius_list(result['Y'], radial_range, ymin=ymin, ymax=ymax)
                result.insert(result.shape[1], 'Radius', radius_list)
                # filter out data out of radial_range
                result = result.loc[
                    (result['Radius'] <= radial_range[1]) & (result['Radius'] >= radial_range[0]),
                    :
                ]

        return result

    def modify_breakends_data(self, data_df):

        def set_site1_coords(data):
            if isinstance(data, pr.PyRanges):
                # without coercing into int32, "intersect" method results in error
                data = data.assign('Chromosome', lambda df: df['CHROM_site1'])
                data = data.assign('Start', lambda df: df['POS0_site1'].astype('int32'))
                data = data.assign('End', lambda df: (df['POS0_site1'] + 1).astype('int32'))
            elif isinstance(data, pd.DataFrame):
                data['Chromosome'] = data['CHROM_site1']
                data['Start'] = data['POS0_site1']
                data['End'] = data['POS0_site1'] + 1
            return data

        def set_site2_coords(data):
            if isinstance(data, pr.PyRanges):
                # without coercing into int32, "intersect" method results in error
                data = data.assign('Chromosome', lambda df: df['CHROM_site2'])
                data = data.assign('Start', lambda df: df['POS0_site2'].astype('int32'))
                data = data.assign('End', lambda df: (df['POS0_site2'] + 1).astype('int32'))
            elif isinstance(data, pd.DataFrame):
                data['Chromosome'] = data['CHROM_site2']
                data['Start'] = data['POS0_site2']
                data['End'] = data['POS0_site2'] + 1
            return data

        # sanity check
        check_column_existence(data_df, ('CHROM_site1', 'POS0_site1', 'CHROM_site2', 'POS0_site2', 'Color'))
        # remove nan rows - those whose bnd2 is unknown (e.g. telomere)
        data_df = data_df.loc[np.logical_not(data_df['POS0_site1'].isna()), :]
        data_df = data_df.loc[np.logical_not(data_df['POS0_site2'].isna()), :]
        # prepare data_df for turning into a pyranges object
        data_df['Chromosome'] = '1'
        data_df['Start'] = 0
        data_df['End'] = 1
        # initialize gr
        data_gr = pr.PyRanges(df=data_df, int64=False)
        # intersect with site1
        data_gr = set_site1_coords(data_gr)
        data_gr = data_gr.intersect(self.sectorlist.gr)
        if data_gr.empty:
            return pd.DataFrame()
        # intersect with site2
        data_gr = set_site2_coords(data_gr)
        data_gr = data_gr.intersect(self.sectorlist.gr)
        if data_gr.empty:
            return pd.DataFrame()
        # run modify_data with site2
        data_df = self.modify_data(data_gr.df, make_radius=False)
        data_df.insert(data_df.shape[1], 'Angle_site2_rad', data_df['Start_rad'])
        data_df = data_df.drop(columns=['Start_rad', 'End_rad', 'ID', 'Start_b', 'End_b'])
        # run modify_data with site1
        data_df = set_site1_coords(data_df)
        data_df = self.modify_data(data_df, make_radius=False)
        data_df.insert(data_df.shape[1], 'Angle_site1_rad', data_df['Start_rad'])
        data_df = data_df.drop(columns=['Start_rad', 'End_rad', 'ID', 'Start_b', 'End_b'])
        # finally drop unused columns
        data_df = data_df.drop(columns=['Chromosome', 'Start', 'End'])

        return data_df

    def draw_frame(self, ax, radial_range, gridline_ys=None, gridlabel_ys=None, gridlabels=None, gridlabel_size=None, ymin=None, ymax=None, linewidth_grid=0.5, linewidth_spine=1, color_spine='black', color_grid='black', where_to_draw_label=None):
        # sanity checks
        if (gridlabel_ys is not None) and (gridlabels is not None):
            if len(gridlabel_ys) != len(gridlabels):
                raise Exception(f'Lengths of "gridlabel_ys" and "gridlabels" must be equal.')

        # determine if grids will be drawn
        draw_grids = (gridline_ys is not None) or (gridlabel_ys is not None)

        # set default gridlabels
        if (gridlabel_ys is not None) and (gridlabels is None):
            gridlabels = [str(x) for x in gridlabel_ys]

        # set default ymin and ymax
        if draw_grids:
            all_ys = list(
                itertools.chain(
                    (list() if gridline_ys is None else gridline_ys),
                    (list() if gridlabel_ys is None else gridlabel_ys)
                )
            )
            if ymin is None:
                ymin = min(all_ys)
            if ymax is None:
                ymax = max(all_ys)

        # set default "where_to_draw_label"
        if where_to_draw_label is None:
            where_to_draw_label = np.repeat(True, len(self.sectorlist))

        # set other parameters
        spine_height = radial_range[1] - radial_range[0]
        spine_bottom = radial_range[0]

        if gridline_ys is None:
            #gridline_radii = None
            gridline_diameters = None
        else:
            #gridline_radii = self.get_radius_list(gridline_ys, radial_range, ymin=ymin, ymax=ymax)
            gridline_diameters = 2 * self.get_radius_list(gridline_ys, radial_range, ymin=ymin, ymax=ymax)

        if gridlabel_ys is None:
            gridlabel_radii = None
            #gridlabel_diameters = None
        else:
            gridlabel_radii = self.get_radius_list(gridlabel_ys, radial_range, ymin=ymin, ymax=ymax)
            #gridlabel_diameters = 2 * gridlabel_radii

        # draw
        for idx, (sector, draw_label) in enumerate(zip(self.sectorlist, where_to_draw_label)):
            # draw spines
            sector_start_angle_rad, sector_end_angle_rad = self.sector_angles_rad[sector.id]
            width = sector_end_angle_rad - sector_start_angle_rad
            ax.bar(x=sector_start_angle_rad, height=spine_height, width=width, bottom=spine_bottom, align='edge', fill=False, edgecolor=color_spine, linewidth=linewidth_spine)

            if draw_grids:
                # set thetas
                if self.theta_direction == 1:
                    theta1 = self.to_canonical_degrees(sector_start_angle_rad)
                    theta2 = self.to_canonical_degrees(sector_end_angle_rad)
                    start_theta = theta1
                elif self.theta_direction == -1:
                    theta1 = self.to_canonical_degrees(sector_end_angle_rad)
                    theta2 = self.to_canonical_degrees(sector_start_angle_rad)
                    start_theta = theta2
                # draw gridlines
                if gridline_diameters is not None:
                    for diameter in gridline_diameters:
                        arc = mpl.patches.Arc(
                            xy=(0, 0), width=diameter, height=diameter,
                            theta1=theta1, theta2=theta2, linewidth=linewidth_grid,
                            edgecolor=color_grid, transform=ax.transData._b,
                            zorder=0, alpha=0.25,
                        )
                        ax.add_patch(arc)
                # draw gridlabels
                if (gridlabel_radii is not None) and draw_label:
                    interspace_rad = self.sector_interspace_rad[idx]
                    text_angle = sector_start_angle_rad - 0.5 * interspace_rad
                    for radius, label in zip(gridlabel_radii, gridlabels):
                        rot = self.get_rotation_upright(start_theta)
                        ax.text(text_angle, radius, label, rotation=rot, ha='center', va='center', fontsize=gridlabel_size)
                        #print(sector_start_angle_rad, radius, label, rot)

    def draw_bar(self, ax, data_df, radial_range):
        check_column_existence(data_df, ('Chromosome', 'Start', 'End', 'Color',))
        # modify data
        data_df = self.modify_data(data_df, make_radius=False)
        # remove rows where Color value is missing
        data_df = data_df.loc[np.logical_not(data_df['Color'].isnull()), :]
        # prepare drawing parameters
        height = radial_range[1] - radial_range[0]
        bottom = radial_range[0]
        # draw
        for idx, row in data_df.iterrows():
            color = row['Color']
            x = row['Start_rad']
            width = row['End_rad'] - row['Start_rad']
            ax.bar(x=x, height=height, width=width, bottom=bottom, align='edge', color=color, linewidth=0)
            
    def draw_scatter(self, ax, data_df, radial_range, ymin=None, ymax=None, size=None, alpha=None):
        check_column_existence(data_df, ('Chromosome', 'Pos0', 'Y', 'Color'))
        # modify data
        data_df['Start'] = data_df['Pos0']
        data_df['End'] = data_df['Pos0'] + 1
        data_df = self.modify_data(data_df, make_radius=True, radial_range=radial_range, ymin=ymin, ymax=ymax)
        # remove rows where Y or Color value is missing
        data_df = data_df.loc[np.logical_not(data_df['Y'].isnull()), :]
        data_df = data_df.loc[np.logical_not(data_df['Color'].isnull()), :]
        # prepare drawing parameters
        x = data_df['Start_rad']
        y = data_df['Radius']
        colors = data_df['Color']
        # draw
        ax.scatter(x=x, y=y, s=size, c=colors, alpha=alpha)

    def draw_arc(self, ax, data_df, radial_range, linewidth=None, ymin=None, ymax=None):
        check_column_existence(data_df, ('Chromosome', 'Start', 'End', 'Y', 'Color'))

        data_df = self.modify_data(data_df, make_radius=True, radial_range=radial_range, ymin=ymin, ymax=ymax)
        # remove rows where Y or Color value is missing
        data_df = data_df.loc[np.logical_not(data_df['Y'].isnull()), :]
        data_df = data_df.loc[np.logical_not(data_df['Color'].isnull()), :]

        # prepare drawing parameters
        if self.theta_direction == 1:
            theta1_list = self.to_canonical_degrees(data_df['Start_rad'])
            theta2_list = self.to_canonical_degrees(data_df['End_rad'])
        elif self.theta_direction == -1:
            theta1_list = self.to_canonical_degrees(data_df['End_rad'])
            theta2_list = self.to_canonical_degrees(data_df['Start_rad'])

        diameter_list = data_df['Radius'] * 2
        color_list = data_df['Color']

        # draw
        for theta1, theta2, diameter, color in zip(
            theta1_list, theta2_list, diameter_list, color_list
        ):
            arc = mpl.patches.Arc(
                xy=(0, 0), width=diameter, height=diameter,
                theta1=theta1, theta2=theta2, linewidth=linewidth,
                edgecolor=color, transform=ax.transData._b,
            )
            ax.add_patch(arc)

    def draw_chord(self, ax, data_df, radial_pos, linewidth=None):

        def make_patch(site1_rad, site2_rad, radial_pos, color, linewidth):
            site1_vert = (site1_rad, radial_pos)
            site2_vert = (site2_rad, radial_pos)
            control_vert = (0, 0)
            verts = [
                site1_vert,
                control_vert,
                site2_vert,
            ]
            codes = [
                mpl.path.Path.MOVETO, 
                mpl.path.Path.CURVE3, 
                mpl.path.Path.CURVE3,
            ]
            patch = mpl.patches.PathPatch(
                mpl.path.Path(verts, codes), fill=False, edgecolor=color,
                linewidth=linewidth,
            )
            return patch

        # main
        check_column_existence(data_df, ('CHROM_site1', 'POS0_site1', 'CHROM_site2', 'POS0_site2', 'Color'))
        data_df = self.modify_breakends_data(data_df)
        for idx, row in data_df.iterrows():
            patch = make_patch(row['Angle_site1_rad'], row['Angle_site2_rad'], radial_pos, row['Color'], linewidth=linewidth)
            ax.add_patch(patch)

    def arghandler_draw_gridlabel(self, arg):
        if arg is True:
            return np.repeat(True, len(self.sectorlist))
        elif arg is False:
            return np.repeat(False, len(self.sectorlist))
        else:
            if len(arg) != len(self.sectorlist):
                raise Exception(f'When "draw_gridlabel" is a sequence, its length must be equal to the number of sectors.')
            if not all(isinstance(x, int) for x in arg):
                raise Exception(f'When "draw_gridlabel" is a sequence, all of its elements must be integers.')
            result = np.repeat(False, len(self.sectorlist))
            result[list(arg)] = True
            return result

    def arghandler_sector_interspace_deg(self, arg):
        if np.isscalar(arg):
            return np.repeat(arg, len(self.sectorlist))
        else:
            if len(arg) != len(self.sectorlist):
                raise Exception(f'When "sector_interspace_deg" is a sequence, its length must be equal to the number of sectors.')
            return np.array(arg)



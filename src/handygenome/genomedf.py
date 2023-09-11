import os
import datetime
import inspect
import itertools
import contextlib

import pandas as pd
from pandas.util import hash_pandas_object
import numpy as np
import pyranges as pr

import handygenome.deco as deco
import handygenome.refgenome.refgenome as refgenome
import handygenome.tools as tools
import handygenome.cnv.rdnacopy as rdnacopy


COMMON_COLUMNS = ['Chromosome', 'Start', 'End']


class GenomeDataFrame:
    chrom_arg_keys = ['chrom', 'chroms', 'chromosome', 'chromosomes']
    start_arg_keys = ['start', 'starts', 'start0', 'start0s']
    end_arg_keys = ['end', 'ends', 'end0', 'end0s']

    MERGETYPES_SIMPLEAGG = ['mean', 'std', 'first', 'last', 'max', 'min', 'median', 'skew']
    MERGETYPES_WEIGHTED_SIMPLEAGG = ['weighted_mean']
    MERGETYPES_BYLENGTH = ['longest', 'shortest']
    AVAILABLE_MERGE_METHODS = (
        MERGETYPES_SIMPLEAGG
        + MERGETYPES_WEIGHTED_SIMPLEAGG
        + MERGETYPES_BYLENGTH
    )

    ########
    # repr #
    ########

    def __repr__(self):
        return f'<{self.__class__.__name__} object>\n{repr(self.df)}'

    ################
    # constructors #
    ################

    def __init__(self, refver):
        self.refver = refgenome.standardize(refver)
        self._df = None
        self._df_lastaccess = datetime.datetime.fromtimestamp(0)
        self._gr = None
        self._gr_lastaccess = datetime.datetime.fromtimestamp(0)
        self._sorted = False
        #self._columns = None

    def spawn(self, frame):
        result = self.__class__(self.refver)
        result.assign_frame(frame)
        return result

    @classmethod
    def from_frame(cls, frame, refver):
        result = cls(refver)
        result.assign_frame(frame)
        return result

    @classmethod
    def from_data(cls, refver, **kwargs):
        # parse arguments
        chrom_key = cls.from_data_argparse(kwargs, cls.chrom_arg_keys)
        start_key = cls.from_data_argparse(kwargs, cls.start_arg_keys)
        end_key = cls.from_data_argparse(kwargs, cls.end_arg_keys)
        if chrom_key is None:
            raise Exception(f'Chromosomes must be specified by one of these keywords (case-insensitive): {cls.chrom_arg_keys}')
        if (start_key is None) and (end_key is not None):
            raise Exception(f'Ends must not be given without starts.')

        # classify keys
        #annot_keys = set(kwargs.keys()).difference({chrom_key, start_key, end_key})
        annot_keys = [x for x in kwargs.keys() if x not in (chrom_key, start_key, end_key)]

        # make source dict
        source_dict = dict()
        source_dict['Chromosome'] = np.atleast_1d(kwargs[chrom_key])

        if start_key is not None:
            source_dict['Start'] = np.atleast_1d(kwargs[start_key])
        else:
            source_dict['Start'] = np.atleast_1d(0)

        if end_key is not None:
            source_dict['End'] = np.atleast_1d(kwargs[end_key])
        else:
            chromdict = refgenome.get_chromdict(refver)
            source_dict['End'] = [chromdict[x] for x in source_dict['Chromosome']]

        for key in annot_keys:
            source_dict[key] = np.atleast_1d(kwargs[key])

        keys, old_values = zip(*list(source_dict.items()))
        new_values = np.broadcast_arrays(*old_values)

        df = pd.DataFrame(dict(zip(keys, new_values)))

        return cls.from_frame(df, refver)

    @classmethod
    def from_records(cls, records, refver, annot_cols=()):
        columns = COMMON_COLUMNS + list(annot_cols)
        kwargs = dict()
        for key, val in itertools.zip_longest(columns, zip(*records), fillvalue=None):
            if (key is None) or (val is None):
                raise Exception(f'Length of annotation columns does not coincide between "records" and "annot_cols"')
            kwargs[key] = val
        return cls.from_data(refver=refver, **kwargs)

    @classmethod
    def all_regions(cls, refver, assembled_only=True):
        chromdict = refgenome.get_chromdict(refver)
        if assembled_only:
            chroms = refgenome.choose_assembled_chroms(chromdict.contigs)
        else:
            chroms = chromdict.contigs
        starts = 0
        ends = tuple(chromdict[x] for x in chroms)
        result = cls.from_data(
            refver=refver,
            chroms=chroms,
            starts=starts,
            ends=ends,
        )
        return result

    @staticmethod
    def from_data_argparse(kwargs, valid_keys):
        selected_keys = [x for x in kwargs.keys() if x.lower() in valid_keys]
        if len(selected_keys) == 0:
            return None
        elif len(selected_keys) == 1:
            return selected_keys[0]
        else:
            raise Exception(f'Among the argument keys, exactly one of these (case-insensitive) must exist: {valid_keys}')

    ###################
    # read/write file #
    ###################

    @classmethod
    def read_tsv(cls, filename, refver, annot_cols=None):
        frame = cls.make_frame_from_tsv(filename, annot_cols=annot_cols)
        return cls.from_frame(frame, refver)

    from_tsv = read_tsv

    def write_tsv(self, filename):
        assert re.search(r'(bed|tsv)(\.gz)?$', filename)
        self.df.to_csv(filename, sep='\t', header=True, index=False)

    to_tsv = write_tsv

    @staticmethod
    def make_frame_from_tsv(filename, annot_cols=None):
        with tools.openfile(filename, 'rt') as infile:
            firstline = next(infile).strip()
        firstline_sp = firstline.split('\t')
        assert len(firstline_sp) >= 3

        has_header = not (firstline_sp[1].isdecimal() and firstline_sp[2].isdecimal())
        if has_header:
            assert (
                (firstline_sp[0] == 'Chromosome')
                and (firstline_sp[1] == 'Start')
                and (firstline_sp[2] == 'End')
            )

        if (not has_header) and (annot_cols is None):
            result = pr.read_bed(filename)
        else:
            if has_header:
                if annot_cols is None:
                    skiprows = None
                    header = 0 
                    names = None
                else:
                    skiprows = 1
                    header = None
                    names = COMMON_COLUMNS + list(annot_cols)
            else:
                skiprows = 0
                header = None
                names = COMMON_COLUMNS + list(annot_cols)

            result = pd.read_csv(
                filename, 
                sep='\t', 
                skiprows=skiprows,
                header=header,
                names=names,
                dtype={'Chromosome': 'string', 'Start': int, 'End': int},
            )

        return result

    ###########################################################
    # assignment of dataframe/pyranges objects data populator #
    ###########################################################

    def assign_df(self, df):
        self.sanitycheck_df(df)
        self._df = self.postprocess_df(df)
        self.update_df_lastaccess()

    def assign_gr(self, gr):
        self.sanitycheck_gr(gr)
        self._gr = self.postprocess_gr(gr)
        self.update_gr_lastaccess()

    def assign_frame(self, frame):
        if isinstance(frame, pd.DataFrame):
            return self.assign_df(frame)
        elif isinstance(frame, pr.PyRanges):
            return self.assign_gr(frame)
        else:
            raise Exception(
                f'Input "frame" object must be either "pandas DataFrame" or "PyRanges" object.'
            )

    @staticmethod
    def sanitycheck_df(df):
        assert isinstance(df, pd.DataFrame)
        if not set(COMMON_COLUMNS).issubset(df.columns):
            raise Exception(f'Required columns: {COMMON_COLUMNS}')
        assert ((df['End'] - df['Start']) > 0).all()

    @staticmethod
    def postprocess_df(df):
        # dtype setting
        df = df.astype({'Chromosome': 'string', 'Start': int, 'End': int})

        # column reordering
        new_cols = COMMON_COLUMNS + [x for x in df.columns if x not in COMMON_COLUMNS] 
            # df.columns.difference(COMMON_COLUMNS).to_list() -> does not preserve ordering
        df = df.loc[:, new_cols]

        # remove index
        df = df.reset_index(drop=True, inplace=False)

        return df

    @staticmethod
    def sanitycheck_gr(gr):
        assert isinstance(gr, pr.PyRanges)

    @staticmethod
    def postprocess_gr(gr):
        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                if gr.dtypes['Start'] == np.int32:
                    gr.Start = gr.Start.to_numpy().astype('int64')
                if gr.dtypes['End'] == np.int32:
                    gr.End = gr.End.to_numpy().astype('int64')

        return gr

    ################################
    # df, gr accessors and setters #
    ################################

    @property
    def df(self):
        if self._df is None:
            self.set_df_from_gr()
        else:
            if self._gr is not None:
                if self.check_gr_is_newer():
                    if not self.check_df_gr_same():
                        self.set_df_from_gr()

        self.update_df_lastaccess()
        return self._df

    @df.setter
    def df(self, value):
        raise Exception(f'Direct assignment of "df" attribute is not allowed.')

    @property
    def gr(self):
        if self._gr is None:
            self.set_gr_from_df()
        else:
            if self._df is not None:
                if not self.check_gr_is_newer():  # df is newer
                    if not self.check_df_gr_same():
                        self.set_gr_from_df()

        self.update_gr_lastaccess()
        return self._gr

    @gr.setter
    def gr(self, value):
        raise Exception(f'Direct assignment of "gr" attribute is not allowed.')

    #####################
    # helpers of df, gr #
    #####################

    def set_gr_from_df(self):
        assert self._df is not None
        self._gr = pr.PyRanges(self._df, int64=True)

    def set_df_from_gr(self):
        assert self._gr is not None
        self._df = self._gr.df.reset_index(drop=True)
            # "df" method of PyRanges object returns a newly created DataFrame
            # (gr.df is gr.df) is False

    def sync_df_gr(self):
        if not self.check_df_gr_same():
            if self.check_gr_is_newer():
                self.set_df_from_gr()
            else:
                self.set_gr_from_df()

    def update_df(self):
        if self.check_gr_is_newer():
            if not self.check_df_gr_same():
                self.set_df_from_gr()

    def update_gr(self):
        if not self.check_gr_is_newer():
            if not self.check_df_gr_same():
                self.set_gr_from_df()

    def update_df_lastaccess(self):
        self._df_lastaccess = datetime.datetime.now()

    def update_gr_lastaccess(self):
        self._gr_lastaccess = datetime.datetime.now()

    def check_gr_is_newer(self):
        return self._gr_lastaccess > self._df_lastaccess

    def check_df_gr_same(self):
        """This does not consider the order of rows because hash is calculated after sorting with genome coordinates."""
        assert self._df is not None
        assert self._gr is not None

        df_for_hash = self._df.loc[:, self._gr.columns]  # make column orders the same
        return (hash_gr(self._gr) == hash_df(df_for_hash)).to_numpy().all()  # It is faster (~2.5 us vs ~13 us) to convert into numpy and then do .all()


    ###################################
    # adaptations of pyranges methods #
    ###################################

    EXCLUDED_PYRANGES_METHODS = [
        'lengths',
    ]

    def _pyranges_method_adapter_alone(methodname):
        def new_method(self, *args, **kwargs):
            result_gr = getattr(self.gr, methodname)(*args, **kwargs)
            return self.spawn(result_gr)

        new_method.__name__ = methodname
        return new_method

    def _pyranges_method_adapter_other(methodname):
        def new_method(self, other, *args, **kwargs):
            if isinstance(other, pr.PyRanges):
                other_gr = other
            elif isinstance(other, GenomeDataFrame):
                other_gr = other.gr
            else:
                raise Exception(
                    f'"other" must be either PyRanges or GenomeDataFrame object. The type of "other" is: {type(other)}'
                )

            result_gr = getattr(self.gr, methodname)(other_gr, *args, **kwargs)
            return self.spawn(result_gr)

        new_method.__name__ = methodname
        return new_method

    for key in set(dir(pr.PyRanges)).difference(EXCLUDED_PYRANGES_METHODS):
        if key.startswith('_'):
            continue
        mthd = getattr(pr.PyRanges, key)
        if not callable(mthd):
            continue
        
        params = inspect.signature(mthd).parameters
        if (
            (len(params) >= 2)
            and (list(params.keys())[1] == 'other')
        ):
            mthd = _pyranges_method_adapter_other(key)
            exec(f'{key} = mthd')
        else:
            mthd = _pyranges_method_adapter_alone(key)
            exec(f'{key} = mthd')

    ###########################################
    # adaptations of pandas DataFrame methods #
    ###########################################

    def __getitem__(self, key):
        return self.df.__getitem__(key)

    def __setitem__(self, key, val):
        if self._df is None:
            self.set_df_from_gr()
            
        self._df.__setitem__(key, val)
        self.assign_df(self._df)

    @property
    def loc(self):
        return GenomeDataFrameLoc(self)

    @property
    def iloc(self):
        return GenomeDataFrameILoc(self)

    @staticmethod
    def dataframe_method_adaptor(methodname):
        original_method = getattr(pd.DataFrame, methodname)
        assert callable(original_method)

        sig = inspect.signature(original_method)
        has_inplace = ('inplace' in sig.parameters.keys())

        if has_inplace:
            def new_method(self, *args, **kwargs):
                ba = sig.bind(self, *args, **kwargs)
                ba.apply_defaults()
                inplace = ba.arguments['inplace']

                self.update_df()
                if inplace:
                    getattr(self._df, methodname)(*args, **kwargs)
                    self.assign_df(self._df)
                    self.update_df_lastaccess()
                else:
                    new_df = getattr(self.df, methodname)(*args, **kwargs)
                    return self.spawn(new_df)

            new_method.__name__ = methodname
            return new_method
        else:
            # assume the original method returns a DataFrame
            def new_method(self, *args, **kwargs):
                new_df = getattr(self.df, methodname)(*args, **kwargs)
                return self.spawn(new_df)

        new_method.__name__ = methodname
        return new_method

    for mthdname in (
        'rename',
    ):
        new_mthd = dataframe_method_adaptor(mthdname)
        exec(f'{mthdname} = new_mthd')

    #################################
    # custom methods and properties #
    #################################

    def get_newer_frame(self):
        if self.check_gr_is_newer():
            assert self._gr is not None
            return self._gr
        else:
            assert self._df is not None
            return self._df

    @property
    def is_empty(self):
        newer_frame = self.get_newer_frame()
        if isinstance(newer_frame, pd.DataFrame):
            return newer_frame.shape[0] == 0
        else:
            return newer_frame.empty

    @property
    def columns(self):
        if self.is_empty:
            return list(COMMON_COLUMNS)  # Chromosome, Start, End
        else:
            return self.get_newer_frame().columns.to_list()

    @property
    def shape(self):
        return self.df.shape

    @property
    def lengths(self):
        return (self['End'] - self['Start']).to_numpy()

    @property
    def annot_cols(self):
        return self.df.columns.drop(COMMON_COLUMNS).to_list()

    @property
    def chromdict(self):
        return refgenome.get_chromdict(self.refver)

    def sort(self):
        """In-place sort, only applied to df, not gr"""
        if self._df is None:
            self.set_df_from_gr()
            do_sort = True
        elif self._gr is None:
            do_sort = (not self._sorted)
        else:
            if self.check_gr_is_newer() and (not self.check_df_gr_same()):
                self.set_df_from_gr()
                do_sort = True
            else:
                do_sort = (not self._sorted)

        if do_sort:
            self_df = self._df
            chrom_indexes = self_df['Chromosome'].apply(self.chromdict.contigs.index)
            sortkey = np.lexsort([self_df['End'], self_df['Start'], chrom_indexes])
            self.assign_df(self_df.iloc[sortkey, :])
            self._sorted = True

    def group_bychrom(self, how=None):
        assert how in (None, 'gr', 'df')

        self_gr = self.gr
        all_chroms = set(self_gr.Chromosome)
        result_grs = {chrom: self_gr[chrom] for chrom in all_chroms}

        if how == 'gr':
            return result_grs
        else:
            if how == 'df':
                return {key: val.df for key, val in result_grs.items()}
            elif how is None:
                return {key: self.spawn(val) for key, val in result_grs.items()}

    def isec_union(self, other):
        self_gr = self.gr
        other_gr = other.gr

        isec = self_gr.intersect(other_gr)
        self_diff_isec = self_gr.subtract(isec)
        other_diff_isec = other_gr.subtract(isec)
        result_gr = pr.concat([isec, self_diff_isec, other_diff_isec])
        return self.spawn(result_gr)

    def subset(self, chrom, start0=None, end0=None):
        """Returns a subset of original rows which overlaps the interval specified by arguments.
        Original rows are returned as is, without intersection with the input interval.
        'start0' and 'end0' arguments are 0-based & half-open.
        """
        getitem_key = (chrom, slice(start0, end0))
        return self.spawn(self.gr[getitem_key])

    def subset_chroms(self, chromlist):
        df = self.df
        selector = df['Chromosome'].isin(set(chromlist))
        return self.spawn(df.loc[selector, :])

    def get_midpoints(self):
        return ((self['Start'] + self['End'] - 1) / 2).to_numpy().astype(int)

    def drop_annots(self):
        return self.spawn(self.df.iloc[:, :3])

    def copy(self):
        return self.spawn(self.df.copy())

    ################
    # segmentation #
    ################

    def make_compact(self):
        gdfs_bychrom = self.group_bychrom(how=None)
        #encoders = dict()
        decoders = dict()
        compact_gdfs_bychrom = dict()
        for chrom, subgdf in gdfs_bychrom.items():
            compact_pos0_list, decoder = compact_coords(subgdf['Start'], subgdf['End'])
            #encoders[chrom] = enc
            decoders[chrom] = decoder
            compact_gdf = subgdf.copy()
            compact_gdf['Start'] = compact_pos0_list
            compact_gdf['End'] = compact_gdf['Start'] + 1
            compact_gdfs_bychrom[chrom] = compact_gdf

        concat_gdf = self.spawn(pd.concat((x.df for x in compact_gdfs_bychrom.values()), axis=0))
        concat_gdf.sort()

        return concat_gdf, compact_gdfs_bychrom, decoders

    def get_segment(
        self, 
        annot_colname,
        N_colname='N',
        smoothing=False, 
        verbose=True, 
        smooth_kwargs=dict(), 
        segment_kwargs=dict(),
    ):
        concat_gdf, compact_gdfs_bychrom, decoders = self.make_compact()
        compact_seg_df = rdnacopy.run_dnacopy_segment(
            chromosomes=concat_gdf['Chromosome'],
            positions=concat_gdf['Start'],
            values=concat_gdf[annot_colname],

            N_colname='N',
            value_colname='value',
            smoothing=smoothing, 
            verbose=verbose, 
            smooth_kwargs=smooth_kwargs, 
            segment_kwargs=segment_kwargs,
        )
        compact_seg_gdf = self.spawn(compact_seg_df)
        bychrom = compact_seg_gdf.group_bychrom(how=None)

        chrom_list = list()
        start0_list = list()
        end0_list = list()
        N_list = list()
        value_list = list()

        for chrom, sub_gdf in bychrom.items():
            chrom_list.extend(sub_gdf['Chromosome'])
            N_list.extend(sub_gdf['N'])
            value_list.extend(sub_gdf['value'])

            noncompact_start0s, noncompact_end0s = uncompact_coords(sub_gdf['Start'], sub_gdf['End'], decoders[chrom])
            start0_list.extend(noncompact_start0s)
            end0_list.extend(noncompact_end0s)
            
        from_data_kwargs = {
            'chrom': chrom_list,
            'start': start0_list,
            'end': end0_list,
            annot_colname: value_list,
            N_colname: N_list,
        }
        result = self.__class__.from_data(self.refver, **from_data_kwargs)
        result.sort()

        return result


    ########
    # join #
    ########

    @deco.get_deco_arg_choices({'how': [None, 'left']})
    def join(
        self, other, 
        other_cols=None,
        how='left', 
        merge=['mean', 'std'],
        #add_std=False, 
        ddof=0,
        overlapping_length=False,
        N_colname='N',
        #std_suffix='_std',
        suffixes=None,
    ):
        other_cols, suffixes, merge = self.join_sanitycheck_arghandler(self, other, other_cols, merge, N_colname, suffixes)

        # make parameters
        left_columns = self.columns
        report_overlap = (
            overlapping_length 
            and set(merge).intersection(self.MERGETYPES_BYLENGTH)
        )

        # run pyranges join
        joined_df = self.gr.join(other.gr[other_cols], how=how, report_overlap=report_overlap).df
        joined_df.reset_index(drop=True, inplace=True)

        # fill missing slots with pd.NA
        unmatched_selector = (joined_df['Start_b'] == -1).to_numpy()
        #joined_df.loc[unmatched_selector, (['Start_b', 'End_b'] + other_cols)] = pd.NA  # this becomes np.nan in float columns
        joined_df.loc[unmatched_selector, ~joined_df.columns.isin(left_columns)] = pd.NA  # this becomes np.nan in float columns

        # merge
        if merge is None:
            result_df = joined_df
        else:
            unmatched_df = joined_df.loc[unmatched_selector, left_columns]

            matched_left = joined_df.loc[~unmatched_selector, left_columns]
            counts, groupkey = get_coord_groupkey(matched_left, self.chromdict)
                # "counts" is aligned with "groupkey"

            merged_left = matched_left.drop_duplicates()
            merged_left.reset_index(drop=True, inplace=True)
            merged_left.insert(merged_left.shape[1], N_colname, counts)

            matched_right = joined_df.loc[~unmatched_selector, ~joined_df.columns.isin(left_columns)]
            matched_right.reset_index(drop=True, inplace=True)
            merged_right_list = [
                self.join_merge_right(matched_right, groupkey, other_cols, mergetype, suffixes, ddof)
                for mergetype in merge
            ]
            merged_right = pd.concat(merged_right_list, axis=1)

            #matched_df = joined_df.loc[~unmatched_selector, :]

            #merged_df = self.join_merge(matched_df, self.chromdict, left_columns, other_cols, merge, add_std, ddof, N_colname, suffix, std_suffix)
            #result_df = pd.concat([unmatched_df, merged_df], axis=0, ignore_index=True)

            result_matched = pd.concat([merged_left, merged_right], axis=1)
            result_df = pd.concat([unmatched_df, result_matched], axis=0)

        # return
        result = self.spawn(result_df)
        result.sort()
        return result

    @classmethod
    def join_sanitycheck_arghandler(cls, self, other, other_cols, merge, N_colname, suffixes):
        # sanitycheck BEFORE args modification #
        assert len(self.columns) == len(set(self.columns))
        assert len(other.columns) == len(set(other.columns))

        assert not check_duplicate_coords(self.df)
        assert not check_duplicate_coords(other.df)
    
        # args modification #
        if other_cols is None:
            other_cols = other.annot_cols

        if merge is not None:
            merge = np.atleast_1d(merge)

        if merge is not None:
            if suffixes is None:
                suffixes = {
                    mergetype: f'_{mergetype}'
                    for mergetype in merge
                }

        # sanitycheck AFTER args modification #

        if merge is not None:
            assert len(merge) > 0
            assert set(merge).issubset(cls.AVAILABLE_MERGE_METHODS)

            assert isinstance(suffixes, dict)
            assert set(suffixes.keys()) == set(merge)
            assert len(tuple(suffixes.values())) == len(set(suffixes.values())), f'Duplicate values in "suffixes" dictionary'

        assert len(other_cols) == len(set(other_cols))
        assert set(other_cols).issubset(other.annot_cols)
        assert not set(self.annot_cols).intersection(other_cols), (
            f'There must be no overlapping annotation columns between two {self.__class__.__name__} objects.'
        )
        assert not {'Start_b', 'End_b'}.intersection(self.annot_cols + other_cols)

        if merge is None:
            cols_tobe_added = other_cols
        else:
            cols_tobe_added = list(itertools.product(other_cols, suffixes.values()))
        assert not set(self.annot_cols).intersection(cols_tobe_added)

        final_cols = self.columns + cols_tobe_added
        assert N_colname not in final_cols

        return other_cols, suffixes, merge

    def join_merge_right(cls, matched_right, groupkey, other_cols, mergetype, suffixes, ddof):
        if mergetype in cls.MERGETYPES_SIMPLEAGG:
            merged_right = cls.join_merge_right_simpleagg(matched_right, groupkey, other_cols, mergetype, ddof)
        elif mergetype in cls.MERGETYPES_WEIGHTED_SIMPLEAGG:
            merged_right = cls.join_merge_right_weighted_simpleagg(matched_right, groupkey, other_cols, mergetype, ddof)
        elif mergetype in cls.MERGETYPES_BYLENGTH:
            merged_right = cls.join_merge_right_bylength(matched_right, groupkey, other_cols, mergetype, ddof)

        merged_right.reset_index(drop=True, inplace=True)
        suffix = suffixes[mergetype]
        merged_right.rename(columns=(lambda x: x + suffix), inplace=True)
        return merged_right

    @classmethod
    def join_merge_right_simpleagg(cls, matched_right, groupkey, other_cols, mergetype, ddof):
        assert mergetype in cls.MERGETYPES_SIMPLEAGG

        groupby = matched_right.groupby(groupkey)[other_cols]

        if mergetype == 'first':
            merged_right = groupby.first()
        elif mergetype == 'last':
            merged_right = groupby.last()
        elif mergetype == 'mean':
            merged_right = groupby.mean()
        elif mergetype == 'std':
            merged_right = groupby.std(ddof=ddof)
        elif mergetype == 'max':
            merged_right = groupby.max()
        elif mergetype == 'min':
            merged_right = groupby.min()
        elif mergetype == 'median':
            merged_right = groupby.median()
        elif mergetype == 'skew':
            merged_right = groupby.skew()
        else:
            raise Exception(f'Unknown mergetype: {mergetype}')

        return merged_right

    @classmethod
    def join_merge_right_weighted_simpleagg(cls, matched_right, groupkey, other_cols, mergetype, ddof):
        assert mergetype in cls.MERGETYPES_WEIGHTED_SIMPLEAGG
        #assert (matched_right.index == pd.RangeIndex(start=0, end=matched_right.shape[0])).all()

        new_matched_right = cls.join_make_weighted_right(matched_right, other_cols)
        groupby = new_matched_right.groupby(groupkey)[other_cols]

        if mergetype == 'weighted_mean':
            merged_right = groupby.sum()
        else:
            raise Exception(f'Unknown mergetype: {mergetype}')

        return merged_right

    @classmethod
    def join_make_weighted_right(cls, matched_right, other_cols):
        if 'Overlap' in matched_right.columns:
            lengths = matched_right['Overlap']
        else:
            lengths = matched_right['End_b'] - matched_right['Start_b']

        lengthsum = lengths.sum()

        src_data = {
            key: ((matched_right[key] * lengths) / lengthsum)
            for key in other_cols
        }
        weighted_matched_right = pd.DataFrame(src_data)

        return weighted_matched_right

    @classmethod
    def join_merge_right_bylength(cls, matched_right, groupkey, other_cols, mergetype, ddof):
        assert mergetype in cls.MERGETYPES_BYLENGTH
        #assert (matched_right.index == pd.RangeIndex(start=0, end=matched_right.shape[0])).all()

        if 'Overlap' in matched_right.columns:
            lengths = matched_right['Overlap']
        else:
            lengths = matched_right['End_b'] - matched_right['Start_b']

        groupby = lengths.reset_index(drop=True).groupby(groupkey)
        if mergetype == 'longest':
            indexes = groupby.idxmax()
        elif mergetype == 'shortest':
            indexes = groupby.idxmin()
        else:
            raise Exception(f'Unknown mergetype: {mergetype}')

        other_cols_idxs = [matched_right.columns.get_loc(x) for x in other_cols]
        merged_right = matched_right.iloc[indexes.to_numpy(), other_cols_idxs]

        return merged_right

    ###########
    # binning #
    ###########

    def get_binsize(self):
        self.sort()

        gdfs_bychrom = self.group_bychrom(how=None)
        equal_length_portions = list()
        last_ones = list()
        for chrom, subgdf in gdfs_bychrom.items():
            lengths = subgdf.lengths()
            last_ones.append(lengths[-1])
            if len(lengths) > 1:
                equal_length_portions.append(lengths[:-1])

        all_lengths = np.unique(np.concatenate(equal_length_portions))
        if all_lengths.shape != (1,):
            raise Exception(f'Bin lengths are not identical.')
        binsize = all_lengths[0]

        if len(last_ones) > 0:
            last_ones = np.asarray(last_ones)
            if not (last_ones <= binsize).all():
                raise Exception(f'Some last interval lengths are greater than the bin length.')

        return binsize

    @staticmethod
    def get_upsizing_groupkey(singlechrom_df, chromlen, new_binsize, old_binsize):
        nrow = singlechrom_df.shape[0]
        assert nrow == np.ceil(chromlen / old_binsize).astype(int)

        n_newbin = np.ceil(chromlen / new_binsize).astype(int)
        olddf_idx_ends = np.fromiter(
            (
                min(int((new_binsize * x) / old_binsize), nrow)
                for x in range(1, n_newbin + 1)
            ),
            dtype=int,
        )
        group_lengths = np.diff(olddf_idx_ends, prepend=0)
        groupkey = np.repeat(range(len(group_lengths)), group_lengths)
        return groupkey, n_newbin

    @staticmethod
    def upsize_singlechrom_df(
        singlechrom_df, 
        new_binsize, 
        groupkey, 
        n_newbin, 
        chromlen, 
        annotcols_numeric,
        annotcols_bool,
        merge_bool,
    ):
        assert merge_bool in ('any', 'all')
        assert len(annotcols_numeric) > 0

        starts = np.arange(n_newbin) * new_binsize
        ends = np.concatenate([starts[1:], [chromlen]])
        result_common = pd.DataFrame(
            {
                'Chromosome': singlechrom_df['Chromosome'][0],
                'Start': starts,
                'End': ends,
            }
        )

        groupby = singlechrom_df.groupby(by=groupkey, sort=False)

        result_annot_mean = groupby[annotcols_numeric].mean()
        result_annot_mean.reset_index(inplace=True, drop=True)

        if len(annotcols_bool) > 0:
            if merge_bool == 'any':
                result_annot_any = groupby[annotcols_bool].any()
            elif merge_bool == 'all':
                result_annot_any = groupby[annotcols_bool].all()
            result_annot_any.reset_index(inplace=True, drop=True)

        if len(annotcols_bool) > 0:
            result_df = pd.concat([result_common, result_annot_mean, result_annot_any], axis=1)
        else:
            result_df = pd.concat([result_common, result_annot_mean], axis=1)

        return result_df

    def upsize_bin(self, size, merge_bool='any'):
        old_binsize = self.get_binsize()
        if size <= old_binsize:
            raise Exception(f'New bin size ({size}) must be greater than the old bin size ({old_binsize})')

        dtypes = self.df.dtypes
        annotcols_numeric = list()
        annotcols_bool = list()
        for x in self.annot_cols:
            dtype = dtypes[x]
            if pd.api.types.is_bool_dtype(dtype):
                annotcols_bool.append(x)
            elif pd.api.types.is_numeric_dtype(dtype):
                annotcols_numeric.append(x)

        dfs_bychrom = self.group_bychrom(how='df')
        newbin_dfs_bychrom = dict()
        for chrom, singlechrom_df in dfs_bychrom.items():
            chromlen = self.chromdict[chrom]
            groupkey, n_newbin = self.get_upsizing_groupkey(
                singlechrom_df=singlechrom_df, 
                chromlen=chromlen, 
                new_binsize=size, 
                old_binsize=old_binsize,
            )
            newbin_dfs_bychrom[chrom] = self.upsize_singlechrom_df(
                singlechrom_df=singlechrom_df, 
                new_binsize=size, 
                groupkey=groupkey, 
                n_newbin=n_newbin, 
                chromlen=chromlen, 
                annotcols_numeric=annotcols_numeric,
                annotcols_bool=annotcols_bool,
                merge_bool=merge_bool,
            )

        result_df = pd.concat(newbin_dfs_bychrom.values(), axis=0)
        result = self.spawn(result_df)
        return result


class GenomeDataFrameLoc:
    def __init__(self, gdf):
        self.gdf = gdf

    def __getitem__(self, key):
        result = self.gdf.df.loc[key]
        if isinstance(result, pd.DataFrame):
            return self.gdf.spawn(result)
        else:
            return result

    def __setitem__(self, key, val):
        self.gdf.df.loc.__setitem__(key, val)


class GenomeDataFrameILoc:
    def __init__(self, gdf):
        self.gdf = gdf

    def __getitem__(self, key):
        result = self.gdf.df.iloc[key]
        if isinstance(result, pd.DataFrame):
            return self.gdf.spawn(result)
        else:
            return result

    def __setitem__(self, key, val):
        self.gdf.df.iloc.__setitem__(key, val)



##########################
# module level functions #
##########################

def concat(gdf_list):
    refvers = set(x.refver for x in gdf_list)
    if len(refvers) != 1:
        raise Exception(f'Reference versions differ')
    refver = refvers.pop()

    new_df = pd.concat([x.df for x in gdf_list], axis=0)
    return GenomeDataFrame.from_frame(new_df, refver)
    

def hash_gr(gr):
    #df_for_hash = pd.concat([getattr(gr, x) for x in gr.columns], axis=1)
    df_for_hash = gr.df
    df_for_hash.sort_values(COMMON_COLUMNS, inplace=True)
    df_for_hash.reset_index(drop=True, inplace=True)
    return hash_pandas_object(df_for_hash, index=False)


def hash_df(df):
    assert set(COMMON_COLUMNS).issubset(df.columns)
    df_for_hash = df.sort_values(COMMON_COLUMNS, inplace=False)
    df_for_hash.reset_index(drop=True, inplace=True)
    return hash_pandas_object(df_for_hash, index=False)


def get_coord_groupkey(df, chromdict):
    """Does not sort before grouping, like itertools.groupby"""
    chrom_indexes = df['Chromosome'].apply(chromdict.contigs.index)
    coord_arr = np.stack(
        [chrom_indexes.to_numpy(), df['Start'].to_numpy(), df['End'].to_numpy()], 
        axis=1,
    )
    _, counts, groupkey = tools.array_grouper(coord_arr, omit_values=True)
    return counts, groupkey


def check_duplicate_coords(df):
    return df.loc[:, ['Chromosome', 'Start', 'End']].duplicated().any()
    

def check_interval_overlap(start0s, end0s):
    start0s = np.asarray(start0s)
    end0s = np.asarray(end0s)
    assert (start0s.ndim == 1) and (end0s.ndim == 1)

    idxs = np.argsort(start0s)
    start0s_sort = start0s[idxs]
    end0s_sort = end0s[idxs]

    return not (start0s_sort[1:] >= end0s_sort[:-1]).all()


def compact_coords(start0s, end0s):
    # sanity check
    if check_interval_overlap(start0s, end0s):
        raise Exception(f'Intervals must not overlap.')

    intervals = list(zip(start0s, end0s))
    #encoder = dict()  # noncompact (start, end) -> compact single pos
    #decoder = dict()  # compact single pos -> noncompact (start, end)
    compact_pos0_list = np.argsort(start0s)
    decoder = dict(zip(compact_pos0_list, intervals))
    #for compact_pos0, intv in zip(compact_pos0_list, intervals):
        #encoder[intv] = compact_pos0
    #    decoder[compact_pos0] = intv
        #compact_pos0_list.append(compact_pos0)

    #return encoder, decoder
    return compact_pos0_list, decoder


def uncompact_coords(compact_start0s, compact_end0s, decoder):
    start0s = list()
    end0s = list()
    for compact_s0, compact_e0 in zip(compact_start0s, compact_end0s):
        s0 = decoder[compact_s0][0]
        start0s.append(s0)

        last_pos0 = compact_e0 - 1
        e0 = decoder[last_pos0][1]
        end0s.append(e0)
        
    return start0s, end0s



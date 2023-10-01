import re
import os
import datetime
import inspect
import itertools
import contextlib
import multiprocessing
import functools

import pandas as pd
import numpy as np
import pyranges as pr

import handygenome.deco as deco
import handygenome.refgenome.refgenome as refgenome
import handygenome.tools as tools
import handygenome.logutils as logutils
import handygenome.cnv.rdnacopy as rdnacopy

import handygenome.genomedf.genomedf_utils as genomedf_utils
import handygenome.genomedf.genomedf_methods as genomedf_methods


##############################
# coordinate argument parser #
##############################

def coord_argparse_OLD(**kwargs):
    all_valid_keys = {
        'chrom': ['chrom', 'chroms', 'chromosome', 'chromosomes'],
        'start': ['start', 'starts', 'start0', 'start0s'],
        'end': ['end', 'ends', 'end0', 'end0s'],
    }

    # parse coordinates from raw arguments
    parsed_keys = dict()
    parsed_values = dict()
    for coordtype, valid_keys in all_valid_keys.items():
        selected_keys = [x for x in kwargs.keys() if x.lower() in valid_keys]
        if len(selected_keys) == 0:
            parsed_keys[coordtype] = None
            parsed_values[coordtype] = None
        elif len(selected_keys) == 1:
            parsed_keys[coordtype] = selected_keys[0]
            parsed_values[coordtype] = kwargs[selected_keys[0]]
            if parsed_values[coordtype] is not None:
                parsed_values[coordtype] = np.atleast_1d(parsed_values[coordtype])
        else:
            raise Exception(
                f'Among the argument keys, exactly one of '
                f'these (case-insensitive) must exist: {valid_keys}'
            )

    # sanitycheck
    if parsed_keys['chroms'] is None:
        raise Exception(
            f'Chromosomes must be specified by one of these '
            f'keywords (case-insensitive): {all_valid_keys["chroms"]}'
        )
    if (parsed_keys['start0s'] is None) and (parsed_keys['end0s'] is not None):
        raise Exception(f'Ends must not be given without starts.')

    # broadcast only not-None values
    not_none_coordtypes, not_none_values = zip(
        *(
            (coordtype, val) for (coordtype, val) in parsed_values.items() 
            if val is not None
        )
    )
    try:
        new_values = np.broadcast_arrays(*not_none_values)
    except ValueError as exc:
        if str(exc).startswith('shape mismatch'):
            myexc = Exception(
                f'"chroms", "start0s", "end0s" must be broadcastable to a single shape.'
            )
            raise myexc from exc
        else:
            raise

    for coordtype, new_val in zip(not_none_coordtypes, new_values):
        parsed_values[coordtype] = new_val

    # make chroms into 'object' type array (by default fixed length unicode)
    parsed_values['chroms'] = parsed_values['chroms'].astype(object)

    # prepare results
    coords = {
        'Chromosome': parsed_values['chroms'],
        'Start': parsed_values['start0s'],
        'End': parsed_values['end0s'],
    }
    coord_keys = set(parsed_keys.keys()).difference({None})
    other_kwargs = dict(
        (k, v) for (k, v) in kwargs.items()
        if k not in coord_keys
    )

    return coords, other_kwargs


COORDARG_MANDATORY_KEYS = ['refver', 'chroms', 'start0s', 'end0s']


def modify_coord_args(kwargs):
    """Modifies kwargs in-place. This is an intended behavior in order to modify "arguments" attribute of BoundArguments object."""
    if not set(COORDARG_MANDATORY_KEYS).issubset(kwargs.keys()):
        raise Exception(f'Mandatory keys: {COORDARG_MANDATORY_KEYS}')

    # sanitycheck
    if kwargs['chroms'] is None:
        if (
            (kwargs['start0s'] is not None) 
            or (kwargs['end0s'] is not None)
        ):
            raise Exception(f'Start and End must not be specified without Chromosomes')
    if (kwargs['start0s'] is None) and (kwargs['end0s'] is not None):
        raise Exception(f'Ends must not be given without starts.')

    # return without any change when chroms is not given
    if kwargs['chroms'] is None:
        return

    # transform chrom values
    kwargs['chroms'] = np.atleast_1d(kwargs['chroms']).astype(object)

    # fill in default start0s and end0s
    if kwargs['start0s'] is None:
        kwargs['start0s'] = 0
    if kwargs['end0s'] is None:
        chromdict = refgenome.get_chromdict(kwargs['refver'])
        kwargs['end0s'] = np.asarray([chromdict[x] for x in kwargs['chroms']])

    # broadcast
    try:
        kwargs['chroms'], kwargs['start0s'], kwargs['end0s'] = np.broadcast_arrays(
            kwargs['chroms'], kwargs['start0s'], kwargs['end0s'],
        )
    except ValueError as exc:
        if str(exc).startswith('shape mismatch'):
            myexc = Exception(
                f'"chroms", "start0s", "end0s" must be broadcastable to a single shape.'
            )
            raise myexc from exc
        else:
            raise


def deco_coordarg(func):
    sig = inspect.signature(func)
    if not set(COORDARG_MANDATORY_KEYS).issubset(sig.parameters.keys()):
        raise Exception(
            f'Parameters of the fuction {func.__name__} must include all of these: {COORDARG_MANDATORY_KEYS}'
        )

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        ba = sig.bind(*args, **kwargs)
        ba.apply_defaults()
        modify_coord_args(ba.arguments)
        return func(*ba.args, **ba.kwargs)

    return wrapper


##############
# exceptions #
##############

class NotBinnedError(Exception):
    pass


#########################
# GenomeDataFrame class #
#########################

class GenomeDataFrame:
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

    def __del__(self):
        del self._df
        del self._gr

    def spawn(self, frame, **kwargs):
        result = self.__class__(self.refver, **kwargs)
        result.assign_frame(frame, self.dtypes.to_dict())
        return result

    def spawn_from_data(self, **kwargs):
        kwargs['refver'] = self.refver
        return self.__class__.from_data(**kwargs)

    @classmethod
    def init_empty(cls, refver, annot_cols=tuple()):
        columns = genomedf_utils.COMMON_COLUMNS + list(annot_cols)
        df = pd.DataFrame([], columns=columns)
        return cls.from_frame(df, refver)

    def spawn_empty(self, annot_cols=tuple()):
        columns = genomedf_utils.COMMON_COLUMNS + list(annot_cols)
        df = pd.DataFrame([], columns=columns)
        return self.spawn(df)

    @classmethod
    def concat(cls, gdf_iterator):
        gdf_list = list(gdf_iterator)
        assert len(gdf_list) > 0

        refvers = set(x.refver for x in gdf_list)
        if len(refvers) != 1:
            raise Exception(f'Reference versions differ: {refvers}')
        refver = refvers.pop()

        new_df = pd.concat([x.df for x in gdf_list], axis=0)
        return cls.from_frame(new_df, refver)

    @classmethod
    def from_frame(cls, frame, refver, dtype=dict()):
        result = cls(refver)
        result.assign_frame(frame, dtype)
        return result

    @classmethod
    @deco_coordarg
    def from_data(cls, refver, chroms, start0s=None, end0s=None, **kwargs):
        # make source dict
        source_dict = dict()
        source_dict['Chromosome'] = chroms

        if start0s is None:
            source_dict['Start'] = np.atleast_1d(0)
        else:
            source_dict['Start'] = start0s

        if end0s is None:
            chromdict = refgenome.get_chromdict(refver)
            source_dict['End'] = [chromdict[x] for x in source_dict['Chromosome']]
        else:
            source_dict['End'] = end0s

        for key, val in kwargs.items():
            source_dict[key] = val

        df = pd.DataFrame(source_dict)

        return cls.from_frame(df, refver)

    @classmethod
    def from_records(cls, records, refver, annot_cols=tuple()):
        columns = genomedf_utils.COMMON_COLUMNS + list(annot_cols)
        df_source = dict()
        for key, val in itertools.zip_longest(columns, zip(*records), fillvalue=None):
            if (key is None) or (val is None):
                raise Exception(f'Length of annotation columns does not coincide between "records" and "annot_cols"')
            df_source[key] = val
        return cls.from_frame(pd.DataFrame(df_source), refver)

    @classmethod
    def all_regions(cls, refver, assembled_only=True):
        chromdict = refgenome.get_chromdict(refver)
        if assembled_only:
            chroms = refgenome.choose_assembled_chroms(chromdict.contigs)
        else:
            chroms = chromdict.contigs
        result = cls.from_data(refver=refver, chroms=chroms)
        return result

    @classmethod
    def from_margins(cls, refver, chrom_left, start0_left, chrom_right, end0_right):
        chromdict = refgenome.get_chromdict(refver)
        cumpos0_left = chromdict.get_cumpos0(chrom_left, start0_left)
        cumpos0_right = chromdict.get_cumpos0(chrom_right, end0_right)
        if cumpos0_left >= cumpos0_right:
            raise Exception(
                f'"left position" comes later than "right_position"; '
                f'chrom_left={chrom_left}, start0_left={start0_left}, '
                f'chrom_right={chrom_right}, end0_right={end0_right}'
            )

        #result = cls(refver)
        chroms = list()
        start0s = list()
        end0s = list()

        if chrom_left == chrom_right:
            chroms.append(chrom_left)
            start0s.append(start0_left)
            end0s.append(end0_right)
        else:
            chrom_left_idx = chromdict.contigs.index(chrom_left)
            chrom_right_idx = chromdict.contigs.index(chrom_right)

            chroms.append(chrom_left)
            start0s.append(start0_left)
            end0s.append(chromdict[chrom_left])

            for chrom in chromdict.contigs[(chrom_left_idx + 1):chrom_right_idx]:
                chroms.append(chrom)
                start0s.append(0)
                end0s.append(chromdict[chrom])

            chroms.append(chrom_right)
            start0s.append(0)
            end0s.append(end0_right)

        result = cls.from_data(
            refver=refver, chroms=chroms, start0s=start0s, end0s=end0s,
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
    def read_tsv(cls, filename, refver, annot_cols=None, dtype=dict()):
        frame = cls.make_frame_from_tsv(
            filename, annot_cols=annot_cols, dtype=dtype,
        )
        return cls.from_frame(frame, refver)

    from_tsv = read_tsv

    def write_tsv(self, filename):
        patstring = r'(bed|tsv)(\.gz)?$'
        assert re.search(patstring, filename), f'File name must match {patstring}'
        self.df.to_csv(filename, sep='\t', header=True, index=False)

    to_tsv = write_tsv

    @classmethod
    def make_frame_from_tsv(cls, filename, annot_cols=None, dtype=dict()):
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
                    names = genomedf_utils.COMMON_COLUMNS + list(annot_cols)
            else:
                skiprows = 0
                header = None
                names = genomedf_utils.COMMON_COLUMNS + list(annot_cols)

            dtype = (
                {'Chromosome': 'string', 'Start': int, 'End': int}
                | dtype
            )
            result = pd.read_csv(
                filename, 
                sep='\t', 
                skiprows=skiprows,
                header=header,
                names=names,
                dtype=dtype,
            )

        return result

    ###########################################################
    # assignment of dataframe/pyranges objects data populator #
    ###########################################################

    def assign_df(self, df, dtype=dict()):
        assert isinstance(df, pd.DataFrame)

        self.sanitycheck_df(df)
        self._df = self.postprocess_df(df, dtype)

        self.update_df_lastaccess()

    def assign_gr(self, gr):
        assert isinstance(gr, pr.PyRanges)
        
        self.sanitycheck_gr(gr)
        self._gr = self.postprocess_gr(gr)

        self.update_gr_lastaccess()

    def assign_frame(self, frame, dtype=dict()):
        if isinstance(frame, pd.DataFrame):
            return self.assign_df(frame, dtype)
        elif isinstance(frame, pr.PyRanges):
            return self.assign_gr(frame)
        else:
            raise Exception(
                f'Input "frame" object must be either "pandas DataFrame" or "PyRanges" object.'
            )

    @staticmethod
    def sanitycheck_common(chromosomes, ends, chromdict):
        if not np.isin(chromosomes, chromdict.keys()).all():
            raise Exception(f'Unknown chromosome in the input frame')

        maxends = chromosomes.apply(lambda x: chromdict[x])
        if not (ends <= maxends).all():
            raise Exception(f'Some of "End" values are greater than the chromosome length.')

    @classmethod
    def sanitycheck_df(cls, df):
        if not set(genomedf_utils.COMMON_COLUMNS).issubset(df.columns):
            raise Exception(f'Required columns: {genomedf_utils.COMMON_COLUMNS}')

        assert (df['End'] > df['Start']).all()  # this works with 0-row DataFrame

        #self.sanitycheck_common(df['Chromosome'], df['End'], self.chromdict)

    @classmethod
    def postprocess_df(cls, df, dtype=dict()):
        # all processes work for 0-row dataframe

        # dtype setting
        astype_arg = genomedf_utils.DEFAULT_DTYPES | dtype
        astype_arg = {
            key: val for key, val in astype_arg.items()
            if key in df.columns
        }
        df = df.astype(astype_arg)

        # column reordering
        new_cols = genomedf_utils.COMMON_COLUMNS + [
            x for x in df.columns if x not in genomedf_utils.COMMON_COLUMNS
        ]
            # df.columns.difference(genomedf_utils.COMMON_COLUMNS).to_list() -> does not preserve ordering
        df = df.loc[:, new_cols]

        # remove index
        df = df.reset_index(drop=True, inplace=False)

        return df

    def sanitycheck_gr(self, gr):
        pass
        #self.sanitycheck_common(gr.Chromosome, gr.End, self.chromdict)

    @staticmethod
    def postprocess_gr(gr):
        if not gr.empty:
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

    @staticmethod
    def make_empty_df():
        result = pd.DataFrame([], columns=genomedf_utils.COMMON_COLUMNS)
        result = result.astype(genomedf_utils.DEFAULT_DTYPES)
        return result

    def set_gr_from_df(self):
        assert self._df is not None
        self._gr = pr.PyRanges(self._df, int64=True)

    def set_df_from_gr(self):
        assert self._gr is not None
        if self._gr.empty:
            if self._df is None:
                self._df = self.make_empty_df()
            else:
                self._df = self.df.iloc[[], :]
        else:
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

        # when gr or df is empty
        gr_is_empty = self._gr.empty
        df_is_empty = (self._df.shape[0] == 0)
        if gr_is_empty or df_is_empty:
            return (gr_is_empty == df_is_empty)

        # both gr and df are not empty
        if set(self._gr.columns) != set(self._df.columns):
            return False
        else:
            df_for_hash = self._df.loc[:, self._gr.columns]  # make column orders the same
            return (
                genomedf_utils.hash_gr(self._gr) 
                == genomedf_utils.hash_df(df_for_hash)
            ).to_numpy().all()  
                # It is faster (~2.5 us vs ~13 us) to convert into numpy and then do .all()

    def get_newer_frame(self):
        if self.check_gr_is_newer():
            assert self._gr is not None
            return self._gr
        else:
            assert self._df is not None
            return self._df

    ###################################
    # adaptations of pyranges methods #
    ###################################

    EXCLUDED_PYRANGES_METHODS = [
        'lengths',
        'assign',
    ]

    def _pyranges_method_adapter_alone(methodname):
        def new_method(self, *args, **kwargs):
            result_gr = getattr(self.gr, methodname)(*args, **kwargs)
            if result_gr.empty:
                return self.__class__.init_empty(
                    refver=self.refver, 
                    annot_cols=self.annot_cols,
                )
            else:
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

            if result_gr.empty:
                return self.__class__.init_empty(
                    refver=self.refver, 
                    annot_cols=self.annot_cols,
                )
            else:
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
        return self.df.__getitem__(key).to_numpy().copy()

    def __setitem__(self, key, val):
        if self._df is None:
            self.set_df_from_gr()
            
        self._df.__setitem__(key, np.asarray(val))
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
        'assign',
    ):
        new_mthd = dataframe_method_adaptor(mthdname)
        exec(f'{mthdname} = new_mthd')

    #################################
    # custom methods and properties #
    #################################

    @property
    def chromdict(self):
        return refgenome.get_chromdict(self.refver)

    ###############
    # coordinates #
    ###############

    @property
    def chromosomes(self):
        return self['Chromosome']

    chroms = chromosomes

    @property
    def chromosome_indexes(self):
        return np.fromiter(
            map(self.chromdict.contigs.index, self.chromosomes),
            dtype=int,
        )

    chromidxs = chromosome_indexes

    @property
    def starts(self):
        return self['Start']

    start0s = starts

    @property
    def ends(self):
        return self['End']

    end0s = ends

    def get_midpoints(self):
        return np.floor((self.start0s + self.end0s - 1) / 2).astype(int)

    @property
    def lengths(self):
        return self.end0s - self.start0s

    @property
    def lengths_cumsum(self):
        return np.cumsum(self.lengths)

    def check_self_overlap(self):
        if self.is_empty:
            return False
        else:
            bychrom = self.group_bychrom(how=None, sort=False)
            return any(
                genomedf_utils.check_interval_overlap(subgdf.start0s, subgdf.end0s)
                for chrom, subgdf in bychrom.items()
            )

    def check_duplicate_coords(self):
        if self.is_empty:
            return False
        else:
            return genomedf_utils.check_duplicate_coords(self.df)

    ##############################
    # basic dataframe attributes #
    ##############################

    @property
    def dtypes(self):
        return self.df.dtypes

    @property
    def shape(self):
        return self.df.shape

    @property
    def nrow(self):
        return self.df.shape[0]

    @property
    def ncol(self):
        return self.df.shape[1]

    @property
    def is_empty(self):
        newer_frame = self.get_newer_frame()
        if isinstance(newer_frame, pd.DataFrame):
            return newer_frame.shape[0] == 0
        else:
            return newer_frame.empty

    def copy(self):
        return self.spawn(self.df.copy())

    ######################
    # annotation columns #
    ######################

    @property
    def columns(self):
        if self.is_empty:
            return list(genomedf_utils.COMMON_COLUMNS)  # Chromosome, Start, End
        else:
            return self.get_newer_frame().columns.to_list()

    @property
    def annot_cols(self):
        return self.df.columns.drop(genomedf_utils.COMMON_COLUMNS).to_list()

    def choose_annots(self, annot_columns):
        annot_columns = np.atleast_1d(annot_columns)
        assert not set(annot_columns).intersection(genomedf_utils.COMMON_COLUMNS)
        return self.loc[:, genomedf_utils.COMMON_COLUMNS + list(annot_columns)]

    def drop_annots(self, columns=None, inplace=False):
        if columns is None:
            if inplace:
                self.assign_df(self.df.iloc[:, :3])
            else:
                return self.iloc[:, :3]
        else:
            columns = np.atleast_1d(columns)
            columns = set(columns).intersection(self.columns)
            #assert set(columns).issubset(self.columns)

            self_df = self.df
            if inplace:
                self_df.drop(columns, axis=1, inplace=True)
                self.assign_df(self_df)
            else:
                return self.spawn(self_df.drop(columns, axis=1, inplace=False))

    ##################
    # set operations #
    ##################

    def isec_union(self, other, nproc=1, drop_annots=False, sort=True):
        assert not self.check_self_overlap()
        assert not other.check_self_overlap()

        isec = self.intersect(other)
        self_diff_isec = self.subtract(isec)
        other_diff_isec = other.subtract(isec)
        if drop_annots:
            isec.drop_annots(inplace=True)
            self_diff_isec.drop_annots(inplace=True)
            other_diff_isec.drop_annots(inplace=True)
        else:
            cols_to_drop = ['Start_b', 'End_b', 'Distance']

            isec = isec.join(
                other, how='left', merge=None, find_nearest=False, nproc=nproc,
            )
            isec.drop_annots(cols_to_drop, inplace=True)

            self_diff_isec = self_diff_isec.join(
                other, how='left', find_nearest=True, merge=None, nproc=nproc,
            )
            self_diff_isec.drop_annots(cols_to_drop, inplace=True)

            other_diff_isec = other_diff_isec.join(
                self, how='left', find_nearest=True, merge=None, nproc=nproc,
            )
            other_diff_isec.drop_annots(cols_to_drop, inplace=True)

        result = self.__class__.concat([isec, self_diff_isec, other_diff_isec])
        if sort:
            result.sort()

        return result

    def subset_pyranges(self, chrom, start0=None, end0=None):
        """Returns a subset of original rows which overlaps the interval 
        specified by arguments. Original rows are returned as is, without 
        intersection with the input interval.
        'start0' and 'end0' arguments are 0-based & half-open.
        """
        getitem_key = (chrom, slice(start0, end0))
        return self.spawn(self.gr[getitem_key])

    def subset(self, chrom, start0=None, end0=None):
        """Returns a subset of original rows which overlaps the interval 
        specified by arguments. Original rows are returned as is, without 
        intersection with the input interval.
        'start0' and 'end0' arguments are 0-based & half-open.
        """
        isec_components = list()

        chrom_selector_idxs = self.get_chrom_selector(chrom, as_indexes=True)
        isec_components.append(chrom_selector_idxs)

        def helper(gdf_coords, query_pos, query_isend):
            relevant_gdf_coords = gdf_coords[chrom_selector_idxs]
            relevant_gdf_coords_selector = (
                (relevant_gdf_coords < query_pos)
                if query_isend else
                (relevant_gdf_coords > query_pos)
            )
            gdf_coords_selector_idxs = chrom_selector_idxs[relevant_gdf_coords_selector]
            return gdf_coords_selector_idxs

        if end0 is not None:
            isec_components.append(helper(self.start0s, end0, True))
        if start0 is not None:
            isec_components.append(helper(self.end0s, start0, False))

        final_selector_idxs = functools.reduce(np.intersect1d, isec_components)
        return self.iloc[final_selector_idxs, :]

    def subset_chroms(self, chromlist):
        return self.iloc[self.get_chrom_selector(chromlist, as_indexes=True), :]

    ########################
    # groupers & iterators #
    ########################

    def group_bychrom(self, how=None, sort=True):
        assert how in (None, 'gr', 'df')
        assert not self.is_empty

        result_asdf = dict(
            (chrom, subdf.reset_index(drop=True)) 
            for chrom, subdf in self.df.groupby(self.chroms, sort=False)
        )
        result_asgdf = {
            key: self.spawn(subdf) 
            for key, subdf in result_asdf.items()
        }
        if sort:
            for val in result_asgdf.values():
                val.sort()

        if how is None:
            return result_asgdf
        elif how == 'df':
            return {key: subgdf.df for key, subgdf in result_asgdf.items()}
        elif how == 'gr':
            return {key: subgdf.gr for key, subgdf in result_asgdf.items()}

    def get_chrom_selector(self, chromlist, as_indexes=False):
        chromlist = np.atleast_1d(chromlist)
        selector = np.isin(self.chromosomes, chromlist)
        if as_indexes:
            return np.nonzero(selector)[0]
        else:
            return selector

    def iter_coords(self):
        return zip(self.chromosomes, self.starts, self.ends)

    def get_coordcols_bychrom(self):
        self.sort()
        result = list()
        for (chrom, subdf) in self.df.iloc[:, :3].groupby(self.chroms, sort=False):
            #chroms = subdf['Chromosome'].to_numpy()
            start0s = subdf['Start'].to_numpy()
            end0s = subdf['End'].to_numpy()
            result.append((chrom, start0s, end0s))
        return result

    #def get_coord_groupkey(self):
    #    return genomedf_utils.get_coord_groupkey(self.df, self.chromdict)

    #########
    # split #
    #########

    @deco.get_deco_num_set_differently(('n', 'width'), 1)
    def equal_nrow_split(self, *, n=None, width=None):
        if n is not None:
            split_nums = tools.get_split_nums_byN(self.nrow, n)
        elif width is not None:
            split_nums = tools.get_split_nums_bywidth(self.nrow, width)

        indexes = np.insert(np.cumsum(split_nums), 0, 0)
        diff = np.diff(indexes)
        grouper = np.repeat(np.arange(len(diff)), diff)
        result = list(
            self.spawn(subdf) 
            for (key, subdf) in self.df.groupby(grouper, sort=False)
        )

        #result = list()
        #for start_idx, end_idx in tools.pairwise(indexes):
        #    result.append(self.iloc[start_idx:end_idx, :])

        return result

    @deco.get_deco_num_set_differently(('n', 'width'), 1)
    def equal_length_split(self, *, n=None, width=None):
        """Does NOT preserve original interval borders"""

        # get result_lengths_cumsum
        total_length = self.lengths_cumsum[-1]
        if n is not None:
            result_lengths_cumsum = np.cumsum(
                tools.get_split_nums_byN(total_length, n)
            )
        elif width is not None:
            result_lengths_cumsum = np.cumsum(
                tools.get_split_nums_bywidth(total_length, width)
            )

        # prepare result values
        result = list()
        for idx in range(len(result_lengths_cumsum)):
            # get start0 and end0 in the single merged coordinate system
            merged_start0 = (
                0 
                if idx == 0 else
                result_lengths_cumsum[idx - 1]
            )
            merged_end0 = result_lengths_cumsum[idx]
            # modify merged coordinates into interval-specific ones
            chrom_start, pos0_start, self_idx_start = self.equal_length_split_modify_coord(
                merged_start0, is_end=False,
            )
            chrom_end, pos0_end, self_idx_end = self.equal_length_split_modify_coord(
                merged_end0, is_end=True,
            )
            # create a GenomicDataFrame corresponding to a split unit
            gdf_data = {
                'chroms': list(),
                'start0s': list(),
                'end0s': list(),
            }

            if self_idx_start == self_idx_end:
                gdf_data['chroms'].append(chrom_start)
                gdf_data['start0s'].append(pos0_start)
                gdf_data['end0s'].append(pos0_end)
            else:
                gdf_data['chroms'].append(chrom_start)
                gdf_data['start0s'].append(pos0_start)
                gdf_data['end0s'].append(self.ends[self_idx_start])

                for self_idx in range(self_idx_start + 1, self_idx_end):
                    gdf_data['chroms'].append(self.chromosomes[self_idx])
                    gdf_data['start0s'].append(self.starts[self_idx])
                    gdf_data['end0s'].append(self.ends[self_idx])

                gdf_data['chroms'].append(chrom_end)
                gdf_data['start0s'].append(self.starts[self_idx_end])
                gdf_data['end0s'].append(pos0_end)

            result.append(self.spawn_from_data(**gdf_data))

        return result

    def equal_length_split_modify_coord(self, merged_pos0, is_end=False):
        lengths_cumsum = self.lengths_cumsum
        def get_interval_values(merged_pos0, idx):
            chrom = self.chromosomes[idx]

            if idx == 0:
                shift_within_interval = merged_pos0
            else:
                shift_within_interval = (
                    merged_pos0 
                    - lengths_cumsum[idx - 1]
                )
            new_pos0 = self.starts[idx] + shift_within_interval

            return chrom, new_pos0

        def handle_current_interval(
            merged_pos0, idx, length_previous, length_current, is_end,
        ):
            if (
                (merged_pos0 > length_previous)
                and (merged_pos0 <= length_current)
            ):
                if merged_pos0 == length_current:
                    if is_end:
                        chrom, new_pos0 = get_interval_values(merged_pos0, idx)
                    else:
                        chrom, new_pos0 = get_interval_values(merged_pos0, idx + 1)
                elif merged_pos0 < length_current:
                    chrom, new_pos0 = get_interval_values(merged_pos0, idx)

                to_break = True
            else:
                chrom = None
                new_pos0 = None
                to_break = False

            return chrom, new_pos0, to_break

        # sanity check
        assert merged_pos0 >= 0, f'"merged_pos0" must be non-negative.'
        assert not (
            (not is_end) 
            and (merged_pos0 >= lengths_cumsum[-1])
        ), (
            f'If ""is_end" is False, "merged_pos0" must be less than '
            f'the total IntervalList length.'
        )
        assert not (is_end and merged_pos0 > lengths_cumsum[-1]), (
            f'If ""is_end" is True, "merged_pos0" must be less than or '
            f'equal to the total IntervalList length.'
        )

        # main
        if merged_pos0 == 0:
            chrom = self.chromosomes[0]
            new_pos0 = self.starts[0]
            self_idx = 0
        else:
            for idx in range(len(lengths_cumsum)):
                length_current = lengths_cumsum[idx]
                length_previous = (
                    0
                    if idx == 0 else
                    lengths_cumsum[idx - 1]
                )
                chrom, new_pos0, to_break = handle_current_interval(
                    merged_pos0, idx, length_previous, length_current, is_end
                )
                self_idx = idx

                if to_break:
                    break

        return chrom, new_pos0, self_idx

    ################
    # segmentation #
    ################

    def make_compact(self):
        compact_pos0_list, decoder = genomedf_utils.compact_coords(self.start0s, self.end0s)
        compact_gdf = self.copy()
        compact_gdf['Start'] = compact_pos0_list
        compact_gdf['End'] = compact_gdf['Start'] + 1

        return compact_gdf, decoder

    @classmethod
    def get_segment_base(
        cls,
        gdf, 
        annot_colname,
        N_colname='N',
        smoothing=False, 
        smooth_kwargs=dict(), 
        segment_kwargs=dict(),
    ):
        # set kwargs
        segment_kwargs = (
            dict()
            | segment_kwargs
        )

        # make compact and run segmentation
        compact_gdf, decoder = gdf.make_compact()
        compact_seg_df = rdnacopy.run_dnacopy_segment(
            chromosomes=compact_gdf.chroms,
            positions=compact_gdf.start0s,
            values=compact_gdf[annot_colname],

            N_colname='N',
            value_colname='value',
            smoothing=smoothing, 
            verbose=False, 
            smooth_kwargs=smooth_kwargs, 
            segment_kwargs=segment_kwargs,
        )
        compact_seg_gdf = gdf.spawn(compact_seg_df)

        # recover noncompact coordinates
        noncompact_start0s, noncompact_end0s = genomedf_utils.uncompact_coords(
            compact_seg_gdf.starts, 
            compact_seg_gdf.ends, 
            decoder,
        )

        # make result
        result = GenomeDataFrame.from_data(
            refver=gdf.refver,
            **{
                'chroms': compact_seg_gdf.chromosomes,
                'start0s': noncompact_start0s,
                'end0s': noncompact_end0s,
                annot_colname: compact_seg_gdf['value'],
                N_colname: compact_seg_gdf['N'],
            }
        )
        return result

    def get_segment(
        self, 
        annot_colname,
        nproc=1,
        split_width=1000,

        N_colname='N',
        smoothing=False, 
        verbose=False, 
        smooth_kwargs=dict(), 
        segment_kwargs=dict(),
    ):
        assert not self.is_empty

        if verbose:
            logutils.log(f'Beginning segmentation')

        #grouper = self.chroms  # this may be more divided by cytoband features
        #grouped_gdfs = list(
        #    self.spawn(val) 
        #    for (key, val) in self.df.groupby(grouper, sort=False)
        #)

        self.sort()
        grouped_bychrom = self.group_bychrom(sort=False)
        grouped_gdfs = list()
        for chrom, subgdf in grouped_bychrom.items():
            grouped_gdfs.extend(subgdf.equal_nrow_split(width=split_width))

        with multiprocessing.Pool(nproc) as pool:
            args = (
                (
                    gdf, 
                    annot_colname,
                    N_colname,
                    smoothing,
                    smooth_kwargs,
                    segment_kwargs,
                )
                for gdf in grouped_gdfs
            )
            pool_result = pool.starmap(GenomeDataFrame.get_segment_base, args)

        result = GenomeDataFrame.concat(pool_result)
        result.sort()

        if verbose:
            logutils.log(f'Finished segmentation')

        return result

    ########
    # join #
    ########

    @deco.get_deco_arg_choices({'how': [None, 'left']})
    def join(
        self, other, 
        right_gdf_cols=None,
        how='left', 
        merge=['mean', 'std'],
        ddof=0,
        overlapping_length=False,
        N_colname='N',
        suffixes=None,
        merge_args=tuple(),
        merge_kwargs=dict(),
        winsorize=None,
        omit_N=False,
        find_nearest=False,

        split_width=3000,
        nproc=1,
    ):
        join_kwargs = dict(
            right_gdf_cols=right_gdf_cols,
            how=how,
            merge=merge,
            ddof=ddof,
            overlapping_length=overlapping_length,
            N_colname=N_colname,
            suffixes=suffixes,
            merge_args=merge_args,
            merge_kwargs=merge_kwargs,
            winsorize=winsorize,
            omit_N=omit_N,
            find_nearest=find_nearest,
        )

        if (
            (nproc == 1)
            or self.is_empty
            or other.is_empty
        ):
            return genomedf_methods.join_base(self, other, **join_kwargs)
        else:
            partial_left_gdf_list = genomedf_methods.split_left_gdf(self, width=split_width)
            with multiprocessing.Pool(nproc) as pool:
                args = (
                    (partial_left_gdf, other, join_kwargs)
                    for partial_left_gdf in partial_left_gdf_list
                )
                mp_result = pool.starmap(genomedf_methods.targetfunc, args)

            return self.__class__.concat(mp_result)

    ###########
    # binning #
    ###########

    def get_binsize(self):
        assert not self.is_empty

        self.sort()

        gdfs_bychrom = self.group_bychrom(how=None, sort=False)
        equal_length_portions = list()
        last_ones = list()
        for chrom, subgdf in gdfs_bychrom.items():
            lengths = subgdf.lengths
            last_ones.append(lengths[-1])
            if len(lengths) > 1:
                equal_length_portions.append(lengths[:-1])

        all_lengths = np.unique(np.concatenate(equal_length_portions))
        if all_lengths.shape != (1,):
            raise NotBinnedError(f'Bin lengths are not identical.')
        binsize = all_lengths[0]

        if len(last_ones) > 0:
            last_ones = np.asarray(last_ones)
            if not (last_ones <= binsize).all():
                raise NotBinnedError(f'Some last interval lengths are greater than the bin length.')

        return binsize

    def check_binned(self, return_binsize=False):
        try:
            binsize = self.get_binsize()
        except NotBinnedError:
            if return_binsize:
                return False, None
            else:
                return False
        else:
            if return_binsize:
                return True, binsize
            else:
                return True

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
        assert not self.is_empty

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

        gdfs_bychrom = self.group_bychrom()
        newbin_dfs_bychrom = dict()
        for chrom, singlechrom_gdf in gdfs_bychrom.items():
            singlechrom_df = singlechrom_gdf.df
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

    ######################
    # other modification #
    ######################

    def sort(self):
        """In-place sort, only applied to df, not gr"""
        if self.is_empty:
            return

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

    def fill_gaps(self, edit_first_last=True):
        self.sort()

        new_start0s = list()
        new_end0s = list()
        for (chrom, old_start0s, old_end0s) in self.get_coordcols_bychrom():
            # fill borders between regions
            new_borders = np.floor(
                0.5 * (old_end0s[:-1] + old_start0s[1:])
            ).astype(int)
            current_new_start0s = np.insert(new_borders, 0, old_start0s[0])
            current_new_end0s = np.append(new_borders, old_end0s[-1])

            # modify start of the first region and end of the last region
            if edit_first_last:
                current_new_start0s[0] = 0
                current_new_end0s[-1] = self.chromdict[chrom]

            # apply
            new_start0s.append(current_new_start0s)
            new_end0s.append(current_new_end0s)

        new_start0s = np.concatenate(new_start0s)
        new_end0s = np.concatenate(new_end0s)

        result = self.spawn_from_data(
            chroms=self.chroms,
            start0s=new_start0s,
            end0s=new_end0s,
            **{key: self[key] for key in self.annot_cols},
        )
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



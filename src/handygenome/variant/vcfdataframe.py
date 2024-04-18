import re
import sys
import os
import itertools
import collections
import logging
import multiprocessing

import numpy as np
import pandas as pd
import pysam

import handygenome.deco as deco
import handygenome.tools as tools
import handygenome.refgenome.refgenome as refgenome
import handygenome.cnv.cnvcall as cnvcall
from handygenome.genomedf.genomedf import GenomeDataFrame
from handygenome.variant.variantplus import VariantPlus, VariantPlusList


class VariantDataFrame(GenomeDataFrame):
    """Mandatory columns: Chromosome, Start, End, REF, ALT1
    """
    COMMON_ALLELE_COLUMNS = ['REF', 'ALT1']
    COMMON_COLUMNS = GenomeDataFrame.COMMON_COLUMNS + COMMON_ALLELE_COLUMNS
    ALLELE_COLNAME_PAT = re.compile(r'(REF|ALT[1-9][0-9]*)')

    DEFAULT_DTYPES = (
        GenomeDataFrame.DEFAULT_DTYPES
        | {
            'REF': 'string',
            re.compile('ALT[1-9][0-9]*'): 'string',
        }
    )

    @property
    def POS(self):
        return self['Start'] + 1

    @staticmethod
    def allele_column_sortkey(colname):
        if colname == 'REF':
            return 0
        else:
            return int(re.sub('^ALT', '', colname))

    @property
    def allele_columns(self):
        result = [
            x for x in self.columns
            if self.__class__.ALLELE_COLNAME_PAT.fullmatch(x)
        ]
        result.sort(key=self.allele_column_sortkey)
        return result

    @property
    def vaf_columns(self):
        result = [x for x in self.columns if x.endswith('_vaf')]
        result.sort(
            key=(lambda x:
                self.allele_column_sortkey(re.sub('_vaf$', '', x))
            )
        )
        return result

    @property
    def nonannot_columns(self):
        return list(GenomeDataFrame.COMMON_COLUMNS) + self.allele_columns

    @property
    def num_alleles(self):
        return len(self.allele_columns)

    ##################
    # set operations #
    ##################

    def get_uids(self):
        allele_values = self[self.allele_columns]
        ref_list = list(allele_values[:, 0])
        alts_list = list()
        for row in allele_values[:, 1:]:
            alts_list.append(
                tuple(row[pd.notna(row)])
            )

        values = list()
        for chrom, pos, ref, alts in zip(
            self.chroms, self.POS, ref_list, alts_list,
        ):
            values.append(
                libvcfspec.make_vcfspec_id(
                    chrom, pos, ref, alts
                )
            )
        return pd.Index(values)

    def get_isec_other_indexes(self, other):
        self_uids = self.get_uids()
        other_uids = other.get_uids()
        joined_uids, self_indexes, other_indexes = self_uids.join(other_uids, how='left', return_indexers=True)
        if other_indexes is None:  # when all elements should be selected, "self_indexes" or "other_indexes" becomes None
            other_indexes = np.arange(len(other_uids))
            
        result = np.where((other_indexes == -1), np.nan, other_indexes)
        return result

    @classmethod
    def from_var_indexed_df(cls, df, refver):
        index_columns = df.index.names
        assert set(['Chromosome', 'Start']).issubset(index_columns)
        non_coord_cols = set(index_columns).difference(['Chromosome', 'Start'])
        assert all(
            cls.ALLELE_COLNAME_PAT.fullmatch(x) 
            for x in non_coord_cols
        )
        assert set(cls.COMMON_ALLELE_COLUMNS).issubset(non_coord_cols)

        new_df = df.reset_index()
        new_df['End'] = new_df['Start'] + new_df['REF'].str.len()

        return cls.from_frame(new_df, refver=refver)

    def variant_join(self, other, how='left', nproc=1):
        self_varidx_df = self.get_var_indexed_df()
        other_varidx_df = other.get_var_indexed_df()
        joined_df = self_varidx_df.join(other_varidx_df, how=how)
        return self.__class__.from_var_indexed_df(joined_df, refver=self.refver)

    @classmethod
    def variant_union(cls, gdfs, nproc=1, num_split=30):
        assert len(gdfs) >= 2
        assert len(set(x.refver for x in gdfs)) == 1
        assert len(set(x.num_alleles for x in gdfs)) == 1

        logutils.log(f'Splitting input gdfs')
        #gdfs_bychrom = [x.group_bychrom(sort=False) for x in gdfs]
        with multiprocessing.Pool(nproc) as pool:
            gdfs_bychrom = pool.map(cls.split_input_gdfs_targetfunc2, gdfs)

        logutils.log(f'Joining split gdfs')
        split_args = list()
        all_chroms = set(
            itertools.chain.from_iterable(x.keys() for x in gdfs_bychrom)
        )
        for chrom in all_chroms:
            sublist = list()
            for x in gdfs_bychrom:
                if chrom in x:
                    sublist.append(x[chrom])
            split_args.append(sublist)

        with multiprocessing.Pool(nproc) as pool:
            mp_result = pool.map(cls.variant_union_base, split_args)

        result = cls.concat(mp_result)
        result.sort()
        return result

    @staticmethod
    def split_input_gdfs_targetfunc2(gdf):
        return gdf.group_bychrom(sort=False)

    @classmethod
    def variant_union_base(cls, gdfs):
        var_indexed_df_list = [x.get_var_indexed_df() for x in gdfs]
        union_df = var_indexed_df_list[0].join(var_indexed_df_list[1:], how='outer')
        result = cls.from_var_indexed_df(union_df, refver=gdfs[0].refver)
        result.sort()
        return result
       

class VCFDataFrame(VariantDataFrame):
    """Support sample-wise annotation columns
    """
    SAMPLEID_SEP = ':::'
    CCF_COLNAME = 'ccf'
    CNm_COLNAME = 'CNm'
    PON_PREFIX = '<PON>'

    ################
    # initializers #
    ################

    @classmethod
    def from_vp(cls, vp, n_allele=2):
        return cls.from_vplist(VariantPlusList([vp]), n_allele=n_allele)

    @classmethod
    @deco.get_deco_atleast1d(['sampleids'], keep_none=True)
    @deco.get_deco_arg_choices({'mode': ['common_only', 'samplewise_1']})
    def from_vplist(
        cls, vplist, 
        chrom=None, start0=None, end0=None,
        n_allele=2,
        mode='samplewise_1',
        sampleids=None,
        exclude_other=False,
    ):
        # make parameters
        all_samples = list(vplist[0].vr.samples.keys())
        vp_iterator = vplist.get_vp_iter_from_vcf(
            chrom=chrom, start0=start0, end0=end0, 
        )
        vr_iterator = (vp.vr for vp in vplist)
        refver = vplist.refver
        vp0 = vplist[0]

        # main
        result = cls.init_common(
            sampleids, all_samples, vp0,
            n_allele, vr_iterator,
            mode, 
            vp_iterator,
            refver,
        )

        return result

    @classmethod
    @deco.get_deco_atleast1d(['sampleids'], keep_none=True)
    @deco.get_deco_arg_choices({'mode': ['common_only', 'samplewise_1']})
    def from_vcfpath(
        cls, vcf_path, 
        chrom=None, start0=None, end0=None,
        n_allele=2,
        mode='samplewise_1',
        sampleids=None,
        exclude_other=False,
    ):
        # make parameters
        in_vcf = pysam.VariantFile(vcf_path)
        all_samples = list(in_vcf.header.samples)
        vplist = VariantPlusList.from_vcf_lazy(vcf_path)
        vp_iterator = vplist.get_vp_iter_from_vcf(
            chrom=chrom, start0=start0, end0=end0, 
        )
        vr_iterator = in_vcf.fetch()
        refver = refgenome.infer_refver_vcfheader(in_vcf.header)
        vp0 = VariantPlus.from_vr(next(in_vcf.fetch()))

        # main
        result = cls.init_common(
            sampleids, all_samples, vp0,
            n_allele, vr_iterator,
            mode, 
            vp_iterator,
            refver,
        )

        return result

    @classmethod
    def init_common(
        cls, 
        sampleids, all_samples, vp0,
        n_allele, vr_iterator,
        mode, 
        vp_iterator,
        refver,
    ):
        # args postprocess
        (
            sampleids, 
            pon_samples, 
            nonpon_samples,
        ) = cls.postprocess_sampleids(sampleids, all_samples, vp0)
        if n_allele is None:
            n_allele = max(len(vr.alleles) for vr in vr_iterator)

        # main
        allele_columns = cls.get_allele_columns(n_allele)
        all_columns = cls.make_columns(
            allele_columns, mode, pon_samples, nonpon_samples,
        )
        data = list()
        for row_dict in cls.iter_row_dict(
            vp_iterator=vp_iterator,
            mode=mode,
            allele_columns=allele_columns,
            sampleids=sampleids,
            pon_samples=pon_samples,
            nonpon_samples=nonpon_samples,
        ):
            data.append(row_dict)
        src_df = pd.DataFrame.from_records(data, columns=all_columns)
        result = cls.from_frame(src_df, refver=refver)

        return result

    def make_rename_dic(cls, all_columns, pon_samples):
        rename_dic = dict()
        for x in all_columns:
            if any(
                x.startswith(f'{pon_sid}{cls.SAMPLEID_SEP}') 
                for pon_sid in pon_samples
            ):
                rename_dic[x] = cls.PON_PREFIX + x

    #######################
    # helpers - immediate #
    #######################

    @classmethod
    def make_columns(cls, allele_columns, mode, pon_samples, nonpon_samples):
        if mode == 'samplewise_1':
            vaf_columns = cls.get_sw_columns_vaf(pon_samples, nonpon_samples, allele_columns)
            other_columns = cls.get_sw_columns_others(nonpon_samples)
        else:
            vaf_columns = list()
            other_columns = list()

        columns = (
            ['Chromosome', 'Start', 'End']
            + allele_columns
            + ['POS', 'ID', 'is_bnd1']
            + vaf_columns
            + other_columns
        )
        return columns

    @classmethod
    def iter_row_dict(
        cls,
        vp_iterator,
        mode,
        allele_columns,
        sampleids,
        pon_samples,
        nonpon_samples,
    ):
        for vp in vp_iterator:
            row_dict = cls.make_row_from_vp(
                vp=vp, 
                allele_columns=allele_columns, 
                sampleids=sampleids, 
                nonpon_samples=nonpon_samples, 
                pon_samples=pon_samples,
                mode=mode,
            )
            yield row_dict

    @classmethod
    def make_sourcedata(
        cls, 
        vp_iterator,
        mode,
        allele_columns,
        sampleids,
        pon_samples,
        nonpon_samples,
    ):
        data = list()
        for vp in vp_iterator:
            row_dict = cls.make_row_from_vp(
                vp=vp, 
                allele_columns=allele_columns, 
                sampleids=sampleids, 
                nonpon_samples=nonpon_samples, 
                pon_samples=pon_samples,
                mode=mode,
            )
            data.append(row_dict)

        return data
       
    @classmethod
    def make_row_from_vp(
        cls, vp, allele_columns, sampleids, nonpon_samples, mode,
        pon_samples,
    ):
        row_dict = dict()
        cls.populate_rowdata_common_basic(vp, row_dict, allele_columns)
        if mode == 'samplewise_1':
            # exclude_other is fixed as False
            cls.populate_rowdata_samplewise_basic(
                vp=vp, 
                row_dict=row_dict,
                allele_columns=allele_columns,
                sampleids=sampleids, 
                pon_samples=pon_samples,
                nonpon_samples=nonpon_samples,
            )
        row_dict['POS'] = row_dict['Start'] + 1
        return row_dict

    ##################
    # helpers - misc #
    ##################

    @classmethod
    def postprocess_sampleids(cls, sampleids, all_samples, vp0):
        if sampleids is None:
            sampleids = list(all_samples)

        assert not any(x.startswith(cls.PON_PREFIX) for x in sampleids), (
            f'All sample IDs must not begin with the preset PON prefix {repr(cls.PON_PREFIX)}'
        )

        pon_samples = [key for key, val in vp0.readstats.items() if val['is_pon']]
        nonpon_samples = [x for x in sampleids if x not in pon_samples]
        #sampleids = [
        #    (f'{cls.PON_PREFIX}{x}' if (x in pon_samples) else x) 
        #    for x in sampleids
        #]

        return sampleids, pon_samples, nonpon_samples

    ###############################
    # helpers - make column names #
    ###############################

    @staticmethod
    def get_allele_columns(n_allele):
        return ['REF'] + [f'ALT{x + 1}' for x in range(n_allele - 1)]

    @classmethod
    def get_sw_columns_vaf(cls, pon_samples, nonpon_samples, allele_columns):
        """
            "sw" means sample-wise

            [
                sid1:::REF, sid1:::ALT1, sid1:::ALT2,
                sid2:::REF, sid2:::ALT1, sid2:::ALT2,
                ...
            ]
        """
        result = list()
        for sid in nonpon_samples:
            for allele_col in allele_columns:
                result.append(f'{sid}{cls.SAMPLEID_SEP}{allele_col}_vaf')
        for sid in pon_samples:
            for allele_col in allele_columns:
                result.append(
                    f'{cls.PON_PREFIX}{sid}{cls.SAMPLEID_SEP}{allele_col}_vaf'
                )
        return result

    @classmethod
    def get_sw_columns_others(cls, nonpon_samples):
        persample_names = [
            'TD', 'RD', 'AD', 'OTHER',
            'CLIPPED_RD', 'CLIPPED_AD',
            'RDF1R2', 'RDF2R1', 'ADF1R2', 'ADF2R1',
            #'SB',
            'mean_mNM_ref', 'mean_mNM_alt', 'mean_NM_ref', 'mean_NM_alt',
            'mean_RP_alt', 'std_RP_alt',
            'mean_MQ_ref', 'mean_MQ_alt',
            'ID_context_ref', 'ID_context_alt',
        ]
        result = list()
        for sid in nonpon_samples:
            for x in persample_names:
                result.append(f'{sid}{cls.SAMPLEID_SEP}{x}')

        return result
        

    ######################
    # helpers - populate #
    ######################

    @classmethod
    def populate_rowdata_common_basic(cls, vp, row_dict, allele_columns):
        if vp.vcfspec.check_is_sv():
            bnds = vp.vcfspec.get_bnds()
            tmpvp_bnd1 = VariantPlus.from_vcfspec(bnds.get_vcfspec_bnd1())
            tmpvp_bnd2 = VariantPlus.from_vcfspec(bnds.get_vcfspec_bnd2())
            uid = vcfspec_bnd1.get_id()
            
            for vp in [tmpvp_bnd1, tmpvp_bnd2]:
                cls.populate_rowdata_common_basic_helper(
                    vp=vp, 
                    row_dict=row_dict,
                    allele_columns=allele_columns,
                    vcfspec_id=uid,
                )
        else:
            cls.populate_rowdata_common_basic_helper(
                vp=vp, 
                row_dict=row_dict,
                allele_columns=allele_columns,
            )

    @staticmethod
    def populate_rowdata_common_basic_helper(
        vp, 
        row_dict,
        allele_columns,
        vcfspec_id=None,
    ):
        # positions
        row_dict['Chromosome'] = vp.vcfspec.chrom
        row_dict['Start'] = vp.vcfspec.start0
        row_dict['End'] = vp.vcfspec.end0
        # alleles 
        for idx, allele_colname in enumerate(allele_columns):
            row_dict[allele_colname] = vp.get_allele(idx)
        # unique ID
        if vcfspec_id is None:
            vcfspec_id = vp.vcfspec.get_id()
        row_dict['ID'] = vcfspec_id
        # is_bnd1
        if vp.vcfspec.check_is_sv():
            row_dict['is_bnd1'] = vp.vcfspec.check_is_bnd1()
        else:
            row_dict['is_bnd1'] = pd.NA

    @classmethod
    def populate_rowdata_samplewise_basic(
        cls, 
        vp, 
        row_dict,
        allele_columns,
        sampleids, 
        pon_samples,
        nonpon_samples,
    ):
        # VAF
        vafs = vp.get_vafs(
            sampleids, n_allele=len(allele_columns), exclude_other=False,
        )
        for sid in nonpon_samples:
            for allele_name, vafval in zip(allele_columns, vafs[sid]):
                key = f'{sid}{cls.SAMPLEID_SEP}{allele_name}_vaf'
                row_dict[key] = vafval
        for sid in pon_samples:
            for allele_name, vafval in zip(allele_columns, vafs[sid]):
                key = f'{cls.PON_PREFIX}{sid}{cls.SAMPLEID_SEP}{allele_name}_vaf'
                row_dict[key] = vafval

        # others
        for sid in nonpon_samples:
            rdst = vp.readstats[sid]
            # TD, RD, AD, OTHER
            row_dict[f'{sid}{cls.SAMPLEID_SEP}TD'] = rdst.get_total_rppcount(exclude_other=False)
            row_dict[f'{sid}{cls.SAMPLEID_SEP}RD'] = rdst['rppcounts'][0]
            row_dict[f'{sid}{cls.SAMPLEID_SEP}AD'] = rdst['rppcounts'][1]
            row_dict[f'{sid}{cls.SAMPLEID_SEP}OTHER'] = rdst['rppcounts'][-1]
            # CLIPPED_RD, AD
            row_dict[f'{sid}{cls.SAMPLEID_SEP}CLIPPED_RD'] = rdst['clipped_rppcounts'][0]
            row_dict[f'{sid}{cls.SAMPLEID_SEP}CLIPPED_AD'] = rdst['clipped_rppcounts'][1]
            # RDF1R2, RDF2R1, ADF1R2, ADF2R1
            row_dict[f'{sid}{cls.SAMPLEID_SEP}RDF1R2'] = rdst['f1r2_rppcounts'][0]
            row_dict[f'{sid}{cls.SAMPLEID_SEP}RDF2R1'] = rdst['f2r1_rppcounts'][0]
            row_dict[f'{sid}{cls.SAMPLEID_SEP}ADF1R2'] = rdst['f1r2_rppcounts'][1]
            row_dict[f'{sid}{cls.SAMPLEID_SEP}ADF2R1'] = rdst['f2r1_rppcounts'][1]
            # mNM, NM
            row_dict[f'{sid}{cls.SAMPLEID_SEP}mean_mNM_ref'] = rdst['mNM'][0]
            row_dict[f'{sid}{cls.SAMPLEID_SEP}mean_mNM_alt'] = rdst['mNM'][1]
            row_dict[f'{sid}{cls.SAMPLEID_SEP}mean_NM_ref'] = rdst['mean_NMs'][0]
            row_dict[f'{sid}{cls.SAMPLEID_SEP}mean_NM_alt'] = rdst['mean_NMs'][1]
            # RP
            row_dict[f'{sid}{cls.SAMPLEID_SEP}mean_RP_alt'] = rdst['mean_varpos_fractions'][1]
            row_dict[f'{sid}{cls.SAMPLEID_SEP}std_RP_alt'] = rdst['std_varpos_fractions'][1]
            # MQ
            row_dict[f'{sid}{cls.SAMPLEID_SEP}mean_MQ_ref'] = rdst['mean_MQs'][0]
            row_dict[f'{sid}{cls.SAMPLEID_SEP}mean_MQ_alt'] = rdst['mean_MQs'][1]
            # ID_context
            row_dict[f'{sid}{cls.SAMPLEID_SEP}ID_context_ref'] = rdst['IDcontext_rppcounts'][0]
            row_dict[f'{sid}{cls.SAMPLEID_SEP}ID_context_alt'] = rdst['IDcontext_rppcounts'][1]

    #######
    # I/O #
    #######

    def write_vcftsv(self, path):
        assert path.endswith('.tsv.gz')

        df_to_write = self.df.copy()
        df_to_write.insert(1, 'POS', self.POS) 
        df_to_write.rename({'Chromosome': 'CHROM'}, axis=1, inplace=True)
        df_to_write.drop(['Start', 'End', 'index'], axis=1, inplace=True)
        df_to_write.to_csv(path, sep='\t', header=True, index=False, na_rep='NA')
            
    ##############
    # properties #
    ##############

    @property
    def nonannot_columns(self):
        return self.colname_getter(self.colnamefilter_nonannot)

    @property
    def info_columns(self):
        return self.colname_getter(self.colnamefilter_info)

    @property
    def format_columns(self):
        return self.colname_getter(self.colnamefilter_format)

    @property
    def samples(self):
        return tools.unique_keeporder(
            [
                x.split(self.__class__.SAMPLEID_SEP)[0] 
                for x in self.format_columns
            ]
        )

    @property
    def format_keys(self):
        return tools.unique_keeporder(
            [
                x.split(self.__class__.SAMPLEID_SEP)[1] 
                for x in self.format_columns
            ]
        )

    #######################
    # column name filters #
    #######################

    def colname_getter(self, colnamefilter):
        return tools.unique_keeporder(filter(colnamefilter, self.columns))

    @classmethod
    def colnamefilter_nonannot(cls, x):
        return (
            (x in cls.COMMON_COLUMNS)
            or (re.fullmatch('ALT[0-9]+', x) is not None)
        )

    @classmethod
    def colnamefilter_info(cls, x):
        return (
            (cls.SAMPLEID_SEP not in x)
            and (not cls.colnamefilter_nonannot(x))
        )

    @classmethod
    def colnamefilter_format(cls, x):
        return (
            (cls.SAMPLEID_SEP in x)
            and (not cls.colnamefilter_nonannot(x))
        )

    ###############
    # sanitycheck #
    ###############
        
    @classmethod
    def sanitycheck_df(cls, df):
        super().sanitycheck_df(df)
        assert all(
            (cls.SAMPLEID_SEP not in x)
            for x in df.columns
            if cls.colnamefilter_nonannot(x)
        )
        assert all(
            (x.count(cls.SAMPLEID_SEP) == 1)
            for x in df.columns
            if cls.colnamefilter_format(x)
        )

    ###########
    # getters #
    ###########

    def get_format_colname(self, sample, key):
        return sample + self.__class__.SAMPLEID_SEP + key

    def get_format(self, samples, keys):
        samples = np.atleast_1d(samples)
        keys = np.atleast_1d(keys)
        selected_colnames = list()
        for s, k in itertools.product(samples, keys):
            colname = s + self.__class__.SAMPLEID_SEP + k
            selected_colnames.append(colname)

        if len(selected_colnames) == 1:
            selected_colnames = selected_colnames[0]
            return self[selected_colnames]
        else:
            return np.stack(
                [self.loc[:, x] for x in selected_colnames],
                axis=1,
            )
            #return self.loc[:, selected_colnames]

    def get_sample_annots(self, sample):
        colnames = [
            x for x in self.format_columns 
            if x.split(self.__class__.SAMPLEID_SEP)[0] == sample
        ]
        return self.loc[:, colnames]

    #######
    # CCF #
    #######

    def add_ccf(
        self, cnv_gdf, cellularity, vcf_sampleid, 
        clonal_pval=cnvcall.DEFAULT_CLONAL_PVALUE,
    ):
        assert 'total_depth' in self.columns
        assert set(['CN', 'B_baf0', 'gCN']).issubset(cnv_gdf.columns)

        CNt_colname = cnv_gdf.get_clonal_CN_colname(germline=False)
        B_colname = cnv_gdf.get_clonal_B_colname('baf0', germline=False)
        CNg_colname = cnv_gdf.get_clonal_CN_colname(germline=True)

        joined_df = self.join(
            cnv_gdf, 
            [CNt_colname, B_colname, CNg_colname], 
            how='left',
            merge='longest', 
            suffixes={'longest': ''},
        )
        valid_flag = (joined_df.get_format(vcf_sampleid, 'ALT1_depth') > 0)

        CNm, ccf = cnvcall.find_ccf(
            vaf=joined_df.get_format(vcf_sampleid, 'ALT1_vaf'),
            total_depth=joined_df['total_depth'],
            alt_depth=np.where(
                valid_flag, 
                joined_df.get_format(vcf_sampleid, 'ALT1_depth'),
                np.nan,
            ),
            CNg=joined_df[CNg_colname],
            CNt=joined_df[CNt_colname],
            Bt=joined_df[B_colname],
            cellularity=cellularity,
            clonal_pval=clonal_pval,
        )

        self[self.__class__.CCF_COLNAME] = ccf
        self[self.__class__.CNm_COLNAME] = CNm
        self[CNt_colname] = joined_df[CNt_colname]
        self[B_colname] = joined_df[B_colname]
        self[CNg_colname] = joined_df[CNg_colname]

    ########
    # draw #
    ########

    def add_vaf_legend_handle(self, handles):
        handles.add_line(
            marker='o', 
            linewidth=0, 
            color='black', 
            markersize=4,
            label='vaf'
        )

    def add_ccf_legend_handle(self, handles):
        handles.add_line(
            marker='x', 
            linewidth=0, 
            color='tab:red', 
            markersize=4,
            label='ccf'
        )

    def draw_vaf(
        self,
        ax=None,
        genomeplotter=None,
        plotdata=None,

        sampleid=None,
        vaf_legend=True,
        ccf_legend=False,

        # fig generation params
        title=None,
        suptitle_kwargs=dict(),
        subplots_kwargs=dict(),

        # drawing kwargs
        plot_kwargs=dict(),

        # axes setting
        setup_axes=True,
        ylabel=False,
        ylabel_prefix='',
        ylabel_kwargs=dict(),
        ymax=False,
        ymin=False,
        yticks=False,
        draw_common_kwargs=dict(),
        rotate_chromlabel=None,

        # multicore plotdata generation
        nproc=1,
        verbose=True,
    ):
        if ymax is False:
            ymax = 1.1
        if ymin is False:
            ymin = -0.1
        if ylabel is False:
            ylabel = 'Mutation'
        if yticks is False:
            yticks = np.round(np.arange(0, 1.2, 0.2), 1)
        if sampleid is None:
            sampleid = self.samples[0]

        # prepare single plotdata
        if genomeplotter is None:
            genomeplotter = self.get_default_genomeplotter()
        if plotdata is None:
            plotdata = genomeplotter.make_plotdata(
                self, 
                log_suffix=' (mutation vaf)',
                nproc=nproc,
                verbose=verbose,
            )

        # draw vaf
        plot_kwargs = (
            {'linewidth': 0, 'alpha': 0.7, 'markersize': 1, 'marker': 'o', 'color': 'black'} 
            | plot_kwargs
        )
        gdraw_result = self.draw_dots(
            y_colname=self.get_format_colname(sampleid, 'ALT1_vaf'),
            ax=ax,
            genomeplotter=genomeplotter,

            plotdata=plotdata,
            plot_kwargs=plot_kwargs,

            setup_axes=False,
            subplots_kwargs=subplots_kwargs,
        )

        ax = gdraw_result.ax
        fig = gdraw_result.fig
        genomeplotter = gdraw_result.genomeplotter

        gdraw_result.set_name('vaf')

        # axes setup
        if setup_axes:
            genomedf_draw.draw_axessetup(
                ax=ax,
                genomeplotter=genomeplotter,

                ylabel=ylabel,
                ylabel_prefix=ylabel_prefix,
                ylabel_kwargs=ylabel_kwargs,
                ymax=ymax,
                ymin=ymin,
                yticks=yticks,

                draw_common_kwargs=draw_common_kwargs,
                rotate_chromlabel=rotate_chromlabel,

                fig=fig,
                title=title, 
                suptitle_kwargs=suptitle_kwargs,
            )

        # make legend - intentionally placed after axes setup to use "bbox_to_anchor"
        handles = plotmisc.LegendHandles()
        if vaf_legend:
            self.add_vaf_legend_handle(handles)
        if ccf_legend:
            self.add_ccf_legend_handle(handles)

        ax.legend(handles=handles, loc='upper right', bbox_to_anchor=(1, 1.3))

        return gdraw_result

    def draw_ccf(
        self,
        ax=None,
        genomeplotter=None,
        plotdata=None,

        sampleid=None,
        vaf_legend=False,
        ccf_legend=True,

        # fig generation params
        title=None,
        suptitle_kwargs=dict(),
        subplots_kwargs=dict(),

        # drawing kwargs
        plot_kwargs=dict(),

        # axes setting
        setup_axes=True,
        ylabel=False,
        ylabel_prefix='',
        ylabel_kwargs=dict(),
        ymax=False,
        ymin=False,
        yticks=False,
        draw_common_kwargs=dict(),
        rotate_chromlabel=None,

        # multicore plotdata generation
        nproc=1,
        verbose=True,
    ):
        if ymax is False:
            ymax = 1.1
        if ymin is False:
            ymin = -0.1
        if ylabel is False:
            ylabel = 'Mutation'
        if yticks is False:
            yticks = np.round(np.arange(0, 1.2, 0.2), 1)
        if sampleid is None:
            sampleid = self.samples[0]

        # prepare single plotdata
        if genomeplotter is None:
            genomeplotter = self.get_default_genomeplotter()
        if plotdata is None:
            plotdata = genomeplotter.make_plotdata(
                self, 
                log_suffix=' (mutation ccf)',
                nproc=nproc,
                verbose=verbose,
            )

        # draw
        plot_kwargs = (
            {'linewidth': 0, 'alpha': 0.7, 'markersize': 1, 'color': 'tab:red', 'marker': 'x'} 
            | plot_kwargs
        )
        gdraw_result = self.draw_dots(
            y_colname=self.__class__.CCF_COLNAME,
            ax=ax,
            genomeplotter=genomeplotter,

            plotdata=plotdata,
            plot_kwargs=plot_kwargs,

            setup_axes=False,
            subplots_kwargs=subplots_kwargs,
        )

        ax = gdraw_result.ax
        fig = gdraw_result.fig
        genomeplotter = gdraw_result.genomeplotter

        gdraw_result.set_name('ccf')

        # axes setup
        if setup_axes:
            genomedf_draw.draw_axessetup(
                ax=ax,
                genomeplotter=genomeplotter,

                ylabel=ylabel,
                ylabel_prefix=ylabel_prefix,
                ylabel_kwargs=ylabel_kwargs,
                ymax=ymax,
                ymin=ymin,
                yticks=yticks,

                draw_common_kwargs=draw_common_kwargs,
                rotate_chromlabel=rotate_chromlabel,

                fig=fig,
                title=title, 
                suptitle_kwargs=suptitle_kwargs,
            )

        # make legend - intentionally placed after axes setup to use "bbox_to_anchor"
        handles = plotmisc.LegendHandles()
        if vaf_legend:
            self.add_vaf_legend_handle(handles)
        if ccf_legend:
            self.add_ccf_legend_handle(handles)

        ax.legend(handles=handles, loc='upper right', bbox_to_anchor=(1, 1.3))

        return gdraw_result



#########################################
# parallelized vaf dataframe generation #
#########################################

@deco.get_deco_nproc_limit(3)
@deco.get_deco_atleast1d(['sampleids'])
def get_vafdf(
    vcf_path, 
    sampleids, 
    n_allele=2,
    nproc=1,
    num_split=200,
    exclude_other=False,
    prop=None,
    vpfilter=None,
    verbose=True,
):
    """Columns: MultiIndex with 2 levels
        level 1: (name:                                              None |                                 sid1 |                                 sid2
        level 2: Chromosome, Start, End, REF, [ALT1, [ALT2, ...]]   REF_vaf, [ALT1_vaf, [ALT2_vaf, ...]]   REF_vaf, [ALT1_vaf, [ALT2_vaf, ...]]
    """
    # get VCF fetch regions for each parallel job
    if verbose:
        logutils.log(f'Extracting vcf position information') 

    #num_split = nproc * 10
    fetchregion_gdf_list = vcfmisc.get_vcf_fetchregions_new(
        vcf_path, 
        n=num_split, 
        nproc=nproc,
    )

    # run multiprocess jobs
    args = (
        (
            fetchregion_gdf,
            vcf_path, 
            sampleids, 
            n_allele, 
            exclude_other,
            prop,
            vpfilter,
        )
        for fetchregion_gdf in fetchregion_gdf_list
    )
    with multiprocessing.Pool(nproc) as pool:
        if verbose:
            logutils.log(f'Running parallel jobs') 
        mp_result = pool.starmap(get_vafdf_targetfunc, args)
        if verbose:
            logutils.log(f'Concatenating split job dataframes') 
        result = VCFDataFrame.concat(itertools.chain.from_iterable(mp_result))

    return result


def get_vafdf_targetfunc(
    fetchregion_gdf,
    vcf_path, 
    sampleids, 
    n_allele, 
    exclude_other,
    prop,
    vpfilter,
):
    vafdf_list = list()
    for chrom, start0, end0 in fetchregion_gdf.iter_coords():
        vplist = VariantPlusList.from_vcf_lazy(
            vcf_path, 
            logging_lineno=None, 
            verbose=False,
            init_all_attrs=False,
            vp_init_params=dict(
                init_readstats=True,
                sampleid_list=sampleids,
            ),
        )
        vafdf = vplist.get_vafdf(
            sampleids=sampleids,
            chrom=chrom, start0=start0, end0=end0,
            n_allele=n_allele,
            exclude_other=exclude_other,
            lazy=True,
            prop=prop,
            vpfilter=vpfilter,
        )
        vafdf_list.append(vafdf)

    return vafdf_list



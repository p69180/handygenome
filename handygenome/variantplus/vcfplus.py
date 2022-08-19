import os
import textwrap
import pprint

import pysam

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
workflow = importlib.import_module(".".join([top_package_name, "workflow"]))
variantplus = importlib.import_module(
    ".".join([top_package_name, "variantplus", "variantplus"])
)
varianthandler = importlib.import_module(
    ".".join([top_package_name, "variantplus", "varianthandler"])
)
initvcf = importlib.import_module(".".join([top_package_name, "vcfeditor", "initvcf"]))
breakends = importlib.import_module(
    ".".join([top_package_name, "variantplus", "breakends"])
)
# annotationdb = importlib.import_module('.'.join([top_package_name, 'annotation', 'annotationdb']))
indexing = importlib.import_module(
    ".".join([top_package_name, "vcfeditor", "indexing"])
)
headerhandler = importlib.import_module(
    ".".join([top_package_name, "vcfeditor", "headerhandler"])
)
libpopfreq = importlib.import_module(
    ".".join([top_package_name, "annotation", "popfreq"])
)
libcosmic = importlib.import_module(
    ".".join([top_package_name, "annotation", "cosmic"])
)


LOGGER = workflow.get_logger(name=__name__)
DEFAULT_LOGGING_LINENO = 10000

GREMLIN_HEADER = [
    "CHR1",
    "POS1",
    "CHR2",
    "POS2",
    "CT",
    "SVTYPE",
    "tumor_var_read",
    "tumor_var_split_read",
    "tumor_var_sa_read",
    "normal_var_read",
    "normal_var_split_read",
    "normal_same_clip",
    "tumor_ref_read_bpt1",
    "tumor_ref_read_bpt2",
    "p_value_ref_var_read",
    "tumor_vaf_bpt1",
    "tumor_vaf_bpt2",
    "tumor_var_mapq_bpt1",
    "tumor_var_mapq_bpt2",
    "tumor_depth_change_bpt1",
    "tumor_depth_change_bpt2",
    "normal_other_var_cluster_bpt1",
    "normal_other_var_cluster_bpt2",
    "normal_depth_bpt1",
    "normal_depth_bpt2",
    "normal_panel",
    "score",
    "mate_chr_diff_bp1",
    "mate_chr_same_bp1",
    "mate_chr_diff_prop_bp1",
    "mate_chr_diff_bp2",
    "mate_chr_same_bp2",
    "mate_chr_diff_prop_bp2",
]


class VcfPlus:
    """
    Attributes:
        vcf_path
        vcf
        refver
        fasta
        chromdict
        vplist
    """

    def __init__(
        self,
        vcf_path,
        fasta=None,
        chromdict=None,
        init_vp=True,
        logger=None,
        init_popfreq=True,
        init_cosmic=True,
        init_transcript=True,
        init_regulatory=True,
        init_motif=False,
        init_repeat=True,
        init_readstats=True,
        logging_lineno=DEFAULT_LOGGING_LINENO,
    ):
        """Args:
            fasta : pysam.FastaFile object
            vcf_path : Path to a vcf file
        """

        def init_vcf(vcf_path):
            self.vcf = pysam.VariantFile(self.vcf_path, "r")
            # annotationdb.add_infometas(self.vcf.header)

        def set_fasta_chromdict(fasta, chromdict):
            if fasta is None:
                self.fasta = pysam.FastaFile(common.DEFAULT_FASTA_PATHS[self.refver])
            else:
                self.fasta = fasta

            if chromdict is None:
                self.chromdict = common.ChromDict(fasta=self.fasta)
            else:
                self.chromdict = chromdict

#        def init_vp_containers(vplist, set_readstats, logging_lineno):
#            self.set_vplist(
#                vplist=vplist,
#                #set_annotdb=set_annotdb,
#                set_readstats=set_readstats,
#                logging_lineno=logging_lineno,
#            )
#            self.vplist_filtered = self.vplist
#            # if set_details:
#            # self.set_ID_attributes()
#            # self.set_more_vp_containers()

        self._vplist_show_len = 30

        self.vcf_path = vcf_path
        self.logger = LOGGER if (logger is None) else logger

        init_vcf(vcf_path)
        self.refver = common.infer_refver_vcfheader(self.vcf.header)
        set_fasta_chromdict(fasta, chromdict)

        # initiation of vp containers 
        if init_vp:
            self.set_vplist(
                vplist=None,
                init_popfreq=init_popfreq,
                init_cosmic=init_cosmic,
                init_transcript=init_transcript,
                init_regulatory=init_regulatory,
                init_motif=init_motif,
                init_repeat=init_repeat,
                init_readstats=init_readstats,
                logging_lineno=logging_lineno,
            )
            self.vplist_filtered = self.vplist
            self.length = len(self.vplist)

    def __repr__(self):
        result = list()
        result.extend([f"<VcfPlus object>", f"- refver: {self.refver}"])

        result.extend(
            [
                f"- vplist:",
                textwrap.indent(
                    pprint.pformat(self.vplist[: self._vplist_show_len]), " " * 4
                ),
            ]
        )
        if len(self.vplist) > self._vplist_show_len:
            result.append(" " * 4 + "...")
        result.append(" " * 4 + f"{len(self.vplist)} VariantPlus objects")

        if self.vplist_filtered is not None:
            result.extend(
                [
                    f"- vplist_filtered:",
                    textwrap.indent(
                        pprint.pformat(self.vplist_filtered[: self._vplist_show_len]),
                        " " * 4,
                    ),
                ]
            )
            if len(self.vplist_filtered) > self._vplist_show_len:
                result.append(" " * 4 + "...")
            result.append(" " * 4 + f"{len(self.vplist_filtered)} VariantPlus objects")

        return "\n".join(result)

    ###########################################################

    def set_header(self, header=None):
        if header is None:
            if self.vcf is None:
                self.header = None
            else:
                self.header = self.vcf.header
        else:
            self.header = header

    def _update_contig(self):
        if not check_having_good_contigs(self.vcf):
            self.vcf = get_updated_pysamvcf(self.vcf, chromdict=self.chromdict)

    def _set_has_good_contigs(self):
        self.has_good_contigs = check_having_good_contigs(self.vcf)

    ###########################################################

    def set_vplist(
        self,
        vplist=None,
        init_popfreq=True,
        init_cosmic=True,
        init_transcript=True,
        init_regulatory=True,
        init_motif=False,
        init_repeat=True,
        init_readstats=True,
        logging_lineno=DEFAULT_LOGGING_LINENO,
    ):
        """SV variant records for bnd2 are not loaded"""

        # construct vplist
        if vplist is None:
            self.vplist = variantplus.VariantPlusList()

            if init_popfreq:
                popfreq_metadata = libpopfreq.PopfreqMetadata.from_vcfheader(self.vcf.header)
            else:
                popfreq_metadata = None

            if init_cosmic:
                cosmic_metadata = libcosmic.CosmicMetadata.from_vcfheader(self.vcf.header)
            else:
                cosmic_metadata = None

            for vr in common.iter_lineno_logging(self.vcf.fetch(), self.logger, logging_lineno):
                if varianthandler.check_SV(vr):
                    vr_svinfo = breakends.get_vr_svinfo_standard_vr(
                        vr, self.fasta, self.chromdict
                    )
                    if not vr_svinfo["is_bnd1"]:
                        continue

                vp = variantplus.VariantPlus.from_vr(
                    vr=vr,
                    refver=self.refver,
                    fasta=self.fasta,
                    chromdict=self.chromdict,
                    init_popfreq=init_popfreq,
                    init_cosmic=init_cosmic,
                    init_transcript=init_transcript,
                    init_regulatory=init_regulatory,
                    init_motif=init_motif,
                    init_repeat=init_repeat,
                    init_readstats=init_readstats,
                    popfreq_metadata=popfreq_metadata, 
                    cosmic_metadata=cosmic_metadata,
                )
                self.vplist.append(vp)
        else:
            self.vplist = vplist

        # set gr
        self.vplist.set_gr()

    def sort_vplist(self):
        self.logger.info("Sorting self.vplist")
        self.vplist.sort(key=variantplus.get_vp_sortkey(self.chromdict))

    def write(
        self, outfile_path, vplist=None, mode_bcftools="z", mode_pysam=None, index=True
    ):
        mode_pysam = common.write_mode_arghandler(mode_bcftools, mode_pysam)

        if vplist is None:
            self.sort_vplist()
            vplist = self.vplist

        with pysam.VariantFile(
            outfile_path, mode=mode_pysam, header=self.vcf.header
        ) as out_vcf:
            for vp in vplist:
                if vp.is_sv:
                    out_vcf.write(vp.vr)
                    out_vcf.write(vp.get_vr_bnd2())
                else:
                    out_vcf.write(vp.vr)

        if index:
            indexing.index_vcf(outfile_path)

    ###########################################################

    def filter_isec(self, gr):
        self.vplist_filtered = self.vplist.get_isec(gr)

    def filter_vplist(self, vpfilter):
        new_vplist_filtered = variantplus.VariantPlusList()
        new_vplist_filtered.extend(filter(vpfilter, self.vplist))
        self.vplist_filtered = new_vplist_filtered

    ###########################################################

    #    def set_ID_attributes(self):
    #        self.set_has_duplicate_pos_alleles()
    #        self.set_has_id()
    #        self.set_has_duplicate_id()
    #        self.set_has_mateid()
    #        self._sanity_check_dupID()
    #
    #    def set_has_duplicate_pos_alleles(self):
    #        self.has_duplicate_pos_alleles = \
    #                varianthandler.check_has_duplicate_pos_alleles(vp.vr for vp in self.vplist)
    #
    #    def set_has_id(self):
    #        self.has_id = varianthandler.check_has_id(vp.vr for vp in self.vplist)
    #
    #    def set_has_duplicate_id(self):
    #        if self.has_id:
    #            self.has_duplicate_id = varianthandler.check_has_duplicate_id(vp.vr for vp in self.vplist)
    #        else:
    #            self.has_duplicate_id = False
    #
    #    def set_has_mateid(self):
    #        if self.has_id:
    #            self.has_mateid = varianthandler.check_has_mateid(vp.vr for vp in self.vplist)
    #        else:
    #            self.has_mateid = False
    #
    #    def _sanity_check_dupID(self):
    #        if self.has_duplicate_id:
    #            raise Exception(f'VCF file {self.vcf_path} has duplicate IDs.')

    ###########################################################

    #    def set_more_vp_containers(self):
    #        self.set_id_dict()
    #        self.set_mateid_dict()
    #        self._sanity_check_mateid()
    #        self.set_vpp_list()
    #
    #    def set_id_dict(self):
    #        self.id_dict = dict()
    #        for vp in self.vplist:
    #            if vp.vr.id is not None:
    #                self.id_dict[vp.vr.id] = vp
    #
    #    def set_mateid_dict(self):
    #        self.mateid_dict = dict()
    #        for vp in self.vplist:
    #            if vp.vr.id is not None:
    #                if vp.check_NA_info('MATEID'):
    #                    self.mateid_dict[vp.vr.id] = None
    #                else:
    #                    self.mateid_dict[vp.vr.id] = vp.get_value_info('MATEID')
    #
    #    def _sanity_check_mateid(self):
    #        all_IDs = set(self.id_dict.keys())
    #        for ID, MATEID in self.mateid_dict.items():
    #            if MATEID is not None:
    #                if MATEID in all_IDs:
    #                    mate_of_mate = self.mateid_dict[MATEID]
    #                    if mate_of_mate != ID:
    #                        raise Exception(f'MATEID validity error. {ID} -> {MATEID} but {MATEID} -> {mate_of_mate}')
    #                else:
    #                    raise Exception(f'MATEID validity error. "{MATEID}", which is the mate of "{ID}", is not in the vcf.')
    #
    #    def set_vpp_list(self):
    #        def get_vpp_from_one_SVvp(vp):
    #            new_vp1 = vp.make_maxspan_vp(pos1 = True)
    #            new_vp2 = vp.make_maxspan_vp(pos2 = True)
    #            return variantpluspair.VariantPlusPair([new_vp1, new_vp2])
    #
    #        self.vpp_list = list()
    #        checked_IDs = set()
    #
    #        for vp in self.vplist:
    #            if vp.vr.id in checked_IDs:
    #                continue
    #            else:
    #                if vp.check_NA_info('MATEID'):
    #                    if varianthandler.check_SV(vp.vr):
    #                        vpp = get_vpp_from_one_SVvp(vp)
    #                    else:
    #                        vpp = variantpluspair.VariantPlusPair([vp])
    #
    #                    self.vpp_list.append(vpp)
    #                    checked_IDs.add(vp.vr.id)
    #                else:
    #                    mateid = self.mateid_dict[vp.vr.id]
    #                    mate_vp = self.id_dict[mateid]
    #
    #                    if vp.bnds.sameseq(mate_vp.bnds):
    #                        vpp = variantpluspair.VariantPlusPair([vp, mate_vp])
    #                        self.vpp_list.append(vpp)
    #                    else:
    #                        vpp1 = get_vpp_from_one_SVvp(vp)
    #                        self.vpp_list.append(vpp1)
    #                        vpp2 = get_vpp_from_one_SVvp(mate_vp)
    #                        self.vpp_list.append(vpp2)
    #
    #                    checked_IDs.add(vp.vr.id)
    #                    checked_IDs.add(mateid)

    # self.vp_pair_list_sorted = self.check_vp_pair_list_is_sorted()

    ###########################################################

    #    def check_vplist_is_sorted(self):
    #        for idx in range(len(self.vplist) - 1):
    #            vp_pre = self.vplist[idx]
    #            vp_post = self.vplist[idx + 1]
    #            if common.compare_coords(
    #                    vp_pre.vcfspec.chrom, vp_pre.vcfspec.pos,
    #                    vp_post.vcfspec.chrom, vp_post.vcfspec.pos,
    #                    self.chromdict) > 0:
    #                return False
    #
    #        return True
    #
    #    def check_vp_pair_list_is_sorted(self):
    #        for idx in range(len(self.vp_pair_list) - 1):
    #            vp_pair_pre = self.vp_pair_list[idx]
    #            vp_pair_post = self.vp_pair_list[idx + 1]
    #            if common.compare_coords(
    #                    vp_pair_pre['vp1'].vr.contig,
    #                    vp_pair_pre['vp1'].vr.pos,
    #                    vp_pair_post['vp1'].vr.contig,
    #                    vp_pair_post['vp1'].vr.pos,
    #                    self.chromdict) > 0:
    #                return False
    #        return True
    #
    #    def sort_vp_pair_list(self):
    #        """
    #        Sort by coordinate of breakend 1
    #        """
    #        self.vp_pair_list.sort(
    #            key = lambda x: varianthandler.vr_sortkey(x[0].vr,
    #                                                           self.chromdict))

    ###########################################################

    def write_gremlin(self, outfile_path):
        NAchar = "."
        with open(outfile_path, "w") as outfile:
            outfile.write("#" + "\t".join(GREMLIN_HEADER) + "\n")
            for pair in self.vp_pair_list:
                vp = pair["vp1"]
                line = list()

                line.append(vp.bnds.chrom_bnd1)  # CHR1
                line.append(str(vp.bnds.pos_bnd1))  # POS1
                line.append(vp.bnds.chrom_bnd2)  # CHR1
                line.append(str(vp.bnds.pos_bnd2))  # POS1

                CT1 = "5" if vp.bnds.endis5_bnd1 else "3"
                CT2 = "5" if vp.bnds.endis5_bnd2 else "3"
                CT = CT1 + "to" + CT2
                line.append(CT)  # CT

                line.append(vp.bnds.SVTYPE)  # SVTYPE

                for field in GREMLIN_HEADER[6:]:
                    if field in vp.vr.info.keys():
                        val = str(vp.vr.info[field])
                    else:
                        val = NAchar
                    line.append(val)

                outfile.write("\t".join(line) + "\n")

    ###########################################################


def check_having_good_contigs(vcf):
    """
    Checks if CHROM names of all variant records are included in metadata
        contig names.
    Returns True if so, or returns False.

    Args:
        vcf: pysam.VariantFile instance
    """

    good = True

    chromdict = common.ChromDict(pysamhdr=vcf.header)
    for vr in vcf.fetch():
        if vr.contig in chromdict.contigs:
            if vr.pos in range(1, chromdict[vr.contig] + 1):
                continue
            else:
                good = False
                break
        else:
            good = False
            break

    return good


######################################


def get_updated_pysamvcf(
    vcf, chromdict=None, samples=None, pysamhdr=None, vcfver=common.DEFAULT_VCFVER
):
    """
    Writes a new temporary vcf file (removed at the end) with the same
        contents as input pysam.VariantFile except contig headers.
    Returns a pysam.VariantFile object generated with the newly written
        vcf file.
    """

    tf_path = workflow.get_tmpfile_path(delete=False)
    write_updated_vcf(tf_path, vcf, chromdict, samples, pysamhdr, vcfver)
    new_vcf = pysam.VariantFile(tf_path, "r")
    os.remove(tf_path)

    return new_vcf


def write_updated_vcf(
    outfile_path,
    vcf,
    chromdict=None,
    samples=None,
    pysamhdr=None,
    vcfver=common.DEFAULT_VCFVER,
    mode_bcftools="z",
    mode_pysam=None,
):
    mode_pysam = common.write_mode_arghandler(mode_bcftools, mode_pysam)
    if pysamhdr is None:
        header = initvcf.create_header(chromdict, samples, vcf.header, vcfver)
    else:
        header = initvcf.create_header(chromdict, samples, pysamhdr, vcfver)

    with pysam.VariantFile(outfile_path, mode=mode_pysam, header=header) as out_vcf:
        for vr in vcf.fetch():
            out_vcf.write(varianthandler.reform_samples(vr, pysamhdr=header))

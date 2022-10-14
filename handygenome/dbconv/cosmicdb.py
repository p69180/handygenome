import gzip
import collections
import functools
import operator
import itertools
import re

import pysam

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))
workflow = importlib.import_module(".".join([top_package_name, "workflow"]))
hgvs = importlib.import_module(".".join([top_package_name, "hgvs"]))
varianthandler = importlib.import_module(".".join([top_package_name, "variant", "varianthandler"]))
infoformat = importlib.import_module(".".join([top_package_name, "variant", "infoformat"]))
annotitem = importlib.import_module(".".join([top_package_name, "annotation", "annotitem"]))
libcosmic = importlib.import_module(".".join([top_package_name, "annotation", "cosmic"]))
indexing = importlib.import_module(".".join([top_package_name, "vcfeditor", "indexing"]))


LOGGER = workflow.get_logger(name='COSMIC converter')


def get_header(infile_path):
    with gzip.open(infile_path, "rt") as infile:
        header = common.get_linesp(next(infile), sep="\t")
    return header


def check_decodable(b):
    try:
        b.decode("utf-8")
    except UnicodeDecodeError:
        return False
    else:
        return True


def get_linedict_byte(byteline, header):
    return dict(zip(
                header, common.get_linesp_byte(byteline)
                ))


def get_infile_looper(infile_path):
    def decoding_error_handler(line_bytes, header, NR):
        try:
            linedict_bytes = get_linedict_byte(line_bytes, header)
            if not check_decodable(bytes([linedict_bytes["Mutation CDS"][-1]])):
                linedict_bytes["Mutation CDS"] = linedict_bytes["Mutation CDS"][:-1]
                # cosmic v95 GRCh38 CosmicMutantExport.tsv.gz file contains 
                # a line where 'Mutation CDS' field value contains a non-ascii byte
            linedict_str = {
                key: val.decode("utf-8") for key, val in linedict_bytes.items()
            }
            return linedict_str
        except:
            print("Unexpected decoding error")
            print(f"{NR:,}th line")
            print(line_bytes)
            raise

    header = get_header(infile_path)
    with gzip.open(infile_path, "rb") as infile:
        _ = next(infile)  # skip the header line
        NR = 1
        for line_bytes in infile:
            NR += 1
            try:
                line_str = line_bytes.decode("utf-8")
            except UnicodeDecodeError:
                linedict_str = decoding_error_handler(line_bytes, header, NR)
            else:
                linedict_str = dict(zip(header, common.get_linesp(line_str)))

            yield linedict_str


def union_all(iter_of_iter):
    return set(itertools.chain.from_iterable(iter_of_iter))


################################
# file type specific functions #
################################

# line_parser
def line_parser_mutantexport(linedict):
    coding_score_key = 'FATHMM score'

    result = dict()

    result['cosv'] = linedict['GENOMIC_MUTATION_ID']
    result['sample_id'] = int(linedict['ID_sample'])
    result['tumor_id'] = int(linedict['ID_tumour'])
    result['site'] = linedict['Primary site']
    result['hgvsg'] = linedict['HGVSG']
    result['hgvsc'] = linedict['HGVSC']
    result['somatic_status'] = linedict['Mutation somatic status']

    if coding_score_key in linedict:
        result['coding_score'] = (None 
                                  if (linedict[coding_score_key] == '') else 
                                  float(linedict[coding_score_key]))
    else:
        result['coding_score'] = None

    result['noncoding_score'] = None

    return result


def line_parser_ncv(linedict):
    coding_score_key = 'FATHMM_MKL_CODING_SCORE'
    noncoding_score_key = 'FATHMM_MKL_NON_CODING_SCORE'

    result = dict()

    result['cosv'] = linedict['GENOMIC_MUTATION_ID']
    result['sample_id'] = int(linedict['ID_SAMPLE'])
    result['tumor_id'] = int(linedict['ID_tumour'])
    result['site'] = linedict['Primary site']
    result['hgvsg'] = linedict['HGVSG']
    result['hgvsc'] = None
    result['somatic_status'] = linedict['Mutation somatic status']

    if coding_score_key in linedict:
        result['coding_score'] = (
                None 
                if (linedict[coding_score_key] == '') else 
                float(linedict[coding_score_key])
                )
    else:
        result['coding_score'] = None

    if noncoding_score_key in linedict:
        result['noncoding_score'] = (
                None 
                if (linedict[noncoding_score_key] == '') else 
                float(linedict[noncoding_score_key])
                )
    else:
        result['noncoding_score'] = None

    return result


################################


def line_filter(parsed_line):
    """exclude if False"""
    if parsed_line['cosv'] == '':
        return False
    elif parsed_line['hgvsg'] == '':
        return False
    else:
        return True


def line_sanity_check(parsed_line):
    pass


def init_cosv_info():
    return {
        #'sampleID': list(),
        'primary_site': dict(),
        'primary_site_somatic': dict(),
        'hgvsg': set(),
        'hgvsc': set(),
        'coding_score': set(),
        'noncoding_score': set(),
    }


def load_datafile(infile_path, line_parser, summary, site_count, logging_lineno=1_000_000):
    """Values of 'Mutation somatic status':
        'Confirmed somatic variant',
        'Not specified',
        'Reported in another cancer sample as somatic',
        'Systematic screen',
        'Variant of unknown origin'
    """
    def update_site_count(site_count, parsed_line, uid):
        site_count.setdefault(parsed_line['site'], set())
        site_count[parsed_line['site']].add(uid)

    def handle_sites(cosv_info, parsed_line, uid):
        cosv_info['primary_site'].setdefault(parsed_line['site'], set())
        cosv_info['primary_site'][parsed_line['site']].add(uid)

        if parsed_line['somatic_status'] == 'Confirmed somatic variant':
            cosv_info['primary_site_somatic'].setdefault(parsed_line['site'], set())
            cosv_info['primary_site_somatic'][parsed_line['site']].add(uid)

    def handle_hgvs(cosv_info, parsed_line):
        cosv_info['hgvsg'].add(parsed_line['hgvsg'])
        if parsed_line['hgvsc'] is not None:
            cosv_info['hgvsc'].add(parsed_line['hgvsc'])

    def handle_scores(cosv_info, parsed_line):
        if parsed_line['coding_score'] is not None:
            cosv_info['coding_score'].add(parsed_line['coding_score'])
        if parsed_line['noncoding_score'] is not None:
            cosv_info['noncoding_score'].add(parsed_line['noncoding_score'])

    # loop over data file
    looper = get_infile_looper(infile_path)
    for linedict in common.iter_lineno_logging(looper, LOGGER, logging_lineno):
        line_sanity_check(linedict)
        parsed_line = line_parser(linedict)
        if not line_filter(parsed_line):
            continue

        cosv = parsed_line['cosv']
        uid = (parsed_line['sample_id'], parsed_line['tumor_id'])

        update_site_count(site_count, parsed_line, uid)

        if cosv not in summary:
            summary[cosv] = init_cosv_info()
        cosv_info = summary[cosv]

        handle_sites(cosv_info, parsed_line, uid)
        handle_hgvs(cosv_info, parsed_line)
        handle_scores(cosv_info, parsed_line)

    # postprocess
    #total_count = sum(site_count.values())
    #site_count['total'] = total_count

    #print(f'The number of COSV IDs: {len(summary)}')  # 15689770
    #print(f'The number of COSV IDs with at least one confirmed somatic: {len(summary["primary_sites_somatic"])}')
        # 14577671


###################################################


class MultipleScoreError(Exception):
    pass


def modify_scores(cosv_info, key):
    if len(cosv_info[key]) == 0:
        cosv_info[key] = None
    elif len(cosv_info[key]) == 1:
        cosv_info[key] = cosv_info[key].pop()
    else:
        raise MultipleScoreError(f'More than one values for this cosv entry: {cosv_info}')


def modify_sitedic(site_dic):
    keys = tuple(site_dic.keys())
    for site in keys:
        site_dic[site] = len(site_dic[site])


def modify_vcfspec(cosv_info, fasta, is_hg19):
    def hgvsg_to_vcfspecs(hgvsg_set, fasta):
        vcfspec_candidates = set()
        for hgvsg in hgvsg_set:
            vcfspec = hgvs.hgvsg_to_vcfspec(hgvsg, fasta, leftmost=True)
            #if vcfspec.check_without_N():
            vcfspec_candidates.add(vcfspec)

        return vcfspec_candidates

    def hgvsc_to_vcfspecs(hgvsc_set, fasta, is_hg19):
        vcfspec_candidates = set()
        for hgvsc in hgvsc_set:
            # Running ensembl variantrecoder with hgvsc sometimes results 
                # in error
                # e.g. (v95 GRCh37 CosmicMutantExport.tsv.gz COSV104545830) 
                # Unable to parse HGVS notation 'ENST00000424344.3:c.-166N>A':  
                # : Reference allele extracted from 
                # ENST00000424344:155227447-155227447 (G) does not match 
                # reference allele given by HGVS notation 
                # ENST00000424344.3:c.-166N>A (N)
            try:
                vcfspec = hgvs.hgvsc_to_vcfspec(hgvsc, is_hg19, fasta, leftmost=True)
            except:
                pass
            else:
                #if vcfspec.check_without_N():
                vcfspec_candidates.add(vcfspec)

        return vcfspec_candidates

    def get_vcfspec(cosv_info, fasta, is_hg19):
        """Invalid COSV ID if linked to multiple vcfspecs"""
        vcfspec_candidates = hgvsg_to_vcfspecs(cosv_info['hgvsg'], fasta)
        if len(vcfspec_candidates) > 1:
            vcfspec = None
        elif len(vcfspec_candidates) == 1:
            vcfspec = vcfspec_candidates.pop()
        elif len(vcfspec_candidates) == 0:
            new_vcfspec_candidates = hgvsc_to_vcfspecs(cosv_info['hgvsc'], fasta, is_hg19)
            if len(new_vcfspec_candidates) == 1:
                vcfspec = new_vcfspec_candidates.pop()
            else:
                vcfspec = None

        return vcfspec

    # main
    cosv_info['vcfspec'] = get_vcfspec(cosv_info, fasta, is_hg19)


@common.get_deco_timestamp('SUMMARY MODIFICATION', LOGGER)
def modify_summary(summary, site_count, fasta, is_hg19, logging_lineno=500_000):
    # summary
    for cosv, cosv_info in common.iter_lineno_logging(summary.items(), LOGGER, logging_lineno):
        # coding & noncoding score
        for key in ('coding_score', 'noncoding_score'):
            try:
                modify_scores(cosv_info, key)
            except MultipleScoreError as e:
                raise Exception(f'Error while modifying "{key}"') from e

        # primary site
        modify_sitedic(cosv_info['primary_site'])
        modify_sitedic(cosv_info['primary_site_somatic'])

        # vcfspec
        modify_vcfspec(cosv_info, fasta, is_hg19)

    # site_count
    modify_sitedic(site_count)


###############################################


@common.get_deco_timestamp('VCFSPEC-TO-COSV MAPPING', LOGGER)
def get_vcfspec_cosv_map(summary, fasta, is_hg19):
    # main
    vcfspec_cosv_map = dict()
    for cosv, cosv_info in summary.items():
        vcfspec = cosv_info['vcfspec']
        if vcfspec is not None:
            vcfspec_cosv_map.setdefault(vcfspec, set())
            vcfspec_cosv_map[vcfspec].add(cosv)

    # deduplicate
    dedup_vcfspec_cosv_map = dict()
    for vcfspec in vcfspec_cosv_map:
        dedup_vcfspec_cosv_map[vcfspec] = vcfspec_cosv_map[vcfspec].pop()

    return dedup_vcfspec_cosv_map


@common.get_deco_timestamp('VCFSPEC SORTING', LOGGER)
def sort_vcfspecs(vcfspec_cosv_map, chromdict):
    return sorted(vcfspec_cosv_map.keys(), key=common.get_vcfspec_sortkey(chromdict))


###############################################


def get_cosmic_metadata(site_count, refver, cosmic_version):
    cosmicmeta = libcosmic.CosmicMetadata()
    cosmicmeta['num_sample_by_site'] = site_count
    cosmicmeta['reference_version'] = refver
    cosmicmeta['cosmic_version'] = cosmic_version

    return cosmicmeta


def get_vcf_header(chromdict, cosmic_metadata):
    header = pysam.VariantHeader()

    for contig, length in chromdict.items():
        if re.fullmatch('(chr)?([0-9]+|X|Y)|chrM|MT', contig) is not None:
            header.contigs.add(contig, length)

    cosmic_metadata.write_header(header)
    libcosmic.CosmicInfoALTlist.add_meta_info(header)
    
    return header


def into_CosmicInfo(cosv, cosv_info):
    cosmicinfo = libcosmic.CosmicInfo()
    cosmicinfo['id'] = cosv
    cosmicinfo['occurrence'] = cosv_info['primary_site']
    cosmicinfo['occurrence_somatic'] = cosv_info['primary_site_somatic']
    cosmicinfo['coding_score'] = cosv_info['coding_score']
    cosmicinfo['noncoding_score'] = cosv_info['noncoding_score']

    cosmicinfolist = libcosmic.CosmicInfoALTlist()
    cosmicinfolist.append(cosmicinfo)

    return cosmicinfolist


@common.get_deco_timestamp('WRITING FINAL OUTPUT FILE', LOGGER)
def write_outfile(outfile_path, summary, site_count, vcfspec_cosv_map, sorted_vcfspecs, chromdict, refver, cosmic_version, logging_lineno=1_000_000):
    cosmic_metadata = get_cosmic_metadata(site_count, refver, cosmic_version)
    header = get_vcf_header(chromdict, cosmic_metadata)
    with pysam.VariantFile(outfile_path, mode='wz', header=header) as out_vcf:
        for vcfspec in common.iter_lineno_logging(sorted_vcfspecs, LOGGER, logging_lineno):
            cosv = vcfspec_cosv_map[vcfspec]
            cosv_info = summary[cosv]
            cosmicinfolist = into_CosmicInfo(cosv, cosv_info)

            vr = header.new_record()
            varianthandler.apply_vcfspec(vr, vcfspec)
            cosmicinfolist.write_info(vr)

            out_vcf.write(vr)

    indexing.index_vcf(outfile_path)


###############################################


@common.get_deco_arg_choices({'refver': ('GRCh37', 'GRCh38')})
def main(infile_path_ncv, infile_path_mutantexport, outfile_path, refver, cosmic_version):
    fasta = pysam.FastaFile(common.DEFAULT_FASTA_PATHS[refver])
    chromdict = common.ChromDict(refver=refver)
    is_hg19 = (refver == 'GRCh37')

    summary = dict()
    site_count = dict()

    LOGGER.info(f'Beginning loading of CosmicMutantExport file')
    load_datafile(infile_path_mutantexport, line_parser_mutantexport, summary, site_count)
    LOGGER.info(f'Finished')

    LOGGER.info(f'Beginning loading of CosmicNCV file')
    load_datafile(infile_path_ncv, line_parser_ncv, summary, site_count)
    LOGGER.info(f'Finished')

    modify_summary(summary, site_count, fasta, is_hg19)
    vcfspec_cosv_map = get_vcfspec_cosv_map(summary, fasta, is_hg19)
    sorted_vcfspecs = sort_vcfspecs(vcfspec_cosv_map, chromdict)

    write_outfile(outfile_path, summary, site_count, vcfspec_cosv_map, sorted_vcfspecs, chromdict, refver, cosmic_version)



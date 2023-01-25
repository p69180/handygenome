import re
import textwrap

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
infoformat = importlib.import_module('.'.join([top_package_name, 'variant', 'infoformat']))
headerhandler = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'headerhandler']))
initvcf = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'initvcf']))
libvcfspec = importlib.import_module('.'.join([top_package_name, 'variant', 'vcfspec']))

#def get_vcfspec(vr):
#    return common.Vcfspec(vr.contig, vr.pos, vr.ref, tuple(vr.alts))


def check_SV(vr):
    alt = vr.alts[0]
    if any(
            (re.fullmatch(f'<{x}(:.+)?>', alt) is not None)
            for x in libvcfspec.SV_ALTS):
        # <DEL>, <INV>, ...
        return True
    elif (
        (common.RE_PATS['alt_bndstring_1'].fullmatch(alt) is not None) or
        (common.RE_PATS['alt_bndstring_2'].fullmatch(alt) is not None)):
        return True
    else:
        return False


def check_cpgmet(vr):
    return vr.alts[0] == f'<{libvcfspec.CPGMET_ALT}>'


# sorting helpers

def vr_sortkey(vr, chromdict):
    return common.coord_sortkey(vr.contig, vr.pos, chromdict)


# applying vcfspec to vr

#def apply_vcfspec(vr, vcfspec):
#    vr.contig = vcfspec.chrom
#    vr.pos = vcfspec.pos
#    vr.ref = vcfspec.ref
#    vr.alts = vcfspec.alts
#    if not vcfspec.components.is_missing:
#        vcfspec.components.write(vr)


# vr reformatting

def merge(vr_list, vcfheader=None):
    """Intended to be run with variant records with identical vcfspecs.
    chrom, pos, id, ref, alt, qual, INFO/END is derived from the first vr.
    FILTER values are union'ed except PASS value.
    INFO and FORMAT values are chosen from the first variant record which
        contains the key.
    """

    if vcfheader is None:
        vcfheader = headerhandler.merge_vcfheaders([x.header for x in vr_list])

    new_vr = vcfheader.new_record()
    _copy_basic_attrs(vr_list[0], new_vr)

    for idx, vr in enumerate(vr_list):
        _copy_filter(vr, new_vr, skip_pass=True)
        _copy_info(vr, new_vr, dont_copy_if_exist=True)
        _copy_format(vr, new_vr, dont_copy_if_exist=True)

    # handle CallerInfo
    from handygenome.annotation.callerinfo import CallerInfo
    import handygenome.annotation.callerinfo as libcallerinfo

    if any(CallerInfo.get_annotkey() in vr.header.info.keys() for vr in vr_list):
        merged_callerinfo = libcallerinfo.get_merged_callerinfo(vr_list)
        merged_callerinfo.write(new_vr)

    return new_vr


def reheader(vr, vcfheader):
    """Return a new variant record created using the template
    header, filled with the contents of the input variant record.

    INFO or FORMAT keys or samples not included in the template header 
    are discarded.
    """

    new_vr = vcfheader.new_record()
    _copy_basic_attrs(vr, new_vr)
    _copy_filter(vr, new_vr, skip_pass=False)
    _copy_info(vr, new_vr, dont_copy_if_exist=False)
    _copy_format(vr, new_vr, dont_copy_if_exist=False)

    # amend missing values in new_vr which are
        # erroneously represented by .................
    new_vr_format_keys = tuple(new_vr.header.formats.keys())
    new_vr_samples = tuple(new_vr.samples.keys())
    for sampleid in new_vr_samples:
        for key in new_vr_format_keys:
            if new_vr.header.formats[key].number == 'G':
                continue
            if infoformat.check_NA_format(vr, sampleid, key):
                infoformat.set_NA_format(new_vr, sampleid, key)

    return new_vr


def rename(vr, samples):
    if len(samples) != len(vr.header.samples):
        raise Exception(
            f'The number of samples must be the same.')

    chromdict = common.ChromDict(vcfheader=vr.header)
    new_header = initvcf.create_header(chromdict=chromdict, samples=samples, vcfheader=vr.header)
    new_vr = new_header.new_record()
    _copy_basic_attrs(vr, new_vr)
    _copy_filter(vr, new_vr, skip_pass=False)
    _copy_info(vr, new_vr, dont_copy_if_exist=False)

    for old_sampleid, new_sampleid in zip(vr.header.samples, samples):
        for key in vr.samples[old_sampleid]:
            if key == 'GT':
                new_vr.samples[new_sampleid].allele_indices = vr.samples[old_sampleid].allele_indices
                new_vr.samples[new_sampleid].phased = vr.samples[old_sampleid].phased
            else:
                infoformat.set_format(
                    vr=new_vr, 
                    sampleid=new_sampleid,
                    key=key, 
                    val=infoformat.get_format(vr, old_sampleid, key),
                )

    return new_vr


# vr reformatting helpers

def _copy_basic_attrs(old_vr, new_vr):
    for key in ('contig', 'pos', 'id', 'qual'):
        setattr(new_vr, key, getattr(old_vr, key))
    if old_vr.alleles is not None:
        # assigning None to vr.alleles raises an Exception
        new_vr.alleles = old_vr.alleles
    if old_vr.stop is not None:
        new_vr.stop = old_vr.stop


def _copy_filter(old_vr, new_vr, skip_pass=True):
    for filter_name in old_vr.filter:
        # check if new vr header defines the key
        if filter_name in new_vr.header.filters:
            if skip_pass:
                # copy filter only if it is not PASS
                if filter_name != 'PASS':
                    new_vr.filter.add(filter_name)
            else:
                new_vr.filter.add(filter_name)


def _copy_info(old_vr, new_vr, dont_copy_if_exist=True):
    for info_key in old_vr.info:
        # check if new vr header defines the key
        if info_key in new_vr.header.info:  
            # skip if Number == 'G'
            if new_vr.header.info[info_key].number == 'G':
                continue
            if dont_copy_if_exist:
                # copy the value only if the value is not already set
                if info_key not in new_vr.info:
                    infoformat.set_info(new_vr, info_key, 
                                        infoformat.get_info(old_vr, info_key))
            else:
                infoformat.set_info(new_vr, info_key, 
                                    infoformat.get_info(old_vr, info_key))


def _copy_format(old_vr, new_vr, dont_copy_if_exist=True):
    def assign_val(old_vr, new_vr, key, sampleid):
        if key == 'GT':
            new_vr.samples[sampleid].allele_indices = old_vr.samples[sampleid].allele_indices
            new_vr.samples[sampleid].phased = old_vr.samples[sampleid].phased
        else:
            infoformat.set_format(new_vr, sampleid, key, infoformat.get_format(old_vr, sampleid, key))

    new_vr_samples = tuple(new_vr.samples)
    for sampleid in old_vr.samples:
        # check if the new vr header defines the sample id of input vr
        if sampleid in new_vr_samples:
            for key in old_vr.samples[sampleid]:
                # check if new vr header defines the key
                if key in new_vr.header.formats:
                    # skip if Number == 'G'
                    if new_vr.header.formats[key].number == 'G':
                        continue

                    if dont_copy_if_exist:
                        # copy the value only if the value is not already set
                        if key not in new_vr.samples[sampleid]:
                            assign_val(old_vr, new_vr, key, sampleid)
                    else:
                        assign_val(old_vr, new_vr, key, sampleid)


# raw VCF string parsers
def parse_raw_info(vr):
    result = dict()
    infosp = str(vr).split('\t')[7].split(';')
    for x in infosp:
        xsp = x.split('=')
        if len(xsp) == 1:
            result[xsp[0]] = vr.info[xsp[0]]
        else:
            result[xsp[0]] = xsp[1]
    return result


# functions dealing with list of vr
    
def check_has_id(vr_list):
    """
    Args:
        vr_list: A list of pysam.VariantRecord objects
    """

    result = True
    for vr in vr_list:
        if vr.id is None:
            result = False
            break

    return result


def check_has_duplicate_pos_alleles(vr_list):
    """
    Args:
        vr_list: A list of pysam.VariantRecord objects
    """

    result = False
    IDset = set()
    for vr in vr_list:
        ID = (vr.contig, vr.pos, vr.alleles)
        if ID in IDset:
            result = True
            break
        else:
            IDset.add(ID)

    return result


def check_has_duplicate_id(vr_list):
    """
    Args:
        vr_list: A list of pysam.VariantRecord objects

    Returns:
        True when two or more vr have no ID.
    """

    has_duplicate_id = False
    IDset = set()
    for vr in vr_list:
        if vr.id in IDset:
            has_duplicate_id = True
            break
        else:
            IDset.add(vr.id)

    return has_duplicate_id


def check_has_mateid(vr_list):
    """
    Args:
        vr_list: A list of pysam.VariantRecord objects
    """

    has_mateid = True
    for vr in vr_list:
        if 'MATEID' not in vr.info.keys():
            has_mateid = False
            break
    
    return has_mateid


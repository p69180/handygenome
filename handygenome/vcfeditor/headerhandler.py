import os
import itertools
import warnings

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


def add_INFOmeta(vcfheader, ID, Type, Number, Description):
    vcfheader.add_meta(key='INFO',
                       items=[('ID', ID), ('Type', Type), ('Number', Number), 
                              ('Description', Description)])


def add_FORMATmeta(vcfheader, ID, Type, Number, Description):
    vcfheader.add_meta(key='FORMAT',
                       items=[('ID', ID), ('Type', Type), ('Number', Number), 
                              ('Description', Description)])


def add_generic(vcfheader, key, value):
    vcfheader.add_meta(key=key, value=value)


def add_structured(vcfheader, key, dic):
    vcfheader.add_meta(key=key, items=dic.items())


##################################


def addmeta_MATEID(vcfheader):
    if 'MATEID' not in vcfheader.info:
        vcfheader.add_meta(
            key='INFO',
            items=[('ID', 'MATEID'), ('Type', 'String'), ('Number', 1),
                   ('Description', 'ID of mate breakends')])


def addmeta_END(vcfheader):
    if 'END' not in vcfheader.info:
        vcfheader.add_meta(
            key='INFO',
            items=[('ID', 'END'), ('Type', 'Integer'), ('Number', 1),
                   ('Description', 'End position on CHROM')])


##################################


def merge_vcfheaders(vcfheader_list):
    """Any conflicting INFO or FORMAT keys (with regard to any of Type, 
    Number, or Description) are discarded.
    """
    
    # make a unified VariantHeader object which contains all INFO/FORMAT keys and samples
    vcfheader_list = list(vcfheader_list)
    merged_header = vcfheader_list[0].copy()
    for vcfheader in vcfheader_list[1:]:
        merged_header.merge(vcfheader)
        for sampleid in vcfheader.samples:
            if sampleid not in merged_header.samples:
                merged_header.add_sample(sampleid)

    # find conflicting keys
    info_metadatas = dict()
    format_metadatas = dict()
    for hdr in vcfheader_list:
        for key, val in hdr.info.items():
            info_metadatas.setdefault(key, set())
            info_metadatas[key].add((val.type, val.number, val.description))
        for key, val in hdr.formats.items():
            format_metadatas.setdefault(key, set())
            format_metadatas[key].add((val.type, val.number, val.description))

    conflicting_keys = {'info': set(), 'format': set()}
    for key, val in info_metadatas.items():
        if len(val) > 1:
            conflicting_keys['info'].add(key)
    for key, val in format_metadatas.items():
        if len(val) > 1:
            conflicting_keys['format'].add(key)

    # remove them
    for key in conflicting_keys['info']:
        merged_header.info.remove_header(key)
    for key in conflicting_keys['format']:
        merged_header.formats.remove_header(key)

    return merged_header, conflicting_keys


def merge_vcfheaders_heavy(vcfheader_list):
    """Any conflicting INFO or FORMAT keys 
    (with regard to any of Type, Number, or Description) are discarded."""
    
    vcfheader_list = list(vcfheader_list)
    result = vcfheader_list[0].copy()
    for vcfheader in vcfheader_list[1:]:
        result.merge(vcfheader)
        for sampleid in vcfheader.samples:
            if sampleid not in result.samples:
                result.add_sample(sampleid)

    # get conflicting keys
    conflicting_keys_sum = {'info': set(), 'format': set()}
    for vcfheader1, vcfheader2 in itertools.combinations(vcfheader_list, 2):
        compatible, conflicting_keys = check_header_compatibility(
            vcfheader1, vcfheader2)
        conflicting_keys_sum['info'].update(conflicting_keys['info'])
        conflicting_keys_sum['format'].update(conflicting_keys['format'])

    # remove them
    for key in conflicting_keys_sum['info']:
        result.info.remove_header(key)
    for key in conflicting_keys_sum['format']:
        result.formats.remove_header(key)

    return result


def check_header_compatibility(vcfheader1, vcfheader2):
    conflicting_keys = {'info': list(), 'format': list()}

    info_keys_1 = set(vcfheader1.info.keys())
    info_keys_2 = set(vcfheader2.info.keys())
    info_keys_isec = info_keys_1.intersection(info_keys_2)
    for key in info_keys_isec:
        meta1 = vcfheader1.info[key]
        meta2 = vcfheader2.info[key]
        for meta_key in ('name', 'type', 'number', 'description'):
            if getattr(meta1, meta_key) != getattr(meta2, meta_key):
                conflicting_keys['info'].append(key)
                break

    format_keys_1 = set(vcfheader1.formats.keys())
    format_keys_2 = set(vcfheader2.formats.keys())
    format_keys_isec = format_keys_1.intersection(format_keys_2)
    for key in format_keys_isec:
        meta1 = vcfheader1.formats[key]
        meta2 = vcfheader2.formats[key]
        for meta_key in ('name', 'type', 'number', 'description'):
            if getattr(meta1, meta_key) != getattr(meta2, meta_key):
                conflicting_keys['format'].append(key)
                break

    compatible = (len(conflicting_keys['info']) == 0 and 
                  len(conflicting_keys['format']) == 0)

    return compatible, conflicting_keys


#def check_header_compatibility_list(pysamhdr_list):
#    compatible_list = list()
#    for pysamhdr1, pysamhdr2 in itertools.combinations(pysamhdr_list, 2):
#        compatible, differing_keys = check_header_compatibility_two(pysamhdr1, pysamhdr2)
#        compatible_list.append(compatible)
#    return all(compatible_list)
#
#
#def check_header_compatibility_two(pysamhdr1, pysamhdr2):
#    differing_keys = dict()
#    differing_keys['info_type'] = list()
#    differing_keys['info_number'] = list()
#    differing_keys['format_type'] = list()
#    differing_keys['format_number'] = list()
#
#    info_keys_1 = set(pysamhdr1.info.keys())
#    info_keys_2 = set(pysamhdr2.info.keys())
#    info_keys_isec = info_keys_1.intersection(info_keys_2)
#    for key in info_keys_isec:
#        meta1 = pysamhdr1.info[key]
#        meta2 = pysamhdr2.info[key]
#        if meta1.type != meta2.type:
#            differing_keys['info_type'].append(key)
#        if meta1.number != meta2.number:
#            differing_keys['info_number'].append(key)
#
#    format_keys_1 = set(pysamhdr1.formats.keys())
#    format_keys_2 = set(pysamhdr2.formats.keys())
#    format_keys_isec = format_keys_1.intersection(format_keys_2)
#    for key in format_keys_isec:
#        meta1 = pysamhdr1.formats[key]
#        meta2 = pysamhdr2.formats[key]
#        if meta1.type != meta2.type:
#            differing_keys['format_type'].append(key)
#        if meta1.number != meta2.number:
#            differing_keys['format_number'].append(key)
#
#    compatible = all(len(x) == 0 for x in differing_keys.values())
#
#    return compatible, differing_keys


##################################


def write_mode_arghandler(mode_bcftools, mode_pysam):
    """
    Args:
        bcftools_mode: bcftools style mode string
        mode_pysam: pysam style mode string

    If 'mode_pysam' is set as a non-None value,
    it overrides 'mode_bcftools'.
    """

    if mode_pysam is None:
        return common.PYSAM_MODE_DICT[mode_bcftools]
    else:
        return mode_pysam


##################################


def close_vcfplus_list(vcfplus_list):
    for vcfp in vcfplus_list:
        vcfp.vcf.close()


import textwrap
import math
import warnings

import pysam


NA_VALUES = ('', '.', None)
NA_VALUES_SET = set(NA_VALUES)


def show_info(vr):
    for key in vr.info:
        print(key, 
              get_value_info(vr, key), 
              sep='\t')


def show_format(vr):
    for sampleid in vr.header.samples:
        print('\n' + sampleid)
        for key in vr.samples[sampleid]:
            print(key, 
                  get_value_format(vr, sampleid, key), 
                  sep='\t')


def get_info(vr, key, collapse_tuple=True):
    """
    Returns:
        NA-equivalent values ("", ".", None) are converted into None
        An empty tuple is converted into None
        If the original value is a length-1 tuple, returns its element, 
            with NA conversion.
        If the original value is a tuple longer than 1, return type is tuple.
        If the original value is a tuple longer than 1 composed of only 
            NA-equivalents, returns a tuple with the same length composed 
            of None.
    """

    if key not in vr.header.info.keys():
        raise KeyError(f'Input key "{key}" is not included in the input '
                       f'variant record INFO metadata.')

    if key in vr.info:
        metadata_number = vr.header.info[key].number
        translated_number = get_translated_number_info(vr, key)
        keytype = vr.header.info[key].type
        return modify_InfoFormatValue(
            key,
            vr.info[key], 
            metadata_number,
            translated_number,
            keytype,
            for_write=False,
            typeconv=False,
            collapse_tuple=collapse_tuple,
            )
    else:
        return None

get_value_info = get_info


def get_format(vr, sampleid, key, collapse_tuple=True):
    """
    Returns:
        NA-equivalent values ("", ".", None) are converted into None
        An empty tuple is converted into None
        If the original value is a length-1 tuple, returns its element, with NA conversion.
        If the original value is a tuple longer than 1, return type is tuple.
        If the original value is a tuple longer than 1 composed of only NA-equivalents, returns a tuple with the same length composed of None.
    """

    # sanity checks
    if key == 'GT':
        raise Exception(f'"GT" is not treated with this function.')

    if key not in vr.header.formats.keys():
        raise KeyError(f'Input key "{key}" is not included in the input variant record FORMAT metadata.')

    if sampleid not in vr.header.samples:
        raise KeyError(f'Input sample ID "{sampleid}" is not included in the input variant record sample name list.')

    if key in vr.samples[sampleid]:
        metadata_number = vr.header.formats[key].number
        translated_number = get_translated_number_format(vr, key)
        keytype = vr.header.formats[key].type
        return modify_InfoFormatValue( 
            key,
            vr.samples[sampleid][key], 
            metadata_number,
            translated_number,
            keytype,
            for_write = False,
            typeconv = False,
            collapse_tuple = collapse_tuple,
        )
    else:
        return None

get_value_format = get_format


def get_genotype(vr, sampleid):
    return vr.samples[sampleid].allele_indices, vr.samples[sampleid].phased


def check_NA_info(vr, key):
    if key in vr.info:
        return check_InfoFormatValue_isNA(vr.info[key])
    else:
        return True


def check_NA_format(vr, sampleid, key):
    if sampleid not in vr.samples.keys():
        return True
    else:
        if key in vr.samples[sampleid]:
            return check_InfoFormatValue_isNA(vr.samples[sampleid][key])
        else:
            return True


def check_InfoFormatValue_isNA(val):
    """
    Args:
        val: A raw value obtained from vr.info[key] 
            or vr.samples[sampleid][key]
    """

    if isinstance(val, (list, tuple)):
        if len(val) == 0:
            """
            Value of FOO1 can be an empty tuple when:
                - Number > 1
                AND
                - INFO column looks like: FOO1;FOO2=1;FOO3=10
            """
            return True
        else:
            return set(val).issubset(NA_VALUES_SET)
    else:
        return val in NA_VALUES_SET


#############################################


def set_info(vr, key, val, typeconv=True):
    # exceptions
    #if vr.header.info[key].number == 'G':
        #warnings.warn(f'Metadata number for input key {key} is "G"; nothing is done')
    #    return

    if key not in vr.header.info:
        raise Exception(f'Key {key} is absent from INFO metadata header.')
    else:
        metadata_number = vr.header.info[key].number
        translated_number = get_translated_number_info(vr, key)
        keytype = vr.header.info[key].type
        if (metadata_number != 1) and (not isinstance(val, (list, tuple))):
            val = (val,)

        modified_val = modify_InfoFormatValue(
            key, 
            val, 
            metadata_number,
            translated_number, 
            keytype, 
            for_write=True,
            typeconv=typeconv,
            collapse_tuple=False,
        )
        vr.info[key] = modified_val

set_value_info = set_info


def set_NA_info(vr, key):
    set_value_info(vr, key, None)


def set_format(vr, sampleid, key, val, typeconv=True):
    # exceptions
    #if vr.header.formats[key].number == 'G':
        #warnings.warn(f'Metadata number for input key {key} is "G"; nothing is done')
    #    return
    if key == 'GT':
        raise Exception(f'"GT" is not treated with this function.')

    # main
    if key not in vr.header.formats:
        raise Exception(f'Key {key} is absent from FORMAT metadata header.')
    else:
        metadata_number = vr.header.formats[key].number
        translated_number = get_translated_number_format(vr, key)
        keytype = vr.header.formats[key].type
        if (metadata_number != 1) and (not isinstance(val, (list, tuple))):
            val = (val,)

        modified_val = modify_InfoFormatValue(
            key,
            val, 
            metadata_number,
            translated_number, 
            keytype, 
            for_write=True,
            typeconv=typeconv,
            collapse_tuple=False,
            )
        vr.samples[sampleid][key] = modified_val


set_value_format = set_format  # alias


def set_genotype(vr, sampleid, gt, phased):
    """
    Args:
        gt: A tuple composed of integers
        phased: True or False
    """

    vr.samples[sampleid].allele_indices = gt
    vr.samples[sampleid].phased = phased


def set_NA_format(vr, sampleid, key):
    if key == 'GT':
        set_GT_NA(vr, sampleid)
    else:
        set_value_format(vr, sampleid, key, None)


def set_GT_NA(vr, sampleid):
    vr.samples[sampleid].allele_indices = None
    vr.samples[sampleid].phased = None


def set_GT(vr, sampleid, allele_indices, phased):
    vr.samples[sampleid].allele_indices = allele_indices
    vr.samples[sampleid].phased = phased


#############################################

"""
def make_vr_for_write(vr):
    new_vr = vr.copy()
    # INFO
    for key, old_val in new_vr.info.items():
        set_value_info(new_vr, key, old_val)
    # FORMAT
    for sampleid in new_vr.header.samples:
        for key, old_val in new_vr.samples[sampleid].items():
            set_value_format(new_vr, sampleid, key, old_val)

    return new_vr
"""


def refine_vr_InfoFormatValue(vr):
    """
    Modifies each INFO or FORMAT value into
    vcf writing-compatible form. (e.g. For Type=String 
    (None,None) into ('.','.'))
    """

    # INFO
    for key, old_val in vr.info.items():
        set_value_info(vr, key, old_val)
    # FORMAT
    for sampleid in vr.header.samples:
        for key, old_val in vr.samples[sampleid].items():
            set_value_format(vr, sampleid, key, old_val)


#############################################

# not used
def get_genotype_count(vr):
    if 'GT' in vr.format.keys():
        GT_count_set = set(len(x['GT']) for x in vr.samples.values())
        if len(GT_count_set) != 1:
            raise Exception(f'Genotype counts are different between samples in this variant record:\n{vr}')
        genotype_count = GT_count_set.pop()
    else:
        genotype_count = None

    return genotype_count


def translate_metadata_number(vr, metadata_number):
    if metadata_number == 'A':
        translated_number = len(vr.alts)
    elif metadata_number == 'R':
        translated_number = len(vr.alts) + 1
    elif metadata_number == 'G':
        # combinations with repetition H(len(vr.alleles), 2)
        translated_number = math.comb(len(vr.alleles) + 2 - 1, 2)
    else:
        translated_number = metadata_number

    return translated_number


def prepare_translated_numbers(vr):
    return {
        'A': len(vr.alts),
        'R': len(vr.alleles),
        'G': math.comb(len(vr.alleles) + 1, 2),
            # combinations with repetition H(len(vr.alleles), 2)
        #'other': metadata_number,
    }


def get_translated_number_info(vr, key):
    if key in vr.header.info.keys():
        return translate_metadata_number(vr, vr.header.info[key].number)
    else:
        return None


def get_translated_number_format(vr, key):
    if key in vr.header.formats.keys():
        return translate_metadata_number(vr, vr.header.formats[key].number)
    else:
        return None


#############################################


def modify_InfoFormatValue(
    key, 
    val, 
    metadata_number,
    translated_number, 
    keytype, 
    for_write=False, 
    typeconv=False, 
    collapse_tuple=True,
):
    """
    Assumption:
        "val" is a tuple or a list if and only if metadata "Number" is 1.

    Returns:
        Type of the return value (when collapse_tuple == False):
            tuple when metadata "Number" is not 1.
            atomic when metadata "Number" is 1.

            * When collapse_tuple == True, if the return value is a 1-length
                tuple, its element is returned instead.

        Treatment of NA values:
            '', '.', and None are all treated as NA, 
                and are converted into None in the return value.

        Exceptional cases:
            If "val" is an empty tuple, it is treated as NA.
            Boolean input is returned unchanged.
            Special cases - value is unchanged for these ones:
                INFO/AS_SB_TABLE (Mutect2)
                FORMAT/GT
                Type=Flag

    Args:
        val: A raw value obtained from vr.info[key] or 
            vr.samples[sampleid][key]
        for_write: If True, modification goes for writing into a vcf file.
        keytype: VCF metadata 'Type' value. Only relevant when 
            for_write == True.
        collapse_tuple: If True, length-1 tuple is converted to its element
    """
    # main
    assert translated_number == '.' or isinstance(translated_number, int)

    stop, return_val = modify_InfoFormatValue_handle_special_cases(key, val, keytype)
    if stop:
        return return_val

    len_val, isNA, unifiedNA = modify_InfoFormatValue_get_params(val, for_write, keytype)
    modify_InfoFormatValue_sanity_check(key, val, metadata_number, translated_number, keytype, len_val, isNA)
    modified_val = modify_InfoFormatValue_get_modified_val(
        val, metadata_number, translated_number, isNA, len_val, unifiedNA,
    )
    if typeconv:
        modified_val = modify_InfoFormatValue_typeconv(modified_val, keytype)
    if collapse_tuple:
        modified_val = modify_InfoFormatValue_collapse_tuple(modified_val)

    return modified_val


def modify_InfoFormatValue_handle_special_cases(key, val, keytype):
    if key == 'AS_SB_TABLE':
        # Exception for Mutect2 INFO annotation AS_SB_TABLE
        # Number == 1 but its value looks like '20,14|4,4'
        stop = True
        return_val = ','.join(val)
    elif keytype == 'Flag' or key == 'GT':
        # Flag or GT should not be modified
        stop = True
        return_val = val
    else:
        stop = False
        return_val = None

    return stop, return_val


def modify_InfoFormatValue_get_params(val, for_write, keytype):
    if isinstance(val, (tuple, list)):
        len_val = 1 if len(val) == 0 else len(val)
    else:
        len_val = 1

    isNA = check_InfoFormatValue_isNA(val)
        # True with zero-length tuple
        # True with a non-length-zero tuple composed of NA equivalents

    if for_write:
        unifiedNA = '.' if keytype == 'String' else None
    else:
        unifiedNA = None

    return len_val, isNA, unifiedNA


def modify_InfoFormatValue_sanity_check(
    key, val, metadata_number, translated_number, keytype, 
    len_val, isNA,
):
    # Check if metadata Number and actual length of the value are the same
    if (
        isinstance(translated_number, int)
        and (len_val != translated_number)
        and (not isNA)
    ):
        raise Exception(
            f'"Number" metadata and actual number of the value '
            f'does not match.\n'
            f'key = {key}, val = {val}, '
            f'translated_number = {translated_number}, '
            f'keytype = {keytype}'
        )
    # Assumed that raw value is not a sequence type if and only if 
    # "Number" is 1.
    if (
        (metadata_number == 1) 
        and isinstance(val, (tuple, list))
    ):
        raise Exception(
            f'Unexpected case: raw value is a tuple or a list'
            f'when metadata Number is 1.'
        )
    if (
        (metadata_number != 1) 
        and (not isinstance(val, (tuple, list)))
    ):
        raise Exception(
            f'Unexpected case: raw value is not a tuple or a list'
            f'when metadata Number is not 1.'
        )


def modify_InfoFormatValue_get_modified_val(
    val,
    metadata_number,
    translated_number,
    isNA,
    len_val,
    unifiedNA,
):
    if isNA:
        if translated_number == '.':
            modified_val = tuple([unifiedNA] * len_val)
                # len_val is 1 when val is an empty tuple
        else:
            if metadata_number == 1:
                modified_val = unifiedNA
            else:
                modified_val = tuple([unifiedNA] * translated_number)
    else:
        if metadata_number == 1:
            # turn NAvalues into unifiedNA
            if val in NA_VALUES:
                modified_val = unifiedNA
            else:
                modified_val = val
        else:
            # turn NAvalues into unifiedNA
            modified_val = tuple(unifiedNA if x in NA_VALUES else x
                                 for x in val)

    return modified_val


def modify_InfoFormatValue_typeconv(modified_val, keytype):
    if keytype == 'Flag':
        return bool(modified_val)
    else:
        if keytype == 'String' or keytype == 'Character':
            converter_func = str
        elif keytype == 'Integer':
            converter_func = int
        elif keytype == 'Float':
            converter_func = float

        def converter(val):
            return val if val in NA_VALUES else converter_func(val)

        if isinstance(modified_val, (tuple, list)):
            new_val = tuple(converter(x) for x in modified_val)
        else:
            new_val = converter(modified_val)

        return new_val


def modify_InfoFormatValue_collapse_tuple(modified_val):
    if isinstance(modified_val, (tuple, list)):
        if len(modified_val) == 1:
            return modified_val[0]
        else:
            return modified_val
    else:
        return modified_val




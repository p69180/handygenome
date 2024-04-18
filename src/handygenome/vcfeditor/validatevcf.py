import os
import re
import multiprocessing
import functools

import pysam

import handygenome.tools as tools
import handygenome.logutils as logutils
import handygenome.vcfeditor.headerhandler as headerhandler
import handygenome.variant.infoformat as infoformat


#def parse_header(vcfheader):
#    header_data = {'info': dict(), 'format': dict()}
#    for key, val in vcfheader.info.items():
#        metadata = {
#            'name': val.name,
#            'number': val.number,
#        }
#        header_data['info'][metadata['name']] = metadata
#    for key, val in vcfheader.formats.items():
#        metadata = {
#            'name': val.name,
#            'number': val.number,
#        }
#        header_data['format'][metadata['name']] = metadata
#    return header_data


def get_annotval_length(val):
    if isinstance(val, tuple):
        return len(val)
    else:
        return 1


def collect_annot_status_helper(key, val, data, translnum_dict, metadata_numdict):
    data.setdefault(key, {'length': set(), 'conflict': False, 'missing': True})
    val_len = get_annotval_length(val)
    data[key]['length'].add(val_len)
    if key in metadata_numdict:
        data[key]['missing'] = False
        try:
            transl_num = translnum_dict[metadata_numdict[key]]
        except KeyError:
            transl_num = metadata_numdict[key]
        if transl_num != val_len:
            data[key]['conflict'] = True


def collect_annot_status_itervcf(vcf_path, linenum=None):
    with pysam.VariantFile(vcf_path) as vcf:
        for idx, vr in enumerate(vcf):
            if linenum is not None:
                if idx % linenum == 0:
                    print(idx, flush=True)
            translnum_dict = infoformat.prepare_translated_numbers(vr)
            info_annots = dict(vr.info.items())
            format_annots = dict()
            for sid, subdata in vr.samples.items():
                format_annots[sid] = dict(subdata.items())
            yield translnum_dict, info_annots, format_annots


def collect_annot_status_unitjob(
    info_data, 
    format_data, 
    metadata_infonums, 
    metadata_formatnums,
    translnum_dict, 
    info_annots, 
    format_annots,
    untouched_keys,
):
    for key, val in info_annots.items():
        collect_annot_status_helper(
            key=key, 
            val=val, 
            data=info_data, 
            translnum_dict=translnum_dict,
            metadata_numdict=metadata_infonums,
        )
    # format
    for sid, subdata in format_annots.items():
        for key, val in subdata.items():
            if key in untouched_keys:
                continue
            collect_annot_status_helper(
                key=key, 
                val=val, 
                data=format_data, 
                translnum_dict=translnum_dict,
                metadata_numdict=metadata_formatnums,
            )


def collect_annot_status(vcf_path, linenum=None):
    info_data = dict()
    format_data = dict()

    with pysam.VariantFile(vcf_path) as vcf:
        metadata_infonums = {
            x.name: x.number 
            for x in vcf.header.info.values()
        }
        metadata_formatnums = {
            x.name: x.number 
            for x in vcf.header.formats.values()
        }


    untouched_keys = ['GT', 'PL']
    for (
        translnum_dict, 
        info_annots, 
        format_annots,
    ) in collect_annot_status_itervcf(vcf_path, linenum=linenum):
        collect_annot_status_unitjob(
            info_data, 
            format_data, 
            metadata_infonums, 
            metadata_formatnums,
            translnum_dict, 
            info_annots, 
            format_annots,
            untouched_keys,
        )

    return info_data, format_data


# not used; no faster than version without mp
def collect_annot_status_parallel(vcf_path, nproc=1):
    manager = multiprocessing.Manager()
    info_data_mngr = manager.dict()
    format_data_mngr = manager.dict()

    with pysam.VariantFile(vcf_path) as vcf:
        metadata_infonums = {
            x.name: x.number 
            for x in vcf.header.info.values()
        }
        metadata_formatnums = {
            x.name: x.number 
            for x in vcf.header.formats.values()
        }

    with multiprocessing.Pool(nproc) as pool:
        args_iter = (
            (
                info_data_mngr, 
                format_data_mngr, 
                metadata_infonums, 
                metadata_formatnums,
                translnum_dict, 
                info_annots, 
                format_annots,
            )
            for (
                translnum_dict, 
                info_annots, 
                format_annots,
            ) in collect_annot_status_itervcf(vcf_path)
        )
        pool.starmap(collect_annot_status_unitjob, args_iter)

    info_data = dict(info_data_mngr)
    format_data = dict(format_data_mngr)

    return info_data, format_data
            

def modify_header(vcfheader, info_data, format_data):
    excl_info = functools.reduce(
        set.union,
        [
            set(vcfheader.info.keys()).difference(info_data.keys()),  
                # keys present in old header but absent from VariantRecord annotations
            set(key for (key, val) in info_data.items() if not val['conflict']),
                # keys with field length conflict
        ],
    )
    excl_format = functools.reduce(
        set.union,
        [
            set(vcfheader.formats.keys()).difference(format_data.keys()),
            set(key for (key, val) in format_data.items() if not val['conflict']),
        ],
    )
    newhdr = headerhandler.copyheader(
        vcfheader, 
        excl_info=excl_info,
        excl_format=excl_format,
    )
    for metakey, status_data in zip(['INFO', 'FORMAT'], [info_data, format_data]):
        for key, subdic in status_data.items():
            if (subdic['missing'] or subdic['conflict']):
                if subdic['missing']:
                    if len(subdic['length']) == 1:
                        number = subdic['length'].pop()
                    else:
                        number = '.'
                    description = '.'
                elif subdic['conflict']:
                    number = '.'
                    description = vcfheader.info[key].description

                newhdr.add_meta(
                    key=metakey,
                    items=[
                        ('ID', key), 
                        ('Type', 'String'), 
                        ('Number', number), 
                        ('Description', description),
                    ],
                )

    return newhdr


def make_metavalid_vcf(vcf_path, out_vcf_path):
    assert out_vcf_path.endswith('.vcf.gz')
    info_data, format_data = collect_annot_status(vcf_path)
    with pysam.VariantFile(vcf_path) as vcf:
        newhdr = modify_header(vcf.header, info_data, format_data)
        with pysam.VariantFile(out_vcf_path, mode='wz', header=newhdr) as out_vcf:
            for vr in vcf:
                out_vcf.write(vr)




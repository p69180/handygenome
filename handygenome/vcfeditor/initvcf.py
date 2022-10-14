import os
import warnings

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
varianthandler = importlib.import_module('.'.join([top_package_name, 'variant', 'varianthandler']))


def write_empty_vcf(outfile_path, chromdict=None, samples=None, pysamhdr=None, 
                    vcfver=common.DEFAULT_VCFVER):
    """
    Args:
        outfile_path: Output vcf file path.
        chromdict: Designed as a mandatory one since vcf header without 
            contig metadata may result in segmentation fault with writing 
            to a file.
        samples: Iterable containing sample names.
        pysamhdr: pysam.VariantHeader object. All metadata from this, 
            except fileformat and contig, is used.
        vcfver: Default 4.3
    """

    with open(outfile_path, 'w') as outfile:
        outfile.write(f'##fileformat=VCFv{vcfver}\n') # fileformat

        if chromdict is not None:
            for chrom, length in chromdict.items(): # contig
                outfile.write(f'##contig=<ID={chrom},length={length}>\n')

        if pysamhdr is not None: # other metadata
            for hdr_rec in pysamhdr.records:
                if chromdict is not None:
                    if (
                            (hdr_rec.type != 'CONTIG') and 
                            (hdr_rec.key != 'fileformat')):
                        outfile.write(str(hdr_rec))
                else:
                    if hdr_rec.key != 'fileformat':
                        outfile.write(str(hdr_rec))

        headerline = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 
                      'INFO']
        if samples is not None:
            headerline.append('FORMAT')
            for sample in samples:
                headerline.append(sample)
        outfile.write('\t'.join(headerline) + '\n')


def create_empty_pysamvcf(chromdict=None, samples=None, pysamhdr=None, 
                          vcfver=common.DEFAULT_VCFVER):
    """
    Args:
        outfile_path: Output vcf file path.
        chromdict: Designed as a mandatory one since vcf header without 
            contig metadata may result in segmentation fault with writing to 
            a file.
        samples: Iterable containing sample names.
        pysamhdr: pysam.VariantHeader object. All metadata from this, except 
            fileformat and contig, is used.
        vcfver: Default 4.3
    """

    tf_path = workflow.get_tmpfile_path(delete=False)
    write_empty_vcf(tf_path, chromdict, samples, pysamhdr, vcfver)
    empty_vcf = pysam.VariantFile(tf_path, mode = 'r')
    os.remove(tf_path)

    return empty_vcf


######################################

def create_header(chromdict, samples=None, vcfheader=None):
    """
    Args:
        chromdict: Designed as a mandatory one since vcf header without 
            contig metadata may result in segmentation fault with writing 
            to a file.
        samples: An iterable containing sample names.
        vcfheader: Another pysam.VariantHeader object. All metadata from 
            this, except fileformat and contig, is used.

    Returns:
        A pysam.VariantHeader object
    """

    result = pysam.VariantHeader()

    for chrom, length in chromdict.items():
        result.contigs.add(chrom, length)

    if samples is not None:
        for sampleid in samples:
            result.samples.add(sampleid)

    if vcfheader is not None:
        for rec in vcfheader.records:
            if (rec.type != 'CONTIG') and (rec.key != 'fileformat'):
                result.add_record(rec)

    return result


#def get_updated_header(vcfhdr, chromdict = None, samples = None):
#    new_header = create_header(chromdict, samples, vcfhdr)
#    for header_rec in vcfhdr.records:
#        if ( header_rec.type != 'CONTIG' ) and ( header_rec.key != 'fileformat' ):
#            new_header.add_record(header_rec)
#
#    return new_header


def create_vr(chromdict, vcfspec=None, end=None, samples=None, 
              vcfheader=None, use_header=False):
    if use_header:
        vr = vcfheader.new_record()
    else:
        hdr = create_header(chromdict, samples, vcfheader)
        vr = hdr.new_record()

    if vcfspec is not None:
        varianthandler.apply_vcfspec(vr, vcfspec)

    if end is not None:
        vr.stop = end

    return vr

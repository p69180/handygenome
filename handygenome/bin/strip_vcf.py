import textwrap

import pysam

import handygenome.workflow as workflow
import handygenome.vcfeditor.indexing as indexing


def argument_parser(cmdargs):
    def sanity_check(args):
        pass

    parser_dict = workflow.init_parser(
        description=f'Makes a new vcf file where all INFO and FORMAT annotations are removed. Relevant header metadata lines are also removed.'
    )

    workflow.add_infile_arg(parser_dict['required'], required=True)
    workflow.add_outfile_arg(parser_dict['required'], required=True)
    workflow.add_outfmt_arg(parser_dict['optional'], required=False, default='z')
    workflow.add_index_arg(parser_dict)

    args = parser_dict['main'].parse_args(cmdargs)
    sanity_check(args)

    return args


def make_new_header(in_vcf):
    new_header = pysam.VariantHeader()
    for chrom, contigitem in in_vcf.header.contigs.items():
        new_header.contigs.add(contigitem.name, contigitem.length)

    return new_header


def write_new_vrs(in_vcf, out_vcf):
    for vr in in_vcf.fetch():
        new_vr = out_vcf.new_record()
        for key in ('contig', 'pos', 'id', 'ref', 'alts'):
            setattr(new_vr, key, getattr(vr, key))
        out_vcf.write(new_vr)


def main(cmdargs):
    args = argument_parser(cmdargs)
    in_vcf = pysam.VariantFile(args.infile_path)
    new_header = make_new_header(in_vcf)
    out_vcf = pysam.VariantFile(args.outfile_path, mode=args.mode_pysam, header=new_header)
    write_new_vrs(in_vcf, out_vcf)
    out_vcf.close()
    in_vcf.close()
        
    if not args.donot_index:
        indexing.index_vcf(args.outfile_path)

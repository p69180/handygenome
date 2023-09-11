import textwrap
import gzip

import pysam

import handygenome.tools as tools
import handygenome.refgenome as refgenome
import handygenome.workflow as workflow
import handygenome.workflow.toolsetup as toolsetup
import handygenome.vcfeditor.indexing as indexing
import handygenome.vcfeditor.initvcf as initvcf
import handygenome.vcfeditor.headerhandler as headerhandler
import handygenome.variant.vcfspec as libvcfspec
import handygenome.variant.varianthandler as varianthandler


def argument_parser(cmdargs):
    def sanity_check(args):
        pass

    parser_dict = workflow.init_parser(
        description='Turns a Bismark *.cov.gz output file into a VCF.')

    workflow.add_infile_arg(
        parser_dict['required'], required=True,
        help=f'A Bismark output file which matched the pattern "*.cov.gz"')
    workflow.add_outfile_arg(parser_dict['required'], required=True, 
                             must_not_exist='ask')
    workflow.add_fasta_arg(parser_dict['required'], required=True)
    workflow.add_outfmt_arg(parser_dict['optional'], required=False)

    workflow.add_logging_args(parser_dict)
    workflow.add_index_arg(parser_dict)

    parser_dict['required'].add_argument(
        '--sampleid', required=True,
        help=f'The sample name for the input file.')

    args = parser_dict['main'].parse_args(cmdargs)
    sanity_check(args)

    return args


def add_meta(vcfheader):
    # ALT meta
    vcfheader.add_meta(
        key='ALT', items=[('ID', libvcfspec.CPGMET_ALT), 
                          ('Description', 'CpG site methylation')])

    # INFO meta
    headerhandler.addmeta_END(vcfheader)

    # FORMAT meta
    vcfheader.add_meta(
        key='FORMAT',
        items=[
            ('ID', 'nread_meth_bisulf_top'), ('Type', 'Integer'), 
            ('Number', 1),
            ('Description', ('The number of reads supporting cytosine '
                             'methylation, '
                             'belonging to the top(plus) strand from a '
                             'bisulfite sequencing result.'))])

    vcfheader.add_meta(
        key='FORMAT',
        items=[
            ('ID', 'nread_total_bisulf_top'), ('Type', 'Integer'), 
            ('Number', 1),
            ('Description', ('The number of all reads, '
                             'belonging to the top(plus) strand from a '
                             'bisulfite sequencing result.'))])

    vcfheader.add_meta(
        key='FORMAT',
        items=[
            ('ID', 'meth_fract_bisulf_top'), ('Type', 'Float'), 
            ('Number', 1),
            ('Description', ('The fraction of reads supporting cytosine '
                             'methylation, '
                             'belonging to the top(plus) strand from a '
                             'bisulfite sequencing result.'))])

    vcfheader.add_meta(
        key='FORMAT',
        items=[
            ('ID', 'nread_meth_bisulf_bot'), ('Type', 'Integer'), 
            ('Number', 1),
            ('Description', ('The number of reads supporting cytosine '
                             'methylation, '
                             'belonging to the bottom(minus) strand from a '
                             'bisulfite sequencing result.'))])

    vcfheader.add_meta(
        key='FORMAT',
        items=[
            ('ID', 'nread_total_bisulf_bot'), ('Type', 'Integer'), 
            ('Number', 1),
            ('Description', ('The number of all reads, '
                             'belonging to the bottom(minus) strand from a '
                             'bisulfite sequencing result.'))])

    vcfheader.add_meta(
        key='FORMAT',
        items=[
            ('ID', 'meth_fract_bisulf_bot'), ('Type', 'Float'), 
            ('Number', 1),
            ('Description', ('The fraction of reads supporting cytosine '
                             'methylation, '
                             'belonging to the bottom(minus) strand from a '
                             'bisulfite sequencing result.'))])


def init_vcfheader(chromdict, sampleid):
    vcfheader = initvcf.create_header(chromdict=chromdict, samples=[sampleid])
    add_meta(vcfheader)

    return vcfheader


def parse_linesp(linesp):
    assert linesp[1] == linesp[2], (
        f'Column 2 and Column 3 are different in this line: {linesp}')

    chrom = linesp[0]
    pos = int(linesp[1])
    meth_fract = 0.01 * float(linesp[3])
    nread_meth = int(linesp[4])
    nread_total = int(linesp[4]) + int(linesp[5])

    linesp_parsed = {'chrom': chrom, 'pos': pos, 'meth_fract': meth_fract,
                     'nread_meth': nread_meth, 'nread_total': nread_total}

    return linesp_parsed


def parse_one_line(infile):
    try:
        line = next(infile)
    except StopIteration:
        return None
    else:
        linesp = tools.get_linesp(line)
        linesp_parsed = parse_linesp(linesp)
        return linesp_parsed


def parse_two_lines(infile):
    linesp_parsed_prev = parse_one_line(infile)
    linesp_parsed = parse_one_line(infile)

    return linesp_parsed_prev, linesp_parsed


def handle_singleton(linesp_parsed, fasta, vcfheader, vr_buffer):
    if linesp_parsed['chrom'] in fasta.references:
        ref = fasta.fetch(linesp_parsed['chrom'], linesp_parsed['pos'] - 2, 
                          linesp_parsed['pos'] + 1)
        if ref[1] == 'C':
            if ref[2] == 'G':
                pos = linesp_parsed['pos']
                vr = vcfheader.new_record(start=(pos - 1))
                vr.contig = linesp_parsed['chrom']
                vr.ref = 'CG'
                vr.alts = ('<CPGMET>',)

                vr.samples[0]['nread_meth_bisulf_top'] = (
                    linesp_parsed['nread_meth'])
                vr.samples[0]['nread_total_bisulf_top'] = (
                    linesp_parsed['nread_total'])
                vr.samples[0]['meth_fract_bisulf_top'] = (
                    linesp_parsed['meth_fract'])
            else:
                vr = None
                #raise Exception(f'Unexpected REF sequence on a singleton line '
                                #f'position: "C" "nonG"\nline contents: {linesp}')
        elif ref[1] == 'G':
            if ref[0] == 'C':
                pos = linesp_parsed['pos'] - 1
                vr = vcfheader.new_record(start=(pos - 1))
                vr.contig = linesp_parsed['chrom']
                vr.ref = 'CG'
                vr.alts = ('<CPGMET>',)

                vr.samples[0]['nread_meth_bisulf_bot'] = (
                    linesp_parsed['nread_meth'])
                vr.samples[0]['nread_total_bisulf_bot'] = (
                    linesp_parsed['nread_total'])
                vr.samples[0]['meth_fract_bisulf_bot'] = (
                    linesp_parsed['meth_fract'])
            else:
                vr = None
                #raise Exception(f'Unexpected REF sequence on a singleton line '
                                #f'position: "nonC" "G"\n'
                                #f'line contents: {linesp_parsed}')
        else:
            vr = None
            #raise Exception(f'Unexpected REF sequence on a singleton line '
                            #f'position: neither C nor G'
                            #f'\nline contents: {linesp_parsed}')
    else:
        vr = None

    if vr is not None:
        vr_buffer.append(vr)


def handle_doublet(linesp_parsed_prev, linesp_parsed, vcfheader, fasta,
                   vr_buffer):
    if linesp_parsed['chrom'] in fasta.references:
        ref = fasta.fetch(linesp_parsed['chrom'], 
                          linesp_parsed_prev['pos'] - 1, linesp_parsed['pos'])
        if ref == 'CG':
            pos = linesp_parsed_prev['pos']
            vr = vcfheader.new_record(start=(pos-1))
            vr.contig = linesp_parsed['chrom']
            vr.ref = 'CG'
            vr.alts = ('<CPGMET>',)
            
            vr.samples[0]['nread_meth_bisulf_top'] = (
                linesp_parsed_prev['nread_meth'])
            vr.samples[0]['nread_total_bisulf_top'] = (
                linesp_parsed_prev['nread_total'])
            vr.samples[0]['meth_fract_bisulf_top'] = (
                linesp_parsed_prev['meth_fract'])

            vr.samples[0]['nread_meth_bisulf_bot'] = (
                linesp_parsed['nread_meth'])
            vr.samples[0]['nread_total_bisulf_bot'] = (
                linesp_parsed['nread_total'])
            vr.samples[0]['meth_fract_bisulf_bot'] = (
                linesp_parsed['meth_fract'])
        else:
            vr = None
    else:
        vr = None

    if vr is not None:
        vr_buffer.append(vr)


def check_linesp_adjacency(linesp_parsed_prev, linesp_parsed):
    return (
        (linesp_parsed_prev['pos'] == linesp_parsed['pos'] - 1) and
        (linesp_parsed_prev['chrom'] == linesp_parsed['chrom']))


def get_linesp_buffer(infile):
    linesp_parsed = parse_one_line(infile)
    if linesp_parsed is None:
        return None
    else:
        linesp_buffer = list()
        linesp_buffer.append(linesp_parsed)
        while True:
            linesp_parsed = parse_one_line(infile)
            if linesp_parsed is None:
                break
            else:
                linesp_buffer.append(linesp_parsed)
                if check_linesp_adjacency(linesp_buffer[-2], linesp_buffer[-1]):
                    break
                else:
                    continue

        return linesp_buffer


def get_vr_buffer(linesp_buffer, fasta, vcfheader):
    if linesp_buffer is None:
        return None
    else:
        vr_buffer = list()

        if len(linesp_buffer) == 1:
            handle_singleton(linesp_buffer[0], fasta, vcfheader, vr_buffer)
        elif len(linesp_buffer) >= 2:
            if check_linesp_adjacency(linesp_buffer[-2], linesp_buffer[-1]):
                for linesp_parsed in linesp_buffer[:-2]:
                    handle_singleton(linesp_parsed, fasta, vcfheader,
                                     vr_buffer)
                vr = handle_doublet(linesp_buffer[-2], linesp_buffer[-1],
                                    vcfheader, fasta, vr_buffer)
            else:
                for linesp_parsed in linesp_buffer:
                    handle_singleton(linesp_parsed, fasta, vcfheader,
                                     vr_buffer)

        return vr_buffer


def write_oufile(infile_path, outfile_path, mode_pysam, vcfheader, fasta, 
                 chromdict, logger):
    N_VR = 0
    out_vr_list = list()
    with gzip.open(infile_path, 'rt') as infile:
        while True:
            linesp_buffer = get_linesp_buffer(infile)
            vr_buffer = get_vr_buffer(linesp_buffer, fasta, vcfheader)

            if vr_buffer is None:
                # arrived at the end of the input file
                break
            else:
                N_VR += len(vr_buffer)
                for vr in vr_buffer:
                    out_vr_list.append(vr)

            if N_VR % 1_000_000 == 0:
                logger.info(
                    f'{N_VR:,} variant records processed. Current '
                    f'raw input line: {linesp_buffer[-1]}')

    # sort
    logger.info('Sorting variant records')
    out_vr_list.sort(key=varianthandler.get_vr_sortkey(chromdict))

    # write
    logger.info('Writing variant records')
    with pysam.VariantFile(outfile_path, mode=mode_pysam,
                           header=vcfheader) as out_vcf:
        for vr in out_vr_list:
            out_vcf.write(vr)


def main(cmdargs):
    args = argument_parser(cmdargs)
    logger = toolsetup.setup_logger(args)

    logger.info('BEGINNING')

    fasta = pysam.FastaFile(args.fasta_path)
    chromdict = refgenome.Chromdict.from_fasta(fasta)
    vcfheader = init_vcfheader(chromdict, args.sampleid)

    write_oufile(args.infile_path, args.outfile_path, args.mode_pysam, 
                 vcfheader, fasta, chromdict, logger)
                
    if not args.donot_index:
        indexing.index_vcf(args.outfile_path)

    logger.info('ALL SUCCESSFULLY FINISHED')


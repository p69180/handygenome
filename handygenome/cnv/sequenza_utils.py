import os
import sys
import re
import subprocess
import multiprocessing
import tempfile

import sequenza.commands

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


DEFAULT_GCWIG_PATHS = common.RefverDict(
    {
        'GRCh37': '/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.gc50Base.txt.gz',
        'GRCh38': '/home/users/pjh/References/reference_genome/GRCh38/GCA_for_alignment_pipelines/no_alt_plus_hs38d1/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.gc50Base.txt.gz',
    }
)
DEFAULT_HOM = 0.9
DEFAULT_HET = 0.25


def run_sequenza_utils(args):
    sys.argv = ['sequenza_utils'] + args
    try:
        result = sequenza.commands.main()
    except SystemExit:
        pass
    #sys.exit(result)


def run_sequenza_utils_subprocess(args):
    pass
    #p = subprocess.run


def fasta_path_to_gcwiggle_path(fasta_path):
    """e.g. ref.fasta -> ref.gc50.wig.gz"""
    return re.sub('\.[^.]+$', f'.gc{window_size}.wig.gz', fasta_path)


def seqzutils_gc_wiggle(fasta_path, outfile_path=None, window_size=50):
    if outfile_path is None:
        outfile_path = fasta_path_to_gcwiggle_path(fasta_path)

    run_sequenza_utils(['gc_wiggle', 
                        '-f', fasta_path, 
                        '-o', outfile_path, 
                        '-w', str(window_size)])


def bam2seqz_sanity_check(force_outfile_name, outfile_path, refver, 
                          fasta_path, gcwiggle_path):
    if not force_outfile_name:
        if not outfile_path.endswith('.seqz.gz'):
            raise Exception(f'Output file name must end with ".seqz.gz".')

    if not ((
                (refver is not None) and 
                (fasta_path is None) and 
                (gcwiggle_path is None)) or
            (
                (refver is None) and
                (fasta_path is not None) and
                (gcwiggle_path is not None))):
        raise Exception(
            f'Argument usage must be one of following:\n'
            f'1) "refver" set, "fasta_path" and "gcwiggle_path" not set;\n'
            f'2) "refver" not set, "fasta_path" and "gcwiggle_path" set.'
            )


@common.get_deco_num_set(('refver', 'fasta_path'), 1)
def seqzutils_bam2seqz(tbam_path, nbam_path, outfile_path, 
                       refver=None, fasta_path=None, gcwiggle_path=None, 
                       force_outfile_name=False,
                       chrom=None, start0=None, end0=None,
                       hom=DEFAULT_HOM, het=DEFAULT_HET):
    def seqz_filter(infile_path, outfile_path, chrom, start0, end0):
        with \
                open(infile_path, 'rt') as infile, \
                common.openfile(outfile_path, 'w') as outfile:
            header_line = next(infile)
            outfile.write(header_line)

            for line in infile:
                linesp = line.split('\t')
                chrom_line = linesp[0]
                pos0_line = int(linesp[1]) - 1
                if region_check(chrom, start0, end0, chrom_line, pos0_line):
                    outfile.write(line)

    def region_check(chrom, start0, end0, chrom_line, pos0_line):
        if chrom is not None:
            if start0 is None or end0 is None:
                return (chrom == chrom_line)
            else:
                return (chrom_line == chrom and 
                        pos0_line >= start0 and
                        pos0_line < end0)
        else:
            return True

    def run_bam2seqz(outfile_path, tbam_path, nbam_path, gcwiggle_path, 
                     fasta_path, hom, het, region):
        args = ['bam2seqz',
                '-t', tbam_path,
                '-n', nbam_path,
                '-gc', gcwiggle_path,
                '-F', fasta_path,
                '-o', outfile_path,
                '--hom', str(hom),
                '--het', str(het),
                '--samtools', common.SAMTOOLS,
                '--tabix', common.TABIX]
        if region is not None:
            args.extend(['-C', region])

        run_sequenza_utils(args)

    def get_fasta_gcwiggle_paths(refver, fasta_path, gcwiggle_path):
        if refver is not None:
            fasta_path = common.DEFAULT_FASTA_PATHS[refver]
            gcwiggle_path = DEFAULT_GCWIG_PATHS[refver]
        else:
            # fasta_path is not None in this case
            if gcwiggle_path is None:
                gcwiggle_path = fasta_path_to_gcwiggle_path(fasta_path)

        return fasta_path, gcwiggle_path

    def get_region(chrom, start0, end0):
        if chrom is not None:
            if (start0 is not None) and (end0 is not None):
                region = f'{chrom}:{start0 + 1}-{end0}'
            else:
                region = chrom
        else:
            region = None

        return region

    # sanity check
    bam2seqz_sanity_check(force_outfile_name, outfile_path, refver, 
                          fasta_path, gcwiggle_path)

    # parameter refinements
    fasta_path, gcwiggle_path = get_fasta_gcwiggle_paths(refver, fasta_path, 
                                                         gcwiggle_path)
    region = get_region(chrom, start0, end0)

    # main
    if region is None:
        run_bam2seqz(outfile_path, tbam_path, nbam_path, gcwiggle_path, 
                     fasta_path, hom, het, region)
    else:
        with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            fifo_path = os.path.join(tmpdir, 'fifo')
            os.mkfifo(fifo_path)

            p_fiforeader = multiprocessing.Process(
                target=seqz_filter,
                args=(fifo_path, outfile_path, chrom, start0, end0))
            p_fifowriter = multiprocessing.Process(
                target=run_bam2seqz,
                args=(fifo_path, tbam_path, nbam_path, gcwiggle_path, 
                      fasta_path, hom, het, region))

            p_fiforeader.start()
            p_fifowriter.start()

            p_fifowriter.join()
            p_fiforeader.join()


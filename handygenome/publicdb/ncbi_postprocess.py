import os

import pysam

import handygenome.tools as tools
import handygenome.fastahandler as fastahandler
import handygenome.logutils as logutils
#import handygenome.refgenome as refgenome


def postprocess_fasta(infile_path, outfile_path):
    logutils.print_timestamp(f'Postprocessing downloaded RefSeq fasta file')
    if tools.check_file_is_plain(infile_path):
        pysam.tabix_compress(infile_path, outfile_path)
        fastahandler.make_index(outfile_path)
    elif tools.check_file_is_gzipped(infile_path):
        tmppath = tools.get_tmpfile_path(infile_path)
        tools.unzip(infile_path, tmppath, rm_src=False)
        pysam.tabix_compress(tmppath, outfile_path)
        fastahandler.make_index(outfile_path)
        os.remove(tmppath)
    else:
        raise Exception(f'Input file is neither plain text nor gzipped.')
        

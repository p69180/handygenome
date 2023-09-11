import os
import json
import hashlib

import handygenome.logutils as logutils
import handygenome.tools as tools
import handygenome.publicdb.ncbi_cache as libncbicache
import handygenome.fastahandler as fastahandler


class PostprocessorBase:
    def __init__(
        self,
        standard_refver,
        edited_path,
        custom_path,
        RefSeq_refver,
        species,
        chrom_converter,
        unedited_file_checksum_path,
        unedited_file_wascustom_path,
        calc_checksum,
    ):
        self.standard_refver = standard_refver
        self.edited_path = edited_path
        self.custom_path = custom_path
        self.RefSeq_refver = RefSeq_refver
        self.species = species
        self.chrom_converter = chrom_converter
        self.unedited_file_checksum_path = unedited_file_checksum_path
        self.unedited_file_wascustom_path = unedited_file_wascustom_path
        self.calc_checksum = calc_checksum
        self.iscustom = (self.custom_path is not None)

    def get_helpmsg(self):
        return f'(refver={self.standard_refver}, path={repr(self.get_unedited_path())})'

    def make_processed_file(self):
        logutils.log(
            f'Postprocessing raw data file {self.get_helpmsg()}', 
            level='info',
        )
        self.postprocess(
            infile_path=self.get_unedited_path(), 
            outfile_path=self.edited_path, 
            chrom_converter=(
                self.chrom_converter
                if self.custom_path is None else
                None
            ),
        )
        logutils.log(f'Finished postprocessing', level='info')

    def prepare_processed_file(self):
        if self.calc_checksum:
            unedited_file_changed = self.compare_unedited_file_checksum()
            if unedited_file_changed:
                self.make_processed_file()
            else:
                assert os.path.exists(self.edited_path), (
                    f'Edited fasta file does not exist: {self.edited_path}'
                )
        else:
            if not os.path.exists(self.edited_path):
                self.make_processed_file()

    def compare_unedited_file_checksum(self, digest='sha256'):
        # load stored checksum which was calculated from previous opening
        if os.path.exists(self.unedited_file_checksum_path):
            with open(self.unedited_file_checksum_path, 'rt') as f:
                stored_checksum = f.read()
        else:
            stored_checksum = None

        if stored_checksum is None:
            wascustom = None
        else:
            if os.path.exists(self.unedited_file_wascustom_path):
                with open(self.unedited_file_wascustom_path, 'rt') as f:
                    wascustom = json.load(f)
            else:
                wascustom = None

        # calculate newly opened fasta checksum ; this takes a few seconds
        unedited_path = self.get_unedited_path()
        logutils.log(
            f'Calculating checksum of the unedited file. It takes a few seconds. {self.get_helpmsg()}', 
            level='info',
            add_locstring=False,
        )
        with open(unedited_path, 'rb') as f:
            new_checksum = hashlib.file_digest(f, digest).hexdigest()  
                # digest() returns bytes object, while hexdigest() returns str object

        # result
        if stored_checksum is None:
            unedited_file_changed = True
        else:
            unedited_file_changed = (stored_checksum != new_checksum) or (wascustom != self.iscustom)

        with open(self.unedited_file_checksum_path, 'wt') as f:
            f.write(new_checksum)
        with open(self.unedited_file_wascustom_path, 'wt') as f:
            json.dump(self.iscustom, f)

        return unedited_file_changed


class FastaProcessor(PostprocessorBase):
    def get_unedited_path(self):
        if self.custom_path is None:
            return libncbicache.get_fasta_path(self.RefSeq_refver, self.species)
        else:
            return self.custom_path

    @staticmethod
    def postprocess(infile_path, outfile_path, chrom_converter=None):
        if os.path.exists(outfile_path):
            os.remove(outfile_path)

        if chrom_converter is None:
            if tools.check_textfile(infile_path):
                os.symlink(infile_path, outfile_path)
            elif tools.check_gzipped(infile_path):
                tools.unzip(infile_path, outfile_path, rm_src=False)
            else:
                raise Exception(f'Input file is neither plain text nor gzipped.')
        else:
            fastahandler.rename_fasta(infile_path, outfile_path, chrom_converter)
            
        fastahandler.make_index(outfile_path)


class GenesetProcessor(PostprocessorBase):
    pass


########################
# fasta postprocessing #
########################

def prepare_processed_fasta(
    standard_refver,
    edited_path,
    custom_path,
    RefSeq_refver,
    species,
    chrom_converter,
    unedited_file_checksum_path,
    unedited_file_wascustom_path,
    calc_checksum,
):
    fasta_processor = FastaProcessor(
        standard_refver=standard_refver,
        edited_path=edited_path,
        custom_path=custom_path,
        RefSeq_refver=RefSeq_refver,
        species=species,
        chrom_converter=chrom_converter,
        unedited_file_checksum_path=unedited_file_checksum_path,
        unedited_file_wascustom_path=unedited_file_wascustom_path,
        calc_checksum=calc_checksum,
    )
    fasta_processor.prepare_processed_file()

###########
# geneset #
###########

def prepare_processed_geneset(self):
    pass




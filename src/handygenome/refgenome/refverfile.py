import os
import pickle

import handygenome.logutils as logutils
import handygenome.tools as tools
import handygenome.refgenome.refgenome as refgenome
import handygenome.fastahandler as fastahandler
import handygenome.publicdb.ncbi as libncbi
import handygenome.publicdb.ncbi_cache as libncbicache
from handygenome.genomedf.genomedf import GenomeDataFrame as GDF


##############
# exceptions #
##############

class PARUnavailableError(Exception):
    pass


############
# N region #
############

def get_Nregion_path(refver):
    return os.path.join(refgenome.get_datadir(refver), 'Nregion.bed.gz')


def write_Nregion_file(refver, path):
    logutils.log(f'Writing a new N region file for refver {refver}')

    fasta = refgenome.get_fasta(refver)
    chroms, start0s, end0s = fastahandler.make_N_region_coords(fasta)
    Nregion_gdf = GDF.from_data(refver=refver, chroms=chroms, start0s=start0s, end0s=end0s)
    Nregion_gdf.write_tsv(path)
    logutils.log(f'Finished writing a new N region file for refver {refver}')


def get_Nregion_gdf(refver):
    Nregion_path = get_Nregion_path(refver)
    if not os.path.exists(Nregion_path):
        write_Nregion_file(refver, Nregion_path)
    return GDF.read_tsv(Nregion_path, refver)


#######
# PAR #
#######

def get_assembly_regions_df(refver, force_download=False):
    RefSeq_refver, species = refgenome.get_RefSeq_standard(refver)
    return libncbi.parse_assembly_regions(
        libncbicache.get_assembly_regions_path(
            RefSeq_refver, species, force_download=force_download,
        )
    )


def get_par_gdf(refver, force_download=False):
    df = get_assembly_regions_df(refver, force_download=force_download)
    par_subdf = df.loc[df['Scaffold-Role'] == 'PAR', :]
    if par_subdf.shape[0] == 0:
        raise PARUnavailableError(f'PAR is not available for this reference version.')

    chrom_converter = refgenome.get_chrom_converter(refver)
    chroms = [chrom_converter[x] for x in par_subdf['Chromosome']]

    # assume the coordinates are 1-based both closed
    start0s = par_subdf['Chromosome-Start'] - 1
    end0s = par_subdf['Chromosome-Stop']

    return GDF.from_data(
        refver=refver,
        chroms=chroms,
        start0s=start0s,
        end0s=end0s,
    )


###############
# GC fraction #
###############

def get_gcfraction_path(refver, binsize):
    refver_datadir = refgenome.get_datadir(refver)
    gcfile_path = os.path.join(refver_datadir, f'gcfraction_bin{binsize}.pickle')
    return gcfile_path


def write_gcfraction(refver, binsize):
    logutils.log(f'Writing a new GC fraction file for refver {refver} and binsize {binsize}')

    gcfile_path = get_gcfraction_path(refver, binsize)
    fasta = refgenome.get_fasta(refver)
    all_regions = GDF.all_regions(refver, assembled_only=False).window(binsize)
    all_regions.sort()

    chroms = all_regions['Chromosome']
    start0s = all_regions['Start']
    end0s = all_regions['End']
    gcs = fastahandler.calculate_gc_fractions(chroms, start0s, end0s, fasta, window=None, as_array=True)

    gc_gdf = GDF.from_data(
        refver=refver, chroms=chroms, start0s=start0s, end0s=end0s, GC=gcs,
    )

    with open(gcfile_path, 'wb') as outfile:
        pickle.dump(gc_gdf, outfile)
    
    
def load_gcfraction(refver, binsize):
    gcfile_path = get_gcfraction_path(refver, binsize)
    if not os.path.exists(gcfile_path):
        write_gcfraction(refver, binsize)

    with open(gcfile_path, 'rb') as infile:
        return pickle.load(infile)
    


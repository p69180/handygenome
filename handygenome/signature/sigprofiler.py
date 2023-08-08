import os
import tempfile
import shutil
import gzip
import contextlib
import logging
import functools

import numpy as np
import pandas as pd
import pysam

from handygenome.sigprofiler_clone.SigProfilerAssignment.single_sample import add_remove_signatures
from handygenome.sigprofiler_clone.SigProfilerMatrixGenerator.scripts.SigProfilerMatrixGeneratorFunc import SigProfilerMatrixGeneratorFunc

import handygenome.refgenome as refgenome
import handygenome.variant.varianthandler as varianthandler
import handygenome.vcfeditor.initvcf as initvcf
import handygenome.signature.misc as signature_misc
import handygenome.signature.signatureresult as signatureresult
import handygenome.deco as deco


LOGGER = logging.getLogger(__name__)


def make_uncompressed_vcf_copy(vcf_path, input_copy_path):
    if vcf_path.endswith('.gz'):
        with \
                gzip.open(vcf_path, 'rb') as in_file, \
                open(input_copy_path, 'wb') as out_file:
            shutil.copyfileobj(in_file, out_file)
    else:
        os.symlink(vcf_path, input_copy_path)


@deco.get_deco_arg_choices({'refver': signature_misc.AVAILABLE_REFVERS})
def get_catalogues(vcf_path, refver, verbose=True):
    sampleid = 'test'
    with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdir:
        # run SigProfilerMatrixGeneratorFunc
        input_copy_path = os.path.join(tmpdir, 'input.vcf')
        make_uncompressed_vcf_copy(vcf_path, input_copy_path)
        matgenfunc_partial = functools.partial(
            SigProfilerMatrixGeneratorFunc,
            project=sampleid, genome=refver, vcfFiles=tmpdir,
            exome=False, bed_file=None, chrom_based=False, 
            plot=False, tsb_stat=False)

        if verbose:
            matgenfunc_partial()
        else:
            with open(os.devnull, 'w') as devnull:
                with \
                        contextlib.redirect_stdout(devnull), \
                        contextlib.redirect_stderr(devnull):
                    matgenfunc_partial()

        # load results
        catalogues = dict()
        for muttype in ('SBS', 'ID', 'DBS'):
            subdir = os.path.join(tmpdir, 'output', muttype)
            for fname in os.listdir(subdir):
                key = fname.split('.')[1]  # e.g. SBS96, ID83
                key_lower = key.lower()
                catalogues[key] = pd.read_csv(os.path.join(subdir, fname), 
                                              sep='\t', index_col=0)
                catalogues[key_lower] = catalogues[key]

    return catalogues


@deco.get_deco_arg_choices({'refver': signature_misc.AVAILABLE_REFVERS})
def get_catalogues_from_vcfspecs(vcfspec_iter, refver, verbose=True):
    with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdir:
        # write a temporary vcf file
        input_vcf_path = os.path.join(tmpdir, 'input.vcf')
        chromdict = refgenome.ChromDict.from_refver(refver)
        header = initvcf.create_header(chromdict)
        with pysam.VariantFile(input_vcf_path, mode='w0', 
                               header=header) as out_vcf:
            for vcfspec in vcfspec_iter:
                vr = header.new_record()
                varianthandler.apply_vcfspec(vr, vcfspec)
                out_vcf.write(vr)
        # create catalogues 
        catalogues = get_catalogues(input_vcf_path, refver, verbose=verbose)

    return catalogues


#############################################


def turn_signames_into_indexes(signames, sigdata):
    result = list()
    sigdata_signames = tuple(sigdata.columns)
    for x in signames:
        if x in sigdata_signames:
            result.append(sigdata_signames.index(x))
        else:
            LOGGER.warning(f'Function "run_assignment": '
                           f'Signature name {x} is not included in the '
                           f'input signature dataset.')

    return result


def add_recommended_sigs(sigdata, background_sigs, permanent_sigs):
    recommended_sigs = ['SBS1', 'SBS5']
    if set(recommended_sigs).issubset(set(sigdata.columns)):
        background_sigs.extend(recommended_sigs)
        permanent_sigs.extend(recommended_sigs)


def run_assignment_sanity_check(catalogue, sigdata):
    if set(catalogue.index) != set(sigdata.index):
        raise Exception(f'The mutation catalogue type of "catalogue" and '
                        f'"sigdata" are different.')


def reorder_catalogue(catalogue, sigdata):
    if tuple(catalogue.index) != tuple(sigdata.index):
        catalogue_reorderd = catalogue.reindex(sigdata.index)
        return catalogue_reordered
    else:
        return catalogue


def run_assignment(catalogue: pd.Series,
                   sigdata: pd.DataFrame,
                   background_sigs=list(), 
                   permanent_sigs=list(),
                   checkrule_sigs=list(), 
                   add_penalty=0.05, 
                   remove_penalty=0.01,
                   checkrule_penalty=1.00,
                   use_recommended_sigs=True,
                   use_connected_sigs=True, 
                   verbose=False):
    """Args:
        catalogue: pandas.Series representing the mutation catalogue of 
            the sample. Its index must indicate mutation contexts 
            (e.g. C[C>T]A) and each value must represent the mutation count.
        sigdata: pandas.DataFrame where rows indicate mutation contexts 
            (e.g. C[C>T]A) and columns indicate signature components.

            If must be that "set(catalogue.index) == set(sigdata.index)".
            If the order is different, "catalogue" is reordered by 
            "pd.Series.reindex" method.

        background_sigs: These signatures are included at the beginning of
            the add-remove cycles.
        permanent_sigs: These signatures are included at the beginning of
            the add-remove cycles. In addition, these are not removed during
            remove cycles.
        add_penalty: A signature component can be added if similarity gain is
            at least this value.
        checkrule_sigs: These signatures get penalty during addition cycles.
        checkrule_penalty: This value is multiplied to the raw similarity
            during checkrule penalty application.
        use_connected_sigs: If True, "background_sigs" is augmented with related 
            signatures as follows:
                ["SBS2", "SBS13"]
                ["SBS7a", "SBS7b", "SBS7c", "SBS7d"]
                ["SBS10a", "SBS10b"]
                ["SBS17a", "SBS17b"]

    Raises:
        When set(catalogue.index) != set(sigdata.index)
    """
    # sanity check
    run_assignment_sanity_check(catalogue, sigdata)
    # reorder "catalogue" if its catalogue key order is different with sigdata
    catalogue = reorder_catalogue(catalogue, sigdata)
    # apply recommended background_sigs and permanent_sigs
    if use_recommended_sigs:
        add_recommended_sigs(sigdata, background_sigs, permanent_sigs)
    # make a temporary file for logging
    fd, logfile_path = tempfile.mkstemp(dir=os.getcwd(), 
                                        prefix='SigProfiler_log.', 
                                        suffix='.txt')
    os.close(fd)
    # convert pd.DataFrame into np.ndarray
    catalogue_array = np.array(catalogue)
    sigdata_array = np.array(sigdata)
    # convert arguments
    background_sigs_indexes = turn_signames_into_indexes(background_sigs,
                                                         sigdata)
    permanent_sigs_indexes = turn_signames_into_indexes(permanent_sigs,
                                                        sigdata)
    checkrule_sigs_indexes = turn_signames_into_indexes(checkrule_sigs,
                                                        sigdata)
    # run main function
    (background_sigs, 
     finalactivities, 
     original_distance, 
     cosine_similarity, 
     kldiv, 
     correlation, 
     cosine_similarity_with_four_signatures) = add_remove_signatures(
        W=sigdata_array, sample=catalogue_array,
        background_sigs=background_sigs_indexes,
        permanent_sigs=permanent_sigs_indexes,
        add_penalty=add_penalty, remove_penalty=remove_penalty,
        check_rule_negatives=checkrule_sigs_indexes,
        checkrule_penalty=checkrule_penalty,
        connected_sigs=use_connected_sigs, verbose=verbose)

    # remove temporary log file
    os.remove(logfile_path)

    # create sigresult
    exposure = pd.Series(dict(zip(sigdata.columns, finalactivities)))
    sigresult = signatureresult.SignatureResult(catalogue=catalogue,
                                                sigdata=sigdata,
                                                exposure=exposure,
                                                cossim=cosine_similarity,
                                                kldiv=kldiv,
                                                correlation=correlation)

    return sigresult




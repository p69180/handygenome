import pysam
import pandas as pd

import handygenome.workflow as workflow
import handygenome.variant.varianthandler as varianthandler
import handygenome.signature.misc as signature_misc
import handygenome.signature.sigprofiler as sigprofiler
import handygenome.signature.cataloguing as cataloguing
import handygenome.signature.plotter_sbs96 as plotter_sbs96
import handygenome.variant.vcfspec as libvcfspec
import handygenome.deco as deco


LOGGER = workflow.get_logger(__name__, level='info')


class SignatureResult:
    def __init__(self, 
                 catalogue: pd.Series, 
                 sigdata: pd.DataFrame, 
                 exposure: pd.Series, 
                 cossim,
                 kldiv=None, correlation=None):
        def sanity_check(catalogue, sigdata, exposure):
            if tuple(catalogue.index) != tuple(sigdata.index):
                raise Exception(
                    f'Index of "catalogue" and index of "sigdata" '
                    f'are different')
            if tuple(exposure.index) != tuple(sigdata.columns):
                raise Exception(
                    f'Index of "exposure" and columns of "sigdata" '
                    f'are different')

        sanity_check(catalogue, sigdata, exposure)

        self.catalogue = catalogue
        self.sigdata = sigdata
        self.exposure = exposure
        self.cossim = cossim
        self.kldiv = kldiv
        self.correlation = correlation

        self.catalogue_type = signature_misc.get_catalogue_type(self.catalogue.index)
        self.nonzero_sigs = tuple(
            self.exposure.index[self.exposure > 0])

    def plot(self, sampleid=None):
        if self.catalogue_type == 'sbs96':
            plotter_sbs96.main(self, sampleid=sampleid)
        else:
            raise Exception(f'Unavailable catalogue type')


@deco.get_deco_arg_choices({'refver': signature_misc.AVAILABLE_REFVERS})
@deco.get_deco_arg_choices({'catalogue_type': ('sbs96', 'id83')})
@deco.get_deco_arg_choices({'cataloguer': ('custom', 'sigprofiler')})
def get_sigresult_from_vcfspecs(vcfspec_iter, refver='GRCh37', 
                                catalogue_type='sbs96',
                                cataloguer='custom',
                                verbose=True,
                                verbose_matgen=True,
                                verbose_assign=False, 
                                background_sigs=list(), 
                                permanent_sigs=list(),
                                checkrule_sigs=list(), 
                                add_penalty=0.05, 
                                remove_penalty=0.01,
                                checkrule_penalty=1.00,
                                use_recommended_sigs=True,
                                use_connected_sigs=True):
    # get catalogue
    if verbose:
        LOGGER.info('Creating mutation catalogue')

    if cataloguer == 'custom':
        if catalogue_type == 'sbs96':
            catalogue = cataloguing.get_sbs96_catalogue_vcfspecs(vcfspec_iter, 
                                                                 refver)
        else:
            raise Exception(f'Custom cataloguer is available only for sbs96.')
    elif cataloguer == 'sigprofiler':
        catalogues = sigprofiler.get_catalogues_from_vcfspecs(
            vcfspec_iter, refver, verbose=verbose_matgen)
        catalogue = catalogues[catalogue_type].iloc[:, 0]

    # loading signature database
    sigdata = signature_misc.load_signature_data(refver, catalogue_type)

    # fitting signature components
    if verbose:
        LOGGER.info('Running SigProfilerAssignment')

    sigresult = sigprofiler.run_assignment(
        catalogue, sigdata,
        background_sigs=background_sigs, 
        permanent_sigs=permanent_sigs,
        checkrule_sigs=checkrule_sigs, 
        add_penalty=add_penalty, 
        remove_penalty=remove_penalty,
        checkrule_penalty=checkrule_penalty,
        use_recommended_sigs=use_recommended_sigs,
        use_connected_sigs=use_connected_sigs,
        verbose=verbose_assign)

    if verbose:
        LOGGER.info('All finished.')

    return sigresult


@deco.get_deco_arg_choices({'refver': signature_misc.AVAILABLE_REFVERS})
def get_sigresult_from_vcfpath(vcf_path, refver='GRCh37', **kwargs):
    vcfspec_iter = (
        libvcfspec.Vcfspec.from_vr(vr)
        #varianthandler.get_vcfspec(vr)
        for vr in pysam.VariantFile(vcf_path).fetch()
    )
    sigresult = get_sigresult_from_vcfspecs(vcfspec_iter, refver, **kwargs)

    return sigresult


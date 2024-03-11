import re

import pysam
import pandas as pd

import handygenome.variant.varianthandler as varianthandler
import handygenome.signature.misc as signature_misc
import handygenome.signature.sigprofiler as sigprofiler
import handygenome.signature.cataloguing as cataloguing
#import handygenome.signature.plotter_sbs96 as plotter_sbs96
import handygenome.signature.plot as sigplot
from handygenome.variant.vcfspec import Vcfspec
import handygenome.deco as deco
import handygenome.logutils as logutils
import handygenome.signature.misc as sigmisc


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

    def get_catalogue_sbs6(self):
        return sigplot.make_sbs6_dict_from_sbs96(self.catalogue)

    def plot(self, **kwargs):
        if self.catalogue_type == 'sbs96':
            return sigplot.draw_onesample(self, **kwargs)
        else:
            raise Exception(f'Unavailable catalogue type')

    def get_artefact_burden(self):
        exposure = sum(
            (self.exposure[key] if (key in self.exposure) else 0)
            for key in sigmisc.ARTEFACTS
        )
        exp_sum = sum(self.exposure)
        fraction = exposure / exp_sum
        return {'exposure': exposure, 'fraction': fraction}


#@deco.get_deco_arg_choices({'refver': signature_misc.AVAILABLE_REFVERS})
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
        logutils.log(f'Creating mutation catalogue', level='info')

    if cataloguer == 'custom':
        if catalogue_type == 'sbs96':
            catalogue = cataloguing.get_sbs96_catalogue_vcfspecs(vcfspec_iter, refver)
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
        logutils.log(f'Running SigProfilerAssignment', level='info')

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
        logutils.log(f'All finished.', level='info')

    return sigresult


#@deco.get_deco_arg_choices({'refver': signature_misc.AVAILABLE_REFVERS})
def get_sigresult_from_vcfpath(vcf_path, refver='GRCh37', **kwargs):
    vcfspec_iter = (
        Vcfspec.from_vr(vr)
        for vr in pysam.VariantFile(vcf_path).fetch()
    )
    sigresult = get_sigresult_from_vcfspecs(vcfspec_iter, refver, **kwargs)

    return sigresult


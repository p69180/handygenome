import re
import itertools

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
hgvs = importlib.import_module('.'.join([top_package_name, 'hgvs']))
infoformat = importlib.import_module('.'.join([top_package_name, 'variant', 'infoformat']))
ensembl_rest = importlib.import_module('.'.join([top_package_name, 'annotation', 'ensembl_rest']))
#annotationdb = importlib.import_module('.'.join([top_package_name, 'annotation', 'annotationdb']))
veplib = importlib.import_module('.'.join([top_package_name, 'annotation', 'veplib']))
ensembl_feature = importlib.import_module('.'.join([top_package_name, 'annotation', 'ensembl_feature']))


VEP_INFOMETA_PAT = re.compile('Consequence annotations from Ensembl VEP. Format: ([^|]+(\|[^|]+)*)')

CMDLINE_VEP_POLYPHEN_SIFT_PAT = re.compile('(.+)\((.+)\)')

REGULATORY_FEATURE_TYPES = {
    'regulatory': {
        'Promoter': 'promoter',
        'Promoter Flanking Region': 'promoter_flank',
        'CTCF Binding Site': 'CTCFBS',
        'TF binding site': 'TFBS',
        'Enhancer': 'enhancer',
        'Open chromatin': 'open_chromatin',
        },
    'overlap': {
        'Predicted promoter': 'promoter',
        'Predicted promoter flanking region': 'promoter_flank',
        'CTCF binding site': 'CTCFBS',
        'Transcription factor binding site': 'TFBS',
        'Predicted enhancer region': 'enhancer',
        'Open chromatin region': 'open_chromatin',
        },
    'vep': {
        'promoter': 'promoter',
        'promoter_flanking_region': 'promoter_flank',
        'CTCF_binding_site': 'CTCFBS',
        'TF_binding_site': 'TFBS',
        'enhancer': 'enhancer',
        'open_chromatin_region': 'open_chromatin',
        },
    }

# reference: https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html
CONSEQUENCE_TRANSCRIPT_CLASSES = {
    'is_splice_region_involved': (
        'splice_acceptor_variant',
        'splice_donor_variant',
        'splice_donor_region_variant',  # added on 220422
        'splice_region_variant',
        'splice_polypyrimidine_tract_variant', # added on 220329
        'splice_donor_5th_base_variant', # added on 220401
        ),
    'is_splice_acceptor_involved': (
        'splice_acceptor_variant',
        ),
    'is_splice_donor_involved': (
        'splice_donor_variant',
        'splice_donor_region_variant',
        'splice_donor_5th_base_variant', # added on 220401
        ),

    'is_5pUTR_involved': (
        '5_prime_UTR_variant',
        ),

    'is_3pUTR_involved': (
        '3_prime_UTR_variant',
        ),

    'is_protein_altering': (
        'stop_gained',
        'frameshift_variant',
        'stop_lost',
        'start_lost',
        'inframe_insertion',
        'inframe_deletion',
        'missense_variant',
        'protein_altering_variant',
        ),
    'is_not_protein_altering': (
        'splice_acceptor_variant',
        'splice_donor_variant',
        'splice_region_variant',
        'start_retained_variant',
        'stop_retained_variant',
        'synonymous_variant',
        'mature_miRNA_variant',
        '5_prime_UTR_variant',
        '3_prime_UTR_variant',
        'non_coding_transcript_exon_variant',
        'intron_variant',
        'NMD_transcript_variant',
        'non_coding_transcript_variant',
        'upstream_gene_variant',
        'downstream_gene_variant',
        'TF_binding_site_variant',
        'regulatory_region_variant',
        'intergenic_variant',
        ),

    'is_synonymous': (
        'synonymous_variant',
        ),
    'is_missense': (
        'missense_variant',
        ),
    'is_frameshift': (
        'frameshift_variant',
        ),
    'is_inframe': (
        'inframe_deletion',
        ),
    'is_stopgain': (
        'stop_gained',
        ),
    'is_stoplost': (
        'stop_lost',
        ),
    'is_startlost': (
        'start_lost',
        ),

    'is_unclassified': (
        'incomplete_terminal_codon_variant',
        'coding_sequence_variant',
        ),
    'is_SV_consequence': (
        'transcript_ablation',
        'transcript_amplification',
        'TFBS_ablation',
        'TFBS_amplification',
        'regulatory_region_ablation',
        'regulatory_region_amplification',
        'feature_elongation',
        'feature_truncation',
        ),
    }


CONSEQUENCE_TRANSCRIPT_VALUES = set(
    itertools.chain.from_iterable(CONSEQUENCE_TRANSCRIPT_CLASSES.values()))


BIOTYPE_TRANSCRIPT_CLASSES = {
    # coding
    'coding': (
        'IG_C_gene',
        'IG_D_gene',
        'IG_J_gene',
        'IG_V_gene',
        'ccds_gene',
        'mRNA',
        'protein_coding',
        ),
    'IG': (
        'IG_C_gene',
        'IG_D_gene',
        'IG_J_gene',
        'IG_V_gene',
       ),
    'TR': (
        'TR_C_gene',
        'TR_D_gene',
        'TR_J_gene',
        'TR_V_gene',
       ),

    # noncoding
    'noncoding': (
        '3prime_overlapping_ncRNA',  # updated on 220421; ENST00000606379
        '3prime_overlapping_ncrna',
        'Mt_rRNA',
        'Mt_tRNA',
        'RNase_MRP_RNA',
        'RNase_P_RNA',
        'Y_RNA',
        'antisense_RNA',
        'antisense',  # updated on 220316; ENST00000442411
        'lncRNA',
        'lincRNA',
        'miRNA',
        'misc_RNA',
        'non_stop_decay',
        'nonsense_mediated_decay',
        'processed_transcript',
        'rRNA',
        'retained_intron',
        'ribozyme',
        'sRNA',
        'scRNA',
        'scaRNA',
        'snRNA',
        'snoRNA',
        'tRNA',
        'telomerase_RNA',
        'vault_RNA',
        'sense_overlapping',
        'sense_intronic',  # updated on 220421; ENST00000428191
        ),
    'rRNA': (
        'Mt_rRNA',
        'rRNA',
        ),
    'tRNA': (
        'Mt_tRNA',
        'tRNA',
        ),
    'miRNA': (
        'miRNA',
        ),
    'lncRNA': (
        'lncRNA',
        'lincRNA',
        ),
    'NMD': (
        'nonsense_mediated_decay',
        ),
    'NSD': (
        'non_stop_decay',
        ),
    'antisense': (
        'antisense_RNA',
        'antisense',
        ),

    # pseudogenes
    'pseudogene': (
        'IG_C_pseudogene',
        'IG_J_pseudogene',
        'IG_V_pseudogene',
        'IG_pseudogene',
        'TR_J_pseudogene',
        'TR_V_pseudogene',
        'ncRNA_pseudogene',
        'polymorphic_pseudogene',
        'processed_pseudogene',
        'pseudogene',
        'rRNA_pseudogene',
        'transcribed_processed_pseudogene',
        'transcribed_pseudogene',
        'transcribed_unitary_pseudogene',
        'transcribed_unprocessed_pseudogene',
        'translated_processed_pseudogene',
        'translated_unprocessed_pseudogene',
        'unitary_pseudogene',
        'unprocessed_pseudogene',
        ),
    'IG_pseudogene': (
        'IG_C_pseudogene',
        'IG_J_pseudogene',
        'IG_V_pseudogene',
        'IG_pseudogene',
        ),
    'TR_pseudogene': (
        'TR_J_pseudogene',
        'TR_V_pseudogene',
        ),
    'processed_pseudogene': (
        'processed_pseudogene',
        'transcribed_processed_pseudogene',
        'transcribed_unprocessed_pseudogene',
        'translated_processed_pseudogene',
        'translated_unprocessed_pseudogene',
        'unprocessed_pseudogene',
        ),
    'unprocessed_pseudogene': (
        'transcribed_unprocessed_pseudogene',
        'translated_unprocessed_pseudogene',
        'unprocessed_pseudogene',
        ),
    'unitary_pseudogene': (
        'transcribed_unitary_pseudogene',
        'unitary_pseudogene',
        ),
    'translated_pseudogene': (
        'translated_processed_pseudogene',
        'translated_unprocessed_pseudogene',
        ),
    'transcribed_pseudogene': (
        'transcribed_processed_pseudogene',
        'transcribed_pseudogene',
        'transcribed_unitary_pseudogene',
        'transcribed_unprocessed_pseudogene',
        ),

    # unknown types
    'unknown_type': (
        'LRG_gene',
        'TEC',
        'aligned_transcript',
        'cdna_update',
        'guide_RNA',
        'other',
        ),
}
BIOTYPE_TRANSCRIPT_VALUES = set( itertools.chain.from_iterable( BIOTYPE_TRANSCRIPT_CLASSES.values() ) )


#######################################################################################


def set_feature_type(annotitem, feature_type):
    assert feature_type in ('transcript', 'regulatory', 'motif', 'repeat')

    if feature_type == 'transcript':
        annotitem['feature_type'] = 'transcript'
        '''
        annotitem['is_transcript'] = True
        annotitem['is_regulatory'] = False
        annotitem['is_motif'] = False
        annotitem['is_repeat'] = False
        '''

        annotitem['feature_type_flags'] = {
            'is_transcript': True,
            'is_regulatory': False,
            'is_motif': False,
            'is_repeat': False,
            }
    elif feature_type == 'regulatory':
        annotitem['feature_type'] = 'regulatory'
        '''
        annotitem['is_transcript'] = False
        annotitem['is_regulatory'] = True
        annotitem['is_motif'] = False
        annotitem['is_repeat'] = False
        '''

        annotitem['feature_type_flags'] = {
            'is_transcript': False,
            'is_regulatory': True,
            'is_motif': False,
            'is_repeat': False,
            }
    elif feature_type == 'motif':
        annotitem['feature_type'] = 'motif'
        '''
        annotitem['is_transcript'] = False
        annotitem['is_regulatory'] = False
        annotitem['is_motif'] = True
        annotitem['is_repeat'] = False
        '''

        annotitem['feature_type_flags'] = {
            'is_transcript': False,
            'is_regulatory': False,
            'is_motif': True,
            'is_repeat': False,
        }
    elif feature_type == 'repeat':
        annotitem['feature_type'] = 'repeat'
        '''
        annotitem['is_transcript'] = False
        annotitem['is_regulatory'] = False
        annotitem['is_motif'] = False
        annotitem['is_repeat'] = True
        '''

        annotitem['feature_type_flags'] = {
            'is_transcript': False,
            'is_regulatory': False,
            'is_motif': False,
            'is_repeat': True,
        }


def set_transcript_subtypes(annotitem, biotype, feature_id):
    """
    Args:
        annotitem: A dictionary
        biotype: Raw biotype string obtained from vep, rest, etc.
    """

    '''
    for key, val in BIOTYPE_TRANSCRIPT_CLASSES.items():
        annotitem[f'is_{key}'] = (biotype in val)
        '''

    if biotype not in BIOTYPE_TRANSCRIPT_VALUES:
        raise Exception(f'Unexpected biotype value; biotype: {biotype}; '
                        f'feature_id: {feature_id}')

    annotitem['transcript_subtype_flags'] = dict()
    for key, val in BIOTYPE_TRANSCRIPT_CLASSES.items():
        annotitem['transcript_subtype_flags'][f'is_{key}'] = (biotype in val)


def set_regulatory_subtypes(annotitem, raw_regulatory_type, option):
    """
    Args:
        annotitem: A dictionary
        raw_regulatory_type: Raw regulatory feature type string obtained 
            from vep, rest, etc.
    """

    assert option in REGULATORY_FEATURE_TYPES.keys()
    assert raw_regulatory_type in REGULATORY_FEATURE_TYPES[option].keys(), (
        f'Unexpected regulatory type "{raw_regulatory_type}", '
        f'option "{option}"')

    '''
    for val in REGULATORY_FEATURE_TYPES[option].values():
        annotitem[f'is_{val}'] = ( REGULATORY_FEATURE_TYPES[option][ raw_regulatory_type ] == val )
        '''

    annotitem['regulatory_subtype_flags'] = dict()
    for val in REGULATORY_FEATURE_TYPES[option].values():
        annotitem['regulatory_subtype_flags'][f'is_{val}'] = (
            REGULATORY_FEATURE_TYPES[option][raw_regulatory_type] == val)


def set_is_forward(annotitem, strand):
    if strand == 1:
        annotitem['is_forward'] = True
    elif strand == -1:
        annotitem['is_forward'] = False
    elif strand in (0, None):
        annotitem['is_forward'] = None
    else:
        raise Exception(f'Unexpected "strand" value: {strand}')


def set_exon_intron_attributes(transcript, exon_value, intron_value):
    def subfun(value):
        if value is None:
            segs = None
        else:
            value_sp = value.split('/')
            value_sp_sp = [int(x) for x in value_sp[0].split('-')]
            if len(value_sp_sp) == 1:
                segs = value_sp_sp
            else:
                segs = list(range(value_sp_sp[0], value_sp_sp[1]+1))

        return segs

    transcript['involved_exons'] = subfun(exon_value)
    transcript['involved_introns'] = subfun(intron_value)
    transcript['is_exon_involved'] = (exon_value is not None)
    transcript['is_intron_involved'] = (intron_value is not None)


def vep_exon_intron_handler(value):
    if value in (None, ''):
        segs = None
    else:
        value_sp = value.split('/')
        value_sp_sp = [int(x) for x in value_sp[0].split('-')]
        if len(value_sp_sp) == 1:
            segs = value_sp_sp
        else:
            segs = list(range(value_sp_sp[0], value_sp_sp[1]+1))

    return segs


def set_consequence_attributes(annotitem, consequences, feature_id):
    """
    Args:
        consequences: a sequence of consequence items
    """

    '''
    for key, val in CONSEQUENCE_TRANSCRIPT_CLASSES.items():
        annotitem[key] = ( len(set(consequences).intersection(val)) != 0 )
        '''

    consequences_set = set(consequences)
    if not consequences_set.issubset(CONSEQUENCE_TRANSCRIPT_VALUES):
        unexpected_consequences = consequences_set.difference(
            CONSEQUENCE_TRANSCRIPT_VALUES)
        raise Exception(f'Unexpected consequence value; '
                        f'consequence: {unexpected_consequences}; '
                        f'feature_id: {feature_id}')

    annotitem['consequence_flags'] = dict()
    for key, val in CONSEQUENCE_TRANSCRIPT_CLASSES.items():
        annotitem['consequence_flags'][key] = (
            len(consequences_set.intersection(val)) != 0)


def set_distance(annotitem, consequences, distance_value):
    assert 'is_forward' in annotitem
    assert distance_value is not None
    assert (
        len(
            {'upstream_gene_variant', 
             'downstream_gene_variant'}.intersection(consequences)) == 1)

    upstream = ('upstream_gene_variant' in consequences)
    downstream = ('downstream_gene_variant' in consequences)
    forward = annotitem['is_forward']

    if (forward and upstream) or ((not forward) and downstream):
        annotitem['distance'] = distance_value
    elif (forward and downstream) or ((not forward) and upstream):
        annotitem['distance'] = -1 * distance_value


def set_codon_frame0(annotitem):
    """0-based. Only for substitution."""

    if annotitem['codon_change'] is None:
        annotitem['codon_frame0'] = None
    else:
        codon_change = annotitem['codon_change']
        if len(codon_change[0]) == len(codon_change[1]):
            if len(codon_change[0]) == 3:
                annotitem['codon_frame0'] = [
                    x[0] for x in enumerate(codon_change[0]) 
                    if x[1].isupper()][0]
            else:
                if ('incomplete_terminal_codon_variant' 
                    in annotitem['consequences']):  # hg19 8:86354359 G>T
                    annotitem['codon_frame0'] = None
                else:
                    raise Exception(f'A feature with non-length-3 '
                                    f'"codon_change" value:\n{annotitem}')
        else:
            annotitem['codon_frame0'] = None


##########################################################


def get_VEPkeys(vcfheader):
    """
    Args:
        vcfheader: pysam.VariantHeader object. Source vcf must be a VEP output.
    
    Returns:
        A list composed of VEP annotation subfield names. An empty list 
            if VEP annotation header is absent.
    """

    if veplib.VEP_INFO_FIELD in vcfheader.info:
        mat = VEP_INFOMETA_PAT.fullmatch(
            vcfheader.info[veplib.VEP_INFO_FIELD].description)
        if mat is None:
            VEPkeys = None
        else:
            VEPkeys = mat.group(1).split('|')
    else:
        VEPkeys = None

    return VEPkeys


def extract_cmdline_vep_annotation(vr):
    if infoformat.check_NA_info(vr, veplib.VEP_INFO_FIELD):
        raw_result = None
    else:
        vepkeys = get_VEPkeys(vr.header)
        raw_result = list()
        for item in vr.info[veplib.VEP_INFO_FIELD]:
            raw_result_item = dict(
                zip(vepkeys, ((None if x == '' else common.str_to_nonstr(x))
                              for x in item.split('|'))))
            if raw_result_item['Feature_type'] is not None:
                raw_result.append(raw_result_item)

    return raw_result


def parse_cmdline_vep(vr):
    """Args:
        vr: pysam.VariantRecord object, derived from a vcf output of cmdline VEP.
    """

    def subfun_common_attributes(result, raw_result):
        def parse_clinvar(item, result):
            if 'VAR_SYNONYMS' in item:
                var_synonyms = dict() 
                for subitem in item['VAR_SYNONYMS'].split('--'):
                    subitem_split = subitem.split('::')
                    var_synonyms[subitem_split[0]] = (
                        subitem_split[1].split('&'))

                if 'ClinVar' in var_synonyms:
                    result['clinvar'] = var_synonyms['ClinVar']
                else:
                    result['clinvar'] = None
            else:
                result['clinvar'] = None

        def parse_pubmed(item, result):
            if 'PUBMED' in item:
                result['pubmed'] = [int(x) for x in item['PUBMED'].split('&')]
            else:
                result['pubmed'] = None

        pass
        '''
        item = raw_result[0]
        if pubmed:
            parse_pubmed(item, result)
        if clinvar:
            parse_clinvar(item, result)
        '''

    def subfun_set_polyphen_sift(transcript, raw_result_item):
        for option in ('polyphen', 'sift'):
            raw_value = (raw_result_item['PolyPhen'] 
                         if (option == 'polyphen') else 
                         raw_result_item['SIFT'])
            if raw_value is None:
                transcript[f'{option}_prediction'] = None
                transcript[f'{option}_score'] = None
            else:
                mat = CMDLINE_VEP_POLYPHEN_SIFT_PAT.fullmatch(raw_value)
                transcript[f'{option}_prediction'] = mat.group(1)
                transcript[f'{option}_score'] = float(mat.group(2))

    def subfun_set_distance(annotitem, raw_result_item, consequences):
        if raw_result_item['DISTANCE'] is None:
            annotitem['distance'] = None
        else:
            set_distance(annotitem, consequences, raw_result_item['DISTANCE'])

        #annotitem['is_non_overlapping'] = (annotitem['distance'] is not None)

    def subfun_transcript(raw_result_item):
        transcript = ensembl_feature.Transcript(is_missing=False)

        transcript['id'] = raw_result_item['Feature']
        transcript['biotype'] = raw_result_item['BIOTYPE']
        transcript['consequences'] = raw_result_item['Consequence'].split('&')
        transcript['gene_id'] = raw_result_item['Gene']
        transcript['gene_name'] = raw_result_item['SYMBOL']
        transcript['is_canonical'] = (raw_result_item['CANONICAL'] == 'YES')
        transcript['mane_select'] = raw_result_item['MANE_SELECT']
        transcript['refseq_id'] = raw_result_item['RefSeq']

        transcript['ccds_id'] = raw_result_item['CCDS']
        transcript['protein_id'] = raw_result_item['ENSP']
        transcript['aa_change'] = (
            None 
            if (raw_result_item['Amino_acids'] is None) else 
            raw_result_item['Amino_acids'].split('/'))
        transcript['codon_change'] = (
            None 
            if (raw_result_item['Codons'] is None) else 
            raw_result_item['Codons'].split('/'))

        #set_codon_frame0(transcript)

        transcript['variant_pos_transcript'] = raw_result_item['cDNA_position']
        transcript['variant_pos_cds'] = raw_result_item['CDS_position']
        transcript['variant_pos_protein'] = raw_result_item['Protein_position']

        transcript['hgvsc'] = raw_result_item['HGVSc']
        transcript['hgvsp'] = raw_result_item['HGVSp']

        set_exon_intron_attributes(transcript, 
                                   exon_value=raw_result_item['EXON'], 
                                   intron_value=raw_result_item['INTRON'])
        set_is_forward(transcript, raw_result_item['STRAND'])
        subfun_set_distance(transcript, raw_result_item, 
                            transcript['consequences'])
        subfun_set_polyphen_sift(transcript, raw_result_item)

        return transcript

    def subfun_regulatory(raw_result_item):
        regulatory = ensembl_feature.Regulatory.init_nonmissing()

        regulatory['id'] = raw_result_item['Feature']
        regulatory['biotype'] = raw_result_item['BIOTYPE']
        regulatory['consequences'] = raw_result_item['Consequence'].split('&')
        set_is_forward(regulatory, None)
        subfun_set_distance(regulatory, raw_result_item, 
                            regulatory['consequences'])

        return regulatory

    def subfun_motif(raw_result_item):
        motif = ensembl_feature.Motif.init_nonmissing()

        motif['id'] = raw_result_item['Feature']
        motif['consequences'] = raw_result_item['Consequence'].split('&')
        motif['matrix_id'] = raw_result_item['MOTIF_NAME']
        motif['TF'] = raw_result_item['TRANSCRIPTION_FACTORS'].split('&')

        set_is_forward(motif, raw_result_item['STRAND'])
        subfun_set_distance(motif, raw_result_item, 
                            motif['consequences'])

        return motif

    def filter_raw_result_item(raw_result_item):
        if raw_result_item['BIOTYPE'] is None:
            return False
        else:
            return True

    # main
    transcriptset = ensembl_feature.TranscriptSet.init_nonmissing()
    regulatoryset = ensembl_feature.RegulatorySet.init_nonmissing()
    motifset = ensembl_feature.MotifSet.init_nonmissing()

    raw_result = extract_cmdline_vep_annotation(vr)
    if raw_result is not None:
        #subfun_common_attributes(parsed, raw_result)
        for raw_result_item in raw_result:
            if not filter_raw_result_item(raw_result_item):
                continue

            featuretype = raw_result_item['Feature_type']
            assert featuretype in ('Transcript', 'RegulatoryFeature', 'MotifFeature'), (
                f'Unexpected "Feature_type":\n{raw_result_item}')

            if featuretype == 'Transcript':
                transcript = subfun_transcript(raw_result_item)
                transcriptset[transcript["id"]] = transcript
            elif featuretype == 'RegulatoryFeature':
                regulatory = subfun_regulatory(raw_result_item)
                regulatoryset[regulatory["id"]] = regulatory
            elif featuretype == 'MotifFeature':
                motif = subfun_motif(raw_result_item)
                motifset[motif["id"]] = motif

    parsed = {'transcript': transcriptset, 'regulatory': regulatoryset,
              'motif': motifset}

    return parsed


##################


def _parse_rest_lookup_transcript_singleitem(raw_result, refver, set_gene_name=True):
    def exon_parser(raw_data):
        result = list()
        for dic in raw_data:
            start0 = dic["start"] - 1
            end0 = dic["end"]
            result.append([start, end])
        return result

    transcript = ensembl_feature.Transcript.init_nonmissing()
    
    set_is_forward(transcript, raw_result['strand'])

    transcript['biotype'] = raw_result['biotype']
    transcript['id'] = raw_result['id']
    transcript['is_canonical'] = (raw_result['is_canonical'] == 1)

    transcript['chrom'] = raw_result['seq_region_name']
    transcript['start0'] = raw_result['start'] - 1
    transcript['start1'] = raw_result['start']
    transcript['end0'] = raw_result['end']
    transcript['end1'] = raw_result['end']

    if 'display_name' not in raw_result.keys():
        raise Exception(f'This rest-lookup transcript result does not have "display_name" key:\n{raw_result}')
    transcript['transcript_name'] = raw_result['display_name']

    transcript['gene_id'] = raw_result['Parent']
    if set_gene_name:
        lookup_gene_result = ensembl_rest.lookup_id(transcript['gene_id'], refver=refver, expand=False)
        if 'display_name' not in lookup_gene_result.keys():
            raise Exception(f'This rest-lookup gene result does not have "display_name" key:\n{lookup_gene_result}')
        transcript['gene_name'] = lookup_gene_result['display_name']

    if 'Exon' in raw_result:
        transcript['exon_borders'] = exon_parser(raw_result['Exon'])
    else:
        transcript['exon_borders'] = None

    return transcript


def parse_rest_lookup_transcript(raw_result, refver, set_gene_name=True):
    transcript = _parse_rest_lookup_transcript_singleitem(
        raw_result, refver=refver, set_gene_name=set_gene_name)

    transcriptset = ensembl_feature.TranscriptSet.init_nonmissing()
    transcriptset[transcript['id']] = transcript
    parsed = {'transcript': transcriptset}

    return parsed


def parse_rest_lookup_transcript_post(raw_result, refver, set_gene_name=True):
    transcriptset = ensembl_feature.TranscriptSet.init_nonmissing()
    for val in raw_result.values():
        transcript = _parse_rest_lookup_transcript_singleitem(
            val, refver=refver, set_gene_name=set_gene_name)
        transcriptset[transcript['id']] = transcript

    parsed = {'transcript': transcriptset}

    return parsed


##################


def parse_rest_regulatory(raw_result):
    regulatory = ensembl_feature.Regulatory.init_nonmissing()

    regulatory['id'] = raw_result['id']
    regulatory['chrom'] = raw_result['seq_region_name']
    regulatory['start0'] = raw_result['start'] - 1
    regulatory['start1'] = raw_result['start']
    regulatory['end0'] = raw_result['end']
    regulatory['end1'] = raw_result['end']
    regulatory['bound_start0'] = raw_result['bound_start'] - 1
    regulatory['bound_start1'] = raw_result['bound_start']
    regulatory['bound_end0'] = raw_result['bound_end']
    regulatory['bound_end1'] = raw_result['bound_end']
    regulatory['biotype'] = raw_result['feature_type']
    regulatory['activity'] = raw_result['activity']

    regulatoryset = ensembl_feature.RegulatorySet.init_nonmissing()
    regulatoryset[regulatory['id']] = regulatory

    parsed = {'regulatory': regulatoryset}
    return parsed


##################


def parse_rest_overlap(raw_result, refver, include_motif_without_evidence=False, set_gene_name=True):
    def subfun_common(annotitem, raw_result_item):
        annotitem['chrom'] = raw_result_item['seq_region_name']
        annotitem['start0'] = raw_result_item['start'] - 1
        annotitem['start1'] = raw_result_item['start']
        annotitem['end0'] = raw_result_item['end']
        annotitem['end1'] = raw_result_item['end']
        set_is_forward(annotitem, raw_result_item['strand'])

    def subfun_transcript(raw_result_item, refver, set_gene_name):
        transcript = ensembl_feature.Transcript.init_nonmissing()
        subfun_common(transcript, raw_result_item)

        transcript['biotype'] = raw_result_item['biotype']
        #result_item['exons'] = None
        transcript['id'] = raw_result_item['id']
        transcript['transcript_name'] = raw_result_item['external_name']
        transcript['gene_id'] = raw_result_item['Parent']
        transcript['is_canonical'] = (raw_result_item['is_canonical'] == 1)

        if set_gene_name:
            transcript['gene_name'] = ensembl_rest.lookup_id(
                transcript['gene_id'], refver=refver, expand=False
            )['display_name']

        if 'ccdsid' in raw_result_item:
            transcript['ccds_id'] = raw_result_item['ccdsid']
        else:
            transcript['ccds_id'] = None

        return transcript

    def subfun_regulatory(raw_result_item):
        regulatory = ensembl_feature.Regulatory.init_nonmissing()
        subfun_common(regulatory, raw_result_item)

        regulatory['id'] = raw_result_item['id']
        regulatory['bound_start0'] = raw_result_item['bound_start'] - 1
        regulatory['bound_start1'] = raw_result_item['bound_start']
        regulatory['bound_end0'] = raw_result_item['bound_end']
        regulatory['bound_end1'] = raw_result_item['bound_end']
        regulatory['biotype'] = raw_result_item['description']

        return regulatory

    def subfun_motif(raw_result_item):
        motif = ensembl_feature.Motif.init_nonmissing()
        subfun_common(motif, raw_result_item)

        if 'epigenomes_with_experimental_evidence' in raw_result_item:
            motif['is_with_evidence'] = True
            motif['evidence_name'] = '&'.join(
                raw_result_item['epigenomes_with_experimental_evidence'].split(',')
            )
        else:
            motif['is_with_evidence'] = False
            motif['evidence_name'] = None

        motif['id'] = raw_result_item['stable_id']
        motif['matrix_id'] = raw_result_item['binding_matrix_stable_id']
        motif['TF'] = raw_result_item['transcription_factor_complex'].split(',')

        return motif

    def subfun_repeat(raw_result_item):
        repeat = ensembl_feature.Repeat.init_nonmissing()
        subfun_common(repeat, raw_result_item)
        repeat['name'] = raw_result_item['description']

        return repeat

    # main
    transcriptset = ensembl_feature.TranscriptSet.init_nonmissing()
    regulatoryset = ensembl_feature.RegulatorySet.init_nonmissing()
    motifset = ensembl_feature.MotifSet.init_nonmissing()
    repeatset = ensembl_feature.RepeatSet.init_nonmissing()

    for raw_result_item in raw_result:
        feature_type = raw_result_item['feature_type']
        if feature_type == 'transcript':
            transcript = subfun_transcript(raw_result_item, refver, set_gene_name)
            transcriptset[transcript["id"]] = transcript
        elif feature_type == 'regulatory':
            regulatory = subfun_regulatory(raw_result_item)
            regulatoryset[regulatory["id"]] = regulatory
        elif feature_type == 'motif':
            motif = subfun_motif(raw_result_item)
            motifset[motif["id"]] = motif
        elif feature_type == 'repeat':
            repeat = subfun_repeat(raw_result_item)
            repeatset[repeat["id"]] = repeat
        else:
            raise Exception(f'Unexpected feature type:\n{raw_result_item}')

    parsed = {
        'transcript': transcriptset,
        'regulatory': regulatoryset,
        'motif': motifset,
        'repeat': repeatset,
    }

    return parsed


##################


def parse_rest_vep(raw_result):
    def subfun_common(result, raw_result, pubmed, clinvar):
        pass
        '''
        null = False

        if 'colocated_variants' in raw_result:
            relevant_item = list(filter(
                        lambda x: ('frequencies' in x), 
                        raw_result['colocated_variants'],
                        ))

            assert len(relevant_item) in (0, 1)
            if len(relevant_item) == 0:
                null = True
            else:
                relevant_item = relevant_item[0]

                # pubmed, clinvar
                if pubmed:
                    result['pubmed'] = relevant_item['pubmed'] if ('pubmed' in relevant_item) else None

                if clinvar:
                    try:
                        result['clinvar'] = relevant_item['var_synonyms']['ClinVar']
                    except KeyError:
                        result['clinvar'] = None

                # population frequencies
                if relevant_item['minor_allele'] == alt:
                    result['AF_1kg'] = relevant_item['minor_allele_freq']
                else:
                    result['AF_1kg'] = None

                if alt in relevant_item['frequencies']:
                    result['AF_gnomAD_exome'] = relevant_item['frequencies'][alt]['gnomad']
                else:
                    result['AF_gnomAD_exome'] = None
        else:
            null = True

        if null:
            if pubmed:
                result['pubmed'] = None
            if clinvar:
                result['clinvar'] = None
            result['AF_1kg'] = None
            result['AF_gnomAD_exome'] = None
        '''

    def subfun_set_distance(annotitem, dic, consequences):
        if 'distance' in dic:
            set_distance(annotitem, consequences, dic['distance'])
        else:
            annotitem['distance'] = None

        #annotitem['is_non_overlapping'] = (annotitem['distance'] is not None)

    def subfun_set_polyphen_sift_cadd(annotitem, dic):
        for key in (
            'polyphen_prediction',
            'polyphen_score',
            'sift_prediction',
            'sift_score',
            'cadd_phred',
            'cadd_raw',
        ):
            if key in dic:
                annotitem[key] = dic[key]
            else:
                annotitem[key] = None

    def subfun_transcript(dic):
        transcript = ensembl_feature.Transcript.init_nonmissing()

        transcript['id'] = dic['transcript_id']
        transcript['biotype'] = dic['biotype']
        transcript['consequences'] = dic['consequence_terms']
        transcript['gene_id'] = dic['gene_id'] 
        transcript['gene_name'] = (dic['gene_symbol'] 
                                  if ('gene_symbol' in dic) else 
                                  None) # 'TEC' biotype transcript does not have a gene symbol
        transcript['is_canonical'] = ('canonical' in dic)
        transcript['mane_select'] = (dic['mane_select'] 
                                    if ('mane_select' in dic) else 
                                    None)

        transcript['ccds_id'] = dic['ccds'] if ('ccds' in dic) else None
        transcript['protein_id'] = (dic['protein_id'] 
                                   if ('protein_id' in dic) else 
                                   None)
        transcript['aa_change'] = (dic['amino_acids'].split('/') 
                                  if ('amino_acids' in dic) else 
                                  None)
        transcript['codon_change'] = (dic['codons'].split('/') 
                                     if ('codons' in dic) else 
                                     None)
        #set_codon_frame0(transcript)

        transcript['variant_start0_transcript'] = (dic['cdna_start'] - 1 
                                                  if ('cdna_start' in dic) 
                                                  else None)
        transcript['variant_start1_transcript'] = (dic['cdna_start'] 
                                                  if ('cdna_start' in dic) 
                                                  else None)
        transcript['variant_end0_transcript'] = (dic['cdna_end'] 
                                                if ('cdna_end' in dic) else 
                                                None)
        transcript['variant_end1_transcript'] = (dic['cdna_end'] 
                                                if ('cdna_end' in dic) else 
                                                None)

        transcript['variant_start0_cds'] = (dic['cds_start'] - 1 
                                           if ('cds_start' in dic) else 
                                           None)
        transcript['variant_start1_cds'] = (dic['cds_start'] 
                                           if ('cds_start' in dic) else 
                                           None)
        transcript['variant_end0_cds'] = (dic['cds_end'] 
                                         if ('cds_end' in dic) else 
                                         None)
        transcript['variant_end1_cds'] = (dic['cds_end'] 
                                         if ('cds_end' in dic) else 
                                         None)

        transcript['variant_start0_protein'] = (dic['protein_start'] - 1 
                                               if ('protein_start' in dic) else 
                                               None)
        transcript['variant_start1_protein'] = (
            dic['protein_start'] 
            if ('protein_start' in dic) else 
            None)
        transcript['variant_end0_protein'] = (
            dic['protein_end'] 
            if ('protein_end' in dic) else 
            None)
        transcript['variant_end1_protein'] = (
            dic['protein_end'] 
            if ('protein_end' in dic) else 
            None)

        transcript['hgvsc'] = dic['hgvsc'] if ('hgvsc' in dic) else None
        transcript['hgvsp'] = dic['hgvsp'] if ('hgvsp' in dic) else None

        set_exon_intron_attributes(
            transcript, 
            exon_value=(dic['exon'] if 'exon' in dic else None), 
            intron_value=(dic['intron'] if 'intron' in dic else None))
        set_is_forward(transcript, dic['strand'])
        subfun_set_distance(transcript, dic, transcript['consequences'])
        subfun_set_polyphen_sift_cadd(transcript, dic)

        return transcript

    def subfun_regulatory(dic):
        regulatory = ensembl_feature.Regulatory.init_nonmissing()

        regulatory['id'] = dic['regulatory_feature_id']
        regulatory['biotype'] = dic['biotype']
        regulatory['consequences'] = dic['consequence_terms']

        set_is_forward(regulatory, None)
        subfun_set_distance(regulatory, dic, regulatory['consequences'])

        return regulatory

    def subfun_motif(dic):
        motif = ensembl_feature.Motif.init_nonmissing()

        motif['id'] = dic['motif_feature_id']
        motif['consequences'] = dic['consequence_terms']
        motif['matrix_id'] = dic['motif_name']
        motif['TF'] = dic['transcription_factors']

        set_is_forward(motif, dic['strand'])
        subfun_set_distance(motif, dic, motif['consequences'])

        return motif

    # main
    transcriptset = ensembl_feature.TranscriptSet.init_nonmissing()
    regulatoryset = ensembl_feature.RegulatorySet.init_nonmissing()
    motifset = ensembl_feature.MotifSet.init_nonmissing()

    #raw_result = raw_result[0]
    # 'intergenic_consequences' is ignored
    if 'transcript_consequences' in raw_result:
        for dic in raw_result['transcript_consequences']:
            transcript = subfun_transcript(dic)
            transcriptset[transcript['id']] = transcript

    if 'regulatory_feature_consequences' in raw_result:
        for dic in raw_result['regulatory_feature_consequences']:
            regulatory = subfun_regulatory(dic)
            regulatoryset[regulatory['id']] = regulatory

    if 'motif_feature_consequences' in raw_result:
        for dic in raw_result['motif_feature_consequences']:
            motif = subfun_motif(dic)
            motifset[motif['id']] = motif

    parsed = {
        'transcript': transcriptset,
        'regulatory': regulatoryset,
        'motif': motifset,
    }

    return parsed


##############################################################


def parse_rest_map(raw_result, adjust_chrom=False, fasta=None):
    assert len(raw_result['mappings']) == 1
    dic = raw_result['mappings'][0]

    if adjust_chrom:
        assert fasta is not None
        chrom = hgvs.modify_chrom(dic['seq_region_name'], fasta)
    else:
        chrom = dic['seq_region_name']

    start1 = dic['start']
    end1 = dic['end']
    is_forward = (dic['strand'] == 1)

    return chrom, start1, end1, is_forward

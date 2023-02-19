import os
import itertools
import pprint
import re

import collections

import handygenome.common as common


DATA_DIR = os.path.join(common.DATA_DIR, 'assembly_reports')
if not os.path.exists(DATA_DIR):
    os.mkdir(DATA_DIR)


ASSEMBLYFILE_URLS = common.RefverDict({
    'NCBI36': 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.12_NCBI36/GCF_000001405.12_NCBI36_assembly_report.txt',
    'GRCh37': 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt',
    'GRCh38': 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt',

    'GRCm38': 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/all_assembly_versions/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_assembly_report.txt',
    'GRCm39': 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/all_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_assembly_report.txt',
})


PATHS = {
    key: os.path.join(DATA_DIR, os.path.basename(val))
    for (key, val) in ASSEMBLYFILE_URLS.items()
}


#INVALID_NAMES = [ 'NT_187507.1' ]

ASSEMBLY_REPORT_KEYDICT = {
    'default': 'Sequence-Name',
    'ucsc': 'UCSC-style-name',
    'genbank': 'GenBank-Accn',
    'refseq': 'RefSeq-Accn',
    'length': 'Sequence-Length',
    'role': 'Sequence-Role',
}


class AssemblySpec:
    """Attributes:
        data: { 
            'length': [249250621, 243199373, ...],
            'role': ['assembled-molecule', ... 'unlocalized-scaffold', ... ],
            'default': ['1', '2', ... 'HSCHR1_RANDOM_CTG5', ... ],
            'ucsc': ['chr1', 'chr2', ... 'chr1_gl000191_random', ... ],
            'nochr_plus_genbank': ['1', '2', ... 'GL000191.1', ... ],
            'genbank': ['CM000663.1', 'CM000664.1', ... ],
            'refseq': ['NC_000001.11', 'NC_000002.12', ... ],
        }
        dicts: {('ucsc', 'genbank'): dict, ... }
        chromdicts: {'ucsc': common.ChromDict object, ... }
        versions: All self.data keys except 'length' and 'role'
    """

    def __init__(self, data):
        """
        Args:
            data: {
                'ucsc': ['chr1', 'chr2', ... ], 
                'refseq': ['NC_000001.11', 'NC_000002.12', ... ],
                ...
                }
        """

        self.data = data
        self.dicts = dict()
        self.chromdicts = dict()

        for key1, key2 in itertools.permutations(data.keys(), 2):
            dic = dict()
            chromdict_data = {'contigs': list(), 'lengths': list()}
            for val1, val2 in zip(self.data[key1], self.data[key2]):
                if val1 is not None:
                    dic[val1] = val2
                    if key2 == 'length':
                        chromdict_data['contigs'].append(val1)
                        chromdict_data['lengths'].append(val2)

            self.dicts[(key1, key2)] = dic
            if key2 == 'length':
                self.chromdicts[key1] = common.ChromDict(custom=chromdict_data)

        self.versions = list(
            set(self.data.keys()).difference(('length', 'role')))

    def __repr__(self):
        tmp_keys = list()
        tmp_vals = list()
        for key, val in self.data.items():
            tmp_keys.append(key)
            tmp_vals.append(val)

        show_result = list()
        for tup in zip(*tmp_vals):
            show_result.append(str(dict(zip(tmp_keys, tup))))

        return '\n'.join(show_result)

    def convert(self, in_chrom, out_version):
        assert out_version in self.versions, (
            f'"out_version" argument must be one of {self.versions}')

        absent = True
        for version, names in self.data.items():
            if version in ('length', 'role'):
                continue
            if in_chrom in names:
                absent = False
                if version == out_version:
                    result = in_chrom
                else:
                    result = self.dicts[(version, out_version)][in_chrom]
                break

        if absent:
            raise Exception(f'Input chrom name "{in_chrom}" is not included '
                            f'in the database.')

        return result


def parse_assembly_report(assembly_report_path):
    data = {'length': list(), 'default': list(), 'ucsc': list(),
            'nochr_plus_genbank': list(), 'role': list(),
            'genbank': list(), 'refseq': list()}

    with open(assembly_report_path, 'r') as infile:
        record_start = False
        for line in infile:
            if not record_start:
                if line.startswith('# Sequence-Name'):
                    keys = common.get_linesp(re.sub('^#\s*', '', line))
                    record_start = True
                    continue
            else:
                linesp = common.get_linesp(line)
                linedict = dict(
                    zip(
                        keys, 
                        [(None if x == 'na' else x) for x in linesp]
                    )
                )

                default = linedict[ASSEMBLY_REPORT_KEYDICT['default']]
                ucsc = linedict[ASSEMBLY_REPORT_KEYDICT['ucsc']]
                genbank = linedict[ASSEMBLY_REPORT_KEYDICT['genbank']]
                refseq = linedict[ASSEMBLY_REPORT_KEYDICT['refseq']]
                length = int(linedict[ASSEMBLY_REPORT_KEYDICT['length']])
                role = linedict[ASSEMBLY_REPORT_KEYDICT['role']]

                #assert genbank is not None

                data['default'].append(default)
                data['ucsc'].append(ucsc)
                data['genbank'].append(genbank)
                data['refseq'].append(refseq)
                data['length'].append(length)
                data['role'].append(role)

                if ucsc is None:
                    data['nochr_plus_genbank'].append(genbank)
                elif ucsc == 'chrM':
                    data['nochr_plus_genbank'].append('MT')
                else:
                    mat = common.RE_PATS['assembled_chromosome'].fullmatch(ucsc)
                    if mat is None:
                        data['nochr_plus_genbank'].append(genbank)
                    else:
                        data['nochr_plus_genbank'].append(mat.group(2))

    return data


def get_assemblyspec_data(refver):
    assemblyfile_url = ASSEMBLYFILE_URLS[refver]
    assembly_report_path = PATHS[refver]
    if not os.path.exists(assembly_report_path):
        common.download(assemblyfile_url, assembly_report_path)

    return parse_assembly_report(assembly_report_path)
    


##############################################################################


SPECS = common.RefverDict(
    {
        refver: AssemblySpec(get_assemblyspec_data(refver))
        for refver in ASSEMBLYFILE_URLS.keys()
    }
)

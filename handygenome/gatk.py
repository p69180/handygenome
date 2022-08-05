import subprocess

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


def run_haplotypecaller(fasta_path, infile_path, outfile_path, tmpdir, 
                        incl_bed_path=None, excl_bed_path=None,
                        mem_gb=6, threads=2, mbq=10):
    """Args:
        threads: value for --native-pair-hmm-threads option. 
            gatk default is 4.
    """
    args = [
        common.GATK, 
        '--java-options', f'-Xmx{mem_gb}G',
        'HaplotypeCaller',
        '--reference', fasta_path,
        '--input', infile_path,
        '--output', outfile_path,
        '--native-pair-hmm-threads', str(threads),
        '--min-base-quality-score', str(mbq),
        '--sample-ploidy', '2', # gatk default
        '--active-probability-threshold', '0.002', # gatk default
        '--standard-min-confidence-threshold-for-calling', '30', # gatk default
        '--dont-use-soft-clipped-bases', 'false',
        '--create-output-variant-index', 'false',
        '--tmp-dir', tmpdir,
        ]

    if incl_bed_path is not None:
        args.extend(['-L', incl_bed_path])
    if excl_bed_path is not None:
        args.extend(['-XL', excl_bed_path])

    p = subprocess.run(args, capture_output=True, text=True, check=True) 
    return p

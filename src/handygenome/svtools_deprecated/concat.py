#!/home/users/pjh/conda_bin/python

import argparse
import pysam

from pprint import pprint


def argument_parser(_args):
	parser = argparse.ArgumentParser(
			description = '''
Receives input from multiple input files, then writes them into a single output file.
Input lines can be sorted with --sort option.
Output file metadata is created by union of all input vcfs.
'''
			) 
	parser.add_argument('-i', nargs = '*', help = 'One or more input vcf files. If not set, reads from stdin.')
	parser.add_argument('-o', required = True, help = 'Output vcf file path')
	parser.add_argument('-O', required = False, default = "z", help = 'Output vcf file format ("v", "z"(default), "u", "b")')
	parser.add_argument('--sort', required = False, help = 'Output sorting method. Default is no sorting, just concatenating input lines. Possible arguments: "pos", "mate"')
	parser.add_argument('--fasta', required = False, help = 'Reference fasta path. Required for sorting.')

	args = parser.parse_args(_args)
	sanity_check_args(args)
	modify_args(args)

	return args


def sanity_check_args(args):
	if args.O not in 'vzub':
		raise Exception('-O option must be one of "v", "z", "u", "b".')
	if args.sort not in (None, 'pos', 'mate'):
		raise Exception('--sort option must be one of "pos" or "mate".')
	if args.sort != None and args.fasta == None:
		raise Exception('--fasta must be set if --sort is used.')


def modify_args(args):
	if args.i == None:
		args.i = ['/dev/stdin']
	args.O = {'v':'wu', 'z':'w', 'u':'wbu', 'b':'wb'}[args.O]


###########################################################

def coord_sortkey(chrom, pos, fasta): 
	'''
	fasta: pysam.FastaFile object
	pos: 1-based
	'''
	return sum(fasta.lengths[ : fasta.references.index(chrom)]) + pos

###########################################################

def check_having_MATEID(args):
	for input_vcf_path in args.i:
		with pysam.VariantFile(input_vcf_path, 'r') as input_vcf:
			for pysamvr in input_vcf.fetch():
				if 'MATEID' not in pysamvr.info.keys():
					return False
	return True


def check_sample_columns(args):
	if len(args.i) == 1:
		return True
	elif len(args.i) > 1:
		samples_list = list()
		for input_vcf_path in args.i:
			with pysam.VariantFile(input_vcf_path, 'r') as input_vcf:
				samples_list.append(list(input_vcf.header.samples))
		for samples in samples_list[1:]:
			if samples != samples_list[0]:
				return False
		return True


def sanity_check(args):
	if not check_sample_columns(args):
		raise Exception('Samples of input vcfs are not identical')
	if args.sort == 'mate':
		if not check_having_MATEID(args):
			raise Exception('"--sort mate" used but one or more input vcf file do not have INFO/MATEID.')


###########################################################

def get_fasta(args):
	if args.fasta == None:
		return None
	else:
		return pysam.FastaFile(args.fasta)


def create_new_header(vcfver = "4.3", fasta = None):
	'''
	chrom_list : list of tuples (chrom, length)
	'''
	import tempfile
	import os

	tf = tempfile.NamedTemporaryFile(mode = 'w', delete = False)
	tf.write(f'##fileformat=VCFv{vcfver}\n')
	if not isinstance(fasta, type(None)):
		for chrom, length in zip(fasta.references, fasta.lengths):
			tf.write(f'##contig=<ID={chrom},length={length}>\n')
	tf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
	tf.close()

	dummy_vcf = pysam.VariantFile(tf.name)
	header = dummy_vcf.header
	dummy_vcf.close()

	os.remove(tf.name)

	return header


def get_merged_metadata_header(args, fasta):
	#records_list = [x.header.records for x in input_vcf_list]
	keys = ('GENERIC', 'CONTIG', 'FILTER', 'INFO', 'FORMAT', 'STRUCTURED')

	record_dict = dict()
	for key in keys:
		record_dict[key] = dict()
		record_dict[key]['keys'] = list()
		record_dict[key]['values'] = list()
	with pysam.VariantFile(args.i[0], 'r') as input_vcf:
		input_sample_list = list(input_vcf.header.samples)

	for input_vcf_path in args.i:
		with pysam.VariantFile(input_vcf_path, 'r') as input_vcf:
			for record in input_vcf.header.records:
				if record.type == 'GENERIC' and record.key == 'fileformat':
					pass
				if record.type == 'GENERIC':
					if record.key not in record_dict[record.type]['keys']:
						record_dict[record.type]['keys'].append(record.key)
						record_dict[record.type]['values'].append(record)
				else:
					if record.get('ID') not in record_dict[record.type]['keys']:
						record_dict[record.type]['keys'].append(record.get('ID'))
						record_dict[record.type]['values'].append(record)

	header = create_new_header(fasta = fasta)
	for key in keys:
		for record in record_dict[key]['values']:
			header.add_record(record)
	for sample in input_sample_list:
		header.add_sample(sample)

	return header


def write_nosort(args, header):
	with pysam.VariantFile(args.o, mode = args.O, header = header) as out_vcf:
		for input_vcf_path in args.i:
			with pysam.VariantFile(input_vcf_path, 'r') as input_vcf:
				for pysamvr in input_vcf.fetch():
					out_vcf.write(pysamvr)


def get_pysamvr_list(args):
	pysamvr_list = list()
	for input_vcf_path in args.i:
		with pysam.VariantFile(input_vcf_path, 'r') as input_vcf:
			pysamvr_list.extend(list(input_vcf.fetch()))

	return pysamvr_list


def write_possort(args, header, fasta):
	with pysam.VariantFile(args.o, mode = args.O, header = header) as out_vcf:
		pysamvr_list = get_pysamvr_list(args)
		pysamvr_list.sort(key = lambda x: coord_sortkey(x.contig, x.pos, fasta))
		for pysamvr in pysamvr_list:
			out_vcf.write(pysamvr)


def write_matesort(args, header, fasta):
	with pysam.VariantFile(args.o, mode = args.O, header = header) as out_vcf:
		pysamvr_dict, mateid_tuple_list = get_mateid_spec(args, fasta)
		for id1, id2 in mateid_tuple_list:
			out_vcf.write(pysamvr_dict[id1])
			out_vcf.write(pysamvr_dict[id2])


def get_mateid_spec(args, fasta):
	'''
	mateid_tuple_list : pysamvr ID tuples (sorted as bnd1, bnd2) are sorted by coordinate of bnd1.
	'''
	pysamvr_list = get_pysamvr_list(args)

	pysamvr_dict = dict()
	for pysamvr in pysamvr_list:
		assert pysamvr.id != None # all pysamvr must have INFO/MATEID by sanity check
		if pysamvr.id in pysamvr_dict:
			raise Exception(f'Duplicate ID detected among all input vcfs: {pysamvr.id}')
		pysamvr_dict[pysamvr.id] = pysamvr

	mateid_dict = dict()
	for pysamvr in pysamvr_list:
		if pysamvr.info['MATEID'] not in pysamvr_dict:
			raise Exception(f'MATEID of this vcf record does not exist in input vcfs:\n{pysamvr}')
		mateid_dict[pysamvr.id] = pysamvr.info['MATEID']

	mateid_tuple_set = set()
	for k,v in mateid_dict.items():
		tup = tuple(sorted(
				(k,v), 
				key = lambda x: coord_sortkey(pysamvr_dict[x].contig, pysamvr_dict[x].pos, fasta),
				))
		mateid_tuple_set.add(tup)
	mateid_tuple_list = sorted(
			mateid_tuple_set, 
			key = lambda x: coord_sortkey(pysamvr_dict[x[0]].contig, pysamvr_dict[x[0]].pos, fasta),
			)

	return pysamvr_dict, mateid_tuple_list


def main(_args = None):
	args = argument_parser(_args)
	sanity_check(args)
	fasta = get_fasta(args)

	header = get_merged_metadata_header(args, fasta)
	if args.sort == None:
		write_nosort(args, header)
	else:
		if args.sort == 'pos':
			write_possort(args, header, fasta)
		elif args.sort == 'mate':
			write_matesort(args, header, fasta)


if __name__ == '__main__':
	main()

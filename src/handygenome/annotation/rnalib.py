import re

import pysam

import importlib
top_package_name = __name__.split('.')[0]
customfile = importlib.import_module('.'.join([top_package_name, 'annotation', 'customfile']))


def get_next_codon_frame(chrom, pos, is5prime, transcript_id, is_forward, tabixfile_geneset):
	next_CDS = get_next_CDS(chrom, pos, transcript_id, tabixfile_geneset, is5prime)
	if is5prime:
		return get_leftmost_codon_frame(next_CDS, is_forward)
	else:
		return get_rightmost_codon_frame(next_CDS, is_forward)


def get_rna_range(chrom, pos, transcript_id, tabixfile_geneset):
	line = customfile.fetch_transcript_tabixline(chrom, pos - 1, pos, 
												 [transcript_id], 
												 tabixfile_geneset)[0]

	return range(line.start, line.end)


def get_next_CDS(chrom, pos, transcript_id, tabixfile_geneset, is5prime):
	"""
	Returns:
		A pysam.libctabixproxies.GTFProxy object representing 
		None if there is no nearby CDS.
	"""

	rna_range = get_rna_range(chrom, pos, transcript_id, tabixfile_geneset)

	if is5prime:
		def filter_other_side(line, pos):
			return line.start > (pos - 1)
	else:
		def filter_other_side(line, pos):
			return (line.end - 1) < (pos - 1)
											
	lines_other_side = list()
	for line in tabixfile_geneset.fetch(chrom, rna_range.start, rna_range.stop):
		if line.feature == 'CDS' and re.search(f'Parent=transcript:{transcript_id}', line.attributes):
			if filter_other_side(line, pos):
				lines_other_side.append(line)

	protein_ids = set(customfile.parse_gff3_attrs(x.attributes)['protein_id'] for x in 
					  lines_other_side)
	if len(protein_ids) > 1:
		raise Exception(f'More than one protein IDs for transcript ID "{transcript_id}"')

	if len(lines_other_side) == 0:
		return None
	else:
		if is5prime:
			return min(lines_other_side, key = lambda x: x.start)
		else:
			return max(lines_other_side, key = lambda x: x.start)


def phase_to_frame0(phase): 
	"""
	phase: As per ensembl gff3 spec. (http://asia.ensembl.org/info/website/upload/gff3.html)
		If the feature is on the reverse strand, when phase is 1, the base on the left of the rightmost base is the first base of a codon.
	frame: the first base of a codon is 0, the second 1, the third 2
	"""

	return -phase % 3


def get_leftmost_codon_frame0(CDS):
	if CDS.strand == '+':
		return phase_to_frame0(CDS.frame)
	else:
		return (phase_to_frame0(CDS.frame) + (CDS.end - CDS.start - 1))%3


def get_rightmost_codon_frame0(CDS):
	if CDS.strand == '+':
		return (phase_to_frame0(CDS.frame) + (CDS.end - CDS.start - 1))%3
	else:
		return phase_to_frame0(CDS.frame)


##################################

# for bnd annotation

def get_next_codon_frame0(transcript_id_list, chrom, pos, is5prime, tabixfile_geneset):
	def get_fetch_range(tabixfile_geneset, chrom, pos):
		range_list = list()
		for tabixline in tabixfile_geneset.fetch(chrom, pos - 1, pos):
			if customfile.check_tabixline_type(tabixline, 'gene'):
				range_list.append(range(tabixline.start, tabixline.end))

		return range(min(x.start for x in range_list),
					 max(x.stop for x in range_list))

	def classify_fetched(tabixfile_geneset, transcript_id_list, chrom, fetch_range):
		fetched_dict = dict()
		for tabixline in tabixfile_geneset.fetch(chrom, fetch_range.start, fetch_range.stop):
			parent = customfile.get_parent(tabixline)
			if parent in transcript_id_list:
				fetched_dict.setdefault(parent, list())
				fetched_dict[parent].append(tabixline)

		return fetched_dict

	def get_next_CDS(pos, transcript_id, fetched_dict, is5prime):
		if is5prime:
			relevant_features = [x for x in fetched_dict[transcript_id]
								 if x.start + 1 > pos]
			next_pos = min(x.start for x in relevant_features) + 1
			next_features = [x for x in relevant_features if x.start == next_pos - 1]
		else:
			relevant_features = [x for x in fetched_dict[transcript_id]
								 if x.end < pos]
			next_pos = max(x.end for x in relevant_features)
			next_features = [x for x in relevant_features if x.end == next_pos]

		next_feature_roles = [x.feature for x in next_features]
		if next_feature_roles.count('exon') != 1:
			raise Exception(f'Unexpected  exon is not included among the next features:\ntranscript_id: {transcript_id}, pos: {pos}, is5prime: {is5prime}\nnext_features: {next_features}')
			
		if 'exon' not in next_feature_roles:
			raise Exception(f'exon is not included among the next features:\ntranscript_id: {transcript_id}, pos: {pos}, is5prime: {is5prime}\nnext_features: {next_features}')

		next_CDS_candidates = [x for x in next_features if x.feature == 'CDS']

		if len(next_CDS_candidates) == 0:
			return None
		elif len(next_CDS_candidates) == 1:
			return next_CDS_candidates[0]
		else:
			raise Exception(f'More than one next CDSs with identical positions:\ntranscript_id="{transcript_id}", pos="{pos}", is5prime="{is5prime}"')

	def get_result(transcript_id_list, pos, fetched_dict, is5prime):
		if is5prime:
			frame0_getter = get_leftmost_codon_frame0
		else:
			frame0_getter = get_rightmost_codon_frame0

		result = dict()

		for transcript_id in transcript_id_list:
			next_CDS = get_next_CDS(pos, transcript_id, fetched_dict, is5prime)
			if next_CDS is None:
				next_codon_frame0 = None
			else:
				next_codon_frame0 = frame0_getter(next_CDS)
			result[transcript_id] = next_codon_frame0

		return result

	# main
	fetch_range = get_fetch_range(tabixfile_geneset, chrom, pos)
	fetched_dict = classify_fetched(tabixfile_geneset, transcript_id_list, chrom, fetch_range)

	return get_result(transcript_id_list, pos, fetched_dict, is5prime)


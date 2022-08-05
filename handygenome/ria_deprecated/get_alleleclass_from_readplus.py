import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
readhandler = importlib.import_module('.'.join([top_package_name, 'readplus', 'readhandler']))


# level 1 functions

def get_pos0_flanking_variant(POS0, REF, ALT, mttype = None):
	'''
	- returns 0-based positions just before/after the bases subject to variation.
	- e.g. if REF == 'ACT' and ALT == 'A', 
	  pos0_before_var is the position of 'A',
	  and pos0_after_var is the position after 'T'.
	- e.g. if REF == 'A' and ALT == 'G', 
	  pos0_before_var is the position before 'A',
	  and pos0_after_var is the position after 'A'.
	- flank length == 1
	'''

	if mttype is None:
		mttype = common.get_mttype(REF, ALT)

	if mttype == 'sv':
		pos0_before_var = None
		pos0_after_var = None
	else:
		if REF[0] == ALT[0]:
			pos0_before_var = POS0
		else:
			pos0_before_var = POS0 - 1

		pos0_after_var = POS0 + len(REF)

	return pos0_before_var, pos0_after_var


# deprecated version #

#def get_pos0_flanking_variant(POS0, REF, ALT, mttype = None):
#	'''
#	- returns 0-based positions just before/after the bases subject to variation.
#	- e.g. if REF == 'ACT' and ALT == 'A', 
#	  pos0_before_var is the position of 'A',
#	  and pos0_after_var is the position after 'T'.
#	- e.g. if REF == 'A' and ALT == 'G', 
#	  pos0_before_var is the position before 'A',
#	  and pos0_after_var is the position after 'A'.
#	- flank length == 1
#	'''
#
#	if mttype == None:
#		mttype = common_funcs.get_mttype(REF, ALT)
#
#	if mttype == 'SV':
#		pos0_before_var = None
#		pos0_after_var = None
#	else:
#		pos0_after_var = POS0 + len(REF)
#	
#		if mttype == 'SNV' or mttype == 'MNV':
#			pos0_before_var = POS0 - 1
#		elif mttype == 'INS' or mttype == 'DEL':
#			pos0_before_var = POS0
#		else: # CINDEL
#			if REF[0] == ALT[0]: 
#				'''
#				Not sure about this but CINDEL REF and ALT may have identical first letter. 
#				In this case, the first position in REF would not be a position
#				undergoing variation.
#				'''
#				pos0_before_var = POS0
#			else:
#				pos0_before_var = POS0 - 1
#
#	return pos0_before_var, pos0_after_var

######################


def generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT):
	if pos0_before_var == None or pos0_after_var == None:
		pos0_before_var, pos0_after_var = get_pos0_flanking_variant(POS0, REF, ALT)
	return pos0_before_var, pos0_after_var


def get_pos0_REFbases(POS0, REF):
	pos0_first_REFbase = POS0
	pos0_last_REFbase = POS0 + len(REF) - 1
	return pos0_first_REFbase, pos0_last_REFbase


def check_read_spans(pairs_dict, pos0):
	return pos0 in pairs_dict['refpos0']


def get_pairs_index(pairs_dict, pos0):
	if check_read_spans(pairs_dict, pos0):
		return pairs_dict['refpos0'].index(pos0)
	else:
		return None


def get_querypos0(pairs_dict, pos0):
	if check_read_spans(pairs_dict, pos0):
		pairs_index = get_pairs_index(pairs_dict, pos0)
		return pairs_dict['querypos0'][pairs_index]
	else:
		return None


def check_read_matches(pairs_dict, pos0):
	# returns True if matches
	# returns False if mismatches or does not span
	if check_read_spans(pairs_dict, pos0):
		pairs_index = get_pairs_index(pairs_dict, pos0)
		return pairs_dict['refseq'][pairs_index].isupper()
	else:
		return False


def check_chrom(read, CHROM):
	return read.reference_name == CHROM


# level 2 functions: utilizes level 1 functions

def check_read_spans_REF_positions(
		read, pairs_dict, CHROM, POS0, REF,
		pos0_first_REFbase = None, 
		pos0_last_REFbase = None,
		):
	if not check_chrom(read, CHROM):
		return False
	if pos0_first_REFbase == None or pos0_last_REFbase == None:
		pos0_first_REFbase, pos0_last_REFbase = get_pos0_REFbases(POS0, REF)

	if \
	check_read_spans(pairs_dict, pos0_first_REFbase) and \
	check_read_spans(pairs_dict, pos0_last_REFbase):
		return True
	else:
		return False


def check_read_partially_spans_REF_positions(
		read, pairs_dict, CHROM, POS0, REF,
		pos0_first_REFbase = None, 
		pos0_last_REFbase = None,
		):
	if not check_chrom(read, CHROM):
		return False
	if pos0_first_REFbase == None or pos0_last_REFbase == None:
		pos0_first_REFbase, pos0_last_REFbase = get_pos0_REFbases(POS0, REF)

	if check_read_spans_REF_positions(
		read, pairs_dict, CHROM, POS0, REF,
		pos0_first_REFbase, pos0_last_REFbase,
		):
		return False
	else:
		for pos0_candidate in range(pos0_first_REFbase, pos0_last_REFbase + 1):
			if check_read_spans(pairs_dict, pos0_candidate):
				return True
		return False


def check_read_flanks_variant(
		read, pairs_dict, CHROM, POS0, REF, ALT,
		pos0_before_var = None, 
		pos0_after_var = None,
		):
	if not check_chrom(read, CHROM):
		return False
	pos0_before_var, pos0_after_var = \
			generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT)

	if \
	check_read_spans(pairs_dict, pos0_before_var) and \
	check_read_spans(pairs_dict, pos0_after_var):
		return True
	else:
		return False


def check_read_matchflanks_variant(
		read, pairs_dict, CHROM, POS0, REF, ALT,
		pos0_before_var = None, 
		pos0_after_var = None,
		):
	if not check_chrom(read, CHROM):
		return False
	pos0_before_var, pos0_after_var = \
			generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT)

	if \
	check_read_matches(pairs_dict, pos0_before_var) and \
	check_read_matches(pairs_dict, pos0_after_var):
		return True
	else:
		return False


def check_read_doesnot_give_information(
		read, CHROM, POS0, REF, ALT,
		pos0_before_var = None, 
		pos0_after_var = None,
		):
	if not check_chrom(read, CHROM):
		return True
	pos0_before_var, pos0_after_var = \
			generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT)

	if \
	read.reference_end - 1 <= pos0_before_var or \
	read.reference_start >= pos0_after_var:
		return True
	else:
		return False


def check_read_gives_ambiguous_information(
		read, pairs_dict, CHROM, POS0, REF, ALT,
		pos0_before_var = None, 
		pos0_after_var = None,
		):
	if not check_chrom(read, CHROM):
		return False
	pos0_before_var, pos0_after_var = \
			generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT)

	if \
	( not check_read_flanks_variant(
		read, pairs_dict, CHROM, POS0, REF, ALT, 
		pos0_before_var = pos0_before_var,
		pos0_after_var = pos0_after_var,
		) ) \
	and \
	( not check_read_doesnot_give_information(
		read, CHROM, POS0, REF, ALT, 
		pos0_before_var = pos0_before_var,
		pos0_after_var = pos0_after_var,
		) ):
		return True
	else:
		return False


# level 3 functions: utilizes level 1 and level 2 functions


# The function below is deprecated.
# It is a version of get_seq_REFrange which returns a valid value
# even if the read does not flank the varying base positions.

#def get_seq_REFrange(read, CHROM, POS0, REF, ALT, pairs_dict = None):
#	'''
#	- Does not require the read to flank the positions of the varying bases,
#	  but the read must span over all REF positions to get a valid result.
#	- If chromosomes do not match, result is None.
#	- If the read does not span over all REF positions, result is None.
#	- Else, if all REF positions are cigar D, result is an empty string.
#	'''
#	if pairs_dict == None:
#		pairs_dict = readhandler.get_pairs_dict(read)
#
#	if not check_read_spans_REF_positions(read, pairs_dict, CHROM, POS0, REF):
#		seq_REFrange = None
#	else:
#		# get sliceidx_1
#		pos0_first_REFbase, pos0_last_REFbase = get_pos0_REFbases(POS0, REF)
#		sliceidx_1 = get_pairs_index(pairs_dict, pos0_first_REFbase)
#
#		# get sliceidx_2
#		if read.reference_end == pos0_last_REFbase + 1: # no downstream flanking
#			sliceidx_2 = get_pairs_index(pairs_dict, pos0_last_REFbase) + 1
#		else:
#			pos0_before_var, pos0_after_var = get_pos0_flanking_variant(POS0, REF, ALT)
#			sliceidx_2 = get_pairs_index(pairs_dict, pos0_after_var)
#
#		querypos0_list = [ 
#				x 
#				for x in pairs_dict['querypos0'][sliceidx_1:sliceidx_2]
#				if x != None
#				]
#
#		if len(querypos0_list) == 0: # all REF positions are cigar D (deletion)
#			seq_REFrange = ''
#		else:
#			querypos0_first = querypos0_list[0]
#			querypos0_last = querypos0_list[-1]
#			seq_REFrange = read.query_sequence[ querypos0_first : querypos0_last + 1 ]
#
#	return seq_REFrange


######################################################################################
# both seq_REFrange and alleleclass are None iff the read does not flank the variant #
######################################################################################

def get_querypos0_REFrange_list_from_read(
		POS0, REF, ALT, 
		pairs_dict = None,
		pos0_before_var = None, 
		pos0_after_var = None,
		):
	'''
	Returns a range object representing indices of read.query_sequence,
	which by slicing results in seq_REFrange.
	'''
	pos0_before_var, pos0_after_var = \
			generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT)

	sliceidx_1 = get_pairs_index(pairs_dict, POS0)
	sliceidx_2 = get_pairs_index(pairs_dict, pos0_after_var)

	querypos0_REFrange_list = [ 
			x 
			for x in pairs_dict['querypos0'][sliceidx_1:sliceidx_2]
			if x != None
			]

	return querypos0_REFrange_list


def get_seq_REFrange_from_read(
		read, CHROM, POS0, REF, ALT, 
		pairs_dict = None,
		pos0_before_var = None, 
		pos0_after_var = None,
		):
	'''
	- If the read does not flank the varying base positions, result is None.
	- The flanking positions need not be cigar M.
	- If chromosomes do not match, result is None.
	- If cigar operations on all REF positions are D, result is an empty string.
	'''
	if pairs_dict == None:
		pairs_dict = readhandler.get_pairs_dict(read)
		if pairs_dict == None:
			return None

	pos0_before_var, pos0_after_var = \
			generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT)

	if not check_read_flanks_variant(
			read, pairs_dict, CHROM, POS0, REF, ALT,
			pos0_before_var = pos0_before_var,
			pos0_after_var = pos0_after_var,
			):
		seq_REFrange = None
	else:
		#pos0_first_REFbase, pos0_last_REFbase = get_pos0_REFbases(POS0, REF)
		sliceidx_1 = get_pairs_index(pairs_dict, POS0)
		sliceidx_2 = get_pairs_index(pairs_dict, pos0_after_var)

		querypos0_list = [ 
				x 
				for x in pairs_dict['querypos0'][sliceidx_1:sliceidx_2]
				if x != None
				]

		if len(querypos0_list) == 0: # all REF positions are cigar D (deletion)
			seq_REFrange = ''
		else:
			seq_REFrange = read.query_sequence[ querypos0_list[0] : querypos0_list[-1] + 1 ]

	return seq_REFrange


def get_alleleclass_from_read(
		read, CHROM, POS0, REF, ALT, 
		pairs_dict = None, 
		pos0_before_var = None, 
		pos0_after_var = None,
		seq_REFrange = None,
		):
	'''
	None : read is not qualified to support any allele
	0 : supports REF
	1 : supports ALT
	2 : supports another allele different from REF or ALT
	'''
	# exclude SV
	mttype = common.get_mttype(REF, ALT)
	if mttype == 'sv':
		alleleclass = None
		return alleleclass

	# get pairs_dict if absent
	if pairs_dict == None:
		pairs_dict = readhandler.get_pairs_dict(read)
		if pairs_dict == None:
			return None

	# get pos0 before/after var if absent
	pos0_before_var, pos0_after_var = \
			generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT)

	# read which does not flank the variant is regarded as not qualified to support any allele (flank length == 1)
	if not check_read_flanks_variant(
			read, pairs_dict, CHROM, POS0, REF, ALT,
			pos0_before_var = pos0_before_var,
			pos0_after_var = pos0_after_var,
			):
		alleleclass = None
	else:
		# if the flanking bases are not cigar M, the read is regarded not to support REF nor ALT
		if not check_read_matchflanks_variant(
				read, pairs_dict, CHROM, POS0, REF, ALT,
				pos0_before_var = pos0_before_var,
				pos0_after_var = pos0_after_var,
				):
			alleleclass = 2
		else:
			# evaluate seq_REFrange if not given as argument
			if seq_REFrange == None:
				seq_REFrange = get_seq_REFrange_from_read(
						read, CHROM, POS0, REF, ALT, 
						pairs_dict = pairs_dict,
						pos0_before_var = pos0_before_var,
						pos0_after_var = pos0_after_var,
						)
			if seq_REFrange == None: # this should not happen
				alleleclass = None
				#raise Exception('seq_REFrange is None after check_read_matchflanks_variant function returned True')
			else:
				if seq_REFrange == REF:
					alleleclass = 0
				elif seq_REFrange == ALT:
					alleleclass = 1
				else:
					alleleclass = 2

	return alleleclass


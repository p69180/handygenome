import importlib
top_package_name = __name__.split('.')[0]
getalleleclass_rp = importlib.import_module('.'.join([top_package_name, 'ria', 'get_alleleclass_from_readplus']))
readhandler = importlib.import_module('.'.join([top_package_name, 'readplus', 'readhandler']))
common = importlib.import_module('.'.join([top_package_name, 'common']))


# seq_REFrange and alleleclass

def add_seq_REFrange_to_readplus(
		rp, CHROM, POS0, REF, ALT,
		pos0_before_var = None,
		pos0_after_var = None,
		):

	pos0_before_var, pos0_after_var = \
			getalleleclass_rp.generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT)

	rp.seq_REFrange = getalleleclass_rp.get_seq_REFrange_from_read(
		rp.read, CHROM, POS0, REF, ALT, 
		pairs_dict = rp.pairs_dict,
		pos0_before_var = pos0_before_var, 
		pos0_after_var = pos0_after_var,
		)


def add_alleleclass_to_readplus(
		rp, CHROM, POS0, REF, ALT, 
		pos0_before_var = None, 
		pos0_after_var = None,
		seq_REFrange = None,
		):

	pos0_before_var, pos0_after_var = \
			getalleleclass_rp.generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT)

	rp.alleleclass = getalleleclass_rp.get_alleleclass_from_read(
		rp.read, CHROM, POS0, REF, ALT, 
		pairs_dict = rp.pairs_dict, 
		pos0_before_var = pos0_before_var, 
		pos0_after_var = pos0_after_var,
		seq_REFrange = seq_REFrange,
		)


def add_seqREFrange_alleleclass_to_readplus(
		rp, CHROM, POS0, REF, ALT,
		pos0_before_var = None, 
		pos0_after_var = None,
		):

	pos0_before_var, pos0_after_var = \
			getalleleclass_rp.generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT)

	add_seq_REFrange_to_readplus(
			rp, CHROM, POS0, REF, ALT,
			pos0_before_var = pos0_before_var,
			pos0_after_var = pos0_after_var,
			)

	add_alleleclass_to_readplus(
		rp, CHROM, POS0, REF, ALT, 
		pos0_before_var = pos0_before_var,
		pos0_after_var = pos0_after_var,
		seq_REFrange = rp.seq_REFrange,
		)


# other attributes

def add_IDCinfo_to_readplus(rp):
	rp.len_ins = rp.cigarstats[1]
	rp.len_del = rp.cigarstats[2]
	rp.len_clip = rp.cigarstats[4]


def add_MM_to_readplus(rp, REF, ALT, mttype=None):
	'''
- MM: Stands for mismatch. Represents only the number of mismatches with cigar M operation.
- NM tag value includes indels as well as base mismatch ; indel length is subtracted
- Works even if the read does not have NM tag
	'''
	if mttype == None:
		mttype = common.get_mttype(REF, ALT)

	if rp.read.has_tag('NM'):
		rp.MM = rp.read.get_tag('NM') - rp.cigarstats[1] - rp.cigarstats[2]
	else:
		rp.MM = readhandler.get_MM(rp.read, rp.pairs_dict)

	# subtract the number of mismatch attributed to the variant itself
	if mttype == 'snv' or mttype == 'mnv':
		rp.MM -= len(REF)


def add_others_to_readplus(rp):
	rp.lowqual = False
	rp.lowqual_cause = list()
	rp.SVread = False
	rp.SVread_cause = list()


# main function

def ria_modify_readplus(
		rp, CHROM, POS0, REF, ALT,
		pos0_before_var = None, 
		pos0_after_var = None,
		):
	if rp.irrelevant:
		raise Exception('Irrelevant readplus object was entered.')

	add_seqREFrange_alleleclass_to_readplus(rp, CHROM, POS0, REF, ALT)
	add_IDCinfo_to_readplus(rp)
	add_MM_to_readplus(rp, REF, ALT)
	add_others_to_readplus(rp)






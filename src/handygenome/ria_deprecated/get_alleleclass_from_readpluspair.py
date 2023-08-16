import itertools
import importlib

top_package_name = __name__.split('.')[0]
module_getalleleclass_rp = importlib.import_module(top_package_name + '.ria.get_alleleclass_from_readplus')

#from .get_alleleclass_from_readplus import get_seq_REFrange_from_read
#from .get_alleleclass_from_readplus import get_alleleclass_from_read
#from .get_alleleclass_from_readplus import check_read_gives_ambiguous_information
#from .get_alleleclass_from_readplus import generate_pos0_flanking_variant_if_missing
#from .get_alleleclass_from_readplus import get_pos0_flanking_variant


def check_containing_ambiguous_read(
		rpp, CHROM, POS0, REF, ALT,
		pos0_before_var = None, 
		pos0_after_var = None,
		):
	'''
	Returns True if one or more of the reads give ambiguous information.
	Being ambiguous means partially overlapping REF positions.
	'''
	pos0_before_var, pos0_after_var = \
			module_getalleleclass_rp.generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT)

	for rp in (rpp.rp1, rpp.rp2):
		if module_getalleleclass_rp.check_read_gives_ambiguous_information(
				rp.read, rp.pairs_dict, CHROM, POS0, REF, ALT,
				pos0_before_var = pos0_before_var,
				pos0_after_var = pos0_after_var,
				):
			return True
	return False


def pick_from_two_values(val_list):
	'''
	- If both None, return None
	- If only one is None, return the non-None value
	- If both non-None and identical, return the value
	- If both non-None and not identical, return None
	'''
	assert len(val_list) == 2

	if val_list[0] == val_list[1]: # both None or both non-None
		result = val_list[0]
	else:
		not_None_list = [ x != None for x in val_list ]
		if sum(not_None_list) == 1:
			result = next(itertools.compress(val_list, not_None_list))
		else: # both not None
			result = None

	return result


def get_seq_REFrange_from_readpluspair(
		rpp, CHROM, POS0, REF, ALT,
		pos0_before_var = None, 
		pos0_after_var = None,
		):

	if rpp.SV_supporting:
		seq_REFrange = None
	else:
		pos0_before_var, pos0_after_var = \
				module_getalleleclass_rp.generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT)
	
		if check_containing_ambiguous_read(
			rpp, CHROM, POS0, REF, ALT,
			pos0_before_var = pos0_before_var,
			pos0_after_var = pos0_after_var,
		):
			seq_REFrange = None
		else:
			seq_REFrange_list = list()
			for rp in (rpp.rp1, rpp.rp2):
				if 'seq_REFrange' not in dir(rp):
					seq_REFrange_list.append(
						module_getalleleclass_rp.get_seq_REFrange_from_read(
								rp.read, CHROM, POS0, REF, ALT,
								pairs_dict = rp.pairs_dict,
								pos0_before_var = pos0_before_var,
								pos0_after_var = pos0_after_var,
						)
					)
				else:
					seq_REFrange_list.append(rp.seq_REFrange)
		
			seq_REFrange = pick_from_two_values(seq_REFrange_list)

	return seq_REFrange


def get_alleleclass_from_readpluspair(
		rpp, CHROM, POS0, REF, ALT,
		pos0_before_var = None, 
		pos0_after_var = None,
		):
	if rpp.SV_supporting:
		alleleclass = None
	else:
		pos0_before_var, pos0_after_var = \
				module_getalleleclass_rp.generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT)
	
		if check_containing_ambiguous_read(
			rpp, CHROM, POS0, REF, ALT,
			pos0_before_var = pos0_before_var,
			pos0_after_var = pos0_after_var,
		):
			alleleclass = None
		else:
			alleleclass_list = list()
			for rp in (rpp.rp1, rpp.rp2):
				if 'alleleclass' not in dir(rp):
					if 'seq_REFrange' not in dir(rp):
						alleleclass_list.append(
							module_getalleleclass_rp.get_alleleclass_from_read(
									rp.read, CHROM, POS0, REF, ALT, 
									pairs_dict = rp.pairs_dict,
									pos0_before_var = pos0_before_var,
									pos0_after_var = pos0_after_var,
									seq_REFrange = None,
							)
						)
					else:
						alleleclass_list.append(
							module_getalleleclass_rp.get_alleleclass_from_read(
									rp.read, CHROM, POS0, REF, ALT, 
									pairs_dict = rp.pairs_dict,
									pos0_before_var = pos0_before_var,
									pos0_after_var = pos0_after_var,
									seq_REFrange = rp.seq_REFrange,
							)
						)
				else:
					alleleclass_list.append(rp.alleleclass)

			alleleclass = pick_from_two_values(alleleclass_list)

	return alleleclass

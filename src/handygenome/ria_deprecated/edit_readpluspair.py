import importlib
import math
import statistics

top_package_name = __name__.split('.')[0]
module_getalleleclass_rp = importlib.import_module(top_package_name + '.ria.get_alleleclass_from_readplus')
module_getalleleclass_rpp = importlib.import_module(top_package_name + '.ria.get_alleleclass_from_readpluspair')
module_edit_readplus = importlib.import_module(top_package_name + '.ria.edit_readplus')


def add_seqREFrange_alleleclass_to_readpluspair(
		rpp, CHROM, POS0, REF, ALT,
		pos0_before_var = None, 
		pos0_after_var = None,
		):
	pos0_before_var, pos0_after_var = \
		module_getalleleclass_rp.generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT)

	rpp.seq_REFrange = module_getalleleclass_rpp.get_seq_REFrange_from_readpluspair(
		rpp, CHROM, POS0, REF, ALT,
		pos0_before_var = pos0_before_var, 
		pos0_after_var = pos0_after_var,
		)
	rpp.alleleclass = module_getalleleclass_rpp.get_alleleclass_from_readpluspair(
		rpp, CHROM, POS0, REF, ALT,
		pos0_before_var = pos0_before_var, 
		pos0_after_var = pos0_after_var,
		)


def add_IDCinfo_to_readpluspair(rpp):
	rpp.len_ins = rpp.rp1.len_ins + rpp.rp2.len_ins
	rpp.len_del = rpp.rp1.len_del + rpp.rp2.len_del
	rpp.len_clip = rpp.rp1.len_clip + rpp.rp2.len_clip


def add_MQ(rpp):
	rpp.MQ = statistics.mean((rpp.rp1.read.mapping_quality, rpp.rp2.read.mapping_quality))


def add_MM(rpp): # mismatch
	rpp.MM = rpp.rp1.MM + rpp.rp2.MM


###

def ria_modify_readpluspair(
		rpp, CHROM, POS0, REF, ALT,
		pos0_before_var = None, 
		pos0_after_var = None,
		):
	pos0_before_var, pos0_after_var = \
		module_getalleleclass_rp.generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT)

	# modify rpp.rp1 and rpp.rp2
	module_edit_readplus.ria_modify_readplus(
			rpp.rp1, CHROM, POS0, REF, ALT, 
			pos0_before_var = pos0_before_var, 
			pos0_after_var = pos0_after_var,
			)
	module_edit_readplus.ria_modify_readplus(
			rpp.rp2, CHROM, POS0, REF, ALT, 
			pos0_before_var = pos0_before_var, 
			pos0_after_var = pos0_after_var,
			)

	# modify rpp itself
	add_seqREFrange_alleleclass_to_readpluspair(
			rpp, CHROM, POS0, REF, ALT,
			pos0_before_var = pos0_before_var, 
			pos0_after_var = pos0_after_var,
			)
	add_IDCinfo_to_readpluspair(rpp)
	add_MQ(rpp)
	add_MM(rpp)


def ria_modify_readpluspairlist(
		rpplist, CHROM, POS0, REF, ALT,
		pos0_before_var = None, 
		pos0_after_var = None,
		):
	pos0_before_var, pos0_after_var = \
		module_getalleleclass_rp.generate_pos0_flanking_variant_if_missing(pos0_before_var, pos0_after_var, POS0, REF, ALT)

	for rpp in rpplist.rpplist:
		if not rpp.irrelevant:
			ria_modify_readpluspair(
					rpp, CHROM, POS0, REF, ALT,
				pos0_before_var = pos0_before_var, 
				pos0_after_var = pos0_after_var,
					)
	

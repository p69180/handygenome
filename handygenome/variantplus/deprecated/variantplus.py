import importlib
import cyvcf2

top_package_name = __name__.split('.')[0]
module_common_params = importlib.import_module(top_package_name + '.common.params')
module_repeats = importlib.import_module(top_package_name + '.variantplus.repeats')
module_caller_format_check = importlib.import_module(top_package_name + '.svtools.caller_format_check')

variant_format_checker_dict, caller_list = moduel_caller_format_check.get_variant_format_checker_dict()


class VariantPlus:
	def __init__(
			self, 
			variant=None, 
			CHROM1=None, 
			POS1=None, 
			REF=None, 
			ALT=None,
			CHROM2=None, 
			POS2=None, 
			CT=None,
			fasta=None,
			):
		'''
		Argument requirements:
		1) 'variant' is supplied
		2) 'CHROM1', 'POS1', 'REF', 'ALT' are supplied
		'''
		set_basic_attributes(self, variant, CHROM1, POS1, REF, ALT, CHROM2, POS2, fasta)
		modify_REF(self)
		self.mttype = get_mttype(self.REF, self.ALT)
		set_repeat_list(self)


	def INFO_has_field(self, field):
		try:
			value = self.variant.INFO[field]
		except KeyError:
			return False
		else:
			if value == '.':
				return False
			else:
				return True


### set_basic_attributes

def set_basic_attributes(vp, variant, CHROM1, POS1, REF, ALT, CHROM2, POS2, fasta):
	if variant == None:
		set_basic_attributes_without_variant(vp, CHROM1, POS1, REF, ALT, CHROM2, POS2, fasta)
	else:
		set_basic_attributes_with_variant(vp, variant, fasta)


def set_basic_attributes_with_variant(vp, variant, fasta): # vp: variantplus
	vp.variant = variant
	vp.fasta = fasta

	vp.CHROM1 = variant.CHROM
	vp.POS1 = variant.POS
	vp.POS01 = vp.POS1 - 1
	vp.REF = variant.REF
	vp.ALT = variant.ALT[0]

	vp.CHROM2 = None
	vp.POS2 = None
	vp.POS02 = None


def set_basic_attributes_without_variant(vp, CHROM1, POS1, REF, ALT, CHROM2, POS2, fasta): # vp: variantplus
	vp.variant = None
	vp.fasta = fasta

	vp.CHROM1 = CHROM1
	vp.POS1 = POS1
	vp.POS01 = vp.POS1 - 1
	vp.REF = variant.REF
	vp.ALT = variant.ALT[0]

	vp.CHROM2 = CHROM2
	vp.POS2 = vp.POS2
	vp.POS02 = vp.POS2 - 1

### modify_REF

def modify_REF(vp):
	assert 'REF' in dir(vp)
	if vp.REF == 'N':
		assert vp.fasta != None
		vp.REF = vp.fasta.fetch(vp.CHROM1, vp.POS01, vp.POS01 + 1)

### set_repeat_list

def set_repeat_list(vp):
	vp.repeat_list = module_repeats.get_repeat_list(vp.CHROM1, vp.POS01, vp.fasta, pre=10, post=10)
	vp.relevant_repeat_list = module_repeats.get_relevant_repeat_list(
			vp.repeat_list, vp.CHROM1, vp.POS01, vp.REF, vp.ALT, mttype = vp.mttype, flanklen = 5,
			)

# get_mttype

def get_mttype(REF, ALT):
	if module_common_params.PAT_NUCLEOBASES.match(ALT): # ALT is [ACGTNacgtn]+
		if len(REF) == len(ALT):
			if len(REF) == 1:
				mttype = 'SNV'
			else:
				mttype = 'MNV'
		else:
			if len(REF) == 1:
				mttype = 'INS'
			elif len(ALT) == 1:
				mttype = 'DEL'
			else:
				mttype = 'CINDEL'
	else:
		mttype = 'SV'

	return mttype

# get_sv_format

def get_sv_format(variant):
	'''
	checks caller format of a variant line of a SV vcf
	'''

	true_formats = list()
	for caller, checker in variant_format_checker_dict.items():
		if checker(variant): # checker(variant) returns a boolean
			true_formats.append(caller)

	if len(true_formats) == 0:
		raise Exception(f'Input variant line format does not conform with any of the callers: {caller_list}')
	elif len(true_formats) > 1:
		raise Exception(f'Input variant line format conforms with more than one of these callers: {caller_list}')

	sv_format = true_formats[0]
	
	return sv_format


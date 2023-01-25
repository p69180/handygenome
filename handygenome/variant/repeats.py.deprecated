import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


class Repeat:
	def __init__(self, CHROM=None, start_localcoord=None, end_localcoord=None, start_refcoord=None, end_refcoord=None, repeat_unit=None, repeat_count=None):
		self.CHROM = CHROM
		self.start_localcoord = start_localcoord
		self.end_localcoord = end_localcoord
		self.start_refcoord = start_refcoord
		self.end_refcoord = end_refcoord
		self.repeat_unit = repeat_unit
		self.repeat_count = repeat_count

# get_repeat_list

# not used #
def check_already_within_repeat(subseq_start, repeat_range_localcoord_list):
	for repeat_range_localcoord in repeat_range_localcoord_list:
		if subseq_start in repeat_range_localcoord:
			return True
	return False
############


def check_already_found_repeat(repeat, repeat_list):
	for known_repeat in repeat_list:
		if repeat.start_localcoord not in range(known_repeat.start_localcoord, known_repeat.end_localcoord):
			continue
		if len(repeat.repeat_unit) % len(known_repeat.repeat_unit) != 0:
			continue
		if known_repeat.repeat_unit * ( len(repeat.repeat_unit) // len(known_repeat.repeat_unit) ) != repeat.repeat_unit:
			continue
		if repeat.start_localcoord not in range(known_repeat.start_localcoord, known_repeat.end_localcoord, len(known_repeat.repeat_unit)):
			continue
		return True

	return False


def get_repeat_count(
		remaining_seq,
		repeat_unit,
		):
	repeat_count_max = len(remaining_seq)//len(repeat_unit)

	repeat_count = 1
	while True:
		repeat_count += 1
		if repeat_unit * repeat_count == remaining_seq[ : len(repeat_unit)*repeat_count ]:
			if repeat_count == repeat_count_max:
				break
			else:
				continue
		else:
			repeat_count -= 1
			break

	return repeat_count


def get_repeat_list_from_seq(seq):
	repeat_list = list() # elements: tuple (repeat_range_refCoord, repeat_unit, repeat_count)

	for subseq_start in range(len(seq)):
		remaining_seq = seq[subseq_start:]
		for subseq_end in range(subseq_start + 1, len(seq) + 1):
			repeat = Repeat(
					start_localcoord = subseq_start,
					repeat_unit = seq[ subseq_start : subseq_end ],
					)
			if check_already_found_repeat(repeat, repeat_list):
				continue

			repeat_count = get_repeat_count(remaining_seq, repeat.repeat_unit)
			if repeat_count >= 2: # this is a true repeat
				repeat.repeat_count = repeat_count
				repeat.end_localcoord = repeat.start_localcoord + len(repeat.repeat_unit) * repeat_count
				repeat_list.append(repeat)

	return repeat_list


def get_repeat_list(CHROM, POS0, fasta, pre=10, post=40):
	seq = fasta.fetch(CHROM, POS0 - pre, POS0 + post + 1)
	pad = POS0 - pre

	repeat_list = get_repeat_list_from_seq(seq)
	for repeat in repeat_list:
		repeat.CHROM = CHROM
		repeat.start_refcoord = repeat.start_localcoord + pad
		repeat.end_refcoord = repeat.end_localcoord + pad
			
	return repeat_list


# get_relevant_repeat_list

def indelseq_in_repeat(indelseq, repeat):
	return indelseq in repeat.repeat_unit * repeat.repeat_count


def check_repeat_is_near_REFpos(repeat, CHROM, REF_start, REF_end, flanklen):
	if repeat.CHROM != CHROM:
		return False
	else:
		padded_REF_start = REF_start - flanklen
		padded_REF_end = REF_end + flanklen
		
		if repeat.start_refcoord > padded_REF_start and repeat.end_refcoord < padded_REF_end:
			return True # repeat near REF
		else:
			return False # repeat not near REF


def get_relevant_repeat_list(
	repeat_list,
	CHROM,
	POS0,
	REF,
	ALT,
	mttype = None,
	flanklen = 5,
):	
	if mttype is None:
		mttype = common.get_mttype(REF, ALT)

	if mttype == 'ins' or mttype == 'del':
		indelseq = common.get_indelseq(REF, ALT)
		REF_start = POS0
		REF_end = POS0 + len(REF)
			
		relevant_repeat_list = list()
		for repeat in repeat_list:
			if \
			indelseq_in_repeat(indelseq, repeat) and \
			check_repeat_is_near_REFpos(repeat, CHROM, REF_start, REF_end, flanklen):
				relevant_repeat_list.append(repeat)
	else:
		relevant_repeat_list = list()
		
	return relevant_repeat_list



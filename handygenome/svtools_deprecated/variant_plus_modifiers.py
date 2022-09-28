import re

svALT_pat_list = [
	re.compile('^(.+)\[(.+):(.+)\[$'),
	re.compile('^(.+)\](.+):(.+)\]$'),
	re.compile('^\[(.+):(.+)\[(.+)$'),
	re.compile('^\](.+):(.+)\](.+)$'),
	]


def get_SVTYPE(vp):
	if vp.INFO_has_field('SVTYPE'):
		return vp.variant.INFO['SVTYPE']


def get_svALT(vp):


def main(vp):
	vp.SVTYPE = get_SVTYPE(vp)

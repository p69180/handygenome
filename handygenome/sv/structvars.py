import pysam
import Bio.Seq

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
initvcf = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'initvcf']))
breakends = importlib.import_module('.'.join([top_package_name, 'sv', 'breakends']))
libvcfspec = importlib.import_module('.'.join([top_package_name, 'variant', 'vcfspec']))


MITOCHONDRIAL = ('chrM', 'MT')


class SimpleStructuralVariant(common.Interval):

	def __init__(self, chrom, start1 = None, end1 = None, start0 = None, end0 = None, fasta = None, refver = None):
		"""
		'chrom' is mandatory
		('start1' and 'end1') or ('start0' and 'end0') must be given(for coordinate setting).
		'fasta' or 'refver' may be given

		'start0' and 'end0' are 0-based half-open system.
		'start1' and 'end1' are 1-based closed system.
		The interval given by [start0, end0) or [start1, end1] itself represents the affected bases. 
			(e.g. The bases on 'start1' or 'end1' are affected)
		The sequence represented by the interval [start0, end0) is affected.
		"""

		super().__init__(chrom, start1, end1, start0, end0)

		if fasta is None:
			self.fasta = pysam.FastaFile(common.DEFAULT_FASTA_PATHS[refver])
		else:
			self.fasta = fasta

	def __repr__(self):
		return f'<{self.__class__.__name__}> {self.chrom}:{self.start1:,}-{self.end1:,}'

	def get_seq(self):
		return self.fasta.fetch(self.chrom, self.start0, self.end0)

	def get_vr_symbolic(self):
		chromdict = common.ChromDict(fasta = self.fasta)
		vr = initvcf.create_vr(chromdict = chromdict)
		vr.contig = self.chrom
		vr.pos = self.start0
		vr.ref = self.fasta.fetch(self.chrom, vr.pos - 1, vr.pos)
		vr.alts = [ self.__class__.ALTsymbol ]
		vr.stop = self.end0

		return vr


class Deletion(SimpleStructuralVariant):

	ALTsymbol = '<DEL>'

	def get_hgvsg(self):
		prefix = 'm' if (self.chrom in MITOCHONDRIAL) else 'g'
		if self.length == 1:
			hgvsg = f'{self.chrom}:{prefix}.{self.start0 + 1}del'
		else:
			hgvsg = f'{self.chrom}:{prefix}.{self.start0 + 1}_{self.end0}del'

		return hgvsg

	def get_vcfspec(self):
		pos = self.start0
		ref = self.fasta.fetch(self.chrom, pos - 1, self.end0)
		alt = ref[0]

		return libvcfspec.Vcfspec(self.chrom, pos, ref, [alt])
			
	def get_bnds(self):
		return breakends.Breakends(
				chrom1 = self.chrom,
				pos1 = self.start0,
				pos1_endis5 = False,
				chrom2 = self.chrom,
				pos2 = self.end0 + 1,
				pos2_endis5 = True,
				inserted_seq = list(),
				fasta = self.fasta,
				)


class TandemDuplication(SimpleStructuralVariant):

	ALTsymbol = '<DUP:TANDEM>'

	def get_hgvsg(self):
		prefix = 'm' if (self.chrom in MITOCHONDRIAL) else 'g'
		if self.length == 1:
			hgvsg = f'{self.chrom}:{prefix}.{self.start0 + 1}dup'
		else:
			hgvsg = f'{self.chrom}:{prefix}.{self.start0 + 1}_{self.end0}dup'

		return hgvsg

	def get_vcfspec(self):
		pos = self.start0
		fetched = self.fasta.fetch(self.chrom, pos - 1, self.end0)
		ref = fetched[0]
		alt = fetched

		return libvcfspec.Vcfspec(self.chrom, pos, ref, [alt])


class Inversion(SimpleStructuralVariant):

	ALTsymbol = '<INV>'

	def get_hgvsg(self):
		prefix = 'm' if (self.chrom in MITOCHONDRIAL) else 'g'
		hgvsg = f'{self.chrom}:{prefix}.{self.start0 + 1}_{self.end0}inv'

		return hgvsg

	def get_vcfspec(self, padded = False):
		if padded:
			pos = self.start0
			fetched = self.fasta.fetch(self.chrom, pos - 1, self.end0)
			ref = fetched
			alt = fetched[0] + Bio.Seq.reverse_complement(fetched[1:])
		else:
			pos = self.start0 + 1
			fetched = self.fasta.fetch(self.chrom, pos - 1, self.end0)
			ref = fetched
			alt = Bio.Seq.reverse_complement(fetched)

		return libvcfspec.Vcfspec(self.chrom, pos, ref, [alt])


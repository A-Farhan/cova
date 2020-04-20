import os,pkg_resources
from unittest import TestCase
from Bio.Data.CodonTable import unambiguous_dna_by_id as codon_table
from Bio import AlignIO
from cova import _utils

DATAPATH = pkg_resources.resource_filename('cova', 'testdata')

class timer(TestCase):
	def test_empty(self):
		self.assertEqual(_utils.timer(),'00:00:0')

class readcsv(TestCase):
	def test_empty(self):
		self.assertRaises(TypeError, _utils.readcsv)
	def test_negative(self):
		self.assertRaises(FileNotFoundError, _utils.readcsv, 'bla')
	def test_positive(self):
		self.assertEqual(first=_utils.readcsv(fl=os.path.join(DATAPATH,'testfile.csv'),header=True), 
			second=(['A','B','C'], [['1','2','3'],['a','b']]))

class splitdata(TestCase):
	data = _utils.readcsv(fl=os.path.join( DATAPATH, 'testfile.csv'), header=True)[1]
	def test_negative(self):
		self.assertRaises(IndexError, _utils.split_data, data=self.data, ix=0, cixs=2)
	def test_positive_int(self):
		self.assertEqual(first=_utils.split_data(self.data,0,1), second={'1': ['2'], 'a': ['b']})
	def test_positive_list(self):
		self.assertEqual(first=_utils.split_data(self.data,0,[0,1]), \
			second={'1': [['1', '2']], 'a': [['a', 'b']]})

class listep(TestCase):
	data = _utils.readcsv(fl=os.path.join( DATAPATH, 'test_listep.csv'))
	data = [ [ int(j) for j in i] for i in data]
	def test_mix(self):
		self.assertEqual(first=_utils.list_ep(self.data[0]), second=[(1,3),(5,5),(7,9)])
	def test_stretch(self):
		self.assertEqual(first=_utils.list_ep(self.data[1]), second=[(10,14)])
	def test_ends(self):
		self.assertEqual(first=_utils.list_ep(self.data[2]), second=[(15,15),(20,20)])
	def test_disorder(self):
		self.assertIsNone(_utils.list_ep(self.data[3]))
	def test_step(self):
		self.assertEqual(first=_utils.list_ep(als=self.data[0],step=2), \
			second=[(1,1),(2,2),(3,7),(8,8),(9,9)])

class n2c(TestCase):
	def test_char_neg(self):
		nseq = 'atgcgtccraaw'
		self.assertRaises(_utils.QualSeqError, _utils.n2c, nseq=nseq)
	def test_cds_neg(self):
		nseq = 'acgt'
		self.assertRaises(_utils.LenSeqError, _utils.n2c, nseq=nseq)
	def test_pos(self):
		nseq = 'atgcgtcca'
		self.assertEqual(first=_utils.n2c(nseq), second=['atg','cgt','cca'])

class nv2av(TestCase):
	cseq = ['ATG','TGC','AAA','TAG']
	
	def test_nonnegp(self):
		try:
			out = _utils.nv2av(p=-3, v='C', seq=self.cseq)
		except ValueError as ex:
			msg = str(ex)
		self.assertEqual(first=msg,second='Position must be a non-negative integer') 
	
	def test_varlcase(self):
		self.assertRaises(_utils.LowercaseSeqError, _utils.nv2av, p=3, v='c', seq=self.cseq)
	def test_vbase(self):
		try:
			out = _utils.nv2av(p=2, v='R', seq=['ATG'])
		except ValueError as ex:
			msg = str(ex)
		self.assertEqual(first=msg,second='v must be an unambiguous DNA base') 
	
	def test_cdslcase(self):
		cseq = [ i.lower() for i in self.cseq]
		self.assertRaises(_utils.LowercaseSeqError, _utils.nv2av, p=3, v='C', seq=cseq)
	def test_cdsbase(self):
		try:
			out = _utils.nv2av(p=2, v='C', seq=['NTG'])
		except ValueError as ex:
		self.assertEqual(first=str(ex),second='seq must be a list of codons') 
	def test_codons(self):
		cseq = [ i for i in self.cseq]
		cseq[0] = 'AT'
		self.assertRaises(_utils.LenSeqError, _utils.nv2av, p=3, v='C', seq=cseq)
	
	def test_vlen(self):
		self.assertRaises(_utils.LenSeqError, _utils.nv2av, p=12, v='C', seq=self.cseq)

class extractcds(TestCase):
	msa = AlignIO.read(handle=os.path.join(DATAPATH,'aln.fna'), format='fasta')
	def test_ends_int(self):
		try:
			out = _utils.extract_cds(begin='a', end=3, msa=self.msa, has_stop=False, outf='cds_aln.fna')
		except TypeError as ex:
			self.assertEqual(first=str(ex), second='end positions must be integers.')
	def test_msa_type(self):
		msa = _utils.MSA(os.path.join(DATAPATH,'aln.fna'))
		try:
			out = _utils.extract_cds(begin=1, end=3, msa=msa, has_stop=False, outf='cds_aln.fna')
		except TypeError as ex:
			self.assertEqual(first=str(ex), second='msa must be a biopython MultipleSeqAlignment.')
	def test_end_after_begin(self):
		try:
			out = _utils.extract_cds(begin=5, end=3, msa=self.msa, has_stop=False, outf='cds_aln.fna')
		except ValueError as ex:
			self.assertEqual(first=str(ex), second='end must be greater than begin.')

class testmsa(TestCase):
	msa = _utils.MSA(os.path.join(DATAPATH,'aln.fna'))
	def test_nofile(self):
		try:
			_utils.MSA('nonexistingfile')
		except IOError as ex:
			self.assertIn(first="failed to open the file", second=str(ex))
	def test_ref_exist(self):
		try:
			msa.limref(ref='4')
		except ValueError as ex:
			self.assertEqual(first=str(ex), second='reference must be present in the MSA.')
	def test_ins_qual(self):
		try:
			msa.ins(ref='2',ambt=-1)
		except ValueError as ex:
			self.assertEqual(first=str(ex), second='threshold must be non-negative.')
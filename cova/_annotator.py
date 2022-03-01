from cova import FEATURETABLE, GENOME, RFS, CDS, PSEQS
from cova import utils
from Bio.Data.CodonTable import unambiguous_dna_by_id as codon_table
import os, sys, pandas, math, multiprocessing, numpy
from time import time

#### Point mutations #######
def ann_pm(vpos,vseq,ft=FEATURETABLE,cdss=CDS,ct=codon_table[1],rfs=RFS):
	"""
	Annotate a point mutation in a given protein with the corresponding amino acid change. 

	Arguments:
	* Mandatory
	- vpos 	- variant's genomic position
	- vseq 	- variant's sequence
	* Optional
	- ft 	- reference feature table with protein ids for index
	- cdss 	- reference CDS biopython sequence records
	- ct 	- biopython codon table
	- rfs 	- dataframe with info on ribosomal slippage cases

	Value:
	A pandas dataframe with 1 row for every variant codon and following columns-
	- position 		- 1-indexed genomic position
	- var_base 		- variant base sequence
	- protein_id 	- product accession of the target protein
	- name 			- common name for the protein
	- ref_codon 	- reference codon at the variant position
	- var_codon 	- variant codon
	- aa_change 	- amino acid substitution caused by the variant
	"""
	# initialize output dataframe
	colnames = ['position','var_base','protein_id','name','ref_codon','var_codon','aa_change']
	out = pandas.DataFrame( columns=colnames)
	# find the affected protein(s)
	prots = ft[ (ft['start'] <= vpos) & (vpos <= ft['end'])].index
	
	# return empty if no protein is affected
	if len(prots) == 0:
		return out
	
	# 0-index of the variant's genome position
	genome_x = vpos-1

	# for every affected protein
	for p in prots:
		# 0-indexed ends of the CDS
		b = ft.loc[p,'start']-1
		e = ft.loc[p,'end']
		# corresponding CDS sequence 
		cds_seq = cdss[p]
		
		# variant's index in the CDS
		if p in rfs.index:
			rfs_pos = rfs.loc[p,'genomic_position']
			rfs_type = rfs.loc[p,'type']

			if genome_x < (rfs_pos-1):
				cds_x = genome_x - b
			else:
				cds_x = genome_x - b - rfs_type

		else:
			cds_x = genome_x - b
		
		# corresponding list of codons
		try:
			codonls = utils.n2c(cds_seq)
		except utils.LenSeqError:
			print("\tinvalid CDS!")
			continue
			
		# list of amino acid variant(s)
		try:
			avs = utils.nv2av(p=cds_x, v=vseq, seq=codonls)
		except (utils.LenSeqError,ValueError):
			print('''Invalid variant {}'''.format(vpos))
			continue

		avs = pandas.DataFrame( [ [vpos,vseq,p,ft.loc[p,'name']]+i for i in avs], columns=colnames)
		out = out.append(avs)

	return out

def ann_pm_apply(vrow,ft=FEATURETABLE,chars='ACGT'):
	# alleles at the position in the samples
	alls = vrow[2:]
	# alternative allelles
	alts = set(alls[ alls != vrow.ref])
	# only retain permissible characters
	alts = { i for i in alts if all(j in chars for j in i)}

	if len(alts) == 0:
		print("No valid alternative alleles!")
		return
	
	alts = pandas.Series(list(alts))
	vrs = alts.apply( lambda x: ann_pm(vrow.pos,x,ft))
	vlist = vrs.tolist()
	
	# if no valid variants
	if len(vlist) == 0:
		return
	
	out = pandas.concat(vlist,ignore_index=True)
	out['ref_base'] = vrow.ref
	# list of genomes
	out['genomes'] = out.var_base.apply( lambda x: alls[alls==x].index)
	# frequency
	out['freq'] = out.genomes.apply(len)
	out['genomes'] = out.genomes.apply( lambda x: ','.join(x))
	return out

def annotate_pm(fin,fout):
	"""
	Annotate point mutations located within protein regions. Identify amino acid changes
	corresponding to nucleotide changes. The output table is directly saved to a file. 

	Arguments:
	- fin 	- full path to input file of point mutations table
	- fout 	- full path to output file of annotated mutations
	
	Value:
	The output table has the following 7 columns.
	- 1 : protein id
	- 2 : protein name, as present in the reference feature table
	- 3 : 1-indexed genomic position
	- 4 : nucleotide variant
	- 5	: reference codon
	- 6 : variant codon
	- 7 : amino acid change
			A string with 3 components:
			- reference aa
			- 1-indexed aa position 
			- variant aa
	"""
	# variant table
	vartab = pandas.read_csv( fin, sep='\t', keep_default_na=False) 
	# perform annotations for every row of the variant table	
	anvs = vartab.apply( axis=1, func=ann_pm_apply)
	anvs = [ i for i in anvs if i is not None]

	# if no valid variants
	if len(anvs) == 0:
		return

	# reconvert to pandas dataframe
	out = pandas.concat( anvs, ignore_index=True)
	# get type of variant ( Syn / Non-syn)
	out['type'] = out.aa_change.apply( lambda x: (x[0] != x[-1]) and 'N' or 'S') 
	# output
	out = out[['protein_id','name','position','ref_base','var_base','ref_codon',
	'var_codon','aa_change','type','freq','genomes']]
	out.to_csv( fout, sep='\t', index=False)

######### Deletions ##############
def ann_del_single(pb,l,fif):
	"""
	Annotate a deletion, specified by its starting position on the genome
	and its length, with a genomic feature.
	"""
	# corresponding amino acid position
	## needs to be adjusted for ribosomal slippage
	p=math.ceil( (pb-fif['start'])/3+1)
	# check if deletion's length is a multiple of 3
	# if not, then its a frameshift
	if l%3 != 0:
		out=fif['name']+'_del'+str(p)+'_'+'frameshift'
	# otherwise, specify length of the corresponding amino acid sequence deletion
	else:
		l=int(l/3)
		aaseq = str(PSEQS[fif.name].seq[p-1:p-1+l])
		out=fif['name']+'_del'+str(p)+'-'+aaseq
	
	return out

def ann_del_apply(bpos,dlen,ft=FEATURETABLE):
	# 1-indexed end position
	epos = bpos+dlen-1
	# extract all features that overlap with these coordinates
	t = ft[ (bpos >= ft.start) & (epos <= ft.end) ]
	# corresponding annotation(s) from these features
	anns = t.apply( axis=1, func= lambda x: ann_del_single(pb=bpos,l=dlen,fif=x))
	# comma-separated string of annotations
	out = ','.join( set(anns.values))
	return out

def annotate_del(fin,fout):
	# load data
	df = pandas.read_csv(fin,sep='\t')
	# get frequencies
	freqs = df.iloc[:,3:].sum(axis=1)
	# extract pos and len column to perform annotation
	df = df[['pos','len']]
	out = df.apply(axis=1,func=lambda x: ann_del_apply(bpos=x['pos'],dlen=x['len']))
	odf = pandas.DataFrame(zip(df['pos'],df['len'],out,freqs), columns=['pos','len','ann','freq'])    
	odf.to_csv(fout,sep='\t',index=False)
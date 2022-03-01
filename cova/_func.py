import os, re, subprocess, matplotlib, seaborn, pandas
from . import utils, FEATURETABLE, GENOME, CODONTABLE, TYPEPOS, SEQTYPES
from time import time
from Bio import SeqIO, AlignIO

def rm_genome_w_stopm(vtab):
	"""Return a dataframe of genomes with nonsense variants given a 
	dataframe of genomes with their shared & unique variants."""
	# dict of genome with their list of variants
	vtab['unique'] = vtab['unique'].fillna('')    
	genome_vrs = vtab[['shared','unique']].apply(axis=1, 
		func=lambda x: ','.join(x.astype(str))).apply( lambda x: x.split(',')) 
	# identify genomes w/ stop mutations
	genome_stopm = genome_vrs.apply(utils.get_stopm)
	genome_stopm = genome_stopm[genome_stopm!='']
	return genome_stopm

def extract_nmsa_prots(fmsa,dr):
	"""
	Extract nucleotide MSA of proteins from Reference-limited Whole-genome Multiple Sequence Alignment.
	
	Arguments:
	- fmsa 	- path to the reference-limited Whole-genome Multiple Sequence Alignment
	- dr 	- path to the output directory of extracted MSAs
	"""
	# biopython multiple sequence alignment
	msa = AlignIO.read(fmsa, 'fasta')
	# proteins' end points
	prot_ends = FEATURETABLE[ ['feature', 'start','end','name']]
	# apply function to extract CDS MSA
	out = prot_ends.apply( axis=1, func=lambda x: 
		utils.extract_cds_msa(p=x.name,pstart=x.start,pstop=x.end,msa=msa)) 
	
	# for every item
	for k,v in out.items():
		# prepare a file path for the output
		fout = os.path.join( dr, prot_ends.loc[k,'name']+'.msa')
		# write the MSA to this path
		AlignIO.write(alignments=[v], handle=fout, format='fasta')

def run_fubar(fmsa,ftree,outdr,prog):
	"""
	Run positive selection analysis using FUBAR from Hyphy.

	Arguments:
	- fmsa 	- full path to input MSA file
	- ftree - full path to input tree file based on the above MSA
	- outdr - full path to the directory to store FUBAR output
	- prog 	- hyphy program, full path may be required
	"""
	## checks
	# msa file is present
	if not os.path.exists(fmsa):
		raise FileNotFoundError('msa file %s must be present.'%fmsa)
	# tree file is present
	if not os.path.exists(ftree):
		raise FileNotFoundError('tree file %s must be present.'%fmsa)
	# program is available
	cmd = [prog,'-h']
	try:
		runstat = subprocess.run(cmd,stdout=subprocess.DEVNULL)
	except FileNotFoundError:
		raise FileNotFoundError("couldn't find the HYPHY program. You may want to check the path.")
	
	# if output file is already present
	indr = os.path.dirname(fmsa)
	f = os.path.basename(fmsa)
	p = f.replace('.msa','') 
	fsout = os.path.join( outdr, p+'.out')
	
	# json output intermediate file
	fjout_int = os.path.join( indr, f +'.FUBAR.json')
	# json output file
	fjout = os.path.join( outdr, p+'.json')
	# cache file
	fcout = os.path.join( outdr, p+'.cache')

	cmd = [ prog, 'fubar', '--alignment', fmsa, '--tree', ftree, '--cache', fcout]
	print("Command: %s"%' '.join(cmd))
	with open(fsout,'w') as flob:
		s = subprocess.run(cmd, stdout=flob)	
	if s.returncode == 0:
		os.rename(src=fjout_int, dst=fjout)

def parse_fubar(indr,frout,fsout):
	"""
	Parse FUBAR output to generate 
	- table of dN/dS rates of proteins
	- table of protein sites under positive selection.
	
	Arguments:
	- indr - full path to directory with the results of FUBAR analysis
	- frout - full path to output rates file
	- fsout - full path to output sites file 
	"""
	## checks
	if not os.path.exists(indr):
		raise IOError("Couldn't find the input directory %s."%indr)
			
	# initialize lists for rates ans sites output tables
	rates_out = []
	sites_out = []
	# list of FUBAR output files, to be used as input here
	infls = [ i for i in os.listdir(indr) if i.endswith('.out')]
	
	# for each of these files
	for f in infls:
		# protein name
		p = f.replace('.out', '')
		# full path to the file
		fin = os.path.join( indr, f)
		
		# extract file's contents
		with open(fin) as flob:
			contents = [ i.strip('\n') for i in flob]
		
		# tree length
		tln = [ i for i in contents if 'Tree length' in i]
		
		# skip the protein, if no change, as estimated by tree length
		if len(tln) == 0:
			print('\t"Tree Length" entry missing! Check the FUBAR output file for details.%s.'%p)
			continue
		else:
			tln = tln[0]
		
		# extract tree length from the string
		t = float( re.split(':\W+', tln)[1])	
		# lines with syn and non-syn rates
		rate_lns = [ i for i in contents if 'synonymous' in i]
		# extract the rates from the lines above and make an entry in the output list
		rates = [p,t] + [ float(i.split(' = ')[1]) for i in rate_lns]
		rates_out.append(rates)
		
		# extract the table of positively selected sites
		sites_tabs = [ re.split('\W*\|\W*',i) for i in contents if '|' in i][2:]
		# process the table and make an entry in the output list
		sites_ls = [ [ p, int(i[1]), float(i[3]), float(i[4]), float(i[5].split(' = ')[1]) ]\
						for i in sites_tabs]
		sites_out.extend(sites_ls)
	
	# further process rates output table	
	rates_out = [ i + [ round(i[-1]-i[-2],3)] for i in rates_out]
	rates_out = sorted( rates_out, key=lambda x: x[-1], reverse=True)
	# convert both to pandas dataframe
	rates_out = pandas.DataFrame(rates_out,columns=['protein', 'exp_subs','syn', 'nonsyn', 'dnds'])
	sites_out = pandas.DataFrame(sites_out,columns=['protein','site','syn', 'nonsyn', 'post_prob'])
	# write both tables to files
	rates_out.to_csv(frout,index=False)
	sites_out.to_csv(fsout,index=False)

def genome_seqtype(msa,fst=None,pos=TYPEPOS,sts=SEQTYPES):
    """
	Return Sequence type of genomes from the input alignment.
	By default, CoVa's sequence types are used. Optionally, a file can be provided for sequence types.
    """    
    # if a sequence type file was provided
    if fst is not None:
    	print("A sequence type file was provided. Not using CoVa's standard sequence types")
    	sts = pandas.read_csv(fst,index_col=0)
    	pos = [ int(i) for i in sts.columns]
    	sts	= sts.apply( axis = 1, func = lambda x: ''.join(x))
    	sts = dict(zip(sts.values,sts.index))  

    # list of column sequences in the aligment at these positions
    pseqs = [ msa[:,p-1] for p in pos]
    # list of row suquences limited to these positions
    id_st = {}
    for x,i in enumerate(zip(*pseqs)):
    	seq = ''.join(i).upper()
    	if seq in sts.keys():
    		v = sts[seq] 
    	else:
    		v = 'U'
    	id_st[ msa[x].id ] = v	
    return id_st

def genome_var(df):
	"""Return a dataframe of genomes and their shared and unique variants given a dataframe of annotated point mutations."""
	# only non-synonymous variants
	nsmv = df[ df.type == 'N']
	# no. of such variants
	N = len(nsmv)
	# join protein's name and amino acid change
	var_seqs = nsmv.apply( axis=1, func=lambda x:  
		pandas.Series( {'var':'_'.join( x[[ 'name','aa_change']]),'seqs':x['genomes'].split(',')}))
	# dict of variants and the list of sequences with the variant
	var_seqs = dict(zip(var_seqs['var'],var_seqs['seqs']))
	# set of all sequences
	allseqs = set( j for i in var_seqs.values() for j in i)  
	# initialize a dataframe with variants for rows and sequences for columns
	var_pa = pandas.DataFrame(index=var_seqs.keys(),columns=allseqs)
	# fill it with P/A of variant in the sequence
	var_pa = var_pa.apply( axis=1, func=lambda x: 
		pandas.Series( int(i in var_seqs[x.name]) for i in allseqs),result_type='broadcast') 
	# dataframe of shared and unique variants
	var_su = var_pa.apply(utils.get_shared_unique_elements,df=var_pa).T
	# dataframe of 
	var_nsu = var_su.apply( axis=1, func=lambda x: 
		pandas.Series([ (len(i)>0 and i.count(',')+1 or 0) for i in x],index=['nshared','nunique'])) 
	var_nsu['total'] = var_nsu.sum(axis=1)
	outdf = pandas.concat([var_nsu,var_su],axis=1)
	return outdf

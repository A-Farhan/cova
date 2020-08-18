import os, re, subprocess, matplotlib, seaborn, pandas
from . import utils, FEATURETABLE, GENOME, PROTNAMES, TYPEPOS, SEQTYPES
from time import time
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data.CodonTable import unambiguous_dna_by_id as codon_table

def rm_genome_w_stopm(fin,fout,finmsa,foutmsa):
	"""
	Remove all genomes with non-sense mutations from an input MSA.
	"""
	# table of variants
	vtab = utils.readcsv(fl=fin,sep='\t',header=True)[1]
	# dict of genome with their list of variants
	genome_vrs = { i[0]:i[4].split(',')+i[5].split(',') for i in vtab}
	# dict of genome with their list of nonsense mutations
	genome_stopm = {}
	
	for k,v in genome_vrs.items():
		stopms = [ i for i in v if bool(re.search('_[A-Z]\d+\*',i))]
		if len(stopms) > 0:
			genome_stopm[k] = ','.join(stopms)

	# list form of the dict of genomes with nonsense mutations 
	out = [ [k,v] for k,v in genome_stopm.items()]
	utils.writecsv(fl=fout, data=out, sep='\t')
	# input alignment 
	aln = AlignIO.read(finmsa,'fasta')
	# list of sequence records with no nonsense mutation
	alnls = [ i for i in aln if i.id not in genome_stopm.keys()]
	SeqIO.write(alnls,foutmsa,'fasta')


def extract_nucmsa_prots(msa, outdr, rfss=13468, prfs='YP_009725307.1', pepswstop = ['methyltransferase']):
	"""
	Extract nucleotide MSA of proteins from Reference-limited Whole-genome Multiple Sequence Alignment.
	
	Arguments:
	- msa 		- reference-limited Whole-genome Multiple Sequence Alignment (Biopython)
	- outdr 	- full path to the output directory of MSAs
	
	Optional arguments:
	- rfss 		- 1-indexed ribosomal frameshifting position 	[ 13468 ]
	- prfs 		- protein affected by the above frameshifting	[ YP_009725307.1' ]
	- pepswstop 	- list of peptides with stop codons 		[ 'methyltransferase' ]
	"""
	## checks
	# msa is a biopython MSA
	if type(msa) is not MultipleSeqAlignment:
		raise TypeError('msa must be a biopython MultipleSeqAlignment.')
	
	# create output directory, if not already present
	if not os.path.exists(outdr):
		os.mkdir(outdr)
		print("%s was not already present. Created now."%outdr)

	# dict of protein ids ( other than of polyproteins) and end points
	pid_ends = { i[10]:[ i[0], int(i[7]), int(i[8]), PROTNAMES[ i[13] ] ] \
	for i in FEATURETABLE if i[10] != '' and 'polyprotein' not in i[13]}

	for k,v in pid_ends.items():
		# output file for subMSA
		fmout = os.path.join( outdr, v[3]+'.msa')
	
		# does the region have stop codon
		if v[0] == 'CDS':
			# by default, CDSs have stop codons included
			has_stop=True
		# some peptides may have stop codon 
		elif v[3] in pepswstop:
			has_stop=True
		# other extracted sequences won't have stop codons
		else:
			has_stop=False

		# for the protein affected by frameshifting
		if k == prfs:
			cds = msa[:,v[1]-1:rfss] + msa[:,rfss-1:v[2]]
			for i in cds:
				i.description = ''
			AlignIO.write(alignments=[cds], handle=fmout, format='fasta')					
		# for regular proteins
		else:
			utils.extract_cds(begin=v[1], end=v[2], msa=msa, has_stop=has_stop, outf=fmout)

def annotate_var(fin,fout,ft=FEATURETABLE,genome=GENOME,codon_table=codon_table[1],\
	rfs_type=-1,rfs_pos=13468,rfs_prot='YP_009725307.1'):
	"""
	Annotate point mutations located within protein regions. Identify amino acid changes
	corresponding to nucleotide changes. The output table is saved to a file and not
	returned. It is a nested list of lists, where each entry is a list of the following
	7 elements.
	- 0 : protein id
	- 1 : protein name, as present in the reference feature table
	- 2 : 1-indexed genomic position
	- 3 : nucleotide variant
	- 4	: reference codon
	- 5 : variant codon
	- 6 : amino acid change
			A string with 3 components:
			- reference aa
			- 1-indexed aa position 
			- variant aa

	Arguments:
	- fin 			- full path to input file of point mutations table
	- fout 			- full path to output file of annotated mutations
	- ft 			- NCBI reference feature table 
	- genome 		- NCBI reference assembly genome
	- codon_table 	- Biopython codon table 	[ Standard]
	
	Additional arguments:
	To account for ribosomal frameshifting present in Coronavirus
	- rfs_type 		- type 				[ -1 ]
	- rfs_pos 		- 1-indexed genomic position 	[ 13468 ]
	- rfs_prot 		- protein id 			[ 'YP_009725307.1' ]
	"""
	# list of stop codons
	stopc = codon_table.stop_codons
	# point mutations variant table
	head, vartab = utils.readcsv(fl=fin,sep='\t',header=True) 
	# list of genomes
	genomes = head[2:]
	# subset feature table to columns of interest
	prot_ftrs = { i[10]:[ int(i[7]), int(i[8]), i[13]] for i in ft if i[0] in ['CDS','mat_peptide'] \
				and 'polyprotein' not in i[13]}

	## dict of proteins variants and genomes
	var_genomes = {}
	# for every row in the variant table
	for row in vartab:
		# 1-indexed variant position
		pos = int(row[0])
		# reference base
		rb = row[1]
		# find the affected protein(s)
		prots = [ k for k,v in prot_ftrs.items() if v[0] <= pos <= v[1] ]
		if len(prots) == 0:
			continue

		# set of tuples( prot, pos, ref, var nuc)
		pprn = set( ( prot, pos, rb, i) for i in row[2:] if i != rb for prot in prots)
		# dict with above tuples as keys and list of genome ids as their values
		entry = { k:[ genomes[x] for x,i in enumerate(row[2:]) if i == k[3] ] for k in pprn}
		var_genomes.update(entry)

	# list of all variant tuples
	vtups = list(var_genomes.keys())
	# dict of protein and their variants
	prot_vrs = utils.split_data(data=vtups, ix=0, cixs=[0,1,2,3])

	# initialize output table
	out = []
	# for every protein and its list of variants
	for prot, vrs in prot_vrs.items():
		# list of relevant features of the protein
		ftrs = prot_ftrs[prot]
		# 0-indexed ends
		b = ftrs[0]-1
		e = ftrs[1]

		# CDS sequence
		if prot == rfs_prot:
			cds_seq = genome[b:rfs_pos] + genome[rfs_pos+rfs_type:e]
		else:
			cds_seq = genome[b:e]
		
		# corresponding list of codons
		try:
			codonls = utils.n2c(cds_seq)
		except utils.LenSeqError:
			print("\tsequence of {} is not a multiple of 3 ( invalid CDS).".format(ftrs[2].upper()))
			continue

		# remove the terminal stop codon, if present
		#if codonls[-1] in stopc:
		#	codonls = codonls[:-1]

		# for every variant position
		for nvp in vrs:
			# 0-index of the genome position
			genome_x = nvp[1]-1
			# corresponding index in the CDS
			if prot == rfs_prot:
				if genome_x < (rfs_pos-1):
					cds_x = genome_x - b
				else:
					cds_x = genome_x - b - rfs_type
			else:
					cds_x = genome_x - b
			
			# list of amino acid variant(s)
			try:
				avs = utils.nv2av(p=cds_x, v=nvp[3], seq=codonls, codon_table=codon_table)
			except utils.LenSeqError:
				print('''Invalid variant {} \n Associated protein's features {}
					'''.format(nvp,prot_ftrs[prot]))
				raise 

			# make an entry in the output table
			entry = [ [ prot, PROTNAMES[ftrs[2]] ] + nvp[1:] + i + [','.join(var_genomes[tuple(nvp)])]\
				for i in avs]
			out.extend(entry)
	
	# add a column for no. of genomes and one to identify mutation type
	out = [ row[:-1] + [ 'N' if row[-2][0] != row[-2][-1] else 'S', row[-1].count(',') + 1, row[-1]] \
	for row in out]

	head = ['protein_id', 'name', 'position', 'ref_base', 'variant_base',\
			'old_codon','new_codon','aa_change', 'type', 'freq', 'genomes']
	utils.writecsv(fl=fout,data=out, header=head,sep='\t')

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
	# write both tables to files
	utils.writecsv(fl=frout, data=rates_out, header=['protein', 'exp_subs','syn', 'nonsyn', 'dnds'])
	utils.writecsv(fl=fsout, data=sites_out, header=['protein','site','syn', 'nonsyn', 'post_prob'])

def genome_seqtype(fin,fst=None,pos=TYPEPOS,sts=SEQTYPES):
    """
	Return Sequence type of genomes from the input alignment.
	By default, CoVa's sequence types are used. Optionally, a file can be provided for sequence types.
    """
    # biopython multiple sequence alignment
    aln = AlignIO.read(fin,'fasta')
    
    # if a sequence type file was provided
    if fst is not None:
    	print("A sequence type file was provided. Not using CoVa's standard sequence types")
    	pos, sts = utils.readcsv(fl=fst,header=True)
    	pos = [ int(i) for i in pos[1:]]
    	sts = { ''.join(i[1:]):i[0] for i in sts}
	
    # list of column sequences in the aligment at these positions
    pseqs = [ aln[:,p-1] for p in pos]
    # list of row suquences limited to these positions
    id_st = {}
    for x,i in enumerate(zip(*pseqs)):
    	seq = ''.join(i).upper()
    	if seq in sts.keys():
    		v = sts[seq] 
    	else:
    		v = 'U'
    	id_st[ aln[x].id ] = v	
    return id_st

def genome_var(fin):
	"""
	Write a table of genomes and their shared and unique annotated point mutations.

	Arguments:
	- fin - full path to the file of annotated variants table
	- fout - full path to output file
	"""
	## checks
	if not os.path.exists(fin):
		raise FileNotFoundError("input file %s must be present."%fin)
	
	# annotation table
	van = pandas.read_csv(fin,sep='\t')
	# only non-synonymous variants
	nsmv = van[ van.type == 'N']
	# no. of such variants
	N = len(nsmv)
	# join protein's name and amino acid change
	allvs = nsmv['name'] + '_' + nsmv['aa_change']
	# dict of genomes and variants
	genome_vrs = {}

	for x in range(N):
		v = allvs.iat[x]
		genomes = str(nsmv['genomes'].iat[x]).split(',')
		for genome in genomes:
			if genome not in genome_vrs.keys():
				genome_vrs[genome] = [v]
			else:
				genome_vrs[genome].append(v)
		
	# for every genome, get number and values of all, shared and unique variants
	genome_vls = {}

	for genome,var in genome_vrs.items():
		# no. of variants
		nvs = len(var)
		
		# set of variants in other genomes
		other_vs = set( i for k,v in genome_vrs.items() if k != genome for i in v)
		# shared variants
		shared_vs = set(var) & other_vs
		nsv = len(shared_vs)
		# unique varuants
		unq_vs = set(var) - other_vs
		nuv = len(unq_vs)
		# make an entry in the dict
		genome_vls[genome] = [ nvs, nsv, nuv, ','.join(shared_vs), ','.join(unq_vs)]

	return genome_vls

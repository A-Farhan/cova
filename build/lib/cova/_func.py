import os, re, subprocess, matplotlib, seaborn, ete3
from . import _utils, FEATURETABLE, GENOME, PROTNAMES, LIBPATH 
from time import time
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data.CodonTable import unambiguous_dna_by_id as codon_table

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

		action = True
		if os.path.exists(fmout):
			response=input('''output file %s already exists. 
Do you want to proceed and rewrite it? [y/n]\n'''%fmout)
			if response == 'n':
				action = False
			elif response == 'y':
				action = True
			else:
				raise ValueError('''Inappropriate input! Please respond with [y] to
proceed and overwrite OR with [n] to skip this protein and retain the MSA file''')
			
		if not action:
			print("okay! Skipping this protein.")
			continue
	
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
			_utils.extract_cds(begin=v[1], end=v[2], msa=msa, has_stop=has_stop, outf=fmout)

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
	- codon_table 		- Biopython codon table 	[ Standard]
	
	Additional arguments:
	To account for ribosomal frameshifting present in Coronavirus
	- rfs_type 		- type 				[ -1 ]
	- rfs_pos 		- 1-indexed genomic position 	[ 13468 ]
	- rfs_prot 		- protein id 			[ 'YP_009725307.1' ]
	"""
	# list of stop codons
	stopc = codon_table.stop_codons
	# data
	head, vt = _utils.readcsv(fl=fin,sep='\t',header=True) # variant table
	# subset feature table to columns of interest
	prot_ftrs = { i[10]:[ int(i[7]), int(i[8]), i[13]] for i in ft if i[0] in ['CDS','mat_peptide'] \
				and 'polyprotein' not in i[13]}

	## dict of proteins and their variants
	prot_vrs = {}
	for prot,ftrs in prot_ftrs.items():
		als = []
		for j in vt:
			# 1-indexed variant position
			pos = int(j[0])
			if ftrs[0] <= pos <= ftrs[1]:
				entry = [ pos, j[1], set(k for k in j[2:] if k != j[1])]
				als.append(entry)
		if len(als) > 0:
			prot_vrs[prot] = als

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
			codonls = _utils.n2c(cds_seq)
		except _utils.LenSeqError:
			print("\tsequence of {} is not a multiple of 3 ( invalid CDS).".format(ftrs[2].upper()))
			continue

		if codonls[-1] in stopc:
			codonls = codonls[:-1]

		# for every variant position
		for nvp in vrs:
			# 0-index of the genome position
			genome_x = nvp[0]-1
			
			# corresponding index in the CDS
			if prot == rfs_prot:
				if genome_x < (rfs_pos-1):
					cds_x = genome_x - b
				else:
					cds_x = genome_x - b - rfs_type
			else:
					cds_x = genome_x - b

			# for every nucleotide variant at this position
			for nv in nvp[2]:
				# amino acid variant
				av = _utils.nv2av(p=cds_x, v=nv, seq=codonls, codon_table=codon_table)
				entry = [ [ prot, ftrs[2], nvp[0], nv] + i for i in av]
				out.extend(entry)
		
	head = ['protein_id', 'name', 'genomic position','variant','old_codon','new_codon','aa_change']
	_utils.writecsv(fl=fout,data=out, header=head,sep='\t')

def plottree(ftree,fplot,fmap=None):
	"""
	Plot phylogeny tree.

	Arguments:
	- ftree - input tree file
	- fplot - output plot file 
	- fmap 	- (optional) file for metadata on sequences, used to annotate tree
	"""
	## checks
	# tree file is present
	if not os.path.exists(ftree):
		raise FileNotFoundError('tree file %s must be present.'%ftree)
	# plot file has png suffix
	if fplot.split('.')[-1] != 'png':
		raise ValueError('output file must have suffix "png".')
	# if map file is provided, it is a valid file path
	if fmap is not None:
	 	if type(fmap) is not str:
	 		raise TypeError('file name must be a string.')
	 	else:
	 		if not os.path.exists(fmap):
	 			raise FileNotFoundError("couldn't find the file %s."%fmap)

	# load tree
	t = ete3.Tree(ftree)
	# list of leaves
	leaves = t.get_leaf_names()

	# create treestyle
	ts = ete3.TreeStyle()
	ts.scale = 100000
	ts.branch_vertical_margin = -1
	
	# If mapping is available, then use it to color leaves and branches
	if fmap is not None:
		dmap = _utils.readcsv(fl=fmap,header=True)[1]
	
		# colors for countries
		isol_country = { i[0]:i[5] for i in dmap}
		countries = sorted( list( set( isol_country.values())))
		nc = len(countries)
		country_colors = _utils.get_N_HexCol(N=nc)
		isol_country_color = {}
		for k,v in isol_country.items():
			x = countries.index(v)
			c = country_colors[x]
			isol_country_color[k] = c

		# colors for months
		isol_month = { i[0]:'-'.join(i[-2].split('-')[:2]) for i in dmap}
		# months
		months = sorted( list( set( isol_month.values())))
		# dict of month names and corresponding key
		month_key = {}
		for x,i in enumerate(months):
			month_key[i] = x+1
		# replace months with key in the above
		isol_mkey = { k:month_key[v] for k,v in isol_month.items()}
		months = sorted( list(set(isol_mkey.values())))
		nm = len(months)
		month_colors = seaborn.color_palette('Blues',n_colors=nm)   
		isol_month_color = {}
		for k,v in isol_mkey.items():
			x = months.index(v)
			c = month_colors[x]
			isol_month_color[k] = matplotlib.colors.to_hex(c)

		for n in t.traverse():
			if n.is_leaf():	
				ns = ete3.NodeStyle()
				ns['fgcolor'] = isol_country_color[n.name]
				ns['size'] = 10
				ns['hz_line_width'] = 4
				ns['hz_line_color'] = isol_month_color[n.name]
				n.set_style(ns)

	# output
	t.render(fplot,h=900,w=1000,units='px',tree_style=ts)

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
	
	action = True
	if os.path.exists(fsout):
		response=input('''output file %s already exists. 
Do you want to proceed and rewrite it? [y/n]\n'''%fsout)
		if response == 'n':
			action = False
		elif response == 'y':
			action = True
		else:
			raise ValueError('''Inappropriate input! Please respond with [y] to
proceed and overwrite OR with [n] to skip this protein and retain the MSA file''')
		
	if not action:
		print("okay! Exiting.")
		return

	# json output intermediate file
	fjout_int = os.path.join( indr, f +'.FUBAR.json')
	# json output file
	fjout = os.path.join( outdr, p+'.json')
	# cache file
	fcout = os.path.join( outdr, p+'.cache')

	cmd = [ prog, 'fubar', 'LIBPATH='+LIBPATH, '--alignment', fmsa, '--tree', ftree, '--cache', fcout]
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

	action = True
	if os.path.exists(frout) or os.path.exists(fsout):
		response = input('''output file(s) already exist. 
Do you want to proceed and rewrite? [y/n]\n''')
		if response == 'n':
			action = False
		elif response == 'y':
			action = True
		else:
			raise ValueError('''Inappropriate input! Please respond with [y] to
proceed and overwrite OR with [n] to skip this protein and retain the MSA file''')
		
	if not action:
		print("okay! Exiting.")
		return
	else:
		print("okay! Proceeding with rewriting")		
	
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
			print("\tno changes in %s."%p)
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
	_utils.writecsv(fl=frout, data=rates_out, header=['protein', 'exp subs','syn', 'nonsyn', 'dnds'])
	_utils.writecsv(fl=fsout, data=sites_out, header=['protein','site','syn', 'nonsyn', 'post. prob'])

def genome_var(fpm,fann,fout):
	"""
	Write a table of genomes and their shared and unique annotated point mutations.

	Arguments:
	- fpm - full path to the file of point mutations table
	- fann - full path to the file of annotated variants table
	- fout - full path to output file
	"""
	## checks
	if not os.path.exists(fpm):
		raise FileNotFoundError("input file %s must be present."%fpm)
	if not os.path.exists(fann):
		raise FileNotFoundError("input file %s must be present."%fann)
	# if output is alaready present
	action = True
	if os.path.exists(fout):
		response=input('''output file %s already exists. 
Do you want to proceed and rewrite it? [y/n]\n'''%fout)
		if response == 'n':
			action = False
		elif response == 'y':
			action = True
		else:
			raise ValueError('''Inappropriate input! Please respond with [y] to
proceed and overwrite OR with [n] to terminate this command''')

	if not action:
		print("okay! Exiting.")
		return
	else:
		print("okay! Proceeding with rewriting")

	head_vtab, vtab = _utils.readcsv( fl=fpm, sep='\t', header=True)
	head_van, van = _utils.readcsv( fl=fann, sep='\t', header=True)
	# list of genome ids
	genomes = head_vtab[2:]
	# create a dict of genomes and all variants
	genome_vrs = { genome:[ ( row[0], row[x+2]) for row in vtab if row[1] != row[x+2] ] \
					for x,genome in enumerate(genomes) }
	# convert annotation table to a dict
	var_ann = { tuple(i[2:4]):PROTNAMES[i[1]] +'_'+ i[-1] for i in van \
				if len(set(re.split('\d+',i[-1]))) > 1}
	# use above to get a dict of genomes and annotated variants
	genome_anv = { k:[ var_ann[i] for i in v if i in var_ann.keys()] for k,v in genome_vrs.items() }
	# pool annotated variants
	all_anv = set( j for i in genome_anv.values() for j in i)
	# for every genome, get number and values of all, shared and unique variants
	genome_vls = {}

	for genome,anv in genome_anv.items():
		# all variants
		n_allv = len(genome_vrs[genome])
		
		# if genome has no variants, skip
		if n_allv == 0:
			print("\t%s has no variants."%genome)
			continue

		# set of variants in other genomes
		other_vs = set(i for k,v in genome_anv.items() if k != genome for i in v)
		# shared variants
		shared_vs = set(anv) & other_vs
		nsv = len(shared_vs)
		# unique varuants
		unq_vs = set(anv) - other_vs
		nuv = len(unq_vs)
		# make an entry in the dict
		genome_vls[genome] = [ n_allv, nsv, nuv, shared_vs, unq_vs]

	# process output table
	out = [ [k] + v[:3] + [ ','.join(v[3]), ','.join(v[4])] for k,v in genome_vls.items() ]
	# write to file
	_utils.writecsv(fl=fout, data=out, sep='\t',\
		header=['genome','#variants','#shared', '#unique', 'shared', 'unique'])
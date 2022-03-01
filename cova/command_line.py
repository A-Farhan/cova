import os, sys, shutil, click, subprocess, time, Bio, cova, pandas

# start clock
start = time.time()

### commands ###
@click.group(chain=True)
@click.version_option()
@click.option('--indr', help='Full path to the working directory', default=os.getcwd(), type=click.Path())
@click.option('--ncpu', help='Number of CPUs to use', default=4, show_default=True, type=int)
@click.option('--ref', help='Reference sequence accession', default=cova.REF, show_default=True)
@click.option('--debug', help='See full traceback for errors.', is_flag=True)
@click.option('--addseq', help='Add new sequences and redo analysis.', is_flag=True)
@click.pass_context
def cli(ctx,indr,ncpu,ref,debug,addseq):
	"""
	Variant analysis using whole-genome Multiple Sequence Alignments.
	By default, it works on the current directory. Directory can be specified using INDR option.
	
	New sequences can be added to update an existing analysis using --addseq. 
	"""
	ctx.ensure_object(dict)
	ctx.obj['DR'] = indr
	ctx.obj['REF'] = ref
	ctx.obj['NCPU'] = str(ncpu)
	ctx.obj['ADDSEQ'] = addseq
	# control traceback
	if debug:
		ctx.obj['DEBUG'] = debug
		click.echo("Debug mode is ON.\n")
	else:
		sys.excepthook = lambda exctype,exc,traceback : print("{}: {}".format(exctype.__name__,exc)) 
	click.echo("CoVa will run in the directory: {}\n".format(indr))

@cli.command()
@click.pass_context
@click.option('--prog', default='mafft', 
	help='''full path to MAFFT program''',show_default=True)
@click.option('--infile', default='genomes.fna', show_default=True)
@click.option('--outfile', default='genome_aln.fna', show_default=True)
@click.option('--mode', type=click.Choice(['standard','fast','ultra'],case_sensitive=False),default='standard',show_default=True)
def msabuild(ctx,prog,infile,outfile,mode):
	"""Build whole-genome MSA.""" 
	fin = os.path.join(ctx.obj['DR'],infile)
	fout = os.path.join(ctx.obj['DR'],outfile)
	
	# throw error if input file missing
	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)

	if cova.utils.outcheck(fout):
		# get number of sequences in the input
		nseq = len(Bio.SeqIO.index( fin, 'fasta'))
		print("%s: Input has %i sequences."%( cova.utils.timer(start), nseq))

		# set path variable to find mafft
		my_env = os.environ
		if 'COVA_BIN_PATH' in my_env.keys():
			my_env['PATH'] = ':'.join([ my_env['COVA_BIN_PATH'], my_env['PATH']])
		
		# create commands for different run modes ( FAST / ULTRA / STANDARD)
		# FAST
		if mode == 'fast':
			cmd = [prog, '--quiet', '--retree', '2', '--thread', ctx.obj['NCPU'], fin]
		# ULTRA
		elif mode == 'ultra':
			# with ULTRA mode, a reference genome file is first placed in the project dir
			fref = os.path.join(ctx.obj['DR'],'ref.fasta')
			Bio.SeqIO.write(cova.GENOME,fref,'fasta')
			cmd = [prog, '--quiet', '--auto', '--thread', ctx.obj['NCPU'], '--keeplength', '--addfragments', fin, fref]
		# STANDARD
		else:
			cmd = [prog, '--quiet', '--nomemsave', '--maxiterate', '5', '--thread', ctx.obj['NCPU'], fin]

		# run the MAFFT command created above
		print("%s: Building MSA from %s in %s mode\n Command: %s,\n Output will be saved to %s"%(\
			cova.utils.timer(start),fin,mode,' '.join(cmd),fout))
		
		with open( fout,'w') as flob:
			s1 = subprocess.run( cmd, stdout=flob, env=my_env)

		# clean-up for ultra mode
		if mode == 'ultra':
			# remove reference file
			if os.path.exists(fref):
				os.remove(fref)
			# remove additional reference seq from MSA
			msa = Bio.AlignIO.read(fout, 'fasta')
			msa = msa[1:,:]
			Bio.AlignIO.write(msa, fout, 'fasta')

	print("%s:\tMSABUILD is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--prog', default='mafft', 
	help='''full path to MAFFT program''',show_default=True)
@click.option('--inmsa', default='genome_aln.fna', show_default=True)
@click.option('--newseq', default='new_seq.fna', show_default=True)
@click.option('--oldcopy', default='old_genome_aln.fna', show_default=True)
def msad(ctx,prog,inmsa,newseq,oldcopy):
	"""Add new sequence(s) to a pre-existing whole-genome MSA.""" 
	fin1 = os.path.join(ctx.obj['DR'],inmsa)
	fin2 = os.path.join(ctx.obj['DR'],newseq)
	fout = fin1
	fcopy = os.path.join(ctx.obj['DR'],oldcopy)
	
	if not os.path.exists(fin1):
		raise FileNotFoundError("couldn't read the input file %s."%fin1)

	if not os.path.exists(fin2):
		raise FileNotFoundError("couldn't read the input file %s."%fin2)

	# set path variable to find mafft
	my_env = os.environ
	if 'COVA_BIN_PATH' in my_env.keys():
		my_env['PATH'] = ':'.join([ my_env['COVA_BIN_PATH'], my_env['PATH']])

	# first, copy the original file
	print("Generating backup for the original MSA")
	shutil.copy(src=fin1, dst=fcopy)
	# then run the addseq command: FFT-NS-2
	cmd = [prog, '--quiet', '--auto', '--thread', ctx.obj['NCPU'], '--addfragments', fin2, fcopy]
	print("%s: Adding sequences from %s to %s,\n Command: %s,\n Output will be saved to %s"%(\
		cova.utils.timer(start),fin2, fcopy,' '.join(cmd),fout))
	
	# rewrite MSA file
	with open( fout,'w') as flob:
		s1 = subprocess.run( cmd, stdout=flob, env=my_env)
	print("%s:\tMSAD is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='genome_aln.fna',show_default=True)
@click.option('--outfile',default='genome_aln_ref.fna',show_default=True)
def msaref(ctx,infile,outfile):
	"""Limit MSA to sites present in the reference."""
	fin = os.path.join(ctx.obj['DR'], infile)
	fout = os.path.join(ctx.obj['DR'],outfile)
	
	# throw error if input missing
	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
	
	# if output doesn't exist	
	if cova.utils.outcheck(fout):
		# cova MSA object
		msa = cova.utils.MSA(fname=fin,ref=ctx.obj['REF'])
		print("{}: Generating reference limited MSA, Output will be saved to {}\n".format(cova.utils.timer(start),fout))
		# use its method to limit MSA to a reference sequence
		out = msa.limref()
		# write output
		Bio.AlignIO.write(alignments=[out], handle=fout, format='fasta')

	print("%s:\tMSAREF is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile', default='genome_aln_ref.fna', show_default=True)
@click.option('--outfile1',default='genome_aln_unq.fna', show_default=True)
@click.option('--outfile2',default='genome_dups.tsv', show_default=True)

def msaunq( ctx, infile, outfile1, outfile2):
	"""Remove duplicate sequences from reference-limited MSA."""
	fin = os.path.join(ctx.obj['DR'], infile)
	fout1 = os.path.join(ctx.obj['DR'],outfile1)
	fout2 = os.path.join(ctx.obj['DR'],outfile2)

	# throw error if input missing
	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
	
	if cova.utils.outcheck(fout1):
		# cova MSA object
		msa = cova.utils.MSA(fname=fin,ref=ctx.obj['REF'])
		print("{}: Removing duplicate sequences from {}\n".format(cova.utils.timer(start),fout1))
		out1, out2 = msa.rmdup()
		# write alignment to output path
		Bio.AlignIO.write(alignments=[out1], handle=fout1, format='fasta')
		# write dataframe of genomes retained and their excluded duplicates
		out2 = pandas.DataFrame(out2)
		out2.to_csv(fout2,sep='\t',index=False,header=False)
	
	print("%s:\tMSAUNQ is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='genome_aln_unq.fna',show_default=True)
@click.option('--typefile',default=None,help='Full path to the file with sequence types definition [Optional]')
@click.option('--outfile',default='genome_types.csv',show_default=True)
def seqtype( ctx, infile, typefile, outfile):
	"""Identify sequence types."""
	fin = os.path.join(ctx.obj['DR'], infile)
	fout = os.path.join(ctx.obj['DR'], outfile)

	if cova.utils.outcheck(fout):
		print("{}: Sequence Typing genomes from {}".format(cova.utils.timer(start),fin))
		# biopython multiple sequence alignment
		msa = Bio.AlignIO.read(fin,'fasta')
		# genomes and their sequence types
		out = [ [k,v] for k,v in cova.genome_seqtype(msa,fst=typefile).items()]
		out = sorted(out, key=lambda x: x[1])
		# covert to df
		out = pandas.DataFrame(out)
		# write table to file
		out.to_csv(fout,header=False,index=False)		

	print("%s:\t SEQTYPE is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='genome_aln_unq.fna',show_default=True)
@click.option('--outfp',default='point_mutations.tsv',show_default=True)
@click.option('--outfd',default='deletions.tsv',show_default=True)
def vcalpd(ctx,infile,outfp,outfd):
	"""Call point mutations / deletions from Reference-limited MSA."""
	fin = os.path.join(ctx.obj['DR'], infile)
	fout1 = os.path.join(ctx.obj['DR'],outfp)
	fout2 = os.path.join(ctx.obj['DR'],outfd)

	# throw error if input missing
	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
	
	msa = cova.utils.MSA(fname=fin,ref=ctx.obj['REF'])
	
	# point mutations	
	if cova.utils.outcheck(fout1):
		print("{}: Calling point mutations from {}\n".format(cova.utils.timer(start),fin))
		tab1 = msa.pointmuts()
		tab1.to_csv(fout1,sep='\t',index=False)
	
	# deletions
	if cova.utils.outcheck(fout2):
		print("{}: Calling deletions from {}\n".format(cova.utils.timer(start),fin))
		tab2 = msa.dels()
		tab2.to_csv(fout2,sep='\t',index=False)
	
	print("%s:\t VCALPD is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile1',default='point_mutations.tsv',show_default=True)
@click.option('--outfile1',default='point_mutations_ann.tsv',show_default=True)
@click.option('--infile2',default='deletions.tsv',show_default=True)
@click.option('--outfile2',default='deletions_ann.tsv',show_default=True)
def anvpd(ctx,infile1,outfile1,infile2,outfile2):
	"""Annotate point mutations and deletions located within protein regions."""
	fin1 = os.path.join(ctx.obj['DR'], infile1)
	fout1 = os.path.join(ctx.obj['DR'],outfile1)
	fin2 = os.path.join(ctx.obj['DR'], infile2)
	fout2 = os.path.join(ctx.obj['DR'],outfile2)

	# throw error if input missing
	if not os.path.exists(fin1):
		raise FileNotFoundError("couldn't read the input file %s."%fin1)

	if not os.path.exists(fin2):
		raise FileNotFoundError("couldn't read the input file %s."%fin2)
	
	# annotate point mutations
	if cova.utils.outcheck(fout1):
		print("%s: Annotating point mutations located within protein regions"%cova.utils.timer(start))
		cova._annotator.annotate_pm(fin1,fout1)
	
	# annotate deletions
	if cova.utils.outcheck(fout2):
		print("%s: Annotating deletions located within protein regions"%cova.utils.timer(start))
		cova._annotator.annotate_del(fin2,fout2)

	print("%s:\t ANVPD is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='point_mutations_ann.tsv',show_default=True)
@click.option('--outfile',default='genome_vars.tsv',show_default=True)
def nsvar( ctx, infile, outfile):
	"""Get shared and unique non-synonymous variants for genomes."""
	fin = os.path.join(ctx.obj['DR'], infile)
	fout = os.path.join(ctx.obj['DR'], outfile)
	
	if cova.utils.outcheck(fout):
		print("%s: Identifying shared and unique variants"%cova.utils.timer(start))
		# dataframe of annotated variants
		van = pandas.read_csv(fin,sep='\t')
		# dataframe of shared and unique variants
		out = cova.genome_var(van)
		# write to output path
		out.to_csv(fout,sep='\t',index_label='id')
		
	print("%s:\t NSVAR is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='genome_vars.tsv',show_default=True)
@click.option('--outfile',default='genome_w_stopm.tsv',show_default=True)
@click.option('--inmsafile',default='genome_aln_unq.fna',show_default=True)
@click.option('--outmsafile',default='genome_aln_sf.fna',show_default=True)
def rmstop(ctx,infile,outfile,inmsafile,outmsafile):
	"""Remove genomes with non-sense mutations."""
	## paths
	fin = os.path.join(ctx.obj['DR'], infile)
	fout = os.path.join(ctx.obj['DR'], outfile)
	finmsa = os.path.join(ctx.obj['DR'], inmsafile)
	foutmsa = os.path.join(ctx.obj['DR'], outmsafile)
	
	# shall we proceed in case output is present?
	if cova.utils.outcheck(fout):
		# table of variants
		vtab = pandas.read_csv(fin,sep='\t',index_col=0)
		print("%s: Identifying and removing genomes with non-sense mutations."%cova.utils.timer(start))
		# dataframe of genomes with nonsense variants
		out = cova.rm_genome_w_stopm(vtab)
		# write dataframe to output path
		out.to_csv(fout, sep='\t',header=False)
		# input alignment 
		aln = Bio.AlignIO.read(finmsa,'fasta')
		# list of sequence records with no nonsense mutation
		alnls = [ i for i in aln if i.id not in out.index]
		# write alignment to output path
		Bio.SeqIO.write(alnls,foutmsa,'fasta')

	print("%s:\t RMSTOP is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='genome_aln_sf.fna',show_default=True)
@click.option('--outdr',default='prots_nmsa',show_default=True)
def msap(ctx,infile,outdr):
	"""Extract nucleotide MSA of proteins from Reference-limited MSA."""
	fin = os.path.join(ctx.obj['DR'], infile)
	dout = os.path.join(ctx.obj['DR'], outdr)

	# throw error if input missing
	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
	
	# create output directory, if not already present
	if not os.path.exists(dout):
		os.mkdir(dout)
		print("%s was not already present. Created now."%outdr)

	print('''%s: Extracting nucleotide MSA of proteins from Reference-limited MSA,
		Outputs will be saved to %s\n.'''%(cova.utils.timer(start),dout))
	cova.extract_nmsa_prots(fmsa=fin, dr=dout)
	
	print("%s:\t MSAP is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='genome_aln.fna',show_default=True)
@click.option('--outfile',default='insertions.tsv',show_default=True)
def vcali(ctx,infile,outfile):
	"""Call insertions from MSA."""
	fin = os.path.join(ctx.obj['DR'], infile)
	fout = os.path.join(ctx.obj['DR'],outfile)

	# throw error if input missing
	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
	
	if cova.utils.outcheck(fout):
		msa = cova.utils.MSA(fname=fin,ref=ctx.obj['REF'])
		print("%s: Calling insertions from MSA , Output will be saved to %s\n"%(cova.utils.timer(start),fout))
		tab = msa.ins()
		tab.to_csv(fout,sep='\t',index_label='pos')
		
	print("%s:\t VCALI is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='genome_aln_sf.fna',show_default=True)
@click.option('--indr',default='prots_nmsa',show_default=True)
@click.option('--window',default=300,show_default=True,type=int)
@click.option('--jump',default=20,show_default=True,type=int)
@click.option('--outfile1',default='divs.csv',show_default=True)
@click.option('--outfile2',default='slide_divs.csv',show_default=True)
@click.option('--slide',is_flag=True,help='Should we calculate sliding diversity?')
def div(ctx,infile,indr,window,jump,outfile1,outfile2,slide):
	"""Compute nucleotide diversity from Reference-limited MSA."""
	fin = os.path.join(ctx.obj['DR'], infile)
	din = os.path.join(ctx.obj['DR'], indr)
	fout1 = os.path.join(ctx.obj['DR'],outfile1)
	fout2 = os.path.join(ctx.obj['DR'],outfile2)
	ncpu = int(ctx.obj['NCPU'])

	# throw error if input missing
	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
	
	# load alignment as cova MSA object
	msa = cova.utils.MSA(fname=fin)

	if cova.utils.outcheck(fout1):
		print('''%s: Computing diversity for whole-genome and peptide-encoding regions.
		Output will be saved to %s\n'''%(cova.utils.timer(start),fout1))
		wndiv = msa.ndiv()
		fpmsas = [ i for i in os.listdir(din) if i.endswith('.msa')]
		pndivs = [ [ i.replace('.msa',''), cova.utils.MSA(os.path.join(din,i)).ndiv()] for i in fpmsas]
		pndivs = [ i for i in pndivs if i[1] is not None]
		pndivs = sorted( pndivs, key=lambda x: x[1], reverse=True)
		out1 = [ ['genome', wndiv] ] + pndivs
		# save dataframe to output file
		out1 = pandas.DataFrame(out1)
		out1.to_csv( fout1, index=False, header=False)
	
	# to compute diversity over sliding window
	if slide:
		if cova.utils.outcheck(fout2):
			print('''%s: Computing diversity within a window sliding over the genome.
			Output will be saved to %s\n'''%(cova.utils.timer(start),fout2))
			out2 = msa.slide_ndiv(window=window,jump=jump,ncpu=ncpu)
			# save dataframe to output file
			out2 = pandas.DataFrame(out2)
			out2.to_csv( fout2, index=False, header=False)
	
	print("%s:\t DIV is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--prog', default='FastTree', 
	help='''full path to FASTTREE program''', show_default=True)
@click.option('--infile', default='genome_aln_sf.fna', show_default=True)
@click.option('--outfile', default='tree.nwk', show_default=True)
@click.option('--boots', help='''# Bootstrap samples''', default=100, show_default=True, type=int)
def tree(ctx,prog,infile,outfile,boots):
	"""Build phyogeny from whole-genome MSA."""
	fin = os.path.join( ctx.obj['DR'], infile)
	fout = os.path.join( ctx.obj['DR'], outfile)
	bs = str(boots)

	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
	
	# building
	if cova.utils.outcheck(fout):
		my_env = os.environ
		my_env['OMP_NUM_THREADS'] = ctx.obj['NCPU']
		
		# set path variable to find FASTTREE
		if 'COVA_BIN_PATH' in my_env.keys():
			my_env['PATH'] = ':'.join([ my_env['COVA_BIN_PATH'], my_env['PATH']])
		
		# create FASTTREE command
		cmd = [prog, '-quiet', '-nt', '-mlnni', '4', '-boot', bs, fin]
		print("%s: Building Phylogeny from %s,\n Command: %s,\n Output will be saved to %s"%(\
			cova.utils.timer(start),fin,' '.join(cmd),fout))
		
		# run the FASTTREE command create above
		with open( fout,'w') as flob:
			s1 = subprocess.run( cmd, stdout=flob, stderr=subprocess.DEVNULL, env=my_env)
		
	print("%s:\t TREE is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--prog', default='hyphy',
	help='''full path to HYPHY program''', show_default=True)
@click.option('--indr', default='prots_nmsa', show_default=True)
@click.option('--tree', default='tree.nwk', show_default=True)
@click.option('--outdr', default='fubar', show_default=True)
@click.option('--outr', default='rates.csv', show_default=True)
@click.option('--outs', default='sites.csv', show_default=True)
@click.option('--excl', default='orf1a,orf1ab', help='comma-separated list of proteins to be excluded from this analysis.',
	show_default=True)
def sel(ctx,prog,tree,indr,outdr,outr,outs,excl):
	"""Identify sites under positive selection."""
	ftree = os.path.join(ctx.obj['DR'], tree)
	din   = os.path.join(ctx.obj['DR'], indr)
	dout  = os.path.join(ctx.obj['DR'],outdr)
	frout = os.path.join(ctx.obj['DR'], outr)
	fsout = os.path.join(ctx.obj['DR'], outs)
	excl = excl.split(',')

	if not os.path.exists(ftree):
		raise FileNotFoundError("couldn't read the input tree file %s."%ftree)
	
	# set path variable to find hyphy
	my_env = os.environ
	if 'COVA_BIN_PATH' in my_env.keys():
		my_env['PATH'] = ':'.join([ my_env['COVA_BIN_PATH'], my_env['PATH']])

	# check if hyphy can find its batch files
	cmd = [prog, 'fubar', '--help']
	s = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, env=my_env)
	if s.returncode != 0:
		raise FileNotFoundError("Hyphy couldn't find where FUBAR is. Check library path")

	# if we shall proceeed in case the output directory already exists
	if cova.utils.outcheck(dout):
		
		# if the directory did exist before, delete it
		if os.path.exists(dout):
			shutil.rmtree(dout)
		
		# make the directory
		os.mkdir(dout)

		print('''{}: Analysing MSAs for positive selection using tree ({})
			Outputs will be saved to {}'''.format(cova.utils.timer(start),ftree,dout))
		
		for i in os.listdir(din):
			
			if i.endswith('.msa') and i.split('.')[0] not in excl:
				cova.run_fubar(fmsa=os.path.join(din,i), ftree=ftree, outdr=dout, prog=prog)
	
	if cova.utils.outcheck(frout):		
		print('''%s: Parsing FUBAR output to generate rates and sites tables'''%cova.utils.timer(start))
		cova.parse_fubar(indr=dout, frout=frout, fsout=fsout)
	
	print("%s:\t SEL is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--excl', default=None, help='list of files to be saved from cleanup [Optional]')

def cleanup(ctx,excl):
	"""Remove I/O files processed by CoVa present in the working directory."""	
	# directory
	dr = ctx.obj['DR']
	# list of files to be removed from this directory
	fls = cova.IOF
	
	# list of files to be saved from this process
	if excl is not None:
		exls = excl.split(',')
	else:
		exls = []

	# for each file in the list
	for f in fls:
		# full path of the file
		fp = os.path.join( dr, f)
		# if this path is to be retained
		if f in exls:
			continue

		# if path already exists
		if os.path.exists(fp):
			
			# and if its a directory, remove the directory
			if os.path.isdir(fp):
				print("\tDeleting dir: %s"%fp)
				shutil.rmtree(fp)
			# else, if it's a file, delete the file
			else:
				print("\tDeleting file: %s"%fp)
				os.remove(fp)

	print("%s:\t CLEANUP is done."%cova.utils.timer(start))

### command to run all other commands
@cli.command()
@click.pass_context
def full(ctx):
	"""
	Run full pipeline.
	Same as running: CoVa msabuild msaref msaunq seqtype vcalpd anvpd nsvar rmstop msap vcali div tree sel
	"""
	# if new sequences are to be added
	if ctx.obj['ADDSEQ']:
		# add new sequences
		ctx.forward(msad)
		# clean up
		ctx.invoke(cleanup,excl='genome_aln.fna')
		print('\t Done with cleanup.')
	else:
		ctx.forward(msabuild)
	
	ctx.forward(msaref)
	ctx.forward(msaunq)
	ctx.forward(seqtype)
	ctx.forward(vcalpd)
	ctx.forward(anvpd)
	ctx.forward(nsvar)
	ctx.forward(rmstop)
	ctx.forward(msap)
	ctx.forward(vcali)
	ctx.forward(div)
	ctx.forward(tree)
	ctx.forward(sel)

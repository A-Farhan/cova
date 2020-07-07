import os, sys, shutil, click, subprocess, time, Bio, cova

# start clock
start = time.time()

### commands ###
@click.group(chain=True)
@click.version_option()
@click.option('--indr', help='Full path to the working directory', default=os.getcwd(), type=click.Path())
@click.option('--ncpu', help='number of CPUs to use', default=4, show_default=True, type=int)
@click.option('--debug', help='See full traceback for errors.', is_flag=True)
@click.option('--addseq', help='Add new sequences and redo analysis.', is_flag=True)
@click.pass_context
def cli(ctx,indr,ncpu,debug,addseq):
	"""
	Variant analysis using whole-genome Multiple Sequence Alignments.
	By default, it works on the current directory. Directory can be specified using INDR option.
	
	New sequences can be added to update an existing analysis using --addseq. 
	"""
	ctx.ensure_object(dict)
	ctx.obj['DR'] = indr
	ctx.obj['REF'] = cova.REF
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
@click.option('--maxseq_accuracy', type=int, default=1000, show_default=True)
def msabuild(ctx,prog,infile,outfile,maxseq_accuracy):
	"""Build whole-genome MSA.""" 
	fin = os.path.join(ctx.obj['DR'],infile)
	fout = os.path.join(ctx.obj['DR'],outfile)
	
	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)

	if cova.utils.outcheck(fout):
		# set path variable to find mafft
		my_env = os.environ
		if 'COVA_BIN_PATH' in my_env.keys():
			my_env['PATH'] = ':'.join([ my_env['COVA_BIN_PATH'], my_env['PATH']])

		# no. of sequences
		nseqs = len(Bio.SeqIO.index(fin,'fasta'))
		click.echo("{}: {} sequences found in the input file.".format(cova.utils.timer(start),nseqs))
		
		if nseqs < maxseq_accuracy:
			cmd = [prog, '--quiet', '--nomemsave', '--maxiterate', '5', '--thread', ctx.obj['NCPU'], fin]
		else:
			click.echo("\ttoo many sequences! MSABUILD will optimize for speed and run FFT-NS-2.")
			cmd = [prog, '--quiet', '--retree', '2', '--thread', ctx.obj['NCPU'], fin]

		print("%s: Building MSA from %s,\n Command: %s,\n Output will be saved to %s"%(\
			cova.utils.timer(start),fin,' '.join(cmd),fout))
		
		with open( fout,'w') as flob:
			s1 = subprocess.run( cmd, stdout=flob, env=my_env)

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
	cmd = [prog, '--quiet', '--retree', '2', '--thread', ctx.obj['NCPU'], '--add', fin2, fcopy]
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
	
	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
		
	if cova.utils.outcheck(fout):
		msa = cova.utils.MSA(fname=fin)
		print("{}: Generating reference limited MSA, Output will be saved to {}\n".format(cova.utils.timer(start),fout))
		out = msa.limref(ref=ctx.obj['REF'])
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

	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
	
	if cova.utils.outcheck(fout1):
		msa = cova.utils.MSA(fname=fin)
		print("{}: Removing duplicate sequences from {}\n".format(cova.utils.timer(start),fout1))
		out1, out2 = msa.rmdup(ref=ctx.obj['REF'])
		Bio.AlignIO.write(alignments=[out1], handle=fout1, format='fasta')
		cova.utils.writecsv(fl=fout2, data=out2, sep='\t')
	
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
		out = [ [k,v] for k,v in cova.genome_seqtype(fin,fst=typefile).items()]
		out = sorted(out, key=lambda x: x[1])
		cova.utils.writecsv(fl=fout,data=out)		

	print("%s:\t SEQTYPE is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='genome_aln_unq.fna',show_default=True)
@click.option('--outfp',default='point_mutations.tsv',show_default=True)
@click.option('--outfd',default='deletions.tsv',show_default=True)
def vcalpd(ctx,infile,outfp,outfd):
	"""Call point mutations/deletions from Reference-limited MSA."""
	fin = os.path.join(ctx.obj['DR'], infile)
	fout1 = os.path.join(ctx.obj['DR'],outfp)
	fout2 = os.path.join(ctx.obj['DR'],outfd)

	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
	
	msa = cova.utils.MSA(fname=fin)
	
	# point mutations	
	if cova.utils.outcheck(fout1):
		print("{}: Calling point mutations from {}\n".format(cova.utils.timer(start),fin))
		tab1,head1 = msa.pointmuts(ref=ctx.obj['REF'],header=True)
		cova.utils.writecsv(fl=fout1, data=tab1, sep='\t', header=head1)
	
	# deletions
	if cova.utils.outcheck(fout2):
		print("{}: Calling deletions from {}\n".format(cova.utils.timer(start),fin))
		tab2,head2 = msa.dels(ref=ctx.obj['REF'],header=True)
		cova.utils.writecsv(fl=fout2, data=tab2, sep='\t', header=head2)
	
	print("%s:\t VCALPD is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='point_mutations.tsv',show_default=True)
@click.option('--outfile',default='prot_point_mutations_ann.tsv',show_default=True)
def annpv(ctx,infile,outfile):
	"""Annotate point mutations located within protein regions."""
	fin = os.path.join(ctx.obj['DR'], infile)
	fout = os.path.join(ctx.obj['DR'],outfile)

	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
	
	if cova.utils.outcheck(fout):
		print("%s: Annotating point mutations located within protein regions"%cova.utils.timer(start))
		cova.annotate_var(fin,fout)
	
	print("%s:\t ANNPV is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='prot_point_mutations_ann.tsv',show_default=True)
@click.option('--outfile',default='genome_vars.tsv',show_default=True)
def nsvar( ctx, infile, outfile):
	"""Get shared and unique non-synonymous variants for genomes."""
	fin = os.path.join(ctx.obj['DR'], infile)
	fout = os.path.join(ctx.obj['DR'], outfile)
	
	if cova.utils.outcheck(fout):
		print("%s: Identifying shared and unique variants"%cova.utils.timer(start))
		out = [ [k]+v for k,v in cova.genome_var(fin).items()]
		header = ['id','#var','#shared','#unique','shared','unique']
		out = sorted(out, key=lambda x: x[1])
		cova.utils.writecsv(fl=fout,data=out,sep='\t',header=header)		

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
		print("%s: Identifying and removing genomes with non-sense mutations."%cova.utils.timer(start))
		cova.rm_genome_w_stopm(fin,fout,finmsa,foutmsa)		

	print("%s:\t RMSTOP is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='genome_aln_sf.fna',show_default=True)
@click.option('--outdr',default='prots_nmsa',show_default=True)
def msap(ctx,infile,outdr):
	"""Extract nucleotide MSA of proteins from Reference-limited MSA."""
	fin = os.path.join(ctx.obj['DR'], infile)
	dout = os.path.join(ctx.obj['DR'], outdr)

	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)

	if cova.utils.outcheck(dout):
		msa = Bio.AlignIO.read(handle=fin, format='fasta')
		print('''%s: Extracting nucleotide MSA of proteins from Reference-limited MSA,
			Outputs will be saved to %s\n.'''%(cova.utils.timer(start),dout))
		cova.extract_nucmsa_prots(msa=msa, outdr=dout)
	
	print("%s:\t MSAP is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='genome_aln.fna',show_default=True)
@click.option('--outfile',default='insertions.tsv',show_default=True)
def vcali(ctx,infile,outfile):
	"""Call insertions from MSA."""
	fin = os.path.join(ctx.obj['DR'], infile)
	fout = os.path.join(ctx.obj['DR'],outfile)

	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
	
	if cova.utils.outcheck(fout):
		msa = cova.utils.MSA(fname=fin)
		print("%s: Calling insertions from MSA , Output will be saved to %s\n"%(cova.utils.timer(start),fout))
		tab, head = msa.ins( ref=ctx.obj['REF'], header=True)
		cova.utils.writecsv(fl=fout, data=tab, sep='\t', header=head)

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

	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
	
	if cova.utils.outcheck(fout1):
		print('''%s: Computing diversity for whole-genome and peptide-encoding regions.
		Output will be saved to %s\n'''%(cova.utils.timer(start),fout1))
		msa = cova.utils.MSA(fname=fin)
		wndiv = msa.ndiv()
		fpmsas = [ i for i in os.listdir(din) if i.endswith('.msa')]
		pndivs = [ [ i.replace('.msa',''), cova.utils.MSA(os.path.join(din,i)).ndiv()] for i in fpmsas]
		pndivs = [ i for i in pndivs if i[1] is not None]
		pndivs = sorted(pndivs, key=lambda x: x[1], reverse=True)
		out1 = [ ['genome', wndiv] ] + pndivs
		cova.utils.writecsv(fl=fout1, data=out1)
	
	if slide:
		if cova.utils.outcheck(fout2):
			print('''%s: Computing diversity within a window sliding over the genome.
			Output will be saved to %s\n'''%(cova.utils.timer(start),fout2))
			out2 = msa.slide_ndiv(window=window,jump=jump,ncpu=ncpu)
			cova.utils.writecsv(fl=fout2, data=out2)
	
	print("%s:\t DIV is done."%cova.utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--prog', default='FastTree', 
	help='''full path to FASTTREE program''', show_default=True)
@click.option('--infile', default='genome_aln_sf.fna', show_default=True)
@click.option('--outfile', default='tree.nwk', show_default=True)
def tree(ctx,prog,infile,outfile):
	"""Build phyogeny from whole-genome MSA."""
	fin = os.path.join(ctx.obj['DR'],infile)
	fout = os.path.join(ctx.obj['DR'],outfile)

	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
	
	# building
	if cova.utils.outcheck(fout):
		my_env = os.environ
		my_env['OMP_NUM_THREADS'] = ctx.obj['NCPU']
		# set path variable to find fasttree
		if 'COVA_BIN_PATH' in my_env.keys():
			my_env['PATH'] = ':'.join([ my_env['COVA_BIN_PATH'], my_env['PATH']])
		cmd = [prog, '-quiet', '-nt', '-mlnni', '4', '-boot', '100', fin]
		print("%s: Building Phylogeny from %s,\n Command: %s,\n Output will be saved to %s"%(\
			cova.utils.timer(start),fin,' '.join(cmd),fout))
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
def sel(ctx,prog,tree,indr,outdr,outr,outs):
	"""Identify sites under positive selection."""
	ftree = os.path.join(ctx.obj['DR'], tree)
	din   = os.path.join(ctx.obj['DR'], indr)
	dout  = os.path.join(ctx.obj['DR'],outdr)
	frout = os.path.join(ctx.obj['DR'], outr)
	fsout = os.path.join(ctx.obj['DR'], outs)
	
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
			
			if i.endswith('.msa'):
				cova.run_fubar(fmsa=os.path.join(din,i), ftree=ftree, outdr=dout, prog=prog)
	
	if cova.utils.outcheck(frout):		
		print('''%s: Parsing FUBAR output to generate rates and sites tables'''%cova.utils.timer(start))
		cova.parse_fubar(indr=dout, frout=frout, fsout=fsout)
	
	print("%s:\t SEL is done."%cova.utils.timer(start))

### command to run all other commands
@cli.command()
@click.pass_context
def full(ctx):
	"""
	Run full pipeline.
	Same as running: CoVa msabuild msaref msaunq seqtype vcalpd annpv nsvar rmstop msap vcali div tree sel
	"""
	# if new sequences are to be added
	if ctx.obj['ADDSEQ']:
		# add new sequences
		ctx.forward(msad)
		# clean up
		dr = ctx.obj['DR']
		fls = ['genome_aln_ref.fna', 'genome_aln_unq.fna','genome_aln_sf.fna','prots_nmsa',\
		'genome_w_stopm.tsv','genome_dups.tsv','genome_types.csv','genome_vars.tsv',\
		'point_mutations.tsv','deletions.tsv','prot_point_mutations_ann.tsv','insertions.tsv',\
		'divs.tsv','tree.nwk','fubar','rates.csv','sites.csv']
		print("%s:Deleting previously generated files"%cova.utils.timer(start))
		
		for f in fls:
			# full path of the file
			fp = os.path.join( dr, f)
			
			# if path already exists
			if os.path.exists(fp):
				
				# and if its a directory, remove the directory
				if os.path.isdir(fp):
					shutil.rmtree(fp)
				# else, if it's a file, delete the file
				else:
					os.remove(fp)
		
		print('\t Done with cleanup.')
	else:
		ctx.forward(msabuild)
	
	ctx.forward(msaref)
	ctx.forward(msaunq)
	ctx.forward(seqtype)
	ctx.forward(vcalpd)
	ctx.forward(annpv)
	ctx.forward(nsvar)
	ctx.forward(rmstop)
	ctx.forward(msap)
	ctx.forward(vcali)
	ctx.forward(div)
	ctx.forward(tree)
	ctx.forward(sel)

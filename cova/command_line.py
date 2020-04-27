import cova, os, sys, shutil, click, subprocess, time
from cova import _utils
from Bio import AlignIO

# start clock
start = time.time()

### commands ###
@click.group(chain=True)
@click.version_option()
@click.option('--indr', help='Full path to the working directory', default=os.getcwd(), type=click.Path())
@click.option('--ref', help='Reference accession', default='NC_045512', show_default=True)
@click.option('--ncpu', help='number of CPUs to use', default=4, show_default=True, type=int)
@click.option('--debug', help='See full traceback for errors.', is_flag=True)
@click.option('--addseq', help='Add new sequences and redo analysis.', is_flag=True)
@click.pass_context
def cli(ctx,indr,ref,ncpu,debug,addseq):
	"""
	Variant analysis using whole-genome Multiple Sequence Alignments.
	By default, it works on the current directory. Directory can be specified using INDR option.
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
def msabuild(ctx,prog,infile,outfile):
	"""Build whole-genome MSA.""" 
	fin = os.path.join(ctx.obj['DR'],infile)
	fout = os.path.join(ctx.obj['DR'],outfile)
	
	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)

	action = True
	if os.path.exists(fout):
		response=input('''output file %s already exists. 
Do you want to proceed and rewrite it? [y/n]\n'''%fout)
		if response == 'n':
			action = False
		elif response == 'y':
			action = True
		else:
			raise ValueError('''Inappropriate input! Please respond with [y] to proceed
and overwrite OR with [n] to terminate this command''')
	
	if action:
		# set path variable to find mafft
		my_env = os.environ
		if 'COVA_BIN_PATH' in my_env.keys():
			my_env['PATH'] = ':'.join([ my_env['COVA_BIN_PATH'], my_env['PATH']])

		cmd = [prog, '--quiet', '--nomemsave', '--maxiterate', '5', '--thread', ctx.obj['NCPU'], fin]
		print("%s: Building MSA from %s,\n Command: %s,\n Output will be saved to %s"%(\
			_utils.timer(start),fin,' '.join(cmd),fout))
		with open( fout,'w') as flob:
			s1 = subprocess.run( cmd, stdout=flob, env=my_env)
	else:
		print("okay! Existing output retained.")
	print("%s:\tMSABUILD is done."%_utils.timer(start))

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
	# then run the addseq command
	cmd = [prog, '--quiet', '--nomemsave', '--maxiterate', '5', '--thread', ctx.obj['NCPU'],\
	 	'--add', fin2, fcopy]
	print("%s: Adding sequences from %s to %s,\n Command: %s,\n Output will be saved to %s"%(\
		_utils.timer(start),fin2, fcopy,' '.join(cmd),fout))
	# rewrite MSA file
	with open( fout,'w') as flob:
		s1 = subprocess.run( cmd, stdout=flob, env=my_env)
	print("%s:\tMSAD is done."%_utils.timer(start))

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
	
	if action:
		msa = _utils.MSA(fname=fin)
		print("%s: Generating reference limited MSA, Output will be saved to %s\n"%(_utils.timer(start),fout))
		out = msa.limref(ref=ctx.obj['REF'])
		AlignIO.write(alignments=[out], handle=fout, format='fasta')
	else:
		print("okay! Existing output retained.")
	print("%s:\tMSAREF is done."%_utils.timer(start))

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
	
	action = True
	if os.path.exists(fout1):
		response=input('''output file %s already exists. 
Do you want to proceed and rewrite it? [y/n]\n'''%fout1)
		if response == 'n':
			action = False
		elif response == 'y':
			action = True
		else:
			raise ValueError('''Inappropriate input! Please respond with [y] to
proceed and overwrite OR with [n] to terminate this command''')
	
	if action:
		msa = _utils.MSA(fname=fin)
		print("%s: Removing duplicate sequences from the MSA\n"%_utils.timer(start))
		out1, out2 = msa.rmdup(ref=ctx.obj['REF'])
		AlignIO.write(alignments=[out1], handle=fout1, format='fasta')
		_utils.writecsv(fl=fout2, data=out2, sep='\t')
	else:
		print("okay! Existing output retained.")
	print("%s:\tMSAUNQ is done."%_utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='genome_aln_ref.fna',show_default=True)
@click.option('--outdr',default='prots_nmsa',show_default=True)
def msap(ctx,infile,outdr):
	"""Extract nucleotide MSA of proteins from Reference-limited MSA."""
	fin = os.path.join(ctx.obj['DR'], infile)
	dout = os.path.join(ctx.obj['DR'], outdr)

	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)

	msa = AlignIO.read(handle=fin, format='fasta')
	print('''%s: Extracting nucleotide MSA of proteins from Reference-limited MSA,
		Outputs will be saved to %s\n.'''%(_utils.timer(start),dout))
	cova.extract_nucmsa_prots(msa=msa, outdr=dout)
	print("%s:\t MSAP is done."%_utils.timer(start))

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
	
	msa = _utils.MSA(fname=fin)
	
	# point mutations
	action = True
	if os.path.exists(fout1):
		response=input('''output file %s already exists. 
Do you want to proceed and rewrite it? [y/n]\n'''%fout1)
		if response == 'n':
			action = False
		elif response == 'y':
			action = True
		else:
			raise ValueError('''Inappropriate input! Please respond with [y] to
proceed and overwrite OR with [n] to terminate this command''')
	
	if action:
		print("%s: Calling point mutations using reference-limited MSA\n"%_utils.timer(start))
		tab1,head1 = msa.pointmuts(ref=ctx.obj['REF'],header=True)
		_utils.writecsv(fl=fout1, data=tab1, sep='\t', header=head1)
	else:
		print("okay! Existing output retained.")
	
	# deletions
	action = True
	if os.path.exists(fout2):
		response=input('''output file %s already exists. 
Do you want to proceed and rewrite it? [y/n]\n'''%fout2)
		if response == 'n':
			action = False
		elif response == 'y':
			action = True
		else:
			raise ValueError('''Inappropriate input! Please respond with [y] to
proceed and overwrite OR with [n] to terminate this command''')

	if action:
		print("%s: Calling deletions using reference-limited MSA\n"%_utils.timer(start))
		tab2,head2 = msa.dels(ref=ctx.obj['REF'],header=True)
		_utils.writecsv(fl=fout2, data=tab2, sep='\t', header=head2)
	else:
		print("okay! Existing output retained.")
	print("%s:\t VCALPD is done."%_utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='point_mutations.tsv',show_default=True)
@click.option('--outfile',default='prot_point_mutations_annotated.tsv',show_default=True)
def annpv(ctx,infile,outfile):
	"""Annotate point mutations located within protein regions."""
	fin = os.path.join(ctx.obj['DR'], infile)
	fout = os.path.join(ctx.obj['DR'],outfile)

	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
	
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
	
	if action:
		print("%s: Annotating point mutations located within protein regions"%_utils.timer(start))
		cova.annotate_var(fin,fout)
	else:
		print("okay! Existing output retained.")
	print("%s:\t ANNPV is done."%_utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='prot_point_mutations_annotated.tsv',show_default=True)
@click.option('--outfile',default='genome_variants.tsv',show_default=True)
def tabvs(ctx,infile,outfile):
	"""Tabulate genomes with their shared and unique variants."""
	fin = os.path.join(ctx.obj['DR'], infile)
	fout = os.path.join(ctx.obj['DR'], outfile)
		
	print("%s: Identifying shared and unique non-synonymous variants"%_utils.timer(start))
	cova.genome_var(fin,fout)
	print("%s:\t TABVS is done."%_utils.timer(start))

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
	
	if action:
		msa = _utils.MSA(fname=fin)
		print("%s: Calling insertions from MSA , Output will be saved to %s\n"%(_utils.timer(start),fout))
		tab, head = msa.ins( ref=ctx.obj['REF'], header=True)
		_utils.writecsv(fl=fout, data=tab, sep='\t', header=head)
	else:
		print("okay! Existing output retained.")
	print("%s:\t VCALI is done."%_utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--infile',default='genome_aln_ref.fna',show_default=True)
@click.option('--indr',default='prots_nmsa',show_default=True)
@click.option('--outfile',default='divs.csv',show_default=True)
def div(ctx,infile,indr,outfile):
	"""Compute nucleotide diversity from Reference-limited MSA."""
	fin = os.path.join(ctx.obj['DR'], infile)
	din = os.path.join(ctx.obj['DR'], indr)
	fout = os.path.join(ctx.obj['DR'],outfile)
	ncpu = int(ctx.obj['NCPU'])
	
	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
	
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
	
	if action:
		print("%s: Computing diversity from MSA, Output will be saved to %s\n"%(_utils.timer(start),fout))
		msa = _utils.MSA(fname=fin)
		wndiv = msa.ndiv(ncpu)
		fpmsas = [ i for i in os.listdir(din) if i.endswith('.msa')]
		pndivs = [ [ i.replace('.msa',''), _utils.MSA(os.path.join(din,i)).ndiv(ncpu)] for i in fpmsas]
		pndivs = [ i for i in pndivs if i[1] is not None]
		pndivs = sorted(pndivs, key=lambda x: x[1], reverse=True)
		out = [ ['genome', wndiv] ] + pndivs
		_utils.writecsv(fl=fout, data=out)
	else:
		print("okay! Existing output retained.")
	print("%s:\t DIV is done."%_utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--prog', default='FastTree', 
	help='''full path to FASTTREE program''', show_default=True)
@click.option('--infile', default='genome_aln_ref.fna', show_default=True)
@click.option('--outfile', default='genomes.nwk', show_default=True)
@click.option('--plotfile', default='genomes_tree.png', show_default=True)
@click.option('--mapfile', default=None, show_default=True)
def tree(ctx,prog,infile,outfile,plotfile,mapfile):
	"""Build phyogeny from whole-genome MSA."""
	fin = os.path.join(ctx.obj['DR'],infile)
	fout = os.path.join(ctx.obj['DR'],outfile)
	fplot = os.path.join(ctx.obj['DR'],plotfile)
	if mapfile is not None:
		fmap = os.path.join(ctx.obj['DR'],mapfile)
	else:
		fmap = None

	if not os.path.exists(fin):
		raise FileNotFoundError("couldn't read the input file %s."%fin)
	
	# building
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
	
	if action:
		my_env = os.environ
		my_env['OMP_NUM_THREADS'] = ctx.obj['NCPU']
		# set path variable to find fasttree
		if 'COVA_BIN_PATH' in my_env.keys():
			my_env['PATH'] = ':'.join([ my_env['COVA_BIN_PATH'], my_env['PATH']])
		cmd = [prog, '-quiet', '-nt', '-mlnni', '4', '-nosupport', fin]
		print("%s: Building Phylogeny from %s,\n Command: %s,\n Output will be saved to %s"%(\
			_utils.timer(start),fin,' '.join(cmd),fout))
		with open( fout,'w') as flob:
			s1 = subprocess.run( cmd, stdout=flob, stderr=subprocess.DEVNULL, env=my_env)
	else:
		print("okay! Existing output retained.")
		
	# plotting	
	action = True
	if os.path.exists(fplot):
		response=input('''output file %s already exists. 
Do you want to proceed and rewrite it? [y/n]\n'''%fplot)
		if response == 'n':
			action = False
		elif response == 'y':
			action = True
		else:
			raise ValueError('''Inappropriate input! Please respond with [y] to
proceed and overwrite OR with [n] to terminate this command''')
	
	if action:
		print("%s: Plotting tree from %s"%(_utils.timer(start),fout))
		cova.plottree(ftree=fout, fmap=fmap, fplot=fplot)	
	else:
		print("okay! Existing output retained.")
	print("%s:\t TREE is done."%_utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--prog', default='hyphy',
	help='''full path to HYPHY program''', show_default=True)
@click.option('--indr', default='prots_nmsa', show_default=True)
@click.option('--tree', default='genomes.nwk', show_default=True)
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
	if not os.path.exists(dout):
		os.mkdir(dout)
	
	# set path variable to find hyphy
	my_env = os.environ
	if 'COVA_BIN_PATH' in my_env.keys():
		my_env['PATH'] = ':'.join([ my_env['COVA_BIN_PATH'], my_env['PATH']])

	# check if hyphy can find its batch files
	cmd = [prog, 'fubar', '--help']
	s = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, env=my_env)
	if s.returncode != 0:
		raise FileNotFoundError("Hyphy couldn't find where FUBAR is. Check library path")

	print('''{}: Analysing MSAs for positive selection using tree ({})
		Outputs will be saved to {}'''.format(_utils.timer(start),ftree,dout))
	for i in os.listdir(din):
		if i.endswith('.msa'):
			cova.run_fubar(fmsa=os.path.join(din,i), ftree=ftree, outdr=dout, prog=prog)
			
	print('''%s: Parsing FUBAR output to generate rates and sites tables'''%_utils.timer(start))
	cova.parse_fubar(indr=dout, frout=frout, fsout=fsout)
	print("%s:\t SEL is done."%_utils.timer(start))

### command to run all other commands
@cli.command()
@click.pass_context
def full(ctx):
	"""Run full pipeline."""
	ctx.forward(msabuild)
	
	# if new sequences are to be added
	if ctx.obj['ADDSEQ']:
		# add new sequences
		ctx.forward(msad)
		# clean up
		dr = ctx.obj['DR']
		fls = ['genome_aln_ref.fna', 'genome_aln_unq.fna','genome_dups.tsv','prots_nmsa',\
		'point_mutations.tsv','deletions.tsv','prot_point_mutations_annotated.tsv','insertions.tsv',\
		'divs.tsv','genomes.nwk','genomes_tree.png','fubar','rates.csv','sites.csv',\
		'genome_variants.tsv']
		print("%s:Deleting pre-existing analysis files"%_utils.timer(start))
		for f in fls:
			fp = os.path.join( dr, f)
			if os.path.exists(fp):
				if os.path.isdir(fp):
					shutil.rmtree(fp)
				else:
					os.remove(fp)
		print('\t Done with cleanup.')
	ctx.forward(msaref)
	ctx.forward(msaunq)
	ctx.forward(msap)
	ctx.forward(vcalpd)
	ctx.forward(annpv)
	ctx.forward(tabvs)
	ctx.forward(vcali)
	ctx.forward(div)
	ctx.forward(tree)
	ctx.forward(sel)
import cova, os, sys, click, subprocess, time
from cova import _utils
from Bio import AlignIO

# start clock
start = time.time()

### commands ###
@click.group(chain=True)
@click.version_option()
@click.option('--debug', is_flag=True)
@click.option('--indr', default=os.getcwd(), type=click.Path())
@click.option('--ref', default='NC_045512', show_default=True)
@click.option('--ncpu', default='4', show_default=True, type=str)
@click.pass_context
def cli(ctx,indr,ref,ncpu,debug):
	"""
	Variant analysis using whole-genome Multiple Sequence Alignments.
	By default, it works on the current directory. Directory can be specified using INDR option.
	"""
	ctx.ensure_object(dict)
	ctx.obj['DR'] = indr
	ctx.obj['REF'] = ref
	ctx.obj['NCPU'] = ncpu
	# control traceback
	if debug:
		ctx.obj['DEBUG'] = debug
		click.echo("Debug mode is ON.\n")
	else:
		sys.excepthook = lambda exctype,exc,traceback : print("{}: {}".format(exctype.__name__,exc)) 
	click.echo("CoVa will run in the directory: {}\n".format(indr))

@cli.command()
@click.pass_context
@click.option('--prog', default=os.path.join(cova.PROGPATH,'mafft'), 
	help='''full path to MAFFT program''',show_default=True)
@click.option('--infile', default='genomes.fna', show_default=True)
@click.option('--outfile', default='genome_aln.fna', show_default=True)
def msabuild(ctx,prog,infile,outfile):
	"""Build whole-genome Multiple Sequence Alignment.""" 
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
		cmd = [prog, '--quiet', '--nomemsave', '--maxiterate', '5', '--thread', ctx.obj['NCPU'], fin]
		print("%s: Building MSA from %s,\n Command: %s,\n Output will be saved to %s"%(\
			_utils.timer(start),fin,' '.join(cmd),fout))
		with open( fout,'w') as flob:
			s1 = subprocess.run( cmd, stdout=flob)
	else:
		print("okay! Existing output retained.")
	print("%s:\tMSABUILD is done."%_utils.timer(start))
			
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

### Remove identical sequences
@cli.command()
@click.pass_context
@click.option('--infile', default='genome_aln_ref.fna', show_default=True)
@click.option('--outfile',default='genome_aln_unq.fna', show_default=True)
def msaunq( ctx, infile, outfile):
	"""Remove duplicate sequences from reference-limited MSA."""
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
		print("%s: Removing duplicate sequences from the MSA\n"%_utils.timer(start))
		out = msa.rmdup(ref=ctx.obj['REF'])
		AlignIO.write(alignments=[out], handle=fout, format='fasta')
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
	"""Call variants ( point mutations/deletions ) from Reference-limited MSA."""
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
@click.option('--infile',default='genome_aln.fna',show_default=True)
@click.option('--outfile',default='insertions.tsv',show_default=True)
def vcali(ctx,infile,outfile):
	"""Call variants ( insertions ) from Reference-limited MSA."""
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
@click.option('--outfile',default='divs.tsv',show_default=True)
def div(ctx,infile,indr,outfile):
	"""Compute whole-genome and gene-wise diversity from Reference-limited MSA."""
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
		_utils.writecsv(fl=fout, data=out, sep='\t')
	else:
		print("okay! Existing output retained.")
	print("%s:\t DIV is done."%_utils.timer(start))

@cli.command()
@click.pass_context
@click.option('--prog', default=os.path.join(cova.PROGPATH,'FastTree'), 
	help='''full path to FASTTREE program''', show_default=True)
@click.option('--infile', default='genome_aln_ref.fna', show_default=True)
@click.option('--outfile', default='genomes.nwk', show_default=True)
@click.option('--plotfile', default='genomes_tree.png', show_default=True)
@click.option('--mapfile', default=None, show_default=True)
def tree(ctx,prog,infile,outfile,plotfile,mapfile):
	"""Build phyogeny from whole-genome Multiple Sequence Alignment."""
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
@click.option('--prog', default=os.path.join(cova.PROGPATH,'hyphy'),
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
	
	# check if hyphy can find its batch files
	cmd = [prog, 'fubar', 'LIBPATH='+cova.LIBPATH, '--help']
	s = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
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

@cli.command()
@click.pass_context
@click.option('--infile1',default='point_mutations.tsv',show_default=True)
@click.option('--infile2',default='prot_point_mutations_annotated.tsv',show_default=True)
@click.option('--outfile',default='genome_variants.tsv',show_default=True)
def tabvs(ctx,infile1,infile2,outfile):
	"""Return table of genomes with their shared and unique non-synonymous changes."""
	fin1 = os.path.join(ctx.obj['DR'], infile1)
	fin2 = os.path.join(ctx.obj['DR'], infile2)
	fout = os.path.join(ctx.obj['DR'], outfile)
		
	print("%s: Identifying shared and unique non-synonymous variants"%_utils.timer(start))
	cova.genome_var(fpm=fin1,fann=fin2,fout=fout)
	print("%s:\t TABVS is done."%_utils.timer(start))

### command to run all other commands
@cli.command()
@click.pass_context
def full(ctx):
	"""Run full pipeline."""
	ctx.forward(msabuild)
	ctx.forward(msaref)
	ctx.forward(msaunq)
	ctx.forward(msap)
	ctx.forward(vcalpd)
	ctx.forward(annpv)
	ctx.forward(vcali)
	ctx.forward(div)
	ctx.forward(tree)
	ctx.forward(sel)
	ctx.forward(tabvs)
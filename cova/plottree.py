import os, sys, click, ete3, seaborn, matplotlib, pandas
from time import time
from . import utils

start = time()

## functions
def typesortingfunc(t):
	if 'ST' not in t:
		return 1000
	else:
		return int(t.replace('ST',''))

def colbyinfo(infodict,sorting_func = None):
	"""Return a dict of colors for labels given a dict of some info on labels."""
	# a list of unique info labels
	if sorting_func is None:
		infotags = sorted( list( set( infodict.values())))
	else:
		infotags = sorted( list( set( infodict.values())), key = sorting_func)
	# no. of such tags
	ntags = len(infotags)
	# a list of these many colors
	colors = utils.get_N_HexCol(N=ntags)
	# dict of infotags and their corresponding color
	infotag_color = { infotags[x]:colors[x] for x in range(ntags)}
	# dict of data labels and their corresponding color
	datalabel_color = { data:infotag_color[info] for data,info in infodict.items()}
	
	return (datalabel_color,infotag_color)

@click.command()
@click.option('--dr', help='full path to input directory', type=click.Path(), default=os.getcwd(), show_default=True)
@click.option('--ftree', help='tree file name', default='tree.nwk', show_default=True)
@click.option('--fplot', help='output plot file, default: <ftree>.png')
@click.option('--fst', help='file with info on sequence types', default='genome_types.csv', show_default=True)
@click.option('--fld', help='file with info on location and date')
@click.option('--typef', help='minimum freq of a type to include in the figure', type=float, default=0.01, show_default=True)
@click.option('--branch_scale', help='Scale used to draw branch lengths', type=int, default=500, show_default=True)
@click.option('--show', help='Show tree instead of saving?', is_flag=True)
@click.option('--branch_support', help='Show branch support?', is_flag=True)
@click.option('--show_legend', help='Include legend?', is_flag=True)
@click.option('--legend_box_size', help='size of the squares drawn in the legend', type=int, default=100, show_default=True)
@click.option('--max_legend_stack', help='maximum size of the squares stack in the legend', type=int, default=500, show_default=True)
@click.option('--legend_font_size', help=' font size for legend text', type=int, default=40, show_default=True)
@click.option('--img_height', help='output image height', type=int, default=2500, show_default=True)
@click.option('--img_dpi', help='dots per inch for output image', type=int, default=600, show_default=True)
@click.option('--typecoldict', help='dict of types and correponding colors', type=str,required=False)

def main_fun(dr,ftree,fplot,fst,fld,typef,branch_scale,branch_support,show_legend,legend_box_size,
	max_legend_stack,legend_font_size,img_height,img_dpi,show,typecoldict):
	## paths
	ftree = os.path.join( dr, ftree)
	if fplot is None:
		fplot = ftree.replace('.nwk','.png')
	fst = os.path.join( dr, fst)
	
	## checks
	# tree file is present
	if not os.path.exists(ftree):
		raise FileNotFoundError('tree file %s must be present.'%ftree)
	
	# should you proceed if the output path already exists
	if not utils.outcheck(fplot):
		return

	# plot file has png suffix
	if fplot.split('.')[-1] != 'png':
		raise ValueError('output file must have suffix "png".')
	
	# if info file ( for sequence types) is provided, it is a valid file path
	if fld:
		click.echo("Location file provided.")
		fld = os.path.join( dr, fld)
		
		if not os.path.exists(fld):
	 		raise FileNotFoundError("couldn't find the file %s."%fld)
	
	else:
		click.echo("No location file! Annotation will only be for sequence types")

	# load tree
	t = ete3.Tree(ftree)
	# list of leaves
	leaves = t.get_leaf_names()

	## create treestyle
	ts = ete3.TreeStyle()
	ts.show_branch_support = branch_support
	ts.mode = "c"
	ts.scale = branch_scale

	### types #####################################
	# table of genomes and their sequence types
	typedata = utils.readcsv(fst)
	# threshold for a type to be shown explicitly in the figure
	th = len(typedata) * typef
	# dict of sequence type with isolates
	type_isols = utils.split_data(data=typedata, ix=1, cixs=0)
	# empty list of types to be removed
	rmkeys = []
	# empty list of such minor isolates
	minors = [] 

	# for every type and its isolates
	for k,v in type_isols.items():
		
		# if the type is unkown
		if k == 'U':
			# skip
			continue
		
		# if no. of isolates for the types are less than the above threshold
		if len(v) < th:
			# minor isolates
			minors.extend(v)
			# excluded type
			rmkeys.append(k)

	# type isolate dict with low represetation types excluded
	type_isols = {k:v for k,v in type_isols.items() if k not in rmkeys}
	# and added back as minors under type 'O'thers
	type_isols['O'] = minors
	# modified table of genome and types
	typedata = [ [i,k] for k,v in type_isols.items() for i in v] 
	# dict of isolate and its type if the isolate is present on the tree
	isol_type = { i[0]:i[1] for i in typedata if i[0] in leaves} 	
	# color representation of types
	isol_type_color, type_color = colbyinfo(infodict=isol_type,sorting_func=typesortingfunc)
	
	# if a color dict was explicitly provided
	if typecoldict is not None:
		tcl = typecoldict.split(',')
		type_color = { tcl[x]:tcl[x+1] for x in range(0,len(tcl),2)}
		isol_type_color = { k:type_color[v] for k,v in isol_type.items() }

	for k,v in isol_type_color.items():
		
		if v == type_color['U']:
			isol_type_color[k] = 'white'

		if 'O' in type_color.keys() and v == type_color['O']:
			isol_type_color[k] = 'grey'

	type_color['U'] = 'white'
	type_color['O'] = 'grey'
	###############################################

	# basic tree style with type annotation
	for n in t.traverse():
		
		# if branch support is less than 0.5, delete the branch
		if n.support < 0.5:
			n.delete()
			continue
		
		n.dist = 0.1
		ns = ete3.NodeStyle()
		if n.is_leaf():	
			ns['size'] = 10
			if n.name in isol_type_color.keys():
				ns['bgcolor'] = isol_type_color[n.name]
			else:
				ns['bgcolor'] = 'grey'
		else:
			ns['size'] = 0
		n.set_style(ns)

	# If mapping is available, then use it to color leaves and branches
	if fld is not None: 
		dmap = pandas.read_csv(fld)
		nrow = len(dmap)
		head = dmap.columns

		# colors for locations
		isol_loc = { dmap.at[x,'accession']:dmap.at[x,'location'] for x in range(nrow)}
		isol_loc_color, loc_color = colbyinfo(infodict=isol_loc)

		# colors for months
		isol_month = { dmap.at[x,'accession']:'-'.join(dmap.at[x,'date'].split('-')[:2]) for x in range(nrow) if dmap.at[x,'date'].count('-') == 2}
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
		
		boxsize = 10*branch_scale/100
		
		for n in t.traverse():
			
			if n.name not in isol_loc_color.keys():
				continue

			if n.is_leaf():	
				rct1 = ete3.RectFace(width=boxsize,height=boxsize,fgcolor='',bgcolor=isol_loc_color[n.name])
				n.add_face( rct1, column=2, position='aligned')
				if n.name in isol_month_color.keys():
					rct2 = ete3.RectFace(width=boxsize,height=boxsize,fgcolor='',bgcolor=isol_month_color[n.name])
					n.add_face( rct2, column=3, position='aligned')
			else:
				n.img_style['size'] = 0
			
		
	### legend ##################################
	if show_legend:
		ts.legend_position = 3
		stack_size = 0
		colx = 0

		for k,v in type_color.items():

			rct = ete3.RectFace(legend_box_size,legend_box_size,'',v)
			rct.margin_left = 10
			rct.margin_right = 10
			txt = ete3.TextFace( k, fsize=legend_font_size)
			txt.margin_left = 10
			txt.margin_right = 10
			
			if stack_size > max_legend_stack:
				stack_size = 0
				colx += 2

			if stack_size == 0:
				rct.margin_top = 20
			
			ts.legend.add_face(rct, column=colx)
			ts.legend.add_face(txt, column=colx+1)
			stack_size += legend_box_size
	###############################################

	## output
	if show is not None:
		t.render(fplot,tree_style=ts,units='px',h=img_height,dpi=img_dpi)
		click.echo("{}: Tree plotting complete. Output was saved in {}".format(utils.timer(start),fplot))
	else:
		t.show(tree_style=ts)
	################
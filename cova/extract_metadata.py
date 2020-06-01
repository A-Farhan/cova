## require
import os, pandas, numpy, click
from . import _utils, REF

## Functions ##
def extract_metadata_gisaid_table(fin,fout=None):
	
	if fout is None:
		fout = os.path.join( os.path.dirname(fin), 'genome_info.csv')
	
	xls = pandas.read_excel(fin, sheet_index = 0, skiprows=[0,1], 
		usecols=['Accession ID', 'Location', 'Collection date'])
	nrows = len(xls)
	ac_info = {}
	
	for x in range(1,nrows):
		row = xls.loc[x]
		ac = row[ 'Accession ID' ].replace('EPI_ISL_','')
		
		try:
			loc = row[ 'Location' ].split( ' / ')[1].strip(' ')
		except IndexError:
			print("Bad format location! {}".format(row['Location']))
			continue

		date = row[ 'Collection date' ]
		ac_info[ac] = [loc,date]
	
	# include reference, if not present alreadt
	if REF not in ac_info.keys():
		ac_info[REF] = ['China','2019-12']

	out = [ [k]+v for k,v in ac_info.items()]
	_utils.writecsv(fl=fout, data=out, header=['accession','location','date'])

def extract_metadata_ncbivirus_table(fin,fout=None):
	
	if fout is None:
		fout = os.path.join( os.path.dirname(fin), 'genome_info.csv')
	
	data = pandas.read_csv(fin,usecols=['Accession', 'Geo_Location', 'Collection_Date'])
	nrows = len(data)
	ac_info = {}

	for x in range(nrows):
		row = data.loc[x]
		ac, loc, date = row
		loc = loc.split(":")[0]
		ac_info[ac] = [loc,date]
		
	# include reference, if not present alreadt
	if REF not in ac_info.keys():
		ac_info[REF] = ['China','2019-12']

	out = [ [k]+v for k,v in ac_info.items()]
	_utils.writecsv(fl=fout, data=out, header=['accession','location','date'])

## Main ##
@click.command()
@click.option('--fin', help='metadata table, either from GISAID(xls) OR NCBI Virus(csv)', 
	type=click.Path(), required=True)
@click.option('--fout', help='output file name', default='genome_info.csv', 
	show_default=True, type=click.Path())
@click.option('--source', help='source of metadata', required=True, 
	show_default=True, type=click.Choice(['ncbivirus', 'gisaid']))
 
def main_fun(fin,fout,source):
	outpath = os.path.join( os.path.dirname(fin), fout)
	click.echo("Output will be saved to path: %s"%outpath)

	if not _utils.outcheck(outpath):
		return

	if source == 'gisaid':
		click.echo('Source is GISAID.')
		extract_metadata_gisaid_table(fin,fout=outpath)

	if source == 'ncbivirus':
		click.echo('Source is NCBI Virus.')
		extract_metadata_ncbivirus_table(fin,fout=outpath)
	print("metadata extraction complete!")
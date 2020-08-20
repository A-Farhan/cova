## require
import os, click
from time import time
from Bio import SeqIO
from . import utils

start = time()

@click.command()
@click.option('--fseq', help='full path to multi-FASTA GISAID file where headers are strain names',
	type=click.Path(), required=True)
@click.option('--fmd', help='full path to GISAID metadata table corresponding to the FASTA file',
	type=click.Path(), required=True)

def main_fun(fseq,fmd):
	"""Replace FASTA headers with corresponding gisaid EPI ISL accessions given metadata."""
	# output file path
	fout = os.path.join( os.path.dirname(fseq), 'sequences_w_gisaid_acc.fasta')
	click.echo("Output will be saved to path: %s"%fout)
	
	# check if we should proceed in case output is present
	if not utils.outcheck(fout):
		return
	
	# sequence records
	recs = SeqIO.parse(fseq,'fasta')
	# metadata
	md = pandas.read_csv(fmd,sep='\t')
	# number of entries with metadata
	nrow = len(md)
	# dict of strain name and gisaid accession
	strain_ac = {md.at[x,'strain']:md.at[x,'gisaid_epi_isl'] for x in range(nrow)}
	# empty list for output records
	outrecs = []

	# for every record
	for rec in recs:
		# replace its ID with gisaid accession
		rec.id = strain_ac[rec.id]
		# scrap its name and description
		rec.name = ''
		rec.description = ''
		# add to output list
		outrecs.append(rec)

	# write output to file
	SeqIO.write(outrecs,fout,'fasta')
	click.echo("%s: Done."%utils.timer(start))
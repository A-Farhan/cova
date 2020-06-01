## require
import os, click
from time import time
from Bio import SeqIO
from . import _utils, GENOME

start = time()

## Functions ##
def rm_lowqualseq(fin,p=1):
	"""
	Remove low quality entries from a FASTA file of nucleotide sequences.
	
	Arguments:
	- fin 	- FASTA file of nucleotide sequences
	- p 	- acceptable percentage of ambiguous characters

	Value: biopython sequence record
	"""
	recs = SeqIO.parse(handle=fin, format='fasta')
	modrecs = [ i for i in recs if _utils.is_seq_qual(rec=i, ambt=p)]
	return modrecs

def trim_header_gisaid(fin):
	"""Trim FASTA headers of sequences downloaded from GISAID."""
	recs = SeqIO.parse(handle=fin, format='fasta')
	outrecs = []
	
	for rec in recs:
		try:
			rec.id = rec.description.split('|')[1].replace('EPI_ISL_','')
		except IndexError:
			print("record: {}".format(rec))
			raise
		rec.description=''
		outrecs.append(rec)
	
	return outrecs

def add_ref(fin):
	"""Add reference sequence to a FASTA file of genomes."""
	recs = list(SeqIO.parse(handle=fin, format='fasta'))
	refrec = GENOME
	refrec.id = refrec.id.split('.')[0]
	refrec.name=''
	refrec.description=''
	ids = [ i.id for i in recs]
	
	if refrec.id not in ids:
		recs.insert(0, refrec)

	return recs

## Main ##
@click.command()
@click.option('--fin', help='FASTA file of genome sequences', type=click.Path(), required=True)
@click.option('--fout', help='output file name', default='genomes.fna', 
	show_default=True, type=click.Path())
@click.option('--tqual', help='percent threshold of ambiguous characters', default=1, 
	show_default=True,type=click.FloatRange(min=0, max=100))
@click.option('--gisaid', help='If sequences were downloaded from GISAID?', is_flag=True)

def main_fun(fin,fout,tqual,gisaid):
	"""Preprocess FASTA file of genomes."""
	# full path to output	
	outpath = os.path.join( os.path.dirname(fin), fout)
	click.echo("Output will be saved to path: %s"%outpath)

	# check if we should proceed in case output is present
	if not _utils.outcheck(outpath):
		return

	# remove low quality genomes
	rec1 = rm_lowqualseq(fin,p=tqual)
	SeqIO.write(sequences=rec1, handle='tmp', format='fasta')
	
	# trim gisaid header	
	if gisaid:
		click.echo('Source is GISAID.')
		rec2 = trim_header_gisaid(fin='tmp')
		SeqIO.write(sequences=rec2, handle='tmp', format='fasta')
	
	# add reference
	rec3 = add_ref(fin='tmp')
	SeqIO.write(sequences=rec3, handle=outpath, format='fasta')
	
	# clean up
	os.remove('tmp')
	click.echo("%s: Processed output was saved in %s."%(_utils.timer(start),outpath))
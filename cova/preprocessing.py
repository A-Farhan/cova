## require
import os, click
from shutil import copyfile
from time import time
from Bio import SeqIO
from . import utils, GENOME

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
	modrecs = [ i for i in recs if utils.is_seq_qual(rec=i, ambt=p)]
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
	ids = [ i.id for i in recs]
	refrec = GENOME
	refrec.id = refrec.id.split('.')[0]
	refrec.name=''
	refrec.description=''
	
	if refrec.id not in ids:
		recs.insert(0, refrec)
	else:
		print('\tReference id "%s" already present in the input.'%refrec.id)

	return recs

## Main ##
@click.command()
@click.option('--fin', help='FASTA file of genome sequences', type=click.Path(), required=True)
@click.option('--fout', help='output file name', default='genomes.fna', 
	show_default=True, type=click.Path())
@click.option('--ftmp', help='temporary intermediate file name', default='tmp.fna', 
	show_default=True, type=click.Path())
@click.option('--tqual', help='percent threshold of ambiguous characters', default=1, 
	show_default=True,type=click.FloatRange(min=0, max=100))
@click.option('--gisaid', help='If sequences were downloaded from GISAID?', is_flag=True)
@click.option('--noref', help='Do not include reference', is_flag=True)

def main_fun(fin,fout,ftmp,tqual,gisaid,noref):
	"""Preprocess FASTA file of genomes."""
	# full path to output	
	outpath = os.path.join( os.path.dirname(fin), fout)
	click.echo("Output will be saved to path: %s"%outpath)

	# check if we should proceed in case output is present
	if not utils.outcheck(outpath):
		return

	# full path to temporary file
	tmppath = os.path.join( os.path.dirname(fin), ftmp)
	click.echo("\tIntermediate file path: %s"%tmppath)
	# remove low quality genomes
	rec1 = rm_lowqualseq(fin,p=tqual)
	SeqIO.write(sequences=rec1, handle=tmppath, format='fasta')
	
	# trim gisaid header	
	if gisaid:
		click.echo('Source is GISAID.')
		rec2 = trim_header_gisaid(fin=tmppath)
		SeqIO.write(sequences=rec2, handle=tmppath, format='fasta')
	
	# add reference
	if noref:
		click.echo("Reference will not be added to the output.")
		copyfile(src=tmppath, dst=outpath)
	else:
		click.echo('Adding Reference to the output..')
		rec3 = add_ref(fin=tmppath)
		SeqIO.write(sequences=rec3, handle=outpath, format='fasta')
	
	# clean up
	os.remove(tmppath)
	click.echo("%s: Processed output was saved in %s."%(utils.timer(start),outpath))
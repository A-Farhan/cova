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
	modrecs = [ i for i in recs if utils.is_seq_qual(seq=i.seq, ambt=p)]
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

def fix_ids(reclist):
	"""Fix sequence IDs."""
	out = []
	
	# for every sequence record
	for i in reclist:
		
		# if it's from NCBI, then the accession would be separated with the rest by |
		if '|' in i.id:
			print("\tSplitting ID on |")
			i.id = i.id.split('|')[0]
		# spaces are removed from the descriptio, if it's gisaid and converted to ID
		elif ' ' in i.description:
			print("\tRemoving spaces from Description")
			i.id = i.description.replace(' ','')
		else:
			pass

		out.append(i)
	return out	

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
	click.echo("Removing low quality genomes at {} % threshold..".format(tqual))
	rec1 = rm_lowqualseq(fin,p=tqual)
	# fix ids: remove spaces
	click.echo("Fixing headers..")
	rec1 = fix_ids(rec1)
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
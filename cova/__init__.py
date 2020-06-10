"""Coronavirus Variant analysis"""
import pkg_resources, os
from Bio import SeqIO
from . import utils

DATAPATH = pkg_resources.resource_filename('cova', 'data')
FEATURETABLE = utils.readcsv(fl = os.path.join( DATAPATH, 'feature_table.txt'),sep='\t')
PROTNAMES = { i[0]:i[1] for i in utils.readcsv(fl=os.path.join( DATAPATH, 'protnames_map.csv'))}
GENOME = SeqIO.read(handle = os.path.join( DATAPATH, 'genome.fna'), format='fasta')
REF = GENOME.id.split('.')[0]
TYPEPOS, SEQTYPES = utils.readcsv(fl = os.path.join( DATAPATH, 'seq_types.csv'), header=True)
TYPEPOS = [ int(i) for i in TYPEPOS[1:]]
SEQTYPES = { ''.join(i[1:]):i[0] for i in SEQTYPES}

from ._func import rm_genome_w_stopm, extract_nucmsa_prots, annotate_var,\
run_fubar, parse_fubar, genome_var, genome_seqtype
from . import _divs

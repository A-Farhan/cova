"""Coronavirus Variant analysis"""
import pkg_resources, os

PROGPATH = pkg_resources.resource_filename('cova', 'softwares/bin')
LIBPATH =  pkg_resources.resource_filename('cova', 'softwares/lib/hyphy')
DATAPATH = pkg_resources.resource_filename('cova', 'data')

from . import _utils
FEATURETABLE = _utils.readcsv(fl = os.path.join( DATAPATH, 'feature_table.txt'),sep='\t')
PROTNAMES = { i[0]:i[1] for i in _utils.readcsv(fl=os.path.join( DATAPATH, 'protnames_map.csv'))}

from Bio import SeqIO
GENOME = SeqIO.read(handle = os.path.join( DATAPATH, 'genome.fna'), format='fasta')

from ._func import extract_nucmsa_prots, annotate_var, plottree, run_fubar, parse_fubar, genome_var

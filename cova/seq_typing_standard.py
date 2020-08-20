## require
import os, sys, pandas, datetime, click
from Bio import AlignIO
from cova import utils, TYPEPOS
from time import time

## Functions ###########
def barcode_genomes(aln,pos,m,alpha='ACGT'):
    """
    Return a dict of a barcoding sequence as key and a list of genomes with that sequence as its value.

    Arguments:
    - aln       - whole-genome reference-limited Multiple Sequence Alignment 
    - pos       - list of barcoding positions
    - m         - minimum no. of genomes to keep a type
    - alpha     - string of acceptable characters in a barcoding sequence
    """ 
    # list of column sequences at barcoding positions
    pseqs = [ aln[:,p-1] for p in pos]
    # list of tuples( genome id, row suquence limited to these positions)
    id_seq = [ ( aln[x].id, ''.join(i).upper()) for x,i in enumerate(zip(*pseqs))]
    # dict of sequence tags and their corresponding list of genome ids
    st_ids = utils.split_data(data=id_seq, ix=1, cixs=0)
    # discard types with non-acceptable characters
    # and with only m examples
    out = { k:v for k,v in st_ids.items() if len(v) >= m and not any(i not in alpha for i in k)}
    return out
########################################

## Main ##
@click.command()
@click.option('--dr', help='full path to working directory', type=click.Path(), required=True)
@click.option('--ming', help='minimum no. of genomes to consider a type', type=click.IntRange(min=1))
@click.option('--ftp', help='full path to the list of barcoding positions', type=click.Path())
@click.option('--fout', help='output file name', default='seqtype_def.csv', 
    show_default=True, type=click.Path())

def main_fun(dr,ming,fout,ftp=None):
    """Generate sequence types."""
    # start timer
    start = time()
    # state arguments
    click.echo("""Arguments will take the following values:\n\tdr - {}\n\tming - {}""".format(dr,ming) )
    ## paths
    faln = os.path.join( dr, 'genome_aln_unq.fna')
    fmd = os.path.join( dr, 'genome_info.csv')
    fout = os.path.join( dr, fout)
    
    # input file checks
    if not os.path.exists(faln):
        raise FileNotFoundError("couldn't find the alignment file %s."%faln)

    if not os.path.exists(fmd):
        raise FileNotFoundError("couldn't find the metadata file %s."%fmd)

    # full path to output
    # check if we should proceed in case output is present
    click.echo("Output will be saved to path: %s"%fout)
    if not utils.outcheck(fout):
        return

    # if a list of barcoding position was not provided, then use the list bundled with cova
    if ftp is None:
        click.echo("No position list provided! Using cova' internal typing positions.")
        tpos = TYPEPOS
    else:
        click.echo("A position list was provided! Loading typing positions from %s"%ftp)
        with open(ftp) as flob:
            tpos = [ int(i.strip('\n')) for i in flob]

    # sort typing positions
    tpos = sorted(list(set(tpos)))
    
    ## data
    # whole-genome Multiple Sequence Alignment
    aln = AlignIO.read( faln,'fasta')
    # list of genome ids in the alignment
    ids = [ i.id for i in aln] 
    # metadata on above genome ids
    metadata = pandas.read_csv(fmd)
    nrow = len(metadata)
    # empty dict of genome ids and their collection date
    seq_date = {}

    # for every genome entry
    for x in range(nrow):
        a = metadata.at[x,'accession']
        d = metadata.at[x,'date']
        # if the date is of appropriate format, capture it
        try:
            seq_date[a] = datetime.datetime.strptime(d,'%Y-%m-%d').date()
        except ValueError:
            #print("Bad format date!")
            # or else, skip to next entry
            continue

    # dict of a barcoding sequence as key and a list of genomes with that sequence as its value
    idst = barcode_genomes(aln,pos=tpos,m=ming)
    click.echo("%s: Barcodes have been generated."%utils.timer(start))
    # empty dict of sequence type and the earliest collection date of a genome of this type
    st_mindate = {}

    # for every entry in the dict of barcodes and their list of genomes
    for k,v in idst.items():
        # list of date objects
        dobs = [ seq_date[i] for i in v if i in seq_date.keys() ]
        # if there is least 1 genome with a valid date
        if len(dobs) != 0:
            # the earliest date
            mindate = min(dobs)
            # make an entry in the above dict
            st_mindate[k] = mindate 

    click.echo("%s: Earliest date for each type has been identified."%utils.timer(start))
    # the above keys sorted in the ascending order of dates
    sorted_sts = sorted( st_mindate, key=lambda x:st_mindate[x])
    # no. of unique sequence tags
    N = len(sorted_sts)
    # dict of sequence type id and its sequences
    stid_seq = pandas.DataFrame([ [ 'ST'+str(x+1)] + list(i) for x,i in enumerate(sorted_sts) ],
        columns = ['type'] + tpos)
    # write output to file
    stid_seq.to_csv(fout, index=False)
    click.echo("%s: Output saved. Done."%utils.timer(start))

if __name__ == '__main__':
    main_fun()
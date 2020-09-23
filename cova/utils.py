## require
import os, re, numpy, csv, colorsys, math
from time import time
from itertools import groupby
from joblib import Parallel, delayed
from collections import Counter
from scipy.special import binom
from Bio import AlignIO
from Bio.Data.CodonTable import unambiguous_dna_by_id as codon_table
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from . import _divs

### Exceptions ###
class LowercaseSeqError(Exception):
    pass
class LenSeqError(Exception):
    pass
class EmptySeqError(Exception):
    pass
class QualSeqError(Exception):
    pass

### Functions ###
def timer(b=0):
    """Print elapsed time in hh:mm:ss format given starting time ( using time module)."""
    if b == 0:
        b=time()
    
    rts = round(time() - b)
    if rts < 60:
        out = '00:00:'+str(rts)
    else:
        rtm_frac,rtm_whole = math.modf(rts/60)
        rtm = rtm_whole
        rtms = round(rtm_frac*60)
        if rtm < 60:
            out = '00:'+str(rtm)+':'+str(rtms)
        else:
            rth_frac,rth_whole = math.modf(rtm/60)
            rth = rth_whole
            rthm = round(rth_frac*60)
            out = str(rth)+':'+str(rthm)+':00'
    return out

def readcsv(fl,sep=',',header=None,comments='#'):
    """
    Return a nested list of lists from a CSV file. 
    
    Arguments:
    - fl       - full path to input file
    - sep      - field seperator                                    [ ',' ]  
    - header   - boolean, if passed, first row is used as a header
    - comments - lines starting with this character are comments    [ '#' ]
    """
    with open(fl) as flob:
        rdob = csv.reader( filter(lambda x: x[0]!=comments, flob), delimiter=sep )
        data = [i for i in rdob]
    if header==True:
        head = data[0]
        data = data[1:]
        return (head,data)
    else:
        return data

def writecsv(fl,data,header=None,sep=',',comments='',mode='w'):
    """
    Write a nested list of lists to a CSV file. 
    
    Arguments:
    - fl       - full path to output file
    - data     - nested list of lists
    - sep      - field seperator                                    [ "," ]  
    - header   - list of column names for the table                 [ None]
    - comments - lines to be written before the actual data         [ '#' ]
    """
    with open(fl,mode) as flob:
        flob.write(comments)
        wrob = csv.writer(flob,delimiter=sep)
        if header != None: wrob.writerow(header)
        for i in data: wrob.writerow(i)

def outcheck(path):
    """Decide whether to proceed or not, if output is already present."""
    proceed = True
    
    if os.path.exists(path):
        response = input('''output path "%s" already exists. 
Do you want to proceed and rewrite? [y/n]\n'''%path)
    
        if response == 'n':
            proceed = False
        elif response == 'y':
            proceed = True
        else:
            raise ValueError('''Inappropriate input! Please respond with [y] to
proceed and overwrite OR with [n] to exit.''')
        
    if not proceed:
        print("okay! Exiting.")
    
    return proceed

def split_data( data, ix, cixs):
    """
    Split a table (nested list of lists) into a dict given a column index
    where each key is a unique entry from the column and its value is the
    list of rows with that entry in the column.    
    
    Arguments:
    - data  - nested list of lists
    - ix    - splitter column index  
    - cixs  - either a single index OR a list of indices of the columns to keep 
    """
    keyfunc = lambda v:v[ix]
    data = sorted( data, key = keyfunc)  
    out = { k:[ i for i in g] for k, g \
            in groupby( data, keyfunc) }
    tp = type(cixs)
    if tp is int:
        out = { k:[ i[cixs] for i in v] for k,v in out.items()}
    else:
        out = { k:[ [ i[x] for x in cixs] for i in v] for k,v in out.items()}
    return out

def list_ep(als,step=1):
    """
    Return a nested list of ends of consecutive runs given an ordered list of numbers.

    Arguments:
    - als   - an ordered list of numbers
    - step  - step size, to define "consecutive runs" [ 1 ]
    """
    # first check if the input list is sorted in ascending order
    try:
        assert als == sorted(als)
    except AssertionError:
        print("input list needs to be in ascending order")
        return
        
    # initialize output list
    out = []
    # no. of elements
    n = len(als)
    # last no. stored, initialize with one number less than the first
    lns = als[0]-1
    # for every number
    for x in range(n):
        # first element
        e1 = als[x]
        # if it is not greater than the last stored no., 
        # then nothing more to do for this no.
        if e1 <= lns:
            continue
        # initial range end-points for this no.
        r1 = e1
        r2 = e1
        # initialize expected no.
        expn = e1
        # compare this no. with all subsequent no.s
        for y in range(x+1,n):
            # second element
            e2 = als[y]
            # expected no. increment by step size in every loop
            expn += step
            # if the next no. matches the expectation
            if e2 == expn:
                # update right end-point of the range with the current no.
                r2 = e2
            # when you reach a no. that doesn't
            else:
                # update the same end-point with the previous no.
                r2 = als[y-1]
                break
        # store current end-points
        out.append((r1,r2))
        # update the last number stored
        lns = out[-1][1]
    return out

def is_seq_qual(rec,ambt,chars='ACGT'):
    """
    Test if sequence is of good quality based on the proportion of ambiguous characters.

    Arguments:
    - rec   - biopython seqrecord
    - ambt  - ambiguous characters threshold (%)
    - chars - string of allowed characters
    """
    # no. of sites in the record
    nsites = len(rec)
    # no. of ambiguous characters
    rec_ambc = sum([ i not in chars for i in rec.seq])
    # no. of gaps ( if any)
    ngaps = rec.seq.count('-')
    # percentage of ambiguous characters
    rec_ambp = rec_ambc/(nsites-ngaps)*100
    # boolean output
    if rec_ambp > ambt:
        return False
    else:
        return True

def n2c(nseq,chars='ACGTU-',thres=90):
    """Return a list of codons given a string of nucleotides."""
    ## checks 
    # length is a multiple of 3
    l = len(nseq)
    
    if l%3 != 0:
        raise LenSeqError('Sequence length must be a multiple of 3')
    
    # sequence is largely of unambiguous nucleotides
    pc = sum(i in chars for i in nseq)/l*100
    
    if pc < thres:
        raise QualSeqError('Sequence has many non-standard characters')
    
    cseq = [ ''.join(nseq[x:x+3]) for x in range(0,l,3)]
    return cseq

def nv2av(p,v,seq,codon_table=codon_table[1]):
    """
    Generate amino acid variant(s) of a given nucleotide variant.
    
    Arguments:
    - p             - 0-indexed starting position of the
                        variant nucleotide in the CDS
    - v             - sequence of the variant, a string 
                        of len >= 1
    - seq           - CDS as a list of codons
    - codon_table   - Biopython codon table                 [ Standard]    
    
    Value:
        a nested list to accomodate changes in the multiple codons
        each entry is a list of 3, with following value at the indices
        0: ancestral codon
        1: variant codon
        2: amino acid variant label, with 3 components
            1: ancestral aa
            2: 1-indexed position on the aa sequence
            3: variant aa
    """
    ## checks
    # position is integer and non-negative
    if type(p) is not int or p < 0:
        raise ValueError('Position must be a non-negative integer')
    # variant is in capital letters
    var_ucase = v.isupper()
    if not var_ucase:
        raise LowercaseSeqError('Variant sequence must be in Uppercase')
    # variant is a string of dna bases
    not_base = any(i not in 'ACGT' for i in v)
    if not_base:
        raise ValueError('v must be an unambiguous DNA base')
    # CDS is in capital letters
    cds_ucase = all( i.isupper() for i in seq)
    if not cds_ucase:
        raise LowercaseSeqError('CDS sequence must be in Uppercase')
    # the CDS is a list of codons
    not_codon = any( any(j not in 'ACGT' for j in i) for i in seq) 
    if not_codon:
        raise ValueError('seq must be a list of codons')
    # length check
    codonlen = all(len(i) == 3 for i in seq)
    if not codonlen:
        raise LenSeqError('All seq elements must be of length 3')
    # the position along with the variant sequence does not exceed cds length
    if p+len(v) > len(seq)*3:
        raise LenSeqError('''Position along with the variant sequence must not exceed CDS length.
            pos = {}, variant = {}, CDS length = {}'''.format(p,v,len(seq)*3)) 

    # codon table attributes
    cotab = codon_table.forward_table
    stopc = codon_table.stop_codons
    
    # variant length
    l = len(v)
    # indices of nucleotide variants
    nixs = list(range(p,p+l))
    # corresponding codon, base indices and variant nucleotides
    cbixs = [ ( math.floor(nixs[x]/3), nixs[x]%3, i) for x,i in enumerate(v) ]
    # split above by codon_indices
    codon_bixs = split_data(data=cbixs, ix=0, cixs=[1,2])

    # initialize output list
    out = []
    # for every modified codon and its modified bases
    for cx,bxs in codon_bixs.items():
        # ancestral codon
        acod = seq[cx]
        # if ancestral codon is a stop codon, then ancestral amino acid is *
        if acod in stopc:
            aa = '*'
        else:
            aa = cotab.get(acod)

        # generate variant codon
        vcod = list(acod)
        for bx in bxs:
            vcod[ bx[0] ] = bx[1]
        vcod = ''.join(vcod)

        # variant amino acid
        if vcod in stopc:
            av = '*'
        else:
            av = cotab.get(vcod)
        
        lab = aa + str(cx+1) + av
        out.append([ acod, vcod, lab])
    return(out)

def extract_cds(begin,end,msa,has_stop,outf):
    """
    Extract a CDS MSA from a whole-genome (biopython) MSA.
    
    Arguments:
    - begin     - 1-indexed starting position of the CDS on the MSA
    - end       - 1-indexed last position of the CDS on the MSA
    - msa       - Whole-genome Biopython MultipleSeqAlignment
    - has_stop  - boolean, if the terminating stop codon is included
    - outf      - full path to the output file for alignment
    """
    ## checks
    # start and end positions are integers
    if type(begin) is not int or type(end) is not int:
        raise TypeError('end positions must be integers.')
    # msa is a biopython alignment
    if type(msa) is not MultipleSeqAlignment:
        raise TypeError('msa must be a biopython MultipleSeqAlignment.')
    # end > begin
    if end <= begin:
        raise ValueError('end must be greater than begin.')
    
    # extraction differs based on whether or not the terminating stop codon is included
    if has_stop:
        cds = msa[ :, begin-1:end-3]
    else:
        cds = msa[ :, begin-1:end]
    # remove description
    for i in cds:
        i.description = ''
    # write alignment to file
    AlignIO.write(alignments=[cds], handle=outf, format='fasta')
    return cds

def genome_seqtype(aln):
    # list of column sequences in the aligment at these positions
    pseqs = [ aln[:,p-1] for p in TYPEPOS]
    # list of row suquences limited to these positions
    id_st = []
    for x,i in enumerate(zip(*pseqs)):
        seq = ''.join(i).upper()
        if seq in SEQTYPES.keys():
            entry = [ aln[x].id, SEQTYPES[seq] ]
            id_st.append(entry)
    return id_st

# (https://stackoverflow.com/questions/876853/generating-color-ranges-in-python)
def get_N_HexCol(N):
    """Return 'N' hexadecimal color codes, uniformly spaced over the entire range."""
    HSV_tuples = [(x/ N, 0.6, 0.6) for x in range(N)]
    hex_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x * 255), colorsys.hsv_to_rgb(*rgb))
        hex_out.append('#%02x%02x%02x' % tuple(rgb))
    return hex_out

##################################################################


### Classes ###
class MSA(object):
    """
    Class for Multiple Sequence Alignments (MSA). 
    Internally, it uses Biopython alignment and numpy array representations. 
    To modify alignments, both rowwise and columnwise, as well as to perform 
    variant calling and diversity computations. 
    
    To initialize, you need to provide full path to an MSA file.

    Attributes:
    - aln       - Biopython MultipleSeqAlignment object
    - array     - Numpy array form of the MSA
    - nseq      - number of sequences in the alignment
    - ncol      - number of columns in the alignment
    - gaps      - list of gapped column indices
    - nogapcol  - number of columns without gaps
    - difs      - list of non-gapped variant columns
    - ndif      - number of variant columns

    Methods:
    - remove_gapcols - return MSA array excluding gapped columns
    - pdistmat       - return pairwise distance matrix
    - ndiv           - compute nucleotide diversity 
    - limref         - limit MSA to sites present in a given reference
    - rmdup          - remove duplicate sequences from the MSA 
    - ins            - find insertion variants relative to a reference
    - dels           - find deletion variants relative to a reference
    - pointmuts      - find point mutations relative to a reference
    """
    def __init__(self,fname):
        """
        Initialize a new MSA object.

        Arguments:
        - fname     - full path to input MSA file 
        """
        try:
            # read
            self.aln = AlignIO.read(handle=fname,format='fasta')
            # if file is empty or in the wrong format, biopython will handle it
            # But it will also pass the file even if it has just a header
            has_empty_rec = any( len(i) == 0 for i in self.aln)
            if has_empty_rec:
                raise EmptySeqError('record(s) are missing sequence(s).')
            # Even if a record has sequence, it is not a valid MSA unless we have at least 2 records
            if len(self.aln) <= 1:
                raise ValueError('need at least 2 sequences for a valid MSA.')
            # list of sequence ids
            self.ids = [ i.id for i in self.aln]
            # ids must be unique
            if len(self.ids) != len(set(self.ids)):
                raise ValueError('sequence IDs must be unique.')

            # capitalize ( in any case)
            for i in self.aln:
                i.seq = i.seq.upper()

            # numpy array form of the MSA
            self.array = numpy.array( [ list(rec) for rec in self.aln])
            # number of sequences and number of sites
            self.nseq, self.ncol = self.array.shape
        
            # indices of columns with gaps
            self.gaps = [ x for x in range(self.ncol) \
                            if '-' in self.array[:,x]]
            # no. of columns without gaps
            self.nogapcol = self.ncol - len(self.gaps)
            # indices of columns with differences but no gaps
            self.difs = [ x for x in range(self.ncol) if \
                            x not in self.gaps and \
                            len( set( self.array[:,x])) > 1]
            # no. of variant sites
            self.ndif = len(self.difs)
        
        except (FileNotFoundError, IsADirectoryError, PermissionError) as ex:
            raise IOError("failed to open the file: %s."%str(ex))
            
    def remove_gapcols(self):   
        """Remove columns with gaps and return array."""
        arr = _divs.rmc_ef(ar = self.array, e='-', f=0)
        return arr

    def pdistmat(self,ncpu=1):
        """
        Return a pairwise distance matrix as a numpy array.
        Distance is the fraction of differing sites between two sequences.

        Optional argument:
        - ncpu  - number of CPUs to be used for parallel computation
        """
        ## checks
        # need to have at least 1-ungapped column
        if self.nogapcol == 0:
            raise EmptySeqError('sequence array needs at least 1-ungapped column.')

        # array after removing columns with gaps
        arr = self.remove_gapcols()
        # number of sequences
        N = self.nseq
        # number of columns Without gaps
        L = self.nogapcol

        ## function to parallelize    
        def fun(x,y,arr=arr,L=L):
            if x == y:
                out=0
            elif x > y:
                out=-1
            else:
                # number of differing elements between two rows of the array
                nd = sum( arr[x,:] != arr[y,:])
                # p-distance: dividing above by no. of sites
                out = round(float(nd)/L,4)
            return out

        # run the above function parallely to get p-distances
        data = Parallel(n_jobs=ncpu) ( delayed(fun) (x,y) for x in range(N) for y in range(N))
        dmat = numpy.array(data, dtype = float).reshape(N,N)     

        return dmat
    
    def limref(self,ref):
        """
        Limit MSA to sites present in a given reference. By excluding columns 
        which have gaps in the row corresponding to the reference.
        The reference must be present in the MSA.
        """
        ## checks
        msa = self.aln
        # ref is present in the MSA
        if ref not in self.ids:
            raise ValueError('reference must be present in the MSA.')

        # reference record
        refrec = [ i for i in msa if i.id == ref][0]

        # check for presence of gaps in the reference record
        # if no gaps, then the alignment is already reference limited
        if '-' not in refrec:
            print("No gaps in the reference. Returning as it is!")
            return msa 

        # list of sequences from the non-gapped columns
        ngcs = [  msa[:,x] for x, i in enumerate(refrec.seq) if i != '-']
        # list of modified sequences for sequence records
        modseqs = [ ''.join(i) for i in zip(*ngcs) ]
        for x,i in enumerate(msa):
           i.seq = Seq(modseqs[x])  
        return msa

    def rmdup(self,ref=None):
        """
        Remove duplicate sequences from the MSA.

        Optional arguments:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
        - ref - reference accession [ None ]

        If reference is not provided, input is the original biopython MSA
        else, input is the MSA returned by self.limref().
        """
        if ref is None:
            old_msa = self.aln
        else:
            old_msa = self.limref(ref)
            
        # list of unq seq records
        unqrec = []
        # set of sequences seen
        seenseqs = set()
        
        for i in old_msa:
            if str(i.seq) not in seenseqs:
                unqrec.append(i)
            seenseqs.add( str(i.seq))
            
        new_msa = MultipleSeqAlignment(records=unqrec)
        
        # nested list of identical records
        genome_dups = []
        for rec in new_msa:
            idenrecs = [ self.ids[x] for x,i in enumerate(old_msa) if i.id != rec.id \
                and i.seq == rec.seq]
            if len(idenrecs) > 0:
                genome_dups.append( [rec.id, ','.join(idenrecs)])
        return ( new_msa, genome_dups)
    
    def ndiv(self,ref=None):
        """
        Compute nucleotide diversity:
            Average pairwise difference of sequences per site
            
            if an MSA has L sites and N sequences, 
            and d_ij represents number of differing sites between two sequences i and j
            then the formula for nucleotide diversity is
                sum(d_ij) / ( choose(N,2) * L )
                    where i : 0 - (N-1)
                      and j : i+1 - N

        It achieves a substantial speed gain by splitting calculation over columns
        as opposed to first calculating pairwise row differences. 
        """
        if ref is None:
            msa = self.aln
        else:
            msa = self.limref()
        return _divs.nucdiv_aln( msa, t=0)

    def slide_ndiv(self,window,jump,ncpu=1, ref=None):
        
        if ref is None:
            aln = self.aln
        else:
            aln = self.rmdup(ref)[0]
        
        out = _divs.nucdiv_slide(msa=aln,window=window,jump=jump,ncpu=ncpu)
        return out

    def ins(self,ref,ambt=10,chars='ACGT-',header=False):
        """
        Identify positions and sequence of insertions across variants relative to a given reference.
        
            Arguments:
        - ref     - reference accession
        - ambt    - ambiguous characters threshold,                 [ 10 (%) ]
                    used both to check sequence quality,
                    and to call insertions
        - chars   - string of expected characters in the msa        [ 'ACGT-']
        - header  - boolean, if a header should also be returned    [ False  ]       
        
        Value:
        a nested list, to be interpreted as a table, where rows are variant positions.
        The first two columns are for reference position and allele respectively,
        and the following columns are for the alleles of other sequences in the MSA.

        if header is True,
            a 2-tuple (Table,Header)
        """
        ## checks
        # ref is present in the MSA
        if ref not in self.ids:
            raise ValueError('reference must be present in the MSA.')
        # threshold is numeric
        if type(ambt) not in (int,float):
            raise TypeError('threshold must be numeric.')
        # threshold is non-negative
        if ambt < 0:
            raise ValueError('threshold must be non-negative.')

        # MSA
        msa = self.aln; nseq = self.nseq; nsites = self.ncol
        # list of genome ids
        genomes = self.ids

        # reference record
        refx, refrec = [ (x, msa[x]) for x,i in enumerate(genomes) if i == ref][0]
        # reference sequence length
        rlen = len(refrec) - refrec.seq.count('-')

        # dict of records with insertions relative to the reference
        rec_ins = {}
        for c,rec in enumerate(msa):
            # skip the reference sequence
            if rec.seq == refrec.seq:
                #print("\tskipping the reference sequence")
                continue

            # indices of sites where rec has a letter but ref has a gap
            ixs = [ x for x,i in enumerate(rec) if refrec[x] == '-' and i != '-' ]
            # skip, if no insertion sites 
            if len(ixs) == 0:
                #print("\tno insertions in %s. Moving on to the next record.."%genomes[c])
                continue

            # convert above list to a list of 2-tuples marking ranges of consecutive stretches
            eps = list_ep(ixs)

            for i in eps:
                # insertion sequence 
                gapseq = rec[i[0]:i[1]+1].seq
                l = len(gapseq)
                # count of ambiguous bases in this sequence
                ambc = sum(j not in chars for j in gapseq)
                ambp = round(ambc/l*100)
                # ignore, if more ambiguous letters than the set threshold
                if ambp > ambt:
                    continue

                # find insertion position in the reference
                # using no. of bases covered on reference before reaching this position on the alignment 
                rx = sum( j != '-' for j in refrec[:i[0]])-1
            
                # ignore end insertions
                if rx == -1 or rx == rlen-1:
                    # print("\tIt's an end-point insertion. Skipping!")
                    continue

                entry = [ (rx, refrec[i[0]-1]), rec[i[0]-1]+str(gapseq)]
            
                if rec.id not in rec_ins.keys():
                    rec_ins[rec.id] = [entry]
                else:
                    rec_ins[rec.id].append(entry)


        # list of all variants
        allv = list( set(i[0] for v in rec_ins.values() for i in v))
        nv = len(allv)

        # initialize output tab with these variants
        vtab = [ [ allv[r][1] for c in range(nseq) ] for r in range(nv)]

        # go through all variants
        for x in range(nv):
            # extract reference position and base
            pos,base = allv[x]
            # head, to be added later to this row
            head = [ pos+1, base]
            # go through dict of genomes and their variants
            for k,v in rec_ins.items():
                # if genome has a variant at this position
                for i in v:
                    if allv[x] in i:
                        # find the genomes.index in the MSA
                        y = genomes.index(k)
                        # make an entry in the variant table
                        vtab[x][y] = i[1]
                        # stop iterating over this genome's list of variants
                        break 
            vtab[x] = head + vtab[x]

        # remove reference from the above table
        tmp = list(zip(*vtab))
        tmp = [ i for x,i in enumerate(tmp) if x != refx+2]
        vtab = list(zip(*tmp))
        vtab = sorted(vtab, key = lambda x: x[0])
        if header is False:
            return vtab
        else:
            # header
            head = ['pos','ref'] + [ i for i in genomes if i != ref]
            return (vtab,head)

    def dels(self,ref,header=False):
        """
        Identify positions and sequence of deletions across variants relative to a given reference.
        
        Arguments:
        - ref       - reference accession
        - header    - boolean, if a header should also be returned [ False ]
        
        Value: 
        a nested list, to be interpreted as a table, where rows are variant positions
        The first two columns are for reference position and allele respectively,
        and the following columns are binary representations of Presence/Absence 
        of the deletion across sequences.

        if header is True,
            a 2-tuple (Table,Header)
        """
        # msa limited to reference positions
        msa = self.limref(ref)
        # no. of isolates
        nseq = self.nseq
        # no. of sites
        nsites = len(msa[0])

        # list of genome ids
        genomes = self.ids
        # reference record
        refrec = [ i for i in msa if i.id == ref][0]
        # reference index
        rx = [ x for x,i in enumerate(msa) if i.id == ref][0]

        # initialize dict of records and their gaps coordinates
        rec_gaps = {}
        for x,rec in enumerate(msa):
            # skip the reference sequence
            if rec.seq == refrec.seq:
                #print("\tskipping the reference sequence")
                continue

            # list of all gap indices
            ixs = [ y for y, i in enumerate(rec.seq) if i == '-']
            # skip, if no gaps
            if len(ixs) == 0:
                #print("\t%s has no deletion. Moving on to the next record"%genomes[x])
                continue
            # convert above list to a list of 2-tuples marking ranges of consecutive stretches
            eps = list_ep(ixs)
            
            # ignore end deletions
            eps = [ i for i in eps if not (i[0] == 0 or i[1] == nsites-1)]
            if len(eps) == 0:
                #print('''Deletions were only present at sequence ends.
                #Disregarded as potential sequencing errors''')
                continue

            rec_gaps[rec.id] = []
            for i in eps:
                # stating and ending indices of the gap
                b = i[0]
                e = i[1]+1
                # sequence deleted from the reference
                delseq = str(refrec.seq[b:e])
                l = len(delseq)

                entry = (b, delseq, l)
                rec_gaps[rec.id].append(entry)
        
        # set of variant positions and their length
        vposlen = set( i for v in rec_gaps.values() for i in v)
        # convert above to a dict
        pos_lens = split_data(data=vposlen,ix=0,cixs=[1,2])
        # only keep positions with unique insertions
        allv = [ (k,v[0][0],v[0][1]) for k,v in pos_lens.items() if len(v) == 1] 
        nv = len(allv)
        
        # initialize output tab with these variants
        vtab = [ [ 0 for c in range(nseq) ] for r in range(nv)]

        # go through all variants
        for x in range(nv):
            # extract reference position, sequence and length
            pos,seq,l = allv[x]
            # head, to be added later to this row
            head = [ pos+1, seq, l]
            # go through dict of genomes and their variants
            for k,v in rec_gaps.items():
                # if genome has a variant at this position
                if allv[x] in v:
                    # find the genomes.index in the MSA
                    y = genomes.index(k)
                    # make an entry in the variant table
                    vtab[x][y] = 1
            vtab[x] = head + vtab[x]

        # remove reference from the above table
        tmp = list(zip(*vtab))
        tmp = [ i for x,i in enumerate(tmp) if x != rx+3]
        vtab = list(zip(*tmp))
        vtab = sorted(vtab, key = lambda x: x[0])
        if header is False:
            return vtab
        else:
            # header
            head = ['pos','ref','len'] + [ i for i in genomes if i != ref]
            return (vtab,head)

    def pointmuts(self,ref,header=False,chars='ACGT'):
        """
        Identify positions and sequence of point mutations across variants relative to a given reference
        
        Arguments:
        - ref       - reference accession
        - header    - boolean, if a header should also be returned  [ False ]
        - chars     - string of allowed variant characters          [ 'ACGT']
        
        Value: 
        a nested list, to be interpreted as a table, where rows are variant positions
        The first two columns are for reference position and allele respectively,
        and the following columns are for the alleles in rest of the sequences.

        if header is True,
            a 2-tuple (Table,Header)
        """

        # msa limited to reference positions
        msa = self.limref(ref)
        # no. of isolates
        nseq = self.nseq
        # no. of sites
        nsites = len(msa[0])

        # list of genome ids
        genomes = self.ids
        # reference record
        refrec = [ i for i in msa if i.id == ref][0]
        # reference index
        rx = [ x for x,i in enumerate(msa) if i.id == ref][0]

        ## dict of records and their point mutations
        rec_pms = {}
        for x,rec in enumerate(msa):
            # skip the reference sequence
            if rec.seq == refrec.seq:
                #print("\tskipping the reference sequence")
                continue

            # list of indices of all variants
            ixs = [ y for y,i in enumerate(rec.seq) if i != refrec[y] and i in chars]    
            # skip, if no variants
            if len(ixs) == 0:
                #print("\tno point mutations in %s. Moving on to the next record.."%genomes[x])
                continue

            # convert above list to a list of 2-tuples marking ranges of consecutive stretches
            eps = list_ep(ixs)
            
            # ignore mutations at either end
            eps = [ i for i in eps if not (i[0] == 0 or i[1] == nsites-1)]
            if len(eps) == 0:
                #print('''Point mutations were only present at sequence ends. Disregarded as
                #potential sequencing errors. Skipping this record!''')
                continue

            rec_pms[rec.id] = eps

        # list of all variant positions
        allv = list( set( i for v in rec_pms.values() for i in v))
        nv = len(allv)

        # initialize output tab with these variants
        vtab = [ [ str(refrec.seq[ allv[r][0]:allv[r][1]+1 ]) for c in range(nseq) ] for r in range(nv)]

        # go through all variants
        for x in range(nv):
            # extract reference position
            pos = allv[x][0]
            # sequence
            rseq = vtab[x][0]

            # head, to be added later to this row
            head = [ pos+1, rseq]
            # go through dict of genomes and their variants
            for k,v in rec_pms.items():
                # if genome has a variant at this position
                if allv[x] in v:
                    # find the genomes.index in the MSA
                    y = genomes.index(k)
                    # make an entry in the variant table
                    vtab[x][y] = str(msa[y][ pos:allv[x][1]+1 ].seq)
            vtab[x] = head + vtab[x]

        # remove reference from the above table
        tmp = list(zip(*vtab))
        tmp = [ i for x,i in enumerate(tmp) if x != rx+2]
        vtab = list(zip(*tmp))
        vtab = sorted(vtab, key = lambda x: x[0])
        
        if header is False:
            return vtab
        else:
            # header
            head = ['pos','ref'] + [ i for i in genomes if i != ref]
            return (vtab,head)
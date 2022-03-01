## require
import os, re, numpy, csv, colorsys, math, pandas, multiprocessing
from time import time
from itertools import groupby
from collections import Counter
from scipy.special import binom
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from . import _divs, RFS, CODONTABLE, GENOME, REF

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
    
    # if there are no elements, terminate
    if n == 0:
        return []
    
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

def extract_dates_gisaid_headers(fin):
    """Return a table of genome accession and their collection dates from gisaid headers of the input FASTA files."""
    headers = list( SeqIO.index( fin, 'fasta').keys())
    tab = pandas.DataFrame([ i.split('|') for i in headers],columns=['accession','gisaid','date'])
    tab = tab[['accession','date']]
    tab.accession = headers
    tab = tab.dropna()
    return tab

def get_recgaps(rec):
    """Return 0-indexed start and end positions of sequence gaps (-)
    in a given biopython sequence record."""
    # no. of sites
    nsites = len(rec.seq)
    # list of all gap indices
    ixs = [ y for y, i in enumerate(rec.seq) if i == '-']
    # convert above list to a list of 2-tuples marking ranges of consecutive stretches
    eps = list_ep(ixs)
    # ignore end deletions
    eps = [ i for i in eps if not (i[0] == 0 or i[1] == nsites-1)]
    return eps

def get_recmups(rec,refrec,chars='ACGT'):
    """
    Return 0-indexed start and end positions of mutated positions
    in a given sequence record.

    Arguments:
    - rec       - biopython sequence record
    - refrec    - reference sequence record
    - chars     - string of allowed characters
    """
    # no. of sites
    nsites = len(rec.seq)
    # list of indices of all variants
    ixs = [ y for y,i in enumerate(rec.seq) if i != refrec[y] and i in chars]    
    # convert above list to a list of 2-tuples marking ranges of consecutive stretches
    eps = list_ep(ixs)    
    # ignore mutations at either end
    eps = [ i for i in eps if not (i[0] == 0 or i[1] == nsites-1)]
    return eps

def get_recins(rec,refrec):
    """
    Return 0-indexed start and end positions of insertions in a given sequence record.

    Arguments:
    - rec       - biopython sequence record from an MSA
    - refrec    - reference sequence record from the same MSA
    """
    # indices of sites where rec has a letter but ref has a gap
    ixs = [ x for x,i in enumerate(rec) if refrec[x] == '-' and i != '-' ]
    # convert above list to a list of 2-tuples marking ranges of consecutive stretches
    eps = list_ep(ixs)
    return eps

def is_seq_qual(seq,ambt,chars='ACGT'):
    """
    Test if sequence is of good quality based on the proportion of ambiguous characters.

    Arguments:
    - seq   - a string of letters
    - ambt  - ambiguous characters threshold (%)
    - chars - string of allowed characters
    """
    # no. of sites in the record
    nsites = len(seq)
    # no. of ambiguous characters
    seq_ambc = sum([ i not in chars for i in seq])
    # no. of gaps ( if any)
    ngaps = seq.count('-')
    
    # percentage of ambiguous characters
    try:
        seq_ambp = seq_ambc/(nsites-ngaps)*100
    except ZeroDivisionError:
        return False

    # boolean output
    if seq_ambp > ambt:
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

def nv2av(p,v,seq):
    """
    Generate amino acid variant(s) of a given nucleotide variant.
    
    Arguments:
    - p             - 0-indexed starting position of the
                        variant nucleotide in the CDS
    - v             - sequence of the variant, a string 
                        of len >= 1
    - seq           - CDS as a list of codons
    
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
    # position is non-negative
    if  p < 0:
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
    cotab = CODONTABLE.forward_table
    stopc = CODONTABLE.stop_codons
    
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

def get_stopm(s):
    """Return nonsense mutations from a list of amino acid changes."""
    stopms = [ i for i in s if bool(re.search('_[A-Z]\d+\*',i))]
    out = ','.join(stopms)
    return out

def extract_cds_msa(p,pstart,pstop,msa):
    """
    Extract a CDS MSA from a whole-genome (biopython) MSA.
    
    Arguments:
    - p         - protein accession
    - pstart    - 1-indexed starting position of the CDS on the MSA
    - pstop     - 1-indexed last position of the CDS on the MSA
    - msa       - Whole-genome Biopython MultipleSeqAlignment
    """
    
    # start and end positions are integers
    if type(pstart) is not int or type(pstop) is not int:
        raise TypeError('end positions must be integers.')
    
    # msa is a biopython alignment
    if type(msa) is not MultipleSeqAlignment:
        raise TypeError('msa must be a biopython MultipleSeqAlignment.')
    
    # pstop > pstart
    if pstop <= pstart:
        raise ValueError('end must be greater than begin.')

    b = pstart-1
    e = pstop
    
    # for regular proteins 
    if p not in RFS.index:
        cds = msa[:,b:e]
    # for proteins affected by ribosomal slippage
    else:
        rfs_pos = RFS.loc[p,'genomic_position']
        rfs_type = RFS.loc[p,'type']
        cds = msa[:,b:rfs_pos] + msa[:,rfs_pos+rfs_type:e]

    # remove stop codon if present
    if cds[0,-3:].seq in CODONTABLE.stop_codons:
        cds = cds[:,:-3]
    
    # remove description
    for i in cds:
        i.description = ''

    return cds

def get_shared_unique_elements(s,df):
    """Return shared and unique elements of a series given a binary dataframe."""
    b = df[s==1].sum(axis=1)>1
    p = b[b==True].index.values
    a = b[b==False].index.values
    out = pandas.Series([','.join(p),','.join(a)],index=['shared','unique'])
    return out

# (https://stackoverflow.com/questions/876853/generating-color-ranges-in-python)
def get_N_HexCol(N):
    """Return 'N' hexadecimal color codes, uniformly spaced over the entire range."""
    HSV_tuples = [(x/ N, 0.6, 0.6) for x in range(N)]
    hex_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x * 255), colorsys.hsv_to_rgb(*rgb))
        hex_out.append('#%02x%02x%02x' % tuple(rgb))
    return hex_out

def parallelize_df(df, func, n_cores=4):
    """Run a function parallely on a pandas dataframe.
    Source: https://towardsdatascience.com/make-your-own-super-pandas-using-multiproc-1c04f41944a1
    """
    df_split = numpy.array_split(df, n_cores)
    
    with multiprocessing.Pool(n_cores) as mp_pool:
        results = mp_pool.map(func, df_split)

    results = [i for i in results if i is not None]
    out = pandas.concat(results,ignore_index=True)
    return out

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
    def __init__(self,fname,ref=None):
        """
        Initialize a new MSA object.

        Arguments:
        - fname     - full path to input MSA file
        - ref       - reference accession included in the MSA ( Optional) 
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
            #if len(self.ids) != len(set(self.ids)):
            #    raise ValueError('sequence IDs must be unique.')

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
            # reference accession
            self.ref = ref
            # if a reference accession was provided
            if self.ref is not None:
                # reference record
                self.refrec = [ i for i in self.aln if i.id == ref][0]
                # variant records
                vrecs = pandas.Series( [ i for i in self.aln if i.seq != self.refrec.seq])
                vrecs.index = [i.id for i in vrecs.values]
                self.variant_recs = vrecs
        
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
        with multiprocessing.Pool(ncpu) as mp_pool:
            data = mp_pool.starmap( fun, [ (x,y) for x in range(N) for y in range(N)])
        dmat = numpy.array(data, dtype = float).reshape(N,N)     

        return dmat
    
    def limref(self):
        """
        Limit MSA to sites present in a given reference. By excluding columns 
        which have gaps in the row corresponding to the reference.
        The reference must be present in the MSA.
        """
        ## checks
        msa = self.aln
        # ref is present in the MSA
        if self.ref not in self.ids:
            raise ValueError('reference must be present in the MSA.')

        # reference record
        refrec = self.refrec

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

    def rmdup(self):
        """
        Remove duplicate sequences from the MSA.
        If reference is not provided, input is the original biopython MSA
        else, input is the MSA returned by self.limref().
        """
        if self.ref is None:
            old_msa = self.aln
        else:
            old_msa = self.limref()
            
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
    
    def ndiv(self):
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
        if self.ref is None:
            msa = self.aln
        else:
            msa = self.limref()
        return _divs.nucdiv_aln( msa, t=0)

    def slide_ndiv(self,window,jump,ncpu=1):
        
        if self.ref is None:
            aln = self.aln
        else:
            aln = self.rmdup()[0]
        
        out = _divs.nucdiv_slide(msa=aln,window=window,jump=jump,ncpu=ncpu)
        return out

    def ins(self,ambt=10,maxf=0.5):
        """
        Identify positions and sequence of insertions across variants relative to a given reference.
        
            Arguments:
        - ambt  - ambiguous characters threshold to keep insertions     [ 10 (%) ]
        - maxf  - maximum fraction of records at which an insertions 
                    fails quality check                                 [ 0.5 ]         
        Value:
        A pandas dataframe with 1 row per variant position and columns for 
        reference position, allele and variant alleles.
        """
        ## checks
        # ref is present in the MSA
        if self.ref not in self.ids:
            raise ValueError('reference must be present in the MSA.')
        # threshold is numeric
        if type(ambt) not in (int,float):
            raise TypeError('threshold must be numeric.')
        # threshold is non-negative
        if ambt < 0:
            raise ValueError('threshold must be non-negative.')

        # reference record
        refrec = self.refrec
        # reference sequence length
        rlen = len(refrec) - refrec.seq.count('-')
        # variant records
        vrecs = self.variant_recs
        # insertion intervals for every variant record
        inxs = vrecs.apply( lambda rec: get_recins(rec,refrec))
        # exclude records w/o insertions
        inxs = inxs[ inxs.apply(len)>0]
        # limit variant records to those w insertions
        vrecs = vrecs[inxs.index]
        # dataframe of ends of all insertions
        allis = pandas.DataFrame( set( j for i in inxs for j in i), columns=['start','end'])
        # corresponding sequences at these positions
        vseqs = allis.apply( axis=1, func=lambda x: vrecs.apply( 
        lambda y: str(y.seq[ x['start']-1:x['end']+1])))
        # find insertion position in the reference
        # using no. of bases covered on reference before reaching this position on the alignment 
        vseqs.index = allis.apply( axis=1, func=lambda x: ( sum( j != '-' for j in refrec[:x['start']])-1))
        # exclude insertions at either ends
        vseqs = vseqs[~vseqs.index.isin([-1,rlen-1])]
        # remove gap characters
        vseqs = vseqs.applymap( lambda x: x.split('-')[0])   
        # apply quality threshold
        b = vseqs.applymap( lambda x: is_seq_qual(x,ambt)) 
        # exclude variants based on majority-consensus
        vseqs = vseqs[ b.apply(axis=1, func=lambda x: sum(x==False))/b.shape[1] < maxf]
        # reference sequence at these positions
        rseq = [ GENOME[x] for x in vseqs.index]
        # output
        out = vseqs
        out.index = out.index+1
        out['ref'] = rseq
        out = out.sort_index()
        out = out[ ['ref']+vrecs.index.tolist()]
        return out

    def dels(self):
        """
        Identify positions and sequence of deletions across variants relative to a given reference.

        Value: 
        a pandas dataframe with 1 row per variant.
        The first three columns are for reference positions, length and sequence 
        respectively, and the following columns have P/A across variant genomes.
        """
        # check that MSA has a reference
        if self.ref is None:
            print('''Deletions can not be identified in the absence of a reference.
                Make sure a value was provided to the <ref> while initializing MSA.''')
            return

        # msa limited to reference positions
        msa = self.limref()
        # reference record
        refrec = self.refrec
        # variant records
        vrecs = self.variant_recs
        # list of gap indices for these genomes
        gapxs = vrecs.apply( lambda rec: get_recgaps(rec))
        # exclude empty entries
        gapxs = gapxs[ gapxs.apply(len)>0]  
        # pooled list of all gap intervals
        allgaps = list( set( j for i in gapxs for j in i))  
        # binary table of gaps P/A
        btab = pandas.DataFrame([ [ g in r for r in gapxs] for g in allgaps],
            columns=gapxs.index).astype(int)
        # starting position, length & sequence of gaps
        pls = pandas.DataFrame( [ [ i[0]+1, i[1]-i[0]+1, str(refrec.seq[i[0]:i[1]+1])]
            for i in allgaps], columns=['pos','len','seq']) 
        out = pls.T.append(btab.T).T
        out = out.sort_values(axis=0,by='pos')   
        return out

    def pointmuts(self,chars='ACGT',ambt=10,maxf=0.1):
        """
        Identify positions and sequence of point mutations across variants relative to a given reference
        
        Arguments:
        - chars - string of allowed variant characters                  [ 'ACGT']
        - ambt  - ambiguous characters threshold to keep a variant      [ 10 (%)]
        - maxf  - maximum fraction of records at which a variant 
                    fails quality check                                 [ 0.1   ]         

        Value: 
        a pandas dataframe with 1 row per variant position and columns for reference position, 
        allele, and for the variant alleles.
        """

        # msa limited to reference positions
        msa = self.limref()
        # reference record
        refrec = self.refrec
        # variant records
        vrecs = self.variant_recs
        # variant position intervals in each record
        pmxs = vrecs.apply( lambda rec: get_recmups(rec,refrec,chars))
        # exclude empty entries
        pmxs = pmxs[ pmxs.apply(len)>0]  
        # limit variant records to those w substitutions
        vrecs = vrecs[pmxs.index]
        # dataframe of ends of all variant positions
        allvs = pandas.DataFrame( set( j for i in pmxs for j in i), columns=['start','end'])
        # corresponding sequences at these positions
        vseqs = allvs.apply( axis=1, func=lambda x: vrecs.apply( 
        lambda y: str(y.seq[ x['start']:x['end']+1])))
        # reference sequence at the same positions
        vseqs['ref'] = allvs.apply( axis=1, func=lambda x: 
            str(refrec.seq[ x['start']:x['end']+1 ]))
        # arrange output
        out = pandas.concat([allvs,vseqs],axis=1)
        out['pos'] = out.start+1
        out = out[ ['pos','ref']+vrecs.index.values.tolist()]
        # apply quality threshold
        b = out.iloc[:,2:].applymap( lambda x: is_seq_qual(x,ambt))
        # filter out variants above maxf poor quality calls
        out = out[ b.apply(axis=1, func=lambda x: sum(x==False))/b.shape[1] < maxf] 
        out = out.sort_values(by='pos')
        return out
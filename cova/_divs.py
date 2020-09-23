import numpy
from scipy import special
from joblib import Parallel, delayed
from collections import Counter
import multiprocessing

def rmc_ef(ar,e,f):   
    """
    Remove columns from the array if an element's frequency is above threshold.

    Arguments:
    - ar    - 2-D numpy array
    - e     - element
    - f     - frequency threshold
    """
    # number of rows and columns
    nrow, ncol = ar.shape
    # frequency in counts
    fc = f*nrow
    # indices of columns where the element exceeds frequency threshold
    ixs = [ x for x in range(ncol) if list(ar[:,x]).count(e) > fc ]
    # delete above columns from the array
    out = numpy.delete( ar, ixs, 1) 
    return out

def pwd_ar1d(ar,alpha):
    """
    Return total pairwise difference of elements in a 1-D numpy array.

    Arguments:
    - ar        - 1-D numpy array
    - alpha     - set of acceptable alphabets, not necessarily characters
    """
    # set of letters
    sl = set(ar) & alpha
    
    # if 0 or only 1 acceptable alphabet is present at the position,
    # then its variation is zero
    if len(sl) <= 1:
        return 0
    
    # list of counts of all acceptable alphabets
    counts = [ numpy.count_nonzero(ar == i) for i in alpha] 
    # formula to get variation for the position
    out = sum([ i*sum(counts[x+1:]) for x,i in enumerate(counts)])
    return out

def nucdiv_aln(aln,alpha={'A','C','G','T'},g='-',t=0.1):
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
    
    Arguments:
    - aln       - biopython MSA object
    - ncpu      - number of CPUs for parallel processing [ 1 ]
    - alpha     - string of acceptable alphabets for the alignment [ {'A','C','G','T'} ]
    - g         - gap string [ '-' ]
    - t         - to remove columns with more gaps than this frequency threshold
    """
    # numpy array form of the MSA
    ar = numpy.array( [ list(rec) for rec in aln]) 
    # delete columns if gap frequency exceeds threshold
    ar = rmc_ef(ar,e=g,f=t)
    # number of sequences and number of sites
    N,L = ar.shape
        
    # no. of pairs
    numpairs = special.binom(N, 2)
    
    # apply the above function to all positions
    if L > 0:
        totd = sum([ pwd_ar1d(ar[:,k],alpha) for k in range(L)])    
        #print("total difference = {}, # pairs = {}, # sites = {}".format(totd,numpairs,L))
        div = round( totd/(numpairs*L), 4)
        return div

def nucdiv_inv(msa,inv):
    """
    Return nucleotide diversity over a range of columns in the MSA.
    
    Arguments:
    - msa - biopython MSA object
    - inv - tuple of start and end indices of the range
    """
    return ( inv[0], inv[1], nucdiv_aln( msa[ :, inv[0]:inv[1] ]))

def nucdiv_slide(msa,window,jump,ncpu):
    """Compute nucleotide diversity in a sliding window over the genome"""
    # no. of columns
    L = msa.get_alignment_length()
    # list of start and end positions with the above 
    # window and jump lengths
    starts = range( 0, L, jump)
    intervals = [ (i, i + window) for i in starts if i + window < L]
    if len(intervals) == 0:
        return

    if ncpu > 1:
        with multiprocessing.Pool(ncpu) as mp_pool:
            divs = mp_pool.starmap(  nucdiv_inv, [ (msa,i) for i in intervals])
    else:
        divs = [ nucdiv_inv(msa,i) for i in intervals]    
    return divs
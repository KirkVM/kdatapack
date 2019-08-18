from collections import Counter
import pandas as pd
import numpy as np
import random
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#builds cdf (codon data frame) using codonfreqtable in this directory on import
#this is an E. coli codon table
cdf=pd.read_excel('codonfreqtable.xlsx')
cdf['GCcontent']=[Counter(cdf['codon'].at[x])['G']+Counter(cdf['codon'].at[x])['C'] for x in cdf.index]
cdf.assign(cumulative_freq_min=0,cumulative_freq_max=0)
for groups in cdf.groupby('aa'):
    fsum=0
    for cdnidx in groups[1].index:
        cdf.loc[cdnidx,'cumulative_freq_min']=fsum
        fsum+=groups[1].loc[cdnidx,'freq']
        cdf.loc[cdnidx,'cumulative_freq_max']=fsum


#def get_new_codon(cdnseq):
#    cdnrow=cdf[cdf['codon']==cdnseq]
#    cdnfreq=cdnrow['freq'].iloc[0]
#    cdnaa=cdnrow['aa'].iloc[0]
#    cdndf=cdf[cdf['aa']==cdnaa]
#    for x in cdndf.index:
#        if cdndf.GCcontent.at[x]<cdnrow.GCcontent.iloc[0] \
#            and cdndf.freq.at[x]>cdnrow.freq.iloc[0]-0.05:
#            newcodon=cdndf.codon.at[x]
#            return newcodon
#    return cdnseq

def get_ratiod_codon(aa):
    '''like get_opt_codon but with aa as argument'''
    cdndf=cdf[cdf['aa']==aa]
    randval=random.random()
    for cdnidx in cdndf.index:
        if randval>=cdndf.loc[cdnidx,'cumulative_freq_min'] and\
            randval<=cdndf.loc[cdnidx,'cumulative_freq_max']:
            return cdndf.loc[cdnidx,'codon']

def get_opt_codon(cdnseqstr):
    ''' returns a new codon with probability according to freq table

    Arguments:
        cdnseq: sequence of starting codon
    Returns:
        newcdnseq: str sequence of new codon
    '''
    cdnrow=cdf[cdf['codon']==cdnseqstr]
    cdnfreq=cdnrow['freq'].iloc[0]
    cdnaa=cdnrow['aa'].iloc[0]

    cdndf=cdf[cdf['aa']==cdnaa]
    randval=random.random()
    for cdnidx in cdndf.index:
        if randval>=cdndf.loc[cdnidx,'cumulative_freq_min'] and\
            randval<=cdndf.loc[cdnidx,'cumulative_freq_max']:
            return cdndf.loc[cdnidx,'codon']
    return cdnseqstr
   
def get_reducedgc_codon(cdnseq):
    ''' returns a new codon with lower GC content iff lower GC AND codon freq not <0.05 than orig

    Arguments:
        cdnseq: sequence of starting codon
    Returns:
        newcdnseq: str sequence of new codon (possibly the starting codon...)
    '''
    cdnrow=cdf[cdf['codon']==cdnseq]
    cdnfreq=cdnrow['freq'].iloc[0]
    cdnaa=cdnrow['aa'].iloc[0]
    cdndf=cdf[cdf['aa']==cdnaa]
    for cdnidx in cdndf.index:
        if cdndf.GCcontent.at[cdnidx]<cdnrow.GCcontent.iloc[0] \
            and cdndf.freq.at[cdnidx]>cdnrow.freq.iloc[0]-0.05:
            newcodon=cdndf.codon.at[cdnidx]
            return newcodon
    return cdnseq

def optimize_codon_frequencies(dnaseq):
    ''' _optimizes_ (via probability) a protein-coding DNA sequence

    Arguments:
        dnaseq: biopython DNA sequence
    Returns:
        newdnaseq: biopython DNA sequence with codon frequencies that correspond to table probabilities
    '''
    newdnastr=''
    for codonposition in range(0,len(dnaseq),3):
        cdnseq=dnaseq[codonposition:codonposition+3]
        newcdnstr=get_opt_codon(str(cdnseq))
        newdnastr+=newcdnstr
    newdnaseq=Seq(newdnastr,alphabet=generic_dna)
    return newdnaseq
 
def reduce_gc_frequency(dnaseq):
    ''' reduces GC-content in a protein-coding DNA sequence

    Arguments:
        dnaseq: biopython DNA sequence
    Returns:
        newdnaseq: biopython DNA sequence with GC content reduced as defined
                    in get_reducedgc_codon
    '''
 
    newdnastr=''
    for codonposition in range(0,len(dnaseq),3):
        cdnseq=dnaseq[codonposition:codonposition+3]
        newcdnstr=get_reducedgc_codon(str(cdnseq))
        newdnastr+=newcdnstr
    newdnaseq=Seq(newdnastr,alphabet=generic_dna)
    return newdnaseq    


def get_avg_freq(dnaseq):
    freqs=[]
    for cdnpos in range(0,len(dnaseq),3):
        cdnseq=str(dnaseq[cdnpos:cdnpos+3])
        cdnrow=cdf[cdf['codon']==cdnseq]
        cdnfreq=cdnrow['freq'].iloc[0]
        freqs.append(cdnfreq)
    return(np.mean(freqs))
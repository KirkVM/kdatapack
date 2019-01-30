import os,sys,shutil
from typing import List
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
import itertools




def get_aligndf(srs:List[SeqRecord],aligner:str='clustal',alnfpath:str=None):
    """Builds dataframe with index & column values set to acc number of each member

    Arguments:
    srs: list of sequence records 

    Keyword Arguments:
    aligner: aligner type ('clustal') (default 'clustal')
    alnfpath: path to store alignment file (default None)

    Returns:
    pandas data frame with index = acc.number & column = acc. number. data type is float
    """
    tempseqfpath="temp_seqfile"
    temppwfpath="temp_pwfile"
    tempalnfpath='temp_align'
    SeqIO.write(srs,tempseqfpath,'fasta')
    #pw-vals I'm getting are smaller than via ebi web interface by ~10%, why? (2/24/17)
    #maybe just weighting?
    #4/17/17-diffs scale with pwdiff- ie are very close for top-similarities=95% vs 94.5% eg...but differ by 10% or more as you get to 30% & below
    cocl=ClustalOmegaCommandline(infile=tempseqfpath,outfile=tempalnfpath,distmat_out=temppwfpath,distmat_full=True,force=True,percentid=True,distmat_full_iter=True,outfmt="st")
    #cocl=ClustalOmegaCommandline(infile=tempseqfpath,outfile=tempalnfpath,iter=0,distmat_out=temppwfpath)#,distmat_full=True,force=True,percentid=True,distmat_full_iter=True,outfmt="st")
    cocl()
    pwfile=open(temppwfpath,'r')
    pwvals_=[]
    for x,l in enumerate(pwfile.readlines()):
        entries_=[x.strip() for x in l.split()[1:]]
        pwvals_.extend(entries_[x+1:])
    #shutil.copy(tempseqfpath,"seqf_{0}.fasta".format(self.klustername))
    os.remove(tempseqfpath)
    pwfile.close()
    if alnfpath is not None:
        shutil.copy(tempalnfpath,alnfpath)
    os.remove(tempalnfpath)
    #pwfile=open(temppwfpath,'r')
    pwdf=pd.read_csv(temppwfpath,sep='\s+',index_col=0,skiprows=1,header=None)
    os.remove(temppwfpath)
    pwdf.columns=list(pwdf.index)
    return pwdf


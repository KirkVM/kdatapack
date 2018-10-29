import kazykools# import pullfamily
import gbkools# import pullgb_fromcazyobjs
from Bio import SeqIO

family="GH43"
czes_=kazykools.pullfamily(family)
aHT=gbkools.pullgb_fromcazyobjs(czes_)
for subfam in aHT.keys():
    SeqIO.write(aHT[subfam],"{0}_{1}seqs.fasta".format(family,subfam),"fasta")

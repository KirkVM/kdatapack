import os,time
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def getgbseqrec(gbid):
	handle=Entrez.efetch(db="protein",id=gbid,rettype="fasta",retmode="xml")
	a=Entrez.read(handle)
	#print a[0]['TSeq_orgname']
	#sr=SeqRecord(Seq(a[0]['TSeq_sequence'],IUPAC.protein),id=sid,name=a[0]['TSeq_defline'])
	sr=SeqRecord(Seq(a[0]['TSeq_sequence'],IUPAC.protein),id=gbid,description=a[0]['TSeq_defline'])
	return sr

#def pullgb_fromidlist(ids_,ofname="cool"):
#	fastasrs_=[]
#	Entrez.email="kirk.vandermeulen@wisc.edu"
#	for x,sid in enumerate(accids_):
#		#if x<6500:continue
#		fastasrs_.append(getseqrec(sid))
#		if x!=0 and x%100==0:
#			print 'accumulated {0} sequences'.format(len(fastasrs_))
#			time.sleep(2)
#		if x!=0 and x%500==0:
#			SeqIO.write(fastasrs_,"{0}{1}new.fasta".format(,ofname,len(fastasrs_)),"fasta")
#			time.sleep(60)
#
#	SeqIO.write(fastasrs_,"{0}.fasta".format(ofname),"fasta")

def pullgb_fromcazyobjs(czes_):
    Entrez.email="kirk.vandermeulen@wisc.edu"
    fastasrs_=[]
    outfamHT={}
    outfamHT['all']=[]
    for x,cze in enumerate(czes_):
        if x%10==0: 
            time.sleep(2)
        if x%25==0: 
            print('grabbed {0} seqs from GB'.format(x))
            time.sleep(1)
            if x%50==0:
                time.sleep(2)
                if x%100==0:
                    time.sleep(10)
        try:
            gbid=cze.gbids_[0]
            print(gbid)
#          if len(cze.gbids_)>1:
#               alt_gbid=cze.gbids[1]
            sr=getgbseqrec(gbid)
            #now we add add'l stuff to id or description
#           full_desc=sr.description
            dbxrefs_=[]
            if len(cze.gbids_)>1:
                #full_desc+="|"+cze.gbids_[1]
                extragbs_=cze.gbids_[1:]
                for gbid in extragbs_:
                    dbxrefs_.append("GB:{0}".format(gbid))
            if len(cze.ecs_)>0:
                for ec in cze.ecs_:
                    dbxrefs_.append("EC:{0}".format(ec))
            if len(cze.pdbids_)>0:
                for pdbid in cze.pdbids_:
                    #full_desc+="|"+pdbid
                    dbxrefs_.append("PDB:{0}".format(pdbid))
            if cze.family!=None:
                #full_desc+="|SUBFAM"+cze.family
                dbxrefs_.append("SUBFAM:{0}".format(cze.family))
#           full_desc=sr.description
            full_desc=sr.description
            for dbxref in dbxrefs_:
                full_desc+="|"+dbxref
            newsr=SeqRecord(sr.seq,id=sr.id,description=full_desc)#sr.description,dbxrefs=dbxrefs_)
            #fastasrs_.append(sr)
            if cze.family is not None:
                if not int(cze.family) in outfamHT.keys(): outfamHT[int(cze.family)]=[]
                outfamHT[int(cze.family)].append(newsr)
            outfamHT['all'].append(newsr)
        except:
            print('could not add {}'.format(cze.gbids_[0]))
    return outfamHT

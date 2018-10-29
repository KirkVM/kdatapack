import operator,pdb,re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from operator import attrgetter

class dMatch:
    def __init__(self,theline):
        self.gh_redict={'name':['[Gg]lyco_hydro','GH\d+','Glyco_hydr'],'desc':['[Gg]lycosyl hydrolase','[Gg]lycosyl-hydrolase'] }
        self.cbm_redict={'name':['CBM']}
        self.maccHT={'GH5':'PF00150','CBM_X2':'PF03442'}
        entries=list(map(lambda x:x.strip(),theline.split()))
        self.dname=entries[0]
        self.motifacc=(entries[1])
        self.dsize=int(entries[2])
        self.accession=entries[3]
        self.psize=int(entries[5])
        self.fevalue=float(entries[6])
        self.fscore=float(entries[7])
        self.fbias=float(entries[8])
        self.matchnum=int(entries[9])
        self.matchtot=int(entries[10])
        self.dcevalue=float(entries[11])
        self.dievalue=float(entries[12])
        self.dscore=float(entries[13])
        self.dbias=float(entries[14])
        self.hmmstart=int(entries[15])
        self.hmmstop=int(entries[16])
        self.alstart=int(entries[17])
        self.alstop=int(entries[18])
        self.envstart=int(entries[19])
        self.envstop=int(entries[20])
        self.acc=float(entries[21])
        self.motifdesc=''
        for e in entries[22:]:
            self.motifdesc+='{} '.format(e)#entries[22:]
        self.motifdesc=self.motifdesc[:-1]
        self.picdstart=None
        self.picdstop=None
        self.modded_start=False #has picdstart position been modified to accommodate adjacent domain?
        self.modded_stop=False #has picdstop position been modified to accommodate adjacent domain?
        self.motifclass=None
        for namere in self.gh_redict['name']:
            if re.search(namere,self.dname):
                self.motifclass="GH"
        for descre in self.gh_redict['desc']:
            if re.search(descre,self.motifdesc):
                self.motifclass="GH"
        for namere in self.cbm_redict['name']:
            if re.search(namere,self.dname):
                self.motifclass="CBM"
#        for descre in self.cbm_redict['desc']:
#            if re.search(descre,self.motifdesc):
#                self.motifclass="CBM"

                


def parse_hmmfile(domtbloutfpath):
    """returns a dictionary with a list for each domain annoation:
    [d_name,hmm_domain_start,hmm_domain_end,protein_domain_start,protein_domain_end,hmm_domain_size,coverage,cvalue]
    lists are sorted according to protein_domain_start"""
    dbcanfile=open(domtbloutfpath,'r')
    dHT={}
    for line in dbcanfile.readlines():
        if line[0]=="#":continue
        hmmermatch=dMatch(line)
        if not hmmermatch.accession in dHT.keys():dHT[hmmermatch.accession]=[]
        dHT[hmmermatch.accession].append(hmmermatch)
    #now sort each by dievalue:
    for acc in dHT:
        dHT[acc]=sorted(dHT[acc],key=lambda x:x.dievalue)
    return dHT

def parse_hmmfile2(domtbloutfpath):
    """returns a dictionary with a list for each domain annoation:
    [d_name,hmm_domain_start,hmm_domain_end,protein_domain_start,protein_domain_end,hmm_domain_size,coverage,cvalue]
    lists are sorted according to protein_domain_start"""
    dbcanfile=open(domtbloutfpath,'r')
    dHT={}
    for line in dbcanfile.readlines():
        if line[0]=="#":continue
        hmmermatch=dMatch(line)
        if not hmmermatch.accession in dHT.keys():dHT[hmmermatch.accession]=[]
        dHT[hmmermatch.accession].append(hmmermatch)
    #now sort each by dievalue:
    return dHT

"""
def cleanannotations2(dmatches_):
    sdmatches_=sorted(dmatches_,key=attrgetter('dievalue'),reverse=False)
    topmatches_=[]
    for dmatch in sdmatches_:
        #print(dmatch.accession,dmatch.dievalue)
        dmatch.picdstart=max(1,dmatch.alstart-dmatch.hmmstart+1)
        dmatch.picdstop=min(dmatch.psize,dmatch.alstop+(dmatch.dsize-dmatch.hmmstop))
        
        addable=True
        modded=False 
        tmindex=0
        for tmindex in range(len(topmatches_)-1,-1,-1):
            if not addable:
                break
            topmatch=topmatches_[tmindex]

            curmacc=dmatch.motifacc.split('.')[0]
            curstart=dmatch.picdstart
            curstop=dmatch.picdstop
            curmotifclass=dmatch.motifclass
            prevmacc=topmatch.motifacc.split('.')[0]
            prevstart=topmatch.picdstart
            prevstop=topmatch.picdstop
            prevmotifclass=topmatch.motifclass
            
            
            
            mergestatus=None
            if (curstart>=prevstart and curstop<prevstop) or (curstart>prevstart and curstop<=prevstop):
                mergestatus="within"
            elif curstart<=prevstart and curstop>=prevstop:
                mergestatus="span"
            elif curstart<=prevstop and curstop>prevstop:
                mergestatus="ohright"
            elif curstart<prevstart and curstop>=prevstart:
                mergestatus="ohleft"
            elif curstart>prevstop or curstop<prevstart:
                mergestatus="clean"
            else:
                exit("this shouldn't happen!")
            
            if mergestatus in ['within','span']:
                addable=False
            elif mergestatus=="ohright":
                overlapsize=prevstop-curstart+1
                freesize=curstop-prevstop

                if curmacc=="PF03442" and prevmacc=="PF00150":
                    if overlapsize<20:
                        dmatch.picdstart=prevstop+1
                        dmatch.modded_start=True
                    else:
                        addable=False
                #trying to add CBM to right of GH
                elif curmotifclass=="CBM" and prevmotifclass=="GH":
                    if overlapsize<10:
                        dmatch.picdstart=prevstop+1
                        dmatch.modded_start=True
                    else:
                        addable=False
                 #trying to add GH to right of CBM
                elif curmotifclass=="GH" and prevmotifclass=="CBM":
                    if overlapsize<10:
                        topmatch.picdstop=curstart-1
                        topmatch.modded_stop=True
                    else:
                        addable=False

            elif mergestatus=="ohleft":
                overlapsize=curstop-prevstart+1
                freesize=prevstart-curstart
                #trying to add GH5 to left of CBM46
                if curmacc=="PF000150" and prevmacc=="PF03442":
                    if overlapsize<20:
                        topmatch.picdstart=curstop+1
                        topmatch.modded_start=True
                    else:
                        addable=False
                #trying to add GH to left of CBM
                elif curmotifclass=="GH" and prevmotifclass=="CBM":
                    if overlapsize<10:
                        topmatch.picdstart=curstop+1
                        topmatch.modded_start=True
                    else:
                        addable=False
                #trying to add CBM to left of GH
                elif curmotifclass=="CBM" and prevmotifclass=="GH":
                    if overlapsize<10:
                        dmatch.picdstop=prevstart-1
                        dmatch.modded_stop=True
                    else:
                        addable=False
 
                       
            #now, if overlap is <4 or <4% for both domains AND not already modded
            #AND good e-values,truncate the larger domain and allow
            if (mergestatus=='ohleft' and dmatch.modded_stop==False and topmatch.modded_start==False) or\
                (mergestatus=='ohright' and dmatch.modded_start==False and topmatch.modded_stop==False):
                if topmatch.dievalue<0.01 and dmatch.dievalue<0.01:
                    #sort the domains so longer is first

                    doms_=sorted([topmatch,dmatch],key=attrgetter('dsize'),reverse=True)
                    if (overlapsize<4) or (float(overlapsize)/doms_[1].dsize<0.04):
                        if doms_[0].picdstart<doms_[1].picdstart: #long one is first
                            doms_[0].picdstop=doms_[1].picdstart-1
                            doms_[0].modded_stop=True
                        else:
                            doms_[0].picdstart=doms_[1].picdstop+1
                            doms_[0].modded_start=True
                    else:
                        addable=False
                else:
                    addable=False
            tmindex+=1
            #now see if necessary truncations make it implausible:
            if (float(dmatch.picdstop-dmatch.picdstart))/(float(dmatch.dsize))<0.5:
                addable=False
        if addable:
            topmatches_.append(dmatch)
    return topmatches_
"""

def cleanannotations_simple(dmatches_):
    sdmatches_=sorted(dmatches_,key=attrgetter('dievalue'),reverse=False)
    topmatches_=[]
    for dmatch in sdmatches_:
        #print(dmatch.accession,dmatch.dievalue)
        dmatch.picdstart=max(1,dmatch.alstart-dmatch.hmmstart+1)
        dmatch.picdstop=min(dmatch.psize,dmatch.alstop+(dmatch.dsize-dmatch.hmmstop))
        
        addable=True
        modded=False 
        tmindex=0
        for tmindex in range(len(topmatches_)-1,-1,-1):
            if not addable:
                break
            topmatch=topmatches_[tmindex]

            curmacc=dmatch.motifacc.split('.')[0]
            curstart=dmatch.picdstart
            curstop=dmatch.picdstop
            curmotifclass=dmatch.motifclass
            prevmacc=topmatch.motifacc.split('.')[0]
            prevstart=topmatch.picdstart
            prevstop=topmatch.picdstop
            prevmotifclass=topmatch.motifclass

            
            mergestatus=None
            if (curstart>=prevstart and curstop<prevstop) or (curstart>prevstart and curstop<=prevstop):
                mergestatus="within"
            elif curstart<=prevstart and curstop>=prevstop:
                mergestatus="span"
            elif curstart<=prevstop and curstop>prevstop:
                mergestatus="ohright"
            elif curstart<prevstart and curstop>=prevstart:
                mergestatus="ohleft"
            elif curstart>prevstop or curstop<prevstart:
                mergestatus="clean"
            else:
                exit("this shouldn't happen!")
            
            if mergestatus in ['within','span']:
                addable=False
            elif mergestatus=="ohright":
                overlapsize=prevstop-curstart+1
                freesize=curstop-prevstop
#                
#                if curmacc=="PF03442" and prevmacc=="PF00150":
#                    if overlapsize<20:
#                        dmatch.picdstart=prevstop+1
#                        dmatch.modded_start=True
#                    else:
#                        addable=False
#                #trying to add CBM to right of GH
#                elif curmotifclass=="CBM" and prevmotifclass=="GH":
#                    if overlapsize<10:
#                        dmatch.picdstart=prevstop+1
#                        dmatch.modded_start=True
#                    else:
#                        addable=False
#                 #trying to add GH to right of CBM
#                elif curmotifclass=="GH" and prevmotifclass=="CBM":
#                    if overlapsize<10:
#                        topmatch.picdstop=curstart-1
#                        topmatch.modded_stop=True
#                    else:
#                        addable=False

            elif mergestatus=="ohleft":
                overlapsize=curstop-prevstart+1
                freesize=prevstart-curstart
#                #trying to add GH5 to left of CBM46
#                if curmacc=="PF000150" and prevmacc=="PF03442":
#                    if overlapsize<20:
#                        topmatch.picdstart=curstop+1
#                        topmatch.modded_start=True
#                    else:
#                        addable=False
#                #trying to add GH to left of CBM
#                elif curmotifclass=="GH" and prevmotifclass=="CBM":
#                    if overlapsize<10:
#                        topmatch.picdstart=curstop+1
#                        topmatch.modded_start=True
#                    else:
#                        addable=False
#                #trying to add CBM to left of GH
#                elif curmotifclass=="CBM" and prevmotifclass=="GH":
#                    if overlapsize<10:
#                        dmatch.picdstop=prevstart-1
#                        dmatch.modded_stop=True
#                    else:
#                        addable=False
 
            if (mergestatus=="ohleft" or mergestatus=="ohright"):
                if overlapsize>=0:# or (float(overlapsize)/doms_[1].dsize<0.04):
                    addable=False
        if addable:
            topmatches_.append(dmatch)
                       
#            #now, if overlap is <4 or <4% for both domains AND not already modded
#            #AND good e-values,truncate the larger domain and allow
#            if (mergestatus=='ohleft' and dmatch.modded_stop==False and topmatch.modded_start==False) or\
#                (mergestatus=='ohright' and dmatch.modded_start==False and topmatch.modded_stop==False):
#                if topmatch.dievalue<0.01 and dmatch.dievalue<0.01:
#                    #sort the domains so longer is first
#
#                    doms_=sorted([topmatch,dmatch],key=attrgetter('dsize'),reverse=True)
#                    if (overlapsize<4) or (float(overlapsize)/doms_[1].dsize<0.04):
#                        if doms_[0].picdstart<doms_[1].picdstart: #long one is first
#                            doms_[0].picdstop=doms_[1].picdstart-1
#                            doms_[0].modded_stop=True
#                        else:
#                            doms_[0].picdstart=doms_[1].picdstop+1
#                            doms_[0].modded_start=True
#                    else:
#                        addable=False
#                else:
#                    addable=False
#            tmindex+=1
#            #now see if necessary truncations make it implausible:
#            if (float(dmatch.picdstop-dmatch.picdstart))/(float(dmatch.dsize))<0.5:
#                addable=False
#        if addable:
#            topmatches_.append(dmatch)
    return topmatches_

#returns a list of domains for each accession
def getdict_fromdomtbloutfile(domtbloutfpath):
    domHT=parse_hmmfile2(domtbloutfpath)
    for acc in domHT:
        domHT[acc]=cleanannotations_simple(domHT[acc])
#    for dm in domHT['AEV59735.1']+domHT['AFH62788.1']:
#        print(dm.accession,dm.dname,dm.motifacc,dm.picdstart,dm.picdstop,dm.motifdesc)
    return domHT


import re,os
from operator import itemgetter,attrgetter
import pdb
qRE=re.compile("Query:\s+(\S+)\s+\[M=(\d+)\].*")

accRE=re.compile("([a-z]+)\|([0-9a-zA-Z_\.]+)\|")
smpl_accRE=re.compile("([0-9A-Z_\.]+)")

fltstrgrp="([0-9e\.\-\+]+)"
smrestr="\s*"
for x in range(8):
    smrestr+=fltstrgrp+"\s+"
smrestr+="(\S+)\s+(\S+)"
smRE=re.compile(smrestr)

dmrestr="\s*(\d+)\s+(\!|\?)\s+"
for x in range(14):
    dmrestr+="\s+(\S+)"
dmRE=re.compile(dmrestr)



def flattenuids(accHT):
    valstr=''
    for k in sorted(accHT.keys()):
        valstr+=accHT[k]+'|'
    return valstr

class DomProfileMatch:
    def __init__(self,l,vals_):
        self.dpm_num=int(vals_[0])
        self.thresh=vals_[1]
        self.score=float(vals_[2])
        self.bias=float(vals_[3])
        self.c_evalue=float(vals_[4])
        self.i_evalue=float(vals_[5])
        self.hmm_start=int(vals_[6])
        self.hmm_stop=int(vals_[7])
        self.hmm_start_state=vals_[8][0]
        self.hmm_stop_state=vals_[8][1]
        self.align_start=int(vals_[9])
        self.align_stop=int(vals_[10])
        self.align_start_state=vals_[11][0]
        self.align_stop_state=vals_[11][1]
        self.env_start=int(vals_[12])
        self.env_stop=int(vals_[13])
        self.env_start_state=vals_[14][0]
        self.env_stop_state=vals_[14][1]
        self.acc=float(vals_[15])
        self.textline=l

class SeqProfileMatch:
    def __init__(self,accHT,pname,pacc,psize,puv=None):#gi=None,ref=None,gb=None,emb=None):
        self.dpms_=[]
        self.accHT=accHT
        self.pname=pname
        self.pacc=pacc
        self.psize=int(psize)
        self.dpm_hdrlines_=None
        if self.accHT is not None:
            self.prot_uval=flattenuids(self.accHT)
        elif puv is not None:
            self.prot_uval=puv
    def addspmvals(self,line,vals_):
        self.evalue=float(vals_[0])
        self.score=float(vals_[1])
        self.bias=float(vals_[2])
        self.bdevalu=float(vals_[3])
        self.bdscore=float(vals_[4])
        self.bdbias=float(vals_[5])
        self.ndomexp=float(vals_[6])
        self.ndomn=int(vals_[7])
        self.textline=line
    def add_dpm(self,hdrlines_,line,vals_):
        self.dpms_.append(DomProfileMatch(line,vals_))
        self.dpms_.sort(key=attrgetter('i_evalue'))
        if self.dpm_hdrlines_ is None:
            self.dpm_hdrlines_=hdrlines_[:]


class HMMERSearchMotif:
    def __init__(self,lines_,pname,psize,pacc,acc_mode="simple"):
        self.lines_=lines_
        self.pname=pname
        self.psize=psize
        self.pacc=None
        self.pacc_version=None
        #dbcan does not have acc_val. 
        if pacc is not None: #pfam. split pacc into label & version
            pfamaccRE=re.compile('(PF\d{5})\.(\d+)')
            pfam_accobj=pfamaccRE.match(pacc)
            if pfam_accobj:
                self.pacc=pfam_accobj.group(1)
                self.pacc_version=pfam_accobj.group(2)
        else: #dbcan/others. set pacc=pname
            self.pacc=pname
        self.hdrlines_=[]
        self.ftrlines_=[]
        self.spms_=[]
        self.readin(acc_mode=acc_mode)
#        print(self.hdrlines_)
#        print(self.pname,len(self.spms_))
#        print(self.ftrlines_)
    def readin(self,acc_mode="simple"):
#        print('reading in')
        accRE=re.compile("([a-z]+)\|([0-9a-zA-Z_\.]+)\|")
        smpl_accRE=re.compile("([0-9A-Z_\.]+)")
        fltstrgrp="([0-9e\.\-\+]+)"
#        smrestr="\s*"
        smrestr="\s*(\d[0-9e\.\-\+]*)\s+"
        for x in range(7):
            smrestr+=fltstrgrp+"\s+"
        smrestr+="(\S+)\s+(\S+)"
        smRE=re.compile(smrestr)
        dmrestr="\s*(\d+)\s+(\!|\?)\s+"
        for x in range(14):
            dmrestr+="\s+(\S+)"
        dmRE=re.compile(dmrestr)
        nhRE=re.compile("\[No hits detected that satisfy reporting thresholds\]")
        in_hdr=True 
        in_domains=False
        in_summary=False
        cur_spm=None
        for x,l in enumerate(self.lines_):
#            if self.pacc=="GH5_7.hmm":
#                pdb.set_trace()
#            print(self.lines_)
            if l[:8]=="Internal": #this means we've passed everything useful
                in_summary=True
                in_domains=False
            entries_=[x.strip() for x in l.split(' ')]
            smobj=smRE.match(l)
            if smobj:
                in_hdr=False
                if acc_mode=="simple":
                    accobj=smpl_accRE.search(smobj.group(9))
                    if accobj:
                        prot_uval=accobj.group(1)
                        spm=SeqProfileMatch(None,self.pname,self.pacc,self.psize,puv=prot_uval)
                        spm.addspmvals(l,smobj.groups()[0:8])
                        self.spms_.append(spm)
                else:
                    accobjs_=accRE.findall(smobj.group(9))
                    accHT={}
                    if len(accobjs_)>0:
                        for accgrp in accobjs_:
                            accHT[accgrp[0]]=accgrp[1]
                        spm=SeqProfileMatch(accHT,self.pname,self.pacc,self.psize)
                        spm.addspmvals(l,smobj.groups()[0:8])
                        self.spms_.append(spm)
            if nhRE.search(l):
                in_hdr=False
            if l[:2]==">>":
                hdrlines_=['\n']
                in_domains=True
                cur_spm=None #reset each time
                dpm_accHT={}
                if acc_mode=="simple":
                    accobj=smpl_accRE.search(l)
                    prot_uval=accobj.group(1)
                    for spm in self.spms_:
                        if spm.prot_uval==prot_uval:
                            cur_spm=spm
                else:
                    accobjs_=accRE.findall(l)
                    for accgrp in accobjs_:
                        dpm_accHT[accgrp[0]]=accgrp[1]
                        prot_uval=flattenuids(dpm_accHT)
                        for spm in self.spms_:
                            if spm.prot_uval==prot_uval:
                                cur_spm=spm
            if in_domains:
                dmobj=dmRE.match(l)
                if dmobj:
                    cur_spm.add_dpm(hdrlines_,l,dmobj.groups())
                else:
                    hdrlines_.append(l)
            if in_hdr:
                self.hdrlines_.append(l)
                #if self.pacc=="GH5_7.hmm":pdb.set_trace()
            if in_summary:
                self.ftrlines_.append(l)
    def get_text(self):
        spmlines_=[]
        dpmlines_=[]
#        if self.pacc=="GH5_7.hmm":pdb.set_trace()
        for spm in self.spms_:
            if len(spm.dpms_)>0:
                spmlines_.append(spm.textline)
                for x,dpm in enumerate(spm.dpms_):
                    if x==0:
                        dpmlines_.extend(spm.dpm_hdrlines_)
                    dpmlines_.append(dpm.textline)
        rtext=""
        for l in self.hdrlines_:
            rtext+=l
        if len(spmlines_)>0: 
            for l in spmlines_:
                rtext+=l
        else:
            rtext+="\n   [No hits detected that satisfy reporting thresholds]\n\n"


        rtext+='\nDomain Annotation for each sequence:'
        if len(dpmlines_)>0:
            for l in dpmlines_:
                rtext+=l
        else:
            rtext+="\n\n   [No targets detected that satisfy reporting thresholds]"
        rtext+='\n\n\n'
        for l in self.ftrlines_:
            rtext+=l
        return rtext
    
class HMMERSearchResParser:
    def __init__(self):#,nrfpath):
        self.motifresults_=[]
        self.spms_=None
        self.motifHT=None
        self.protHT=None
    def file_read(self,nrfpath):
        nrfile=open(nrfpath,'r')
        prof_ghname=None
        prof_size=None
        prof_acc=None
        qRE=re.compile("Query:\s+(\S+)\s+\[M=(\d+)\].*")
        paccRE=re.compile("Accession:\s+(\S+)")
        pls_=[]
        for x,l in enumerate(nrfile.readlines()):
            qobj=qRE.match(l)
            if qobj:
                prof_ghname=qobj.group(1)
                prof_size=qobj.group(2)
                pls_=[l]
            else:
                pls_.append(l)
            paccobj=paccRE.match(l)
            if paccobj:
                prof_acc=paccobj.group(1)
            if l[:2]=='//':
#                print(prof_ghname,prof_size,prof_acc,len(pls_))
#                print(pls_)
                curmotifobj=HMMERSearchMotif(pls_,prof_ghname,prof_size,prof_acc)
                self.motifresults_.append(curmotifobj)
                pls_=[]
                prof_ghname=None
                prof_size=None
                prof_acc=None
#            print(x,l)
        nrfile.close()

        for x in range(len(self.motifresults_)-1,-1,-1):
            for y in range(x-1,-1,-1):
                if self.motifresults_[x].pacc==self.motifresults_[y].pacc:
                    self.motifresults_[y].spms_.extend(self.motifresults_[x].spms_)#.append(spm)
                    self.motifresults_.pop(x)
                    break


    def get_spm_list(self):
        self.spms_=[]
        for mr in self.motifresults_:
            for spm in mr.spms_:
                self.spms_.append(spm)
        return self.spms_

    def get_motifHT(self):
        self.motifHT={}
        for mr in self.motifresults_:
            #if mr.pacc not in self.motifHT.keys():#merge during read so not needed
            self.motifHT[mr.pacc]=[]
            for spm in mr.spms_:
                self.motifHT[mr.pacc].append(spm)
        return self.motifHT
    
    def get_protHT(self):
        self.protHT={}
        #if self.spms_ is None: 
        self.get_spm_list()
        for spm in self.spms_:
            curp_uval=spm.prot_uval
            if curp_uval not in self.protHT.keys():
                self.protHT[curp_uval]=[]
            self.protHT[curp_uval].append(spm)
        return self.protHT 

    def write_searchres(self,ofpath,pmotifaccs_):
        ofile=open(ofpath,'w')
        ofile.write('#kirk re-write\n')
        for motifresult in self.motifresults_:
            if motifresult.pacc in pmotifaccs_:
                motiftext=motifresult.get_text()
                ofile.write(motiftext)
        ofile.close()


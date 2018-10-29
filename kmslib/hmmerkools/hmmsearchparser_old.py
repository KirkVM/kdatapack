import re,os
from operator import itemgetter,attrgetter

qRE=re.compile("Query:\s+(\S+)\s+\[M=(\d+)\].*")
paccRE=re.compile("Accession:\s+(\S+)")

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
    def __init__(self,vals_):
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

class SeqProfileMatch:
    def __init__(self,accHT,pname,pacc,psize,puv=None):#gi=None,ref=None,gb=None,emb=None):
        self.dpms_=[]
        self.accHT=accHT
        self.pname=pname
        self.pacc=pacc
        self.psize=int(psize)
        if self.accHT is not None:
            self.prot_uval=flattenuids(self.accHT)
        elif puv is not None:
            self.prot_uval=puv
    def addspmvals(self,vals_):
        self.evalue=float(vals_[0])
        self.score=float(vals_[1])
        self.bias=float(vals_[2])
        self.bdevalu=float(vals_[3])
        self.bdscore=float(vals_[4])
        self.bdbias=float(vals_[5])
        self.ndomexp=float(vals_[6])
        self.ndomn=int(vals_[7])
    def add_dpm(self,vals_):
        self.dpms_.append(DomProfileMatch(vals_))
        self.dpms_.sort(key=attrgetter('i_evalue'))



def parse_profmatches(lines_,pname,pacc,psize,acc_mode="simple"):
    in_domains=False
    in_summary=False
    cur_spm=None
    spms_=[]
    for l in lines_:
        if l[:9]=="Internal":
            in_summary=True
            in_domains=False
        
        entries_=[x.strip() for x in l.split(' ')]
        smobj=smRE.match(l)
        if smobj:
            if acc_mode=="simple":
                accobj=smpl_accRE.search(smobj.group(9))
                if accobj:
                    prot_uval=accobj.group(1)
                    spm=SeqProfileMatch(None,pname,pacc,psize,puv=prot_uval)
                    spm.addspmvals(smobj.groups()[0:8])
                    spms_.append(spm)
            else:
                accobjs_=accRE.findall(smobj.group(9))
                accHT={}
                if len(accobjs_)>0:
                    for accgrp in accobjs_:
                        accHT[accgrp[0]]=accgrp[1]
                    spm=SeqProfileMatch(accHT,pname,pacc,psize)
                    spm.addspmvals(smobj.groups()[0:8])
                    spms_.append(spm)
        if l[:2]==">>":
            #print(l)
            in_domains=True
            cur_spm=None #reset each time
            dpm_accHT={}
            if acc_mode=="simple":
                accobj=smpl_accRE.search(l)
                prot_uval=accobj.group(1)
                for spm in spms_:
                    if spm.prot_uval==prot_uval:
                        cur_spm=spm
            else:
                accobjs_=accRE.findall(l)
                for accgrp in accobjs_:
                    dpm_accHT[accgrp[0]]=accgrp[1]
                    prot_uval=flattenuids(dpm_accHT)
                    for spm in spms_:
                        if spm.prot_uval==prot_uval:
                            cur_spm=spm

        if in_domains:
            dmobj=dmRE.match(l)
            if dmobj:
                cur_spm.add_dpm(dmobj.groups())
    return spms_

def parse_hmmsearchres(nrfpath,acc_mode="simple"):
#first pass looks for domain of interest
    nrfile=open(nrfpath,'r')
    prof_ghname=None
    prof_size=None
    prof_acc=None
    pls_=[]
    motifHT={}
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
            spms_=parse_profmatches(pls_,prof_ghname,prof_acc,prof_size,acc_mode=acc_mode)
            motifHT[prof_acc]=spms_
    return motifHT

def parse_hmmsearchres_toprotdict(nrfpath,acc_mode="simple"):
    motifHT=parse_hmmsearchres(nrfpath,acc_mode=acc_mode)
    protHT={}
    for pfam_motif in motifHT.keys():
        cur_spm=motifHT[pfam_motif]
        curp_uval=cur_spm.prot_uval
        if cur_spm.prot_uval not in protHT.keys():
            protHT[curp_uval]=[]
        protHT[curp_uval].append(cur_spm)
    return protHT

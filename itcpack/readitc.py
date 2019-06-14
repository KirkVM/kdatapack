import re,sys
import numpy as np
from dataclasses import dataclass
from typing import Iterable

@dataclass  
class injectspec:
    injvol: float
    injtime: float
    injspacing: float
    injfilter: float
    lnum:int
@dataclass
class expdetails:
    syringe_lconc:float
    cell_mstartconc: float
    cell_volume: float
@dataclass
class injection:
    injum:float
    injvol:float
    injtime:float
    ltoti:float
    ltotf:float
    mtoti:float
    mtotf:float
    seconds_total:Iterable
    seconds:Iterable
    power:Iterable

#def get_ftots(ltoti,mtoti,expd,injspec):
#   moddy=expd.cell_volume-0.5*injspec.injvol)/(expd.cell_volume+0.5*injspec)
#		Variable moddy=(ts.vo-0.5*ts.inject)/(ts.vo+0.5*ts.inject)
#		ts.mactot[counts]=ts.mactot[counts-1]*moddy
#		ts.ligtot[counts]=(ts.vo*ts.ligtot[counts-1]+ts.inject*ts.syringe)/ts.vo
#		ts.ligtot[counts]*=moddy

def readitc(fpathstr):
    with open(fpathstr,'r') as f:
        lines=f.readlines()
    injectspecs=[]
    injRE=re.compile('\$\s+([\d.]+)\s+,\s+([\d.]+)\s+,\s+([\d.]+)\s+,\s+([\d.]+)\s+')
    for lnum,l in enumerate(lines):
        injobj=injRE.match(l)
        if not injobj: continue
        injectspecs.append(injectspec(*[float(x) for x in injobj.groups()],lnum))
        assert (len(injectspecs)==lnum-injectspecs[0].lnum+1),'missed an injection row!'
    expdRE=re.compile('#\s+([\d.]+)\s+')
    expdvals=[]
    for lnum,l in enumerate(lines):
        if lnum<injectspecs[-1].lnum:continue
        expdobj=expdRE.match(l)
        if expdobj:
            expdvals.append(expdobj.group(1))
    expd=expdetails(*[float(x) for x in [expdvals[1],expdvals[2],expdvals[3]] ])
    inj_initRE=re.compile('@(\d+)(,([\d.]+),([\d.]+)){0,1}')
    injection_primers=[]
    for lnum,l in enumerate(lines):
        iiobj=inj_initRE.match(l)
        if not iiobj:continue
        injnum,inj_startlnum,inj_stoplnum,injvol,injtime=int(iiobj.group(1)),lnum+1,None,None,None
        if injnum!=0:
            injection_primers[-1][2]=lnum
            injvol=float(iiobj.group(3))
            injtime=float(iiobj.group(4))
        injection_primers.append([injnum,inj_startlnum,inj_stoplnum,injvol,injtime])
    injection_primers[-1][2]=lnum+1
    print(expd)
    all_seconds=[]
    all_power=[]
    ltoti=0
    mtoti=expd.cell_mstartconc*1e-3
    syrconc=expd.syringe_lconc*1e-3
    vo=expd.cell_volume*1e-3
    injdataRE=re.compile('([\d.-]+),([\d.-]+),([\d.-]+),([\d.-]+),([\d.-]+),([\d.-]+),([\d.-]+)\s+')
    injections=[]
    for ip in injection_primers:
        cur_seconds=[]
        cur_power=[]
        injlines=lines[ip[1]:ip[2]]
        for injline in injlines:
            injdataobj=injdataRE.match(injline)
            cur_seconds.append(float(injdataobj.group(1)))
            cur_power.append(float(injdataobj.group(2)))
        all_seconds.extend(cur_seconds)
        all_power.extend(cur_power)
        if ip[0]!=0:        
            moddy=(vo-0.5*ip[3]*1e-6)/(vo+0.5*ip[3]*1e-6)
            mtotf=mtoti*moddy
            ltotf=moddy*(vo*ltoti + ip[3]*1e-6*syrconc)/vo
            new_injection=injection(ip[0],ip[3]*1e-6,ip[4],ltoti,ltotf,mtoti,mtotf,\
                          np.array(cur_seconds),np.array(cur_seconds)-cur_seconds[0],np.array(cur_power))
            injections.append(new_injection)
            ltoti=ltotf
            mtoti=mtotf
    for injecty in injections:
        print(injecty)
#        all_seconds.append[ip]
#        if ip[0]!=0
#			Variable moddy=(ts.vo-0.5*ts.inject)/(ts.vo+0.5*ts.inject)
#			ts.mactot[counts]=ts.mactot[counts-1]*moddy
#			ts.ligtot[counts]=(ts.vo*ts.ligtot[counts-1]+ts.inject*ts.syringe)/ts.vo
#			ts.ligtot[counts]*=moddy
			
#	bTerm=-1-activeLigtoti/(activeMactoti)-1/(ts.coefw[0]*(activeMactoti))
#	squareTerm=bTerm^2 - 4*activeLigtoti/(activeMactoti)
#	thetai=0.5*(-bTerm - sqrt(squareTerm))
#	bTerm=-1-activeLigtotf/(activeMactotf)-1/(ts.coefw[0]*(activeMactotf))
#	squareTerm=bTerm^2 - 4*activeLigtotf/(activeMactotf)
#	thetaf=0.5*(-bTerm - sqrt(squareTerm))

#	ts.myHeat = (ts.coefw[1]/(ts.syringe*ts.inject))*ts.vo*(thetaf*(activeMactotf) - thetai*(activeMactoti))
#	ts.myHeat +=0.5*(ts.coefw[1]/(ts.syringe*ts.inject))*ts.inject*(thetaf*(activeMactotf)+thetai*(activeMactoti))
#	ts.myBaseline=ts.coefw[3]*ts.ligtotf/ts.mactotf+ts.coefw[4]

#        #weird taht inject specs does not go past line 50 in header section?...so this can fail
#        if injnum!=0:
#            assert ( float(iiobj.group(3))== injectspecs[injnum-1].injvol)

#        print(iiobj.groups())
#    print(injectspecs)
#    print(expd)
            
        
        
        
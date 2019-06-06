import re,sys
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
    print(injectspecs)
    print(expd)
            
        
        
        
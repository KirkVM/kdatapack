import re,sys
import numpy as np
from dataclasses import dataclass
from typing import Iterable
from pathlib import Path
from itcpack import fititc

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
    injnum:float
    injvol:float
    injtime:float
    seconds_total:Iterable
    seconds:Iterable
    power:Iterable

def readitc(fpathstr):
    """
    returns list of injections as [@dataclass_injection]

    Arguments:
        fpathstr- path to itc file
    
    Returns:
        list of injections-injnum,injvol,injtime,ltoti,ltotf,mtoti,mtotf,seconds_totaal,seconds,power
    """
    with open(fpathstr,'r') as f:
        lines=f.readlines()
    #scan through file-- looking for injection details in header section
    injectspecs=[]
    injRE=re.compile('\$\s+([\d.]+)\s+,\s+([\d.]+)\s+,\s+([\d.]+)\s+,\s+([\d.]+)\s+')
    for lnum,l in enumerate(lines):
        injobj=injRE.match(l)
        if not injobj: continue
        injectspecs.append(injectspec(*[float(x) for x in injobj.groups()],lnum))
        assert (len(injectspecs)==lnum-injectspecs[0].lnum+1),'missed an injection row!'
    #scan through file-- looking for expdetails in header section (currently just syrconc,mconc,cellvol)
    expdRE=re.compile('#\s+([\d.]+)\s+')
    expdvals=[]
    for lnum,l in enumerate(lines):
        if lnum<injectspecs[-1].lnum:continue
        expdobj=expdRE.match(l)
        if expdobj:
            expdvals.append(expdobj.group(1))
    expd=expdetails(*[float(x)*1e-3 for x in [expdvals[1],expdvals[2],expdvals[3]] ])
    #now scan through file and find lines indicating a new injection (start with '@')
    inj_initRE=re.compile('@(\d+)(,([\d.]+),([\d.]+)){0,1}')
    injection_primers=[]
    for lnum,l in enumerate(lines):
        iiobj=inj_initRE.match(l)
        if not iiobj:continue
        injnum,inj_startlnum,inj_stoplnum,injvol,injtime=int(iiobj.group(1)),lnum+1,None,None,None
        if injnum!=0:
            injection_primers[-1][2]=lnum
            injvol=float(iiobj.group(3))*1e-6
            injtime=float(iiobj.group(4))
        injection_primers.append([injnum,inj_startlnum,inj_stoplnum,injvol,injtime])
    injection_primers[-1][2]=lnum+1
    all_seconds=[]
    all_power=[]
    ltoti=0
    mtoti=expd.cell_mstartconc
    syrconc=expd.syringe_lconc
    vo=expd.cell_volume
    injdataRE=re.compile('([\d.-]+),([\d.-]+),([\d.-]+)(,([\d.-]+),([\d.-]+),([\d.-]+),([\d.-]+)\s+)*')
    injections=[]
    #injection primers gives line locations for injection data
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
            moddy=(vo-0.5*ip[3])/(vo+0.5*ip[3])
            mtotf=mtoti*moddy
            ltotf=moddy*(vo*ltoti + ip[3]*syrconc)/vo
            new_injection=injection(ip[0],ip[3],ip[4],\
                          np.array(cur_seconds),np.array(cur_seconds)-cur_seconds[0],np.array(cur_power))
            injections.append(new_injection)
            ltoti=ltotf
            mtoti=mtotf
        else:
            new_injection=injection(0,None,None,\
                          np.array(cur_seconds),np.array(cur_seconds)-cur_seconds[0],np.array(cur_power))
            injections.append(new_injection)

    #itcd=fititc.ITCDataset(expd,injections)
    return expd,injections

#def loaditc(fpathstr,fdir=None,cachedfpathdir='same',cachedfpathname='saved_itc.thermodynamics'):
def loaditc(fpathstr,fdir=None,cachedfpathdir='itc_saved_data',reset_titration=False):#,cachedfpathname='saved_itc.thermodynamics'):
    '''runs readitc and returns an ITCDataset'''
    if fdir is not None:
        fpath=Path(fdir)
        fpath=fpath / fpathstr
    else:
        fpath=Path(fpathstr)
    expdeets,injections=readitc(fpath)
    #now look for an existing
    #filename=fpath.name

    #if cachedfpathdir=='same':
    itcd=fititc.ITCDataset(expdeets,injections,fpath.name)
    cachedfpath=fpath.parent / cachedfpathdir / f'{fpath.stem}.json'
    if cachedfpath.exists() and reset_titration==False:
        itcd.extract_peaks()
        itcd.update_with_stored_vals()
        itcd.make_titrationdf()
    else:
        itcd.get_xspower()
        itcd.extract_peaks()
        itcd.make_titrationdf()
#    itcd.create_titration_data()
    return itcd

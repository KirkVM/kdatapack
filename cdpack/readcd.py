import os,re
import pandas as pd  
import numpy as np
from pathlib import Path

#def readwlscan(fldrpath,fbasename):
def read_wlscan(fpath):
    wlineRE=re.compile('(\d{3}\.\d{3})\s+([0-9\.-]{3,10})\s+([0-9\.]{3,10})\s+([0-9\.]{3,10})')
    wlengths=[]
    cdsignals=[]
    with open(fpath,'r') as f:
        lines=f.readlines()
    for line in lines:
        wline_obj=wlineRE.match(line)
        if wline_obj:
            wlengths.append(int(float(wline_obj.groups()[0])))
            cdsignals.append(float(wline_obj.groups()[1]))
    return wlengths,cdsignals

def read_tempscan(fpath):
    wlineRE=re.compile('([0-9\.]{3,10})\s+([0-9\.-]{3,10})\s+([0-9\.]{3,10})\s+([0-9\.]{3,10})')
    temps=[]
    cdsignals=[]
    with open(fpath,'r') as f:
        lines=f.readlines()
    for line in lines:
        wline_obj=wlineRE.match(line)
        if wline_obj:
            temps.append(float(wline_obj.groups()[0]))
            cdsignals.append(float(wline_obj.groups()[1]))
    return temps,cdsignals


def get_cell_label(fpath,cellnum):
    rotorRE=re.compile('Rotor #(\d) : (\S+)')#{3}\.\d{3})\s+([0-9\.-]{3,10})\s+([0-9\.]{3,10})\s+([0-9\.]{3,10})')
    with open(fpath,'r') as f:
        lines=f.readlines()
    for line in lines:
        rotor_obj=rotorRE.match(line)
        if rotor_obj:
            rnum=int(rotor_obj.groups()[0])
            if rnum==cellnum:
                return rotor_obj.groups()[1]
 

def readset(fldrp,fbasename,xvar='wl'):
    assert(xvar in ['wl','temp']),'xvar should be wl or temp'
    base_fldrpath=Path(fldrp)
    fldrfpaths=base_fldrpath.glob('*.dat')
    cellRE=re.compile('Cell(\d)')
    fpaths=[]
    cdsig_dict={}
    for curfpath in fldrfpaths:
        if curfpath.stem.split(' ')[0]==fbasename:
            curcellnum=int(cellRE.search(curfpath.stem).groups()[0])
            if xvar=='wl':
                xvals,cdsigs=read_wlscan(curfpath)
                clabel=get_cell_label(curfpath,curcellnum)
                cdsig_dict[clabel]=cdsigs
            elif xvar=='temp':
                xvals,cdsigs=read_tempscan(curfpath)
                clabel=get_cell_label(curfpath,curcellnum)
                cdsig_dict[clabel]=cdsigs
                cdsig_dict[clabel+'_temp']=xvals
    #return cdsig_dict
    #cdsig_dict[xvar]=xvals
    longest_array=max([len(cdsig_dict[x]) for x in cdsig_dict])
    for colname in cdsig_dict.keys():
        curlist=cdsig_dict[colname]
        curlist.extend([np.nan for _ in range(longest_array-len(curlist))])
    setdf=pd.DataFrame.from_dict(cdsig_dict)
    return setdf


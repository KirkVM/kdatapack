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
 

def readset(fldrp,fbasename):
    base_fldrpath=Path(fldrp)
    fldrfpaths=base_fldrpath.glob('*.dat')
    cellRE=re.compile('Cell(\d)')
    fpaths=[]
    cdsig_dict={}
    for curfpath in fldrfpaths:
        if curfpath.stem.split(' ')[0]==fbasename:
            curcellnum=int(cellRE.search(curfpath.stem).groups()[0])
            wls,cdsigs=read_wlscan(curfpath)
            clabel=get_cell_label(curfpath,curcellnum)
            cdsig_dict[clabel]=cdsigs
    cdsig_dict['wl']=wls

    setdf=pd.DataFrame.from_dict(cdsig_dict)
    return setdf
#            fpaths.append(curfpath)

import sklearn
import numpy as np
from dataclasses import dataclass
from sklearn import linear_model
from sklearn.preprocessing import PolynomialFeatures
#import scipy,math,sklearn,pickle
#import pandas as pd
#from scipy import optimize
#from typing import Iterable
#from dfitlib import fitmodels

import bqplot
from bqplot import pyplot as plt
from bqplot import Figure
from bqplot import LinearScale
from bqplot import Axis

def make_guided_figure(injpeak,injnum=1):
    fsx=LinearScale()
    fsy=LinearScale()
    xax=Axis(label='seconds',scale=fsx)
    yax=Axis(label='uWatt',scale=fsy,orientation='vertical')
    flpr=bqplot.marks.Lines(x=injpeak.seconds,y=injpeak.xs_power,scales={'x':fsx,'y':fsy},colors=['black'])
    blpts=bqplot.marks.Scatter(x=injpeak.guided_bp_seconds,y=injpeak.guided_bp_power,scales={'x':fsx,'y':fsy})
    rngpts=bqplot.marks.Scatter(x=[injpeak.guided_intstart,injpeak.guided_intstop],\
                                y=[np.quantile(injpeak.final_xs_powerbl,0.25),np.quantile(injpeak.final_xs_powerbl,0.25)],\
                                scales={'x':fsx,'y':fsy},
                                colors=['gold'],marker='rectangle',default_skew=0.99,default_size=400)
    bls=[]
    for blseg in injpeak.blsegs:
        newbl=bqplot.marks.Lines(x=blseg.seconds,y=blseg.blvals,scales={'x':fsx,'y':fsy})
        bls.append(newbl)
    injpeak.plot_bls=bls
    fig=Figure(marks=[flpr,blpts,rngpts,*bls],axes=[xax,yax],min_aspect_ratio=1.2,\
                fig_margin={'top':60, 'bottom':60, 'left':60, 'right':0})
    fig.title=f'Injection {injnum}'
    blpts.enable_move=True
    blpts.restrict_y=True
    blpts.on_drag_end(injpeak.bl_guided_callback)
    rngpts.enable_move=True
    rngpts.restrict_x=True
    rngpts.on_drag_end(injpeak.rng_guided_callback)
    return fig

def make_guided_figure_fixed(injpeak,injnum=1):
    fsx=LinearScale()
    fsy=LinearScale()
    xax=Axis(label='seconds',scale=fsx)
    yax=Axis(label='uWatt',scale=fsy,orientation='vertical')
    fig=Figure(axes=[xax,yax],min_aspect_ratio=1.2,\
                fig_margin={'top':60, 'bottom':60, 'left':60, 'right':0})
#
#    fig=Figure()#axes=[xax,yax],min_aspect_ratio=1.2,\
    return fig


def create_interp_line(x0,x1,y0,y1,seconds):
    poly=PolynomialFeatures(1)
    xtarr=poly.fit_transform(np.expand_dims([x0,x1],1))
    clf = sklearn.linear_model.LinearRegression()
    clf.fit(xtarr,np.expand_dims([y0,y1],1))
    sec_tarr=poly.fit_transform(np.expand_dims(seconds,1))
    return clf.predict(sec_tarr)[:,0]

@dataclass
class Blsegment:
    x0: float
    x1: float
    y0: float
    y1: float
    seconds: list
    blvals: list


class ITCInjectionPeak:
    def __init__(self,seconds,power,smooth_powerbl,xs_power,xsbl,injnum,injvol):
        self.seconds=seconds
        self.power=power
        self.smooth_powerbl=smooth_powerbl
        self.xs_power=xs_power
        self.init_xs_powerbl=xsbl
        self.injnum=injnum
        self.injvol=injvol

        self.init_intstart=self.seconds[0]
        self.init_intstop=self.seconds[-1]
        self.xs_heat=0.0
        self.final_xs_powerbl=self.init_xs_powerbl
        self.final_intstart=self.init_intstart
        self.final_intstop=self.init_intstop

        self.guided_intstart=self.final_intstart
        self.guided_intstop=self.final_intstop
        self.guided_bp_indices=[] #guided_breakpoint_indices...
        self.guided_bp_seconds=[]
        self.guided_bp_power=[]
        self.blsegs=[]
        self.create_guided_bl()
        self.plot_bls=[]
        self.calc_xs_heat()
    
    def save_guided_bl(self):
        self.final_xs_powerbl=list(self.blsegs[0].blvals[:])
        for blseg in self.blsegs[1:]:
            self.final_xs_powerbl.extend(blseg.blvals[1:])
        self.final_xs_powerbl=np.array(self.final_xs_powerbl)
        #self.final_intstart=self.guided_intstart
        #self.final_intstop=self.guided_intstop
        self.calc_xs_heat()

    def create_guided_bl(self):
        steps50sec_idx=[int(x) for x in np.linspace(0,len(self.seconds)-1,int(self.seconds[-1]/50))]#guide pts ~every 50seconds
        if len(steps50sec_idx)<=2: #make sure there are at least 3 guide points
            steps50sec_idx=[int(x) for x in np.linspace(0,len(self.seconds)-1,3)] #the ensuring override
        steps8_idx=[int(x) for x in np.linspace(0,len(self.seconds)-1,8)] #8 evenly-spaced guide points
        self.guided_bp_indices=[max(x,y) for x,y in zip(steps50sec_idx[:10],steps8_idx)]#choose option with bigger step size
        self.guided_bp_seconds=[self.seconds[x] for x in self.guided_bp_indices]
        self.guided_bp_power=[self.final_xs_powerbl[x] for x in self.guided_bp_indices]
        self.blsegs=[]
        for curidx in range(1,len(self.guided_bp_indices)):
            x0=self.seconds[self.guided_bp_indices[curidx-1]]
            x1=self.seconds[self.guided_bp_indices[curidx]]
            y0=self.final_xs_powerbl[self.guided_bp_indices[curidx-1]]
            y1=self.final_xs_powerbl[self.guided_bp_indices[curidx]]
            segseconds=self.seconds[self.guided_bp_indices[curidx-1]:self.guided_bp_indices[curidx]+1]
            segline=create_interp_line(x0,x1,y0,y1,segseconds)
#            import pdb;pdb.set_trace()
            self.blsegs.append(Blsegment(x0,x1,y0,y1,segseconds,segline))            
    def calc_xs_heat(self):
        start_index=list(self.seconds).index(self.final_intstart)
        stop_index=list(self.seconds).index(self.final_intstop)
#        if self.injnum==1:
#            import pdb;pdb.set_trace()
        xs_power_forint=np.sum(self.xs_power[start_index:stop_index+1]-\
                        self.final_xs_powerbl[start_index:stop_index+1])
        self.xs_heat=1e-6*xs_power_forint*((self.final_intstop-self.final_intstart)/(stop_index-start_index))
    
    def bl_guided_callback(self,name,value):
        xval=value['point']['x']
        yval=value['point']['y']
        for blseg in self.blsegs:
            if np.isclose(blseg.x0,xval):
                blseg.y0=value['point']['y']
                blseg.blvals=create_interp_line(blseg.x0,blseg.x1,blseg.y0,blseg.y1,blseg.seconds)
            if np.isclose(blseg.x1,xval):
                blseg.y1=value['point']['y']
                blseg.blvals=create_interp_line(blseg.x0,blseg.x1,blseg.y0,blseg.y1,blseg.seconds)
        for blseg,fig_bl in zip(self.blsegs,self.plot_bls):
            fig_bl.y=blseg.blvals
        newbp_pwr=[]
        for bpidx in range(len(self.guided_bp_seconds)):
            if np.isclose(self.guided_bp_seconds[bpidx],xval):
                newbp_pwr.append(yval)
            else:
                newbp_pwr.append(self.guided_bp_power[bpidx])
        self.guided_bp_power=newbp_pwr
        #self.save_guided_bl()

    def rng_guided_callback(self,name,value):
#        xval=value['point']['x']
        int_start_time=min(name.x)
        start_diffs=[abs(x-int_start_time) for x in self.seconds]
        start_idx=start_diffs.index(min(start_diffs))
        self.final_intstart=self.seconds[start_idx]

        int_stop_time=max(name.x)
        stop_diffs=[abs(x-int_stop_time) for x in self.seconds]
        stop_idx=stop_diffs.index(min(stop_diffs))
        self.final_intstop=self.seconds[stop_idx]
        #self.save_guided_bl()
#        ,self.final_intstop=name.x
#        self.cool=name
#        self.dumb=value
#        self.final_intstart=value['point']
#        self.save_guided_bl()
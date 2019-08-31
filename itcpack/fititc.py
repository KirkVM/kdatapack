import scipy,math,sklearn,pickle
import numpy as np
import pandas as pd
from scipy import optimize
from dataclasses import dataclass
from typing import Iterable
from sklearn import linear_model
from sklearn.preprocessing import PolynomialFeatures
from pathlib import Path
from dfitlib import fitmodels
from itcpack import itcpeaks

import ipywidgets
from ipywidgets import TwoByTwoLayout,GridspecLayout,Button
#def make_makefixy(itcd):
#    gs = GridspecLayout(16, 20)
#    button=Button()
#    button.on_click(on_button_clicked)
#    gs[0,0]=Button()
#    gs[1:,:]=myfig2
#gs

def linear_lossfunc(tup,xs,ys):
    inty,slopey=tup
    estimates= slopey*xs+inty
    delta = estimates - ys
    return np.dot(delta, delta) 

def itc_heats_prediction(Ka,DelH,ltoti,ltotf,mtoti,mtotf,syrconc,injvol,vo):
    bTerm=-1-ltoti/mtoti- 1/(Ka*mtoti)
    squareTerm=np.power(bTerm,2) - 4*ltoti/mtoti
    thetai=0.5*(-bTerm - math.sqrt(squareTerm))
    try:
        bTerm=-1-ltotf/mtotf- 1/(Ka*mtotf)
        squareTerm=np.power(bTerm,2) - 4*ltotf/mtotf
        thetaf=0.5*(-bTerm - math.sqrt(squareTerm))
    except:
        import pdb;pdb.set_trace()

    heat=(DelH/(syrconc*injvol))*vo*(thetaf*mtotf -thetai*mtoti) \
         + 0.5*(DelH/(syrconc*injvol))*injvol*(thetaf*mtotf+thetai*mtoti)
    return heat

def itc_l2lossfunc(kdatup,ltotis,ltotfs,mtotis,mtotfs,injvols,syrconc,vo,ys):#ltotis,ltotfs,mtotis,mtotfs):
    Ka,DelH,Mact=kdatup
    ltoti_active=ltotis
    ltotf_active=ltotfs
    mtoti_active=mtotis*Mact
    mtotf_active=mtotfs*Mact
    syrconc_active=syrconc 
    predicted_heats=[]
    for injidx in range(len(ltotis)):
        predicted_heat=itc_heats_prediction(Ka,DelH,ltoti_active[injidx],ltotf_active[injidx],\
                        mtoti_active[injidx],mtotf_active[injidx],syrconc_active,injvols[injidx],vo)
        predicted_heats.append(predicted_heat)
    delta = [predheat - y for predheat,y in zip(predicted_heats,ys)]

    return np.dot(delta,delta)

#from sklearn.linear_model import Ridge

#def polybl_fitandsubtract(fulltrace,xblpoints,yblpoints):
#button = widgets.Button(description="Click Me!")
output = ipywidgets.Output()

#display(button, output)

def on_button_clicked(b):
    b.description='hi'
    with output:
        print("Button clicked.")

#button.on_click(on_button_clicked)
#def on_button_clicked(b):

class ITCDataset:
    def __init__(self,expdetails,injections,fname):
        self.expdetails=expdetails
        self.syrconc=expdetails.syringe_lconc
        self.vo=expdetails.cell_volume
        self.mtot0=expdetails.cell_mstartconc
        self.injection_details=injections #this includes pre-inject baseline
        self.fname=fname
        self.injection_peaks=[] #this will not include pre-inject baseline
        self.tracedf=None
        self.titrationdf=None
        self.mact=1.0
        self.lact=1.0

        self.fit_Ka=None
        self.fit_DelH=None
        self.fit_stoich=None
        self.build_tracedf()
        #fields for figure...
        self.gs=None
        self.fbutton=None
        self.peak_figs=None
        self.peak_figs_idx=-1
        self.get_xspower()

    def build_interactive(self,start_injnum=1):
        self.peak_figs=[itcpeaks.make_guided_figure(self.injection_peaks[x],injnum=x+1) for x in range(len(self.injection_peaks))]
        self.gs = GridspecLayout(16, 20)
        self.rbutton=Button(description='<')
        self.rbutton.on_click(self.on_rbutton_clicked)
        self.gs[0,0]=self.rbutton

        self.fbutton=Button(description='>')
        self.fbutton.on_click(self.on_fbutton_clicked)
        self.gs[0,1]=self.fbutton
        
        self.svbutton=Button(description='save changes')
        self.svbutton.on_click(self.on_svbutton_clicked)
        self.gs[0,8]=self.svbutton
#        self.gs[1,0]=output
        self.gs[1:,:]=self.peak_figs[start_injnum-1]
        self.peak_figs_idx=start_injnum-1
        return self.gs
    
    def on_rbutton_clicked(self,b):
        if self.peak_figs_idx>0:
            self.peak_figs_idx-=1
            self.gs[1:,:]=self.peak_figs[self.peak_figs_idx]

    def on_fbutton_clicked(self,b):
        if self.peak_figs_idx<len(self.injection_peaks):
            self.peak_figs_idx+=1
            self.gs[1:,:]=self.peak_figs[self.peak_figs_idx]

    def on_svbutton_clicked(self,b):
        for injpk in self.injection_peaks:
            injpk.save_guided_bl()
        self.titrationdf.xs_heat.loc[:]=[x.xs_heat for x in self.injection_peaks]
        #self.titrationdf=self.titrationdf.assign(ndh=self.titrationdf.xs_heat/(self.syrconc*self.titrationdf.injvol))
        self.titrationdf.ndh.loc[:]=self.titrationdf.xs_heat/(self.syrconc*self.titrationdf.injvol)
        #self.update_heats([x.heat for x in self.injection_peaks])

    def build_tracedf(self):
        tracedfs=[]
        for injd in self.injection_details:
            numpnts=len(injd.seconds)
            injdict={'injnum':[injd.injnum for _ in range(numpnts)],
                     'injvol':[injd.injvol for _ in range(numpnts)],
                     'seconds':injd.seconds,'seconds_total':injd.seconds_total,
                     'power':injd.power,'xs_power':injd.power,
                     'spl_power':[0 for _ in range(numpnts)],
                     'spl_power_weight':[0 for _ in range(numpnts)],
                     'smooth_powerbl':[0 for _ in range(numpnts)] }
            tracedfs.append(pd.DataFrame.from_dict(injdict))
        self.tracedf=pd.concat(tracedfs,ignore_index=True)

    def get_xspower(self):
        prev_last5avg=self.tracedf.power.iloc[0]
        for inj in self.tracedf.groupby('injnum'):
            start_iloc=inj[1].index[0]
            stop_iloc=inj[1].index[-1]
            numpnts=inj[1].shape[0]
            if np.isnan(inj[1].injvol.iloc[0]):
                self.tracedf.spl_power.iloc[start_iloc:stop_iloc+1]=self.tracedf.power.iloc[start_iloc:stop_iloc+1]
                self.tracedf.spl_power_weight.iloc[start_iloc:stop_iloc+1]=\
                            [1/((numpnts)**1) for _ in range(numpnts)]
                prev_last5avg=np.mean(inj[1].power.iloc[-5:])
            else:
                last5avg=inj[1].power.iloc[-5:].mean()
                self.tracedf.spl_power.iloc[start_iloc:stop_iloc-4]= \
                                            np.array([last5avg*       \
                                            (x-self.tracedf.seconds_total.iloc[start_iloc])/ \
                                            (self.tracedf.seconds_total.iloc[stop_iloc]-self.tracedf.seconds_total.iloc[start_iloc]) \
                                            + prev_last5avg*    \
                                            (1- (x-self.tracedf.seconds_total.iloc[start_iloc])/ \
                                            (self.tracedf.seconds_total.iloc[stop_iloc]-self.tracedf.seconds_total.iloc[start_iloc])) \
                                            for x in self.tracedf.seconds_total.iloc[start_iloc:stop_iloc-4].values ])
                self.tracedf.spl_power_weight.iloc[start_iloc:stop_iloc-4]=[1/((numpnts-5)**3) for _ in range(numpnts-5)]
                self.tracedf.spl_power.iloc[stop_iloc-4:stop_iloc+1]= inj[1].power.iloc[-5:]
                self.tracedf.spl_power_weight.iloc[stop_iloc-4:stop_iloc+1]=[1 for _ in range(5)]
                prev_last5avg=last5avg
       
        poly=PolynomialFeatures(5)
        xtarr=poly.fit_transform(np.expand_dims(self.tracedf.seconds_total.values,1))
        clf = sklearn.linear_model.Ridge(alpha=0,fit_intercept=False)
        yblcenter=self.tracedf.power.mean()
        yfit=self.tracedf.spl_power-yblcenter
        clf.fit(xtarr,np.expand_dims(yfit,1),sample_weight=self.tracedf.spl_power_weight)
        prediction=clf.predict(xtarr)+yblcenter
        self.tracedf=self.tracedf.assign(smooth_powerbl=prediction)
        self.tracedf.xs_power=self.tracedf.power-self.tracedf.smooth_powerbl

    def create_titration_dataset(self):
#        trcgrps=self.tracedf.groupby('injnum')
        mtotis=[];mtotfs=[];ltotis=[];ltotfs=[];injvols=[];xs_heats=[]
        injdidx=0
        mtoti=self.mtot0
        ltoti=0.0
        prev_last5avg=self.tracedf.xs_power.iloc[0]
        for tracegrp in self.tracedf.groupby('injnum'):
            last5avg=tracegrp[1].xs_power.iloc[-5:].mean()
            cur_injdetail=self.injection_details[injdidx]
            assert cur_injdetail.injnum==tracegrp[1].injnum.iloc[0]
            if np.isnan(tracegrp[1].injvol.iloc[0]) == False:
                #get linear baseline approx to start with...
                xpoints=[tracegrp[1].seconds.iloc[0],tracegrp[1].seconds.iloc[-1]]
                ypoints=[prev_last5avg,last5avg]
                poly=PolynomialFeatures(1)
                xtarr=poly.fit_transform(np.expand_dims(xpoints,1))
                clf = sklearn.linear_model.Ridge(alpha=0,fit_intercept=True)
                clf.fit(xtarr,np.expand_dims(ypoints,1))
                start_xsbl=clf.predict(poly.fit_transform(np.expand_dims(tracegrp[1].seconds.values,1)))
                #done getting approx
                self.injection_peaks.append(itcpeaks.ITCInjectionPeak(tracegrp[1].seconds.values,tracegrp[1].power.values,\
                                                    tracegrp[1].smooth_powerbl.values,tracegrp[1].xs_power.values,start_xsbl[:,0]))
                moddy=(self.vo-0.5*cur_injdetail.injvol)/(self.vo+0.5*cur_injdetail.injvol)
                mtotf=mtoti*moddy
                ltotf=moddy*(self.vo*ltoti + cur_injdetail.injvol*self.syrconc)/self.vo
                mtotis.append(mtoti);mtotfs.append(mtotf)
                ltotis.append(ltoti);ltotfs.append(ltotf)
                xs_heats.append(self.injection_peaks[-1].xs_heat)
                injvols.append(cur_injdetail.injvol)
                mtoti=mtotf
                ltoti=ltotf
            injdidx+=1
            prev_last5avg=last5avg
        self.titrationdf=pd.DataFrame.from_dict({'mtoti':mtotis,'mtotf':mtotfs,'ltoti':ltotis,'ltotf':ltotfs,
                                                'xs_heat':xs_heats,'injvol':injvols})
        self.titrationdf=self.titrationdf.assign(lmratio=self.titrationdf.ltotf/self.titrationdf.mtotf,
                                                 ndh=self.titrationdf.xs_heat/(self.syrconc*self.titrationdf.injvol))
#    def update_heats(self,xs_heats):
#        self.titrationdf.xs_heats=xs_heats
#        self.titrationdf=self.titrationdf.assign(ndh=self.titrationdf.xs_heat/(self.syrconc*self.titrationdf.injvol))
##       self.lmratio=self.ltotfs/self.mtotfs
#        self.ndh_heats=self.xs_heats/(self.syrconc*self.injvols)#np.array([calc_xsheat(x) for x in self.injections[1:]])
    def convenience_fit(self,Ka=1e5,DelH=-4000,Mact=1.0):
        fitvals=scipy.optimize.minimize(itc_l2lossfunc,(Ka,DelH,Mact),args=\
                (self.ltotis[1:],self.ltotfs[1:],self.mtotis[1:],self.mtotfs[1:],self.injvols[1:],self.syrconc,self.vo,\
                 self.ndh_heats[1:]))
        print(fitvals)
        self.fitKa,self.fitDelH,self.fitMact=fitvals.x
    def create_fit_plot(self,numpnts=1000):
        fit_ndhs=[]
        for injidx in range(len(self.ltotis)):
            fit_ndh=itc_heats_prediction(self.fitKa,self.fitDelH,self.ltotis[injidx],self.ltotfs[injidx],\
                        self.mtotis[injidx]*self.fitMact,self.mtotfs[injidx]*self.fitMact,self.syrconc,self.injvols[injidx],self.vo)
            fit_ndhs.append(fit_ndh)
        self.fit_ndhs=fit_ndhs

#    def store_file(self):
#        '''stores key data about titration'''
#        pkl_tracedf=self.tracedf.pickle
#        pkl_df=self.tracedf.pickle
#        smooth_powerbl=self.smooth_powerbl
#    
#    def pickle_itcd(self,fpathstr='itc_saved_files'):
#        pklpath=Path(fpathstr) 
#        pklpath = pklpath / f'pkl_{self.fname}'
#        with open(pklpath,'w') as f:
#            pickle.dump(self,f)
#
    
    def store_itcdata(self,fdirpathstr='itc_saved_data'):
        datapath=Path(fdirpathstr) 
        datapath = datapath / f'data_{self.fname}'
        #dump key data to json...
#        with open(datapath,'w') as f:
#            pickle.dump(self,f)
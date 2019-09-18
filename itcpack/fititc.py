import scipy,sklearn,json,math
import numpy as np
import pandas as pd
from scipy import optimize
from dataclasses import dataclass
from typing import Iterable
from sklearn import linear_model
from sklearn.preprocessing import PolynomialFeatures
from pathlib import Path
from lmfit import Parameters,minimize
from dfitlib import fitmodels
from itcpack import itcpeaks

import ipywidgets
from ipywidgets import TwoByTwoLayout,GridspecLayout,Button

def linear_lossfunc(tup,xs,ys):
    inty,slopey=tup
    estimates= slopey*xs+inty
    delta = estimates - ys
    return np.dot(delta, delta) 

def itc_heats_prediction(Ka,DelH,ltoti,ltotf,mtoti,mtotf,injvol,syrconc,vo):
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
                        mtoti_active[injidx],mtotf_active[injidx],injvols[injidx],syrconc_active,vo)
        predicted_heats.append(predicted_heat)
    delta = [predheat - y for predheat,y in zip(predicted_heats,ys)]

    return np.dot(delta,delta)

def itc_titration_residual(params,ltotis,ltotfs,mtotis,mtotfs,injvols,syrconc,vo,ys):
    logKa=params['logKa']
    DelH=params['DelH']
    Mact=params['Mact']
    DilHeat=params['DilHeat']
    ltoti_active=ltotis
    ltotf_active=ltotfs
    mtoti_active=mtotis*Mact
    mtotf_active=mtotfs*Mact
    syrconc=syrconc 
    predicted_heats=[]

    for injidx in range(len(ltotis)):
        predicted_heat=itc_heats_prediction(np.power(10,logKa),DelH,ltoti_active[injidx],ltotf_active[injidx],\
                        mtoti_active[injidx],mtotf_active[injidx],injvols[injidx],syrconc,vo)
        predicted_heats.append(predicted_heat+DilHeat)
#    print(ys-predicted_heats)
    return ys-predicted_heats
    #delta = [predheat - y for predheat,y in zip(predicted_heats,ys)]
    #
    #return np.dot(delta,delta)



class ITCDataset:
    def __init__(self,expdetails,injections,fname,fdirpathstr='itc_saved_data'):#,reset_titration=True):
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
        self.fit_DilHeat=None
        self.fit_ndhs=None
        self.fit_params=None
        self.build_tracedf()
        #fields for figure...
        self.gs=None
        self.fbutton=None
        self.peak_figs=None
        self.peak_figs_idx=-1
#        if reset_titration:
#            self.get_xspower()

        self.stored_file_path=None
        jfname=Path(self.fname).stem
        self.stored_file_path=Path(fdirpathstr) / f'{jfname}.json'

    def adjust_peaks(self,start_injnum=1):
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
        self.titrationdf.ndh.loc[:]=self.titrationdf.xs_heat/(self.syrconc*self.lact*self.titrationdf.injvol)
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

    def make_titrationdf(self):
        mtotis=[];mtotfs=[];ltotis=[];ltotfs=[];injvols=[];xs_heats=[]
        injdidx=0
        injnum=1
        mtoti=self.mtot0*self.mact
        syr_active=self.syrconc*self.lact
        ltoti=0.0
        
        for ival,injpeak in enumerate(self.injection_peaks):
            assert(ival+1==injpeak.injnum)
            moddy=(self.vo-0.5*injpeak.injvol)/(self.vo+0.5*injpeak.injvol)
            mtotf=mtoti*moddy
            ltotf=moddy*(self.vo*ltoti + injpeak.injvol*syr_active)/self.vo
            mtotis.append(mtoti);mtotfs.append(mtotf)
            ltotis.append(ltoti);ltotfs.append(ltotf)
            injvols.append(injpeak.injvol)
            mtoti=mtotf
            ltoti=ltotf

        xs_heats=[x.xs_heat for x in self.injection_peaks]#.append(self.injection_peaks[-1].xs_heat)
        self.titrationdf=pd.DataFrame.from_dict({'mtoti':mtotis,'mtotf':mtotfs,'ltoti':ltotis,'ltotf':ltotfs,
                                                'xs_heat':xs_heats,'injvol':injvols})
        self.titrationdf=self.titrationdf.assign(lmratio=self.titrationdf.ltotf/self.titrationdf.mtotf,
                                                 ndh=self.titrationdf.xs_heat/(syr_active*self.titrationdf.injvol))
#        self.titrationdf.ndh.loc[:]=self.titrationdf.xs_heat/(syr_active*self.titrationdf.injvol)

    def extract_peaks(self):
        injidx=0
        prev_last5avg=self.tracedf.xs_power.iloc[0]
        for tracegrp in self.tracedf.groupby('injnum'):
            last5avg=tracegrp[1].xs_power.iloc[-5:].mean()
            cur_injdetail=self.injection_details[injidx]
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
                curinjpeak=itcpeaks.ITCInjectionPeak(tracegrp[1].seconds.values,tracegrp[1].power.values,\
                                        tracegrp[1].smooth_powerbl.values,tracegrp[1].xs_power.values,\
                                        start_xsbl[:,0],cur_injdetail.injnum,cur_injdetail.injvol)
                self.injection_peaks.append(curinjpeak)
            prev_last5avg=last5avg
            injidx+=1

    def superfit(self,logKa=1e5,DelH=-4000,Mact=1.0,DilHeat=-500,pts='all'):
        params = Parameters()
        params.add('logKa', value=logKa,min=3,max=7)
        params.add('DelH', value=DelH)
        params.add('Mact', value=Mact,min=0.5,max=2.5)
        params.add('DilHeat', value=DilHeat,min=-2000,max=2000)
#        params.add('DilHeat', value=0,vary=False)#,min=-1000,max=1000)
#        if pts=='all':
#            itc_args=[self.ltotis[1:],self.ltotfs[1:],self.mtotis[1:],self.mtotfs[1:],\
#                [x.injvol for x in self.injection_details[1:]],self.syrconc*self.lact,self.vo,\
#                self.titrationdf.xs_heat.values[1:]]
#        out = minimize(itc_titration_residual, params, args=(*itc_args))
        fit_out = minimize(itc_titration_residual, params, \
               args=(self.titrationdf.ltoti.values[1:],self.titrationdf.ltotf.values[1:],\
                self.titrationdf.mtoti.values[1:],self.titrationdf.mtotf.values[1:],\
                self.titrationdf.injvol.values[1:],self.syrconc*self.lact,self.vo,\
                self.titrationdf.ndh.values[1:]))
        self.fit_params=fit_out.params
        self.fit_logKa=fit_out.params['logKa'].value
        self.fit_DelH=fit_out.params['DelH'].value
        self.fit_stoich=fit_out.params['Mact'].value
        self.fit_DilHeat=fit_out.params['DilHeat'].value
#        print(params['Mact'])
        return fit_out

    def create_fit_plot(self,numpnts=1000):
        fit_ndhs=[]
        for injidx in range(self.titrationdf.shape[0]):
            fit_ndh=itc_heats_prediction(np.power(10,self.fit_logKa),self.fit_DelH,self.titrationdf.ltoti.values[injidx],\
                    self.titrationdf.ltotf.values[injidx],self.titrationdf.mtoti.values[injidx]*self.fit_stoich,\
                    self.titrationdf.mtotf.values[injidx]*self.fit_stoich,self.titrationdf.injvol.values[injidx],\
                    self.syrconc,self.vo)
            fit_ndhs.append(fit_ndh+self.fit_DilHeat)
        self.fit_ndhs=fit_ndhs
#        self.fit_lmratios=self.titrationdf.ltotf.values/(self.titrationdf.mtotf.values*self.fit_stoich)

    
    def store_itcdata(self,fdirpathstr='itc_saved_data'):
        store_dict={}

        store_dict['syrconc']=self.syrconc
        store_dict['vo']=self.vo
        store_dict['mtot0']=self.mtot0
        store_dict['mact']=self.mact
        store_dict['lact']=self.lact
        #store_dict['fit_Ka']=self.fit_Ka
        #store_dict['fit_DelH']=self.fit_DelH
        #store_dict['fit_stoich']=self.fit_stoich
        store_dict['xs_power']=self.tracedf.xs_power.values.tolist()
        store_dict['smooth_powerbl']=self.tracedf.xs_power.values.tolist()

        store_dict['injection_peaks']={}
        for injpeak in self.injection_peaks:
            store_dict['injection_peaks'][injpeak.injnum]={}
            store_dict['injection_peaks'][injpeak.injnum]['final_intstart']=injpeak.final_intstart
            store_dict['injection_peaks'][injpeak.injnum]['final_intstop']=injpeak.final_intstop
            store_dict['injection_peaks'][injpeak.injnum]['final_xs_powerbl']=injpeak.final_xs_powerbl.tolist()
            store_dict['injection_peaks'][injpeak.injnum]['xs_power']=injpeak.xs_power.tolist()
            store_dict['injection_peaks'][injpeak.injnum]['seconds']=injpeak.seconds.tolist()
            store_dict['injection_peaks'][injpeak.injnum]['guided_bp_seconds']=injpeak.guided_bp_seconds
            store_dict['injection_peaks'][injpeak.injnum]['guided_bp_power']=injpeak.guided_bp_power
            store_dict['injection_peaks'][injpeak.injnum]['guided_bp_indices']=injpeak.guided_bp_indices
            store_dict['injection_peaks'][injpeak.injnum]['blsegs']=[]
            for blseg in injpeak.blsegs:
                curseg={}
                curseg['x0']=blseg.x0
                curseg['x1']=blseg.x1
                curseg['y0']=blseg.y0
                curseg['y1']=blseg.y1
                curseg['seconds']=blseg.seconds.tolist()
                curseg['blvals']=blseg.blvals.tolist()
                store_dict['injection_peaks'][injpeak.injnum]['blsegs'].append(curseg)
                #import pdb;pdb.set_trace()
        with open(self.stored_file_path,'w') as f:
            json.dump(store_dict,f)

    def update_with_stored_vals(self):
        store_dict={}
        with open(self.stored_file_path,'r') as f:
            store_dict=json.load(f)
        self.tracedf.smooth_powerbl=np.array(store_dict['smooth_powerbl'])
        self.tracedf.xs_power=np.array(store_dict['xs_power'])
        self.syrconc=store_dict['syrconc']
        self.vo=store_dict['vo']
        self.mtot0=store_dict['mtot0']
        self.mact=store_dict['mact']
        self.lact=store_dict['lact']
        #self.fit_Ka=store_dict['fit_Ka']
        #self.fit_DeH=store_dict['fit_DelH']
        #self.fit_stoich=store_dict['fit_stoich']
        for injpeak in self.injection_peaks:
            stored_peak_dict=store_dict['injection_peaks'][str(injpeak.injnum)]
            injpeak.final_intstart=stored_peak_dict['final_intstart']
            injpeak.final_intstop=stored_peak_dict['final_intstop']
            injpeak.guided_intstart=stored_peak_dict['final_intstart']
            injpeak.guided_intstop=stored_peak_dict['final_intstop']
            injpeak.final_xs_powerbl=np.array(stored_peak_dict['final_xs_powerbl'])
            injpeak.xs_power=np.array(stored_peak_dict['xs_power'])
            injpeak.seconds=np.array(stored_peak_dict['seconds'])
#            injpeak.create_guided_bl()
            cursegs=[]
            for sblseg in stored_peak_dict['blsegs']:
                bl_x0=sblseg['x0']
                bl_x1=sblseg['x1']
                bl_y0=sblseg['y0']
                bl_y1=sblseg['y1']
                bl_seconds=np.array(sblseg['seconds'])#.tolist()#=blseg.seconds.tolist()
                bl_blvals=np.array(sblseg['blvals'])
                curblseg=itcpeaks.Blsegment(bl_x0,bl_x1,bl_y0,bl_y1,bl_seconds,bl_blvals)
                cursegs.append(curblseg)
                injpeak.guided_bp_seconds=stored_peak_dict['guided_bp_seconds']
                injpeak.guided_bp_power=stored_peak_dict['guided_bp_power']
                injpeak.guided_bp_indices=stored_peak_dict['guided_bp_indices']
            injpeak.blsegs=cursegs
            injpeak.save_guided_bl()
            injpeak.calc_xs_heat()
        self.make_titrationdf()

import scipy,math
from scipy import optimize
import numpy as np
from dataclasses import dataclass
from typing import Iterable
from dfitlib import fitmodels

def linear_lossfunc(tup,xs,ys):
    inty,slopey=tup
    estimates= slopey*xs+inty
    delta = estimates - ys
    return np.dot(delta, delta) 

def calc_xsheat(inj):
    '''subtract baseline as linear fit to injection.power[0]-avg(injection.power[-6:-1])'''
    fpoints=np.array([inj.power[0],np.average(inj.power[-6:-1])])
    blinesecs=np.array([inj.seconds[0],inj.seconds[-1]])
    fstats=scipy.optimize.minimize(linear_lossfunc,(fpoints[0],0),args=(blinesecs,fpoints))
    powerbline=fstats.x[0]+fstats.x[1]*np.array(inj.seconds)
    xspower=np.array(inj.power)-powerbline
    xsheat=1e-6*np.sum(xspower*(inj.seconds[-1]/len(inj.seconds)))
    return xsheat

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
    

class ITCDataset:
    def __init__(self,expdetails,injections):
        self.expdetails=expdetails
        self.injections=injections
        self.syrconc=expdetails.syringe_lconc
        self.vo=expdetails.cell_volume
        self.xsheats=None
        self.ltotis=None
        self.ltotfs=None
        self.mtotis=None
        self.mtotfs=None
        self.injvols=None
        self.lmratio=None
        self.fitKa=None
        self.fitDelH=None
        self.fitMact=None
        self.fit_ndhs=None
    def create_titration_dataset(self):
        self.xs_heats=np.array([calc_xsheat(x) for x in self.injections[1:]])
        self.ltotis=np.array([x.ltoti for x in self.injections[1:]])
        self.ltotfs=np.array([x.ltotf for x in self.injections[1:]])
        self.mtotis=np.array([x.mtoti for x in self.injections[1:]])
        self.mtotfs=np.array([x.mtotf for x in self.injections[1:]])
        self.injvols=np.array([x.injvol for x in self.injections[1:]])
        self.lmratio=self.ltotfs/self.mtotfs
        self.ndh_heats=self.xs_heats/(self.syrconc*self.injvols)#np.array([calc_xsheat(x) for x in self.injections[1:]])
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




#	bTerm=-1-activeLigtoti/(activeMactoti)-1/(ts.coefw[0]*(activeMactoti))
#	squareTerm=bTerm^2 - 4*activeLigtoti/(activeMactoti)
#	thetai=0.5*(-bTerm - sqrt(squareTerm))
#	bTerm=-1-activeLigtotf/(activeMactotf)-1/(ts.coefw[0]*(activeMactotf))
#	squareTerm=bTerm^2 - 4*activeLigtotf/(activeMactotf)
#	thetaf=0.5*(-bTerm - sqrt(squareTerm))
#	ts.myHeat = (ts.coefw[1]/(ts.syringe*ts.inject))*ts.vo*(thetaf*(activeMactotf) - thetai*(activeMactoti))
#	ts.myHeat +=0.5*(ts.coefw[1]/(ts.syringe*ts.inject))*ts.inject*(thetaf*(activeMactotf)+thetai*(activeMactoti))
#	ts.myBaseline=ts.coefw[3]*ts.ligtotf/ts.mactotf+ts.coefw[4]


#	ts.myHeat = (ts.coefw[1]/(ts.syringe*ts.inject))*ts.vo*(thetaf*(activeMactotf) - thetai*(activeMactoti))
#	ts.myHeat +=0.5*(ts.coefw[1]/(ts.syringe*ts.inject))*ts.inject*(thetaf*(activeMactotf)+thetai*(activeMactoti))
#	ts.myBaseline=ts.coefw[3]*ts.ligtotf/ts.mactotf+ts.coefw[4]

#        #weird taht inject specs does not go past line 50 in header section?...so this can fail
#        if injnum!=0:
#            assert ( float(iiobj.group(3))== injectspecs[injnum-1].injvol)

#        print(iiobj.groups())
#    print(injectspecs)
#    print(expd)
            
        
  
    
    #def get_ftots(ltoti,mtoti,expd,injspec):
#   moddy=expd.cell_volume-0.5*injspec.injvol)/(expd.cell_volume+0.5*injspec)
#		Variable moddy=(ts.vo-0.5*ts.inject)/(ts.vo+0.5*ts.inject)
#		ts.mactot[counts]=ts.mactot[counts-1]*moddy
#		ts.ligtot[counts]=(ts.vo*ts.ligtot[counts-1]+ts.inject*ts.syringe)/ts.vo
#		ts.ligtot[counts]*=moddy



#def do11fit(injections):

#        all_seconds.append[ip]


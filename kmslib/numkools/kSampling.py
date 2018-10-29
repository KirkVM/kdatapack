import numpy as np
import random

#def get_resample_old(thelist,wtdA_=None):
#    '''takes an array or list and optionally a weighting array (wtdA_=...)
#    and returns a resampled with replacement version with the selections
#    weighted accordingly'''
#    if wtdA_ is None:wtdA_=np.ones((len(thelist)))
#    wtdSumA_=np.zeros((len(thelist)))
#    for x in range(len(wtdA_)):
#        wtdSumA_[x:]+=wtdA_[x]
#    newsample_=[]
#    for x in range(len(thelist)):
#        randval=random.random()*wtdSumA_[len(wtdSumA_)-1]       
#        for y in range(len(wtdSumA_)):
#            if randval<wtdSumA_[y]:
#                newsample_.append(thelist[y])
#                break
##    return np.array(newsample_)
#    return newsample_

import pdb
def get_resample(vals_,wtdA_=None,valsrvs_=None):
    """Returns resampled array/list. 
    
    Args:
        vals: list of vals/entries. Agnostic on contents
    Optional Args and Defaults
        wtdA_=None. weights selection using linear scaling according to vals
        valsrvs_=None. the rv_continuous probability distribution for each corresponding value
            in vals. simply calls the rv .cdf method to sample from that distribution
    """
    thelist=vals_[:] #start by making a copy of the list. 
    if valsrvs_ is not None:
        thelist=[]
        for x0,litems in enumerate(vals_):
            rvs=valsrvs_[x0]
            if type(litems)==list or type(rvs)==list:
                assert(len(litems)==len(rvs))
            elif type(litems)==tuple or type(rvs)==tuple:
                assert(len(litems)==len(rvs))
            else:
                #litemslist=[litems]
                #rvslist=[rvs]
                litems=[litems]
                rvs=[rvs]
            newentrys=[] 
            for x1,litem in enumerate(litems):
                #rv=rvslist[x1]
                rv=rvs[x1]

                #selection should be 0.01-0.99 (obv should use a different random module)
                randval=max(0.01,float(int(random.random()*100))/100)
                sampledval=rv.ppf(randval)
#                print(randval,sampledval,rv.kwds)
                newentrys.append(sampledval)
            #if each entry is a single item, revert to a single value-
            if len(newentrys)==1: newentrys=newentrys[0]
            thelist.append(newentrys)
#        pdb.set_trace()
    if wtdA_ is None:wtdA_=np.ones((len(thelist)))
    wtdSumA_=np.zeros((len(thelist)))
    for x in range(len(wtdA_)):
        wtdSumA_[x:]+=wtdA_[x]
    newsample_=[]
    for x in range(len(thelist)):
        randval=random.random()*wtdSumA_[len(wtdSumA_)-1]       
        for y in range(len(wtdSumA_)):
            if randval<wtdSumA_[y]:
                newsample_.append(thelist[y])
                break
#    return np.array(newsample_)
    return newsample_

def get_resample_weird(vals_,wtdA_=None,valsrvs_=None):
    """Returns resampled array/list. 
    
    Args:
        vals: list of vals/entries. Agnostic on contents
    Optional Args and Defaults
        wtdA_=None. weights selection using linear scaling according to vals
        valsrvs_=None. the rv_continuous probability distribution for each corresponding value
            in vals. simply calls the rv .cdf method to sample from that distribution
    """
    thelist=vals_[:] #start by making a copy of the list. 
    if valsrvs_ is not None:
        thelist=[]
        for x0,litems in enumerate(vals_):
            rvs=valsrvs_[x0]
            if type(litems)==list or type(rvs)==list:
                assert(len(litems)==len(rvs))
            elif type(litems)==tuple or type(rvs)==tuple:
                assert(len(litems)==len(rvs))
            else:
                #litemslist=[litems]
                #rvslist=[rvs]
                litems=[litems]
                rvs=[rvs]
            newentrys=[] 
            for x1,litem in enumerate(litems):
                #rv=rvslist[x1]
                rv=rvs[x1]
                if rv is not None:
                    #selection should be 0.01-0.99 (obv should use a different random module)
                    randval=max(0.01,float(int(random.random()*100))/100)
                    sampledval=rv.ppf(randval)
#                   print(randval,sampledval,rv.kwds)
                    newentrys.append(sampledval)
                else:
                    newentrys.append(litem)
            #if each entry is a single item, revert to a single value-
            if len(newentrys)==1: newentrys=newentrys[0]
            thelist.append(newentrys)
#        pdb.set_trace()
    if wtdA_ is None:wtdA_=np.ones((len(thelist)))
    wtdSumA_=np.zeros((len(thelist)))
    for x in range(len(wtdA_)):
        wtdSumA_[x:]+=wtdA_[x]
    newsample_=[]
    for x in range(len(thelist)):
        randval=random.random()*wtdSumA_[len(wtdSumA_)-1]       
        for y in range(len(wtdSumA_)):
            if randval<wtdSumA_[y]:
                newsample_.append(thelist[y])
                break
#    return np.array(newsample_)
    return newsample_
def wtdresample(npA_,wtdA_=None,iters=100):
    samples_=[]
    for x in range(iters):
        sample=get_resample(npA_,wtdA_=wtdA_)
        samples_.append(sample)
    return samples_

def get_bootstrapstats(npA_,iters=100):
    averagesA_=np.zeros((iters))
    for x in range(iters):
        sampledA_=sklearn.utils.resample(npA_)
        averagesA_[x]=np.average(sampledA_)
    #lb=sorted(averagesA_)[int(iters/3)]
    #ub=sorted(averagesA_)[int(iters*2/3)]
    #print lb,ub
    #sample_avg=np.average(averagesA_)
    #uncertainty=max(sample_avg-lb,ub-sample_avg)
    #return sample_avg,uncertainty
    return np.average(averagesA_),np.std(averagesA_)


class KBSParams:
    def __init__(self):
        self.paramHT={}
        self.xconfA_=None
        self.ylconfA_=None
        self.yuconfA_=None
        self.param_sampleAHT={}
    def add_param(self,pname,lc,uc):
        self.paramHT[pname]=[lc,uc]
    def add_param_sampleAHT(self,pname,sampleA_):
        self.param_sampleAHT[pname]=sampleA_
    def add_ciplot(self,xcA_,ylcA_,yucA_):
        self.xconfA_=xcA_ 
        self.ylconfA_=ylcA_ 
        self.yuconfA_=yucA_ 

def get_linearKBSP(intsA_,slopesA_,lconf=10,uconf=90,getciplot=True,ciplotlen=100,
        cix_start_stop=None):
    kbsObj=KBSParams()
    minslope=np.percentile(slopesA_,10)
    maxslope=np.percentile(slopesA_,90)
    kbsObj.add_param('slope',minslope,maxslope)
    kbsObj.add_param_sampleAHT('slope',slopesA_)
    minint=np.percentile(intsA_,10)
    maxint=np.percentile(intsA_,90)
    kbsObj.add_param('int',minint,maxint)
    kbsObj.add_param_sampleAHT('int',intsA_)
    if getciplot:
        boundxA_=np.linspace(cix_start_stop[0],cix_start_stop[1],ciplotlen)

        #add array to bsys_ for each bootstrap iter
        bsys_=[]
        for a,b in zip(intsA_,slopesA_):
            curYA_=np.zeros(len(boundxA_))
            for x in range(len(boundxA_)):
                curYA_[x]=a+boundxA_[x]*b
            bsys_.append(curYA_)
            
        # for each position in plot,
        # make an array then add percentile to corresponding plot
        bsmaxYA_=np.zeros((len(boundxA_)))
        bsminYA_=np.zeros((len(boundxA_)))
        for x in range(len(boundxA_)):
            curYA_=np.array([f[x] for f in bsys_])
            bsmaxYA_[x]=np.percentile(curYA_,uconf)
            bsminYA_[x]=np.percentile(curYA_,lconf)
        kbsObj.add_ciplot(boundxA_,bsminYA_,bsmaxYA_)

    return kbsObj

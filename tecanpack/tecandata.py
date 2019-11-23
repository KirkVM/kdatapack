import sklearn
import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression,RANSACRegressor

def cleandf(df,overs='nan'):
    cleandf=df.copy()
    #if df has 'OVER' in measurements, it will not be type float
    if df.measurement.dtype != float:
        if overs=='nan':
            cleandf=cleandf.replace(to_replace={'measurement':{'OVER':np.nan}})
    return cleandf

def build_simple_yldconverter(calint,calslope):#make_cylinder_volume_func(r):
    def yield_converter(measurement):
        return (measurement-calint)/calslope 
    return yield_converter

def build_linear_converter(int,slope):
    def linear_converter(value):
        return (value-int)/slope
    return linear_converter
#def build_simple_ebgsubtracter()

def set_relative_dlnfactor(dfrow):
#    if dfrow.experimenter.split(' ')[0]=="Evan"
    if dfrow.detection=='BCA':
        default_predvlp_dlnfactor=5/105
    elif dfrow.detection=='DNS':
        default_predvlp_dlnfactor=1/3
    return dfrow.predvlp_dlnfactor/default_predvlp_dlnfactor

def assign_rxn_types(df):
    df=df.assign(eblank_status=lambda x: (x.ename.isnull()==False) & (x.sname.isnull()) & (x.standardname.isnull()),
                 sblank_status=lambda x: (x.sname.isnull()==False) & (x.ename.isnull()) & (x.standardname.isnull()),
                 calstd_status=lambda x: (x.standardname.isnull()==False) & (x.sname.isnull()) & (x.ename.isnull()),
                 ebg_status=lambda x: (x.sname.isnull()) & (x.standardname.isnull()) & (x.ename=='wge'),
                 rxn_status=lambda x: (x.sname.isnull()==False) & (x.sconc>0) & (x.ename.isnull()==False) & (x.econc>0))
                 #sbg_status=lambda x: (x.sname.isnull()) & (x.standardname.isnull()) & (x.sname=='?'))
#    df=df.assign()
#    df=df.assign(econc_mgmL_adj=lambda x:x.econc_mgmL*x.predvlp_dlnfactor/x.base_detection_factor)
#ebdf_bca['econc_mgmL_adj']=ebdf_bca.apply(lambda x:x.econc_mgmL*x.predvlp_dlnfactor/(5/105),axis=1)

#        get_adjusted_econc(self.expdf,self.base_predvlp_dlnfactor)
    return df

def assign_dlnadjusted_values(df):
    df=df.assign(econc_mgmL_reldvlpadj=lambda x:x.econc_mgmL*x.relative_dlnfactor)
    df=df.assign(econc_molar_reldvlpadj=lambda x:x.econc_molar*x.relative_dlnfactor)
    df=df.assign(econc_reldvlpadj=lambda x:x.econc*x.relative_dlnfactor)
    df=df.assign(sconc_reldvlpadj=lambda x:x.sconc*x.relative_dlnfactor)
    df=df.assign(standardconc_reldvlpadj=lambda x:x.standardconc*x.relative_dlnfactor)
    return df

def add_interpreted_columns(df):
    df=assign_rxn_types(df)
    df.loc[:,'relative_dlnfactor']=df.copy().apply(set_relative_dlnfactor,axis=1)
    df=assign_dlnadjusted_values(df)
    return df

def calc_rxn_yields(df):
    eblankdf=df[df.eblank_status]
    sblankdf=df[df.sblank_status]
    ebgdf=df[df.ebg_status]
    ####if purified enzyme
    #step 1: remove enzyme-only background (by default, use globally-averaged one based on mgmL and known dilution scheme)
    #step 2: remove substrate-only background (by default use averaged from that day. correct for dilution scheme)
    return df
#def get_eblankdf_bca(df):

def calc_bca_ebgconverter(df):
    ebdf_bca=df[(df['eblank_status']) & (df['detection']=='BCA')]
    poly=PolynomialFeatures(1)
    xtarr=poly.fit_transform(np.expand_dims(ebdf_bca.econc_mgmL_reldvlpadj,1))
    rlf=LinearRegression(fit_intercept=False)
    ransac=RANSACRegressor(base_estimator=rlf,min_samples=0.9,residual_threshold=0.5)
    #ransac.fit(xtarr,ebdf_bca.measurement)
    ransac.fit(xtarr,ebdf_bca.re_yield)
    #these lines determine ids of outliers
    maskdf=ebdf_bca.assign(inliers=ransac.inlier_mask_)
    maskdf[['ename','econc','expdate']][maskdf.inliers==False]
    #now return a conversion function built with this model
    print(f'using default ebg conversion of int={ransac.estimator_.coef_[0]} and slope={ransac.estimator_.coef_[1]}')
    return build_linear_converter(*ransac.estimator_.coef_)
        #calculate bca_eblank_default by
        #   (1) grouping into dilution schemes
        #   (2) get correct to get into one curve with 5/105 and 1
        #   (3) fit line to get avg re_yield vs mgmL_conc equation

#def get_sblavg_ydict(df):
#    sbldf=df[df['sblank_status']]
#    sbldict={}
#    sblgrps=sbldf.groupby('sname')
#    for
def calc_re_yield(dfrow,measure_field='measurement',bca_converter=None,dns_converter=None):
    if dfrow.calstd_status==True:
        return np.nan
#        self.alldf.re_yield=self.alldf.measurement.where((self.alldf.detection=='DNS') & (self.alldf.calstd_status==False)).\
    assert(dfrow.detection in ['DNS','BCA'])
    if dfrow.detection=='DNS':
        return(dns_converter(dfrow[measure_field]))
    elif dfrow.detection=='BCA':
        return(bca_converter(dfrow[measure_field]))
#        self.alldf['re_yield']=self.alldf.apply(calc_re_yield,bca_converter=BCA_YldConverter,dns_converter=DNS_YldConverter)

def calc_sbg_yield(dfrow,slopedict):#,sname,sctrl_status):
    if dfrow.sname is not None:
        return slopedict[dfrow.sname]*dfrow.sconc
    else:
        return np.nan

def calc_ebg_yield(dfrow,eslope):#,sname,sctrl_status):
    if dfrow.econc is not None:
        if dfrow.detection=='BCA':
            return dfrow.econc_mgmL*eslope
        else:
            return 0
    else:
        return np.nan

class CAZyExperiment:
    def __init__(self,selected_df,alldf=None,expname=None):
        self.expdf=selected_df.copy()
        self.expdf=assign_rxn_types(self.expdf)
        #self.expdf=assign_rxn_types(self.expdf)
        self.expdf=add_interpreted_columns(self.expdf)#assign_rxn_types(self.expdf)
        if alldf is not None:
            #self.alldf=assign_rxn_types(alldf.copy())
            self.alldf=add_interpreted_columns(alldf.copy())#assign_rxn_types(self.expdf)
        if expname is not None:
            self.expname=expname
        else:
            self.expname=str(self.expdf.expdate.unique()[0])
        self.default_calibration_standard={'DNS':{'slope':0.45,'int':0.05}}
        self.default_calibration_standard.update({'BCA':{'slope':3.36,'int':0.23}})
        self.base_bca_predvlp_dlnfactor=5/105
        #default ctrl slopes are in units of  (mg/mL re_yield equivalents)/(mg/mL e or s)
        self.default_sctrl_slopes={'cmc':0.012,'corn stover':0.003,'galactomannan':0.002,'lichenan':0.010,
                                   'mannan':0.006,'pasc':0.002,'poplar':0.000,'sorghum':0.000,'switchgrass':0.000,
                                   'xylan':0.012,'xyloglucan':0.005}
        self.default_ectrl_slope=0.36 #

    def calc_yields(self,colname='re_yield',subset_selections={},calibration_standard='default',\
                    s_background='default',e_background='default'):
        assert (calibration_standard in ['default']),\
            "calibration_standard must ==['default']"
        if calibration_standard=='default':
            DNS_YldConverter=build_simple_yldconverter(self.default_calibration_standard['DNS']['int'],\
                                                        self.default_calibration_standard['DNS']['slope'])
            BCA_YldConverter=build_simple_yldconverter(self.default_calibration_standard['BCA']['int'],\
                                                        self.default_calibration_standard['BCA']['slope'])
        #now get yields and stuff 
        self.expdf=cleandf(self.expdf)
        self.expdf['re_yield']=self.expdf.apply(calc_re_yield,measure_field='measurement',\
                                            bca_converter=BCA_YldConverter,dns_converter=DNS_YldConverter,axis=1)
        self.alldf=cleandf(self.alldf)
        self.alldf['re_yield']=self.alldf.apply(calc_re_yield,measure_field='measurement',\
                                            bca_converter=BCA_YldConverter,dns_converter=DNS_YldConverter,axis=1)
        
        #get the econc - yield default converter  
        #out of laziness (or thoroughness), currently calculating this each time instead of using default values
        #(like I am for sbldict)     
        self.alldf=cleandf(self.alldf)
        if e_background=='default':
            self.ebl_slope=self.default_ectrl_slope#=0.37 #
        else:
            self.ebg_converter=calc_bca_ebgconverter(self.alldf)
        #now here could get a dictionary of default substrate baselines
        #this would need to change if use diff't than default converters or change default converters
#        self.sbldict=get_sblavg_dict(self.alldf)
        if s_background=='default':
            self.sbldict=self.default_sctrl_slopes
        #COULD RUN CHECK HERE THAT SBLANK CONTROL MATCHES DEFAULT VALS

        self.expdf['measurement_1X']=self.expdf.apply(lambda x:x.measurement/x.relative_dlnfactor,axis=1)
        #self.expdf['re_yield_1X']=self.expdf.apply(calc_re_yield,measure_field='measurement_1X',
        #                            bca_converter=BCA_YldConverter,dns_converter=DNS_YldConverter,axis=1)
        self.expdf['re_yield_1X']=self.expdf.apply(lambda x:x.re_yield/x.relative_dlnfactor,axis=1)
        self.expdf=self.expdf.assign(rxn_yield=np.nan,sbg_yield=np.nan,ebg_yield=np.nan)
#        self.expdf.sbg_yield=self.expdf.apply(lambda x:self.sbldict[x.sname]*x.sconc,axis=1).where(self.expdf.sconc.notnull())

#        self.expdf.sbg_yield=self.expdf.apply(lambda x:self.sbldict[x.sname]*x.sconc,axis=1).where(self.expdf.sconc.notnull())
#        self.expdf.sbg_yield= self.expdf.where(self.expdf.sconc.notnull()).apply(lambda x:self.sbldict[x.sname]*x.sconc if np.isnan(x.sname) else np.nan,axis=1)
        self.expdf.sbg_yield= self.expdf.apply(calc_sbg_yield,args=(self.sbldict,),axis=1)
        self.expdf.ebg_yield= self.expdf.apply(calc_ebg_yield,args=(self.ebl_slope,),axis=1)
        self.expdf.rxn_yield=self.expdf.where(self.expdf.rxn_status).apply(lambda x:x.re_yield_1X-x.sbg_yield-x.ebg_yield,axis=1)


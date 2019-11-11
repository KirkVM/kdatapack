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
    ransac.fit(xtarr,ebdf_bca.measurement)
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
def calc_re_yield(dfrow,bca_converter=None,dns_converter=None):
    if dfrow.calstd_status==True:
        return np.nan
#        self.alldf.re_yield=self.alldf.measurement.where((self.alldf.detection=='DNS') & (self.alldf.calstd_status==False)).\
    assert(dfrow.detection in ['DNS','BCA'])
    if dfrow.detection=='DNS':
        return(dns_converter(dfrow.measurement))
    elif dfrow.detection=='BCA':
        return(bca_converter(dfrow.measurement))
#        self.alldf['re_yield']=self.alldf.apply(calc_re_yield,bca_converter=BCA_YldConverter,dns_converter=DNS_YldConverter)

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
        self.default_calibration_standard={'DNS':{'slope':0.45,'int':0.04}}
        self.default_calibration_standard.update({'BCA':{'slope':2.8,'int':0.17}})
        self.base_bca_predvlp_dlnfactor=5/105

    def calc_yields(self,colname='re_yield',subset_selections={},calibration_standard='default',\
                    s_background='defaualt',e_background='default'):
        assert (calibration_standard in ['default']),\
            "calibration_standard must ==['default']"
        if calibration_standard=='default':
            DNS_YldConverter=build_simple_yldconverter(self.default_calibration_standard['DNS']['int'],\
                                                        self.default_calibration_standard['DNS']['slope'])
            BCA_YldConverter=build_simple_yldconverter(self.default_calibration_standard['BCA']['int'],\
                                                        self.default_calibration_standard['BCA']['slope'])
        #get the econc - yield default converter       
        self.alldf=cleandf(self.alldf)
        self.ebg_converter=calc_bca_ebgconverter(self.alldf)
        #now get a dictionary of default substrate baselines
#        self.sbldict=get_sblavg_dict(self.alldf)

        #now get yields and stuff 
        self.expdf=cleandf(self.expdf)
        self.expdf=self.expdf.assign(re_yield=np.nan)
        #create re_yield column for every entry that isn't a calibration standard
        self.expdf.re_yield=self.expdf.measurement.where((self.expdf.detection=='DNS') & (self.expdf.calstd_status==False)).\
                            apply(DNS_YldConverter)
        self.expdf.re_yield=self.expdf.measurement.where((self.expdf.detection=='BCA') & (self.expdf.calstd_status==False)).\
                            apply(BCA_YldConverter)

        self.alldf['re_yield']=self.alldf.apply(calc_re_yield,bca_converter=BCA_YldConverter,dns_converter=DNS_YldConverter,axis=1)
        
#        self.alldf=self.alldf.assign(re_yield=np.nan)
#        self.alldf.re_yield=self.alldf.measurement.where((self.alldf.detection=='DNS') & (self.alldf.calstd_status==False)).\
#                            apply(DNS_YldConverter)
#        self.alldf.re_yield=self.alldf.measurement.where((self.alldf.detection=='BCA') & (self.alldf.calstd_status==False)).\
#                            apply(BCA_YldConverter)

#        self.expdf=calc_rxn_yields(self.expdf)
import numpy as np
def cleandf(df,overs='nan'):
    cleandf=df.copy()
    if overs=='nan':
        cleandf=cleandf.replace(to_replace={'measurement':{'OVER':np.nan}})
    return cleandf

def build_simple_yldconverter(calint,calslope):#make_cylinder_volume_func(r):
    def yield_converter(measurement):
        return (measurement-calint)/calslope 
    return yield_converter

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
    df=df.assign(econc_mgmL_adjusted=lambda x:x.econc_mgmL*x.relative_dlnfactor)
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
        self.default_calibration_standard={'DNS':{'slope':1.0,'int':0.05}}
        self.default_calibration_standard.update({'BCA':{'slope':1.99,'int':0.15}})
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
        self.expdf=cleandf(self.expdf)
        self.expdf=self.expdf.assign(re_yield=np.nan)
        #create re_yield column for every entry that isn't a calibration standard
        self.expdf.re_yield=self.expdf.measurement.where((self.expdf.detection=='DNS') & (self.expdf.calstd_status==False)).\
                            apply(DNS_YldConverter)
        self.expdf.re_yield=self.expdf.measurement.where((self.expdf.detection=='BCA') & (self.expdf.calstd_status==False)).\
                            apply(BCA_YldConverter)
        #get default bca e_blank:
#        get_dilution_adjusted_values(self.expdf,self.base_bca_predvlp_dlnfactor)
#        get_dilution_adjusted_values(self.alldf,self.base_bca_predvlp_dlnfactor)


#        from sklearn.linear_model import LinearRegression,Lasso,RANSACRegressor

#lf=LinearRegression(fit_intercept=False)
#lasso_f=Lasso()
#lf.fit(xtarr,ebdf.measurement)
#lasso_f.fit(xtarr,ebdf.measurement)
#rlf=LinearRegression(fit_intercept=False)
#ransac=RANSACRegressor(base_estimator=rlf,min_samples=0.9,residual_threshold=0.5)
#ransac.fit(xtarr,ebdf.measurement)
        #calculate bca_eblank_default by
        #   (1) grouping into dilution schemes
        #   (2) get correct to get into one curve with 5/105 and 1
        #   (3) fit line to get avg re_yield vs mgmL_conc equation

        #self.all_eblank_df=get_eblankdf_bca(alldf,converter=BCA_YldConverter)
        self.expdf=calc_rxn_yields(self.expdf)

        #now do some background correction
        #every rxn requires subtraction of eblank,sblank (and potentially ebg and sbg)
#        self.expdf.re_yield_eblank=
#        self.expdf.re_yield_eblank_type=

        # 
#WHAT IS re_yield after correction only for e-only control (0 if DNS)
#        self.expdf=re_yield_ectrlsub
#WHAT IS re_yield after correction only for s-only control
#        self.expdf=re_yield_sctrlsub
#WHAT IS re_yield after correction only for enzyme-solution background (ie wge).  (0 if purified enzyme)
#        self.expdf=re_yield_bgctrlsub



#    """returns a df with average value for each unique sname,sconc,ename,econc"""
#    if controldf is not None:
#        calint,calslope=get_calibration_standard(controldf,'BCA')
#        calibration_standard['BCA']['int']=calint
#        calibration_standard['BCA']['slope']=calslope
#    if bl_subtract=='eavg':
#        dfpctrl=df.dropna(subset=['econc'])
#        dfpctrl=dfpctrl[dfpctrl.sconc.isnull()]
#        pctrl_avg=dfpctrl.measurement.mean()
#        pbl_subtract=pctrl_avg-calibration_standard['BCA']['int']
#        print(f'subtracting {pbl_subtract} from e wells')
#        df=df.assign(measurement_bg=df.measurement)
#        for idx in df[df.econc.notnull()].index:
#            df.loc[idx,'measurement_bg']=df.measurement.loc[idx]-pbl_subtract
##        df.loc[:,'measurement_bg']=[df.measurement.loc[x]-pbl_subtract if x in df.econc.notnull().index]# for x in df.index else df.measurement.loc[x]]
#        df.loc[:,'wyield']=df.measurement_bg.apply(calc_well_yield,\
#                       args=((calibration_standard['BCA']['int'],calibration_standard['BCA']['slope'])) )
#    elif type(bl_subtract)==float:
#        df=df.assign(measurement_bg=df.measurement)
#        for idx in df[df.econc.notnull()].index:
#            df.loc[idx,'measurement_bg']=df.measurement.loc[idx]-bl_subtract
#        #df.loc[:,'measurement_bg']=df.measurement-bl_subtract#(pctrl_avg-calibration_standard['BCA']['int'])
#        df.loc[:,'wyield']=df.measurement_bg.apply(calc_well_yield,\
#                       args=((calibration_standard['BCA']['int'],calibration_standard['BCA']['slope'])) )
#
#    else:
#        df.loc[:,'wyield']=df.measurement.apply(calc_well_yield,\
#                       args=((calibration_standard['BCA']['int'],calibration_standard['BCA']['slope'])) )
#    return df
##
#
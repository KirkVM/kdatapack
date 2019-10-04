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

def assign_rxn_types(df):
    df=df.assign(ectrl_status=lambda x: (x.ename.isnull()==False) & (x.sname.isnull()) & (x.standardname.isnull()),
                 sctrl_status=lambda x: (x.sname.isnull()==False) & (x.ename.isnull()) & (x.standardname.isnull()),
                 calstd_status=lambda x: (x.standardname.isnull()==False) & (x.sname.isnull()) & (x.ename.isnull()),
                 bgctrl_status=lambda x: (x.sname.isnull()) & (x.standardname.isnull()) & (x.standardname.isnull()),
                 rxn_status=lambda x: (x.sname.isnull()==False) & (x.sconc>0) & (x.ename.isnull()==False) & (x.econc>0))
    return df

class CAZyExperiment:
    def __init__(self,selected_df,alldf=None,expname=None):
        self.expdf=selected_df.copy()
        self.expdf=assign_rxn_types(self.expdf)
        if alldf is not None:
            self.alldf=assign_rxn_types(alldf.copy())
        if expname is not None:
            self.expname=expname
        else:
            self.expname=str(self.expdf.expdate.unique()[0])
        self.default_calibration_standard={'DNS':{'slope':1.0,'int':0.05}}
        self.default_calibration_standard.update({'BCA':{'slope':1.99,'int':0.15}})

    def calc_yields(self,colname='re_yield',subset_selections={},calibration_standard='default'):
        assert (calibration_standard in ['default']),\
            "calibration_standard must ==['default']"
        if calibration_standard=='default':
            DNS_YldConverter=build_simple_yldconverter(self.default_calibration_standard['DNS']['int'],\
                                                        self.default_calibration_standard['DNS']['slope'])
            BCA_YldConverter=build_simple_yldconverter(self.default_calibration_standard['BCA']['int'],\
                                                        self.default_calibration_standard['BCA']['slope'])
#        self.expdf=expdf.measurement.apply(DNS_YldConverter) for x in expdf.where(expdf.detection=='DNS')
        self.expdf=cleandf(self.expdf)
        #self.expdf['re_yield']=self.expdf.apply(lambda x:DNS_YldConverter(x.measurement) if x.detection=='DNS' else
        #                                        BCA_YldConverter(x.measurement),axis='columns')
        self.expdf=self.expdf.assign(re_yield=np.nan)
        self.expdf.re_yield=self.expdf.measurement.where((self.expdf.detection=='DNS') & (self.expdf.calstd_status==False)).\
                            apply(DNS_YldConverter)
        self.expdf.re_yield=self.expdf.measurement.where((self.expdf.detection=='BCA') & (self.expdf.calstd_status==False)).\
                            apply(BCA_YldConverter)
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
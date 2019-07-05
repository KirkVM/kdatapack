import scipy
from scipy import optimize
import numpy as np
from dfitlib import loss_funcs
def calculate_activity_plate(df,std_intercept=.07,std_slope=0.5,relative_dilfactor=1.0,
                            standardname='glucose',standardconc_units='mgmL',
                            rxnvol_units='uL',rxntime_units='hours',econc_units='mgmL'):
    if standardconc_units=='mgmL' and standardname=='glucose':
#        df['yield']=[std_slope*(df.measurement.at[x]-std_intercept) for x in df.index]
        df['yield_mgmL']=[std_slope*(df.measurement.at[x]-std_intercept)*relative_dilfactor \
                            for x in df.index]
        df['yield_uM']=[df.yield_mgmL.at[x]*1e6/180.16 for x in df.index]
        if rxnvol_units=='uL':
            df['yield_umoles']=[df.yield_uM.at[x]*df.rxnvol.at[x]*1e-6 for x in df.index]
        if rxntime_units=='hours':
            df['activity_U']=[df.yield_umoles.at[x]/(df.rxntime.at[x]*60) for x in df.index]
        if econc_units=='mgmL' and rxnvol_units=='uL':
            df['spactivity_Umg']=[df.activity_U.at[x]/(df.econc.at[x]*df.rxnvol.at[x]*1e-3) for x in df.index]
#    return df
#            #for x in df.index]
#        if rxntime
#    df['sp_activity']=[]
#def get_calibration
#class Calibrated_Standard:
#    def __init__(self,df):
#        for expdate in df.expdate.unique():

def get_calibration_standard(df,detection):
    #going to start by splitting to individual dsets then returning average slope,int
    df=df.dropna(subset=['standardname']) #redundant if df is actually controldf
    df=df[df.detection==detection.upper()]
    ints=[]
    slopes=[]
    for expdate in df.expdate.unique():
        minidf=df[df.expdate==expdate]
        minidf=minidf.dropna(subset=['measurement'])
        res=scipy.optimize.minimize(loss_funcs.linearfit_l2,(0,1.1),args=(minidf.standardconc.values,minidf.measurement.values))
        ints.append(res.x[0])
        slopes.append(res.x[1])
    return np.mean(ints),np.mean(slopes)

def calc_well_yield(measuredval,int,slope):
    return (measuredval-int)/slope

def add_wyield(df,controldf=None):
    """returns a df with average value for each unique sname,sconc,ename,econc"""
    calibration_standard={'DNS':{'slope':1.0,'int':0.05}}
    calibration_standard.update({'BCA':{'slope':1.99,'int':0.15}})
    if controldf is not None:
        calint,calslope=get_calibration_standard(controldf,'BCA')
        calibration_standard['BCA']['int']=calint
        calibration_standard['BCA']['slope']=calslope
    df.loc[:,'wyield']=df.measurement.apply(calc_well_yield,\
                       args=((calibration_standard['BCA']['int'],calibration_standard['BCA']['slope'])) )
    return df
#    else:
#        calibration_int,calibration_slope=0.15,1.8
#    rxnyields=(rxndf.measurement-0.15)/1.99
#
#    print(calibration_int,calibration_slope)    
##    groups=df.groupby()


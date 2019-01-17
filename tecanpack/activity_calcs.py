
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

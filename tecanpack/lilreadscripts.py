import datetime
from pathlib import Path
import pandas as pd
from . import tecanio

def readnatetempdep():
    infldr_path=Path('./data')
    infname='1-23_40-50_08.16.2018.xlsx'
    ifpath=infldr_path / infname
    shname='40C'

    #set up wells to read
    wellreads_=[]
    rows=['A','B','C','D','E','F','G','H']
    cols=['1','2','3','4','5','6','7','8','9','10','11','12']
    enames=[]
    rxnphs=[]
    for repnum in range(3):
        enames.extend([str(x) for x in range(1,24)])
        enames.append('wge')
        rxnphs.extend([5,4,5,5,5,5,4,5,5,5,4,4])
        rxnphs.extend([6 for x in range(12)]) #coincidence that last 12 (incl. ctrl) were pH6
    enames.extend([None for x in range(24)])
    rxnphs.extend([6 for x in range(12)]) #I think
    rxnphs.extend([None for x in range(12)]) #I think
    #econcs=[None for x in range(72)]
    snames=['xylan' for x in range(48)]
    snames.extend([None for x in range(24)])
    snames.extend(['xylan' for x in range(12)])
    snames.extend([None for x in range(12)])
    standardconcs=[None for x in range(84)]
    standardconcs.extend([2.5,2.5,2.5,1.0,1.0,1.0,0.5,0.5,0.5,0.0,0.0,0.0])
    predvlp_dlnfactors=[2.0 if (x is not None) else 1.0 for x in enames]
    predvlp_dlnfactors[18]=10.0
    predvlp_dlnfactors[42]=10.0
#    postdvlp_dlnfactors=[1.0 for x in range(96)]
#    predvlp_dlnfactors=[1.0 for x in range(96)]
    defaultsdict={'sname':'xylan','sconc':10,'sconc_units':'mgmL',
                  'expdate':datetime.date(2018,8,16),'expfname':infname,
                  'expsheet':shname,'detection':'BCA','rxnph':6.0,'rxntemp':40,
                  'standardconc_units':'mM','standardname':'glucose'}
#                  'predvlp_dlnfactor':2.0}

    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,
                                      standardconcs=standardconcs,rxnphs=rxnphs,
                                      predvlp_dlnfactors=predvlp_dlnfactors)
    t40df=tecanio.platedf_to_datadf(platedf,wellreads_)
    t40df.measurement.iat[17]*=5
    t40df.measurement.iat[41]*=5
    shname='50C'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=50
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,
                                      standardconcs=standardconcs,rxnphs=rxnphs,
                                      predvlp_dlnfactors=predvlp_dlnfactors)
    t50df=tecanio.platedf_to_datadf(platedf,wellreads_)
    t50df.measurement.iat[17]*=5
    t50df.measurement.iat[41]*=5
    t50df.drop(index=[35,47],inplace=True)
    tdf_1thru23=pd.concat([t40df,t50df],ignore_index=True)

    infname='24-46_40-50_08.23.2018.xlsx'
    ifpath=infldr_path / infname
    shname='40C'
    defaultsdict['expsheet']=shname
    defaultsdict['expdate']=datetime.date(2018,8,23)
    defaultsdict['expfname']=infname
    defaultsdict['rxntemp']=40
    enames=[]
    rxnphs=[]
    for repnum in range(3):
        enames.extend([str(x) for x in range(24,47)])
        enames.append('wge')
    enames.extend([None for x in range(24)])
    rxnphs.extend([6 for x in range(84)])
    rxnphs.extend([None for x in range(84)])
    predvlp_dlnfactors=[2.0 if (x is not None) else 1.0 for x in enames]
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),
                          index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,
                                      standardconcs=standardconcs,rxnphs=rxnphs,
                                      predvlp_dlnfactors=predvlp_dlnfactors)
    t40df=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='50C'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=50
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,
                                      standardconcs=standardconcs,rxnphs=rxnphs,
                                      predvlp_dlnfactors=predvlp_dlnfactors)
    t50df=tecanio.platedf_to_datadf(platedf,wellreads_)
    tdf_24thru46=pd.concat([t40df,t50df],ignore_index=True)

    infname='47-51_40-50_08.30.2018.xlsx'
    ifpath=infldr_path / infname
    shname='40_50'
    defaultsdict['expsheet']=shname
    defaultsdict['expdate']=datetime.date(2018,8,30)
    defaultsdict['expfname']=infname
#    defaultsdict['rxntemp']=40
    enames=[]
    snames=[]
    rxntemps=[]
    standardconcs=[]
    rxnphs=[]
    for tempval in [40,50]:
        rxntemps.extend([tempval for x in range(48)])
        for repnum in range(4):
            enames.extend([str(x) for x in range(47,52)])
            enames.append('wge')
            rxnphs.extend([4,6,6,6,4,6])
        enames.extend([None for x in range(24)])
        rxnphs.extend([6 for x in range(12)])
        rxnphs.extend([None for x in range(12)])
        snames.extend(['xylan' for x in range(12)])
        snames.extend([None for x in range(12)])
        snames.extend(['xylan' for x in range(12)])
        snames.extend([None for x in range(12)])
        standardconcs.extend([None for x in range(36)])
        standardconcs.extend([2.5,2.5,2.5,1.0,1.0,1.0,0.5,0.5,0.5,0.0,0.0,0.0])
    predvlp_dlnfactors=[2.0 if (x is not None) else 1.0 for x in enames]
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),
                          index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,
                                      standardconcs=standardconcs,rxnphs=rxnphs,rxntemps=rxntemps,
                                      predvlp_dlnfactors=predvlp_dlnfactors)
    t4050df=tecanio.platedf_to_datadf(platedf,wellreads_)
    t4050df.drop(index=60,inplace=True)
    tdf=pd.concat([tdf_1thru23,tdf_24thru46,t4050df],ignore_index=True)
    return tdf
    #now change settings to sheet2-ph5
#    shname='5_30'
#    defaultsdict['expsheet']=shname
#    defaultsdict['rxnph']=5.0
#    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
#    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs)
#    ph5df=tecanio.platedf_to_datadf(platedf,wellreads_)
def readnatephdep():
    infldr_path=Path('./data')
    infname='1-23_pH at 30_08.15.2018.xlsx'
    ifpath=infldr_path / infname
    shname='4_30'

    #set up wells to read
    wellreads_=[]
    rows=['A','B','C','D','E','F','G','H']
    cols=['1','2','3','4','5','6','7','8','9','10','11','12']
    enames=[]
    for repnum in range(3):
        enames.extend([str(x) for x in range(1,24)])
        enames.append('wge')
    enames.extend([None for x in range(24)])
    #econcs=[None for x in range(72)]
    snames=['xylan' for x in range(48)]
    snames.extend([None for x in range(24)])
    snames.extend(['xylan' for x in range(12)])
    snames.extend([None for x in range(12)])
    standardconcs=[None for x in range(84)]
    standardconcs.extend([2.5,2.5,2.5,1.0,1.0,1.0,0.5,0.5,0.5,0.0,0.0,0.0])
    predvlp_dlnfactors=[2.0 if (x is not None) else 1.0 for x in enames]
    predvlp_dlnfactors[18]=20.0
    predvlp_dlnfactors[42]=20.0

    defaultsdict={'sname':'xylan','sconc':10,'sconc_units':'mgmL',
                  'expdate':datetime.date(2018,8,15),'expfname':infname,
                  'expsheet':shname,'detection':'BCA','rxnph':4.0,'rxntemp':30,
                  'standardconc_units':'mM','standardname':'glucose'}

    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,
                                    standardconcs=standardconcs,predvlp_dlnfactors=predvlp_dlnfactors)
    ph4df=tecanio.platedf_to_datadf(platedf,wellreads_)
    ph4df.measurement.iat[17]*=10
    ph4df.measurement.iat[41]*=10
    #now change settings to sheet2-ph5
    shname='5_30'
    defaultsdict['expsheet']=shname
    defaultsdict['rxnph']=5.0
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,
                                    standardconcs=standardconcs,predvlp_dlnfactors=predvlp_dlnfactors)
    ph5df=tecanio.platedf_to_datadf(platedf,wellreads_)
    ph5df.measurement.iat[17]*=10
    ph5df.measurement.iat[41]*=10
    #now change settings to sheet3-ph6
    shname='6_30'
    defaultsdict['expsheet']=shname
    defaultsdict['rxnph']=6.0
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,
                                    standardconcs=standardconcs,predvlp_dlnfactors=predvlp_dlnfactors)
    ph6df=tecanio.platedf_to_datadf(platedf,wellreads_)
    ph6df.measurement.iat[17]*=10
    ph6df.measurement.iat[41]*=10
    #now change settings to sheet4-ph7
    shname='7_30'
    defaultsdict['expsheet']=shname
    defaultsdict['rxnph']=7.0
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,
                                    standardconcs=standardconcs,predvlp_dlnfactors=predvlp_dlnfactors)
    ph7df=tecanio.platedf_to_datadf(platedf,wellreads_)
    ph7df.measurement.iat[17]*=10
    ph7df.measurement.iat[41]*=10
    phdf_1thru23=pd.concat([ph4df,ph5df,ph6df,ph7df],ignore_index=True)

    #now switch to 24-46 file
    infname="24-46_pH at 30_08.22.2018.xlsx"
    ifpath=infldr_path / infname
    defaultsdict['expfname']=infname
    #redo enames only
    enames=[]
    for repnum in range(3):
        enames.extend([str(x) for x in range(24,47)])
        enames.append('wge')
    enames.extend([None for x in range(24)])
    predvlp_dlnfactors=[2.0 if (x is not None) else 1.0 for x in enames]

    #sheet1-ph4
    shname='4_30'
    defaultsdict['expsheet']=shname
    defaultsdict['rxnph']=4.0
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,
                                    standardconcs=standardconcs,predvlp_dlnfactors=predvlp_dlnfactors)
    ph4df=tecanio.platedf_to_datadf(platedf,wellreads_)
    #now change settings to sheet2-ph5
    shname='5_30'
    defaultsdict['expsheet']=shname
    defaultsdict['rxnph']=5.0
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,
                                    standardconcs=standardconcs,predvlp_dlnfactors=predvlp_dlnfactors)
    ph5df=tecanio.platedf_to_datadf(platedf,wellreads_)
    #now change settings to sheet3-ph6
    shname='6_30'
    defaultsdict['expsheet']=shname
    defaultsdict['rxnph']=6.0
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,
                                    standardconcs=standardconcs,predvlp_dlnfactors=predvlp_dlnfactors)
    ph6df=tecanio.platedf_to_datadf(platedf,wellreads_)
    #now change settings to sheet4-ph7
    shname='7_30'
    defaultsdict['expsheet']=shname
    defaultsdict['rxnph']=7.0
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,
                                    standardconcs=standardconcs,predvlp_dlnfactors=predvlp_dlnfactors)
    ph7df=tecanio.platedf_to_datadf(platedf,wellreads_)
    phdf_1thru46=pd.concat([phdf_1thru23,ph4df,ph5df,ph6df,ph7df],ignore_index=True)
#    phdf_1thru23=pd.concat([ph4df,ph5df,ph6df,ph7df])


    infname='47-51_pH at 30_08.29.2018.xlsx'
    ifpath=infldr_path / infname
    defaultsdict['expfname']=infname
    enames=[]
    for repnum in range(8): #2 for each row
        enames.extend([str(x) for x in range(47,52)])
        enames.append('wge')
    enames.extend([None for x in range(48)])
    snames=[]
    snames.extend(['xylan' for x in range(24)])
    snames.extend([None for x in range(24)])
    snames.extend(['xylan' for x in range(24)])
    snames.extend([None for x in range(24)])

    standardconcs=[None for x in range(72)]
    for repnum in range(2): #2 for each row
        standardconcs.extend([2.5,2.5,2.5,1.0,1.0,1.0,0.5,0.5,0.5,0.0,0.0,0.0])
    rxnphs=[]
    for repnum in range(4): #2 for each row
        rxnphs.extend([4.0 for x in range(12)])
        rxnphs.extend([5.0 for x in range(12)])
    predvlp_dlnfactors=[2.0 if (x is not None) else 1.0 for x in enames]
    shname='4-5_30'
    defaultsdict['expsheet']=shname
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,rxnphs=rxnphs,
                                    standardconcs=standardconcs,predvlp_dlnfactors=predvlp_dlnfactors)
    ph45df=tecanio.platedf_to_datadf(platedf,wellreads_)
    #
    rxnphs=[]
    for repnum in range(4): #2 for each row
        rxnphs.extend([6.0 for x in range(12)])
        rxnphs.extend([7.0 for x in range(12)])
    shname='6-7_30'
    defaultsdict['expsheet']=shname
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,rxnphs=rxnphs,
                                    standardconcs=standardconcs,predvlp_dlnfactors=predvlp_dlnfactors)
    ph67df=tecanio.platedf_to_datadf(platedf,wellreads_)

    phdf=pd.concat([phdf_1thru46,ph45df,ph67df],ignore_index=True)
    return phdf


def read_all_natedata():
    phdf=readnatephdep()
    tdf=readnatetempdep()
    gh_8_43_df=pd.concat([phdf,tdf],ignore_index=True)
    return gh_8_43_df

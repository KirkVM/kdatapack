import datetime,os,itertools,uuid
from pathlib import Path
import pandas as pd
from . import tecanio

def readevanxg():
    infldr_path=Path(os.environ['HOME'])#/'Box Sync'/'/data')
    infldr_path=infldr_path/ 'Box Sync' /'Data'/'GH5SequenceActivity'/'ActivityDataFiles'
    infname='Batch 1 XGase.xlsx'
    ifpath=infldr_path / infname
    shname='XG'

    wellreads_=[]
    ecycle=itertools.cycle(['ADB80100.1','ADX05734.1','ADI13081.1','AFC30482.1','AFC68970.1','BAJ22272.1',\
            'ADD62401.1','BAB04322.1','ADD61809.1','ADU14847.1','ABD79589.1','ACR12247.1',\
            'ADB44000.1','BAL90206.1','BAE44526.1','AAR65335.1','wge','AEQ16450.1',\
            'ADL53038.1','WP_009983134.1','AAA23221.1','ADK55024.1','ABX76045.1','CAA38693.1'])

    ec_cycle=itertools.cycle([0.001904408,0.000624629,0.005507954,0.006292763,\
                              0.002755381,0.006124524,0.003913738,0.002554118,\
                              0.003227074,0.004569871,0.003543612,0.006437947,\
                              0.00365951,0.004098801,0.004815998,0.005191109,\
                              None,0.025,0.000737335,0.002505248,0.00658687,\
                              0.004138057,0.005165561,0.002537294])



    enames=[next(ecycle) for x in range(72)]
    enames.extend([None for x in range(12)])
    econcs=[next(ec_cycle) for x in range(72)]
    econcs.extend([None for x in range(12)])
    snames=['xyloglucan' for x in range(72)]
    snames.extend([None for x in range(12)])
    standardconcs=[None for x in range(72)]
    for x in range(3):
        standardconcs.extend([0.0,0.25,0.5,1.0])
    egskippy=['E10','E11','E12','F10','F11','F12','G10','G11','G12','H10','H11','H12']

    defaultsdict={'sname':'xylan','sconc':10,'sconc_units':'mgmL','econc_units':'uM',
                  'buffername':'NaPhos','bufferconc':0.1,'bufferconc_units':'M',
                  'expdate':datetime.date(2018,10,24),'expfname':infname,
                  'expsheet':shname,'detection':'DNS','rxnph':6.0,'rxntemp':30,
                  'standardconc_units':'mgmL','standardname':'glucose',
                  'ename':None,'rxntime':3.0,'rxntime_units':'hours',
                  'rxnvol':50,'rxnvol_units':'uL','predvlp_dlnfactor':3.0,
                  'experimenter':'Evan Glasgow','biorepid':uuid.uuid4().hex}#,'postdvlp_dlnfactor':1.0}
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,snames=snames,enames=enames,skipwells=egskippy,\
                                      econcs=econcs,standardconcs=standardconcs,itermethod='bycol_wrap')
    xg1df=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='PASC'
    defaultsdict['expsheet']=shname
    snames=['PASC' for x in range(72)]
    snames.extend([None for x in range(12)])
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,snames=snames,enames=enames,skipwells=egskippy,\
                                      econcs=econcs,standardconcs=standardconcs,itermethod='bycol_wrap')
    p1df=tecanio.platedf_to_datadf(platedf,wellreads_)

#
    #2nd experiment set:
    infname='Batch 2 XGase.xlsx'
    ifpath=infldr_path / infname
    shname='XG'
    defaultsdict['expdate']:datetime.date(2018,10,26)
    defaultsdict['expsheet']:shname
    defaultsdict['expfname']=infname
    defaultsdict['biorepid']=uuid.uuid4().hex

    ecycle=itertools.cycle(['AEV87015.1','ADY35478.1','ABW39344.1','CBL35246.1','AEE44521.1',\
                            'ACZ29221.1','BAL62675.1','ACI18413.1','ADU86909.1','AAA23233.1',\
                            'AEV59736.1','ADL34447.1','CAL91973.1','BAC57893.1','ACE84905.1',\
                            'ADQ18612.1','wge','AAF00074.2','WP_005355457.1','WP_010248927.1',\
                            'ADZ19876.1','AAC37035.1','ADZ22510.1','AAB38548.1'])

    ec_cycle=itertools.cycle([0.001537515,0.001535438,0.002164158,0.000413355,\
                              0.001646774,0.001856983,0.004649103,0.002569101,\
                              0.001534607,0.001878585,0.001631403,0.006155733,\
                              0.001989921,0.004575987,0.002728975,0.00163348,\
                              None,0.000982913,0.005837869,0.003215446,\
                              0.005434686,0.005855934,0.004744237,0.001549147])

    enames=[next(ecycle) for x in range(72)]
    enames.extend([None for x in range(12)])
    econcs=[next(ec_cycle) for x in range(72)]
    econcs.extend([None for x in range(12)])
    snames=['xyloglucan' for x in range(72)]
    snames.extend([None for x in range(12)])

    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,snames=snames,enames=enames,skipwells=egskippy,\
                                      econcs=econcs,standardconcs=standardconcs,itermethod='bycol_wrap')
    xg2df=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='PASC'
    defaultsdict['expsheet']=shname
    snames=['PASC' for x in range(72)]
    snames.extend([None for x in range(12)])
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,snames=snames,enames=enames,skipwells=egskippy,\
                                      econcs=econcs,standardconcs=standardconcs,itermethod='bycol_wrap')
    p2df=tecanio.platedf_to_datadf(platedf,wellreads_)
    egxgpdf1018=pd.concat([xg1df,p1df,xg2df,p2df],ignore_index=True)
    return egxgpdf1018

###TODO: add rxnvol, owner to all except Evan datas
####add pconcs to others
def readkirkxg():
    infldr_path=Path(os.environ['HOME'])#/'Box Sync'/'/data')
    infldr_path=infldr_path/ 'Box Sync' /'Data'/'ssEnzymeAssays'
    infname='061617_ZackTest2.xlsx'
    ifpath=infldr_path / infname
    shname='Control'

    #control plate
    #set up wells to read
    wellreads_=[]
    rows=['A','B','C','D','E','F','G','H']
    cols=['1','2','3','4','5','6','7','8','9','10','11','12']
    snames=['xylan' for x in range(36)]
    snames.extend(['xyloglucan' for x in range(36)])
    snames.extend(['cmc' for x in range(24)])

    defaultsdict={'sname':'xylan','sconc':10,'sconc_units':'mgmL',
                  'buffername':'NaPhos','bufferconc':0.1,'bufferconc_units':'M',
                  'expdate':datetime.date(2017,6,16),'expfname':infname,
                  'expsheet':shname,'detection':'DNS','rxnph':6.0,'rxntemp':50,
                  'standardconc_units':'mgmL','standardname':'glucose',
                  'ename':None,'rxntime':4.0,'rxntime_units':'hours',
                  'predvlp_dlnfactor':3.0,'postdvlp_dlnfactor':5.44444,
                  'experimenter':'Kirk Vander Meulen','biorepid':uuid.uuid4().hex}

    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,snames=snames,prows=rows)
    c50df=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='A50Incubator'
    defaultsdict['expsheet']=shname
    rows=['A','B','C','D','E','F','G','H']
#    ecycle=itertools.cycle(['ACE86198.1','AFC30484.1','WP_009984919.1','AFL98798.1','ACZ54907.1','ADZ84561.1'\
#            'ADU23111.1','AEV46039.1','ADU86903.1','ADL33047.1','WP_009984467.1','wge'])
    ecycle=itertools.cycle(['wge','WP_009984467.1','ADL33047.1','ADU86903.1','AEV46039.1',\
                            'ADU23111.1','ADZ84561.1','ACZ54907.1','AFL98798.1',\
                            'WP_009984919.1','AFC30484.1','ACE86198.1'])
    ec_cycle=itertools.cycle([0.75*x if x is not None else None for x in\
                        [None,.0118,.0067,.0074,.0093,.0073,.0078,.0088,.0093,.0060,.0113,.0109]])
                        #ad-hoc adjustment b/c gel loading may be underestimate
    enames=[next(ecycle) for x in range(72)]
    enames.extend([None for x in range(24)])
    econcs=[next(ec_cycle) for x in range(72)]
    econcs.extend([None for x in range(24)])
    del defaultsdict['ename']
    snames=['xylan' for x in range(36)]
    snames.extend(['xyloglucan' for x in range(36)])
    snames.extend(['cmc' for x in range(24)])
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,snames=snames,enames=enames,econcs=econcs)
    h50df=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='A50'
    defaultsdict['expsheet']=shname
    enames=[]
    econcs=[]
    snames=[]
    standardconcs=[]
    for rnum in range(2):
        enames.extend([next(ecycle) for x in range(36)])
        enames.extend([None for x in range(12)])
        econcs.extend([next(ec_cycle) for x in range(36)])
        econcs.extend([None for x in range(12)])
        snames.extend(['xyloglucan' for x in range(12)])
        snames.extend(['xylan' for x in range(12)])
        snames.extend(['cmc' for x in range(12)])
        snames.extend([None for x in range(12)])
        standardconcs.extend([None for x in range(36)])
        standardconcs.extend([0,.5,1,1.5,2,2.5,3,3.5,4,0,.5,1])
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,snames=snames,enames=enames,standardconcs=standardconcs,econcs=econcs)
    inc50df=tecanio.platedf_to_datadf(platedf,wellreads_)


    shname='A30'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=30
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,snames=snames,enames=enames,standardconcs=standardconcs,econcs=econcs)
    inc30df=tecanio.platedf_to_datadf(platedf,wellreads_)
    kdf061617=pd.concat([c50df,h50df,inc50df,inc30df],ignore_index=True)




    subfldrname='061917Translation'
    infname='062117_phs_4hScreen.xlsx'
    ifpath=infldr_path / subfldrname / infname
    shname='phcontrols'
    defaultsdict['expsheet']=shname
    rows=['B','C','D','E','F','G']
    defaultsdict['rxntemp']=50
    defaultsdict['ename']=None
    defaultsdict['expdate']=datetime.date(2017,6,21)
    defaultsdict['expfname']=infname
    defaultsdict['biorepid']=uuid.uuid4().hex
    snames=['xyloglucan' for x in range(24)]
    snames.extend(['cmc' for x in range(24)])
    snames.extend(['xylan' for x in range(24)])
    rxnphs=[]
    for repnum in range(3):
        rxnphs.extend([4 for x in range(12)])
        rxnphs.extend([8 for x in range(12)])
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,prows=rows,snames=snames,rxnphs=rxnphs)
    cphdf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='ph_XG_CMC'
    defaultsdict['expsheet']=shname
    ecycle=itertools.cycle(['AGM71677.1','CBL34359.1','ACZ98591.1','AAA23224.1','EIM57503.1',\
                            'ABS61403.1','WP_004082283.1','BAL45505.1','ADK97059.1',\
                            'CBL16523.1','ADD61911.1','wge'])

    ec_cycle=itertools.cycle([.01016,.0076,.0103,.001,.0055,.001,.001,.001,.001,.001,.0019,None])
    enames=[next(ecycle) for x in range(96)]
    econcs=[next(ec_cycle) for x in range(96)]
    snames=['xyloglucan' for x in range(48)]
    snames.extend(['cmc' for x in range(48)])
    rxnphs=[]
    for repnum in range(2):
        rxnphs.extend([4 for x in range(12)])
        rxnphs.extend([5 for x in range(12)])
        rxnphs.extend([7 for x in range(12)])
        rxnphs.extend([8 for x in range(12)])
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,rxnphs=rxnphs,econcs=econcs)
    phvdf=tecanio.platedf_to_datadf(platedf,wellreads_)

    infname='062117_305070_4hScreen.xlsx'
    ifpath=infldr_path / subfldrname / infname
    shname='30C'
    defaultsdict['expsheet']=shname
    rows=['B','C','D','E','F']
    defaultsdict['rxntemp']=30
    defaultsdict['expfname']=infname
    enames=[next(ecycle) for x in range(36)]
    enames.extend([None for x in range(24)])
    econcs=[next(ec_cycle) for x in range(36)]
    econcs.extend([None for x in range(24)])
    snames=['xyloglucan' for x in range(12)]
    snames.extend(['xylan' for x in range(12)])
    snames.extend(['cmc' for x in range(12)])
    snames.extend([None for x in range(24)])
    standardconcs=[None for x in range(36)]
    standardconcs.extend([0,.5,1,1.5,2,2.5,3,3.5,4,0,.5,1])
    standardconcs.extend([0,.5,1,1.5,2,2.5,3,3.5,4,0,.5,1])
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,prows=rows,econcs=econcs)
    t30h4df=tecanio.platedf_to_datadf(platedf,wellreads_)


    shname='50C'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=50
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,prows=rows,econcs=econcs)
    t50h4df=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='70C'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=70
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,prows=rows,econcs=econcs)
    t70h4df=tecanio.platedf_to_datadf(platedf,wellreads_)

    infname='DNS_0_19_A_G10.xlsx'
    ifpath=infldr_path / subfldrname / infname
    shname='30degreeplate'
    defaultsdict['expsheet']=shname
    defaultsdict['expfname']=infname
    defaultsdict['rxntime']=20
    defaultsdict['rxntemp']=30
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,prows=rows,econcs=econcs)
    t30ondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='50degreeplate'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=50
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,prows=rows,econcs=econcs)
    t50ondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='70degreeplate'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=70
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,prows=rows,econcs=econcs)
    t70ondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='ControlPlate'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=50
    snames=['xyloglucan' for x in range(36)]
    snames.extend(['xylan' for x in range(36)])
    snames.extend(['cmc' for x in range(24)])
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,snames=snames)
    c50ondf=tecanio.platedf_to_datadf(platedf,wellreads_)
    #cphdf,phvdf,t30h4df,t50h4df,t70h4df,t30ondf,t50ondf,t70ondf,c50ondf
    kdf062117=pd.concat([cphdf,phvdf,t30h4df,t50h4df,t70h4df,t30ondf,t50ondf,t70ondf,c50ondf],ignore_index=True)

    #######################################################
    ############### 7/12/17 translation####################
    #######################################################
    subfldrname='071217Translation'
    infname='DNS_112G_7_14_17.xlsx'
    ifpath=infldr_path / subfldrname / infname
    shname='Control Plate'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=50
    defaultsdict['rxntime']=20
    defaultsdict['rxnph']=6
    defaultsdict['ename']=None
    defaultsdict['expdate']=datetime.date(2017,7,13)
    defaultsdict['expfname']=infname
    defaultsdict['biorepid']=uuid.uuid4().hex
    snames=[None for x in range(12)]
    snames.extend(['xylan' for x in range(24)])
    snames.extend(['cmc' for x in range(24)])
    snames.extend(['xyloglucan' for x in range(24)])
    snames.extend([None for x in range(12)])
    standardconcs=[0,.5,1,1.5,2,2.5,3,3.5,4,4,5,6]
    standardconcs.extend([None for x in range(72)])
    standardconcs.extend([0,.5,1,1.5,2,2.5,3,3.5,4,4,5,6])
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,snames=snames,standardconcs=standardconcs)
    cdf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='30'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=30
    defaultsdict['rxntime']=20
    ecycle=itertools.cycle(['CAN91059.1','ADX05688.1','AGB29444.1','ABD81896.1','ABX76046.1',\
                            'ADU21423.1','ACZ98591.1','AEV67086.1','ADD82896.1','BAA92146.1',\
                            'CBL16772.1','wge','CBX94583.1','AFC28691.1','AET58208.1','ADL44291.1',\
                            'CBL17363.1','AEV98717.1','ABQ03808.1','ADJ49096.1','ABW39343.1',\
                            'ABW39339.1','ABW39338.1','wge'])
    ec_cycle=itertools.cycle([.0057,.0118,.0058,.0045,.0051,.0036,.0050,.0055,.0052,.005,.0061,None,\
                            .0056,.006,.0048,.0057,.0046,.0048,.006,.0115,.0089,.007,.0106,None])
    enames=[None for x in range(12)]
    enames.extend([next(ecycle) for x in range(72)])
    enames.extend([None for x in range(12)])
    econcs=[None for x in range(12)]
    econcs.extend([next(ec_cycle) for x in range(72)])
    econcs.extend([None for x in range(12)])
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t30hondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='50'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=50
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t50hondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='70'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=70
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t70hondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='pH 4'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=50
    defaultsdict['rxnph']=4
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t50ph4hondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='pH 5'
    defaultsdict['expsheet']=shname
    defaultsdict['rxnph']=5
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t50ph5hondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='pH 7'
    defaultsdict['expsheet']=shname
    defaultsdict['rxnph']=7
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t50ph7hondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='pH 8'
    defaultsdict['expsheet']=shname
    defaultsdict['rxnph']=8
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t50ph8hondf=tecanio.platedf_to_datadf(platedf,wellreads_)


    infname='071317_305070_4hScreen.xlsx'
    ifpath=infldr_path / subfldrname / infname
    defaultsdict['expfname']=infname
    shname='30'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=30
    defaultsdict['rxntime']=4
    defaultsdict['rxnph']=6
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t30h4df=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='50'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=50
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t50h4df=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='70'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=70
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t70h4df=tecanio.platedf_to_datadf(platedf,wellreads_)

    kdf071217=pd.concat([cdf,t30hondf,t50hondf,t70hondf,t50ph4hondf,t50ph5hondf,t50ph7hondf,t50ph8hondf,t30h4df,t50h4df,t70h4df],ignore_index=True)
    #
    #######################################################
    ############### 7/19/17 translation####################
    #######################################################
    subfldrname='071917Translation'
    infname='072117_onScreens.xlsx'
    ifpath=infldr_path / subfldrname / infname
    shname='pH4ctrl'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=50
    defaultsdict['rxntime']=20
    defaultsdict['rxnph']=4
    defaultsdict['ename']=None
    defaultsdict['expdate']=datetime.date(2017,7,20)
    defaultsdict['expfname']=infname
    defaultsdict['biorepid']=uuid.uuid4().hex
    #snames, standardconcs same as 7/12 but keeping code here for clarity
    snames=[None for x in range(12)]
    snames.extend(['xylan' for x in range(24)])
    snames.extend(['cmc' for x in range(24)])
    snames.extend(['xyloglucan' for x in range(24)])
    snames.extend([None for x in range(12)])
    standardconcs=[0,.5,1,1.5,2,2.5,3,3.5,4,4,5,6]
    standardconcs.extend([None for x in range(72)])
    standardconcs.extend([0,.5,1,1.5,2,2.5,3,3.5,4,4,5,6])
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,snames=snames,standardconcs=standardconcs)
    c4df=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='pH6ctrl'
    defaultsdict['expsheet']=shname
    defaultsdict['rxnph']=6
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,snames=snames,standardconcs=standardconcs)
    c6df=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='pH8ctrl'
    defaultsdict['expsheet']=shname
    defaultsdict['rxnph']=8
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,snames=snames,standardconcs=standardconcs)
    c8df=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='30C'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=30
    defaultsdict['rxntime']=20
    defaultsdict['rxnph']=6
    ecycle=itertools.cycle(['CAJ70720.1','AEV59735.1','AAC19169.1','AAC97596.1','CAJ19139.1',\
                            'CBX94583.1','AEV98717.1','AFY52522.1','ADX05688.1','ADL44291.1',\
                            'AAA23220.1','wge','ABD81896.1','CAH69214.1','CBL35063.1','AAA23224.1',\
                            'ACE86198.1','BAL45505.1','CBL16523.1','AFJ44730.1','CBL35076.1',\
                            'AFC30484.1','WP_009984919.1','wge'])
    ec_cycle=itertools.cycle([.005,.0089,.0047,.0056,.0023,.0058,.0053,.0043,.0047,.0041,.0078,None,\
                             .0065,.0054,.0056,.0063,.0067,.0051,.0034,.0050,.0137,.0101,.0101,None])
    enames=[None for x in range(12)]
    enames.extend([next(ecycle) for x in range(72)])
    enames.extend([None for x in range(12)])
    econcs=[None for x in range(12)]
    econcs.extend([next(ec_cycle) for x in range(72)])
    econcs.extend([None for x in range(12)])
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t30hondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='50C'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=50
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t50hondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='70C'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=70
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t70hondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='pH4'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=50
    defaultsdict['rxnph']=4
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t50ph4hondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='pH5'
    defaultsdict['expsheet']=shname
    defaultsdict['rxnph']=5
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t50ph5hondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='pH7'
    defaultsdict['expsheet']=shname
    defaultsdict['rxnph']=7
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t50ph7hondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='pH8'
    defaultsdict['expsheet']=shname
    defaultsdict['rxnph']=8
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t50ph8hondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    infname='072017_305070_4hScreen.xlsx'
    ifpath=infldr_path / subfldrname / infname
    defaultsdict['expfname']=infname
    shname='30'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=30
    defaultsdict['rxntime']=4
    defaultsdict['rxnph']=6
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t30h4df=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='50'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=50
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t50h4df=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='70'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=70
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,econcs=econcs)
    t70h4df=tecanio.platedf_to_datadf(platedf,wellreads_)

    kdf071917=pd.concat([c4df,c6df,c8df,t30hondf,t50hondf,t70hondf,\
                        t50ph4hondf,t50ph5hondf,t50ph7hondf,t50ph8hondf,t30h4df,t50h4df,t70h4df],\
                        ignore_index=True)
    #c4df,c6df,c8df,t30hondf,t50hondf,t70hondf,t50ph4hondf,t50ph7ondf,t50ph8ondf,t30h4df,t50h4df,t70h4df

    #######################################################
    ############### 7/25/17 translation####################
    #######################################################
    subfldrname='072517Translation'
    infname='DNS_112A_A3_H1_112B_G2_H11_overnight.xlsx'
    ifpath=infldr_path / subfldrname / infname
    shname='pH 5&7'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=50
    defaultsdict['rxntime']=20
    defaultsdict['expdate']=datetime.date(2017,7,26)
    defaultsdict['expfname']=infname
    defaultsdict['biorepid']=uuid.uuid4().hex
    #snames, standardconcs same as 7/12 but keeping code here for clarity

    ecycle=itertools.cycle(['ACA61140.1','ADX05742.1','ADR64666.1','BAL45505.1','CCK24845.1',\
                            'AEB43836.1','ACU73736.1','ADL52801.1','AEV98714.1','CCG00653.1',\
                            'ACX74396.1','wge'])
    ec_cycle=itertools.cycle([.001,.0089,.001,.0085,.010,.008,.011,.008,.0045,.009,.008,None])
    enames=[None for x in range(12)]
    enames.extend([next(ecycle) for x in range(72)])
    enames.extend([None for x in range(12)])
    econcs=[None for x in range(12)]
    econcs.extend([next(ec_cycle) for x in range(72)])
    econcs.extend([None for x in range(12)])
    snames=[None for x in range(12)]
    snames.extend(['xylan' for x in range(12)])
    snames.extend(['cmc' for x in range(12)])
    snames.extend(['xyloglucan' for x in range(12)])
    snames.extend(['xylan' for x in range(12)])
    snames.extend(['cmc' for x in range(12)])
    snames.extend(['xyloglucan' for x in range(12)])
    snames.extend([None for x in range(12)])
    rxnphs=[5 for x in range(48)]
    rxnphs.extend([7 for x in range(48)])
    standardconcs=[0,.5,1,1.5,2,2.5,3,3.5,4,4,5,6]
    standardconcs.extend([None for x in range(72)])
    standardconcs.extend([0,.5,1,1.5,2,2.5,3,3.5,4,4,5,6])
    rows=['A','B','C','D','E','F','H']
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,rxnphs=rxnphs,snames=snames,standardconcs=standardconcs,econcs=econcs)
    ph57df=tecanio.platedf_to_datadf(platedf,wellreads_)
    #because ran out of protein in row G:
    ph57df.drop(index=[72,73,74,75,76,77,78,79,80,81,82,83],inplace=True)

    shname='pH 4&8'
    defaultsdict['expsheet']=shname
    rxnphs=[8 for x in range(48)]
    rxnphs.extend([8 for x in range(48)])
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,rxnphs=rxnphs,snames=snames,standardconcs=standardconcs,econcs=econcs)
    ph48df=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='30'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=30
    defaultsdict['rxnph']=6
    rows=['B','C','D','E','F']
    standardconcs=[0,.5,1,1.5,2,2.5,3,3.5,4,4,5,6]
    standardconcs.extend([None for x in range(36)])
    standardconcs.extend([0,.5,1,1.5,2,2.5,3,3.5,4,4,5,6])
    snames=[None for x in range(12)]
    snames.extend(['xylan' for x in range(12)])
    snames.extend(['cmc' for x in range(12)])
    snames.extend(['xyloglucan' for x in range(12)])
    snames.extend([None for x in range(12)])
    enames=[None for x in range(12)]
    enames.extend([next(ecycle) for x in range(36)])
    enames.extend([None for x in range(12)])
    econcs=[None for x in range(12)]
    econcs.extend([next(ec_cycle) for x in range(36)])
    econcs.extend([None for x in range(12)])
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,prows=rows,econcs=econcs)
    t30hondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='50'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=50
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,prows=rows,econcs=econcs)
    t50hondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='70'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=70
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,prows=rows,econcs=econcs)
    t70hondf=tecanio.platedf_to_datadf(platedf,wellreads_)

    infname='DNS_112A_A3_H1_112B_G2_H11.xlsx'
    ifpath=infldr_path / subfldrname / infname
    defaultsdict['expfname']=infname
    shname='30'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=30
    defaultsdict['rxntime']=4
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,prows=rows,econcs=econcs)
    t30h4df=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='50'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=50
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,prows=rows,econcs=econcs)
    t50h4df=tecanio.platedf_to_datadf(platedf,wellreads_)

    shname='70'
    defaultsdict['expsheet']=shname
    defaultsdict['rxntemp']=70
    platedf=pd.read_excel(ifpath,skiprows=list(range(22)),nrows=8,usecols=range(13),index_col=0,sheet_name=shname)
    wellreads_=tecanio.describe_wells(defaultsdict,enames=enames,snames=snames,standardconcs=standardconcs,prows=rows,econcs=econcs)
    t70h4df=tecanio.platedf_to_datadf(platedf,wellreads_)

    kdf072517=pd.concat([ph57df,ph48df,t30hondf,t50hondf,t70hondf,t30h4df,t50h4df,t70h4df],ignore_index=True)
    kdf=pd.concat([kdf061617,kdf062117,kdf071217,kdf071917,kdf072517],ignore_index=True)

    return kdf
#    return kdf061617,kdf062117,kdf071217,kdf071917,kdf072517

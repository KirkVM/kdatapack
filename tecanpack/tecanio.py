import datetime,itertools,re,uuid
import pandas as pd
import dataclasses
from dataclasses import dataclass

@dataclass
class WellReading:
    '''Base class for 96-well plate reading'''
    measurement: float = None
    expdate: datetime.date = None
    expfname: str = None
    expsheet: str = None
    volume: float = None
    detection: str = None #dns, bca, etc
    wellid: str = None #A1-H12
    predvlp_dlnfactor: float = 1.0
    postdvlp_dlnfactor: float = 1.0
    experimenter: str = None
    absorbance_wl: int = None #3/5/19---add these as defaults according to dns/bca/etc
    excitation_wl: int = None #3/5/19---add these as defaults according to dns/bca/etc
    src_vesselid: str = None #reflects whether this measurement is unique to a source vessel
                             #this will by default get a new uuid.uuid4().hex
                             #somehow incorporate multiple reads from same well intelligence
#somehow incorporate batches/freshness?
#add
@dataclass
class WellRxnReading(WellReading):
    '''Base class for well involving reaction or incubation.
    Overrides WellReading and adds rxn and solultion details'''
    rxntime: float = None
    rxntime_units: str = None
    rxnvessel: str = None
    rxnvessel_wellid: str = None #in case different from wellid of reading
    incubator: str = None
    rxntemp: float = None
    rxnph: float = None
    rxnvol: float = None
    rxnvol_units: str = None
    buffername: str = None
    bufferconc: float = None
    bufferconc_units: str = None
    solutionstr: str = None
    biorepid: str = None #uuid (hex-rep) reflecting uniquely prepared biosamples

@dataclass
class WellEnzymeSubstrateReading(WellRxnReading):
    '''Enzyme-Substrate reaction.
    Overrides WellRxnReading, adding enzyme and substrate names & concs'''
    ename: str = None
    econc: float = None
    econc_units: str = None
    sname: str = None
    sconc: float = None
    sconc_units: str = None
    etype: str = None #'cellfree","purified"
    ename_type: str = "acc" #default is that the name=accession code

@dataclass
class WellEnzymeControl(WellRxnReading):
    '''Enzyme-only control well.
    Overrides WellRxnReading, adding enzyme name & conc'''
    ename: str = None
    econc: float = None
    econc_units: str = None
    etype: str = None #'cellfree","purified"
    ename_type: str = "acc" #default is that the name=accession code

@dataclass
class WellSubstrateControl(WellRxnReading):
    '''Substrate-only control well.
    Overrides WellRxnReading, adding substrate name & conc'''
    sname: str = None
    sconc: float = None
    sconc_units: str = None

@dataclass
class WellStandard(WellRxnReading):
    '''Standard well for calibration of some reaction.
    Overrides WellRxnReading, adding standard name & conc'''
    standardname: str = None
    standardconc: float = None
    standardconc_units: float = None


def platedf_to_datadf(platedf,wreadings_):
    """
    Takes 'dumb' 96well plate dataframe + list of WellReadings, returns a tidy dataframe

    Arguments:
    platedf: raw plate data (matrix with row letter as index, column number as column heading)
    wreadings_: list of WellReadings (or subclass of WellReading)

    Returns:
    dataframe containing all possible WellReading fields, filled in where appropriate
    """
    #get all possible fields related to WellReading
    allfields=[dataclasses.fields(x) for x in
                [WellReading,WellRxnReading,WellEnzymeSubstrateReading,
                WellEnzymeControl,WellSubstrateControl,WellStandard]]
    allfnames=[x.name for x in itertools.chain.from_iterable(allfields)]

    rowcolRE=re.compile('([A-H])(\d{1,2})')
    #set up target dataframe columns
    welldatadf=pd.DataFrame(columns=list(set(allfnames)))
    #build temporary dataframe from wreadings_ then add to target dataframe
    newrows_=[]
    for wread in wreadings_:
        update_dict=dataclasses.asdict(wread)
        newrows_.append(pd.Series(update_dict))
    newdf=pd.DataFrame(newrows_)
    welldatadf=welldatadf.append(newdf,ignore_index=True,sort=False)
    #add measurements from platedf according to row-col
    for ridx in welldatadf.index:
        rowstr,colstr=rowcolRE.match(welldatadf.wellid.at[ridx]).groups()
        welldatadf['measurement'].at[ridx]=platedf[int(colstr)].at[rowstr]
    return welldatadf


def describe_wells(defaultsdict,
                   prows=['A','B','C','D','E','F','G','H'],
                   pcols=[1,2,3,4,5,6,7,8,9,10,11,12],
                   skipwells=[],itermethod='byrow_wrap',
                   enames=None,econcs=None,snames=None,sconcs=None,
                   standardnames=None,standardconcs=None,rxnphs=None,
                   rxntemps=None,rxntimes=None,buffernames=None,
                   predvlp_dlnfactors=None,postdvlp_dlnfactors=None,
                   repeat_measurements=None):
    """
    Takes lists of well settings and returns a list of WellReadings or subclasses as appropriate

    defaultsdict={'sname':'xylan','sconc':10,'sconc_units':'mgmL',
                  'expdate':datetime.date(2018,8,15),'expfname':infname,
                  'expsheet':shname,'detection':'BCA','rxnph':4.0,'rxntemp':30,
                  'standardconc_units':'mM','standardname':'glucose'}
    Arguments:
    defaultsdict: various fields relevant to WellReadings and subclasses that apply to all wells

    Keyword arguments:
    prows: list of row letters in the list (default ['A'-'H'])
    pcols: list of integer column values (default [1-12])
    skipwells: list of str rowcol values to not include (eg. ['A11','C7']) (default [])
    itermethod: direction/rule to iterate through rows/cols in parallel with passing over list
                (default 'byrow_wrap' = 'A1','A2',...,'B1','B2',...'H12')
    enames: list of str enzyme name vals in each well (default = None)
    econcs: list of float enzyme conc vals in each well (default = None)
    snames: list of str substrate name vals in each well (default = None)
    sconcs: list of float substrate conc vals in each well (default = None)
    standardnames: list of str standard name vals in each well (default = None)
    standardconcs: list of float standard conc vals in each well (default = None)
    rxnphs: list of float rxn ph vals (default = None)
    rxntemps: list of float rxn temp vals (default = None)
    rxntimes: list of float rxn time vals (default = None)
    buffernames: list of str buffer names (default = None)
    predvlp_dlnfactors: list of dilution factor before develop step (default = None)
    postdvlp_dlnfactors: list of dilution factor post develop step (default = None)

    Returns:
    List of WellReadings or appropriate subclass
    """
    wellreads_=[] #this is the list we'll build
    if itermethod not in ['byrow_wrap','bycol_wrap']:
        print('itermethod not implemented yet..')
        return None
    #defaulting to byrow_wrap-
    if itermethod=='byrow_wrap':
        allwids_=[''.join(x) for x in itertools.product(prows,[str(x) for x in pcols])]
    elif itermethod=='bycol_wrap':
        allwids_=[''.join(x[::-1]) for x in itertools.product([str(x) for x in pcols],prows)]
    wids_=[x for x in allwids_ if x not in skipwells]
    #assign src_vesselids
    svrepsdict={}
    if repeat_measurements is not None:
        for rmwellids in repeat_measurements:
            newsvid=uuid.uuid4().hex
            svrepsdict.update({x:newsvid for x in rmwellids})
    #iterate through wids_, build a WellReading for each--
    #use defaultsdict but override with one of optional list arguments if present
    #build a well_dict containing all values for each well, then use it to build WellReading object
    for widx,wid in enumerate(wids_):
        well_dict=defaultsdict.copy()
        well_dict['wellid']=wid
        if wid in svrepsdict.keys():
            well_dict['src_vesselid']=svrepsdict[wid]
        else:
            well_dict['src_vesselid']=uuid.uuid4().hex
        if enames is not None:
            well_dict['ename']=enames[widx]
        if econcs is not None:
            well_dict['econc']=econcs[widx]
        if snames is not None:
            well_dict['sname']=snames[widx]
        if sconcs is not None:
            well_dict['sconc']=sconcs[widx]
        if standardconcs is not None:
            well_dict['standardconc']=standardconcs[widx]
        if rxnphs is not None:
            well_dict['rxnph']=rxnphs[widx]
        if rxntemps is not None:
            well_dict['rxntemp']=rxntemps[widx]
        if rxntimes is not None:
            well_dict['rxntime']=rxntimes[widx]
        if buffernames is not None:
            well_dict['buffernames']=buffernames[widx]
        if predvlp_dlnfactors is not None:
            well_dict['predvlp_dlnfactor']=predvlp_dlnfactors[widx]
        if postdvlp_dlnfactors is not None:
            well_dict['postdvlp_dlnfactor']=postdvlp_dlnfactors[widx]
        #now determine what type of WellReading subclass this is based on keys in well_dict
        if well_dict['ename'] is not None and well_dict['sname'] is not None:
            reading=WellEnzymeSubstrateReading(**{x:well_dict[x] for x in
                             [y.name for y in dataclasses.fields(WellEnzymeSubstrateReading)]
                             if x in well_dict.keys()})
        elif well_dict['ename'] is not None and well_dict['sname'] is None:
            reading=WellEnzymeControl(**{x:well_dict[x] for x in
                             [y.name for y in dataclasses.fields(WellEnzymeControl)]
                             if x in well_dict.keys()})
        elif well_dict['sname'] is not None and well_dict['ename'] is None:
            reading=WellSubstrateControl(**{x:well_dict[x] for x in
                             [y.name for y in dataclasses.fields(WellSubstrateControl)]
                             if x in well_dict.keys()})
        else:
            if well_dict['standardconc'] is not None:
                reading=WellStandard(**{x:well_dict[x] for x in
                             [y.name for y in dataclasses.fields(WellStandard)]
                             if x in well_dict.keys()})
            else:
                print('well not classified')
        wellreads_.append(reading)
        #should check and correct for dilution factor here.........

    return wellreads_

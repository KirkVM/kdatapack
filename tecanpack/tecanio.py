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

def read_tecan_excel(ifpath,shname,plate_header_row=23,plate_header_col="A"):
    skiprows=list(range(plate_header_row-1))
    if plate_header_col=='A':
        usecols=range(13)
    else:
        raise NotImplementedError("column offset other than A not yet implemented")
    platedf=pd.read_excel(ifpath,skiprows=skiprows,nrows=8,usecols=usecols,index_col=0,sheet_name=shname)
    #check structure
    assert (list(platedf.index)==['A','B','C','D','E','F','G','H'] and \
            list(platedf.columns)!=['1','2','3','4','5','6','7','8','9','10','11','12']), \
            "excel df read-in indexing incorrect. Is row/col set correctly?"
    return platedf

class TecanWell:
    def __init__(self,wellreading):
        self.wellreading=wellreading

class TecanPlate:
    def __init__(self,defaultsdict):
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
        """
        #wellids=[x+y for x in ['A','B','C','D','E','F','G','H']
        #     for y in ['1','2','3','4','5','6','7','8','9','10','11','12']]
        #get all possible fields related to WellReading

        #previous...
        self.all_settings=[dataclasses.fields(x) for x in
                   [WellReading,WellRxnReading,WellEnzymeSubstrateReading,
                    WellEnzymeControl,WellSubstrateControl,WellStandard]]
        self.defaultsdict={x:None for x in self.all_settings}

        #new...        
#        self.defaultsdict={x.name:x.default for x in dataclasses.fields(y) for y in 
#                            [WellReading,WellRxnReading,WellEnzymeSubstrateReading,
#                            WellEnzymeControl,WellSubstrateControl,WellStandard]}
#        self.all_settings=self.defaultsdict.keys()

        self.defaultsdict.update(defaultsdict)
        self.plateid=uuid.uuid4().hex
        self.active_wids=None
        self.welldict={}#dictionary with wellreadings values
        self.exceldf=None#df pulled straight from excel sheet
        #go ahead and set up welldatadf here...
        self.welldatadf=pd.DataFrame(columns=list(set(self.all_settings)))
    def read_excel_sheet(self,ifpath,shname):
        self.exceldf=read_tecan_excel(ifpath,shname)
        #automatically create if create_layout has already run through
        if len(self.welldict)>0:
            self.make_welldatadf()

    def create_layout(self,active_rows=['A','B','C','D','E','F','G','H'],
                      active_cols=[1,2,3,4,5,6,7,8,9,10,11,12],
                      skipwells=[],repeat_measurements=[],itermethod='by_row',
                      enames=None,econcs=None,snames=None,sconcs=None,standardnames=None,standardconcs=None,
                      rxnphs=None,rxntemps=None,rxntimes=None,buffernames=None,
                      predvlp_dlnfactors=None,postdvlp_dlnfactors=None):
        if itermethod not in ['by_row','by_column']:
            raise NotImplementedError('Only byrow_wrap and bycol_wrap iterations across plate are implemneted')
        if itermethod=='by_row':
            allwids_=[''.join(x) for x in itertools.product(active_rows,[str(x) for x in active_cols])]
        elif itermethod=='by_column':
            allwids_=[''.join(x[::-1]) for x in itertools.product([str(x) for x in active_cols],active_rows)]
        self.active_wids=[x for x in allwids_ if x not in skipwells]
        #now let's check
        for x in [enames,econcs,snames,sconcs,standardnames,standardnames,standardconcs,
                  rxnphs,rxntemps,rxntimes,buffernames,predvlp_dlnfactors,postdvlp_dlnfactors]:
            if x is not None:
                assert (len(x)==len(self.active_wids)), \
                   "length of a varied condition (enames,rxnphs,etc) does not match # of active wells"
        #assign src_vesselids
        svrepsdict={}
        if len(repeat_measurements)>0:
            for rmwellids in repeat_measurements:
                newsvid=uuid.uuid4().hex
                svrepsdict.update({x:newsvid for x in rmwellids})
    #iterate through wids_, build a WellReading for each--
    #use defaultsdict but override with one of optional list arguments if present
    #build a well_dict containing all values for each well, then use it to build WellReading object
        for widx,wid in enumerate(self.active_wids):
            well_dict=self.defaultsdict.copy()
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
            self.welldict[wid]=TecanWell(reading)
            #automatically create if read_excel_sheet has already run through
            if self.exceldf is not None:
                self.make_welldatadf()

    def make_welldatadf(self):
        assert (self.active_wids is not None), 'plate layout not yet created (call .create_layout())'
        assert (self.welldict[x] is not None for x in self.active_wids), 'not all active_wids created. (call .create_layout())'
        assert (self.exceldf is not None), 'read in the excel file before extracting df (call .read_excel_sheet())'
        #build temporary dataframe from self.welldict then add to target dataframe
        dfrows_=[]
        for wid in self.active_wids:
            dict2add=dataclasses.asdict(self.welldict[wid].wellreading)
            dfrows_.append(pd.Series(dict2add))
        newdf=pd.DataFrame(dfrows_)
        rowcolRE=re.compile('([A-H])(\d{1,2})')
        self.welldatadf=self.welldatadf.append(newdf,ignore_index=True,sort=False)
        #for rcidx in self.welldict.keys():
        #for rcidx in self.welldict.keys():
        for dfidx in self.welldatadf.index:
            #rowstr,colstr=rowcolRE.match(rcidx).groups()
            rowstr,colstr=rowcolRE.match(self.welldatadf.wellid.at[dfidx]).groups()
#            self.welldict[rcidx]['measurement'].at[rcidx]=self.exceldf[int(colstr)].at[rowstr]
            xlval=self.exceldf[int(colstr)].at[rowstr]
            self.welldict[rowstr+colstr].measurement=xlval
            self.welldatadf['measurement'].at[dfidx]=xlval#self.welldatadf.append(newdf,ignore_index=True,sort=False)

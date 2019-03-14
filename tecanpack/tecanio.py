import datetime,itertools,re,uuid
import pandas as pd
import dataclasses
from dataclasses import dataclass

@dataclass
class WellReading:
    '''96-well plate reading'''
    measurement: float = None
    expdate: datetime.date = None
    expfname: str = None
    expsheet: str = None
    msrvolume: float = None
    detection: str = None #dns, bca, etc
    wellid: str = None #A1-H12
    predvlp_dlnfactor: float = None
    postdvlp_dlnfactor: float = None
    experimenter: str = None
    absorbance_wl: int = None #3/5/19---add these as defaults according to dns/bca/etc
    excitation_wl: int = None #3/5/19---add these as defaults according to dns/bca/etc
    plateid: str = None 
#somehow incorporate batches/freshness?
    '''Settings for reaction/incubation. All about the setup pre-developing''' 
    rxntime: float = None
    rxntime_units: str = None
    rxnvessel: str = None
    rxnvesselid: str = None #this will by default get a new uuid.uuid4().hex
    rxnvessel_wellid: str = None #in case different from wellid of reading
    rxn_incubator: str = None
    rxn_rpm: int = None
    rxntemp: float = None
    rxnph: float = None
    rxnvol: float = None
    rxnvol_units: str = None
    rxn_description_string: str = None
    buffername: str = None
    bufferconc: float = None
    bufferconc_units: str = None
    solution_description_string: str = None
    biorepid: str = None #uuid (hex-rep) reflecting uniquely prepared biosamples
    #substrate settings
    sname: str = None
    sconc: float = None
    sconc_units: str = None
    sbarcode: str = None
    sprepdate: datetime.date = None
    sstock_conc: float = None
    sstock_conc_units: str = None
    s_description_string: str = None
    #enzyme settings
    ename: str = None
    econc: float = None
    econc_units: str = None
    ebarcode: str = None
    eprepdate: datetime.date = None
    epreptype: str = None #'cellfree","purified"
    enametype: str = "acc" #default is that the name=accession code
    evariant: str = None #ko, ea, +cbmx2,etc
    ethawdate: datetime.date = None
    estock_conc: float = None
    estock_conc_units: str = None
    e_description_string: str = None
    #standard settings
    standard_name: str = None
    standard_conc: float = None
    standard_conc_units: str = None
    standard_stock_conc: float = None
    standard_stock_conc_units: float = None
    standard_description_string: str = None

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

def make_plates_rxn_duplicates(plate1,plate2):
    assert (len(set(plate1.welldict.keys()).difference(plate2.welldict.keys()))==0),\
        "cannot make entire plate rxn duplicates. don't have same wells"
    for dfidx in plate2.welldatadf.index:
        wid=plate2.welldatadf.loc[dfidx].wellid
        rxnvesselid2keep=plate1.welldict[wid].wellreading.rxnvesselid
        plate2.welldict[wid].wellreading.rxnvesselid=rxnvesselid2keep #set TecanWell
        plate2.welldatadf.loc[dfidx,'rxnvesselid']=rxnvesselid2keep #set df also

class TecanWell:
    def __init__(self,well_settings_dict,use_experimenter_defaults=True):
        #can add some error handling here...
        self.well_settings_dict=well_settings_dict
        if use_experimenter_defaults:
            self.set_experimenter_defaults()
        
        self.wellreading=WellReading(**{x:well_settings_dict[x] for x in well_settings_dict.keys()})
           #now determine what type of WellReading subclass this is based on keys in well_dict
    def set_experimenter_defaults(self):
        #this could be redone with a dictionary much nicer...
        assert (self.well_settings_dict['experimenter'].upper() in ['KIRK VANDER MEULEN','ERIC HENNEMAN']),\
                f"experimenter {self.well_settings_dict['experimenter']} has no default settings"
        if 'DNS'==self.well_settings_dict['detection'].upper():
            if 'absorbance_wl' not in self.well_settings_dict.keys():
                self.well_settings_dict['absorbance_wl']=540
            if 'predvlp_dlnfactor' not in self.well_settings_dict.keys():
                if self.well_settings_dict['experimenter'].upper() in ['KIRK VANDER MEULEN','ERIC HENNEMAN']:
                    self.well_settings_dict['predvlp_dlnfactor']=1/3
            if 'postdvlp_dlnfactor' not in self.well_settings_dict.keys():
                if self.well_settings_dict['experimenter'].upper() in ['KIRK VANDER MEULEN','ERIC HENNEMAN']:
                    self.well_settings_dict['postdvlp_dlnfactor']=36/196
            if 'msrvolume' not in self.well_settings_dict.keys():
                if self.well_settings_dict['experimenter'].upper() in ['KIRK VANDER MEULEN','ERIC HENNEMAN']:
                    self.well_settings_dict['msrvolume']=196
        elif 'BCA'==self.well_settings_dict['detection'].upper():
            if 'absorbance_wl' not in self.well_settings_dict.keys():
                self.well_settings_dict['absorbance_wl']=562
            if 'predvlp_dlnfactor' not in self.well_settings_dict.keys():
                if self.well_settings_dict['experimenter'].upper() in ['KIRK VANDER MEULEN','ERIC HENNEMAN']:
                    self.well_settings_dict['predvlp_dlnfactor']=5/105
            if 'postdvlp_dlnfactor' not in self.well_settings_dict.keys():
                if self.well_settings_dict['experimenter'].upper() in ['KIRK VANDER MEULEN','ERIC HENNEMAN']:
                    self.well_settings_dict['postdvlp_dlnfactor']=1
            if 'msrvolume' not in self.well_settings_dict.keys():
                if self.well_settings_dict['experimenter'].upper() in ['KIRK VANDER MEULEN','ERIC HENNEMAN']:
                    self.well_settings_dict['msrvolume']=90
        if 'HEIDOLPH'==self.well_settings_dict['rxn_incubator'].upper():
            self.well_settings_dict['rxn_rpm']=900
        elif 'THERMOCYCLER'==self.well_settings_dict['rxn_incubator'].upper():
            self.well_settings_dict['rxn_rpm']=0


 
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
        """
       #wellids=[x+y for x in ['A','B','C','D','E','F','G','H']
        #     for y in ['1','2','3','4','5','6','7','8','9','10','11','12']]
        #get all possible fields related to WellReading

        #previous...
        self.defaultsdict=defaultsdict
        self.plateid=uuid.uuid4().hex
        self.defaultsdict['plateid']=self.plateid
        self.active_wids=None
        self.welldict={}#dictionary with wellreadings values
        self.exceldf=None#df pulled straight from excel sheet
        #go ahead and set up welldatadf here...
        self.welldatadf=None
    def read_excel_sheet(self,ifpath,shname,plate_header_row=23,plate_header_col='A'):
        """
        Reads excel sheet into the class

        Arguments:
            ifpath: full path and filename of excel file
            shname: sheet name in excel file 
        
        Keyword arguments:
            kwargs to pass.... (needs implementation) 
        """
        self.exceldf=read_tecan_excel(ifpath,shname,plate_header_row=plate_header_row,
                                        plate_header_col=plate_header_col)
        #automatically create if create_layout has already run through
        if len(self.welldict)>0:
            self.make_welldatadf()

    def create_layout(self,active_rows=['A','B','C','D','E','F','G','H'],
                      active_cols=[1,2,3,4,5,6,7,8,9,10,11,12],
                      skipwells=[],repeat_measurements=[],itermethod='by_row',
                      enames=None,econcs=None,snames=None,sconcs=None,standardnames=None,standardconcs=None,
                      rxnphs=None,rxntemps=None,rxntimes=None,buffernames=None,
                      predvlp_dlnfactors=None,postdvlp_dlnfactors=None):
        """
        Keyword arguments:
        active_rows: list of row letters in the list (default ['A'-'H'])
        active_cols: list of integer column values (default [1-12])
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
 
        if itermethod not in ['by_row','by_column']:
            raise NotImplementedError('Only byrow_wrap and bycol_wrap iterations across plate are implemneted')
        if itermethod=='by_row':
            allwids_=[''.join(x) for x in itertools.product(active_rows,[str(x) for x in active_cols])]
        elif itermethod=='by_column':
            allwids_=[''.join(x[::-1]) for x in itertools.product([str(x) for x in active_cols],active_rows)]
        self.active_wids=[x for x in allwids_ if x not in skipwells]
        #now let's check
        varied_dict={'ename':enames,'econc':econcs,'sname':snames,'sconc':sconcs,
                     'standardname':standardnames,'standardconc':standardconcs,
                     'rxnph':rxnphs,'rxntemp':rxntemps,'rxntime':rxntimes,'buffername':buffernames,
                     'predvlp_dlnfactor':predvlp_dlnfactors,'postdvlp_dlnfactor':postdvlp_dlnfactors}
        for ckey in list(varied_dict.keys()):
            if varied_dict[ckey] is None:
                del varied_dict[ckey]
            else:
                assert (len(varied_dict[ckey])==len(self.active_wids)), \
                   f"length of a varied condition {ckey} does not match # of active wells"

        #assign src_vesselids
        rvrepsdict={}
        if len(repeat_measurements)>0:
            for rmwellids in repeat_measurements:
                newrvid=uuid.uuid4().hex
                rvrepsdict.update({x:newrvid for x in rmwellids})
    #iterate through wids_, build a WellReading for each--
    #use defaultsdict but override with one of optional list arguments if present
    #build a well_dict containing all values for each well, then use it to build WellReading object
        for widx,wid in enumerate(self.active_wids):
            well_settings_dict=self.defaultsdict.copy()
            well_settings_dict['wellid']=wid
            if wid in rvrepsdict.keys():
                well_settings_dict['rxnvesselid']=rvrepsdict[wid]
            else:
                well_settings_dict['rxnvesselid']=uuid.uuid4().hex
            well_settings_dict.update({ckey:varied_dict[ckey][widx] for ckey in varied_dict})
            self.welldict[wid]=TecanWell(well_settings_dict)
           #should add exception handling here
           #automatically create if read_excel_sheet has already run through
        if self.exceldf is not None:
            self.make_welldatadf()

    def make_welldatadf(self):
        assert (self.active_wids is not None), 'plate layout not yet created (call .create_layout())'
        assert (self.welldict[x] is not None for x in self.active_wids), 'not all active_wids created. (call .create_layout())'
        assert (self.exceldf is not None), 'read in the excel file before extracting df (call .read_excel_sheet())'
        #build temporary dataframe from self.welldict then add to target dataframe
        self.welldatadf=pd.DataFrame(columns=[x.name for x in dataclasses.fields(WellReading)])
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

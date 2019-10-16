import datetime,itertools,re,uuid
import pandas as pd
import dataclasses
from dataclasses import dataclass
from bokeh.plotting import figure,show,reset_output,ColumnDataSource,output_notebook
from bokeh.models.glyphs import Text
from bokeh.models import FixedTicker,HoverTool,Range1d

class TecanSet:
    '''Utility class returned by load_tecandata. list of plates (.plates) is only property

    Methods:
        get_df(): returns merged dataframe
        get_plate(): returns a particular plate (typical use is diagnosing a problem with loadscript)
    '''
    def __init__(self,plates):
        '''TecanSet constructor, requires list of plates''' 
        self.plates=plates
    def get_df(self):
        '''returns pandas dataframe for all included plates'''
        return pd.concat([x.welldatadf for x in self.plates],ignore_index=True)
    def get_plate(self,plateid=None):
        '''returns a single plate from TecanSet
        
        Keyword arguments:
        plateid: plateid value associated with plate of interest 
                (this is currently a required keyword argument)
        '''
        for plate in self.plates:
            if plate.plateid==plateid:
                rplate=plate
                break
        return rplate    


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
    econc_molar: float = None
    econc_mgmL: float = None
    emw: float = None #molecular weight in Da
    ebarcode: str = None
    eprepdate: datetime.date = None
    epreptype: str = None #"cellfree","ecoli_purified","cellfree_purified"
    enametype: str = None #("acc","jgi_shorthand","kirk_shorthand","evan_shorthand","nate_shorthand")#default is that the name=accession code
    evariant: str = None #ko, ea, +cbmx2,etc
    ethawdate: datetime.date = None
    estock_conc: float = None
    estock_conc_units: str = None
    e_description_string: str = None
    egbacc: str = None
    ename_shorthand : str = None
    epurified_status : bool = None

    #standard settings
    standardname: str = None
    standardconc: float = None
    standardconc_units: str = None
    standardstock_conc: float = None
    standardstock_conc_units: float = None
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
            list(platedf.columns)==[1,2,3,4,5,6,7,8,9,10,11,12]), \
            "excel df read-in indexing incorrect. Is row/col set correctly?"
    return platedf

def make_plates_rxn_duplicates(plate1,plate2,allow_subplate=False):
    '''Set a plate's rxnvesselid values to those of previous plate
    Iterates over new plate (plate2) and copies from previous based on wellid
    
    Arguments:
    plate1: existing plate to copy rxnvesselid values from
    plate2: plate that has been added
    '''
    if allow_subplate==False:
        assert (len(set(plate1.welldict.keys()).difference(plate2.welldict.keys()))==0),\
             "cannot make entire plate rxn duplicates. don't have same wells"
    else:
        assert (len(set(plate2.welldict.keys()).difference(plate1.welldict.keys()))==0),\
             "cannot make entire plate rxn duplicates. don't have same wells"
    for dfidx in plate2.welldatadf.index:
        wid=plate2.welldatadf.loc[dfidx].wellid
        rxnvesselid2keep=plate1.welldict[wid].wellreading.rxnvesselid
        plate2.welldict[wid].wellreading.rxnvesselid=rxnvesselid2keep #set TecanWell
        plate2.welldatadf.loc[dfidx,'rxnvesselid']=rxnvesselid2keep #set df also

class TecanWell:
    def __init__(self,well_settings_dict,use_experimenter_defaults=False):
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
                #print('CAREFUL-- using experimenter defaults for predvlp_dlnfactor')
                if self.well_settings_dict['experimenter'].upper() in ['KIRK VANDER MEULEN','ERIC HENNEMAN']:
                    self.well_settings_dict['predvlp_dlnfactor']=1/3
            if 'postdvlp_dlnfactor' not in self.well_settings_dict.keys():
                #print('CAREFUL-- using experimenter defaults for postdvlp_dlnfactor')
                if self.well_settings_dict['experimenter'].upper() in ['KIRK VANDER MEULEN','ERIC HENNEMAN']:
                    self.well_settings_dict['postdvlp_dlnfactor']=36/196
            if 'msrvolume' not in self.well_settings_dict.keys():
                if self.well_settings_dict['experimenter'].upper() in ['KIRK VANDER MEULEN','ERIC HENNEMAN']:
                    self.well_settings_dict['msrvolume']=196
        elif 'BCA'==self.well_settings_dict['detection'].upper():
            if 'absorbance_wl' not in self.well_settings_dict.keys():
                self.well_settings_dict['absorbance_wl']=562
            if 'predvlp_dlnfactor' not in self.well_settings_dict.keys():
                #print('CAREFUL-- using experimenter defaults for predvlp_dlnfactor')
                if self.well_settings_dict['experimenter'].upper() in ['KIRK VANDER MEULEN','ERIC HENNEMAN']:
                    self.well_settings_dict['predvlp_dlnfactor']=5/105
            if 'postdvlp_dlnfactor' not in self.well_settings_dict.keys():
                #print('CAREFUL-- using experimenter defaults for postdvlp_dlnfactor')
                if self.well_settings_dict['experimenter'].upper() in ['KIRK VANDER MEULEN','ERIC HENNEMAN']:
                    self.well_settings_dict['postdvlp_dlnfactor']=1
            if 'msrvolume' not in self.well_settings_dict.keys():
                if self.well_settings_dict['experimenter'].upper() in ['KIRK VANDER MEULEN','ERIC HENNEMAN']:
                    self.well_settings_dict['msrvolume']=90
        if 'rxn_incubator' in self.well_settings_dict.keys():
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
        #previous...
        self.defaultsdict=defaultsdict
        self.plateid=uuid.uuid4().hex
        self.defaultsdict['plateid']=self.plateid
        self.active_wids=None
        self.welldict={}#dictionary with wellreadings values
        self.exceldf=None#df pulled straight from excel sheet
        #go ahead and set up welldatadf here...
        self.welldatadf=None
        self.ifpath=None
        self.expsheet=None
    def read_excel_sheet(self,ifpath,shname,plate_header_row=23,plate_header_col='A'):
        """
        Reads excel sheet into the class

        Arguments:
            ifpath: full path and filename of excel file
            shname: sheet name in excel file 
        
        Keyword arguments:
            kwargs to pass.... (needs implementation) 
        """
        self.ifpath=ifpath
        self.expsheet=shname
        self.exceldf=read_tecan_excel(ifpath,shname,plate_header_row=plate_header_row,
                                        plate_header_col=plate_header_col)
        #automatically create if create_layout has already run through
        if len(self.welldict)>0:
            self.make_welldatadf()

    def create_layout(self,active_rows=['A','B','C','D','E','F','G','H'],
                      active_cols=[1,2,3,4,5,6,7,8,9,10,11,12],
                      skipwells=[],repeat_measurements=[],itermethod='by_row',
                      enames=None,econcs=None,snames=None,sconcs=None,standardnames=None,standardconcs=None,
                      rxnphs=None,rxntemps=None,rxntimes=None,buffernames=None,rxnvols=None,
                      predvlp_dlnfactors=None,postdvlp_dlnfactors=None,sbarcodes=None):
        """
        Keyword arguments:
        active_rows: list of row letters in the list (default ['A'-'H'])
        active_cols: list of integer column values (default [1-12])
        skipwells: list of str rowcol values to not include (eg. ['A11','C7']) (default [])
        repeat_measurements: str rowcol values that represent repeat measurements of same assay/tube
                            list of lists- eg [['A1','B1'],['H7','G7','B3]]
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
        rxnvols: list of float rxnvols (default = None)
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
                     'rxnvol':rxnvols,
                     'predvlp_dlnfactor':predvlp_dlnfactors,'postdvlp_dlnfactor':postdvlp_dlnfactors,
                     'sbarcode':sbarcodes}
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
        for dfidx in self.welldatadf.index:
            rowstr,colstr=rowcolRE.match(self.welldatadf.wellid.at[dfidx]).groups()
            xlval=self.exceldf[int(colstr)].at[rowstr]
            self.welldict[rowstr+colstr].measurement=xlval
            self.welldatadf['measurement'].at[dfidx]=xlval#self.welldatadf.append(newdf,ignore_index=True,sort=False)
    def assess_cleanup_disambiguate(self,enamedf=None,snamedf=None):
        '''this runs after welldatadf to add molecular weight info etc'''
        #PART 1: CHECKS
        #start by ensuring that all necessary fields are present & logic is correct
        possible_epreptype_values=['cellfree','ecoli_purified','cellfree_purified']
        possible_enametype_values=['gbacc','jgi_shorthand','kirk_shorthand','evan_shorthand','nate_shorthand']
        possible_evariant_values=['KO','EA','CBMX2']
        possible_econc_units_values=['uM','mM','mgmL']
        assert(self.welldatadf.epreptype.dropna().isin(possible_epreptype_values).all()),\
                f"epreptype value must be one of {possible_epreptype_values}. Plate {self.ifpath.name}-{self.expsheet}"
        assert(self.welldatadf.enametype.dropna().isin(possible_enametype_values).all()),\
                f"enametype value must be one of {possible_enametype_values}. Plate {self.ifpath.name}-{self.expsheet}"
        assert(self.welldatadf.evariant.dropna().isin(possible_evariant_values).all()),\
                f"evariant value must be one of {possible_evariant_values}. Plate {self.ifpath.name}-{self.expsheet}"
        assert(self.welldatadf.econc_units.dropna().isin(possible_econc_units_values).all()),\
                f"econc_units value must be one of {possible_econc_units_values}. Plate {self.ifpath.name}-{self.expsheet}"
        #make sure proper fields present in all cases
        #two-way requirements---
        for rtype in [['sconc','sname'],['econc','ename'],['standardname','standardconc']]:
            assert(self.welldatadf[self.welldatadf[rtype[0]].isnull()][rtype[1]].dropna().shape[0]==0),\
                    f"{rtype[1]} value is set for a well with no {rtype[0]}. Plate {self.ifpath.name}-{self.expsheet}"
            assert(self.welldatadf[self.welldatadf[rtype[1]].isnull()][rtype[0]].dropna().shape[0]==0),\
                    f"{rtype[0]} value is set for a well with no {rtype[1]}. Plate {self.ifpath.name}-{self.expsheet}"
        #one-way requirements---
        for rtype in [['ename','epreptype'],['ename','enametype'],['econc','econc_units']]:
            assert(self.welldatadf[self.welldatadf[rtype[0]].notna()][rtype[1]].notna().all()),\
                f"{rtype[1]} is required for all wells with a {rtype[0]}. Plate {self.ifpath.name}-{self.expsheet}"

        #PART 2: FIXED LOGIC - add purified status to enzymes...currently only purified if epreptype=="ecoli_purified"
        self.welldatadf=self.welldatadf.assign(epurified_status=lambda x:x.epreptype=='ecoli_purified')
        sname=[]
        for snidx in self.welldatadf.index:
            cur_sname=self.welldatadf.loc[snidx,'sname']
            if cur_sname is None:
                sname.append(None)
            else:
                sname.append(cur_sname.lower())
        self.welldatadf=self.welldatadf.assign(sname=sname)
#        self.welldatadf.sname=self.welldatadf.sname.apply(str.lower)
        #PART 3: ENAMEDF
        if enamedf is not None:
            egbacc=[]
            emw=[]
            econc_molar=[]
            econc_mgmL=[]
            for dfidx in self.welldatadf.index:
                cur_ename=self.welldatadf.loc[dfidx,'ename']
                if cur_ename is None: 
                    egbacc.append(None)
                    emw.append(None)
                    econc_molar.append(None)
                    econc_mgmL.append(None)
                    continue
                cur_enametype=self.welldatadf.loc[dfidx,'enametype']
                #enamerow=enamedf[enamedf[cur_enametype].str.match(cur_ename,case=False)]
                enamerow=enamedf[enamedf[cur_enametype].str.lower()==cur_ename.lower()]
                assert(enamerow.shape[0]==1),'should be 1ao1 ename match'
                egbacc.append(enamerow['gbacc'].values[0])
                emw.append(enamerow['mw'].values[0])
                #now update conc's using that molecular weight and econc + econc_units
                cur_econc=self.welldatadf.loc[dfidx,'econc']
                cur_econc_units=self.welldatadf.loc[dfidx,'econc_units']
                if cur_econc_units=='uM':
                    econc_molar.append(cur_econc*1e-6)
                    econc_mgmL.append(cur_econc*1e-6*enamerow['mw'].values[0])
                if cur_econc_units=='mM':
                    econc_molar.append(cur_econc*1e-3)
                    econc_mgmL.append(cur_econc*1e-3*enamerow['mw'].values[0])
                if cur_econc_units=='mgmL':
                    if enamerow['mw'].notna().all():
                        econc_molar.append(cur_econc/enamerow['mw'].values[0])
                        econc_mgmL.append(cur_econc)#*1e-3*enamerow['mw'].values[0])
                    else:
                        econc_molar.append(None)
                        econc_mgmL.append(None)#*1e-3*enamerow['mw'].values[0])
            self.welldatadf=self.welldatadf.assign(egbacc=egbacc,emw=emw,econc_mgmL=econc_mgmL,econc_molar=econc_molar)
        #PART 4: SNAMEDF
        if snamedf is not None:
            all_possible_snames=[x.lower() for x in list(snamedf.allnames) for x in x]
            assert set(list(self.welldatadf.sname.dropna().str.lower().unique())).issubset(all_possible_snames),\
                   "one more more sname values doesn't exist in current sname file"
            sname=[]
            for dfidx in self.welldatadf.index:
                cur_sname=self.welldatadf.loc[dfidx,'sname']
                if cur_sname is None: 
                    sname.append(None)
                else:
                    for snidx in snamedf.index:
                        cur_allnames=snamedf.loc[snidx,'allnames']
                        if cur_sname in [x.lower() for x in cur_allnames]:
                            sname.append(snamedf.loc[snidx,'sname'])

            self.welldatadf['full_sname']=self.welldatadf.sname
            self.welldatadf=self.welldatadf.assign(sname=sname)

#    
#        assert()
#                f". Plate {self.ifpath.name}-{self.expsheet}"
#        assert()
#                f". Plate {self.ifpath.name}-{self.expsheet}"
#        assert()
#                f". Plate {self.ifpath.name}-{self.expsheet}"
#        assert()
#                f". Plate {self.ifpath.name}-{self.expsheet}"
        


#        if edb is not None:
#            conn=sqlite3.connect(edb)
#            self.e_accession=....

    def view_plate(self,inline_jupyter=True):
        reset_output()
        hover_tool=HoverTool(names=['stuff'],tooltips=[('Well','@wellid'),('Enzyme','@enzyme'),('Substrate','@substrate'),('Standard','@standard')]
                    )
        pscale=90
        p=figure(plot_width=12*pscale,plot_height=8*pscale,tools=[hover_tool])
        rows=['A','B','C','D','E','F','G','H']
        cols=[1,2,3,4,5,6,7,8,9,10,11,12]
        ecpairs=[]
        for rnum,rname in enumerate(rows):
            for cnum,cname in enumerate(cols):
                curwell=f'{rname}{cname}'
                xpos0,xpos1=cnum*pscale,(cnum+1)*pscale
                ypos0,ypos1=(8-rnum)*pscale,(8-rnum+1)*pscale
                if curwell in self.welldatadf.wellid.values:
                    welldf=self.welldatadf[self.welldatadf.wellid==curwell]
                    cds=get_plthovers(welldf)
                    ecpair=(welldf.ename.values[0],welldf.sname.values[0])
                    if welldf.ename.values[0] is not None and welldf.sname.values[0] is not None:
                        rcolor='purple'
                    elif welldf.ename.values[0] is not None:# and welldf.sname.values[0] is not None:
                        rcolor='blue'
                    elif welldf.sname.values[0] is not None:# and welldf.sname.values[0] is not None:
                        rcolor='red'
                    elif welldf.standardname.values[0] is not None:
                        rcolor='yellow'
                    if ecpair not in ecpairs:
                        ecpairs.append(ecpair)
                    ecidx=ecpairs.index(ecpair)+1
                    ecidxcds=ColumnDataSource(data={'ecidx':[str(ecidx)]})
                    gt=Text(x=xpos0-0.3*pscale,y=ypos0+0.05*pscale,text='ecidx')
                    p.add_glyph(ecidxcds,gt)#x=xpos0-pscale*0.05,y=ypos0-0.05,text=str(ecidx))
                    p.rect(x=xpos0-pscale*0.05,y=ypos0-0.05*pscale,color=rcolor,width=0.9*pscale,height=0.9*pscale,alpha=0.1,name='stuff',source=cds)#  rnum*100,width=100,y=cnum*100,width=100
                else:
                    p.rect(x=xpos0-pscale*0.05,y=ypos0-0.05*pscale,width=0.9*pscale,height=0.9*pscale,alpha=0,line_width=1,line_alpha=1,line_color='gray')#  rnum*100,width=100,y=cnum*100,width=100
        p.grid.grid_line_width=0
        p.xaxis.axis_line_width=0
        p.yaxis.axis_line_width=0
        xaxis_ticks_tuple=[(int(pscale*x),cols[x]) for x in range(len(cols))]
        p.xaxis.ticker=FixedTicker(ticks=[x[0] for x in xaxis_ticks_tuple])
        p.xaxis.major_label_overrides={x:str(y) for x,y in xaxis_ticks_tuple }
        p.xaxis.major_tick_in=0
        p.xaxis.major_tick_out=0
        p.xaxis.major_label_text_font_size=f'{int(14*pscale/100)}pt'

        yaxis_ticks_tuple=[(int(pscale*(8-x)),rows[x]) for x in range(len(rows))]
        p.yaxis.ticker=FixedTicker(ticks=[x[0] for x in yaxis_ticks_tuple])
        p.yaxis.major_label_overrides={x:str(y) for x,y in yaxis_ticks_tuple }
        p.yaxis.major_tick_in=0
        p.yaxis.major_tick_out=0
        p.yaxis.major_label_text_font_size=f'{int(16*pscale/100)}pt'
        #p.yaxis.major_label_standoff=int(-30*pscale/100)#-20#f'{int(14*pscale/100)}pt'
        p.y_range=Range1d(0.4*pscale,pscale*(9-0.4))
        p.x_range=Range1d(-(0.6*pscale),pscale*(12-0.4))
        if inline_jupyter:
            output_notebook()
        show(p)

def get_plthovers(welldfrow):
    cdict={}
    cdict['wellid']=welldfrow['wellid']
    if welldfrow.ename is not None:
        cdict['enzyme']=[f'{welldfrow["ename"].values[0]}: {welldfrow["econc"].values[0]} {welldfrow["econc_units"].values[0]}']
    else:
        cdict['enzyme']=['none']
    if welldfrow.sname is not None:
        cdict['substrate']=[f'{welldfrow["sname"].values[0]}: {welldfrow["sconc"].values[0]} {welldfrow["sconc_units"].values[0]}']
    else:
        cdict['substrate']=['none']
    if welldfrow.standardname.values[0] is not None:
        cdict['standard']=[f'{welldfrow["standardname"].values[0]}: {welldfrow["standardconc"].values[0]} {welldfrow["standardconc_units"].values[0]}']
    else:
        cdict['standard']=['none']
    return ColumnDataSource(data=cdict)
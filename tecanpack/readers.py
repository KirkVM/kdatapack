import sqlite3,atexit
import yaml,os,uuid,datetime
from pathlib import Path
import importlib.util
import pandas as pd

def determinepath(curpath,cfdict):
    if cfdict['filepath_prefix_type']=='relative':
        thepath=curpath
    elif cfdict['filepath_prefix_type']=='envar':
        thepath=Path(os.environ[cfdict['filepath_prefix']])
    elif cfdict['filepath_prefix_type']=='absolute':
        thepath=Path(cfdict['filepath_prefix'])
    thepath= thepath / Path(cfdict['filepath_ending'])
    return thepath

def read_load_confile(cfpathstr,refresh_all):
    cfp=Path(cfpathstr)
    cfdir=cfp.parent
    with open(cfp,'r') as f:
        configgy=yaml.load(f)
    #logdb first:
    logdb_dict=configgy['logdb']
    logdbpath=determinepath(cfdir,logdb_dict)
#    ldbname=logdb_dict['dbname']
    
    pkldf_dict=configgy['pkldfs']
    pkldf_fldrpath=determinepath(cfdir,pkldf_dict)

    load_scripts=[]
    for ldentry in configgy['load_scripts']:
        load_dict=ldentry['load_script']
        lspath=determinepath(cfdir,load_dict)
        ls_funcname=load_dict['funcname']
        load_scripts.append([lspath,ls_funcname])

    conn=sqlite3.connect(logdbpath)
    atexit.register(conn.close)
    conn.row_factory=sqlite3.Row
    c=conn.cursor()

    if refresh_all:
        c.execute('''DROP TABLE PLATES''')
        for pkldf in os.listdir(pkldf_fldrpath):
            os.remove(os.path.join(pkldf_fldrpath,pkldf))
    c.execute('''CREATE TABLE IF NOT EXISTS PLATES (plateid text, filename text, sheetname text,
                expdate text, experimenter text, scriptname text, pkl_ctime text, pkldfname text, pklplate glob)''')

    dfs=[]
    for ls in load_scripts:
        load_required=False
        lsentry=str(ls[0].name)+':'+ls[1]
        seltuple=(lsentry,)
        try:
            c.execute('''SELECT * FROM lscripts WHERE scriptname=?''',seltuple)
            prev_entry=c.fetchone()
            prev_name=prev_entry['scriptname']
            prev_load_time=prev_entry['pkl_ctime']
            prev_pklname=prev_entry['pkldfname']
            prev_picklepath=pkldf_fldrpath / prev_pklname
            #compare create time with load_script's mod time
            ls_mtime=os.path.getmtime(ls[0])
            if ls_mtime>float(prev_load_time):
                load_required=True
        except:
            print(f'reloading df using {lsentry}')
            load_required=True
        if load_required:
            spec=importlib.util.spec_from_file_location(str(ls[0]),ls[0])#'kdfs/dload_dir/load_scripts/testscript.py','kdfs/dload_dir/load_scripts/testscript.py')#
            themodule=importlib.util.module_from_spec(spec)
            spec.loader.exec_module(themodule)
            m2c=getattr(themodule,ls[1])
            curdf=m2c()
            pklname=f'pkldf-{uuid.uuid4().hex}'
            pklpath=pkldf_fldrpath / pklname
            curdf.to_pickle(pklpath)
            curdftime=os.path.getmtime(pklpath)
            intuple=(lsentry,curdftime,pklname)
            c.execute('''INSERT INTO lscripts VALUES (?,?,?)''',intuple)
        else:
         #   curdf=pd.module_from_spec
            curdf=pd.read_pickle(prev_picklepath)
        dfs.append(curdf)


#    edf=pd.concat([dns_std_df,bca_std_df,df0305dns,df0305bca],ignore_index=True)
            #c.execute('''CREATE)
            #print(curdf.shape)
    conn.commit()
    conn.close()
    return dfs
 

def load_tecandata(cfpathstr,refresh_all=False):
    '''reads excel files using load scripts/functions as specified in a config file

    Arguments:
        cfpathstr: path to yaml config file
    Keyword Arguments:
        refresh_all: whether to reset all files for a fresh re-read (default False)
    '''
    dfs=read_load_confile(cfpath,refresh_all)
    mergedf=pd.concat([x for x in dfs],ignore_index=True)
    return mergedf

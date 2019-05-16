import sqlite3,atexit,pickle
import yaml,os,uuid,datetime
from pathlib import Path
import importlib.util
import pandas as pd
from . import tecanio

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
    
    pkl_dict=configgy['pklplates']
    pkl_fldrpath=determinepath(cfdir,pkl_dict)

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
        c.execute('''DROP TABLE LSCRIPTS''')
        for pklf in os.listdir(pkl_fldrpath):
            os.remove(os.path.join(pkl_fldrpath,pklf))
    c.execute('''CREATE TABLE IF NOT EXISTS PLATES (plateid text, filename text, sheetname text,scriptname text, pklplate glob)''')
    c.execute('''CREATE TABLE IF NOT EXISTS LSCRIPTS (scriptname text, pkl_ctime text, pklfname text)''')

    outplates=[]
    for ls in load_scripts:
        load_required=False
        lsentry=str(ls[0].name)+':'+ls[1]
        try:
            c.execute('''SELECT * FROM LSCRIPTS WHERE scriptname=(?)''',(lsentry,))
            prev_entry=c.fetchone()
            prev_load_time=prev_entry['pkl_ctime']
            prev_pklfname=prev_entry['pklfname']
            prev_picklefpath=pkl_fldrpath / prev_pklfname
            #compare create time with load_script's mod time
            ls_mtime=os.path.getmtime(ls[0])
            if ls_mtime>float(prev_load_time):
                load_required=True
                c.execute('''DELETE FROM LSCRIPTS WHERE scriptname=(?)''',(lsentry,))
                c.execute('''DELETE FROM PLATES WHERE scriptname=(?)''',(lsentry,))
                os.remove(prev_picklefpath)
                #os.remove(prev_platepath???)
        except:
            print(f'reloading plates using {lsentry}')
            load_required=True
        if load_required:
            spec=importlib.util.spec_from_file_location(str(ls[0]),ls[0])#'kdfs/dload_dir/load_scripts/testscript.py','kdfs/dload_dir/load_scripts/testscript.py')#
            themodule=importlib.util.module_from_spec(spec)
            spec.loader.exec_module(themodule)
            m2c=getattr(themodule,ls[1])
            plates=m2c()
            for plate in plates:
                newtuple=(plate.plateid,plate.ifpath.name,plate.expsheet,lsentry,pickle.dumps(plate))
                c.execute('''INSERT INTO PLATES VALUES (?,?,?,?,?)''',newtuple)
            pklfname=f'pltscached-{uuid.uuid4().hex}'
            pklpath=pkl_fldrpath / pklfname
            with open(pklpath,'wb') as f:
                pickle.dump(plates,f)
            curlstime=os.path.getmtime(pklpath)
            newlstuple=(lsentry,curlstime,pklfname)
            c.execute('''INSERT INTO LSCRIPTS VALUES (?,?,?)''',newtuple)
        conn.commit()
        #now read in the plates--
        c.execute('''SELECT * FROM LSCRIPTS WHERE scriptname=(?)''',(lsentry,))
        lsrow=c.fetchone()
        picklefpath=pkl_fldrpath / lsrow['pklfname']
        outplates.extend(pickle.load(open(picklefpath,'rb')))
    conn.close()
    return outplates


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

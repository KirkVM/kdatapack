import sqlite3,atexit,pickle
import yaml,os,uuid,datetime
from pathlib import Path
import importlib.util
import pandas as pd
from . import tecanio
from .tecanio import TecanSet

#calls to 'determinepath' require current running path...
#for now, using hack/shortcut that current path is same location as the config files
def determinepath(curpath,cfdict):
    if cfdict['filepath_prefix_type']=='relative':
        thepath=curpath
    elif cfdict['filepath_prefix_type']=='envar':
        thepath=Path(os.environ[cfdict['filepath_prefix']])
    elif cfdict['filepath_prefix_type']=='absolute':
        thepath=Path(cfdict['filepath_prefix'])
    thepath= thepath / Path(cfdict['filepath_ending'])
    return thepath

def loadplates_from_confiles(cfpathstr,settingsfpathstr,username,refresh_all):
    '''Loads plates using a yaml settings file

    Arguments:
        cfpathstr: path to load_scripts yaml settings file
        settingsfpathstr: path to yaml file containing information about datadirectory for a user
        username: name to refer to in settingsfpathstr
        refresh_all: whether to delete all cached (pickled) plates & tables in DB and start over
    '''
    cfp=Path(cfpathstr)
    cfdir=cfp.parent
    with open(cfp,'r') as f:
        configgy=yaml.safe_load(f)
    #logdb first:
    logdb_dict=configgy['logdb']
    logdbpath=determinepath(cfdir,logdb_dict)
    
    pkl_dict=configgy['pklplates']
    pkl_fldrpath=determinepath(cfdir,pkl_dict)

    #now get user/device settings from settingsfpath
    sfp=Path(settingsfpathstr)
    sfdir=sfp.parent
    with open(sfp,'r') as f:
        usersdict=yaml.safe_load(f)
    datapathdict={}
    for name in usersdict:
        if name==username:
            datapathdict[name]=determinepath(sfdir,usersdict[name])
    assert (username.lower() in [x.lower() for x in datapathdict.keys()] ),\
        f"invalid user name {username}. Accepts one of {[x.lower() for x in datapathdict.keys()]}"

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
        c.execute('''DROP TABLE IF EXISTS PLATES''')
        c.execute('''DROP TABLE IF EXISTS LSCRIPTS''')
        if pkl_fldrpath.exists():
            for pklf in os.listdir(pkl_fldrpath):
                os.remove(os.path.join(pkl_fldrpath,pklf))
        print('**reloading all files**')
    c.execute('''CREATE TABLE IF NOT EXISTS PLATES (plateid text, filename text, sheetname text,scriptname text, pklplate glob)''')
    c.execute('''CREATE TABLE IF NOT EXISTS LSCRIPTS (scriptname text, pkl_ctime text, pklfname text)''')

    ##TODO- move all DB commands to end so that they only execute if all functions successful
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
                print(f'reloading plates using {lsentry} due to loadscript update')
                #os.remove(prev_platepath???)
        except:
            print(f'loading plates from file using {lsentry} (no existing cache)')
            load_required=True
#            conn.close()
        if load_required:
            modulepath=Path(ls[0])
            spec=importlib.util.spec_from_file_location(modulepath.name,modulepath)#'kdfs/dload_dir/load_scripts/testscript.py','kdfs/dload_dir/load_scripts/testscript.py')#
            themodule=importlib.util.module_from_spec(spec)
            spec.loader.exec_module(themodule)
            m2c=getattr(themodule,ls[1])
            try:
                plates=m2c(datapathdict[username])
            except:
                conn.close()
                exit(f'Error occurred in load script {lsentry}. Quitting.')
            for plate in plates:
                newpltuple=(plate.plateid,plate.ifpath.name,plate.expsheet,lsentry,pickle.dumps(plate))
                c.execute('''INSERT INTO PLATES VALUES (?,?,?,?,?)''',newpltuple)
            pklfname=f'pltscached-{uuid.uuid4().hex}'
            pklpath=pkl_fldrpath / pklfname
            if pkl_fldrpath.exists()==False:
                os.mkdir(pkl_fldrpath)
            with open(pklpath,'wb') as f:
                pickle.dump(plates,f)
            curlstime=os.path.getmtime(pklpath)
            newlstuple=(lsentry,curlstime,pklfname)
            c.execute('''INSERT INTO LSCRIPTS VALUES (?,?,?)''',newlstuple)
        conn.commit()
        #now read in the plates--
        c.execute('''SELECT * FROM LSCRIPTS WHERE scriptname=(?)''',(lsentry,))
        lsrow=c.fetchone()
        picklefpath=pkl_fldrpath / lsrow['pklfname']
        outplates.extend(pickle.load(open(picklefpath,'rb')))
    conn.close()
    return outplates


def load_tecandata(cfpathstr,username,refresh_all=False):
    '''reads excel files using load scripts/functions as specified in a config file

    Arguments:
        cfpathstr: path to yaml config file
    Keyword Arguments:
        refresh_all: whether to reset all files for a fresh re-read (default False)
    '''
    plates=loadplates_from_confiles(cfpathstr,'user_settings.yml',username,refresh_all)
    tpset=TecanSet(plates)
    return tpset
#    mergedf=pd.concat([x for x in dfs],ignore_index=True)
#    return mergedf

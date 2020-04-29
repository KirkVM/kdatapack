import requests,bs4

import re
from datetime import date
def get_headers(row):
    categories=[]
    cat_indices={}
    for x,header in enumerate(row.find_all('th')):
        catname=header.text
        categories.append(catname)
        if catname=='Entry':
            cat_indices[x]='uparcid'
        if catname=='Organisms':
            cat_indices[x]='organisms'
        if catname=='UniProtKB':
            cat_indices[x]='uprotids'
        if catname=='First seen':
            cat_indices[x]='first_seen'
        if catname=='Last seen':
            cat_indices[x]='last_seen'
        if catname=='Length':
            cat_indices[x]='length'
    assert(len(set(['uparcid','organisms','uprotids','first_seen','last_seen','length']).difference\
           (cat_indices.values()))==0),f"missing a category. Id are: {cat_indices.values()}"
    return cat_indices


def get_data(row,cat_indices):
    datadict={}
    for didx,data in enumerate(row.find_all('td')):
        if didx not in cat_indices.keys():continue
        datastr=data.get_text(separator=';',strip=True)
        if cat_indices[didx] in ['first_seen','last_seen']:
            thedates=[]
            for val in datastr.split(';'):
                thedates.append(date(*[int(x) for x in val.split('-')]))
            datadict[cat_indices[didx]]=thedates          
        else:
            datadict[cat_indices[didx]]=datastr.split(';')
    return datadict

def get_upids(strid):
    requestURL = f"https://uniprot.org/uniparc/?query={strid}"
    r = requests.get(requestURL)
    soup=bs4.BeautifulSoup(r.content,'lxml')
    soupbody=soup.body
   
    crt=None
    cr=None
    divs=soupbody.find_all('div')
    for div in divs:
        #print(div.attrs)
        if 'class' in div.attrs:
            #print(div.attrs['class'])
            if 'results' in div.attrs['class']:
                cr=div
    for child in cr.children:
        if child.name=='table':
            crt=child
        
    rows=crt.find_all('tr')
    for row in rows:
        if len(list(row.find_all('th')))>1:
            headerdict=get_headers(row)
        if len(list(row.find_all('td')))>1:
            datadict=get_data(row,headerdict)
    return datadict
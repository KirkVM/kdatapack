import bs4,requests,urllib
from Bio import SeqIO
class Uniprot_ProteinEntry:
    def __init__(self):
        self.record=None #the unmodified Biopython-parsed xml record
        self.id=None #
        self.seq=None #a Bio.Seq Seq instance
        self.name=None
        self.description=None
        self.cazy_motifs=[]
        self.pfam_motifs=[]
        self.pdbs=[]
        self.ecs=[]
        self.dois=[]
        self.activities=[]
        self.otherstuff=None
        self.journal_references=None
        self.other_references=None
        self.kinetic_params=None
        self.exp_evidences=None

#class Activity:
#    def __init__(self,evidx,)

class JrnlReference:
    def __init__(self,refidx,reftype,refdate):#,pmid=None,doi=None):
        self.refidx=refidx
        self.reftype=reftype
        self.refdate=refdate
        self.pmid=None
        self.doi=None

class KineticParam:
    def __init__(self,text,evidxs,ktype_label,km_mined,kcat_mined,vmax_mined):#infodict):
        self.text=text
        self.evidence_index=evidxs #list of integers
        self.ktype_label=ktype_label
        self.km_mined=km_mined
        self.vmax_mined=vmax_mined
        self.kcat_mined=kcat_mined
        self.journal_reference=None #add this later if possible
        self.doi=[]
    def __str__(self):
        retstr=f'ktype_field: {self.ktype_label}\n'
        retstr+=f'km= {self.km_mined or self.ktype_label.lower()=="km"},'
        retstr+=f'kcat= {self.kcat_mined or self.ktype_label.lower()=="kcat"},'
        retstr+=f'vmax= {self.vmax_mined or self.ktype_label.lower()=="vmax"}\n'
        retstr+=f'description: {self.text}\n'
        if self.doi is not None:
            retstr+=f'article DOIs: {self.doi}\n'
        return retstr
#        return str(self.evidence_number)
#        print(self.evidence_number)

class EvidenceSource:
    def __init__(self,evidx,pmid=None,doi=None):
        self.evidence_index=evidx
        self.pmid=pmid
        self.doi=doi

def parse_protein_record(record):
    upe=Uniprot_ProteinEntry()
    upe.record=record #this will keep unmodified record
    upe.id=record.id
    upe.seq=record.seq
    upe.description=record.description
    upe.name=record.name
    for dbxref in record.dbxrefs:
        dbname,dbval=dbxref.split(':',maxsplit=1)
        if dbname.lower()=='cazy':
            upe.cazy_motifs.append(dbval)
        if dbname.lower()=='ec':
            upe.ecs.append(dbval)
    return upe

#description='Mannan endo-1,4-beta-mannosidase', dbxrefs=['BRENDA:3.2.1.78', 'BioCyc:RMAR518766:G1GGK-18-MONOMER', 'CAZy:CBM35', 'CAZy:GH26', 'DOI:10.1007/s002530000351', 'DOI:10.4056/sigs.46736', 'EC:3.2.1.78', 'EMBL:CP001807', 'EMBL:X90947', 'EnsemblBacteria:ACY46925', 'GO:GO:0006080', 'GO:GO:0016985', 'GO:GO:0030246', 'Gene3D:2.60.120.260', 'HOGENOM:CLU_016930_1_1_10', 'InterPro:IPR000805', 'InterPro:IPR005084', 'InterPro:IPR008979', 'InterPro:IPR017853', 'InterPro:IPR022790', 'InterPro:IPR026444', 'KEGG:rmr:Rmar_0016', 'KO:K01218', 'NCBI Taxonomy:518766', 'OMA:WGWYLID', 'OrthoDB:669354at2', 'PANTHER:PTHR40079', 'PIR:T10748', 'PRINTS:PR00739', 'PROSITE:PS51175', 'PROSITE:PS51764', 'Pfam:PF02156', 'Proteomes:UP000002221', 'PubMed:10919332', 'PubMed:21304669', 'RefSeq:WP_012842537.1', 'SMR:P49425', 'STRING:518766.Rmar_0016', 'SUPFAM:SSF49785', 'SUPFAM:SSF51445', 'Swiss-Prot:D0MK12', 'Swiss-Prot:MANA_RHOM4', 'Swiss-Prot:P49425', 'TIGRFAMs:TIGR04183', 'eggNOG:COG4124', 'eggNOG:ENOG4105D3P'])


def parse_protein_xml(accession):
    requestURL = f"https://uniprot.org/uniprot/{accession}.xml"
    r=requests.get(requestURL)
    #soup=bs4.BeautifulSoup(xmlpage.read(),'lxml')
    soup=bs4.BeautifulSoup(r.content,'lxml')
    refxmls=soup.find_all('reference')
    litrefs=[]
    orefs=[]
    for ref in refxmls:
        rtype=ref.citation['type']
        if rtype!='journal article': #'submission' seems only other alternative
            orefs.append(ref)
            continue
        rkey=ref['key']
        rdate=ref.citation['date'] #this just gives year though??
        litref=JrnlReference(rkey,rtype,rdate)
        dbrefs=ref.citation.find_all('dbreference')
        for dbref in dbrefs:
            dbref_type=dbref['type']
            assert(dbref_type.lower() in ['pubmed','doi']),f'what type of ref is this? {ref.citation}'
            if dbref_type.lower()=='pubmed':
                litref.pmid=dbref['id']
            elif dbref_type.lower()=='doi':
                litref.doi=dbref['id']
        litrefs.append(litref)
    
    annotations=soup.find_all('comment')
    functions=[]
    activities=[]
    kps=[]
    for annotation in annotations:
        if annotation['type']=='function':
            for t in annotation.find_all('text'):
                funky={'evidence':t['evidence']}
                funky['text']=t.text
                functions.append(funky)
        if annotation['type']=='catalytic activity':
            for r in annotation.find_all('reaction'):
                reakt={'evidence':r['evidence']}
                reakt['text']=r.text
                reakt['dbref']=[]
                for dbref in r.find_all('dbreference'):
                    reakt['dbref'].append(dbref.attrs)
                activities.append(reakt)
            j=annotation
        if annotation['type']=='biophysicochemical properties':
            kins=annotation.find_all('kinetics') 
            assert(len(kins)==1),'more than 1 kinetics field?'
            kinkids=[x for x in kins[0].children if x.name!=None]
            for kinkid in kinkids:
                kinevidxs=[int(x) for x in kinkid['evidence'].split()] #can be more than one!
                if kinkid.name!='text':
                    ktype_lbl=kinkid.name
                else:
                    ktype_lbl='unknown'
                kintext=kinkid.text
                #now let's mine the text
                kcat_mined='kcat' in kinkid.text.lower().split()
                km_mined='km' in kinkid.text.lower().split()
                vmax_mined='vmax' in kinkid.text.lower().split()
                kp=KineticParam(kintext,kinevidxs,ktype_lbl,km_mined,kcat_mined,vmax_mined)#infodict):P10477
                kps.append(kp)
    
    evidences=soup.find_all('evidence')
    exp_evidences=[]
    for evidence in evidences:
        expev={'idx':None,'pmid':None,'doi':None}
        if evidence['type']=='ECO:0000269':
            expev['idx']=int(evidence['key'])
            sources=evidence.find_all('source')
            assert(len(sources)==1),"more than 1 source"
            #print(sources[0])#.dbreference.attr)
            dbrefs=sources[0].find_all('dbreference')
            assert(len(dbrefs)==1),"more than 1 dbref"
            dbattrdict=dbrefs[0].attrs
            if dbattrdict['type'].lower()=='pubmed':
                expev['pmid']=dbattrdict['id']
            elif dbattrdict['type'].lower()=='doi':
                expev['doi']=dbattrdict['id']
            ev=EvidenceSource(expev['idx'],pmid=expev['pmid'],doi=expev['doi'])
            exp_evidences.append(ev)
#    evs=[]
#    for expev in exp_evidences:
#        ev=EvidenceSource(expev['idx'],pmid=expev['pmid'],doi=expev['doi'])
#        evs.append(ev)
    #now use SeqIO utility for most stuff---
    #I don't understand how to use requests to get biopython parseable---
    xmlpage=urllib.request.urlopen(requestURL)    
    bio_record= SeqIO.read(xmlpage, "uniprot-xml")
    upe=parse_protein_record(bio_record)
    upe.journal_refs=litrefs
    upe.other_references=orefs
    upe.kinetic_params=kps
    upe.exp_evidences=exp_evidences
    for exp_ev in upe.exp_evidences:
        for jref in upe.journal_refs:
            if exp_ev.pmid is not None:
                if exp_ev.pmid==jref.pmid:
                    exp_ev.doi=jref.doi
    for kp in upe.kinetic_params:
        for exp_ev in upe.exp_evidences:
            if exp_ev.evidence_index in kp.evidence_index:
                kp.doi.append(exp_ev.doi)
#                    print('match')
#            print(exp_ev.evidence_index)
    return upe
#b=example_ref.citation.type
#for dbref in example_ref.find_all('dbreference'):
#    print(dbref['type'],dbref['id'])



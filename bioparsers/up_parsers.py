import bs4,requests,urllib
from Bio import SeqIO
class Uniprot_ProteinEntry:
    def __init__(self,requestURL):
        self.requestURL=requestURL
        self.biopython_uprecord=None #the unmodified Biopython-parsed xml record
        self.uniprot_id=None #
        self.seq=None #a Bio.Seq Seq instance
        self.name=None
        self.description=None
        self.cazy_motifs=[]
        self.ecs=[]
        #not implemented
        #self.pfamsomething
        #above is currently filled in by biopython
        #below is filled in by my parsing
        self.functions=None
        self.activities=None
        self.journal_references=None
        self.other_references=None
        self.kinetic_params=None
        self.exp_evidences=None
        self.pfam_matches=None
        self.dbseqs=None
        self.soup=None #as utility attaching extracted page to this field

        #currently not implemented
        self.pdbs=[] #this is easily parse-able from bipython_uprecord
        self.dois=[] #this is available within journal_references
        self.alternative_ids=[] #are there other entries redundant w/ this
        self.otherstuff=None

    def biopython_parse(self):
        xmlpage=urllib.request.urlopen(self.requestURL)    
        self.biopython_uprecord= SeqIO.read(xmlpage, "uniprot-xml")
        self.uniprot_id=self.biopython_uprecord.id
        self.seq=self.biopython_uprecord.seq
        self.description=self.biopython_uprecord.description
        self.name=self.biopython_uprecord.name
        for dbxref in self.biopython_uprecord.dbxrefs:
            dbname,dbval=dbxref.split(':',maxsplit=1)
            if dbname.lower()=='cazy':
                self.cazy_motifs.append(dbval)
            if dbname.lower()=='ec':
                self.ecs.append(dbval)
    def acc_info_string(self):
        retstr=self.uniprot_id
        for d in self.dbseqs:
            if d.type.lower() in ['embl','genbank','ddbj'] and d.accession is not None:
                retstr+=f', {d.accession}'
        return retstr

#class Activity:
#    def __init__(self,evidx,)

class JrnlReference:
    def __init__(self,refxml):
        self.refidx=refxml['key']
        self.reftype=refxml.citation['type']
        self.refdate=refxml.citation['date']
        self.reftitle=refxml.citation.title.text
        self.pmid=None
        self.doi=None
        dbrefs=refxml.citation.find_all('dbreference')
        for dbref in dbrefs:
            dbref_type=dbref['type']
            assert(dbref_type.lower() in ['pubmed','doi']),f'what type of ref is this? {ref.citation}'
            if dbref_type.lower()=='pubmed':
                self.pmid=dbref['id']
            elif dbref_type.lower()=='doi':
                self.doi=dbref['id']
        self.function_scope=False
        self.activity_scope=False
        self.biophysicochemical_scope=False
        self.specificity_scope=False
        self.xtal_scope=False
        self.mechanism_scope=False
        scopes=refxml.find_all('scope')
        for scope in scopes:
            if scope.text.lower().strip()=='function':
                self.function_scope=True
            elif scope.text.lower().strip()=='catalytic activity':
                self.activity_scope=True
            elif 'biophysicochemical' in scope.text.lower():
                self.biophysicochemical_scope=True
            elif scope.text.lower().strip()=='substrate specificity':
                self.specificity_scope=True
            elif 'x-ray' in scope.text.lower():
                self.xtal_scope=True
            elif scope.text.lower().strip()=='reaction mechanism':
                self.mechanism_scope=True


class KineticParam:
    def __init__(self,kinxml):#text,evidxs,ktype_label,km_mined,kcat_mined,vmax_mined):#infodict):
        self.evidence_index=[int(x) for x in kinxml['evidence'].split()] #can be more than one!
        self.ktype_label='unknown'
        if kinxml.name!='text':
            self.ktype_label=kinxml.name
        self.text=kinxml.text
        #now let's mine the text
        self.kcat_mined='kcat' in kinxml.text.lower().split()
        self.km_mined='km' in kinxml.text.lower().split()
        self.vmax_mined='vmax' in kinxml.text.lower().split()
        self.journal_reference=None
        self.dois=[]

    def __str__(self):
        retstr=f'ktype_field: {self.ktype_label}\n'
        retstr+=f'km= {self.km_mined or self.ktype_label.lower()=="km"},'
        retstr+=f'kcat= {self.kcat_mined or self.ktype_label.lower()=="kcat"},'
        retstr+=f'vmax= {self.vmax_mined or self.ktype_label.lower()=="vmax"}\n'
        retstr+=f'description: {self.text}\n'
        if self.dois is not None:
            retstr+=f'article DOIs: {self.dois}\n'
        return retstr
#        return str(self.evidence_number)
#        print(self.evidence_number)

class EnzFunction:
    def __init__(self,funcxml):
        self.evidence_index=[]
        allchildren=list(set([x.name for x in funcxml.children]))
        allchildren.remove(None)
        assert(len(allchildren)==1), "should only be one field (and next step enforces that it is text)"
        functexts=funcxml.find_all('text')
        assert (len(functexts)==1),"more than 1 reaction here"
        func_textxml=functexts[0]
        if 'evidence' in list(func_textxml.attrs.keys()):# not None:
            self.evidence_index=[int(x) for x in func_textxml['evidence'].split()] #can be more than one!
        self.text=func_textxml.text

class EnzActivity:
    def __init__(self,funcxml):
        self.evidence_index=[]
        reacts=funcxml.find_all('reaction')
        assert (len(reacts)==1),"more than 1 reaction here"
        reactxml=reacts[0]
        if 'evidence' in list(reactxml.attrs.keys()):# not None:
            self.evidence_index=[int(x) for x in reactxml['evidence'].split()] #can be more than one!
#        self.evidence_index=[int(x) for x in reactxml['evidence'].split()] #can be more than one!
        self.evidence_sources=[]
        self.text=reactxml.text
        self.ec=None
        self.dois=[]
        #self.?=None #are there other relevant dbreference types?
        for dbref in reactxml.find_all ('dbreference'):
            if dbref['type'].lower()=='ec':
                assert (self.ec is None),'2 ecs in a single reaction???'
                self.ec=dbref['id']
    def __str__(self):
        retstr=f'Activity ec: {self.ec}\n'
        retstr+=f'Details: {self.text}'
        if self.dois is not None:
            retstr+=f'article DOIs: {self.dois}\n'
        return retstr
#
class EvidenceSource:
    def __init__(self,evxml):
        self.evidence_index=int(evxml.attrs['key'])
        self.doi=None
        self.pmid=None
        self.osrc_dict={}
        sources=evxml.find_all('source')
        assert (len(sources)==1),"more than 1 source here"
        srcxml=sources[0]

        dbrefs=srcxml.find_all('dbreference')
        assert (len(dbrefs)==1),f"{self.evidence_index},{len(dbrefs)},more than 1 source here"
        dbxml=dbrefs[0]
        dbattrdict=dbxml.attrs
        if dbattrdict['type'].lower()=='pubmed':
            self.pmid=dbattrdict['id']
        elif dbattrdict['type'].lower()=='doi':
            self.doi=dbattrdict['id']
        else:
            self.osrc_dict[dbattrdict['type']]=dbattrdict['id']

#description='Mannan endo-1,4-beta-mannosidase', dbxrefs=['BRENDA:3.2.1.78', 'BioCyc:RMAR518766:G1GGK-18-MONOMER', 'CAZy:CBM35', 'CAZy:GH26', 'DOI:10.1007/s002530000351', 'DOI:10.4056/sigs.46736', 'EC:3.2.1.78', 'EMBL:CP001807', 'EMBL:X90947', 'EnsemblBacteria:ACY46925', 'GO:GO:0006080', 'GO:GO:0016985', 'GO:GO:0030246', 'Gene3D:2.60.120.260', 'HOGENOM:CLU_016930_1_1_10', 'InterPro:IPR000805', 'InterPro:IPR005084', 'InterPro:IPR008979', 'InterPro:IPR017853', 'InterPro:IPR022790', 'InterPro:IPR026444', 'KEGG:rmr:Rmar_0016', 'KO:K01218', 'NCBI Taxonomy:518766', 'OMA:WGWYLID', 'OrthoDB:669354at2', 'PANTHER:PTHR40079', 'PIR:T10748', 'PRINTS:PR00739', 'PROSITE:PS51175', 'PROSITE:PS51764', 'Pfam:PF02156', 'Proteomes:UP000002221', 'PubMed:10919332', 'PubMed:21304669', 'RefSeq:WP_012842537.1', 'SMR:P49425', 'STRING:518766.Rmar_0016', 'SUPFAM:SSF49785', 'SUPFAM:SSF51445', 'Swiss-Prot:D0MK12', 'Swiss-Prot:MANA_RHOM4', 'Swiss-Prot:P49425', 'TIGRFAMs:TIGR04183', 'eggNOG:COG4124', 'eggNOG:ENOG4105D3P'])
class PfamMatch:
    def __init__(self,pxml):
        self.accession=pxml.attrs['id']
        self.name=None
        self.num_matches=1
        for propxml in pxml.find_all('property'):
            if propxml.attrs['type'].lower()=='entry name':
                self.name=propxml.attrs['value']
            if propxml.attrs['type'].lower()=='match status':
                self.num_matches=int(propxml.attrs['value'])
    def __str__(self):
        retstr=self.accession
        retstr+=f' ({self.name})'
        if self.num_matches>1:
            retstr+=f' (x{self.num_matches})'
        return retstr

class DBSequence:
    def __init__(self,sxml):
        self.type=sxml.attrs['type'] #EMBL, Genbank,...
        self.id=sxml.attrs['id']
        self.accession=None
        self.otherinfo_dicts=[]
        for propxml in sxml.find_all('property'):
            if propxml.attrs['type'].lower()=='protein sequence id':
                self.accession=propxml.attrs['value']
            else:
                self.otherinfo_dicts.append(propxml.attrs) 


def parse_uniprot_xml(accession):
    requestURL = f"https://uniprot.org/uniprot/{accession}.xml"
    #use SeqIO utility for most stuff---
    #I don't understand how to use requests to get biopython parseable---
    upe=Uniprot_ProteinEntry(requestURL)
    upe.biopython_parse()

    #currently calling requests twice (once in UProtEntry Biopython call, once here)
    #soup=bs4.BeautifulSoup(xmlpage.read(),'lxml')
    r=requests.get(requestURL)
    soup=bs4.BeautifulSoup(r.content,'lxml')
    upe.soup=soup #keep starting soup object attached to the object for convenience

    #########Parse references##############
    refxmls=soup.find_all('reference')
    litrefs=[]
    orefs=[]
    for refxml in refxmls:
        rtype=refxml.citation['type']
        if rtype!='journal article': #'submission' seems only other alternative
            orefs.append(refxml)
            continue
        litrefs.append(JrnlReference(refxml))
    
    #########Parse annotations################
    annotations=soup.find_all('comment')
    kps=[]
    functions=[]
    activities=[]
    for annotation in annotations:
        ##########biophys....
        if annotation['type']=='biophysicochemical properties':
            #currently only pulling kinetics. Should be just one (potentially multiple KM/Kcat entries w/in)
            kins=annotation.find_all('kinetics') 
            if len(kins)==0: #property is something diff
                continue
            assert(len(kins)==1),'more than 1 kinetics field?'
            kinxmls=[x for x in kins[0].children if x.name!=None]
            for kinxml in kinxmls:
                kps.append(KineticParam(kinxml))
        ######function########
        if annotation['type']=='function':
            functions.append(EnzFunction(annotation))
        if annotation['type']=='catalytic activity':
           activities.append(EnzActivity(annotation))
      #add here what about catalytic activity or substrate specificty....
    #check parsing vs biopython
    assert (len(set([x.ec for x in activities]).difference(upe.ecs))==0),'ec conflict?'
    #######Parse evidence sources#########
    allsrc_evidences=[]
    evxmls=soup.find_all('evidence')
    ev_srcxmls= [x  for x in evxmls if 'dbreference' in list(y.name for y in x.descendants)]
    for evxml in ev_srcxmls:
        #print(evxml)
        allsrc_evidences.append(EvidenceSource(evxml))
    ######Parse domain annotations###########
    pfams=[]
    pfamxmls=[x for x in soup.find_all('dbreference') if \
                (x.attrs['type'].lower()=='pfam' and 'property' in [y.name for y in x.children])]

    for pxml in pfamxmls:
        pfams.append(PfamMatch(pxml))
    #####Parse alternative accession codes for sequence######
    dbseqs=[]
    dbseqxmls=[x for x in soup.find_all('dbreference') if x.attrs['type'].lower() in ['embl','genbank','ddbj']]
    for dbseqxml in dbseqxmls:
        dbseqs.append(DBSequence(dbseqxml))

    #add my parsed data to the upe
    upe.journal_references=litrefs
    upe.other_references=orefs
    upe.functions=functions
    upe.activities=activities
    upe.kinetic_params=kps
    upe.exp_evidences=allsrc_evidences
    upe.pfam_matches=pfams
    upe.dbseqs=dbseqs
    #now connect journals to dois
    for exp_ev in upe.exp_evidences:
        for jref in upe.journal_references:
            if exp_ev.pmid is not None:
                if exp_ev.pmid==jref.pmid:
                    exp_ev.doi=jref.doi
    #now connect kinetic data to journals+doi
    for kp in upe.kinetic_params:
        for exp_ev in upe.exp_evidences:
            if exp_ev.evidence_index in kp.evidence_index:
                if exp_ev.doi is not None:
                    kp.dois.append(exp_ev.doi)
    for ca in upe.activities:
        for exp_ev in upe.exp_evidences:
            if exp_ev.evidence_index in ca.evidence_index:
                if exp_ev.doi is not None:
                    ca.dois.append(exp_ev.doi) #this is a convenience shortcut. full evsrc linked in next line
                ca.evidence_sources.append(exp_ev)
    return upe
#b=example_ref.citation.type
#for dbref in example_ref.find_all('dbreference'):
#    print(dbref['type'],dbref['id'])



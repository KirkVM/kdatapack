from bokeh.plotting import curdoc,ColumnDataSource,figure #show,figure,ColumnDataSource,curdoc
from bokeh.layouts import row,column
from bokeh.models import MultiLine,Quad,HoverTool,ResetTool,BoxZoomTool,PanTool,Range1d,Text,TapTool,Label#,Line
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import DataTable,TableColumn,MultiSelect,AutocompleteInput
from bokeh.models.widgets.buttons import Button
#from bokeh.models import HoverTool,ResetTool,BoxZoomTool,PanTool,Line,Range1d
#from bokeh.layouts import row
#from bokeh.models.glyphs import Text
import sys,pickle,ete3,collections,os
from ete3 import Tree
from phylotreepack import phylotree

stuffcode="""
    var data=cds.data;
    var dtdata=dtcds.data;
    const inds = cds.selected.indices;

    var addedaccs=[];
    var rmaccs=[];
    var curaccs=dtdata['accs'].slice(0);
    for (var i = 0; i < inds.length; i++) {
        if (data['frame_lws'][inds[i]]==0.5) {
            data['frame_lws'][inds[i]]=3;
            data['rpnl_tfsize'][inds[i]]="10pt";
            data['qhatch'][inds[i]]="diagonal_cross";
            addedaccs.push(data['gbacc'][inds[i]]);
        }else{
            data['frame_lws'][inds[i]]=0.5;
            data['rpnl_tfsize'][inds[i]]="4pt";
            data['qhatch'][inds[i]]="blank";
            rmaccs.push(data['gbacc'][inds[i]]);
        }
    } 
    for (var i=0; i<addedaccs.length;i++){
        curaccs.push(addedaccs[i]);
    }
    for (var i=0; i<rmaccs.length;i++){
        for (var j=0;j<curaccs.length;j++){
            if (curaccs[j]==rmaccs[i]) {
                curaccs.splice(j,1);
            }
        }
    }
    dtcds.data={'accs':curaccs};
    dtcds.change.emit();
    cds.change.emit();
"""

def tablecallback(attr,old,new):
    for acc in dtcds.data['accs']:
        tax_sking=' '
        tax_phylum=' '
        tax_class=' '
        tax_genus=' '
        tax_sp=' '
        md_ecs=' '
        md_pdbs=' '
        for x,cdacc in enumerate(ptree.leaf_cds.data['gbacc']):
            if cdacc==acc:
                tax_sking=ptree.leaf_cds.data['superkingdom'][x]
                tax_phylum=ptree.leaf_cds.data['phylum'][x]
                tax_class=ptree.leaf_cds.data['class'][x]
                tax_genus=ptree.leaf_cds.data['genus'][x]
                tax_sp=ptree.leaf_cds.data['species'][x]

        for x,ecacc in enumerate(ec_cds.data['gbacc']):
            if ecacc==acc:
                md_ecs=ec_cds.data['ecs'][x]
        for x,mdacc in enumerate(pdb_cds.data['gbacc']):
            if mdacc==acc:
                md_pdbs=pdb_cds.data['pdbids'][x]
        print(f'{acc}\t{tax_sking}\t{tax_phylum}\t{tax_class}\t{tax_genus}\t{tax_sp}\t{md_ecs}\t{md_pdbs}')



def pmsfunc():
    qclist=ptree.leaf_cds.data['qcolor'][:]
    qfalist=ptree.leaf_cds.data['qfillalpha'][:]
#    print(type(list(gms.value)))
    for x,sp in enumerate(ptree.leaf_cds.data['phylum']):
        if sp in list(pms.value):
            qclist[x]='Blue'
            qfalist[x]=0.5
        else:
            qclist[x]=None
            qfalist[x]=0
    ptree.leaf_cds.data['qcolor']=qclist#['Blue' for _ in range(len(leaf_source.data['species']))]
    ptree.leaf_cds.data['qfillalpha']=qfalist

def cmsfunc():
    qclist=ptree.leaf_cds.data['qcolor'][:]
    qfalist=ptree.leaf_cds.data['qfillalpha'][:]
    for x,sp in enumerate(ptree.leaf_cds.data['class']):
        if sp in list(cms.value):
            qclist[x]='Blue'
            qfalist[x]=0.5
        else:
            qclist[x]=None
            qfalist[x]=0
    ptree.leaf_cds.data['qcolor']=qclist#['Blue' for _ in range(len(leaf_source.data['species']))]
    ptree.leaf_cds.data['qfillalpha']=qfalist

def gmsfunc():
    qclist=ptree.leaf_cds.data['qcolor'][:]
    qfalist=ptree.leaf_cds.data['qfillalpha'][:]
    for x,sp in enumerate(ptree.leaf_cds.data['genus']):
        if sp in list(gms.value):
            qclist[x]='Blue'
            qfalist[x]=0.5
        else:
            qclist[x]=None
            qfalist[x]=0
    ptree.leaf_cds.data['qcolor']=qclist#['Blue' for _ in range(len(leaf_source.data['species']))]
    ptree.leaf_cds.data['qfillalpha']=qfalist

def sfmsfunc():
    qclist=ptree.leaf_cds.data['qcolor'][:]
    qfalist=ptree.leaf_cds.data['qfillalpha'][:]
    for x,sp in enumerate(ptree.leaf_cds.data['subfamstr']):
        if sp in list(sfms.value):
            qclist[x]='Green'
            qfalist[x]=0.5
        else:
            qclist[x]=None
            qfalist[x]=0
    ptree.leaf_cds.data['qcolor']=qclist#['Blue' for _ in range(len(leaf_source.data['species']))]
    ptree.leaf_cds.data['qfillalpha']=qfalist

def kmsfunc():
    qclist=ptree.leaf_cds.data['qcolor'][:]
    qfalist=ptree.leaf_cds.data['qfillalpha'][:]
    for x,sp in enumerate(ptree.leaf_cds.data['superkingdom']):
        if sp in list(kms.value):
            qclist[x]='Blue'
            qfalist[x]=0.5
        else:
            qclist[x]=None
            qfalist[x]=0
    ptree.leaf_cds.data['qcolor']=qclist#['Blue' for _ in range(len(leaf_source.data['species']))]
    ptree.leaf_cds.data['qfillalpha']=qfalist



def acfunc():
    qclist=ptree.leaf_cds.data['qcolor'][:]
    qfalist=ptree.leaf_cds.data['qfillalpha'][:]
    for x,sp in enumerate(ptree.leaf_cds.data['species']):
        if sp==acw.value:
            qclist[x]='Blue'
            qfalist[x]=0.5
        else:
            qclist[x]=None
            qfalist[x]=0
    ptree.leaf_cds.data['qcolor']=qclist#['Blue' for _ in range(len(leaf_source.data['species']))]
    ptree.leaf_cds.data['qfillalpha']=qfalist#[1.0 for _ in range(len(leaf_source.data['species']))]

def accacfunc():
    qclist=ptree.leaf_cds.data['qcolor'][:]
    qfalist=ptree.leaf_cds.data['qfillalpha'][:]
    for x,sp in enumerate(ptree.leaf_cds.data['gbacc']):
        if sp==accacw.value:
            qclist[x]='Green'
            qfalist[x]=0.5
#        else:
#            qclist[x]=None
#            qfalist[x]=0
    ptree.leaf_cds.data['qcolor']=qclist#['Blue' for _ in range(len(leaf_source.data['species']))]
    ptree.leaf_cds.data['qfillalpha']=qfalist#[1.0 for _ in range(len(leaf_source.data['species']))]


def tcbfunc():
    tc=ptree.calc_treecoverage(dtcds.data['accs'])
    tccds.data={'tc':[tc]}
#    print(tc)
#    print(dtcds.data['accs'])

#    print(ptree.leaf_cds.data['gbacc'])
#    selectedaccs=[ptree.leaf_cds.data['gbacc'][x] for x in \
#                   range(len(ptree.leaf_cds.data['gbacc'])) if ptree.leaf_cds.data['frame_lws'][x]==3]
#    print(dtbl.columns)#data['accs'])
#    print(dtbl.columns[0])#data['accs'])
#    print(dtbl.columns[0]['accs'])#data['accs'])
#    print(ptree.leaf_cds.data['frame_lws'])
#    print(dtbl['accs'])


def psbfunc():
    qclist=ptree.leaf_cds.data['qcolor'][:]
    qfalist=ptree.leaf_cds.data['qfillalpha'][:]
    for x,gbacc in enumerate(ptree.leaf_cds.data['gbacc']):
        if gbacc in preselections:
            qclist[x]='Gray'
            qfalist[x]=0.5
        else:
            qclist[x]=None
            qfalist[x]=0
    ptree.leaf_cds.data['qcolor']=qclist#['Blue' for _ in range(len(leaf_source.data['species']))]
    ptree.leaf_cds.data['qfillalpha']=qfalist#[1.0 for _ in range(len(leaf_source.data['species']))]


print(sys.argv)
ghfam=sys.argv[1]
pkltreename=sys.argv[2]
#with open('st.pkl','rb') as f:
treefpath=os.path.join(os.environ['SCIENCEDIR'],'GHSeqs',ghfam.upper(),pkltreename)
dbfpath=os.path.join(os.environ['SCIENCEDIR'],'GHSeqs',ghfam.upper(),f'{ghfam.upper()}DB.sql')
ncfpath=os.path.join(os.environ['SCIENCEDIR'],'GHSeqs',ghfam.upper(),f'{ghfam.lower()}ds.nc')
print(treefpath)
with open(treefpath,'rb') as f:
    etetree=pickle.load(f)

hover_tool = HoverTool(names=['leaf_node'],tooltips=[    ("GB acc", "@gbacc"),('sf','@subfam'),\
                        ('phylum','@phylum'),('class','@class'),('species','@species')])
hover_tool2 = HoverTool(names=['metadata'],tooltips=[    ("GB acc", "@gbacc"),("ECs", "@ecs"),('PDBs','@pdbids')])
pheight=1400
p1 = figure(plot_width=850,plot_height=pheight,tools=[hover_tool,ResetTool(),BoxZoomTool(),PanTool()])#plot_width=1100, plot_height=700,
    #pheight=int(len(ptree.get_leaves())*1.5)
p2 = figure(plot_width=60,plot_height=pheight,tools=[hover_tool2],x_range=Range1d(0,2),y_range=p1.y_range)
ptree=phylotree.PhyloTree(etetree)
ptree.update_leafcdsdict_fromdb(dbfpath)
ptree.update_leafcdsdict_fromxr(ncfpath)
ptree.bokehdraw(plot=p1,rotation=0)

fglyph=MultiLine(xs='frame_xs',ys='frame_ys',line_width='frame_lws')
p1.add_glyph(ptree.leaf_cds,fglyph)#,name='leaf_node')
p1.add_glyph(ptree.internal_cds,fglyph)#,name='internal_node')
qglyph=p1.quad(left='nodebox_lefts',right='nodebox_rights',bottom='nodebox_bottoms',top='nodebox_tops',fill_color='qcolor',line_alpha=0,\
                        fill_alpha='qfillalpha',hatch_pattern='qhatch',source=ptree.leaf_cds,name='leaf_node')#,fill_alpha=0,line_alpha=0)                
#qglyph=Quad(left='nodebox_lefts',right='nodebox_rights',bottom='nodebox_bottoms',top='nodebox_tops',fill_color='qcolor',line_alpha=0,\
#                        hatch_pattern='qhatch')#,fill_alpha=0,line_alpha=0)                
#p1.add_glyph(ptree.leaf_cds,qglyph,name='leaf_node')#'leaf_node')#self.ntype)
    


p1.x_range=Range1d(0,ptree.branch_edgecoords.boundbox.xmax,bounds='auto')

ptree.leaf_cds.add(['normal' for _ in ptree.leaf_cds.data['gbacc']],'rpnl_tfstyle')
ptree.leaf_cds.add(['4pt' for _ in ptree.leaf_cds.data['gbacc']],'rpnl_tfsize')
#ptree.leaf_cds.add(['blank' for _ in ptree.leaf_cds.data['gbacc']],'rpnl_tfsize')

textglyph=Text(x=1,y='nodebox_bottoms',text='gbacc',text_font_size='rpnl_tfsize')#,text_font_style='rpnl_tfstyle')#,text_font='Arial')#,text_align='center')#,name='leaf_node')
p2.add_glyph(ptree.leaf_cds,textglyph)
p2selquad=Quad(left=0,right=1,bottom='nodebox_bottoms',top='nodebox_tops',hatch_pattern="qhatch",hatch_alpha=0.5,fill_color=None,line_color=None)
p2.add_glyph(ptree.leaf_cds,p2selquad)

leafdict=ptree.leaf_cds.data.copy()
eckeepers=[x is not None for x in leafdict['ecs']]
ecdict={}
for key in leafdict:
    ecdict[key]=[]
    for x,keepvalue in enumerate(eckeepers):
        if keepvalue:
            ecdict[key].append(leafdict[key][x])
ec_cds=ColumnDataSource(ecdict)
ecqglyph=Quad(left=0,right=2,fill_color='blue',fill_alpha=0.5,line_alpha=0,bottom='nodebox_bottoms',top='nodebox_tops',hatch_pattern="qhatch",hatch_alpha=0.5)                
p2.add_glyph(ec_cds,ecqglyph,name='metadata')#'leaf_node')#self.ntype)

pdbkeepers=[x is not None for x in leafdict['pdbids']]
pdbdict={}
for key in leafdict:
    pdbdict[key]=[]
    for x,keepvalue in enumerate(pdbkeepers):
        if keepvalue:
            pdbdict[key].append(leafdict[key][x])
pdb_cds=ColumnDataSource(pdbdict)
pdbqglyph=Quad(left=1,right=3,fill_color='red',fill_alpha=0.5,line_alpha=0,bottom='nodebox_bottoms',top='nodebox_tops',hatch_pattern='qhatch',hatch_alpha=0.5)                
p2.add_glyph(pdb_cds,pdbqglyph,name='metadata')#'leaf_node')#self.ntype)

#stuff=CustomJS(args={'cds':ptree.leaf_cds,'dtcds':dtcds,'dtbl':dtbl},code=stuffcode)
#stuff=CustomJS(args={'cds':ptree.leaf_cds},code=stuffcode)
#p1.add_tools(TapTool(callback=stuff,names=['leaf_node']))


#    selections=new['1d']['indices']


setupdict={'accs':[]}
dtcds=ColumnDataSource(data=setupdict)
#    dtcds.on_change('data',on_change_data_source)
dtcolumns=[TableColumn(field='accs',title='GBAccession')]
dtbl=DataTable(source=dtcds,columns=dtcolumns,width=200,height=80,editable=True,selectable=True)
dtbl.source.selected.on_change('indices',tablecallback)
#layout=column(row(dtbl,gms,acw),row(p1,p2))

jscallback=CustomJS(args={'cds':ptree.leaf_cds,'dtcds':dtcds},code=stuffcode)
p1.add_tools(TapTool(callback=jscallback))#,names=['leaf_node']))

#cooldude=CustomJS(args={'cds':leaf_source},code=coolcode)

pms=MultiSelect(title='Select phylum',\
                options=[x[0] for x in collections.Counter(ptree.leaf_cds.data['phylum']).most_common()],\
                width=200,height=70)
pms.on_change('value',lambda attr,old,new:pmsfunc())

cms=MultiSelect(title='Select class',
                options=[x[0] for x in collections.Counter(ptree.leaf_cds.data['class']).most_common()],\
                width=200,height=70)
cms.on_change('value',lambda attr,old,new:cmsfunc())

gms=MultiSelect(title='Select genus',\
                options=[x[0] for x in collections.Counter(ptree.leaf_cds.data['genus']).most_common()],\
                width=200,height=70)
gms.on_change('value',lambda attr,old,new:gmsfunc())

kms=MultiSelect(title='Kingdom',\
                options=[x[0] for x in collections.Counter(ptree.leaf_cds.data['superkingdom']).most_common()],\
                width=200,height=70)
kms.on_change('value',lambda attr,old,new:kmsfunc())

sfms=MultiSelect(title='Subfams',\
                options=[x[0] for x in collections.Counter(ptree.leaf_cds.data['subfamstr']).most_common()],\
                width=150,height=70)
sfms.on_change('value',lambda attr,old,new:sfmsfunc())

tcb=Button(label='CalcTC',width=80,height=40)
tcb.on_click(tcbfunc)
#splist=list(set(ptree.leaf_cds.data['species']))
#splist.sorted(key=ptree.leaf_cds.data['species'])
#acw=AutocompleteInput(title='Organism name',completions=list(set(ptree.leaf_cds.data['species'])),width=200,height=50)
acw=AutocompleteInput(title='Organism name',\
    completions=[x[0] for x in collections.Counter(ptree.leaf_cds.data['species']).most_common()],\
    width=200,height=50)
acw.on_change('value',lambda attr,old,new:acfunc())

accacw=AutocompleteInput(title='Accession',\
    completions=ptree.leaf_cds.data['gbacc'][:],width=200,height=50)
accacw.on_change('value',lambda attr,old,new:accacfunc())


if len(sys.argv)>3:
    preselfpath=os.path.join(os.environ['SCIENCEDIR'],'GHSeqs',ghfam.upper(),sys.argv[3])
    with open(preselfpath,'r') as f:
        preselections=[x.strip() for x in f.readlines()]
        print(preselections)
else:
    preselections=[]

psb=Button(label='Show Presels',width=80,height=40)
psb.on_click(psbfunc)

#tlbl=Label(0,0,text=str(dtcds.data.accs))

tlsetupdict={'tl':[ptree.calc_treelength()]}
tlcds=ColumnDataSource(data=tlsetupdict)
tlcolumns=[TableColumn(field='tl',title='TLength')]
tltbl=DataTable(source=tlcds,columns=tlcolumns,width=100,height=100)#,editable=True,selectable=True)

tcsetupdict={'tc':[0.0]}
tccds=ColumnDataSource(data=tcsetupdict)
tccolumns=[TableColumn(field='tc',title='TCoverage')]
tctbl=DataTable(source=tccds,columns=tccolumns,width=100,height=100)#,editable=True,selectable=True)
#tctbl.source.selected.on_change('indices',tablecallback)



p2.x_range=Range1d(0,3)#ptree.branch_edgecoords.boundbox.xmax,bounds='auto')
p1.toolbar_location = 'above'
p2.toolbar_location = 'above'
p1.xaxis.visible=False
p1.xgrid.visible=False
p1.yaxis.visible=False
p1.ygrid.visible=False
p2.xaxis.visible=False
p2.xgrid.visible=False
p2.yaxis.visible=False
p2.ygrid.visible=False


#layout=column(row(column(dtbl,sfms),column(pms,cms),column(gms,acw),column(row(psb,tcb),row(tltbl,tctbl))),row(p1,p2))
layout=column(row(column(dtbl,sfms),column(pms,cms),column(gms,kms,acw),column(row(psb,tcb),row(tltbl,tctbl),row(accacw))),row(p1,p2))
#layout=row(p1,p2)
curdoc().add_root(layout)
#superq=p1.quad(left='nodebox_lefts',right='nodebox_rights',bottom='nodebox_bottoms',top='nodebox_tops',fill_alpha='qfillalpha',fill_color='qcolor',line_alpha=0,\
#                        hatch_pattern='qhatch',source=leaf_source,name='leaf_node')#,fill_alpha=0,line_alpha=0) )

#with open('st.pkl','rb') as f:
#    etetree=pickle.load(f)
#layout,p1,p2,ptree,_=bokehapp_layout.build_layout(etetree,'junk')
#p1.js_on_event('tap',callback)
##output_notebook()
#show(row(p1,p2))
#curdoc().add_root(layout)
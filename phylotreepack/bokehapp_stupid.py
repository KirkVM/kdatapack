from bokeh.plotting import reset_output,show,figure,ColumnDataSource,curdoc
from bokeh.models import HoverTool,ResetTool,BoxZoomTool,PanTool,Line,Range1d,TapTool
from bokeh.layouts import row,column
from bokeh.models.glyphs import Text,Quad,MultiLine
from bokeh.models.widgets import DataTable,TableColumn
from bokeh.models.widgets.inputs import MultiSelect,Select,AutocompleteInput,TextInput

from bokeh.models.callbacks import CustomJS
from bokeh.io import output_notebook

#from bokeh.plotting import figure,show,ColumnDataSource,output_file,curdoc,save,reset_output
#from bokeh.models import Whisker,Band,NormalHead,VeeHead,TeeHead,Range1d,Line,Label,Plot,Oval
#from bokeh.models import HoverTool,ResetTool,BoxZoomTool,LinearAxis,WheelPanTool,PanTool,Legend,LegendItem
##from bokeh.models.Scroll import %%!
#from bokeh.models.widgets import Panel, Tabs
#from bokeh.layouts import row,column
#from bokeh.io import export_png,export_svgs
#from bokeh.layouts import gridplot
import pickle

from bokeh.plotting import reset_output,show,figure,ColumnDataSource,curdoc
from bokeh.models.callbacks import CustomJS
#from bokeh.models import HoverTool,ResetTool,BoxZoomTool,PanTool,Line,Range1d
#from bokeh.layouts import row
#from bokeh.models.glyphs import Text
import sys,pickle,ete3
from ete3 import Tree
from phylotreepack import phylotree

stuffcode="""
    var data=cds.data;
    var dtdata=dtcds.data;
    const inds = cds.selected.indices;

    for (var i = 0; i < inds.length; i++) {
        if (data['frame_lws'][inds[i]]==0.5) {
            data['frame_lws'][inds[i]]=3;
            data['rpnl_tfsize'][inds[i]]="10pt";
            data['qhatch'][inds[i]]="diagonal_cross";
        }else{
            data['frame_lws'][inds[i]]=0.5;
            data['rpnl_tfsize'][inds[i]]="4pt";
            data['qhatch'][inds[i]]="blank";
        }
        dtdata['accs'][0]="cool";///.push("junk");///data['gbacc'][inds[i]])
    } 
    data.change.emit();
///    dtdata.change.emit();
    dtcds.change.emit();
    dtbl.change.emit();
"""

coolcode="""
    var data=cds.data;
    for (var i=0;i<data['phylum'].length;i++){
        for (j=0;j<cb_obj.value.length;j++){
            if (data['phylum'][i]==cb_obj.value[j]){
                data['qcolor'][i]="blue";
                data['qfillalpha'][i]=0.5;
            }else{
                data['qcolor'][i]=null;
                data['qfillalpha'][i]=0;
            }
        }
    }
    data.change.emit();
    cds.change.emit();
"""
accode="""
    var data=cds.data;
    for (var i=0;i<data['species'].length;i++){
        if (data['species'][i]==cb_obj.value){
            data['qcolor'][i]="green";
            data['qfillalpha'][i]=0.5;
        }else{
            data['qcolor'][i]=null;
            data['qfillalpha'][i]=0;
        }
    }
    data.change.emit();
    cds.change.emit();
"""

#with open('pt.pkl','rb') as f:
with open('st.pkl','rb') as f:
    etetree=pickle.load(f)

#def get_source_data(currentsource,fchange,vals):
def get_updated_values(src):
    src2=ColumnDataSource()#data=src.data)
    for k in src.data:
        src2.add(src.data[k],k)
#    del src
    return src2

hover_tool = HoverTool(names=['leaf_node'],tooltips=[    ("GB acc", "@gbacc"),('sf','@subfam'),\
                        ('phylum','@phylum'),('class','@class'),('species','@species')])
hover_tool2 = HoverTool(names=['metadata'],tooltips=[    ("GB acc", "@gbacc"),("ECs", "@ecs"),('PDBs','@pdbids')])
pheight=1100
p1 = figure(plot_width=850,plot_height=pheight,tools=[hover_tool,ResetTool(),BoxZoomTool(),PanTool()])#plot_width=1100, plot_height=700,
    #pheight=int(len(ptree.get_leaves())*1.5)
p2 = figure(plot_width=60,plot_height=pheight,tools=[hover_tool2],x_range=Range1d(0,2),y_range=p1.y_range)
ptree=phylotree.PhyloTree(etetree)
ptree.update_leafcdsdict_fromdb("GH5/GH5DB.sql")
ptree.update_leafcdsdict_fromxr("GH5/gh5ds.nc")
ptree.bokehdraw(plot=p1,rotation=0)

internal_source=ColumnDataSource(data=ptree.internal_cds.data)
leaf_source=ColumnDataSource(data=ptree.leaf_cds.data)
fglyph=MultiLine(xs='frame_xs',ys='frame_ys',line_width='frame_lws')
p1.add_glyph(leaf_source,fglyph)#,name='leaf_node')
p1.add_glyph(internal_source,fglyph)#,name='internal_node')

#qglyph=Quad(left='nodebox_lefts',right='nodebox_rights',bottom='nodebox_bottoms',top='nodebox_tops',fill_color='qcolor',line_alpha=0,\
#                        hatch_pattern='qhatch')#,fill_alpha=0,line_alpha=0)                
superq=p1.quad(left='nodebox_lefts',right='nodebox_rights',bottom='nodebox_bottoms',top='nodebox_tops',fill_alpha='qfillalpha',fill_color='qcolor',line_alpha=0,\
                        hatch_pattern='qhatch',source=leaf_source,name='leaf_node')#,fill_alpha=0,line_alpha=0) )
#p1.add_glyph(leaf_source,qglyph,name='leaf_node')#'leaf_node')#self.ntype)

p1.x_range=Range1d(0,ptree.branch_edgecoords.boundbox.xmax,bounds='auto')
leaf_source.add(['normal' for _ in ptree.leaf_cds.data['gbacc']],'rpnl_tfstyle')
leaf_source.add(['4pt' for _ in ptree.leaf_cds.data['gbacc']],'rpnl_tfsize')

textglyph=Text(x=1,y='nodebox_bottoms',text='gbacc',text_font_size='rpnl_tfsize')#,text_font_style='rpnl_tfstyle')#,text_font='Arial')#,text_align='center')#,name='leaf_node')
p2.add_glyph(leaf_source,textglyph)
p2selquad=Quad(left=0,right=1,bottom='nodebox_bottoms',top='nodebox_tops',hatch_pattern="qhatch",hatch_alpha=0.5,fill_color=None,line_color=None)
p2.add_glyph(leaf_source,p2selquad)

leafdict=leaf_source.data.copy()
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
#
#ecdict={}
#ecaccs=ptree.leaf_cds

#    ptree.qglyph.subscribed_events=['selected','tap']
#    taptool=p1.select(type=TapTool)
#    taptool.callback=callback

#    p3=figure(plotwidth=200,plot_height=pheight/10)
setupdict={'accs':['silly1']}
dtcds=ColumnDataSource(data=setupdict)
#    dtcds.on_change('data',on_change_data_source)
dtcolumns=[TableColumn(field='accs',title='GBAccession')]
dtbl=DataTable(source=dtcds,columns=dtcolumns,width=200,height=50,editable=True,selectable=True)



#def stupid(attr,old,new):
#    global dtcds
#    dtcds.data['accs'].append('fdsfsfsdf')
#    ptree.leaf_cds['frame_lws'][0]=3
#    print('cool')

#ptree.qglyph.subscribed_events=['tap']
#p1.js_on_event('tap',CustomJS(args={'cds':ptree.leaf_cds,'dtcds':dtcds,'dtbl':dtbl},code="""
#    var data=cds.data;
#    var dtdata=dtcds.data;
#    const inds = cds.selected.indices;
#
#    for (var i = 0; i < inds.length; i++) {
#        if (data['frame_lws'][inds[i]]==1) {
#            data['frame_lws'][inds[i]]=3
#        }else{
#            data['frame_lws'][inds[i]]=1
#        }
#        dtdata['accs'][0]="cool";///.push("junk");///data['gbacc'][inds[i]])
#    } 
#    data.change.emit();
#    dtcds.change.emit();
#    dtbl.change.emit();
#"""))


#ptree.leaf_cds.data.on_change('frame_lws',stupid)
#ptree.fglyph.on_change('line_width',stupid)
#CustomJS(args=dict(dtcds=dtcds,dtbl=dtbl),code="""
#    var dtdata=dtcds.data;
#    dtdata['accs'][0]='stupid';
#    dtcds.change.emit();
#    dtbl.change.emit();
#"""))

stuff=CustomJS(args={'cds':leaf_source,'dtcds':dtcds,'dtbl':dtbl},code=stuffcode)
p1.add_tools(TapTool(callback=stuff,names=['leaf_node']))

cooldude=CustomJS(args={'cds':leaf_source},code=coolcode)
gms=MultiSelect(title='Select phylum',callback=cooldude,options=list(set(leaf_source.data['phylum'])),width=200,height=70)

acdude=CustomJS(args={'cds':leaf_source},code=accode)
acw=AutocompleteInput(title='Organism name',callback=acdude,completions=list(set(leaf_source.data['species'])),width=200,height=50)
#acw=TextInput(title='Organism name',callback=acdude,width=200,height=50)
#import numpy as np
#def adjuster(attr,old,new):
#    global leaf_source
#    global p1,superq
#    newsource=ColumnDataSource(data=leaf_source.data)
#    for x,phylum in enumerate(leaf_source.data['phylum']):
#        if phylum in new:
#            leaf_source.data['qcolor'][x]='blue'
#            leaf_source.data['qfillalpha'][x]=0.5
##            leaf_source.data['nodebox_lefts'][x]=np.nan#[np.nan for _ in range(len(leaf_source.data))]#leaf_source.data['nodebox_lefts'][x]-0.1
#        else:
#            leaf_source.data['qcolor'][x]=None
#            leaf_source.data['qfillalpha'][x]=0
##            leaf_source.data['nodebox_lefts'][x]=np.nan#[np.nan for _ in range(len(leaf_source.data))]#leaf_source.data['nodebox_lefts'][x]-0.1
#            leaf_source.data['nodebox_lefts'][x]-=0.1
    #    print(leaf_source.data['qcolor'])
#    newsource=get_updated_values(leaf_source)
#    superq.source=newsource#leaf_source.data['qcolor'][:]
#    newsource.data['qcolor']=leaf_source.data['qcolor'][:]
#    newsource.data['qfillalpha']=leaf_source.data['qcolor'][:]
#    superq=p1.quad(left='nodebox_lefts',right='nodebox_rights',bottom='nodebox_bottoms',top='nodebox_tops',fill_alpha='qfillalpha',fill_color='qcolor',line_alpha=0,\
#                        hatch_pattern='qhatch',source=newsource,name='leaf_node')#,fill_alpha=0,line_alpha=0) )
#    leaf_source.data=newsource.data.copy()
#
##
#    hover_tool = HoverTool(names=['leaf_node'],tooltips=[    ("GB acc", "@gbacc"),('sf','@subfam'),\
#                        ('phylum','@phylum'),('class','@class'),('species','@species')])
#    p1.tools=[hover_tool,ResetTool(),BoxZoomTool(),PanTool()]#plot_width=1100, plot_height=700,
##    p1.toolbar_location = 'above'
#    replot_superq(superq,leaf_source)
#    superq.source=leaf_source
#    p1.add_glyph(leaf_source,superq,name='leaf_node')#'leaf_node')#self.ntype)
#
#gms=MultiSelect(title='Select phylum',options=list(set(leaf_source.data['phylum'])),width=200,height=70)
#gms.on_change('value',adjuster)




#    p1.js_on_event('tap',callback)

#    tap_tool=TapTool(names=['leaf_node'],callback=stuff)

    #p2.text(x=0,y='nodebox_bottoms',text='gbacc')
    #for leaf in ptree.get_leaves():
    #l=Line(x='junk',y='stupid')
    #p1.add_glyph(cds,l)
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
layout=column(row(dtbl,gms,acw),row(p1,p2))
#layout=column(row(dtbl),row(p1,p2))
#output_notebook()
#show(layout)
curdoc().add_root(layout)
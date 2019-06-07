from bokeh.plotting import reset_output,show,figure,ColumnDataSource,curdoc
from bokeh.models import HoverTool,ResetTool,BoxZoomTool,PanTool,Line,Range1d,TapTool
from bokeh.layouts import row,column
from bokeh.models.glyphs import Text
from bokeh.models.widgets import DataTable,TableColumn

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
from . import phylotree


#stuff=CustomJS(args={sourcecode="""
#console.log('Tap event')
#
#""")

with open('st.pkl','rb') as f:
    etetree=pickle.load(f)

hover_tool = HoverTool(names=['leaf_node'],tooltips=[    ("GB acc", "@gbacc"),('sf','@subfam')])#, ("species","@species"),
pheight=1100
p1 = figure(plot_width=850,plot_height=pheight,tools=[hover_tool,ResetTool(),BoxZoomTool(),PanTool()])#plot_width=1100, plot_height=700,
    #pheight=int(len(ptree.get_leaves())*1.5)
p2 = figure(plot_width=60,plot_height=pheight,tools=[],x_range=Range1d(0,2),y_range=p1.y_range)
ptree=phylotree.PhyloTree(etetree)
ptree.bokehdraw(plot=p1,rotation=0)
p1.x_range=Range1d(0,ptree.branch_edgecoords.boundbox.xmax,bounds='auto')
textglyph=Text(x=0,y='nodebox_bottoms',text='gbacc',text_font='Arial',text_font_size='7pt')#,text_align='center')#,name='leaf_node')
p2.add_glyph(ptree.leaf_cds,textglyph)
#    ptree.qglyph.subscribed_events=['selected','tap']
#    taptool=p1.select(type=TapTool)
#    taptool.callback=callback

#    p3=figure(plotwidth=200,plot_height=pheight/10)
setupdict={'accs':['silly1']}
dtcds=ColumnDataSource(data=setupdict)
#    dtcds.on_change('data',on_change_data_source)
dtcolumns=[TableColumn(field='accs',title='GBAccession')]
dtbl=DataTable(source=dtcds,columns=dtcolumns,width=200,height=50,editable=True,selectable=True)

ptree.qglyph.js_on_event('tap',CustomJS(args=dict(dtcds=dtcds,dtbl=dtbl),code="""
    var dtdata=dtcds.data;
    dtdata['accs'][0]='stupid';
    dtcds.change.emit();
    dtbl.change.emit();
"""))
#    stuff=CustomJS(args={'cds':ptree.leaf_cds,'dtcds':dtcds,'dtbl':dtbl},code=stuffcode)
#    p1.add_tools(TapTool(callback=stuff,names=['leaf_node']))

#    p1.js_on_event('tap',callback)

#    tap_tool=TapTool(names=['leaf_node'],callback=stuff)

    #p2.text(x=0,y='nodebox_bottoms',text='gbacc')
    #for leaf in ptree.get_leaves():
    #l=Line(x='junk',y='stupid')
    #p1.add_glyph(cds,l)
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
layout=column(dtbl,row(p1,p2))
return layout
#output_notebook()
#show(layout)

#output_notebook()
#show(row(p1,p2))
from bokeh.plotting import reset_output,show,figure,ColumnDataSource,curdoc
from bokeh.models import HoverTool,ResetTool,BoxZoomTool,PanTool,Line,Range1d,TapTool
from bokeh.layouts import row,column
from bokeh.models.glyphs import Text

from bokeh.models.callbacks import CustomJS

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
stuffcode="""
    var data=cds.data;
    const inds = cds.selected.indices;

    for (var i = 0; i < inds.length; i++) {
        if (data['frame_lws'][inds[i]]==1) {
            data['frame_lws'][inds[i]]=3
        }else{
            data['frame_lws'][inds[i]]=1
        }
///        data['frame_lws'][inds[i]] = 4
///        ym += d['y'][inds[i]]
    } 
///    data['frame_lws'][0]=4;
///    data['gbacc'][0]='dfasf';
    data.change.emit();
"""
#callback = CustomJS(code="""
#// the event that triggered the callback is cb_obj:
#// The event type determines the relevant attributes
#console.log('Tap event occurred at x-position: ' + cb_obj.x)
#""")
def callback(attr,old,new):
    print('junk')

def build_layout(etetree,dbpathstr):
    hover_tool = HoverTool(names=['leaf_node'],tooltips=[    ("GB acc", "@gbacc"),('sf','@subfam')])#, ("species","@species"),
    pheight=1100
    p1 = figure(plot_width=850,plot_height=pheight,tools=[hover_tool,ResetTool(),BoxZoomTool(),PanTool()])#plot_width=1100, plot_height=700,
    #pheight=int(len(ptree.get_leaves())*1.5)
    p2 = figure(plot_width=100,plot_height=pheight,tools=[],x_range=Range1d(0,2),y_range=p1.y_range)
    ptree=phylotree.PhyloTree(etetree)
    ptree.bokehdraw(plot=p1,rotation=0)
    p1.x_range=Range1d(0,ptree.branch_edgecoords.boundbox.xmax,bounds='auto')
    textglyph=Text(x=0,y='nodebox_bottoms',text='gbacc',text_font='Arial',text_font_size='7pt')#,name='leaf_node')
    p2.add_glyph(ptree.leaf_cds,textglyph)
#    ptree.qglyph.subscribed_events=['selected','tap']
#    ptree.qglyph.on_change('tap',callback) 
#    taptool=p1.select(type=TapTool)
#    taptool.callback=callback

    stuff=CustomJS(args={'cds':ptree.leaf_cds},code=stuffcode)
    p1.add_tools(TapTool(callback=stuff,names=['leaf_node']))

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
    layout=row(p1,p2)
    return layout,p1,p2,ptree

#output_notebook()
#show(row(p1,p2))
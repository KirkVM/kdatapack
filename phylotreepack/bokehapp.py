from bokeh.plotting import reset_output,show,figure,ColumnDataSource,curdoc
#from bokeh.models import HoverTool,ResetTool,BoxZoomTool,PanTool,Line,Range1d
#from bokeh.layouts import row
#from bokeh.models.glyphs import Text
import sys,pickle,ete3
from ete3 import Tree
from phylotreepack import bokehapp_layout
with open('st.pkl','rb') as f:
    etetree=pickle.load(f)
layout=bokehapp_layout.build_layout(etetree,'junk')
##output_notebook()
#show(row(p1,p2))
curdoc().add_root(layout)
from bokeh.plotting import reset_output,show,figure,ColumnDataSource,curdoc
from bokeh.models import HoverTool,ResetTool,BoxZoomTool,PanTool,Line
from phylotreepack import phylotree
import pickle
from ete3 import Tree

def update_leafcds_fromdb(ptree,dbcursor,fields=['pdbids','ecs','subfam','extragbs'],searchby='gbacc'):
    assert (searchby=='gbacc')
    for lnode in ptree.get_leaves():
        lname=lnode.name
        dbcursor.execute('''SELECT * FROM CAZYSEQDATA WHERE acc=(?)''',(lname,))
        row=dbcursor.fetchone()
        if row is None: continue
        for field in fields:
            if row[field] is None: continue
            lnode.decoration_dict[field]=row[field]

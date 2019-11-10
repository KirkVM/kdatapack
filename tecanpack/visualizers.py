from bokeh.plotting import show,output_notebook,figure
from bokeh.models import ColumnDataSource,HoverTool,Range1d
from bokeh import palettes
def plot_controls(df,detection):
    df=df.dropna(subset=['standardname'])
    df=df[df.detection==detection.upper()]
    df['strdate']=df.expdate.apply(str)
#    sdf_dict={}
    cds_dict={}
    for expdate in df.expdate.unique():
        #sdf_dict[expdate]=df[df.expdate==expdate]
        cds_dict[expdate]=ColumnDataSource.from_df(df[df.expdate==expdate])#sdf_dict[expdate])

    p=figure()
    ht=HoverTool(tooltips=[('date','@strdate'),('plateid','@plateid')])
    for dnum,expdate in enumerate(cds_dict):
        p.circle('standardconc','measurement',legend='strdate',size=10,color=palettes.Category10[10][dnum],source=cds_dict[expdate])
    p.tools=[ht]
    output_notebook()
    show(p)
#for expdate in bcasdf_dict:
#    
#    plt.plot('measurement','standardconc',data=bcasdf_dict[expdate])
#plt.legend(bcasdf_dict.keys())
#s1.count()
#s1data=s1.to_dict(orient='list')
def panel_plot(df):
    p=figure()
    output_notebook()
    show(p)

from bokeh.models import BoxZoomTool,ResetTool
def plot_eblanks(df,xcolname='econc_mgmL',ycolname='measurement'):
    df['strdate']=df.expdate.apply(str)
#    sdf_dict={}
    cds_dict={}
    for egroup in df.groupby(by=['ename','expdate']):
        grpname=''.join([str(x)+' ' for x in egroup[0]])[:-1]
        cds_dict[grpname]=ColumnDataSource.from_df(egroup[1])#df[df.expdate==expdate])#sdf_dict[expdate])

    p=figure(plot_width=800)
    ht=HoverTool(tooltips=[('enzyme','@ename'),('M','@econc_molar'),('date','@strdate'),('plateid','@plateid'),\
                            ('thing','@predvlp_dlnfactor')])
    for dnum,grpname in enumerate(cds_dict):
        didx=dnum%10
        p.circle(xcolname,ycolname,legend=grpname,size=10,color=palettes.Category20[20][didx],source=cds_dict[grpname])
    p.tools=[ht,BoxZoomTool(),ResetTool()]
#    p.tools=[ht,"box_zoom"]
    p.x_range=Range1d(-0.1,4.5)
    output_notebook()
    show(p)
#
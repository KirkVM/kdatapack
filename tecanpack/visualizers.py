from bokeh.plotting import show,output_notebook,figure
from bokeh.models import ColumnDataSource,HoverTool
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

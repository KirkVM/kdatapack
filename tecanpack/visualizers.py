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

import panel as pn
import datetime
class ActivityPanel:
    def __init__(self,initdf):
        self.initdf=initdf
        self.df=initdf.copy()
        self.tabtups=None
        self.dtabs=None
        self.sels=None
        self.selcol=None
        self.p=None
        self.test='cool'
        self.panel_layout=None
    def date_tabs(self):
        expdates=[str(x) for x in self.df.expdate.unique()]
        tvals=[(str(x),pn.Row(pn.widgets.StaticText(value=str(expdates[x])),width=300 ) ) for x in range(len(expdates))]
        tabtups=[]
        expdates=[x for x in self.df.expdate.unique()]
        info_fields={'expdate':'date','ename':'enzymes','sname':'substrates','experimenter':'person'}
        for expnum,expdate in enumerate(expdates):
            ifields=[]
            for ifkey in info_fields.keys():
        
                vals=[str(x) for x in self.df[self.df.expdate==expdate][ifkey].unique()]
                if len(vals)==0:
                    rowvalstr='-'
                else:
                    rowvalstr=''.join(f'{x},' for x in vals[:-1])
                    rowvalstr+=f'{vals[-1]}'
                r_stext=pn.widgets.StaticText(value=f'{info_fields[ifkey]}: {rowvalstr}')
                ifields.append(r_stext)
            tabtups.append(( str(expnum),pn.Column(*ifields,width=200)  ))#r1_stext,r2_stext,width=300)   ))
        self.dtabs=pn.Tabs(*tabtups)
    def selector_windows(self):
        self.sels=[]
        info_fields={'expdate':'date','ename':'enzymes','sname':'substrates','experimenter':'person'}
        for ifkey,ifname in info_fields.items():
            sel_options=['All']
            sel_options.extend([str(x) for x in self.df[ifkey].unique()])
            sel=pn.widgets.Select(name=ifname,options=sel_options)
            sel.param.watch(self.selector_callback,['value'],onlychanged=True)
            self.sels.append(sel)
        self.selcol=pn.Column(*self.sels,width=200)
    def theplot(self):
        self.p=figure()
        self.p.plot_width=750
        self.p.plot_height=450
  #      output_notebook()
#        show(p)
    def selector_callback(self,*events):
        info_fields={'expdate':'date','ename':'enzymes','sname':'substrates','experimenter':'person'}
        del self.df
        self.df=self.initdf
#        for ifkey,ifname in info_fields.items():
#                if event.obj.name==ifname:

        for sel in self.sels:
            if sel.value!='All':
                for ifkey,ifname in info_fields.items():
                    if ifname==sel.name:
                        if ifkey=='expdate':
                            selectorvalue=datetime.date(*[int(x) for x in sel.value.split('-')])
                        else:
                            selectorvalue=sel.value
                        self.df=self.df[self.df[ifkey]==selectorvalue]
        del self.dtabs
        self.date_tabs()
        return self.panel_layout
#            self.test=
#        for event in events:
#            self.test=event#[event.name,event.new,event.options,event.value]
#            for ifkey,ifname in info_fields.items():
#                if event.obj.name==ifname:
#                    self.df=self.df[self.df[ifkey]==event.new]




#
#                if event.obj.value=='All':
#                    self.df=self.initdf

    def visualize(self):
        self.panel_layout=pn.Row(self.selcol,self.p,self.dtabs)
        return self.panel_layout

#    def selector_callback(self,*events):
#        pass
#        s
#    tabs=pn.Tabs(*tabtups)#,sizing_mode='stretch_b
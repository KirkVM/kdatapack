import datetime
import param
import panel as pn
import pandas as pd
from bokeh.plotting import show,output_notebook,figure
from bokeh.models import ColumnDataSource,HoverTool

from tecanpack import readers
#plateset=readers.load_tecandata('allfiles.yml','kirk')#refresh_all=True)
#alldf=plateset.get_df()


alldf=pd.read_pickle('kirkdf.pkl')
class ActivityPanel(param.Parameterized):
    enames=param.Selector(objects=[str(x) for x in alldf['ename'].unique()])
    snames=param.Selector(objects=[str(x) for x in alldf['sname'].unique()])
    expdates=param.Selector(objects=[str(x) for x in alldf['expdate'].unique()])
    filter_enames=param.Boolean(False)
    filter_snames=param.Boolean(False)
    filter_expdates=param.Boolean(False)
    expdf=None
    plotdf=None
    seldf=alldf.copy()
    stupid=0
#    enames.param.watch(self.selector_callback,['value'],onlychanged=True)
#    pname_dict={'enzymes':'ename','exp_date':'expdate'}
#    cool='funny'
#    def selector_callback(self,*events):
#        self.stupid=events[0]
#        for event in events:
#            self.stupid=event#[event.name,event.new,event.options,event.value]
#            if event.obj.name=='enames':
#                self.seldf=self.seldf[self.seldf['ename']==event.new]
#        self.param.enames.objects=[str(x) for x in self.seldf['ename'].unique()]
#        self.param.snames.objects=[str(x) for x in self.seldf['sname'].unique()]
#        self.param.expdates.objects=[str(x) for x in self.seldf['expdate'].unique()]
#            for ifkey,ifname in info_fields.items():
#                if event.obj.name==ifname:
#                    self.df=self.df[self.df[ifkey]==event.new]
    @param.depends('enames',watch=True)
    def update_seldf(self):
        self.seldf=self.seldf[self.seldf['ename']==self.enames]
        self.param.enames.objects=[str(x) for x in self.seldf['ename'].unique()]
        self.param.snames.objects=[str(x) for x in self.seldf['sname'].unique()]
        self.param.expdates.objects=[str(x) for x in self.seldf['expdate'].unique()]


#    @param.depends('enames','snames','expdates',watch=True)
    def sel_viewer(self):
#        self.seldf=self.seldf[self.seldf['ename']==self.enames]
       #gspec = pn.GridSpec(sizing_mode='stretch_both', max_height=800,max_width=200,height_policy='max')
        #gspec = pn.GridSpec(sizing_mode='scale_both')#max_width=50,max_height=20,width_policy='auto')#='stretch_width')#width=200)#,align='center')#sizing_mode='stretch_both', max_height=800,max_width=200,height_policy='max')
        gspec = pn.GridSpec(height=500,width=180,max_width=10,max_height=20,width_policy='max',height_policy='max')#max_width=50,max_height=20,width_policy='auto')#='stretch_width')#width=200)#,align='center')#sizing_mode='stretch_both', max_height=800,max_width=200,height_policy='max')
#        gspec = pn.GridSpec(max_width=50,width_policy='max')#max_width=50,max_height=20,width_policy='auto')#='stretch_width')#width=200)#,align='center')#sizing_mode='stretch_both', max_height=800,max_width=200,height_policy='max')
#        self.param.enames.objects=[str(x) for x in self.seldf['ename'].unique()]
#        self.param.snames.objects=[str(x) for x in self.seldf['sname'].unique()]
#        self.param.expdates.objects=[str(x) for x in self.seldf['expdate'].unique()]
        gspec[0:1,0:1] = pn.Column(pn.Param(self.param.filter_enames,
                        widgets={'filter_enames':{'type':pn.widgets.Checkbox,'name':'filter'}}),
                        sizing_mode='stretch_both',margin=(10,3,0,0))
#        gspec[0:1,0:1] = pn.Column(self.param.filter_enames,sizing_mode='stretch_both',margin=(10,3,0,0))
        gspec[1:2,0:1] = pn.Column(pn.Param(self.param.filter_snames,
                        widgets={'filter_snames':{'type':pn.widgets.Checkbox,'name':'filter'}}),
                        #widgets={'filter_snames':{'type':pn.widgets.Toggle}}),
                        sizing_mode='stretch_both',margin=(10,3,0,0))
        gspec[2:3,0:1] = pn.Column(pn.Param(self.param.filter_expdates,
                        widgets={'filter_expdates':{'type':pn.widgets.Checkbox,'name':'filter'}}),
                        #widgets={'filter_expdates':{'type':pn.widgets.Toggle}}),
                        sizing_mode='stretch_both',margin=(10,3,0,0))
        gspec[0:1,1:5]= pn.Column(pn.Param(self.param.enames,
                         widgets={'enames':{'type':pn.widgets.Select,'name':'Enzyme'}}),
#                         widgets={'enames':{'type':pn.widgets.Select,'name':'enames'}}),
                        sizing_mode='stretch_both')
#        gspec[0:1,1:5]= pn.Column(self.param.enames,sizing_mode='stretch_both')
        gspec[1:2,1:5]= pn.Column(pn.Param(self.param.snames,
                         widgets={'snames':{'type':pn.widgets.Select,'name':'Substrate'}}),
                         #widgets={'snames':{'type':pn.widgets.Select}}),
                        sizing_mode='stretch_both')
        gspec[2:3,1:5]= pn.Column(pn.Param(self.param.expdates,
                         widgets={'expdates':{'type':pn.widgets.Select,'name':'Exp Date'}}),
                         #widgets={'expdates':{'type':pn.widgets.Select}}),
                         sizing_mode='stretch_both')
        gspec[3:10,:5]=pn.Spacer(margin=0)
        return gspec

#    @param.depends('enames','snames','expdates')#,watch=True)
    def exp_viewer(self):
        tabtups=[]
        expdf=self.seldf.copy()
        if self.filter_enames:
            expdf=expdf[expdf['ename']==self.enames] 
            self.stupid='dumb'
        if self.filter_snames:
            expdf=expdf[expdf['sname']==self.snames] 
            self.stupid='rad'
        if self.filter_expdates:
            dt_expdate=datetime.date(*[int(x) for x in self.expdates.split('-')])
            expdf=expdf[expdf['expdate']==dt_expdate] 

        expdates=[x for x in expdf.expdate.unique()]
        for expnum,expdate in enumerate(expdates):
            ifields=[]
            for pname in ['ename','sname','expdate']:
                vals=[str(x) for x in expdf[expdf.expdate==expdate][pname].unique()]
                if len(vals)==0:
                    rowvalstr='-'
                else:
                    rowvalstr=''.join(f'{x},' for x in vals[:-1])
                    rowvalstr+=f'{vals[-1]}'
                r_stext=pn.widgets.StaticText(value=f'{pname}: {rowvalstr}')
                ifields.append(r_stext)
            tabtups.append(( str(expnum),pn.Column(*ifields,width=200)  ))#r1_stext,r2_stext,width=300)   ))
        dtabs=pn.Tabs(*tabtups)
#        self.expdf=expdf
        return dtabs

    def theplot(self):
        self.p=figure()
        self.p.plot_width=750
        self.p.plot_height=450
    
    def visualize(self):
        self.panel_layout=pn.Row(self.sel_viewer,self.exp_viewer)
        return self.panel_layout
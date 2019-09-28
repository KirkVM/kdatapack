import datetime
import param
import panel as pn
import pandas as pd
from bokeh.plotting import show,output_notebook,figure
from bokeh.models import ColumnDataSource,HoverTool
from functools import partial

from tecanpack import readers
#plateset=readers.load_tecandata('allfiles.yml','kirk')#refresh_all=True)
#alldf=plateset.get_df()
alldf=pd.read_pickle('kirkdf.pkl')

class ActivityPanel(param.Parameterized):
    initdf=param.DataFrame(alldf.copy())
    seldf=param.DataFrame(alldf.copy())
    viewdf=param.DataFrame()
    exploredf=param.DataFrame()
    enames=param.Selector(objects=[str(x) for x in alldf['ename'].unique()])
    snames=param.Selector(objects=[str(x) for x in alldf['sname'].unique()])
    expdates=param.Selector(objects=[str(x) for x in alldf['expdate'].unique()])
    filter_enames=param.Boolean(False)
    filter_snames=param.Boolean(False)
    filter_expdates=param.Boolean(False)

    expview_dates=param.List(['stuff'])
    expview_enames=param.List(['dummy'])
    expview_snames=param.List(['dummy'])
    tabbies=param.Dynamic()
#    dtabs=pn.Tabs()
    cool='dumb'

    def get_plotdf(self,event,**kwargs):
        filter_num=0
        for k,v in kwargs.items():
            if filter_num==0:
                self.exploredf=self.initdf[self.initdf[k]==v]
            else:
                self.exploredf=self.exploredf[self.exploredf[k]==v]
            filter_num+=1


    def get_tabs(self,ptuples):
        panels=[]
        self.viewdf=self.initdf[self.initdf.expdate.isin(self.seldf.expdate.unique())]
        tabwidgs=[]#[] for _ in range(len(ptuples))]
        for dt,dtdf in self.viewdf.groupby('expdate'):
            dscrptext=pn.widgets.StaticText(value=f'Date: {str(dt)}\n',width=200)
            for ptidx,ptuple in enumerate(ptuples):
                pvals=list(dtdf[ptuple[0]].unique())
                pstr_list=''.join(f'{x},' for x in pvals[:-1])
                pstr_list+=f'{pvals[-1]}'
                dscrptext.value+=f'{ptuple[1]}: {pstr_list}'
                if ptidx+1<len(ptuples): dscrptext.value+='\n'
            explore_button=pn.widgets.Button(name=f'button{dt}',button_type='primary',width=100)
            action_with_arg = partial(self.get_plotdf,expdate=dt)
            explore_button.on_click(action_with_arg)
            tabwidgs.append([dscrptext,explore_button])
        return [pn.Column(*tabwidgs[x],name=str(x)) for x in range(len(tabwidgs))]


    @param.depends('enames',watch=True)
    def update_seldf(self):
        self.seldf=self.seldf[self.seldf['ename']==self.enames]
        self.param.enames.objects=[str(x) for x in self.seldf['ename'].unique()]
        self.param.snames.objects=[str(x) for x in self.seldf['sname'].unique()]
        self.param.expdates.objects=[str(x) for x in self.seldf['expdate'].unique()]

#        self.expview_dates=param.List([])
#        self.expview_enames=param.List([])
#        self.expview_snames=param.List([])
    
#    @param.depends('seldf',watch=True)
#    def view_updater(self):
#        dategrps=self.seldf.groupby('expdate')
#        tabnum=1
#        for dt,dtdf in dategrps:
#            for pname in ['ename','sname','expdate']:
#                vals=[str(x) for x in dtdf[pname].unique()]
#                if pname=='expdate':
#                    self.expview_dates.append(vals)
#                if pname=='ename':
#                    self.expview_enames.append(vals)
#                if pname=='sname':
#                    self.expview_snames.append(vals)


    @param.depends('enames',watch=True)
    def exp_updater(self):
        #self.tabbies=self.tabtext(self,'ename')
        self.tabbies=self.get_tabs([['ename','Enzymes'],['sname','Substrates'],['econc','EnzConc'],\
                                  ['sconc','SConc'],['rxnph','pH'],['rxntemp','temp']])
        self.dtabs.active=0
        self.dtabs.objects=self.tabbies

    def sel_view(self):
        gspec = pn.GridSpec(height=500,width=180,max_width=10,max_height=20,width_policy='max',height_policy='max')#max_width=50,max_height=20,width_policy='auto')#='stretch_width')#width=200)#,align='center')#sizing_mode='stretch_both', max_height=800,max_width=200,height_policy='max')
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
 
#    def exp_view(self):
#        self.tabbies=self.tabtext(self,[['ename','Enzymes']])
#        self.dtabs=pn.Tabs(*self.tabbies)
#        return self.dtabs
    def exp_view(self):
        self.tabbies=self.get_tabs([['ename','Enzymes'],['sname','Substrates'],['econc','EnzConc'],\
                                  ['sconc','SConc'],['rxnph','pH'],['rxntemp','temp']])
        self.dtabs=pn.Tabs(*self.tabbies)
        return self.dtabs


    def theplot(self):
        self.p=figure()
        self.p.plot_width=750
        self.p.plot_height=450
    
    def visualize(self):
        self.panel_layout=pn.Row(self.sel_view(),self.exp_view())
        return self.panel_layout
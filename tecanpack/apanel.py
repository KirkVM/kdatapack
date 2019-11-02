import datetime,collections
import param
import panel as pn
import pandas as pd
from bokeh.plotting import show,output_notebook,figure
from bokeh.models import ColumnDataSource,HoverTool
from functools import partial

from tecanpack import readers
from tecanpack import tecandata
#plateset=readers.load_tecandata('allfiles.yml','kirk')#refresh_all=True)
#alldf=plateset.get_df()
alldf=pd.read_pickle('kirkdf.pkl')

class ActivityPanel(param.Parameterized):
    initdf=param.DataFrame(alldf.copy())
    seldf=param.DataFrame(alldf.copy())
    viewdf=param.DataFrame()
    exploredf=param.DataFrame()
#    enames=param.Selector(objects=[str(x) for x in alldf['ename'].unique()])
    enames=param.Selector(objects= [x[0] for x in collections.Counter(
                                    [str(x) for x in alldf['ename'].unique()]).most_common()])
    snames=param.Selector(objects= [x[0] for x in collections.Counter(
                                    [str(x) for x in alldf['sname'].unique()]).most_common()])
    #snames=param.Selector(objects=[str(x) for x in alldf['sname'].unique()])
    expdates=param.Selector(objects=[str(x) for x in alldf['expdate'].unique()])
    filter_enames=param.Boolean(False)
    filter_snames=param.Boolean(False)
    filter_expdates=param.Boolean(False)
    prev_filter_enames=False
    prev_filter_snames=False
    prev_filter_expdates=False
    expview_dates=param.List(['stuff'])
    expview_enames=param.List(['dummy'])
    expview_snames=param.List(['dummy'])
    tabbies=param.Dynamic()

    plotview_exps=param.Selector()
    plotview_enames=param.Selector()
    plotview_snames=param.Selector()
    plotview_xvariable=param.Selector()
    plotview_yvariable=param.Selector()
#    dtabs=pn.Tabs()
#    cool='dumb'

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
            explore_button=pn.widgets.Button(name='explore',button_type='primary',width=100)
            action_with_arg = partial(self.get_plotdf,expdate=dt)
            explore_button.on_click(action_with_arg)
            tabwidgs.append([dscrptext,explore_button])
        return [pn.Column(*tabwidgs[x],name=str(x)) for x in range(len(tabwidgs))]


    @param.depends('filter_enames','filter_snames','filter_expdates',watch=True)
    def update_seldf_filters(self):
        #did filter_enames change to False?
        an_action=False
        if self.filter_enames==False and self.prev_filter_enames:
            an_action=True
            add_df=self.initdf.copy()
            if self.filter_snames:
                add_df=add_df[add_df['sname']==self.snames]
            if self.filter_expdates:
                add_df=add_df[add_df['expdate']==datetime.date(*[int(x) for x in self.expdates.split('-')])]
            #self.seldf=self.seldf.merge(add_df,how='left',left_index=True,right_index=True)
            add_df=add_df[add_df.index.isin(self.seldf.index)==False]
            self.seldf=pd.concat([self.seldf,add_df])#,verify_integrity=True)
        #did filter_snames change to False?
        if self.filter_snames==False and self.prev_filter_snames:
            an_action=True
            add_df=self.initdf.copy()
            if self.filter_enames:
                add_df=add_df[add_df['ename']==self.enames]
            if self.filter_expdates:
                add_df=add_df[add_df['expdate']==datetime.date(*[int(x) for x in self.expdates.split('-')])]
            add_df=add_df[add_df.index.isin(self.seldf.index)==False]
            self.seldf=pd.concat([self.seldf,add_df])#,verify_integrity=True)
        #did filter_expdates change to False?
        if self.filter_expdates==False and self.prev_filter_expdates:
            add_df=self.initdf.copy()
            if self.filter_enames:
                add_df=add_df[add_df['ename']==self.enames]
            if self.filter_snames:
                add_df=add_df[add_df['sname']==self.snames]
            add_df=add_df[add_df.index.isin(self.seldf.index)==False]
            self.seldf=pd.concat([self.seldf,add_df])#,verify_integrity=True)

        self.param.enames.objects=[x[0] for x in collections.Counter(
                                    [str(x) for x in self.seldf['ename'].unique()]).most_common()]
        self.param.snames.objects=[x[0] for x in collections.Counter(
                                    [str(x) for x in self.seldf['sname'].unique()]).most_common()]
        self.param.expdates.objects=[str(x) for x in self.seldf['expdate'].unique()]

        self.prev_filter_enames=self.filter_enames
        self.prev_filter_snames=self.filter_snames
        self.prev_filter_expdates=self.filter_expdates

    @param.depends('enames','snames','expdates',watch=True)
    def update_seldf(self):
        if self.filter_enames:
            self.seldf=self.seldf[self.seldf['ename']==self.enames]
        if self.filter_snames:
            self.seldf=self.seldf[self.seldf['sname']==self.snames]
        if self.filter_expdates:
            self.seldf=self.seldf[self.seldf['expdate']==datetime.date(*[int(x) for x in self.expdates.split('-')])]

        self.param.enames.objects=[x[0] for x in collections.Counter(
                                    [str(x) for x in self.seldf['ename'].unique()]).most_common()]
        self.param.snames.objects=[x[0] for x in collections.Counter(
                                    [str(x) for x in self.seldf['sname'].unique()]).most_common()]
        self.param.expdates.objects=[str(x) for x in self.seldf['expdate'].unique()]


    @param.depends('seldf',watch=True)
    def exp_updater(self):
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
        self.dtabs=pn.Tabs(*self.tabbies,max_width=200)
        return self.dtabs

    def theplot(self):
        self.p=figure()
        self.p.plot_width=750
        self.p.plot_height=500
        return self.p

    def plot_view(self):
        gspec = pn.GridSpec(height=200,width=800,max_width=10,max_height=20,width_policy='max',height_policy='max')#max_width=50,max_height=20,width_policy='auto')#='stretch_width')#width=200)#,align='center')#sizing_mode='stretch_both', max_height=800,max_width=200,height_policy='max')
        gspec[1:2,0:1]=pn.Column(pn.Param(self.param.plotview_exps,
                        widgets={'plotview_exps':{'type':pn.widgets.Select,'name':'Active Experiment'}}),
                        sizing_mode='stretch_both')#,margin=(10,3,0,0))
        gspec[0:1,1:2]=pn.Column(pn.Param(self.param.plotview_enames,
                        widgets={'plotview_enames':{'type':pn.widgets.Select,'name':'Enzyme'}}),
                        sizing_mode='stretch_both')#,margin=(10,3,0,0))
        gspec[1:2,1:2]=pn.Column(pn.Param(self.param.plotview_snames,
                        widgets={'plotview_snames':{'type':pn.widgets.Select,'name':'Substrate'}}),
                        sizing_mode='stretch_both')#,margin=(10,3,0,0))

        return gspec
#        pn.Row(pn.Column(pn.Column(pn.Param(self.param.plotview_exps,
                        

#    plotview_exps=param.List()
#    plotview_enames=param.Selector()
#    plotview_snames=param.Selector()
#    plotview_xvariable=param.Selector()
#        ))))
#        gspec[0:1,0:1] = 
#        pn.Column(pn.Param(self.param.filter_enames,
#                        widgets={'filter_enames':{'type':pn.widgets.Checkbox,'name':'filter'}}),
#                        sizing_mode='stretch_both',margin=(10,3,0,0))
#
    def visualize(self):
        gspec = pn.GridSpec(height=500,width=800,max_width=10,max_height=20,width_policy='max',height_policy='max')#max_width=50,max_height=20,width_policy='auto')#='stretch_width')#width=200)#,align='center')#sizing_mode='stretch_both', max_height=800,max_width=200,height_policy='max')
        gspec[0:22,0:25]=self.theplot()
        gspec[0:22,25:30]=self.sel_view()
        gspec[0:22,30:35]=self.exp_view()
        gspec[22:30,0:35]=self.plot_view()#pltpanel.PlotViewPanel()
        return gspec
#        self.panel_layout=pn.Row(self.sel_view(),self.exp_view())
#        return self.panel_layout

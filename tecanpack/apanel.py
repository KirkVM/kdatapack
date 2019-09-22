import datetime
import param
import panel as pn
import pandas as pd
from tecanpack import readers
#plateset=readers.load_tecandata('allfiles.yml','kirk')#refresh_all=True)
#alldf=plateset.get_df()
alldf=pd.read_pickle('kirkdf.pkl')
class ActivityPanel2(param.Parameterized):
#    stuff=['All']+[str(x) for x in alldf['expdate'].unique()]
    exp_date=param.Selector(objects=[str(x) for x in alldf['expdate'].unique()])#,default='All')
    #exp_date=param.Selector(objects=['All']+[str(x) for x in alldf['expdate'].unique()],default='All')
    enzymes=param.Selector(objects=[str(x) for x in alldf['ename'].unique()])
    filter_exp_date=param.Boolean(False)
    filter_enzymes=param.Boolean(False)
#    junk=['a','b','c']
    pname_dict={'enzymes':'ename','exp_date':'expdate'}
    #def __init__(self,initdf,**params):
    #    #super(ActivityPanel2, self).__init__(**params)
    #    self.df=initdf.copy()
    ##expdate=param.Selector(objects=['things','stuff'])
    def sel_viewer(self):
#        return pn.GridSpec(pn.Column(
#                    pn.Param(self.param.filter_exp_date,
#                        widgets={'filter_exp_date':{'type':pn.widgets.Toggle,'name':'filter'}}),
#                    pn.Param(self.param.exp_date,
#                        widgets={'exp_date':{'type':pn.widgets.Select,'name':'Exp Date'}}),
#                    width=50),
#                    
#                    pn.Column(
#                    pn.Param(self.param.filter_enzymes,
#                        widgets={'filter_enzymes':{'type':pn.widgets.Toggle,'name':'filter'}}),
#                    pn.Param(self.param.enzymes,
##                        widgets={'enzymes':{'type':pn.widgets.Select,'name':'Enzyme'}})
#                    )
#
#
#
#                )
#
        #gspec = pn.GridSpec(sizing_mode='stretch_both', max_height=800,max_width=200,height_policy='max')
        #gspec = pn.GridSpec(sizing_mode='scale_both')#max_width=50,max_height=20,width_policy='auto')#='stretch_width')#width=200)#,align='center')#sizing_mode='stretch_both', max_height=800,max_width=200,height_policy='max')
        gspec = pn.GridSpec(height=500,width=250,max_width=10,max_height=20,width_policy='max',height_policy='max')#max_width=50,max_height=20,width_policy='auto')#='stretch_width')#width=200)#,align='center')#sizing_mode='stretch_both', max_height=800,max_width=200,height_policy='max')
#        gspec = pn.GridSpec(max_width=50,width_policy='max')#max_width=50,max_height=20,width_policy='auto')#='stretch_width')#width=200)#,align='center')#sizing_mode='stretch_both', max_height=800,max_width=200,height_policy='max')
        #gspec = pn.GridSpec(height_policy='fixed',width_policy='fit')#max_width=50,max_height=20,width_policy='auto')#='stretch_width')#width=200)#,align='center')#sizing_mode='stretch_both', max_height=800,max_width=200,height_policy='max')
        #gspec = pn.GridSpec(height=800,height_policy='min')

        gspec[0:1,0:1] = pn.Column(pn.Param(self.param.filter_exp_date,
                        widgets={'filter_exp_date':{'type':pn.widgets.Toggle,'name':'filter'}}),
                        sizing_mode='stretch_both',margin=(10,3,0,0))
        gspec[1:2,0:1] = pn.Column(pn.Param(self.param.filter_enzymes,
                        widgets={'filter_enzymes':{'type':pn.widgets.Toggle,'name':'filter'}}),
                        sizing_mode='stretch_both',margin=(10,3,0,0))
                        #sizing_mode='stretch_both')
        gspec[0:1,1:5]= pn.Column(pn.Param(self.param.exp_date,
                         widgets={'exp_date':{'type':pn.widgets.Select,'name':'Exp Date'}}),
                         sizing_mode='stretch_both')
        gspec[1:2,1:5]= pn.Column(pn.Param(self.param.enzymes,
                         widgets={'enzymes':{'type':pn.widgets.Select,'name':'Enzymes'}}),
                        sizing_mode='stretch_both')
        gspec[2:10,:5]=pn.Spacer(margin=0)
        return gspec


#        return pn.Row(
#                pn.Column(
#                    pn.Param(self.param.filter_exp_date,
#                        widgets={'filter_exp_date':{'type':pn.widgets.Toggle,'name':'filter'}}),
#                    pn.Param(self.param.filter_enzymes,
#                        widgets={'filter_enzymes':{'type':pn.widgets.Toggle,'name':'filter'}}),
#                    sizing_mode='stretch_both',width=20
#                ),
#                pn.Column(
#                    pn.Param(self.param.exp_date,
#                         widgets={'exp_date':{'type':pn.widgets.Select,'name':'Exp Date'}}),
#                    pn.Param(self.param.enzymes,
#                         widgets={'enzymes':{'type':pn.widgets.Select,'name':'Enzymes'}}),
#                )
#                )

#        gspec[:2,1]=pn.Spacer(margin=0,width=2)
#        gspec[2:8,:6]=pn.Spacer(margin=0)
#        .Spacer(background='#FF0000')
#        return gspec
#gspec[0, :3] = pn.Spacer(background='#FF0000')
#gspec[1:3, 0] = pn.Spacer(background='#0000FF')
#gspec[1:3, 1:3] = fig
#gspec[3:5, 0] = hv.Curve([1, 2, 3])
#gspec[3:5, 1] = 'https://upload.wikimedia.org/wikipedia/commons/4/47/PNG_transparency_demonstration_1.png'
#gspec[4:5, 2] = pn.Column(
#    pn.widgets.FloatSlider(),
#    pn.widgets.ColorPicker(),
#    pn.widgets.Toggle(name='Toggle Me!'))
#
#        pn.Param(self.param.exp_date,widgets={'exp_date':pn.widgets.RadioButtonGroup})
#        return pn.Column(self.param.exp_date)
#        return pn.Pane(pn.Column(pn.Param(self.param.exp_date,widgets={'exp_date':pn.widgets.RadioButtonGroup})))
#        return pn.Column(pn.Param(self.param.exp_date,widgets={'exp_date':pn.widgets.RadioButtonGroup}))
#pn.Param(ActivityPanel2.param, widgets={'exp_date': pn.widgets.RadioButtonGroup})#,'button_type':'success'}})
#        return pn.Column(pn.Row(self.param.filter_exp_date,self.param.exp_date),pn.Row(self.param.enzymes))
        #return pn.Column(self.exp_date)#,self.enzymes)
#        return pn.widgets.RadioBoxGroup(self.paraexp_date.#,self.enzymes)
#        return pn.Row(self.param.exp_date)
    @param.depends('enzymes','exp_date',watch=True)
    def exp_viewer(self):
#        expdates=[str(x) for x in self.df.expdate.unique()]
#        tvals=[(str(x),pn.Row(pn.widgets.StaticText(value=str(expdates[x])),width=300 ) ) for x in range(len(expdates))]
        tabtups=[]
        expdf=alldf.copy()
        if self.enzymes!='All':
            expdf=expdf[expdf['ename']==self.enzymes] 
        if self.exp_date!='All':
            dt_expdate=datetime.date(*[int(x) for x in self.exp_date.split('-')])
            expdf=expdf[expdf['expdate']==dt_expdate] 

#        expdf=alldf[(alldf['ename']==self.enzymes) & (alldf['expdate']==dt_expdate)]
        expdates=[x for x in expdf.expdate.unique()]

        for expnum,expdate in enumerate(expdates):
            ifields=[]
            for pname in ['enzymes','exp_date']:
                vals=[str(x) for x in expdf[expdf.expdate==expdate][self.pname_dict[pname]].unique()]
                if len(vals)==0:
                    rowvalstr='-'
                else:
                    rowvalstr=''.join(f'{x},' for x in vals[:-1])
                    rowvalstr+=f'{vals[-1]}'
                r_stext=pn.widgets.StaticText(value=f'{pname}: {rowvalstr}')
                ifields.append(r_stext)
            tabtups.append(( str(expnum),pn.Column(*ifields,width=200)  ))#r1_stext,r2_stext,width=300)   ))
        dtabs=pn.Tabs(*tabtups)
        return(dtabs)
#
pn.Param(ActivityPanel2.param, widgets={'exp_date': pn.widgets.RadioButtonGroup})#,'button_type':'success'}})
#    'enzymes': pn.widgets.DiscretePlayer})
#
#
#
#
#
#            for ifkey in info_fields.keys():
#        
#                vals=[str(x) for x in self.df[self.df.expdate==expdate][ifkey].unique()]
#                if len(vals)==0:
#                    rowvalstr='-'
#                else:
#                    rowvalstr=''.join(f'{x},' for x in vals[:-1])
#                    rowvalstr+=f'{vals[-1]}'
#                r_stext=pn.widgets.StaticText(value=f'{info_fields[ifkey]}: {rowvalstr}')
#                ifields.append(r_stext)
#            tabtups.append(( str(expnum),pn.Column(*ifields,width=200)  ))#r1_stext,r2_stext,width=300)   ))
#        self.dtabs=pn.Tabs(*tabtups)
# 
#        for expnum,expdate in enumerate(expdates):
#            ifields=[]
#            vals=[str(x) for x in expdf[expdf.expdate==expdate]['ename'].unique()]
#            if len(vals)==0:
#                rowvalstr='-'
#            else:
#                rowvalstr=''.join(f'{x},' for x in vals[:-1])
#                rowvalstr+=f'{vals[-1]}'
#                r_stext=pn.widgets.StaticText(value=f'enzymes: {rowvalstr}')
##        return r_stext
#                ifields.append(r_stext)
#            tabtups.append(( str(expnum),pn.Column(*ifields,width=200)  ))#r1_stext,r2_stext,width=300)   ))
#        dtabs=pn.Tabs(*tabtups)
#        return(dtabs)
#


#        return pn.Pane(dtabs)
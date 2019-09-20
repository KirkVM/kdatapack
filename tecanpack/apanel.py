import datetime
import param
import panel as pn
from tecanpack import readers
plateset=readers.load_tecandata('allfiles.yml','kirk')#refresh_all=True)
alldf=plateset.get_df()
class ActivityPanel2(param.Parameterized):
    stuff=['All']+[str(x) for x in alldf['expdate'].unique()]
    exp_date=param.Selector(objects=stuff)#,default='All')
    #exp_date=param.Selector(objects=['All']+[str(x) for x in alldf['expdate'].unique()],default='All')
    enzymes=param.Selector(objects=[str(x) for x in alldf['ename'].unique()])
#    junk=['a','b','c']
    pname_dict={'enzymes':'ename','exp_date':'expdate'}
    #def __init__(self,initdf,**params):
    #    #super(ActivityPanel2, self).__init__(**params)
    #    self.df=initdf.copy()
    ##expdate=param.Selector(objects=['things','stuff'])
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
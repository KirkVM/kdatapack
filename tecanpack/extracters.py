def getexperimentdf(df,**kwargs):
    for kwarg,value in kwargs.items():
        df=df.dropna(subset=[kwarg])
        if type(value)==str:
            df=df[df[kwarg].apply(str.lower)==value.lower()]
        else:
            df=df[df[kwarg]==value]
    return df

def getcontroldf(df,detection,**kwargs):
    df=df.dropna(subset=['standardname'])
    df=df[df.detection=='BCA']
    if kwargs is not None:
        df=getexperimentdf(df,**kwargs)
    return df
 

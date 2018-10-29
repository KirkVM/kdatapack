import random
import pandas as pd

def get_random_replicates(df,xvarname,replace=True,seed=None):
    sampledfs_=[]
    for x,gdf in df.groupby(xvarname):
        sampledf=gdf.sample(frac=1,replace=replace,random_state=seed)
        sampledfs_.append(sampledf)
    newdf=pd.concat(sampledfs_)
    return newdf

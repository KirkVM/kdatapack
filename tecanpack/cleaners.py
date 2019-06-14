import numpy as np
def cleandf(df,overs='nan'):
    cleandf=df.copy()
    if overs=='nan':
        cleandf=cleandf.replace(to_replace={'measurement':{'OVER':np.nan}})
    return cleandf

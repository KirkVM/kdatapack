from ghseqdb import seqdbutils #requires that biodb package is installed
def remove_similars(etetree,cutoff,ghdbpathstr=None,pwidpathstr=None,keepaccs=[],keeporgs=[]):
    #first delete iff one has no activity data and no gene data
    t=etetree.copy()
    assert (ghdbpathstr is not None), "no ghdb"
    assert (pwidpathstr is not None), "no pwidfile"
    pwdf=pd.read_csv(pwidpathstr,sep='\s+',index_col=0,skiprows=1,header=None)
    pwdf.columns=list(pwdf.index)
    pwdf.rename(lambda x:x.split('.')[0],inplace=True,axis='index')
    pwdf.rename(lambda x:x.split('.')[0],inplace=True,axis='columns')

    #print(list(pwdf.index)[0])
    conn=seqdbutils.gracefuldbopen(ghdbpathstr) 
    c=conn.cursor()
    for l in t.get_leaves():
        l.name=l.name.split('.')[0] #make sure all of these get fixed
        c.execute('''SELECT * FROM CAZYSEQDATA WHERE acc=(?)''',(l.name,)) 
        row=c.fetchone()
        if row is None:
            print(f'no db entry for {l.name}')
            continue
        if row['pdbids'] is not None:
            l.add_feature('pdbids',row['pdbids'])
        if row['ecs'] is not None:
            l.add_feature('ecs',row['ecs']) 
    conn.close()
    print('tree decoration complete')
    for rownum in range(len(pwdf.index)-1,-1,-1):
        if rownum%100==0:
            print(rownum,pwdf.shape[0])
        for colnum in range(rownum-1,-1,-1):
            
            if pwdf.iloc[rownum,colnum]>cutoff:
#                ival=pwdf.index[rownum]
                cval=pwdf.index[colnum]
                if 'pdbids' in (t&cval).features:
                    continue
                if 'ecs' in (t&cval).features:
                    continue
                if cval in keepaccs:
                    continue
                #now add something to check if it's within some distance of key acc?
                pwdf.drop(cval,axis='columns',inplace=True)
                pwdf.drop(cval,axis='index',inplace=True) 
            if rownum>=pwdf.shape[0]:
                break
    t.prune(list(pwdf.index))
    return pwdf,t
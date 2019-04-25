
def calc_treelength(t):
    tlen_traverse=0
    for n in t.traverse():
        tlen_traverse+=n.dist
    return tlen_traverse

#def calc_treecoverage(t,accs_):
#    '''calculates tree length covered by a list of leaf nodes'''
#    tlen=0
#    tlen_list=0
#    for n in t.traverse():
#        tlen+=n.dist
#        if n.is_leaf():
#            if n.name in accs_:
#                tlen_list+=n.dist
#        else:
#            curleafnames_=list(map(lambda x:x.name,n.get_leaves()))
#            if len(set(accs_).intersection(curleafnames_))>0:
#                tlen_list+=n.dist
#    return tlen_list/tlen

def calc_treecoverage(t,thenames):
    '''calculates tree length covered by a list of leaf node names'''
    tlen=0
    tlen_list=0
    for n in t.traverse():
        tlen+=n.dist
        if n.is_leaf():
            if n.name in thenames:
                tlen_list+=n.dist
        else:
            if len(  list(set(thenames).intersection(n.get_leaf_names()))  )>0:
                tlen_list+=n.dist
    if tlen==0:
        print('WARNING: tree length==0. Returning coverage as fraction of leaf nodes')
        return len(thenames)/len(t.get_leaf_names())
    return tlen_list/tlen


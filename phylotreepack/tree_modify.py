import os,yaml,argparse
import pandas as pd
import numpy as np
from ete3 import Tree

def clusterstats(t,pwdf):
    tlength=0.0
    for n in t.traverse():
        tlength+=n.dist
    lnames=t.get_leaf_names()
    numleaves=len(lnames)
    sub_pwdf=pwdf.loc[lnames,lnames]
    avgpwid=sub_pwdf.mean().mean()
    return numleaves,tlength,avgpwid



def cladecluster_bysimilarity(treefpath,pwiddf_fpath,outgroup,cluster_minpwid,
                     cluster_minsize):
    loneraccs_=[]
    print(treefpath)
    t=Tree(treefpath)
    t.set_outgroup(t&outgroup)
    print(len(t.get_leaf_names()))
    pwid_df=pd.read_pickle(pwiddf_fpath)
#    tovisit_.append(rnode)
    tovisit_=[t]
    while(len(tovisit_)>0):
        node=tovisit_.pop()
        numleaves,tlength,avgpwid=clusterstats(node,pwid_df)
        print(numleaves,tlength,avgpwid)
        if numleaves>cluster_minsize and avgpwid>cluster_minpwid:
            print('building a new leaf group')
            groupaccs_=node.get_leaf_names()
            print('group size= {}'.format(len(groupaccs_)))
            node.add_feature("accs",groupaccs_)
            nc=node.children[:]
            for x,c in enumerate(nc):#node.children:
                node.remove_child(c)
        else:
            tovisit_.extend(node.children)

    clusternodes_=[]
    orphanleaves_=[]
    totalaccs=0
    for node in t.traverse():
        if node.is_leaf():
            try:
                numaccs=len(set(node.accs))
                clusternodes_.append(node)
                totalaccs+=numaccs
            except:
                orphanleaves_.append(node)
    for orphanleaf in orphanleaves_:
        orphacc=orphanleaf.name
        bestcluster=None
        bestcluster_pwid=0.0
        for clusternode in clusternodes_:
            clusteraccs_=clusternode.accs
            sub_pwdf=pwid_df.loc[orphacc,clusteraccs_]
            avgpwid=sub_pwdf.mean()
            if avgpwid>bestcluster_pwid:
                bestcluster_pwid=avgpwid
                bestcluster=clusternode
        bestcluster.accs.append(orphacc)
        orphanleaf.delete()#detach()
#    for node in t.traverse():
#        if node.is_leaf():
#            try:
#                print(len(node.accs))
#            except:
#                node.detach()
#
    for node in t.traverse():
        #print(node.is_root(),len(node.children),node.is_leaf())
        try:
            print(len(node.accs))
            accstr='cluster|'+node.accs[0]
            for acc in node.accs[1:]:accstr+='|'+acc
            node.add_feature('name',accstr)
        except:
            pass
    newtreefpath=os.path.join(os.path.split(treefpath)[0],"newtree.nw")
    print(newtreefpath)
    t.write(outfile=newtreefpath,features=['name'])

def ete_cluster_bysize(t:Tree,cluster_maxsize:int=50,cluster_minsize:int=5,\
                            collapse_sisters:bool=True,cleanup_merge=True,outgroup:str=None):
    """takes ete3 tree, then groups sets of leaf nodes and truncates tree at common ancestor.
    merge_sister_orphans options can be useful just for visualization purposes if tree has many polytomies
    
    Arguments:
    t: ete3 tree object
    
    Keyword Arguments:
    cluster_maxsize: maximum size for each cluster (default=50)
    cluster_minsize: minimum size for each cluster (default=5)
    collapse_sisters: whether to merge leaves or groups<threshold into a single branch (default True)
    cleanup_merge: whether to perform final step that clusters all subtrees even if size<cluster_minsize
    outgroup: name of outgroup node (default=None)
    
    Returns: 
    ete tree file with added feature 'subtrees', 'cluster_numleaves', and 'accs' for collapsed nodes, 
    which correspond to a list of child nodes, number of leaves, and all (exapanded) leaf names in the cluster.
    Collapsed nodes names: 'm_<numsubtrees>_<numleaves>', 's_<numsubtrees>_<numleaves>', 'c_<numsubtrees>_<numleaves>'
    """
    t=t.copy()
    if outgroup is not None:
        t.set_outgroup(t&outgroup)

    #Breadth-First Tree Traversal, stop when no.leaves<cluster_maxsize
    orphans=[]
    tovisit_=[t]
    print(f'--starting number of leaves: {len(t.get_leaf_names())}--')
    cluster_merges=[]
    while(len(tovisit_)>0):
        node=tovisit_.pop()
        lnames=node.get_leaf_names()
        numleaves=len(lnames)
        if numleaves<cluster_maxsize: 
            groupaccs_=node.get_leaf_names()
            if len(groupaccs_)>cluster_minsize:
                node.add_feature('cluster_numleaves',len(node.get_leaf_names()))
                node.add_feature('subtrees',[nc.detach() for nc in node.get_children()])#
                node.name=f'm_{len(node.subtrees)}_{node.cluster_numleaves}'
                cluster_merges.append(len(groupaccs_))
            else:
                orphans.append(node)
        else:
            tovisit_.extend(node.children)
    print(f'cluster collapse sizes: {cluster_merges}')

    #merge sister orphans
    sister_merges=[]
    if collapse_sisters:
        while(len(orphans)>0):
            cur_orphan=orphans.pop()
            pnode=cur_orphan.up
            sisters=cur_orphan.get_sisters()
            sis_orphs=[]
            for orphos in range(len(orphans)-1,-1,-1):
                if orphans[orphos] in sisters:
                    sis_orphs.append(orphans.pop(orphos))
            if len(sis_orphs)>1:
                size_of_merged=len(cur_orphan.get_leaf_names())+np.sum([len(x.get_leaf_names()) for x in sis_orphs])
                if size_of_merged>cluster_minsize:
                    newnode=pnode.add_child(dist=0)
                    newnode.add_feature('cluster_numleaves',size_of_merged)
                    newnode.add_feature('subtrees',[cur_orphan.detach()])
                    for so in sis_orphs:
                        newnode.subtrees.append(so.detach())
                    newnode.name=f's_{len(newnode.subtrees)}_{newnode.cluster_numleaves}'
                    sister_merges.append(size_of_merged)
    print(f'sisters collapse sizes: {sister_merges}')
    #now cleanup_merge
    if cleanup_merge:
        visited=set()
        while len(visited)<len(t.get_leaves()):
            for n in set(t.get_leaves()).difference(visited):
                visited.add(n)
                nosubtrees='subtrees' not in n.features
                addn=None
                #if no subtrees, see how far we can climb
                while nosubtrees:
                    #if we've climbed once, this is a cluster
                    if len(n.get_descendants())>0:
                        addn=n
                    n=n.up
                    #continue only if no descendants have subtree
                    subtree_status=['subtrees' in x.features for x in n.traverse()]
                    nosubtrees=True not in subtree_status
                #clusterify the node, then add it and all sub-leaves to the visited set
                if addn is not None:
                    addn.add_feature('cluster_numleaves',len(addn.get_leaf_names()))
                    addn.add_feature('subtrees',[nc.detach() for nc in addn.get_children()])#
                    addn.name=f'c_{len(addn.subtrees)}_{addn.cluster_numleaves}'
                    visited=visited.union([x for x in addn.get_leaves()])
                    break #break to reset leaf candidates with updated visited set as filter

    #now at end add accs feature (a list of acc under each)
    for lnode in t.get_leaves():
        lnode.add_feature('accs',[])
        if 'subtrees' in lnode.features:
            for st in lnode.subtrees:
                lnode.accs.extend(st.get_leaf_names())#lnode.accs=[*y for y in [x.get_leaf_names() for x in lnode.subtrees]]
        else:
            lnode.accs.append(lnode.name)

    #final consistency check and readout
    num_leaves=0
    for lnode in t.get_leaves():
        if 'subtrees' in lnode.features:
            num_leaves+=np.sum([len(x.get_leaf_names()) for x in lnode.subtrees])
            #num_leaves+=lnode.cluster_numleaves#np.sum([len(x.get_leaf_names()) for x in lnode.subtrees])
        else:
            num_leaves+=1
    print(f'--total leaves at end: {num_leaves}--')
    return t
    #add merge with nephew orphan?

def expand_eteclustertree(ct:Tree,delete_cluster_names=True,delete_cluster_features=True):
    """expands cluster tree to original topology. 
    
    Collapsed sisters are inferred from a cluster node with dist=0 to parent.
    Returned tree will have some differences in node names from parent?
    Arguments:
    ct: the cluster tree

    Keyword Arguments:
    delete_cluster_names: whether to delete cluster names (default=True)
    delete_cluster_features: whether to delete cluster features accs,cluster_numleaves,subtrees (default=True)

    Returns:
    ete tree expanded so that leaves represent single enzymes
    """
    ct=ct.copy()
    for lnode in ct.get_leaves():
        if 'subtrees' in lnode.features:
            #special handling for collapsed sisters
            if lnode.dist==0: 
                print(f'sister node detected for {lnode.name}')
                for st in lnode.subtrees:
                    lnode.up.add_child(st)
                lnode.detach()
            else:
                for st in lnode.subtrees:
                    lnode.add_child(st)
                if delete_cluster_names:
                    lnode.name=''
                if delete_cluster_features:
                    lnode.del_feature('subtrees')
                    lnode.del_feature('cluster_numleaves')
                    lnode.del_feature('accs')
        #need to do a check... is this also proper handling collapsed sisters?
    return ct

def write_clustertree_tonewick(ct:Tree,otfpath:str='clustertree.nw'):
    """redefines node names property according to child leaf accesion codes and writes out as 
    newick tree
    
    Arguments:
    ctree: ete tree clustered using ete_clustertree_bysize, or ...? (contains 'subtrees' feature)
    Keyword Arguments:
    otfpath: path of output tree (default='clustertree.nw')

    Returns:
    int value number of nodes
    """
    ct=ct.copy()
    for lnode in ct.get_leaves():
        if 'subtrees' in lnode.features:
            lnode.name=''
            for st in lnode.subtrees:
                for acc in st.get_leaf_names():
                    lnode.name+=f'{acc}|'
        else:
            lnode.name+='|'
    ct.write(outfile=otfpath,features=['name'])
    print(f'clusterified newick tree written to {otfpath}')

def read_clustertree_fromnewick(treefpath:str):
    """reads in clustertree as defined in write_clustertree_fromnewick

    Arguments:
    treefpath: newick tree with leaf clusternames = including |-delimited accs

    Returns:
    ete3.Tree with same names and added feature value of accs = list of subtree accs
    """
    ctree=Tree(treefpath)
    for lnode in ctree.get_leaves():
        accnames=lnode.name
        lnode.add_feature('accs',[x for x in accnames.strip('|').split('|')])
    return ctree


"""
def calc_leafweights(t):
        seqwts_=[]
        accs_=[]
        for l in t.get_leaves():
                curseqname=l.name
                divisor=1.0
                seqwt=0.0
                node=l
                while not node.is_root():
                        seqwt+=t.get_distance(node,node.up)/divisor
                        node=node.up
                        divisor*=len(node.children)
                seqwts_.append(seqwt)
                accs_.append(l.name)
        #print len(accs_),np.sum(np.array(seqwts_))
        seqwtsA_=np.array(seqwts_)
        wtsum=np.sum(seqwtsA_)
        seqwtsA_*=(len(accs_)/wtsum)
        wtHT={}
        for a,w in zip(accs_,seqwtsA_):
                wtHT[a]=w
        return wtHT

def calc_treelength(t):
        #childnodes_=t.children[:] #here's the manual BFS implementation
        #tlen=0
        #for cn in childnodes_:
        #       childnodes_.extend(cn.children)
        #       tlen+=cn.dist
        tlen_traverse=0
        for n in t.traverse():
                tlen_traverse+=n.dist
        return tlen_traverse

def bottomup_prune(t,accs_):
    '''returns a tree pruned to contain only leaves and nodes
    with names indicated in accs_; also removes parent nodes that dont
    have any children in that list. deepcopied&pruned tree is returned'''
    t=deepcopy(t)
    for node in t.traverse(strategy="postorder"):
        leafnodes_=node.get_leaves()
        lnames_=list(map(lambda x:x.name,leafnodes_))
        lnames_.append(node.name)
        if len(set(lnames_).intersection(accs_))==0:
            node.detach()        
    return t

def topdown_prune(t,accs_):
    '''returns a tree pruned to contain only leaves and nodes
    with names indicated in accs_; also removes parent nodes that dont
    have any children in that list. deepcopied&pruned tree is returned'''
    t=deepcopy(t)
    for node in t.traverse(strategy="preorder"):
        leafnodes_=node.get_leaves()
        lnames_=list(map(lambda x:x.name,leafnodes_))
        if len(set(lnames_).intersection(accs_))==0:
            node.detach()        
    return t

def bottomup_prunev2(t,accs_):
    '''returns a tree pruned to contain only leaves and nodes
    with names indicated in accs_; also removes parent nodes that dont
    have any children in that list. deepcopied&pruned tree is returned'''
    t=deepcopy(t)
    ca=t.get_common_ancestor(accs_)
    for node in t.traverse(strategy="postorder"):
        leafnodes_=node.get_leaves()
        lnames_=list(map(lambda x:x.name,leafnodes_))
        lnames_.append(node.name)
        if len(set(lnames_).intersection(accs_))==0:
#            print(node.get_ancestors())
            if ca in node.get_ancestors():
                print('cant prune')
            #if ca not in node.get_ancestors():
            else:
                node.detach()        
#                print('pruned')
    return t

def calc_treecoverage(t,accs_):
        '''calculates tree length covered by a list of leaf nodes'''
        tlen=0
        tlen_list=0
        for n in t.traverse():
                tlen+=n.dist
                if n.is_leaf():
                        #print n.name,accs_
                        if n.name in accs_:
                                tlen_list+=n.dist
                else:
                        curleafnames_=list(map(lambda x:x.name,n.get_leaves()))
                        if len(set(accs_).intersection(curleafnames_))>0:
                                tlen_list+=n.dist
        return tlen_list/tlen

def calc_treecoverage_v2(t,accs_):
        tlen=calc_treelength(t)
        pt=bottomup_prune(t)
        tlen_list=calc_treecoverage(pt)
        return tlen_list/tlen
"""
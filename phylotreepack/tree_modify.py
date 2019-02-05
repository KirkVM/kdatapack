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

def ete_cluster_bysize(treefpath:str,cluster_maxsize:int=50,cluster_minsize:int=5,\
                            collapse_sisters:bool=True,cleanup_merge=True,outgroup:str=None):
    """reads in a newick tree file, groups sets of leaf nodes and truncates tree at common ancestor.
    merge_sister_orphans options can be useful just for visualization purposes if tree has many polytomies
    
    Arguments:
    treefpath: str representing path to newick tree file
    
    Keyword Arguments:
    cluster_maxsize: maximum size for each cluster (default=50)
    cluster_minsize: minimum size for each cluster (default=5)
    collapse_sisters: whether to merge leaves or groups<threshold into a single branch (default True)
    cleanup_merge: whether to perform final step that clusters all subtrees even if size<cluster_minsize
    outgroup: node (default=None)
    
    Returns: 
    ete tree file with added feature 'subtrees' and 'cluster_numleaves' for collapsed nodes, 
    which correspond to a list of child nodes and number of leaves in the cluster
    """
    print(f'reading newick file {treefpath}')
    t=Tree(treefpath)
    if outgroup is not None:
        t.set_outgroup(t&outgroup)

    #Breadth-First Tree Traversal, stop when no.leaves<cluster_maxsize
    clusters=[]
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
#                print(f'building a new leaf group: {len(groupaccs_)} lnodes')
                node.add_feature('cluster_numleaves',len(node.get_leaf_names()))
                node.add_feature('subtrees',[nc.detach() for nc in node.get_children()])#
                node.name=f'cool {node.cluster_numleaves} {node.get_distance(node.up)}'
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
                    newnode.name=f'mcool {np.sum([len(x) for x in newnode.subtrees])}'
                    sister_merges.append(size_of_merged)
    print(f'sisters collapse sizes: {sister_merges}')
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
                    addn.name=f'fresh {addn.cluster_numleaves} {addn.get_distance(addn.up)}'
                    visited=visited.union([x for x in addn.get_leaves()])
                    break #break to reset leaf candidates with updated visited set as filter

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

if __name__=="__main__":
    local_configfile=open('config.yml','r')
    defaultsHT=yaml.load(local_configfile)
    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--seqfldrpath',default=defaultsHT['seqfldr_basepath'],
                        help="path to sequence families base directory\n"
                        "current default=%(default)s")
    parser.add_argument('--ghfam',default=defaultsHT['ghfam'],
                        help="gh family (will look for this subfolder)\n"
                        "current default=%(default)s")
    parser.add_argument('--ghfamsubfldr',default=defaultsHT['mergefldr'],
                        help="subfolder within ghfam folder\n"
                        "current default=%(default)s")
    parser.add_argument('--treename',default=defaultsHT['filtertree'],
                        help="name of tree to assist in final selection\n"
                        "current default=%(default)s")
    parser.add_argument('--pwid_df',default=defaultsHT['pwid_df'],
                        help="name of (integer) pwid dataframe\n"
                        "current default=%(default)s")
    parser.add_argument('--outgroup',default=defaultsHT['tree_outgroup'],
                        help="accession code of tree outgroup\n"
                        "current default=%(default)s")
    parser.add_argument('--cluster_minpwid',type=int,default=20,
                        help="int value pwid threshold for cluster\n"
                        "current default=%(default)s")
    parser.add_argument('--cluster_minsize',type=int,default=50,
                        help="int value minimum size of cluster\n"
                        "current default=%(default)s")
    parser.add_argument('-f','--force',action='store_true',default=False,
                        help="overwrite one or both output files\n"
                        "default=%(default)s")

    args=parser.parse_args()
    treefpath=os.path.join(args.seqfldrpath,args.ghfam,args.ghfamsubfldr,
                        args.treename)
    pwiddf_fpath=os.path.join(args.seqfldrpath,args.ghfam,args.ghfamsubfldr,
                        args.pwid_df)
    make_cluster_tree(treefpath,pwiddf_fpath,args.outgroup,
                            args.cluster_minpwid,args.cluster_minsize)
                            #force=args.force)

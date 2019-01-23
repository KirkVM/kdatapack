import os,yaml,argparse
import pandas as pd
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

def ete_cladecluster_bysize(treefpath:str,cluster_maxsize:int,merge_sister_orphans:bool=False,orphan_size:int=5,\
                              outgroup:str=None):
    """reads in a newick tree file, groups sets of leaf nodes and truncates tree at common ancestor.
    merge_sister_orphans options can be useful just for visualization purposes if tree has many polytomies
    
    Arguments:
    treefpath: str representing path to newick tree file
    cluster_maxsize: maximum size for each cluster
    
    Keyword Arguments:
    merge_orphans: whether to merge single nodes into nearest cluster (default True)
    orphan_size: minimum cluster size that would be merged (default=2)
    
    Returns: 
    ete tree file with feature 'name' for leaf nodes, set to "cluster|acc1|acc2|..."
    """
    print(f'reading newick file {treefpath}')
    t=Tree(treefpath)
    if outgroup is not None:
        t.set_outgroup(t&outgroup)
    #print(len(t.get_leaf_names()))
    #pwid_df=pd.read_pickle(pwiddf_fpath)
#    tovisit_.append(rnode)

    #Breadth-First Tree Traversal, stop when no.leaves<cluster_maxsize
    clusters=[]
    orphans=[]
    tovisit_=[t]
    print(len(t.get_leaf_names()))
    newnodes=0
    while(len(tovisit_)>0):
        node=tovisit_.pop()
        lnames=node.get_leaf_names()
        numleaves=len(lnames)
        if numleaves<cluster_maxsize: 
            groupaccs_=node.get_leaf_names()
#            print(f'building a new leaf group: {len(groupaccs_)} lnodes')
            node.add_feature("accs",groupaccs_)
#            node.add_feature('subtree',node.detach)
            node.name=f'cool {len(groupaccs_)}, {node.get_distance(node.up)}'
            newnodes+=len(groupaccs_) 
            if len(groupaccs_)>orphan_size:
                clusters.append(node)
            else:
                orphans.append(node)
            nc=node.children[:]
#            for c in nc:#node.children:
#                c.detach()
#            node.name='cool'
        else:
            tovisit_.extend(node.children)
#    for orph in orphans:
#        print(orph)
    #for n in t.traverse():
    print(newnodes)
#    print(t)
    #lengthier part is merging orphan clusters into larger
    if merge_sister_orphans:
    #first sort orphans so smallest group is last-
    orphans.sort(key=lambda x:len(x.accs),reverse=True)#)#attrgetter('accs'))
    #try to merge sister orphans first
    while(len(orphans)>0 ):
        cur_orphan=orphans.pop()
        groupaccs_=cur_orphan.accs
#        print(len(cur_orphan.get_sisters()))
        sisters=[n for n in cur_orphan.get_sisters() if (n.is_leaf() and n in orphans)]
        sisters.sort(key=lambda x:len(x.accs))
        if len(sisters)>0:
            mergenode=sisters[0] 
            mergenode.accs.extend(cur_orphan.accs)
            mergenode.name=f'mcool {len(mergenode.accs)}, {node.get_distance(node.up)}'
            cur_orphan.detach()
        for orphnum,orph in enumerate(orphans):
            if len(orph.accs)>orphan_size:
                clusters.append(orphans.pop(orphnum))
        orphans.sort(key=lambda x:len(x.accs),reverse=True)#)#attrgetter('accs'))
#        print(len(sisters))
    t.ladderize()
    print(t)
    return
    #add merge with nephew orphan?
    print(len(t.get_leaf_names()))



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

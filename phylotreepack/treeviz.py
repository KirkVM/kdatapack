import matplotlib.pyplot as plt
import numpy as np

def depth_sort(trees):
    dtzips=[]
    for t in trees:
        leaves=t.get_leaves()
        max_depth=max([len(l.get_ancestors()) for l in leaves])
        dtzips.append([max_depth,t])
    dtzips.sort(key=lambda x:x[0])#,reverse=True)
    sorted_trees=[x[1] for x in dtzips]#.sort(key=lambda x:x[0])]
    return sorted_trees

def plot_branch(tree,xcoord,ycoord):
    if tree.is_leaf():
        plt.plot([xcoord,xcoord+tree.dist],[ycoord,ycoord],'g-')
        return ycoord,ycoord-0.1
    else:
        xcoord=xcoord+tree.dist
        branches=tree.children
        sorted_branches=depth_sort(branches)
    ycoord_ascends=[]
    for x,branch in enumerate(sorted_branches):
        #ycoord0s.append(ycoord)
        ycoord_ascend,ycoord_descend=plot_branch(branch,xcoord,ycoord)
        ycoord=ycoord_descend
        ycoord_ascends.append(ycoord_ascend)
    plt.plot([xcoord for x in range(len(ycoord_ascends))],ycoord_ascends,lw=5.0)
    ycoord_ascend=np.mean(ycoord_ascends)
    plt.plot([xcoord-tree.dist,xcoord],[ycoord_ascend,ycoord_ascend],lw=3.0)
    return ycoord_ascend,ycoord_descend

def mpl_tree_render(tree):
    plt.figure(figsize=(15,45))
    plot_branch(tree,0,0)
#plt.axis((-0.3,40,-9,0.3))
#plt.xticks(ticks=[0,1,2])
#plt.axes.tick_params(axis='x',)
#plot_branch(mytree,0,0)
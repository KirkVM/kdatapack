import matplotlib.pyplot as plt
import numpy as np
from typing import Iterable
from ete3 import Tree
from matplotlib import markers
from matplotlib.path import Path

def align_marker(marker, halign='center', valign='middle',):
    """
    create markers with specified alignment.

    Parameters
    ----------

    marker : a valid marker specification.
      See mpl.markers

    halign : string, float {'left', 'center', 'right'}
      Specifies the horizontal alignment of the marker. *float* values
      specify the alignment in units of the markersize/2 (0 is 'center',
      -1 is 'right', 1 is 'left').

    valign : string, float {'top', 'middle', 'bottom'}
      Specifies the vertical alignment of the marker. *float* values
      specify the alignment in units of the markersize/2 (0 is 'middle',
      -1 is 'top', 1 is 'bottom').

    Returns
    -------

    marker_array : numpy.ndarray
      A Nx2 array that specifies the marker path relative to the
      plot target point at (0, 0).

    Notes
    -----
    The mark_array can be passed directly to ax.plot and ax.scatter, e.g.::

        ax.plot(1, 1, marker=align_marker('>', 'left'))

    -----
    swiped from StackOverflow-
    https://stackoverflow.com/questions/26686722/align-matplotlib-scatter-marker-left-and-or-right
    """

    if isinstance(halign, (str)):
        halign = {'right': -1.,
                  'middle': 0.,
                  'center': 0.,
                  'left': 1.,
                  }[halign]

    if isinstance(valign, (str)):
        valign = {'top': -1.,
                  'middle': 0.,
                  'center': 0.,
                  'bottom': 1.,
                  }[valign]

    # Define the base marker
    bm = markers.MarkerStyle(marker)

    # Get the marker path and apply the marker transform to get the
    # actual marker vertices (they should all be in a unit-square
    # centered at (0, 0))
    m_arr = bm.get_path().transformed(bm.get_transform()).vertices

    # Shift the marker vertices for the specified alignment.
    m_arr[:, 0] += halign / 2
    m_arr[:, 1] += valign / 2

    return Path(m_arr, bm.get_path().codes)


def depth_sort(trees):
    dtzips=[]
    for t in trees:
        leaves=t.get_leaves()
        max_depth=max([len(l.get_ancestors()) for l in leaves])
        max_dist=max([t.get_distance(l) for l in leaves])
        dtzips.append([max_depth,max_dist,t])
#    print(dtzips)
    dtzips.sort(key=lambda x:x[1])#,reverse=True)
    sorted_trees=[x[2] for x in dtzips]#.sort(key=lambda x:x[0])]
    return sorted_trees

def mpl_plot_branch(tree,xcoord,ycoord,sepsize=0.1):
    if tree.is_leaf():
        plt.plot([xcoord,xcoord+tree.dist],[ycoord,ycoord],lw=3.0,color='black')#g-')
#        plt.plot([xcoord+tree.dist],[ycoord],'<',ms=50*tree.cluster_relsize)
        plt.plot([xcoord+tree.dist],[ycoord],marker=align_marker('<', halign='left'),
           clip_on=False, color='k',ms=50*tree.cluster_relsize)#, transform=plt.gca().get_xaxis_transform())
        return ycoord,ycoord-sepsize
    else:
        xcoord=xcoord+tree.dist
        branches=tree.children
        sorted_branches=depth_sort(branches)
    ycoord_ascends=[]
    for branch in sorted_branches:
        ycoord_ascend,ycoord_descend=mpl_plot_branch(branch,xcoord,ycoord)
        ycoord=ycoord_descend
        ycoord_ascends.append(ycoord_ascend)
    plt.plot([xcoord for x in range(len(ycoord_ascends))],ycoord_ascends,lw=3.0,color='black')
    ycoord_ascend=np.mean(ycoord_ascends)
    plt.plot([xcoord-tree.dist,xcoord],[ycoord_ascend,ycoord_ascend],lw=3.0,color='black')
    return ycoord_ascend,ycoord_descend

class EteMplTree:
    def __init__(self,tree:Tree):
        self.tree=tree
        self.orientation='left'
        self.dashed=True
        self.cluster_property='accs'
        self.cluster_viz='triangle'
        #self.figsize=(15,45)
        self.set_cluster_size()
    def render(self,figsize:Iterable=None):
        if figsize is not None:
            plt.figure(figsize=figsize)
        else:
            plt.figure(figsize=self.figsize)
#        else:
#            plt.figure()
        mpl_plot_branch(self.tree,0,0)
    def set_cluster_size(self):
        if self.cluster_property=='accs':
            cluster_sizes=[len(l.accs) for l in self.tree.get_leaves()]
            for l in self.tree.get_leaves():
                l.add_feature('cluster_relsize',len(l.accs)/max(cluster_sizes))

#def mpl_tree_render(tree:Tree,orientation:str='left',dashed:bool=):
#    plt.figure(figsize=(15,45))
#    plot_branch(tree,0,0)
#plt.axis((-0.3,40,-9,0.3))
#plt.xticks(ticks=[0,1,2])
#plt.axes.tick_params(axis='x',)
#plot_branch(mytree,0,0)
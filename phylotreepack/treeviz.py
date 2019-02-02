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

def mpl_plot_branch(branch,xcoord,ycoord,sepsize=0.1):
    if branch.is_leaf():
        #branch.add_feature('stem_coords',[[xcoord,xcoord+branch.dist],[ycoord,ycoord]])
        branch.add_feature('stem_coord_offset',np.array([ycoord,ycoord]))#[[xcoord,xcoord+branch.dist],[ycoord,ycoord]])
        branch.add_feature('stem_coord_span',np.array([xcoord,xcoord+branch.dist]))
        return ycoord,ycoord-sepsize
    else:
        xcoord=xcoord+branch.dist
        sub_branches=branch.children
        sorted_sub_branches=depth_sort(sub_branches)
    ycoord_ascends=[]
    for sub_branch in sorted_sub_branches:
        ycoord_ascend,ycoord_descend=mpl_plot_branch(sub_branch,xcoord,ycoord)
        ycoord=ycoord_descend
        ycoord_ascends.append(ycoord_ascend)
    #branch.add_feature('base_coords',[[xcoord for x in range(len(ycoord_ascends))],ycoord_ascends])
    branch.add_feature('base_coord_offset', np.array([xcoord for x in range(len(ycoord_ascends))]))
    branch.add_feature('base_coord_span',np.array(ycoord_ascends))#[[xcoord for x in range(len(ycoord_ascends))],ycoord_ascends])
    ycoord_ascend=np.mean(ycoord_ascends)
    #branch.add_feature('stem_coords',[[xcoord-branch.dist,xcoord],[ycoord_ascend,ycoord_ascend]])
    branch.add_feature('stem_coord_offset',np.array([ycoord_ascend,ycoord_ascend]))
    branch.add_feature('stem_coord_span',np.array([xcoord-branch.dist,xcoord]))#[[xcoord-branch.dist,xcoord],[ycoord_ascend,ycoord_ascend]])
    return ycoord_ascend,ycoord_descend

class EteMplTree:
    def __init__(self,tree:Tree):
        self.tree=tree.copy()
        self.orientation='left'
        self.dashed=True
        self.cluster_property='accs'
        self.cluster_viz='triangle'
        self.figsize=(15,45)
        self.set_cluster_size()
        self.plot_coords=[[np.inf,np.inf],[-np.inf,-np.inf]]
#        self.leaf_node_edges=[]
#        self.leaf_node_dots=[]
#        self.leaf_coords=[]
#        self.leaf_begin=None
#        self.leaf_end=None
        self.scale=1.0
    def render(self,figsize:Iterable=None,orientation='left'):
#        if figsize is not None:
#            plt.figure(figsize=figsize)
#        else:
#            plt.figure(figsize=self.figsize)
#        else:
#            plt.figure()
        plt.figure(figsize=(15,45))
        mpl_plot_branch(self.tree,0,0)
        self.get_plot_coords()
        print(self.plot_coords)
        self.plot_it()
    def set_cluster_size(self):
        if self.cluster_property=='accs':
            cluster_sizes=[len(l.accs) for l in self.tree.get_leaves()]
            for l in self.tree.get_leaves():
                l.add_feature('cluster_relsize',len(l.accs)/max(cluster_sizes))
    def get_plot_coords(self):
        for node in self.tree.traverse():
            stem_coords=[[np.nan],[np.nan]]
            base_coords=[[np.nan],[np.nan]]
            if self.orientation in ['left','right']:
                stem_coords[0]= self.scale*node.stem_coord_span
                stem_coords[1]= self.scale*node.stem_coord_offset
            else:
                stem_coords[1]= self.scale*node.stem_coord_span
                stem_coords[0]= self.scale*node.stem_coord_offset
            node.add_feature('stem_plot_coords',stem_coords)
            if node.is_leaf() is False:
                if self.orientation in ['left','right']:
                    base_coords[0]= self.scale*node.base_coord_offset
                    base_coords[1]= self.scale*node.base_coord_span
                else:
                    base_coords[1]= self.scale*node.base_coord_offset
                    base_coords[0]= self.scale*node.base_coord_span
                node.add_feature('base_plot_coords',base_coords)
#            print(stem_coords[0])
            self.plot_coords[0]=[min(*stem_coords[0],*base_coords[0],self.plot_coords[0][0]),\
                                 min(*stem_coords[1],*base_coords[1],self.plot_coords[0][1])]
            self.plot_coords[1]=[max(*stem_coords[0],*base_coords[0],self.plot_coords[1][0]),\
                                max(*stem_coords[1],*base_coords[1],self.plot_coords[1][1])]
    def plot_it(self):
        for node in self.tree.traverse():
            plt.plot(node.stem_plot_coords[0],node.stem_plot_coords[1],lw=3.0,color='black')
            if node.is_leaf() is False:
                plt.plot(node.base_plot_coords[0],node.base_plot_coords[1],lw=3.0,color='black')
            if node.is_leaf() and self.cluster_property=='accs':
                plt.plot(node.stem_plot_coords[0][-1],node.stem_plot_coords[1][-1],marker=align_marker('<', halign='left'),\
                 clip_on=False, color='k',ms=50*node.cluster_relsize)#, transform=plt.gca().get_xaxis_transform())


#        if self.dashed:
#            for lcoord in self.leaf_coords:
#                plt.plot('')
          


                #if 'stem' in node.features:
                #    plt.plot(node.stem[0],node.stem[1],lw=3.0,color='black')

#        plt.plot([xcoord,xcoord+tree.dist],[ycoord,ycoord],lw=3.0,color='black')#g-')
        
 
#def mpl_tree_render(tree:Tree,orientation:str='left',dashed:bool=):
#    plt.figure(figsize=(15,45))
#    plot_branch(tree,0,0)
#plt.axis((-0.3,40,-9,0.3))
#plt.xticks(ticks=[0,1,2])
#plt.axes.tick_params(axis='x',)
#plot_branch(mytree,0,0)
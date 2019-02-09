import matplotlib.pyplot as plt
import numpy as np
from typing import Iterable
from ete3 import Tree
from matplotlib import markers
from matplotlib.path import Path
from operator import attrgetter

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
    """sorts trees according to ranked criteria: 1.#nodes to deepest leaf, 2.dist. to deepest leaf

    Arguments:
        [ete trees] 
    Returns:
        [sorted ete trees], shortest first"""
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

def mpl_plot_branch(branch,xcoord,ycoord,sepsize):
    """Function called by EteMplTree to determine how to plot tree. 
    Recursively calls, determining 'stem_coord_offset' for leaves (equivalent to y-value 
    in a standard left-oriented graph) and 'stem_coord_span' (x-offset and length in a left-oriented graph)
    for leaves in descent into call stack. On return, intermediate node values are set:
    'stem_coord_offset' and 'stem_coord_span' as well as 'base_coord_offset' and 'base_coord_span', with
    latter 2 referring to x-offset and y-offset/length for lines connecting branches in left-oriented graph

    Arguments:
        branch: the ete3 tree
        xcoord: current xcoord value (in left-oriented graph) in descent into call stack (float,first call xcoord should=0)
        ycoord: current ycoord value (in left-orientd graph) in descent into call stack (float,first call ycoord should=0)
        sepsize: coordinate step size to take between nodes (float, default=0.1)
    Returns:
        ete3 tree with feature values ('stem_coord_offset',etc) set to enable plotting a tree
    """
    if branch.is_leaf():
        branch.add_feature('stem_coord_offset',np.array([ycoord,ycoord]))#[[xcoord,xcoord+branch.dist],[ycoord,ycoord]])
        branch.add_feature('stem_coord_span',np.array([xcoord,xcoord+branch.dist]))
        return ycoord,ycoord-sepsize
    else:
        xcoord=xcoord+branch.dist
        sub_branches=branch.children
        sorted_sub_branches=depth_sort(sub_branches)
    ycoord_ascends=[]
    for sub_branch in sorted_sub_branches:
        ycoord_ascend,ycoord_descend=mpl_plot_branch(sub_branch,xcoord,ycoord,sepsize)
        ycoord=ycoord_descend
        ycoord_ascends.append(ycoord_ascend)
    branch.add_feature('base_coord_offset', np.array([xcoord for x in range(len(ycoord_ascends))]))
    branch.add_feature('base_coord_span',np.array(ycoord_ascends))#[[xcoord for x in range(len(ycoord_ascends))],ycoord_ascends])
    ycoord_ascend=np.mean(ycoord_ascends)
    branch.add_feature('stem_coord_offset',np.array([ycoord_ascend,ycoord_ascend]))
    branch.add_feature('stem_coord_span',np.array([xcoord-branch.dist,xcoord]))#[[xcoord-branch.dist,xcoord],[ycoord_ascend,ycoord_ascend]])
    return ycoord_ascend,ycoord_descend
from matplotlib.axes import Axes
from matplotlib.transforms import IdentityTransform
class EteMplTree:
    """Class for plotting an ete3 tree using matplotlib rendering

    Attributes:
        orientation ('left','right','top','bottom'): tree orientation (str, default 'left')
        dashed_leaves: whether to add dashes from leaves to edge of graph (default True)
        cluster_feature: ete node feature to use to define cluster size (str, default 'accs')
        cluster_viz: symbol to use as cluster indicator (default 'triangle')---not yet impld
        ordered_leaves: list of leaves from graph start (in ladderized form) to end
                        (calculated through calling .render() method)
        plot_coords: [xmin,xmax],[ymin,ymax] values for graph
                        (calculated through calling .render() method)
    """
    def __init__(self,tree:Tree):
        """ constructor for EteMplTree.
        Calls: self.cluster_size() to add feature .cluster_relsize to each leaf node

        Arguments:
            tree: an ete3 tree instance
            [Also currently makes some assumed settings that are configurable public properties-
            -orientation,cluster_feature,cluster_viz,scale]
        """

        self.tree=tree.copy()
        self.orientation='left'
        self.scale=1.0
        self.dashed_leaves=True
        self.cluster_viz='triangle'
        self.cviz_symboldict={'left':'<','right':'>','top':'^','bottom':'v'}
        self.cviz_hadict={'left':'left','right':'right','top':'middle','bottom':'middle'}
        self.cviz_vadict={'left':'middle','right':'middle','top':'top','bottom':'bottom'}
        self.initial_leafspacing=0.1
#        self.figsize=(15,45)
        
        self.cluster_feature='accs'
        self.set_cluster_size()

        self.plot_coords=[[np.inf,-np.inf],[np.inf,-np.inf]]
        self.ordered_leaves=None
        self.decorated_plot_coords=None

#        self.
    def render(self,orientation=None,scale=None,figsize:Iterable=(20,20),initial_leafspacing:float=None,ax=None):
        """ method to generate figure.
        Calls: (1)mpl_plot_branch to set node attributes,
               (2)self.get_plot_coords() to set plot_coords attribute
               (3)self.plot_it() to make the actual plot

        Keyword arguments:
            orientation: 'left','right','bottom', or 'top' (default is to use object property .orientation)
            scale: float setting figure scaling (default is to use object property .scale)
            figsize: size-2 tuple with hxw (default= (20,20))

        """
        if initial_leafspacing is None:initial_leafspacing=self.initial_leafspacing
        if orientation is None: orientation=self.orientation
        if scale is None: scale=self.scale
        leafspacing=initial_leafspacing*scale

        mpl_plot_branch(self.tree,0,0,initial_leafspacing) #builds tree into ete3 tree
        self.set_plot_coords(orientation,scale) #adjusts these using scale/orientation
        self.ordered_leaves=sorted(self.tree.get_leaves(),key=lambda x:x.stem_coord_offset[0],reverse=True)

        if ax is None:
            plt.figure(figsize=figsize)
            ax=plt.gca()

        self.plot_it(ax,orientation,leafspacing) #self.decorated_plot_coords gets set in here
        #return silly
        self.trim_plotspace(ax,orientation,leafspacing) #and used here
    def set_cluster_size(self):
        """ 
        Method to add .cluster_relsize feature to each node in the .tree
        Currently only supports if .cluster_feature=='accs'
        """
 
        if self.cluster_feature=='accs':
            cluster_sizes=[len(l.accs) for l in self.tree.get_leaves()]
            for l in self.tree.get_leaves():
                l.add_feature('cluster_relsize',np.log(len(l.accs)+1)/np.log(max(cluster_sizes)+1))
    def set_plot_coords(self,orientation,scale):
        """ 
        Method to calculate plot coordinates. Calculates then adds feature 'base_plot_coordinates' 
        (internal nodes only) and 'stem_plot_coordinates' (all nodes) to all ete3 tree nodes for plotting. 
        Also generates overall .plot_coordinates

        Arguments:
            orientation: 'left','right','bottom' or 'top' (default is set in constructor (probably 'left'))
            scale: float (default is set in constructor (probably 1.0))
        """
        self.plot_coords=[[np.inf,-np.inf],[np.inf,-np.inf]]
        for node in self.tree.traverse():
            stem_coords=[[np.nan],[np.nan]]
            base_coords=[[np.nan],[np.nan]]
            if orientation in ['left','right']:
                stem_coords[0]= scale*node.stem_coord_span
                stem_coords[1]= scale*node.stem_coord_offset
            elif orientation in ['bottom','top']:
                stem_coords[1]= scale*node.stem_coord_span
                stem_coords[0]= scale*node.stem_coord_offset
            if orientation in ['right','top']:
                stem_coords[0]*=-1
                stem_coords[1]*=-1
            node.add_feature('stem_plot_coords',stem_coords)
            if node.is_leaf() is False:
                if orientation in ['left','right']:
                    base_coords[0]= scale*node.base_coord_offset
                    base_coords[1]= scale*node.base_coord_span
                elif orientation in ['bottom','top']:
                    base_coords[1]= scale*node.base_coord_offset
                    base_coords[0]= scale*node.base_coord_span
                if orientation in ['right','top']:
                    base_coords[0]*=-1
                    base_coords[1]*=-1
                node.add_feature('base_plot_coords',base_coords)
    def plot_it(self,ax,orientation,leafspacing):
        """ method to generate figure. 
        Currently called from within render() method as that first calls necessary get_plot_coordinates() method
        """
        for node in self.tree.traverse():
            ax.plot(node.stem_plot_coords[0],node.stem_plot_coords[1],lw=3.0,color='black')
            if node.is_leaf() is False:
                ax.plot(node.base_plot_coords[0],node.base_plot_coords[1],lw=3.0,color='black')
        self.plot_coords=[[ax.dataLim.x0,ax.dataLim.x1],[ax.dataLim.y0,ax.dataLim.y1]]
        self.decorated_plot_coords=self.plot_coords[:]
        print(self.plot_coords)
#        Axes.convert_xunits(0.1)
        if self.cluster_feature=='accs':
            dummy=ax.transData.transform((leafspacing,leafspacing))#leafspacing
            if orientation in ['left','right']:
                fig_leafspacing=dummy[1]
            elif orientation in ['bottom','top']:
                fig_leafspacing=dummy[0]#-dummy[1]
#                dummy=ax.transData.transform((leafspacing,0))#leafspacing
            for node in self.tree.get_leaves():
                cool=ax.plot(node.stem_plot_coords[0][-1],node.stem_plot_coords[1][-1],\
                    marker=align_marker(self.cviz_symboldict[orientation],halign=self.cviz_hadict[orientation],\
                    valign=self.cviz_vadict[orientation]),\
                     clip_on=False, color='k',ms=np.sqrt(fig_leafspacing)*node.cluster_relsize)
                c=cool[0].get_tightbbox(ax)
                inv = ax.transData.inverted()
                silly=inv.transform(c.corners())
                symbol_xvals=[x[0] for x in silly]
                symbol_yvals=[x[1] for x in silly]
                self.decorated_plot_coords[0]=[min(self.decorated_plot_coords[0][0],*symbol_xvals),
                                               max(self.decorated_plot_coords[0][1],*symbol_xvals)]
                self.decorated_plot_coords[1]=[min(self.decorated_plot_coords[1][0],*symbol_yvals),
                                               max(self.decorated_plot_coords[1][1],*symbol_yvals)]
        print(self.decorated_plot_coords)
        if self.dashed_leaves:
            for lnode in self.tree.get_leaves():
                if orientation=='left':
                    xs=[lnode.stem_plot_coords[0][-1],self.decorated_plot_coords[0][1]]
                    ys=[lnode.stem_plot_coords[1][-1],lnode.stem_plot_coords[1][-1]]
                if orientation=='right':
                    xs=[lnode.stem_plot_coords[0][-1],self.decorated_plot_coords[0][0]]
                    ys=[lnode.stem_plot_coords[1][-1],lnode.stem_plot_coords[1][-1]]
                if orientation=='bottom':
                    xs=[lnode.stem_plot_coords[0][-1],lnode.stem_plot_coords[0][-1]]
                    ys=[lnode.stem_plot_coords[1][-1],self.decorated_plot_coords[1][1]]
                if orientation=='top':
                    xs=[lnode.stem_plot_coords[0][-1],lnode.stem_plot_coords[0][-1]]
                    ys=[lnode.stem_plot_coords[1][-1],self.decorated_plot_coords[1][0]]
                ax.plot(xs,ys,lw=1,ls='--',zorder=0,color='gray')
            plt.plot()
        #ax.set_xticks([])
        ax.set_yticks([])

    def trim_plotspace(self,ax,orientation,leafspacing):
        #this can be improved! should adjust scale using cluster glyph locations...
        cur_mins=ax.transData.transform([x[0] for x in self.decorated_plot_coords])
        cur_maxes=ax.transData.transform([x[1] for x in self.decorated_plot_coords])
        print('.',cur_mins,cur_maxes)
        inv=ax.transData.inverted()
#        cur_mins[0]-=2
#        cur_maxes[0]+=2
#        cur_mins=inv.transform([x for x in cur_mins])
#        cur_maxes=inv.transform([x for x in cur_maxes])
#        print(self.decorated_plot_coords)
#        print(cur_mins)
        if orientation in ['left','right']:
#            ax.set_xlim(self.plot_coords[0][0]- \
#                    0.01*(self.plot_coords[0][1]-self.plot_coords[0][0]),\
#                    self.plot_coords[0][1]+ \
#                    0.01*(self.plot_coords[0][1]-self.plot_coords[0][0]))
            #ax.set_xlim(self.decorated_plot_coords[0],ax.IdentityTransform)#,ax.transData)#self.decorated_plot_coords[0][0],self.decorated_plot_coords[0][1]
            cur_mins[0]-=2
            cur_maxes[0]+=2
            cur_mins=inv.transform([x for x in cur_mins])
            cur_maxes=inv.transform([x for x in cur_maxes])
            ax.set_xlim(cur_mins[0],cur_maxes[0])#,transform=IdentityTransform)#,ax.transData)#self.decorated_plot_coords[0][0],self.decorated_plot_coords[0][1]
#                    0.01*(self.plot_coords[0][1]-self.plot_coords[0][0]),\
#                    self.plot_coords[0][1]+ \
#                    0.01*(self.plot_coords[0][1]-self.plot_coords[0][0]))
            ax.set_ylim(self.plot_coords[1][0]-0.5*leafspacing,self.plot_coords[1][1]+0.5*leafspacing)
        elif orientation in ['top','bottom']:
            #ax.set_ylim(self.plot_coords[1][0]- \
            #        0.01*(self.plot_coords[1][1]-self.plot_coords[1][0]),\
            #        self.plot_coords[1][1]+ \
            #        0.01*(self.plot_coords[1][1]-self.plot_coords[1][0]))
            cur_mins[1]-=2
            cur_maxes[1]+=2
            cur_mins=inv.transform([x for x in cur_mins])
            cur_maxes=inv.transform([x for x in cur_maxes])
            ax.set_ylim(cur_mins[1],cur_maxes[1])#,transform=IdentityTransform)#,ax.transData)#self.decorated_plot_coords[0][0],self.decorated_plot_coords[0][1]
            ax.set_xlim(self.plot_coords[0][0]-0.5*leafspacing,self.plot_coords[0][1]+0.5*leafspacing)
                
            

#plt.gca().set_xlim(0.0,3.3)
#plt.gca().spines['left'].set_visible(False)
         
#def mpl_tree_render(tree:Tree,orientation:str='left',dashed:bool=):
#    plt.figure(figsize=(15,45))
#    plot_branch(tree,0,0)
#plt.axis((-0.3,40,-9,0.3))
#plt.xticks(ticks=[0,1,2])
#plt.axes.tick_params(axis='x',)
#plot_branch(mytree,0,0)
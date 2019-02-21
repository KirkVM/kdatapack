import itertools,functools,collections
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from typing import Iterable
from ete3 import Tree
from matplotlib import markers
from matplotlib.path import Path
from operator import attrgetter
from dataclasses import dataclass
#from functools import reduce
#from itertools import product

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

def ete_set_branch_coordinates(branch,xcoord,ycoord,sepsize):
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
        ycoord_ascend,ycoord_descend=ete_set_branch_coordinates(sub_branch,xcoord,ycoord,sepsize)
        ycoord=ycoord_descend
        ycoord_ascends.append(ycoord_ascend)
    branch.add_feature('base_coord_offset', np.array([xcoord for x in range(len(ycoord_ascends))]))
    branch.add_feature('base_coord_span',np.array(ycoord_ascends))#[[xcoord for x in range(len(ycoord_ascends))],ycoord_ascends])
    ycoord_ascend=np.mean(ycoord_ascends)
    branch.add_feature('stem_coord_offset',np.array([ycoord_ascend,ycoord_ascend]))
    branch.add_feature('stem_coord_span',np.array([xcoord-branch.dist,xcoord]))#[[xcoord-branch.dist,xcoord],[ycoord_ascend,ycoord_ascend]])
    return ycoord_ascend,ycoord_descend
#from matplotlib.axes import Axes
#from matplotlib.transforms import IdentityTransform

@dataclass(repr=True)
class FigureCoordBox:
    xmin:float=None
    xmax:float=None
    ymin:float=None
    ymax:float=None
    def copy(self):
        return FigureCoordBox(self.xmin,self.xmax,self.ymin,self.ymax)

class NodeLayout:
    def __init__(self,ax=None,steml2d:mpl.lines.Line2D=None,basel2d:mpl.lines.Line2D=None,
                      orientation:str=None,leafspacing:float=None):
        self.ax=ax
        self.inverter=self.ax.transData.inverted()
        self.orientation=orientation
        figure_coords=steml2d.get_window_extent(self.ax)
        if orientation in ['left','right']:
            self.treebox_figureheight=np.abs(figure_coords.y1-figure_coords.y0)
        elif orientation in ['top','bottom']:
            self.treebox_figureheight=np.abs(figure_coords.x1-figure_coords.x0)
        self.stem_xycoords=self.inverter.transform(figure_coords)
        sycoord=np.array([self.stem_xycoords[0][1],self.stem_xycoords[1][1]])
        sxcoord=np.array([self.stem_xycoords[0][0],self.stem_xycoords[1][0]])
        if orientation in ['left','right']:
            stembox_y=np.concatenate(( np.array([x-0.5*leafspacing for x in sycoord]), np.array([x+0.5*leafspacing for x in sycoord]) ))
            stembox_x=np.concatenate((sxcoord,sxcoord))
        elif orientation in ['bottom','top']:
            stembox_x=np.concatenate(( np.array([x-0.5*leafspacing for x in sxcoord]), np.array([x+0.5*leafspacing for x in sxcoord]) ))
            stembox_y=np.concatenate((sycoord,sycoord))
        self.treebox=self.get_boundbox(stembox_x,stembox_y)
        self.branch_glyphs=[] #glyphs appended to end of stem
        self.align_glyphs=[] #glyphs appended to edge of tree (defined by all leaf nodes)
        self.branchbox=self.treebox.copy() #try this for now
        self.alignbox=self.treebox.copy()
        self.prep_next_branch_glyph()
       
    def add_branch_glyph_box(self,xypairs):
        branch_xarray=np.concatenate((np.array([self.branchbox.xmin,self.branchbox.xmax]),
                                       np.array([x[0] for x in xypairs])))
        branch_yarray=np.concatenate((np.array([self.branchbox.ymin,self.branchbox.ymax]),
                                       np.array([x[1] for x in xypairs])))
        self.branchbox=self.get_boundbox(branch_xarray,branch_yarray)
        self.prep_next_branch_glyph()

    def prep_next_branch_glyph(self):
        if self.orientation=='left':
            self.nextbranch_x=self.branchbox.xmax
        if self.orientation=='right':
            self.nextbranch_x=self.branchbox.xmin
        if self.orientation in ['left','right']:
            self.nextbranch_y=0.5*(self.branchbox.ymax+self.branchbox.ymin)

        if self.orientation=='top':
            self.nextbranch_y=self.branchbox.ymin
        if self.orientation=='bottom':
            self.nextbranch_y=self.branchbox.ymax
        if self.orientation in ['top','bottom']:
            self.nextbranch_x=0.5*(self.branchbox.xmax+self.branchbox.xmin)
  
    def get_boundbox(self,xarray,yarray):
        xmin=min(xarray)
        xmax=max(xarray)
        ymin=min(yarray)
        ymax=max(yarray)
        return FigureCoordBox(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)


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
        self.cviz_hadict={'left':'left','right':'right','top':'center','bottom':'center'}
        self.cviz_vadict={'left':'center','right':'center','top':'top','bottom':'bottom'}
        self.tree_lw=3.0
        self.tree_color='black'
        self.initial_leafspacing=0.1
        self.create_leaf_names=False
        self.draw_leaf_names=False
        
        self.cluster_feature='accs'
        self.set_cluster_size()

        self.plot_coords=[[np.inf,-np.inf],[np.inf,-np.inf]]
        self.ordered_leaves=None
        self.decorated_plot_coords=None

    def set_cluster_size(self):
        """ 
        Method to add .cluster_relsize feature to each node in the .tree
        Currently only supports if .cluster_feature=='accs'
        """
 
        if self.cluster_feature=='accs':
            cluster_sizes=[len(l.accs) for l in self.tree.get_leaves()]
            for l in self.tree.get_leaves():
                l.add_feature('cluster_relsize',np.log(len(l.accs)+1)/np.log(max(cluster_sizes)+1))

    def render(self,orientation=None,scale=None,figsize:Iterable=(20,20),initial_leafspacing:float=None,\
               tree_lw:float=None,tree_color:str=None,create_leaf_names:bool=None,draw_leaf_names:bool=None,\
                ax=None):
        """ method to generate figure.
        Calls: (1)ete_set_branch_coordinates to set node attributes,
               (2)self.get_plot_coords() to set plot_coords attribute
               (3)self.plot_it() to make the actual plot

        Keyword arguments:
            orientation: 'left','right','bottom', or 'top' (default is to use object property .orientation)
            scale: float setting figure scaling (default is to use object property .scale)
            figsize: size-2 tuple with hxw (default= (20,20))

        """
        #get all settings, falling back to class
        if initial_leafspacing is None:initial_leafspacing=self.initial_leafspacing
        if orientation is None: orientation=self.orientation
        if scale is None: scale=self.scale
        if tree_lw is None: tree_lw=self.tree_lw
        if tree_color is None: tree_color=self.tree_color
        if create_leaf_names is None: create_leaf_names=self.create_leaf_names
        if draw_leaf_names is None: draw_leaf_names=self.draw_leaf_names
        leafspacing=initial_leafspacing*scale
        #get branch relative positioning
        ete_set_branch_coordinates(self.tree,0,0,initial_leafspacing) #builds tree into ete3 tree
        if ax is None:
            plt.figure(figsize=figsize)
            ax=plt.gca()
        ax.get_figure().canvas.draw()
        #uses defaults and tree's coordinates set above to plot the tree
        self._plotit(ax,orientation,scale,leafspacing,tree_lw,tree_color,create_leaf_names,draw_leaf_names)

        #now get various utility attributes 
        self.ordered_leaves=sorted(self.tree.get_leaves(),key=lambda x:x.stem_coord_offset[0],reverse=True)


    def _plotit(self,ax,orientation,scale,leafspacing,tree_lw,tree_color,create_leaf_names,draw_leaf_names):
        """ method to generate figure. 
        Currently called from within render() method as that first calls necessary get_plot_coordinates() method
        """
        self.set_plot_coords(orientation,scale) #adjusts these using scale/orientation

        #plot the tree and create NodeLayout feature
        l2r=list(self.tree.traverse())[::-1] #sorted deepest to ro
        for node in l2r:
            stemplot=ax.plot(node.stem_plot_coords[0],node.stem_plot_coords[1],lw=tree_lw,color=tree_color)
            stemplot_l2d=stemplot[0]
            if node.is_leaf() is False:
                baseplot=ax.plot(node.base_plot_coords[0],node.base_plot_coords[1],lw=tree_lw,color=tree_color)
                baseplot_l2d=baseplot[0]
            else:
                baseplot_l2d=None
            node.add_feature('node_layout',NodeLayout(ax=ax,steml2d=stemplot_l2d,basel2d=baseplot_l2d,\
                              orientation=orientation,leafspacing=leafspacing))
        self.tree_plot_coords=[[ax.dataLim.x0,ax.dataLim.x1],[ax.dataLim.y0,ax.dataLim.y1]]

        fig_leafspacing=None
        converter=ax.transData.transform([leafspacing,leafspacing])
        if orientation in ['left','right']:
            fig_leafspacing=converter[1]
        elif orientation in ['top','bottom']:
            fig_leafspacing=converter[0]
        self.add_cluster_visualization(ax,orientation,leafspacing,fig_leafspacing)
        if draw_leaf_names:
            self.draw_lnames(ax,orientation,leafspacing,create_leaf_names,fig_leafspacing)
        if self.dashed_leaves:
            self.add_leaf_dashextensions(ax,orientation)
        self.trim_plotspace(ax,orientation,self.tree_plot_coords)#,self.decorated_plot_coords)  
        ax.set_xticks([])
        ax.set_yticks([])

    def draw_lnames(self,ax,orientation,leafspacing,create_leaf_names,fig_leafspacing):
        inverter=ax.transData.inverted()
        if create_leaf_names:
            name_iter=itertools.product('ABCDEFGHIJKLMNOPQRSRTUVWXYZ','ABCDEFGHIJKLMNOPQRSRTUVWXYZ')
            for lnode in self.tree.get_leaves():
                lnode.name=functools.reduce(lambda x,y:x+y,next(name_iter))
        texty_dict={}
        for lnode in self.tree.get_leaves():
            nl=lnode.node_layout
            #if orientation in ['left','right']:
            bbox_props=dict(boxstyle='square',fc='white',lw=0,alpha=0)
            #    texty=ax.text(nl.nextbranch_x,nl.nextbranch_y,lnode.name,bbox=bbox_props,\
            #                    ha='left',va='center',size=0.33*np.sqrt(fig_leafspacing))
            texty=ax.text(nl.nextbranch_x,nl.nextbranch_y,lnode.name,bbox=bbox_props,\
                           ha=self.cviz_hadict[orientation],va=self.cviz_vadict[orientation],size=0.5*np.sqrt(fig_leafspacing))
            texty_dict[nl]=texty
        #the .draw() here is required to get box locations, so run once then do adds to each box
        ax.get_figure().canvas.draw() 
        for nl in texty_dict.keys():
            tbox_ext=texty_dict[nl].get_bbox_patch().get_window_extent()
            twindow=inverter.transform(tbox_ext)
            nl.add_branch_glyph_box(twindow)

    def add_cluster_visualization(self,ax,orientation,leafspacing,fig_leafspacing):
        inverter=ax.transData.inverted()
        if self.cluster_feature=='accs':
            for node in self.tree.get_leaves():
                nl=node.node_layout
                cmarker=ax.plot(nl.nextbranch_x,nl.nextbranch_y,\
                    marker=align_marker(self.cviz_symboldict[orientation],halign=self.cviz_hadict[orientation],\
                    valign=self.cviz_vadict[orientation]),\
                     clip_on=False, color='k',ms=np.sqrt(fig_leafspacing)*node.cluster_relsize)
                tbbox=cmarker[0].get_tightbbox(ax)
                dabounds=inverter.transform(tbbox.corners())
                nl.add_branch_glyph_box(dabounds)

    def add_leaf_dashextensions(self,ax,orientation):
        inverter=ax.transData.inverted()
        nlxs=[lnode.node_layout.nextbranch_x for lnode in self.tree.get_leaves()]
        nlys=[lnode.node_layout.nextbranch_y for lnode in self.tree.get_leaves()]
        for lnode in self.tree.get_leaves():
            nl=lnode.node_layout
            if orientation=='left':
                xs=[nl.nextbranch_x,max(nlxs)]#stem_plot_coords[0][-1],plotcoords[0][1]]
                ys=[nl.nextbranch_y,nl.nextbranch_y]#stem_plot_coords[1][-1],lnode.stem_plot_coords[1][-1]]
            elif orientation=='right':
                xs=[nl.nextbranch_x,min(nlxs)]#stem_plot_coords[0][-1],plotcoords[0][1]]
                ys=[nl.nextbranch_y,nl.nextbranch_y]#stem_plot_coords[1][-1],lnode.stem_plot_coords[1][-1]]
            elif orientation=='top':
                xs=[nl.nextbranch_x,nl.nextbranch_x]#stem_plot_coords[0][-1],plotcoords[0][1]]
                ys=[nl.nextbranch_y,min(nlys)]#stem_plot_coords[1][-1],lnode.stem_plot_coords[1][-1]]
            elif orientation=='bottom':
                xs=[nl.nextbranch_x,nl.nextbranch_x]#stem_plot_coords[0][-1],plotcoords[0][1]]
                ys=[nl.nextbranch_y,max(nlys)]#stem_plot_coords[1][-1],lnode.stem_plot_coords[1][-1]]
            dashy=ax.plot(xs,ys,lw=1,ls='--',zorder=0,color='gray')
            tbbox=dashy[0].get_tightbbox(ax)
            dabounds=inverter.transform(tbbox.corners())
            nl.add_branch_glyph_box(dabounds)


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

    def trim_plotspace(self,ax,orientation,tree_plot_coords):
        xmin=tree_plot_coords[0][0]
        xmax=tree_plot_coords[0][1]
        ymin=tree_plot_coords[1][0]
        ymax=tree_plot_coords[1][1]
        for lnode in self.tree.get_leaves():
            nl=lnode.node_layout
            if orientation in ['left','right']:
                xmin=min(xmin,nl.branchbox.xmin)
                xmax=max(xmax,nl.branchbox.xmax)
                ymin=min(ymin,nl.treebox.ymin)
                ymax=max(ymax,nl.treebox.ymax)
            elif orientation in ['top','bottom']:
                xmin=min(xmin,nl.treebox.xmin)
                xmax=max(xmax,nl.treebox.xmax)
                ymin=min(ymin,nl.branchbox.ymin)
                ymax=max(ymax,nl.branchbox.ymax)
        coords=[[xmin,ymin],[xmax,ymax]]
        fig_coords=ax.transData.transform(coords)
        #fig_coords[0][0]-=2
        #fig_coords[0][1]-=2
        #fig_coords[1][0]+=2
        #fig_coords[1][1]+=2
        fig_coords=ax.transData.inverted().transform(fig_coords)
        ax.set_xlim(fig_coords[0][0],fig_coords[1][0])#cur_mins[0],cur_maxes[0])#,transform=IdentityTransform)#,ax.transData)#self.decorated_plot_coords[0][0],self.decorated_plot_coords[0][1]
        ax.set_ylim(fig_coords[0][1],fig_coords[1][1])#cur_mins[0],cur_maxes[0])#,transform=IdentityTransform)#,ax.transData)#self.decorated_plot_coords[0][0],self.decorated_plot_coords[0][1]
#        ax.set_ylim(ymin,ymax)#tree_plot_coords[1][0]-0.5*leafspacing,tree_plot_coords[1][1]+0.5*leafspacing)
#        xmin=ax.transData.transform()
#        ax.set_xlim(xmin,xmax)#cur_mins[0],cur_maxes[0])#,transform=IdentityTransform)#,ax.transData)#self.decorated_plot_coords[0][0],self.decorated_plot_coords[0][1]
#        ax.set_ylim(ymin,ymax)#tree_plot_coords[1][0]-0.5*leafspacing,tree_plot_coords[1][1]+0.5*leafspacing)


#plt.gca().set_xlim(0.0,3.3)
#plt.gca().spines['left'].set_visible(False)
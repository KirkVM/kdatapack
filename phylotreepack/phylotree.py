import dask
from ete3 import Tree
from dataclasses import dataclass
import numpy as np

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

@dataclass(repr=True)
class FigureCoordBox:
    xmin:float=None
    xmax:float=None
    ymin:float=None
    ymax:float=None
    def copy(self):
        return FigureCoordBox(self.xmin,self.xmax,self.ymin,self.ymax)

@dataclass
class FrameCoords:
    stem_start:(float,float)=(None,None)
    stem_end:(float,float)=(None,None)
    base_start:(float,float)=(None,None)
    base_end:(float,float)=(None,None)
    orientation:str='left'

def set_branch_coordinates(branch,xcoord,ycoord,sepsize):
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
        branch.framecoords=FrameCoords(stem_start=(xcoord,ycoord),stem_end=(xcoord+branch.dist,ycoord)) 
        return ycoord,ycoord-sepsize
    else:
        xcoord=xcoord+branch.dist
        sub_branches=branch.children
        #sorted_sub_branches=depth_sort(sub_branches)
    ycoord_ascends=[]
    for sub_branch in sub_branches:#sorted_sub_branches:
        ycoord_ascend,ycoord_descend=set_branch_coordinates(sub_branch,xcoord,ycoord,sepsize)
        ycoord=ycoord_descend
        ycoord_ascends.append(ycoord_ascend)

    branch.framecoords=FrameCoords(stem_start=(xcoord-branch.dist,np.mean(ycoord_ascends)),
                                   stem_end=(xcoord,np.mean(ycoord_ascends)),
                                    base_start=(xcoord,min(ycoord_ascends)),
                                    base_end=(xcoord,max(ycoord_ascends)))
    return np.mean(ycoord_ascends),ycoord_descend


class PTNodeGlyph:
    def __init__(self,glyph):
        self.boundbox=None

class PhyloTree:
    def __init__(self,etenode,depth=0):
        if depth==0:
            self.etenode=etenode.copy() #make a copy of top node
        else:
            self.etenode=etenode
        self.depth=depth
        self.dist=etenode.dist
        self.etenode.add_feature('ptnode',self)
        self.framecoords=None
        self.glyphs=None
        self.branchbox=None
        self.alignbox=None
        etechildren=self.etenode.children
        self.children=[PhyloTree(etekid,depth=self.depth+1) for etekid in etechildren] 
    def is_leaf(self):
        if len(self.children)==0: 
            return True 
        else:
            return False
    def get_leaves(self):
        eteleaves=self.etenode.get_leaves()
        return [eteleaf.ptnode for eteleaf in eteleaves]
    def draw(self,sepsize=0.1,ax=None):
        set_branch_coordinates(self,0,0,sepsize)
#        drawtasks=[]
        for etenode in self.etenode.traverse():
            ptn=etenode.ptnode
#            xs=[ptn.framecoords.stem_start[0],ptn.framecoords.stem_end[0]]
#            ys=[ptn.framecoords.stem_start[1],ptn.framecoords.stem_end[1]] 
#            drawtasks.append(dask.delayed(ax.plot,pure=True)(xs,ys))
            ax.plot([ptn.framecoords.stem_start[0],ptn.framecoords.stem_end[0]],[ptn.framecoords.stem_start[1],ptn.framecoords.stem_end[1]]) 
#        dask.compute(drawtasks)

#class PhyloTree:
#    """Class for plotting an ete3 tree using matplotlib rendering
#
#    Attributes:
#        orientation ('left','right','top','bottom'): tree orientation (str, default 'left')
#        dashed_leaves: whether to add dashes from leaves to edge of graph (default True)
#        cluster_feature: ete node feature to use to define cluster size (str, default 'accs')
#        cluster_viz: symbol to use as cluster indicator (default 'triangle')---not yet impld
#        ordered_leaves: list of leaves from graph start (in ladderized form) to end
#                        (calculated through calling .render() method)
#        plot_coords: [xmin,xmax],[ymin,ymax] values for graph
#                        (calculated through calling .render() method)
#    """
#    def __init__(self,tree:Tree):#,cluster_feature='accs'):
#        """ constructor for EteMplTree.
#        Calls: self.cluster_size() to add feature .cluster_relsize to each leaf node
#
#        Arguments:
#            tree: an ete3 tree instance
#            [Also currently makes some assumed settings that are configurable public properties-
#            -orientation,cluster_feature,cluster_viz,scale]
#        
#        Keyword Arguments:
#            cluster_feature: what feature to use as indication of cluster size: None, 'accs' (default='accs')
#        """
#
#        self.tree=tree.copy()
#        self.orientation='left'
#        self.scale=1.0
#        self.dashed_leaves=True
#        self.cluster_viz='triangle'
#        self.cviz_symboldict={'left':'<','right':'>','top':'^','bottom':'v'}
#        self.cviz_hadict={'left':'left','right':'right','top':'center','bottom':'center'}
#        self.cviz_vadict={'left':'center','right':'center','top':'top','bottom':'bottom'}
#        self.tree_lw=3.0
#        self.tree_color='black'
#        self.initial_leafspacing=0.1

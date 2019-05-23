#import dask
from ete3 import Tree
from dataclasses import dataclass
import numpy as np
import itertools
from matplotlib.path import Path
import matplotlib.patches as patches

from bokeh.models import ColumnDataSource
from bokeh.models.glyphs import Line,MultiLine

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

def trotater(rotation,*xyvals):
    rotate_angle=(rotation/360)*2*np.pi
    rmatrix=np.array(([np.cos(rotate_angle),-np.sin(rotate_angle)],\
                        [np.sin(rotate_angle),np.cos(rotate_angle)]))
    retxys=[]
    for xyval in xyvals:
        if xyval[0] is None:
            retxys.append(None)
        else:
            retxys.append(  np.array(xyval)@rmatrix )
    return retxys

class PTCoords:
    def __init__(self,framecoords,rotation):
        self.framecoords=framecoords
        self.rotation=None
        self.stem_xys=[] #for matplotlib xys
        self.base_xys=[]
        self.mpl_path_verts=[]
        self.mpl_path_codes=[]
        self.stem_xvals=[] #for bokeh lines
        self.stem_yvals=[]
        self.base_xvals=[]
        self.base_yvals=[]
        self.xmin=None
        self.ymin=None
        self.xmax=None
        self.ymax=None
        self.framepath=None
        self.set_coords(rotation)
    def set_coords(self,rotation):
        if self.rotation==rotation:
            return
        if self.framecoords.stem_start[0] is not None:
            self.stem_xys=trotater(rotation,self.framecoords.stem_start,self.framecoords.stem_end)
        if self.framecoords.base_start[0] is not None:
            self.base_xys=trotater(rotation,self.framecoords.base_start,self.framecoords.base_end)
        self.stem_xvals = [x[0] for x in self.stem_xys]
        self.stem_yvals = [x[1] for x in self.stem_xys]
        self.base_xvals = [x[0] for x in self.base_xys]
        self.base_yvals = [x[1] for x in self.base_xys]
        self.mpl_path_verts=self.stem_xys+self.base_xys
        code_cycle=itertools.cycle([Path.MOVETO,Path.LINETO])
        self.mpl_path_codes=[next(code_cycle) for x in range(len(self.mpl_path_verts))]
#        self.framepath=Path(self.mpl_path_verts,self.mpl_path_codes)
        self.xmin=min([x[0] for x in self.stem_xys+self.base_xys])
        self.xmax=max([x[0] for x in self.stem_xys+self.base_xys])
        self.ymin=min([x[1] for x in self.stem_xys+self.base_xys])
        self.ymax=max([x[1] for x in self.stem_xys+self.base_xys])
        self.rotation=rotation #now it's set so can skip next time

@dataclass
class FrameCoords:
    '''stores coordinates of tree bases,stems in left-to-right rectangular form'''
    stem_start:(float,float)=(None,None)
    stem_end:(float,float)=(None,None)
    base_start:(float,float)=(None,None)
    base_end:(float,float)=(None,None)


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
    #import pdb;pdb.set_trace()
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
        self.name=etenode.name
        self.depth=depth
        self.dist=etenode.dist
        self.etenode.add_feature('ptnode',self)
        self.framecoords=None
        self.ptcoords=None
        self.glyphs=None
        self.branchbox=None
        self.alignbox=None
        self.cds=None
        etechildren=self.etenode.children
        self.children=[PhyloTree(etekid,depth=self.depth+1) for etekid in etechildren] 
        self.set_coords(0.1)
        
        self.cds_dict={}#'name':self.name}
        if self.is_leaf():
            self.cds_dict.update({'gbacc':[self.name,self.name]})
    def set_coords(self,sepsize):
        set_branch_coordinates(self,0,0,sepsize)
        for etenode in self.etenode.traverse():
            ptn=etenode.ptnode
            ptn.ptcoords=PTCoords(ptn.framecoords,rotation=0)

    def is_leaf(self):
        if len(self.children)==0: 
            return True 
        else:
            return False
    def get_leaves(self):
        eteleaves=self.etenode.get_leaves()
        return [eteleaf.ptnode for eteleaf in eteleaves]
    def mpldraw(self,sepsize=0.1,ax=None,rotation=0):
#        set_branch_coordinates(self,0,0,sepsize)
        for etenode in self.etenode.traverse():
            ptn=etenode.ptnode
            ptn.ptcoords.set_coords(rotation)
            framepathobj=Path(ptn.ptcoords.mpl_path_verts,ptn.ptcoords.mpl_path_codes)
            patch=patches.PathPatch(framepathobj)
            ax.add_patch(patch)
        all_ptnc=[x.ptnode.ptcoords for x in self.etenode.traverse()]
        xmin=min([x.xmin for x in all_ptnc])
        xmax=max([x.xmax for x in all_ptnc])
        ymin=min([x.ymin for x in all_ptnc])
        ymax=max([x.ymax for x in all_ptnc])
        ax.set_ylim((ymin,ymax))
        ax.set_xlim((xmin,xmax))

    def bokehdraw(self,sepsize=0.2,plot=None,rotation=0):
#        set_branch_coordinates(self,0,0,sepsize)
        for etenode in self.etenode.traverse():
            ptn=etenode.ptnode
            ptn.ptcoords.set_coords(rotation)
            ptn.cds=ColumnDataSource(ptn.cds_dict)
            if ptn.is_leaf():
                ntype='leaf_node'
            else:
                ntype='internal_node'
#            source=ColumnDataSource({'xs':[ptn.ptcoords.stem_xvals,ptn.ptcoords.base_xvals],\
#                                    'ys':[ptn.ptcoords.stem_yvals,ptn.ptcoords.base_yvals]})
            if len(ptn.ptcoords.stem_xvals)==2:
                ptn.cds.add(ptn.ptcoords.stem_xvals,'stem_xs')
                ptn.cds.add(ptn.ptcoords.stem_yvals,'stem_ys')
                sglyph=Line(x='stem_xs',y='stem_ys')
#``                plot.add_glyph(ptn.cds,sglyph)
                plot.add_glyph(ptn.cds,sglyph,name=ntype)
            if len(ptn.ptcoords.base_xvals)==2:
                ptn.cds.add(ptn.ptcoords.base_xvals,'base_xs')
                ptn.cds.add(ptn.ptcoords.base_yvals,'base_ys')
                bglyph=Line(x='base_xs',y='base_ys')
                plot.add_glyph(ptn.cds,bglyph,name=ntype)
#            fglyph=MultiLine(xs='xs',ys='ys')
#            plot.add_glyph(source,fglyph)
#            cds=ColumnDataSource(dict(x=))


#
#                
#                [ptn.framecoords.stem_start[0],ptn.framecoords.stem_end[0]],[ptn.framecoords.stem_start[1],ptn.framecoords.stem_end[1]]) 
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

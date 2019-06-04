#import dask
import itertools,math,os,sqlite3
from ete3 import Tree
from dataclasses import dataclass
import numpy as np
from matplotlib.path import Path
import matplotlib.patches as patches

from bokeh.models import ColumnDataSource,Label
from bokeh.models.glyphs import Line,MultiLine,Quad,Patch

@dataclass(repr=True)
class FigBoundBox:
    xmin:float=None
    xmax:float=None
    ymin:float=None
    ymax:float=None
    def copy(self):
        return FigBoundBox(self.xmin,self.xmax,self.ymin,self.ymax)

class PathXYCoords:
    def __init__(self,xvals=[],yvals=[],xyvals=[],rotation=0):
        self.xyvals=[]
        self.xvals=[]
        self.yvals=[]
        if len(xyvals)>0:
            assert (len(xvals)==0 and len(yvals)==0), \
                "if xyvals provided, cannot provide xvals or yvals"
            self.xyvals=[(xy[0],xy[1]) for xy in xyvals] #make sure it's a list of tuples
            self.xvals=[xy[0] for xy in self.xyvals] 
            self.yvals=[xy[1] for xy in self.xyvals] 
        if len(xvals)>0:
            assert ((len(xvals)==len(yvals)) and len(xyvals)==0),\
                "xvals and yvals lengths do not match"
            self.xvals=xvals
            self.yvals=yvals
            self.xyvals=[(x,y) for x,y in zip(self.xvals,self.yvals)]

        self.boundbox=FigBoundBox(xmin=min(self.xvals),xmax=max(self.xvals),\
                                  ymin=min(self.yvals),ymax=max(self.yvals) )
        self.rotation=rotation 
        self.orderxy_clockwise()
    def orderxy_clockwise(self):
        centroid_x=np.sum(self.xvals)/len(self.xvals)
        centroid_y=np.sum(self.yvals)/len(self.yvals)
        xy_sorted=sorted(self.xyvals,key=lambda x:math.atan2((x[1]-centroid_y),(x[0]-centroid_x)))
        minx_idx=np.argmin([x[0] for x in xy_sorted])
        self.xyvals=xy_sorted[minx_idx:]+xy_sorted[:minx_idx]
    def copy(self):
        return PathXYCoords(xyvals=self.xyvals,rotation=self.rotation)
    def apply_rotation(self,rotation):
        if rotation==self.rotation:
            return
        rotate_angle=(rotation/360)*2*np.pi
        rmatrix=np.array(([np.cos(rotate_angle),-np.sin(rotate_angle)],\
                        [np.sin(rotate_angle),np.cos(rotate_angle)]))
        rotxys=[]
        for xy in self.xyvals:
            rotxys.append( np.array(xy)@rmatrix )
        self.xyvals=rotxys
        self.xvals=[xy[0] for xy in self.xyvals] 
        self.yvals=[xy[1] for xy in self.xyvals] 
        self.boundbox=FigBoundBox(xmin=min(self.xvals),xmax=max(self.xvals),\
                                  ymin=min(self.yvals),ymax=max(self.yvals) )
        self.rotation=rotation 
        self.orderxy_clockwise()

class FrameCoords:
    '''stores coordinates of tree bases,stems in left-to-right rectangular form'''
    def __init__(self,stem=None,base=None,rotation=0):
        self.stem=stem
        self.base=base
        self.rotation=rotation
    def copy(self):
        if ((self.base is not None) and (self.stem is not None)):
            return FrameCoords(stem=self.stem.copy(),base=self.base.copy())
        elif self.base is None:
            return FrameCoords(stem=self.stem.copy())
        else:
            return FrameCoords(base=self.base.copy())
    def apply_rotation(self,rotation):
        if self.stem is not None:
            self.stem.apply_rotation(rotation)
        if self.base is not None:
            self.base.apply_rotation(rotation)
        self.rotation=rotation

def set_branch_coordinates(branch,xcoord,ycoord,sepsize,tdistmax):
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
        stempath=PathXYCoords( xvals=[xcoord,xcoord+branch.dist],yvals=[ycoord,ycoord] )
        branch.r0framecoords=FrameCoords(stem=stempath)
        #nodebox_xmax=max(xcoord+branch.dist,xcoord+tdistmax*0.025)
        #nodebox_xmax=max(xcoord+branch.dist,1.025*tdistmax)
        nodebox_xmax=1.025*tdistmax
        branch.r0node_edgecoords=PathXYCoords(xyvals=[(xcoord,ycoord-0.5*sepsize),(nodebox_xmax,ycoord-0.5*sepsize),
                                            (nodebox_xmax,ycoord+0.5*sepsize),(xcoord,ycoord+0.5*sepsize)])
        branch.r0branch_edgecoords=branch.r0node_edgecoords.copy()
        return ycoord,ycoord-sepsize
    else:
        xcoord=xcoord+branch.dist
        sub_branches=branch.children
        #sorted_sub_branches=depth_sort(sub_branches)
    ycoord_ascends=[]
    for sub_branch in sub_branches:#sorted_sub_branches:
        ycoord_ascend,ycoord_descend=set_branch_coordinates(sub_branch,xcoord,ycoord,sepsize,tdistmax)
        ycoord=ycoord_descend
        ycoord_ascends.append(ycoord_ascend)

    branchstempath=PathXYCoords( xvals=[xcoord-branch.dist,xcoord],
                                 yvals=[np.mean(ycoord_ascends),np.mean(ycoord_ascends)] )
    branchbasepath=PathXYCoords( xvals=[xcoord,xcoord], 
                                 yvals=[min(ycoord_ascends),max(ycoord_ascends)] )
    branch.r0framecoords=FrameCoords(stem=branchstempath,base=branchbasepath) 
    node_edgecoords=PathXYCoords(xyvals=[(xcoord-branch.dist,min(ycoord_ascends)),(xcoord,min(ycoord_ascends)),
                                                (xcoord,max(ycoord_ascends)),(xcoord-branch.dist,max(ycoord_ascends)) ])
    branch.r0node_edgecoords=node_edgecoords
    branch_ymin=min(node_edgecoords.boundbox.ymin,*[subx.r0branch_edgecoords.boundbox.ymin for subx in branch.children])
    branch_ymax=max(node_edgecoords.boundbox.ymax,*[subx.r0branch_edgecoords.boundbox.ymax for subx in branch.children])
    branch_xmin=min(node_edgecoords.boundbox.xmin,*[subx.r0branch_edgecoords.boundbox.xmin for subx in branch.children])
    branch_xmax=max(node_edgecoords.boundbox.xmax,*[subx.r0branch_edgecoords.boundbox.xmax for subx in branch.children])
    branch.r0branch_edgecoords=PathXYCoords(xyvals=[(branch_xmin,branch_ymin),(branch_xmax,branch_ymin),
                                                (branch_xmax,branch_ymax),(branch_xmin,branch_ymax)] )   

    return np.mean(ycoord_ascends),ycoord_descend

class PTGlyph:
    def __init__(self,name,type,position='edge'):
        self.name=name
        self.position=position
        self.boundbox=None
    def get_rendering(self,plot,ptree,rotation):
#        align=TextAlign(15,20,5)
        label=Label(x=ptree.ptcoords.stem_xys[1][0],y=ptree.ptcoords.stem_xys[1][1],text=ptree.decoration_dict[self.name],text_align='center',text_baseline='middle')
        plot.add_layout(label)
#        plot.label(x=ptree.ptcoords.stem_xys[1][0],y=ptree.ptcoords.stem_xys[1][1],text=ptree.decoration_dict[name])
#        pass

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
        etechildren=self.etenode.children
        #recursively add children
        self.children=[PhyloTree(etekid,depth=self.depth+1) for etekid in etechildren] 
        
        if self.is_leaf():
            self.ntype='leaf_node' 
        else:
            self.ntype='internal_node'
        self.leaf_cds=None
        self.internal_cds=None
        self.leaf_cds_dict={}#'name':self.name}
        self.internal_cds_dict={}#'name':self.name}

        #################
        #these will stay as left-to-right-oriented ('rotation0'='r0') base tree
        self.r0framecoords=None
        self.r0branch_edgecoords=None
        self.r0node_edgecoords=None
        self.r0leaf_dashcoords=None
        #and these are the rotatable
        self.framecoords=None
        self.branch_edgecoords=None
        self.node_edgecoords=None
        self.leaf_dashcoords=None
        if self.depth==0:
            #self.set_r0coords(0.1) #values are set in here...
            set_branch_coordinates(self,0,0,0.1,self.etenode.get_farthest_leaf()[1]+self.dist)
            #create plot details with rotation=0
            rhe=self.r0branch_edgecoords.boundbox.xmax
            for ptn in self.traverse():
                if ptn.is_leaf():
                    ptn.r0leaf_dashcoords=PathXYCoords(xyvals=[ (ptn.r0framecoords.stem.boundbox.xmax,ptn.r0framecoords.stem.boundbox.ymax),
                                                        (rhe,ptn.r0framecoords.stem.boundbox.ymax)] )
                ptn.framecoords=ptn.r0framecoords.copy()
                ptn.branch_edgecoords=ptn.r0branch_edgecoords.copy()
                ptn.node_edgecoords=ptn.r0node_edgecoords.copy()
                if ptn.is_leaf():
                    ptn.leaf_dashcoords=ptn.r0leaf_dashcoords.copy()
            self.leaf_cds_dict['gbacc']=[ptn.name for ptn in self.traverse() if ptn.ntype=='leaf_node']#append(ptn.name#[ptn.name for ptn in self.traverse()]
#            self.cds_dict['ntype']=[ptn.ntype for ptn in self.traverse()]
        ####################################################
##        self.branch_glyphcoords=[] #deprecate?
        #these are coordinates to plot        
        #self.ptgcoords=[]
#        self.ptannotations=[]
#        self.ptglyphs=[] #deprecate?
#
#        self.alignbox=None
##        
#        self.decoration_dict={}
#        if self.is_leaf():
#            #self.cds_dict.update({'gbacc':[self.name,self.name]})
#            self.cds_dict.update({'gbacc':[self.name]})#,self.name]})

#    def set_r0coords(self,sepsize):
#        set_branch_coordinates(self,0,0,sepsize)
#        for ptn in self.traverse():
#            ptn.ptcoords=PTCoords(ptn.r0framecoords,rotation=0)


    ######\/##\/##ETE WRAPPERS#\/#\/#####
    def is_leaf(self): #~wrapper
        return len(self.children)==0
    def get_leaves(self): #ETEWRAPPER
        eteleaves=self.etenode.get_leaves()
        return [eteleaf.ptnode for eteleaf in eteleaves]
    def get_leaf_names(self): #ETE WRAPPER
        return self.etenode.get_leaf_names()
    def traverse(self): #ETE WRAPPER
        '''yields a generator over all nodes by wrapping ete's traverse() method'''
        for x in self.etenode.traverse():
            yield x.ptnode
    #####^^^^^ETE WRAPPERS^^^^########

    def update_leafcdsdict_fromdb(self,dbpathstr,fields=['pdbids','ecs','subfam','extragbs'],searchby='gbacc'):
        assert(os.path.exists(dbpathstr))
        conn=sqlite3.connect('GH5/GH5DB.sql')
        conn.row_factory=sqlite3.Row
        dbcursor=conn.cursor()
 
        assert (searchby=='gbacc')
        for field in fields:
            self.leaf_cds_dict[field]=[]
        for gbacc in self.leaf_cds_dict['gbacc']:
            dbcursor.execute('''SELECT * FROM CAZYSEQDATA WHERE acc=(?)''',(gbacc,))
            row=dbcursor.fetchone()
            for field in fields:
                if row is None and row[field] is not None:
                    self.leaf_cds_dict[field].append(None)
                else:
                    self.leaf_cds_dict[field].append(row[field])
        conn.close()
#                if row[field] is None: continue
#                
#                lnode.decoration_dict[field]=row[field]
#    
#    def add_leaf_decoration(self,keyname):
#        for lnode in self.get_leaves():
#            if keyname in [x.name for x in lnode.ptglyphs]: continue #it's already there
#            lnode.ptglyphs.append(PTGlyph(keyname,'annotation'))
    def apply_rotation(self,rotation):
        if rotation!=self.framecoords.rotation:
            self.framecoords=self.r0framecoords.copy()
            self.framecoords.apply_rotation(rotation)
        if rotation!=self.branch_edgecoords.rotation:
            self.branch_edgecoords=self.r0branch_edgecoords.copy() 
            self.branch_edgecoords.apply_rotation(rotation)
        if rotation!=self.node_edgecoords.rotation:
            self.node_edgecoords=self.r0node_edgecoords.copy() 
            self.node_edgecoords.apply_rotation(rotation)

    def mpldraw(self,sepsize=0.1,ax=None,rotation=0):
        for ptn in self.traverse():
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
        for ptn in self.traverse():
            ptn.apply_rotation(rotation)
        leaf_frame_xs=[]
        leaf_frame_ys=[]
        int_frame_xs=[]
        int_frame_ys=[]
        self.leaf_cds=ColumnDataSource(self.leaf_cds_dict)
        self.internal_cds=ColumnDataSource(self.internal_cds_dict)
        #plot the frame as one MultiLine
        for ptn in self.traverse():
            frame_xs=[]
            frame_ys=[]
            if ptn.framecoords.stem is not None: #has stem?
                frame_xs.append(ptn.framecoords.stem.xvals)
                frame_ys.append(ptn.framecoords.stem.yvals)
            if ptn.framecoords.base is not None:
                frame_xs.append(ptn.framecoords.base.xvals)
                frame_ys.append(ptn.framecoords.base.yvals)
            if ptn.ntype=='leaf_node':
                leaf_frame_xs.extend(frame_xs)
                leaf_frame_ys.extend(frame_ys)
            else:
                int_frame_xs.extend(frame_xs)
                int_frame_ys.extend(frame_ys)
        self.leaf_cds.add(leaf_frame_xs,'frame_xs')
        self.leaf_cds.add(leaf_frame_ys,'frame_ys')
        self.internal_cds.add(int_frame_xs,'frame_xs')
        self.internal_cds.add(int_frame_ys,'frame_ys')
        fglyph=MultiLine(xs='frame_xs',ys='frame_ys')
        plot.add_glyph(self.leaf_cds,fglyph)#,name='leaf_node')
        plot.add_glyph(self.internal_cds,fglyph,name='internal_node')

        qlefts=[]
        qrights=[]
        qtops=[]
        qbottoms=[]
        for ptn in self.traverse():
            if ptn.ntype=='leaf_node':
                qlefts.append(ptn.node_edgecoords.boundbox.xmin)
                qrights.append(ptn.node_edgecoords.boundbox.xmax)
                qtops.append(ptn.node_edgecoords.boundbox.ymax)
                qbottoms.append(ptn.node_edgecoords.boundbox.ymin)

        self.leaf_cds.add(qlefts,'nodebox_lefts')
        self.leaf_cds.add(qrights,'nodebox_rights')
        self.leaf_cds.add(qtops,'nodebox_tops')
        self.leaf_cds.add(qbottoms,'nodebox_bottoms')
        qglyph=Quad(left='nodebox_lefts',right='nodebox_rights',bottom='nodebox_bottoms',top='nodebox_tops',fill_color='#b3de69',fill_alpha=0,line_alpha=0)                
        plot.add_glyph(self.leaf_cds,qglyph,name='leaf_node')#'leaf_node')#self.ntype)
#        print(len(self.cds.data['frame_xs']))#,len(self.cds.data['gbacc']))
#            
#            for x in ptn.decoration_dict:
#                ptn.cds_dict[x]=[ptn.decoration_dict[x]]
#            ptn.cds=ColumnDataSource(ptn.cds_dict)
#            xs2plot=[]
#            ys2plot=[]
#            if ptn.framecoords.stem is not None: #has stem?
#                xs2plot.append(ptn.framecoords.stem.xvals)
#                ys2plot.append(ptn.framecoords.stem.yvals)
#            if ptn.framecoords.base is not None:
#                xs2plot.append(ptn.framecoords.base.xvals)
#                ys2plot.append(ptn.framecoords.base.yvals)
#            ptn.cds.add(xs2plot,'frame_xs')
#            ptn.cds.add(ys2plot,'frame_ys')
#            fglyph=MultiLine(xs='frame_xs',ys='frame_ys')
#            plot.add_glyph(ptn.cds,fglyph,name=ptn.ntype)#'leaf_node')#self.ntype)
#            if ptn.ntype=='leaf_node':
#                ptn.cds.add([ptn.node_edgecoords.boundbox.xmin],'left')
#                ptn.cds.add([ptn.node_edgecoords.boundbox.xmax],'right')
#                ptn.cds.add([ptn.node_edgecoords.boundbox.ymin],'bottom')
#                ptn.cds.add([ptn.node_edgecoords.boundbox.ymax],'top')
#                qglyph=Quad(left='left',right='right',bottom='bottom',top='top',fill_color='#b3de69',fill_alpha=0,line_alpha=0)                
#                plot.add_glyph(ptn.cds,qglyph,name=ptn.ntype)
#
#                ptn.cds.add([ptn.leaf_dashcoords.xvals],'dash_xs')
#                ptn.cds.add([ptn.leaf_dashcoords.yvals],'dash_ys')
#                dglyph=MultiLine(xs='dash_xs',ys='dash_ys')
#                plot.add_glyph(ptn.cds,dglyph,name=ptn.ntype)
#        for lnode in self.get_leaves():
#            for ptg in lnode.ptglyphs:
#                ptg.get_rendering(plot,lnode,rotation)
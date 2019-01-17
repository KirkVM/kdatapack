from scipy.optimize import curve_fit
import numpy as np
from .fitstats import FitStats
from typing import List,Iterable,Callable

class NLRegression:
    """Wraps scipy curve_fit and stores fit stats and visualization data"""
    def __init__(self,fitfunc:callable,xvals:List[float],yvals:List[float]):
        """"Builds NLRegression"""
        self.xvals=xvals
        self.yvals=yvals
        self.fitfunc=fitfunc

        self.guess_params=None
        self.yvals_sigma=None
        self.popt=None
        self.pcov=None
        self.fit_r2=None
        self.fitstats=FitStats()
        self.ellipsoid_coords=None
        self.ellipsoid_nstd=None
    def do_fit(self,guess_params: List[float] = None, yvals_sigma: List[float] = None):
        """performs fit using scipy curve_fit, and generates fit stats
        such as fit_r2, fit_sigmas"""
        self.guess_params=guess_params
        self.yvals_sigma=yvals_sigma
        self.popt,self.pcov=curve_fit(self.fitfunc,self.xvals,self.yvals,p0=self.guess_params,
                                    sigma=self.yvals_sigma)
        sstot=np.sum([np.power((y-np.mean(self.yvals)),2) for y in self.yvals])
        ssres=np.sum([np.power((y-self.fitfunc(x,*self.popt)),2)\
                        for x,y in zip(self.xvals,self.yvals)])
        fit_r2=1.0-ssres/sstot
        self.fitstats.r2=fit_r2
        fit_sigmas=[np.sqrt(self.pcov[x,x]) for x in range(np.shape(self.pcov)[0]) ]
        self.fitstats.param_sigma=fit_sigmas
    def get_error_ellipsoid(self,nstd:float=1.0,plotlen:int=100):
        """ returns array with x,y,(opt z) coordinates for error ellipse

            Keyword Arguments:
            nstd -- number of standard deviations for contour (default 1.0)
            plotlen -- number of points for each dimension (default 100)
            Returns:
            array with cartesian plot coordinates either 2d or 3d
        """
        ndim=len(self.popt) #for now
        if ndim>3:
            #could implement ability to sub-select...?
            print("max 3d ellipsoid please")
            return None,None,None
        elif ndim==1:
            print("only 1 variable, you don't need this!")
            return None,None,None
        U,s,Vh=np.linalg.svd(self.pcov)
        if ndim==3:
            w,l,h=nstd*np.sqrt(s)     #radius on the x-axis
            u = np.linspace(0, 2 * np.pi, plotlen)
            v = np.linspace(0, np.pi, plotlen)
            #set up ellipsoid
            x = w * np.outer(np.cos(u), np.sin(v))
            y = l * np.outer(np.sin(u), np.sin(v))
            z = h * np.outer(np.ones(np.size(u)), np.cos(v))
        #apply rotation in U
            for i in range(len(x)):
                for j in range(len(x)):
                    [x[i,j],y[i,j],z[i,j]] = U@([x[i,j],y[i,j],z[i,j]])# + center
        #shift results using popt
            x+=self.popt[0]
            y+=self.popt[1]
            z+=self.popt[2]
            self.ellipsoid_coords=[x,y,z]
            self.ellipsoid_nstd=nstd
            return x,y,z
        elif ndim==2:
            w,l=nstd*np.sqrt(s)
            u = np.linspace(0, 2 * np.pi, 100)
            x = w*np.cos(u)
            y = l*np.sin(u)
            for i in range(len(x)):
                [x[i],y[i]] = U@([x[i],y[i]])# + center
            x+=self.popt[0]
            y+=self.popt[1]
            self.ellipsoid_coords=[x,y]
            self.ellipsoid_nstd=nstd
            return x,y

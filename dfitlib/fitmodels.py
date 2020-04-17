import math
import numpy as np
from typing import Callable,Sequence

###############MODELS################################
def linear_model(xvals:np.ndarray,*coefs:float):
    '''returns linear function of x values. 
    
    arguments:
        xvals: ndarray preprocessed and of shape (xlen,D)
        coefs: float parameter estimates. num coefs should==D
    
    returns:
        dot product xvals@coefs
    
    example:
        xvals w/ shape 4,3: [1,1,1],[1,2,4],[1,3,9],[1,10,-5] (4 datapoints w/ 3 dimensions)
        coefs: [0,1,2]
        xvals@coefs=0*[1 1 1 1]+ 1*[1 2 3 10], 2*[1 4 9 -5] = [3, 10, 21, 0] ie [3 10 21 0]
    '''
    assert(xvals.shape[1]==len(coefs)), "need # of coefs to match # data points (ie 2nd dim of xvals)"
    
    yvals=xvals@coefs
    return yvals

def logistic_simple(xvals:np.ndarray,amp:float,offset:float,gamma:float,theta:float):
    '''the sorta-intuitive version of the logistic equation: 
            offset + amp*(1/1+exp(-gamma(x-theta)))
    
    arguments:
        xvals: 1-d np array (length n)
        amp: amplitude term 
        offset: y-offset in data
        gamma: 'sharpness' of transition
        theta: lag
    
    returns: length-n array of corresponding dependent (ie y) values
    '''
    assert type(xvals)==np.ndarray
    logistic_core=1/(1+np.exp(-gamma*(xvals-theta)))
    return offset+amp*logistic_core

def michaelismenten(xvals=np.ndarray,km=float,vmax=float):
    return vmax*xvals/(km+xvals)

def poisson_distribution(xvals=np.ndarray,lam=float):
    return np.power(lam,xvals)*np.exp(-1*lam)/np.array([math.factorial(x) for x in xvals])


#####lmfit wrappers###############
def lmlogistic_simple(pars,modelfunc,xvals,yvals):
    amp=pars['amp']
    offset=pars['offset']
    gamma=pars['gamma']
    theta=pars['theta']
    
    predicted=logistic_simple(xvals,amp,offset,gamma,theta)
    return yvals-predicted

###########scipy.optimize.minimize wrappers###########
def l2_minimizer(modelargs:Sequence,modelfunc:Callable,indep_vals:np.ndarray,dep_vals:np.ndarray):
    '''abstracted minimizer to be called from scipy.optimize.minimize (or from others?) arguments: modelargs: tuple/list of parameters that will be passed to the model modelfunc: the function to be called that model the data
                    -modelfunc requires signature (indepvals,arg1,arg2,arg3,etc)
        indep_vals: array of indep_vals
        dep_vals: array of dep_vals

    example usage:
        scipy.optimize.minimize(l2_minimizer,(modelfunc_guesses),args=(modelfunc,indep_values,dep_values))
    '''
    model_estimates=modelfunc(indep_vals,*modelargs)
    delta=model_estimates-dep_vals
    return np.dot(delta,delta)

def l1_minimizer(modelargs:Sequence,modelfunc:Callable,indep_vals:np.ndarray,dep_vals:np.ndarray):
    model_estimates=modelfunc(indep_vals,*modelargs)
#    return abs(np.sum(model_estimates-dep_vals))#delta,delta)
    normsum=np.sum([abs(x-y) for x,y in zip(model_estimates,dep_vals)])
    return normsum
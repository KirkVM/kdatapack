import numpy as np  
from dfitlib import fitmodels

#def l2_loss(tup, xs, ys):
#    k, m = tup
#    delta = model(xs, k, m) - ys
#    return numpy.dot(delta, delta)


def linearfit_l2(ptuple,xs, ys):
    intercept=ptuple[0]
    slopes=ptuple[1:]
    estimates=fitmodels.linear_model(xs,intercept,*slopes) 
    delta=estimates-ys
    return np.dot(delta,delta)



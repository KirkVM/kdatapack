import numpy as np  
from dfitlib import fitmodels

#def l2_loss(tup, xs, ys):
#    k, m = tup
#    delta = model(xs, k, m) - ys
#    return numpy.dot(delta, delta)


def linearfit_l2(xs, ys,intercept,*lin_coefficients):
    estimate=fitmodels.linear_model(xs,intercept,lin_coefficients) 
    delta=xs-ys
    return np.dot(delta,delta)


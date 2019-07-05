import numpy as np

#def linear_model(xs,intercept,*lin_coefficients):
#    '''returns linear function of x values
#    
#    Arguments:
#    xs: values of independent variable (shape=(x,ndimensions))
#    intercept: value of offset term
#    lin_coefficients: variable number of coefficients (num must match ndimensions)
#    '''
#    print(xs.shape)
#    print(len(xs.shape))
#    print(len(lin_coefficients))
#    if len(xs.shape)==1:
#        assert (len(lin_coefficients)==1)
#    else:
#        assert (xs.shape[0]==len(lin_coefficients))
#        if xs.shape[1]==len(lin_coefficients):
#            print('xs are expected to be in [X0s,X1s,X2s,...]. Bc dims are same length, please confirm this')
#    #return intercept+np.dot(lin_coefficients,xs)
#    #estimates=xs.copy()
#    estimates=np.array([xidx*xcoef for xidx,xcoef in zip(xs,lin_coefficients)]) 
#    return intercept
#    return np.dot(lin_coefficients,xs)


def linear_model(xvals,intercept,*slopes):
    '''returns linear function of x values
    '''
    rval=np.zeros(len(xvals))
    rval+=intercept
    if len(slopes)==1:
        slopes=slopes[0]
    else:
        slopes=np.array(slopes)
    newsum=np.dot(xvals,slopes)
    
    return rval+newsum
#    return slope*xs+intercept
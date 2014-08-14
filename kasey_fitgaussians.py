# -*- coding: utf-8 -*-
from scipy.optimize import leastsq
from pylab import exp,shape,size,zeros

def gaussians(params,x,numgaussians=1):
    """Defines multiple Gaussians with characteristics defined by 'params':
    params = yoffset,ymax,stddev,x0
    For multiple Gaussians, string together parameters into one long list (so we can fit using leastsq).
    Important: for multiple Gaussians, there is only one yoffset, and this will be the first element of the list;
    each individual Gaussian then has three parameters."""
    freeparams = 3 # number free parameters per gaussian, not counting yoffset
    if numgaussians==1:
        if size(params) != freeparams+1: raise NameError('Incorrect number of parameters supplied to function gaussians.')
        yoffset,ymax,stddev,x0=params
        return yoffset+ymax*exp(-((x-x0)/stddev)**2/2.0)
    else:
        if numgaussians != (size(params)-1)/freeparams: 
            raise NameError('Incorrect number of parameters supplied to function gaussians.')
        yoffset=params[0]
        total=zeros(len(x))
        for ii in range(numgaussians):
            ymax=params[freeparams*ii+1]
            stddev=params[freeparams*ii+2]
            x0=params[freeparams*ii+3]
            total = total + ymax*exp(-((x-x0)/stddev)**2/2.0)
        return total+yoffset
        
    
def fitgaussians(x, y, initparams, numgaussians=1):
    """Simultaneously fit multiple gaussians to a single data set (e.g. a spectrum).
    Returns bestfit, bestfitparameters."""
    def residuals(p, y, x):
        return y-gaussians(p,x,numgaussians)
    bestparams = leastsq(residuals, initparams, args=(y, x))[0]
    bestfit = gaussians(bestparams,x,numgaussians)
    return bestfit, bestparams
    
def splitparams(params,freeparams=3):
    """Say you've got a long parameter list specifying several Gaussians. This function splits it into
    a separate list for each Gaussian and then returns a list of those lists.
    yoffset (first element of parameter list) is not returned with any of the Gaussians."""
    splitapart = [ [0.0]*freeparams ]*((size(params)-1)/freeparams)
    for ii in range((size(params)-1)/freeparams):
        ymax=params[freeparams*ii+1]
        stddev=params[freeparams*ii+2]
        x0=params[freeparams*ii+3]
        splitapart[ii]=[ymax,stddev,x0]
    return splitapart
    
def absparams(params):
    # cant' just do abs(params) because it's a list and not an array...
    for ii in range(1,len(params)):
        params[ii] = abs(params[ii])
    return params
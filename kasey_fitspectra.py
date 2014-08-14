# -*- coding: utf-8 -*-
####
#### This is a newer version of kasey_fitgaussians.py, so use this.
#### This adds ability to fit Lorentzians as well as better parameter handling (dictionaries rather than lists).
####
####
from __future__ import division
from scipy.optimize import leastsq, curve_fit
import numpy as np
import matplotlib.pyplot as plt

def absparams(params):
    """can't just do abs(params) because it's a list and not an array..."""
    for ii in range(1,len(params)):
        params[ii] = np.abs(params[ii])

def blank_yoffset(params):
    """Set the yoffset of each of the functions to zero (usually done for plotting purposes)."""
    for peak in params:
        if type(peak)==dict: 
            peak['yoffset']=0.0
        else:
            peak[0]=0.0
            
def get_initparams_from_mouse( halfwidth=5.0 ):
    """
    Assuming you've already plotted the data you want to fit, this function makes a list 
    of initparams from your mouse clicks on the peaks you want to fit.
    """
    print 'Picking peaks only (not peaks and widths; see code)'
    peaks_or_peakswidths='p' #peaks_or_peakswidths = raw_input('Pick only peaks (p) or peaks and widths (pw)')
    if peaks_or_peakswidths=='p': 
        if True:
            # explitly use "connect" to get points
            peaks = [] # this will be a list of dictionaries of xloc and height of each peak
            def click(event):
                print 'you clicked', event.xdata
                peaks.append({'xloc':event.xdata, 'ymax':event.ydata})
            cid = plt.connect('button_press_event', click)
            raw_input('Click on the peaks you want to fit, then press Enter to continue.\n')
            plt.disconnect(cid)
        else:
            # use ginput
            print ('Click on the peaks you want to fit, then press Enter to continue.\n')
            points = plt.ginput(-1)
            peaks = [{'xloc':p[0], 'ymax':p[1]} for p in points]
    elif peaks_or_peakswidths=='pw':
        print 'This function is not yet implemented'
        peaks = []
    else:
        print 'Incorrect input.'
        peaks = []

    initparams=[0.0] # this is the initial guess at yoffset
    for peak in peaks:
        initparams.append(peak['ymax'])
        initparams.append( halfwidth )
        initparams.append(peak['xloc'])
    return initparams

def easyfitgaussians( datax, datay, halfwidth=5.0 ):
    """
    This function assumes that you have already plotted the data that you want
    to fit. You click on the peaks that you want to fit, then hit return and
    it tries to fit the data with that many gaussians.
    """
    initparams = get_initparams_from_mouse( halfwidth )
    ybest, pbest, std_err = fitgaussians(datax,datay,initparams)
    absparams( pbest )
    
    return ybest, params_to_dicts( splitparams(pbest), splitparams(std_err) )


def easyfitlorentzians( datax, datay, halfwidth=5.0 ):
    """
    This function assumes that you have already plotted the data that you want
    to fit. You click on the peaks that you want to fit, then hit return and
    it tries to fit the data with that many lorentzians.
    usage:
    ybest, bestparams = easyfitlorentzians( datax, datay, halfwidth=5.0 )
    where
    ybest is the line of best fit and bestparams is a dictionary of fit
    parameters and their 95% confidence interval (standard error).
    """
    initparams = get_initparams_from_mouse( halfwidth )
    ybest, pbest, std_err = fitlorentzians( datax, datay, initparams )
    absparams(pbest)
    dict_of_bestparams = params_to_dicts( splitparams(pbest), splitparams(std_err) )
    for peak in dict_of_bestparams:
        peak['Q'] = peak['x0']/2/peak['halfwidth']
        # error in Q is returned as a min/max tuple of the standard error
        peak['Q_err'] = np.asarray( (peak['x0']-peak['x0_err'])/2/(peak['halfwidth']+peak['halfwidth_err']),
                                 (peak['x0']+peak['x0_err'])/2/(peak['halfwidth']-peak['halfwidth_err']) )
        peak['Q_err'] = np.abs( peak['Q'] - peak['Q_err'] )
        
    return ybest, dict_of_bestparams


def fitgaussians(x, y, initparams):
    """Simultaneously fit multiple gaussians to a single data set (e.g. a spectrum).
    Returns bestfit, bestfitparameters."""
    if type(initparams) == dict:
        initparams = [initparams['yoffset'], initparams['ymax'],
                      initparams['halfwidth'], initparams['x0']]
        
    #def residuals(p, y, x):
    #    return y-gaussians(p,x)
    #bestparams = leastsq(residuals, initparams, args=(y, x))[0]


    bestparams, pcov = curve_fit( gaussians, x, y, p0=initparams )
    bestfit = gaussians( x, *bestparams )
    std_err = np.sqrt( np.diag(pcov) ) # is this true?

    return bestfit, bestparams, std_err
    
def fitlorentzians(x, y, initparams):
    """
    Simultaneously fit multiple lorentzians to a single data set (e.g. a spectrum).
    Returns bestfit, bestfitparameters,
    params = yoffset,ymax,halfwidth,x0
    """
    if type(initparams) == dict:
        initparams = [initparams['yoffset'], initparams['ymax'],
                      initparams['halfwidth'], initparams['x0']]
    elif type(initparams) == list and type(initparams[0]) == dict:
        initparams = params_to_lists(initparams)
        
    #def residuals(p, y, x):
    #    return y-lorentzians(p,x)
    #bestparams = leastsq(residuals, initparams, args=(y, x))[0]

    bestparams, pcov = curve_fit( lorentzians, x, y, p0=initparams )
    bestfit = lorentzians(x, *bestparams)
    std_err = np.sqrt( np.diag(pcov) ) # is this true?

    return bestfit, bestparams, std_err
    
def gaussians( x, *args, **kwargs ):
    """
    Defines multiple Gaussians with characteristics defined by 'params':
    params = yoffset,ymax,stddev,x0
    For multiple Gaussians, string together parameters into one long list (so we can fit using leastsq).
    Important: for multiple Gaussians, there is only one yoffset, and this will be the first 
    element of the list; each individual Gaussian then has three parameters.
    For a single Gaussian, params can be a dict containing the keys yoffset, ymax, x0, and halfwidth.
    """
    if len(kwargs) > 0:
    #if type(params)==dict:
        return kwargs['yoffset']+kwargs['ymax']/(1+((x-kwargs['x0'])/kwargs['halfwidth'])**2)
    else:
        freeparams = 3 # number free parameters per lorentzian, not counting yoffset
        if np.mod(len(args)-1,freeparams) != 0: raise NameError('Incorrect number of parameters supplied to function gaussians.')
        total=np.zeros(len(x))
        for ii in range(1,len(args),freeparams):
            ymax=args[ii]
            stddev=args[ii+1]
            x0=args[ii+2]
            total = total + ymax*np.exp(-((x-x0)/stddev)**2/2.0)
        yoffset=args[0]
        return yoffset+total
            
def lorentzians( x, *args, **kwargs ):
    """Defines multiple Lorentzians with characteristics defined by 'params':
    params = yoffset,ymax,halfwidth,x0
    For multiple Lorentzians, string together parameters into one long list (so we can fit using leastsq).
    Important: for multiple Lorentzians, there is only one yoffset, and this will be the first element of the list;
    each individual Lorentzian then has three parameters.
    For a single Lorentzian, params can be a dict containing the keys yoffset, ymax, x0, and halfwidth."""

    if len(kwargs) == 1:
        if 'params' in kwargs.keys():
            params = kwargs['params']
            return params['yoffset']+params['ymax']/(1+((x-params['x0'])/params['halfwidth'])**2)
        else:
            raise KeyError, "Single keyword argument must be 'params'."
    elif len(kwargs) > 1:
        return kwargs['yoffset']+kwargs['ymax']/(1+((x-kwargs['x0'])/kwargs['halfwidth'])**2)
    else:
        freeparams = 3 # number free parameters per lorentzian, not counting yoffset
        if np.mod(len(args)-1,freeparams) != 0: raise NameError('Incorrect number of parameters supplied to function lorentzians.')
        total=np.zeros(len(x))
        for ii in range(1,len(args),freeparams):
            ymax=args[ii]
            halfwidth=args[ii+1]
            x0=args[ii+2]
            total = total + ymax/(1+((x-x0)/halfwidth)**2)
        yoffset=args[0]
        return yoffset+total
            
def params_to_dicts( params, *args, **kwargs ):
    """It's a pain to remember which element in the parameter list corresponds
    to x0, yoffset, etc., so this function converts the list (or list of lists)
    to a list of dictionaries labeled with halfwidth (same as stddev for
    gaussians), height, x0, and yoffset (which is zero by default for
    multiple gaussians or lorentzians). This is meant to be called after the
    function splitparams."""
    
    std_err = None
    if 'std_err' in kwargs.keys():
        """ std_err from the fit has also been passed.
            like 'params', the errors should also be a list of lists
            (one for each peak)
        """
        std_err = kwargs['std_err']
        if size(std_err) != np.size(params):
            """ the length of the std_err list doesn't match the parameters,
                so we can't match them up properly
            """
            raise ValueError('the length of the list of std_err does not match length of list of parameters')
    elif len(args) > 0 and np.size(args[0]) == np.size(params):
        std_err = args[0]
    
    if len(np.shape(params)) == 1: 
        # make sure it's a list of lists, even for a single function
        params = list(params) 
        if std_err is not None: std_err = list(std_err)

    dicts = []
    for i,peak in enumerate(params):
        if std_err is not None:
            dicts.append({'yoffset':peak[0],   'yoffset_err':std_err[i][0],
                          'ymax':peak[1],      'ymax_err':std_err[i][1], 
                          'halfwidth':peak[2], 'halfwidth_err':std_err[i][2],
                          'x0':peak[3],        'x0_err':std_err[i][3]})
        else:
            dicts.append({'yoffset':peak[0], 'ymax':peak[1], 'halfwidth':peak[2], 'x0':peak[3]})
    
    return dicts
    
def params_to_lists(params):
    """Dictionaries are more convenient for storing and working with the parameters of the gaussians and lorentzians,
    but leastsq wants the initial parameters as a list (as far as I can tell...). This will take a list of dictionaries
    and convert it to a single list according to the following order:
    yoffset, ymax, halfwidth, x0 (repeating ymax, halfwidth, and x0 for any additional functions)"""
    if type(params) != list: raise TypeError('Incorrect data type: function params_to_list needs a list of dictionaries.')
    listofparams = [params[0]['yoffset']] # yoffset should be the same for all functions, so just pass it from the first one.
    for peak in params:
        listofparams.append(peak['ymax'])
        listofparams.append(peak['halfwidth'])
        listofparams.append(peak['x0'])
    return listofparams
    
def splitparams(params,freeparams=3):
    """Say you've got a long parameter list specifying several Gaussians. This function splits it into
    a separate list for each Gaussian and then returns a list of those lists.
    The same yoffset (first element of parameter list) is passed to all of the functions, so each list looks like:
    [yoffset,ymax,halfwidth,x0]."""
    splitapart=[]
    yoffset=params[0]
    for ii in range( int((np.size(params)-1)/freeparams) ):
        ymax=params[freeparams*ii+1]
        halfwidth=params[freeparams*ii+2]
        x0=params[freeparams*ii+3]
        splitapart.append([yoffset,ymax,halfwidth,x0])
    return splitapart
    
    

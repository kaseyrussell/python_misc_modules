# -*- coding: utf-8 -*-
import numpy as np
import glob
import os.path
import WinspecUtils as spe
import matplotlib.pyplot as plt


def background_and_spike_correct( fname_spectrum, fname_bg, threshold=3.0, filetype='txt', cts_per_sec=False ):
    if ( fname_spectrum.split('.')[-1] in ['spe', 'SPE'] ) or ( filetype == 'spe' ):
        wavelength, data = spe.getSpectrum( fname_spectrum, cts_per_sec )
        wavelength, bg = spe.getSpectrum( fname_bg, cts_per_sec )
        data -= bg
        diffx, diffy = despike( wavelength, data, threshold )
        return diffx, diffy
    else:
        data = loadtxt( fname_spectrum )
        bg = loadtxt( fname_bg )
        bgcorrected = data[:,2]-bg[:,2]
        diffx, diffy = despike( data[:,0], bgcorrected, threshold )
        return diffx, diffy
    
def centers_to_corners(x):
    """
    This function is for pcolor/pcolormesh plotting. For example, you take spectra over
    a grid of places on your sample, and you want to make a colormap of the intensity versus
    position. The positions you have are actually the CENTER points of the polygons that you
    want to make your colorplot out of, and functions like pcolormesh assume that you pass
    them the CORNERS of the polygons. (They thus want you to pass them one more row and column
    of location data than intensity/height data, and if they're the same size they'll drop
    the last row and column of intensity/height data.)
    
    This function takes a 1D ordered vector or list of locations and gives you back a 1D array
    of the locations in between the points you gave (with points outside the original end points):
    
    # i.e. we want the x points transformed to X as: [and note len(X)=len(x)+1]
    # X0 = -(x1-x0)/2
    # X1 = 0 +(x1-x0)/2
    # X2 = x1 + (x2-x1)/2
    # ...
    # X(N-1) = x(N-1) + (xN-x(N-1))/2
    # XN = xN + (xN-x(N-1))/2
    
    usage:
    X = centers_to_corners(x)
    """

    X = np.zeros(len(x)+1)

    if len(x) == 1:
        # we've only got one point along this dimension
        # (e.g. our raster scan is a line of points)
        # so we'll just make it 200nm wide by default so that
        # we have something to plot until we make it to the next line
        X[0] = x[0]
        X[1] = x[0] + 200.0 # units are in nm
        return X
    
    for i in range(len(x)+1):
        if 0<i<len(x):
            X[i] = x[i-1]+(x[i]-x[i-1])/2.0
        elif i==0:
            X[i] = x[i]-(x[i+1]-x[i])/2.0
        elif i==len(x):
            X[i] = x[i-1]+(x[i-1]-x[i-2])/2.0
        else:
            print ("how the hell did this happen?")
    return X
    
def despike(datax,datay,threshold=3.0,interpolate=False):
    """
    Remove spikes larger than threshold*(standard deviation of diffs in datay).
    interpolate=False (default) just removes the points;
    interpolate=True will just make them equal to the point to the left 
    (if you need to keep array length constant, say).
    """
    d = np.diff(datay)
    spikes = np.find(abs(d)>threshold*std(d)) 
    if (len(spikes)>0.1*len(d)): 
        print 'Did not remove spikes because it wanted to remove too many.'
        return datax,datay # if we're trying to remove a lot of points, just send it back...
    spikes=np.delete(spikes,find(diff(spikes)==1)+1) # if spike is one point, don't delete point after as well.
    if interpolate==False:
        datax = np.delete(datax,spikes+1)
        datay = np.delete(datay,spikes+1)
    else:
        # don't actually 'interpolate', just set it equal to the previous point
        # (actual interpolation could get messy in the case of a two-point spike, for example)
        for i in spikes+1: datay[i] = datay[i-3]
    if (np.any(np.abs(diff(datay))>threshold*np.std(d))): datax,datay = despike(datax,datay,threshold,interpolate)
    return datax, datay

def parse_fname_for_wirenum_and_angle(fname_parallel,fname_perpendicular):
    """
    Try to get the angle of the wire (returns angle parallel to axis) as well as the wire number:
    wire, angle = parse_fname_for_wirenum_and_angle(fname_parallel,fname_perpendicular)
    By using both the parallel and perpendicular filenames, I hope to make the angle algorithm more robust.
    """
    import string # has string constants
    delimiter = '_' # this is typically what I use...
    xpol=fname_perpendicular.split(delimiter)
    for ii,section in enumerate(fname_parallel.split(delimiter)):
        if section != xpol[ii]: break # first non-equal string segment is usually where I have angle info...
    wire=xpol[ii-1] # the wire info is usually the field just before the angle info.
    nums=[]
    for char in section:
        if char in string.digits: nums.append(char)
    if nums==[]: raise NameError('Function parse_fname_for_wirenum_and_angle could not successfully find angle.')
    return wire, float(''.join(nums)) # return the angle as a floating point number
    
def plot_raster( fname_base, interval='full', plot_relative_position=True, fignum=1, program='python', flipx=False, flipy=True, normalize=True ):
    """ 
    Plot a raster scan acquired using either the Excel program or (in the future) the python program.
    The core of this function is the same either way; the only reason for including the "program" flag
    is that I didn't want to be tied to the horiz, vert encoding process for the python program.
    """
    plt.figure( fignum )

    if '*.SPE' in fname_base:
        files = glob.glob( fname_base )
    else:
        files = glob.glob( fname_base + '*.SPE' ) # get a list of all spe files

    if len(files) == 0:
        # nothing to plot
        return None

    files.sort()
    # example file name:
    # PbS10_wire06_vert9083_horiz12796_1.txt

    # first just get the range of horizontal and vertical locations:
    vlist = []
    hlist = []
    if program == 'python':
        for fullfname in files:
            fname = str( os.path.splitext( os.path.basename( fullfname ) )[0] )
            horiz = float( fname.split('x')[1].split('_')[0] ) # this is an absolute position in nanometers
            hlist.append( horiz )
            vert = float( fname.split('y')[1] )
            vlist.append( vert )
    elif program == 'excel':
        raise ValueError( 'now only supporting python-generated data' )
        for fname in files:
            vert = float( fname.split('vert')[1].split('_')[0] ) # this is an absolute position in nanometers
            vlist.append( vert )
            horiz = float( fname.split('horiz')[1].split('.')[0] )
            hlist.append( horiz )

    vlist = list( set(vlist) ) # remove duplicate elements
    vlist.sort()
    vcen = centers_to_corners( vlist )

    hlist = list( set(hlist) )
    hlist.sort()
    hcen = centers_to_corners( hlist )

    if plot_relative_position:
        " plot relative distance, not absolute piezo position"
        vcen = vcen-vcen.min()
        hcen = hcen-hcen.min()

    # these step across in horizontal positions for each vertical position
    # (i.e. they go along a row before switching rows) dist. are in nm
    lum = np.zeros([ len(vlist), len(hlist) ]) # initialize matrix for holding luminescence data
    fname_matrix = [ [0]*len(hlist) for i in range(len(vlist)) ]  # initialize matrix for holding names of files
    for fullfname in files:
        fname = str( os.path.splitext( os.path.basename( fullfname ) )[0] )
        horiz = float( fname.split('x')[1].split('_')[0] ) # this is an absolute position in nanometers
        vert = float( fname.split('y')[1] )
        wavelen, d = spe.getSpectrum( fullfname )
        lum[ vlist.index(vert), hlist.index(horiz) ] = ( np.sum(d) if type(interval) is str else np.sum(d[interval]) )
        fname_matrix[vlist.index(vert)][hlist.index(horiz)] = fullfname
        
    xlabel('Position ($\mu m$)')
    ylabel('Position ($\mu m$)')
    if flipx:
        hpoints = (-hcen+hcen.max())/1000
    else:
        hpoints = hcen/1000

    if flipy:
        vpoints =( -vcen+vcen.max())/1000
    else:
        vpoints = vcen/1000


    if normalize:
        plot = pcolormesh( hpoints, vpoints, lum, cmap='copper' )
    else:
        from matplotlib.colors import NoNorm 
        plot = plt.pcolormesh( hpoints, vpoints, lum/lum.max(), cmap='copper', norm=NoNorm() )
        
    #pcolormesh( hcen/1000, vcen/1000, (lum-xpol_min)/(lum.max()-xpol_min), cmap='copper', norm=NoNorm() )

    plt.axis([ hpoints.min(), hpoints.max(), vpoints.min(), vpoints.max() ])

    ax = plt.gca()
    ax.set_aspect('equal')
    ax.axis('image') # not sure what this does, but it keeps ginput from changing the axes after a click...

    plt.show()
    return plot, (fname_matrix, vcen/1000, hcen/1000)
    
def process_polarized_data(fname_parallel, fname_perpendicular, fname_bg='none', 
        threshold=3.0, filetype='txt', cts_per_sec=False ):
    """
    If a background spectrum is supplied (for example, with glued spectra), this function first subtracts the
    background trace from both the parallel and perpendicular spectra, and then (regardless) subtracts the 
    perpendicular spectrum from the parallel one and despikes the result before returning it.
    It is set up to handle spe data (with flag filetype='spe') but wavelength info is not yet processed unless
    the dataset was glued."""
    if ( fname_parallel.split('.')[-1] in ['spe', 'SPE'] ) or ( filetype == 'spe' ):
        wavelength, parallel = spe.getSpectrum( fname_parallel, cts_per_sec )
        wavelength, perp = spe.getSpectrum( fname_perpendicular, cts_per_sec )
        ppol_cor = parallel
        xpol_cor = perp
        if fname_bg != 'none': 
            """ actually, I think this is pointless since we just subtract xpol from ppol anyway... """
            wavelength, bg = spe.getSpectrum( fname_bg, cts_per_sec )
            ppol_cor -= bg
            xpol_cor -= bg
        diffx, diffy = despike( wavelength, ppol_cor-xpol_cor, threshold )
        return diffx, diffy
    else:
        parallel = loadtxt(fname_parallel)
        perp = loadtxt(fname_perpendicular)
        wavelength = parallel[:,0]
        ppol_cor = parallel[:,2]
        xpol_cor = perp[:,2]
        if fname_bg != 'none': 
            bg = np.loadtxt(fname_bg)
            ppol_cor -= bg[:,2]
            xpol_cor -= bg[:,2]
        diffx, diffy = despike( wavelength, ppol_cor-xpol_cor, threshold )
        return diffx, diffy
        
def pt_to_inches( size ):
    """
    convert pt to inches: inches_per_pt = 1.0/72.27.
    takes a tuple, returns a tuple
    usage:
    (width, height) = pt_to_inches( (width, height) )
    """
    return ( size[0]/72.27, size[1]/72.27 )


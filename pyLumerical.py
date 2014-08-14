# -*- coding: utf-8 -*-
from __future__ import division
import kasey_utils as kc
import numpy as np
import matplotlib.pyplot as plt

"""
Various functions to parse Lumerical output using python.
KJR June 2010 (start)
"""

def loadmat_timeseries(fname, pos=True, keys=None):
    """
    Loads a 1D timeseries dataset that was exported from Lumerical in matlab format.
    We assume that the dataset contains the keys:
    t, Ex, Ey, Ez

    Function calculates fft, returns wavelength in nm, E2(wavelength)
    If pos==True (default), then return only the positive-frequency values.

    keys (optional) is a dict() with keys t, Ex, Ey, Ez for specifying the names
    of the corresponding keys in the matlab dataset.
    """
    from scipy.io import loadmat

    d = loadmat(fname)
    if keys is None:
        t_key  = 't'
        Ex_key = 'Ex'
        Ey_key = 'Ey'
        Ez_key = 'Ez'
    else:
        t_key  = keys['t']
        Ex_key = keys['Ex']
        Ey_key = keys['Ey']
        Ez_key = keys['Ez']

    assert d.has_key(t_key)
    assert d.has_key(Ex_key)
    assert d.has_key(Ey_key)
    assert d.has_key(Ez_key)
   
    t  = d[t_key]
    Ex = d[Ex_key][0,0,0]
    Ey = d[Ey_key][0,0,0]
    Ez = d[Ez_key][0,0,0]

    fs  = np.fft.fftfreq(t.size, d=(t[1]-t[0])[0]) #/2/np.pi
    fs  = np.fft.fftshift(fs)
    Exw = np.fft.fft(Ex)    # fft each component separately
    Eyw = np.fft.fft(Ey)
    Ezw = np.fft.fft(Ez)
    E2w = np.fft.fftshift(np.abs(Exw)**2+np.abs(Eyw)**2+np.abs(Ezw)**2)  # combine field components
    if pos:
        start = np.where(fs>0)[0][0]
        return 3.0e8/fs[start:]*1e9, E2w[start:]
    else:
        return 3.0e8/fs*1e9, E2w

def load_bandstructure( fname ):
    """
    This function loads in a *.txt file that was manually saved after running a bandstructure sweep.
    This will only work if you used the following formatting style to save the data to text file:

        fname = "bands.txt";
        write(fname,"#nx:"+num2str(length(kindices)));
        write(fname,"#ny:"+num2str(length(f)));

        write(fname,"## k-indices");
        write(fname,num2str(kindices));

        write(fname,"## f (Hz)");
        write(fname,num2str(f));

        write(fname,"## E-field intensity matrix");
        write(fname,num2str(transpose(fs_all)));
    
    That is, the first two lines need to specify the number of columns and rows of the data set;
    the order data is written is: column values (usually wavevector),
    row values (usually frequency), and amplitude matrix values (e.g. electric field amplitude);
    and each vector or matrix of data is preceded by a line that starts with '#'.
    """
    with open( fname ) as f:
        print "file:", fname
        
        """ first line will tell you number of columns. """
        num_columns = int(f.readline().split('\n')[0].split(':')[-1])
        print "Num. cols:", num_columns
        
        """ second line will tell you number of rows. """
        num_rows = int(f.readline().split('\n')[0].split(':')[-1])
        print "Num. rows:", num_rows
        
        """ next important line should start with '#', then we load x-data """
        line = f.readline()
        while not line.startswith('#'):
            f.readline()
        
        d1 = np.zeros( num_columns )
        for i in range(num_columns):
            d1[i] = float( f.readline().split('\n')[0] )
        
        """ next important line should start with '#', then we load y-data """
        line = f.readline()
        while not line.startswith('#'):
            f.readline()
        
        d2 = np.zeros( num_rows )
        for i in range(num_rows):
            d2[i] = float( f.readline().split('\n')[0] )
        
        """ next important line should start with '#', then we load matrix """
        line = f.readline()
        while not line.startswith('#'):
            f.readline()
        
        """ then we start with the matrix of field data: """
        fields = np.zeros( [num_columns, num_rows] )
        for i in range( num_columns ):
            line = f.readline().split('\n')[0]
            fields[i,:] = np.array( [ float(number) for number in line.split() ] )
        
    return d1, d2, fields


def load1D_txt( fname ):
    """
    Parse a txt exported from the Lumerical GUI. Works generally. Somehow np.loadtxt
    is unhappy with the txt files generated from the newest Lumerical version, so
    do this all by hand...
    """
    data1 = []
    data2 = []
    with open( fname, 'r' ) as fobj:
        """ First 3 lines are header."""
        print fobj.readline()
        print fobj.readline()
        print fobj.readline()
        for line in fobj:   
            d1, d2 = line.split(',')
            data1.append(float(d1))
            data2.append(float(d2))
        
    return np.array(data1), np.array(data2)

def load1D_timeseries( fname ):
    """
    UPDATE: I think this needs to be edited to work with the new Lumerical version.
    Say you exported data from a time monitor to a txt file and now want to load it in.
    This function gets rid of the header for you.
    usage:
    array([t, f]) = load1D_timeseries( filename )
    """
    with open( fname, 'r' ) as fileobject:
        line = fileobject.readline()
        # assume there's a header
        skiplines = 1
        for line in fileobject:
            if line[0] not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                skiplines += 1
            else:
                break

        return np.loadtxt( fname, delimiter=',', skiprows=skiplines,dtype='float' )

def load2D( fname ):
    """
    Say you used a frequency-domain field profile monitor and
    now you want a nice colorplot of your modes. Lumerical does
    this quickly, but it's hard to get it to fix the scale
    so that the width and length are proportional.
    
    This function assumes that you exported the data as a txt file
    straight from the CAD program.
    
        An aside:
        
        You can also apparently save the data to a .ldf file
        and then run lum2mat on it to convert it to a matlab file.
        (However, this didn't seem to work for me, and the image I plotted
        before exporting the data looked different than what the CAD program
        plotted automatically... ???)
        Like so: ::

            # get data from the simulation to be saved
            mname="xy_monitor";             # monitor name
            x=getdata(mname,"x");           # position vectors associated with Ex fields
            y=getdata(mname,"y");           # position vectors associated with Ex fields
            Ex=getdata(mname,"Ex");         # Ex fields at monitor
            f=getdata(m,"f");               # frequency
            fi=find(f,334.432e12);          # index of desired frequency (334.432 THz)
            onefreq=pinch(Ey,4,fi);         # select the desired frequency
            onefreq=pinch(onefreq);         # remove the singleton dimension (cuz it was 2D to start)
            
            # save variables x, y, Ex, T to a data file
            filename="results.ldf";  # set filename
            savedata(filename, x,y,onefreq);

            #convert the Lumerical Data File to a Matlab file
            lum2mat(filename);
    """
    with open( fname ) as f:
        print "file:", fname
        
        """ first line will tell you what monitor you're loading from, etc. """
        print f.readline().split('\n')[0]
        
        """ second line should be blank: """
        if f.readline() != '\n':
            raise("Shit. Problem loading second line of file within pyLumerical module.")

        """ third line will tell you what dimension the first array will represent (x or y or z)
            as well as its dimensions (nm vs um, etc.)."""
        print "Dimension 1 is:", f.readline().split('\n')[0]
        
        """ and fourth line tells you how big this first dimension is (array length) """
        d1_length = int(f.readline().split(',')[0].split('(')[1])
        print "length is:", d1_length
        d1 = np.zeros(d1_length)
        for i in range(d1_length):
            d1[i] = float( f.readline().split('\n')[0] )
        
        """ after loading the first dimension, it skips a line """
        f.readline()
        
        """ then we start with the second dimension: """
        print "Dimension 2 is:", f.readline().split('\n')[0]
        d2_length = int(f.readline().split(',')[0].split('(')[1])
        print "length is:", d2_length
        d2 = np.zeros(d2_length)
        for i in range(d2_length):
            d2[i] = float( f.readline().split('\n')[0] )

        """ after loading the second dimension, it skips a line """
        f.readline()
        
        """ then we start with the matrix of field data: """
        print "Loading matrix for:", f.readline().split('\n')[0]
        fields = np.zeros( [d1_length, d2_length] )
        for i in range( d1_length ):
            line = f.readline().split('\n')[0]
            if line[0] == ' ': line = line[1:] # cut off leading blank
            fields[i,:] = np.array( [ float(number) for number in line.split(' ') ] )
        
    return d1, d2, fields

def plot2Dmonitor( fname, **kwargs ):
    """
    make a nice pcolormesh plot that has same units in both axes
    usage: plot2Dmonitor( fname ) where fname is a txt file exported manually from CAD.
    Various keyword arguments are supported (in addition to any :
        rotate:    True/False for whether colormap should be rotated 90deg
        offset:    [x,y] tuple specifying lateral shift of origin for plotting
        ax:        axes instance for plotting. If not specified, gca() is used.
        take_sqrt: True/False. If True, take np.sqrt() before plotting (formerly E2)
        logscale:  True/False. If True, take np.log10() before plotting
        flipvert:  True/False. If True, multiply the locations of the vertical axis by -1 to flip it.
        
    """
    d1, d2, fields = load2D( fname )
    D1 = kc.centers_to_corners( d1 )
    D2 = kc.centers_to_corners( d2 )
    
    if 'E2' in kwargs.keys():
        raise ValueError, "E2 is deprecated and was replaced with take_sqrt (which acts equivalently)."

    if 'rotate' in kwargs.keys():
        rotate = kwargs['rotate']
        del kwargs['rotate']
    else:
        rotate = False
    
    if 'offset' in kwargs.keys():
        offset = kwargs['offset']
        del kwargs['offset']
    else:
        offset = [0,0]
    
    if 'ax' in kwargs.keys():
        ax = kwargs['ax']
        del kwargs['ax']
    else:
        ax=plt.gca()
        
    if 'flipvert' in kwargs.keys():
        flipvert = kwargs['flipvert']
        del kwargs['flipvert']
    else:
        flipvert=False
        
    if 'take_sqrt' in kwargs.keys():
        if kwargs['take_sqrt']:
            fields = np.sqrt(fields) # you exported E-squared (E intensity) but you want abs(E).
        del kwargs['take_sqrt']
        
    if 'logscale' in kwargs.keys():
        if kwargs['logscale']:
            fields = np.log10(fields)
        del kwargs['logscale']
        
    if rotate==True:
        if flipvert: D2 *= -1
        plot = ax.pcolormesh( D1+offset[0], D2+offset[1], np.transpose( fields ), **kwargs )
        ax.axis([ D1.min()+offset[0], D1.max()+offset[0], D2.min()+offset[1], D2.max()+offset[1] ])
    else:
        if flipvert: D1 *= -1
        plot = ax.pcolormesh( D2+offset[0], D1+offset[1], fields, **kwargs )
        ax.axis([ D2.min()+offset[0], D2.max()+offset[0], D1.min()+offset[1], D1.max()+offset[1] ])

    ax.set_aspect('equal')
    #show()

    return plot, D1, D2, fields

def plot_bandstructure( fname, rotate=True, cmap='jet', logscale=False, ax=None, E2=False ):
    """
    make a nice pcolormesh plot of the bandstructure
    usage:
    
        plot_bandstructure( fname )
        
    where fname is a txt file exported manually using a lumerical script. See
    'load_bandstructure' for proper formatting of this txt file.
    """
    d1, d2, fields = load_bandstructure( fname )
    D1 = kc.centers_to_corners( d1 )
    D2 = kc.centers_to_corners( d2 )
    
    if ax==None: ax=plt.gca()
    
    if logscale: fields = np.log(fields)

    if rotate==True:
        """ Note: the /1e12 converts Hz to THz """
        plot = ax.pcolormesh( D1, D2/1e12, np.transpose( fields ), cmap=cmap )
        ax.axis([ D1.min(), D1.max(), D2.min()/1e12, D2.max()/1e12 ])
    else:
        plot = ax.pcolormesh( D2, D1, fields, cmap=cmap )
        ax.axis([ D2.min(), D2.max(), D1.min(), D1.max() ])

    #show()

    return plot, D1, D2, fields



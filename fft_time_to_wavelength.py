# -*- coding: utf-8 -*-
from scipy import fft
from scipy.fftpack import fftfreq

# assume that we're trying to find the wavelength in air, which is what we
# would see with a spectrometer.
def fftwavelen( t, x, order=False ):
    n=1 # this will be length of vector that we fft
    while n<len(x): # fft works best on data sets of size that is a power of 2
        n = 2*n  
    n=n*8 # pad w/ extra zeros to inc. res. in freq. space
        # Lumerical seems to do a factor of 8 (at least in one instance that I checked)
    Xfft = abs(fft(x,n=n)[:n/2-1])**2 # "amplitude" according to Lumerical
    
    dt = t[2]-t[1] # assume you have at least 3 elements and all spacing uniform
    freq=fftfreq(n,dt)
    c=2.99792458e8 # speed of light in m/s
    l=c/freq[:n/2-1] # wavelength
    
    if order==True:
        return l[::-1], Xfft[::-1] # reverse arays so that they are in order of increasing wavelength
    else:
        return l, Xfft
    

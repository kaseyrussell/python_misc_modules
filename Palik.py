# I copied (by hand!) (the long-wavelength part of) the optical constants from Johnson and Christy
# PRB v6, pp 4370 (1972)

from constants import hbar, ee, me, A, e0, c
from pylab import *
import PalikTable
from scipy.interpolate import interp1d

data = PalikTable.table1

freq = data[:,0] # Hz
energy = freq * 2*pi*hbar
wavelength = c/freq # cm

Al_ep = data[:,1] + 1.0j * data[:,2]

interpolate_Al = interp1d( energy, Al_ep )
def ep_interp_Al( E ):
    """
    The data is kinda sparse, so interpolate between points using spline interpolation (via interp1d).
    This is just a wrapper around the function returned by interp1d, with bounds checking:
    it throws an error if you're outside the range of the dataset (no extrapolation, only interpolation)
    usage:
    ep = ep_interp_Al( E )
    where E is the energy (in eV), i.e. hbar times the radial frequency you're interested in
    """
    if any(E < energy[0]) or any(E > energy[-1]):
        raise ValueError( "Your energy is outside the range of the JC data set." )
    if size(E)==1: return interpolate_Al(E).astype(complex)
    return interpolate_Al( E )

def wavelen_to_energy( wavelen ):
    """ simply convert wavelength in nm to energy in eV
    """
    wavelen *= 1.0e-7 # convert nm to cm
    return 2*pi*c*hbar/wavelen

if __name__ == '__main__':

    if False:
        # dielectric constants versus energy
        plot( energy, Al_ep.real, 'r' )
        plot( energy, Al_ep.imag, 'b' )

    if False:
        # versus wavelength
        plot( wavelength*1e7, Al_ep.real, 'r-*' )
        plot( wavelength*1e7, Al_ep.imag, 'b-*' )

    if False:
        # versus wavelength, separate plots
        figure(1)
        plot( wavelength*1e7, Al_ep.real, 'r-*' )
        figure(2)
        plot( wavelength*1e7, Al_ep.imag, 'b-*' )

    if True:
        # versus wavelength, separate plots, with interpolated line
        wi = linspace(wavelength[1], wavelength[-2], 1000)*1e7
        epi = ep_interp_Al(wavelen_to_energy(wi.copy()))
        figure(1)
        plot( wi, epi.real, '-k' )
        plot( wavelength*1e7, Al_ep.real, 'ro' )
        figure(2)
        plot( wi, epi.imag, '-k' )
        plot( wavelength*1e7, Al_ep.imag, 'bo' )

    show()


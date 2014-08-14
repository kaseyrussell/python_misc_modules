# I copied (by hand!) (the long-wavelength part of) the optical constants from Johnson and Christy
# PRB v6, pp 4370 (1972)

from constants import hbar, ee, me, A, e0, c
from pylab import *
import JohnsonChristyTable
from scipy.interpolate import interp1d

#data = loadtxt( 'johnsonchristy.txt' )
data = array( JohnsonChristyTable.table1 )

energy = data[:,0]
freq = energy/hbar/2/pi # Hz
wavelength = c/freq # cm

Cu_n = data[:,1]
Cu_k = data[:,2]
Ag_n = data[:,3]
Ag_k = data[:,4]
Au_n = data[:,5]
Au_k = data[:,6]
dn = data[:,7]
dk = data[:,8]

# they also have a table of optical masses and scattering times
# (though the scattering times appear to be in inverse angular frequency, so I multiply by 2pi to get
# good fits between Drude and the tabulated data)
m_Cu = 1.49
tau_Cu = 6.9e-15*2*pi
Cu_ep = (Cu_n + 1.0j*Cu_k)**2 # the complex dielectric constant is simply the square of (n + ik)

m_Ag = 0.96
tau_Ag = 31.0e-15*2*pi
Ag_ep = (Ag_n + 1.0j*Ag_k)**2 # the complex dielectric constant is simply the square of (n + ik)

m_Au = 0.99
tau_Au = 9.3e-15*2*pi
Au_ep = (Au_n + 1.0j*Au_k)**2 # the complex dielectric constant is simply the square of (n + ik)



Ag_Density = 10.49 # g*cm^-3, mass density
Ag_AtomicMass = 107.8682 # g/mol, standard atomic weight
Ag_ElectronDensity = A*Ag_Density/Ag_AtomicMass*1e6 # number of atoms per cubic m., and one electron per atom is assumed
Ag_wp = sqrt(Ag_ElectronDensity * ee**2 / (m_Ag * me * e0))/2/pi # Hz (NOT RADIAL FREQUENCY!!!)

def ep_Ag_Drude( freq, temp=300.0 ):
    """
    usage:
    ep = ep_Ag_Drude( freq, temp=300.0 )
    where frequency freq is in Hz (not radial frequency)
    """
    if temp==300.0:
        return 1-Ag_wp**2/(freq**2 + 1j*freq/tau_Ag)
    elif temp==20.0:
        resistanceratio = 0.00257827
        return 1-Ag_wp**2/(freq**2 + 1j*freq*resistanceratio/tau_Ag)
    else:
        raise ValueError('Only temperatures of 300K (default) and 20K are supported.')
        

interpolate_Ag = interp1d( energy, Ag_ep )
def ep_interp_Ag( E ):
    """
    The JC data is kinda sparse, so interpolate between points using spline interpolation (via interp1d).
    This is just a wrapper around the function returned by interp1d, with bounds checking:
    it throws an error if you're outside the range of the dataset (no extrapolation, only interpolation)
    usage:
    ep = ep_interp_Ag( E )
    where E is the energy (in eV), i.e. hbar times the radial frequency you're interested in
    """
    if any(E < energy[0]) or any(E > energy[-1]):
        raise ValueError( "Your energy is outside the range of the JC data set." )
    if size(E)==1: return interpolate_Ag(E).astype(complex)
    return interpolate_Ag( E )

def wavelen_to_energy( wavelen ):
    """ simply convert wavelength in nm to energy in eV
    """
    wavelen *= 1.0e-7 # convert nm to cm
    return 2*pi*c*hbar/wavelen

if __name__ == '__main__':

    if False:
        # dielectric constants versus energy
        plot( energy, Ag_ep.real/30, 'r' )
        plot( energy, Ag_ep.imag, 'b' )

    if True:
        # versus wavelength
        plot( wavelength*1e7, Ag_ep.real, 'r-*' )
        plot( wavelength*1e7, Ag_ep.imag, 'b-*' )

    if False:
        # imaginary part divided by wavelength versus square of wavelength (see Fig. 6 in J&C)
        # this compares well, meaning I think I entered the data correctly...
        plot( wavelength**2 * 1e8, Ag_ep.imag/wavelength/1e4, 'b-*' )
        ax = gca()
        ax.set_ylim(0,16)
    
    if False:
        # plot the fit from the mass and scattering time given by JC:
        ep_Drude = ep_Ag_Drude( freq, temp=300.0 )

        fig=figure(1)
        fig.clf()
        ax1 = fig.add_subplot(111)
        ax1.set_title('Ag dielectric properties')
        ax1.plot( freq/1e12, Ag_ep.real, 'ro', label='Real' )
        ax1.plot( freq/1e12, ep_Drude.real, 'r-' )
        ax1.set_xlabel( 'Frequency (THz)' )
        ax1.set_ylabel( 'Real epsilon' )
        handles, labels = ax1.get_legend_handles_labels()
               
        ax2 = twinx()
        ax2.plot( freq/1e12, Ag_ep.imag, 'bo', label='Imag' )
        ax2.plot( freq/1e12, ep_Drude.imag, 'b-' )
        ax2.set_ylabel( 'Imaginary epsilon' )
        handles2, labels2 = ax2.get_legend_handles_labels()
        handles.append( handles2[0] )
        labels.append( labels2[0] )
        
        legend( handles, labels, loc=4 )
        
    show()


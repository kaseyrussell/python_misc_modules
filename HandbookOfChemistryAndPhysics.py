# here's some data a copied out of the handbook
from pylab import * # i do this so I can make it a numpy array

# temperature dependence of resistivity of metals:
# these are from on/about pp 12-40 (online 2010)
# here's the disclaimer from those tables:
##    The data refer to polycrystalline samples. The number of
##    significant figures indicates the accuracy of the values. However, at
##    low temperatures (especially below 50 K) the electrical resistivity
##    is extremely sensitive to sample purity. Thus the low-temperature
##    values refer to samples of specified purity and treatment. The references
##    should be consulted for further information on this point,
##    as well as for values at additional temperatures.

T = array([1, 10, 20, 40, 60, 80, 100, 150, 200, 273, 293, 298, 300], dtype=float) # (K)

rho_Au = array([
    0.0220,
    0.0226,
    0.035,
    0.141,
    0.308,
    0.481,
    0.650,
    1.061,
    1.462,
    2.051,
    2.214,
    2.255,
    2.271 ], dtype=float) * 1.0e-8 # (Ohm*m)

sigma_Au = 1.0/rho_Au
ratio_Au = rho_Au/rho_Au[-1] # normalize by 300K value

rho_Ag = array([ 
    0.00100,
    0.00115,
    0.0042,
    0.0539,
    0.162,
    0.289,
    0.418,
    0.726,
    1.029,
    1.467,
    1.587,
    1.617,
    1.629, ], dtype=float) * 1.0e-8 # (Ohm*m)

sigma_Ag = 1.0/rho_Ag
ratio_Ag = rho_Ag/rho_Ag[-1] # normalize by 300K value

rho_Al = array([
    0.000100,
    0.000193,
    0.000755,
    0.0181,
    0.0959,
    0.245,
    0.442,
    1.006,
    1.587,
    2.417,
    2.650,
    2.709,
    2.733], dtype=float) * 1.0e-8 # (Ohm*m)

sigma_Al = 1.0/rho_Al
ratio_Al = rho_Al/rho_Al[-1] # normalize by 300K value

if __name__ == '__main__':
    # plot the constants:
    figure(1)
    clf()
    semilogy( T, rho_Ag, 'bo-', label='Ag' )
    semilogy( T, rho_Au, 'ro-', label='Au' )
    semilogy( T, rho_Al, 'go-', label='Al' )
    ax = gca()
    xlabel( "Temperature (K)" )
    ylabel( "Resistivity (Ohm m)" )
    legend( loc=4 )

    figure(2)
    clf()
    semilogy( T, ratio_Ag, 'bo-', label='Ag' )
    semilogy( T, ratio_Au, 'ro-', label='Au' )
    semilogy( T, ratio_Al, 'go-', label='Al' )
    ax = gca()
    xlabel( "Temperature (K)" )
    ylabel( "Resistivity Relative to 300K Value" )
    legend( loc=4 )

    show()


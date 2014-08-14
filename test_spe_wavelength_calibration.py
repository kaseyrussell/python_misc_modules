# -*- coding: utf-8 -*-
from pylab import *
import piUtils as spe


if False:
    #test glued files
    datapath = '/home/kc/lab/DATA/2010/2010_feb25/'
    txtfile = 'glued_wire25_ppol353deg_1.txt'
    dtxt = loadtxt( datapath + txtfile )
    plot( dtxt[:,0], dtxt[:,2], 'b' )

    spefile = 'glued_wire25_ppol353deg.SPE'
    wavelength, lum = spe.getSpectrum (datapath + spefile )
    plot( wavelength, lum, 'g' )

if True:
    clf()
    # test non-glued files (harder!)
    datapath = '/home/kc/lab/DATA/2010/2010_feb25/'
    txtfile = 'wire27_crazystep_40_1.txt'
    dtxt = loadtxt( datapath + txtfile )
    wtxt = dtxt[:,0] 
    #plot( wtxt, dtxt[:,2], 'b' )

    spefile = 'wire27_crazystep_40.SPE'
    spe_struct = spe.readSpe( datapath + spefile )
    luminescence = transpose( spe_struct['data'][0] )
    xcal = spe_struct['XCALIB']
    wavelength, lum = spe.getSpectrum( datapath + spefile )
    #plot( wavelength, lum, 'g' )

    # these two are identical to less than 10^-3 nm.
    # In fact, wtxt == wavelength.round(3), so using the SPE file
    # is actually more accurate (although our spectrometer calibration
    # isn't anywhere close to accurate enough to make this matter).

    datapath = '/home/kc/lab/DATA/2010/2010_jan17/'
    txtfile = 'AgNWonAg_1sec633nm_1.txt'
    dtxt = loadtxt( datapath + txtfile )
    wtxt = dtxt[:,0] 
    plot( wtxt, dtxt[:,2], 'b' )

    spefile = 'AgNWonAg_1sec633nm.SPE'
    spe_struct = spe.readSpe( datapath + spefile )
    luminescence = transpose( spe_struct['data'][0] )
    xcal = spe_struct['XCALIB']
    wavelength, lum = spe.getSpectrum( datapath + spefile )
    plot( wavelength, lum, 'g' )


show()
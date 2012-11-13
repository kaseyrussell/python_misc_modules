
# routines for fitting my capacitance-voltage data. 
# I wrote this when trying to put together the third(!) revision of the PRB
# KJR 2010

from pylab import *
from scipy.optimize import leastsq
from scipy.optimize import fmin as fmin_plain
from scipy.optimize import fmin_bfgs
from scipy.optimize import curve_fit
import cPickle
import kasey_fitspectra as kcfit

datadir = '/home/kc/lab/DATA'
diameter=200.0e-6     # diameter in m
A = pi*(diameter/2)**2 # area in sq. m
Cbp = 11.0            # bonding pad capacitance in pF
Cgeom=80.0            # appx Ctj, pF
scaleC = 1.0
q = 1.60217646e-19

### here are functions to fit the TDOS (from Cq) versus energy to a set of Lorentzians.
###
def Cq_to_DOS( Cq ):
    return Cq/A/q

def LorentziansDOS( E, gd, halfwidth=5.0 ):
    """ Fit the DOS to a sum of Lorentzians.
        Usage:
        fit, params = LorentziansDOS( E, gd, halfwidth=5.0 )
        where:
        E is the Fermi energy
        fit is the best fit to gd
        params is a dict containing the parameters that gave the best fit
    """
    fit, params = kcfit.easyfitlorentzians( E, gd, halfwidth )
    return fit, params

### here are a few routines for fitting to the complex impedance
###
def Zdevice( p, pfixed, f ):
    """ equivalent circuit of device impedance:
        Z = Zdevice( p, pfixed, f )
        where
        pfixed = [Cg, Cbp, Rs, Rp]
        p = [Ctj, Cq, tau]
    """   
    # use petroff's treatment of the equivalent circuit
    if len(pfixed) == 4:
        Cg, Cbp, Rs, Rp = pfixed
        Ctj, Cq, tau = p
    elif len(pfixed) == 5:
        Ctj, Cg, Cbp, Rs, Rp = pfixed
        Cq, tau = p
    
    w = 2*pi*f
    s = 1e12 # a scaling factor to convert pF to F.  we want all elements of p to have similar order of mag.
    Zqstar = 1/(1j*w*Cq/s/(1+1j*w*tau))
    Ztj = 1/(1j*w*Ctj/s)
    Zg = 1/(1j*w*Cg/s)
    Zideal = Zg + Ztj*Zqstar/(Ztj+Zqstar) # impedance of ideal device (no Rs, Rp or Cbp)
    return Rs + 1/(1/Zideal + 1/(1e9*Rp) + 1j*w*Cbp/s) # impedance of total device, with Rs feeding Cbp and device
#    return 1/(1j*w*Cbp/s + 1/(Rs + 1/(1/Zideal + 1/(1e9*Rp))) )# impedance of total device, with Cbp in parallel with all

def Zfit_to_G( tau, Cq ):
    """ Find the tunneling conductance from MY equivalent circuit model (not Ashoori)
        Usage:
        G = Zfit_to_G( tau, Cq )
    """
    return Cq/tau*1e-12

def Zidealize( Z, Rs, Rp, Cbp, f ):
    """ "clean" the measured complex impedance Z to remove the effects of Rs, Rp, and Cbp.
        Usage:
        Zideal = Zidealize( Z, Rs, Rp, Cbp, f )
        where Rs is in Ohms, Rp in GOhms, and Cbp in pF
    """
    w = 2*pi*f
    s = 1e12 # a scaling factor to convert pF to F.  we want all elements of p to have similar order of mag.
    return 1/(1/(Z-Rs) - 1/(1e9*Rp) - 1j*w*Cbp/s)

def sumresidualsZ( p, pfixed, data, f ):
    """ calculate residuals and sum to package into a single value for optimizers like fmin.
    """
    d = data - Zdevice( p, pfixed, f )
    d *= f # weighting; will get squared in next step
    return sum((d*conj(d)).real)

def fitZ( f, Zdata, p0, pfixed ):
    """ find the capacitance and loss of a set of data (in f,Z form, where Z is complex)
        and fit this with Ashoori's model to find fpeak, Clow, Chigh, and the conductivity G
        Ctj, Cq, tau = fitZ( f, Zdata, p0, pfixed ) 
        where
        p0 = [Ctj0, Cq0, tau0]     # initial guesses
        pfixed = [Cg, Cbp, Rs, Rp] # fixed parameters
    """
    lsq = fmin_plain( sumresidualsZ, p0, args=(pfixed, Zdata, f), ftol=1e-8, xtol=1e-8, maxfun=2000, disp=False )
    if len(pfixed) == 4:
        Ctj, Cq, tau = abs(lsq)
        return Ctj, Cq, tau
    elif len(pfixed) == 5:
        Cq, tau = abs(lsq)
        return Cq, tau
    else:
        print "shit! wrong number fixed params to fitZ"
        return None
    

def findCg_SampleS( Vdc ):
    #   THIS IS NOT Cq!!! we're trying to fit holding Cg constant, but its value varies with Vdc.
    #   if we hold Cq constant, we can see that Cg looks like it varies like 1/(1+sqrt(V0-V)), as it should by Poisson.
    #   see pp. 90 of NB 6 for more info.
    #   this formula is from a fit in the pyplot_outp.py script in this directory (apr0208 python_corrected)
    e0 = 8.854e-12
    ep = 12*e0
    q=1.60217e-19
    Nd = 3e24    # 3D contact doping [m**3]
    V0 = 0.3     # bias where the 3D contact adjacent to the gate barrier starts to deplete
    wg = 38.0e-9 # gate barrier thickness w/o depletion [m] (include 1/2 of QW?)
    x=1/sqrt((V0-Vdc)*2*ep/q/Nd);
    p0 = 3.4021218081662091e-09
    p1 = -1.291285781497841e-17
    y = p0 + p1*x
    return A*ep/(y+wg)*1e12 # return Cg in pF

def findCg_SampleNU( Vdc ):
    #   THIS IS NOT Cq!!! we're trying to fit holding Cg constant, but its value varies with Vdc.
    #   if we hold Cq constant, we can see that Cg looks like it varies like 1/(1+sqrt(V0-V)), as it should by Poisson.
    #   see pp. 90 of NB 6 for more info.
    #   I stole these params from the 2009/feb1009 dir, which may or may not be right...
    e0 = 8.854e-12
    ep = 12*e0
    q=1.60217e-19
    Nd = 3e24    # 3D contact doping [m**3]
    V0 = 0.2     # bias where the 3D contact adjacent to the gate barrier starts to deplete
    wg = 35.0e-9 # gate barrier thickness w/o depletion [m] (include 1/2 of QW?)
    x=1/sqrt((V0-Vdc)*2*ep/q/Nd);
    p0 = 3.4021218081662091e-09
    p1 = -1.291285781497841e-17
    y = p0 + p1*x
    return A*ep/(y+wg)*1e12 # return Cg in pF

### Following few are routines for fitting capacitance and loss and finding conductance (as in Ashoori)
### I'm not sure this is a good way to do it. maybe better to fit to complex Z and then cast to C,D?
###
def ZtoCD( f, Z ):
    """ returns C in pF and Dielectric loss tangent, D
    """
    Y = 1/Z # admittance
    C = 1.0/f/2.0/pi * abs(Y.imag) * 1.0e12
    
    D = Y.real/abs(Y.imag)  # Loss
    return C, D

    
def CDtoZ( f, C, D ):
    """ Input C (pF) and D (unitless loss tangent) and the frequency f and get back the complex impedance Z:
        Usage:
        Z = CDtoZ( f, C, D )
    """
    B = 2*pi*f*C*1e-12
    G = D*B
    Y = G + 1j*B
    return 1/Y


def CD( p, f ):
    """ calculate capacitance and loss:
        C, D = CD( p, f )
        where 
        p = [fpeak, Clow, Chigh]
    """
    fpeak = p[0]
    Clow = p[1]
    Chigh = p[2]
    C = Chigh*Clow*(1+(f/fpeak)**2)/(Chigh + Clow*(f/fpeak)**2)
    D = 2*sqrt(Chigh/Clow)*(Clow/Chigh-1)*f/fpeak/(1+(f/fpeak)**2)
#    D = sqrt(Chigh/Clow)*(Clow/Chigh-1)*f/fpeak/(1+(f/fpeak)**2) # drop the two...
    return C, D


def sumresidualsCD( p, data, f ):
    """ calculate residuals and sum to package into a single value for optimizers like fmin.
    """
    C, D = CD( p, f )
    
    d = data - CDtoZ( f, C, D )
    d *= f**1.5 # weighting; will get squared in next step
    return sum( (d*conj(d)).real )
#    Z = CDtoZ( f, C, D )
#    d = abs( log(data.real)-log(Z.real) ) + abs( log(abs(data.imag))-log(abs(Z.imag)) )
#    d *= f # weighting; will get squared in next step
#    return sum( (d*conj(d)).real )


def fitCD( f, Z, p0, Cgeom=Cgeom ):
    """ find the capacitance and loss of a set of data (in f,Z form, where Z is complex)
        and fit this with Ashoori's model to find fpeak, Clow, Chigh, and the conductivity G
        fpeak, Clow, Chigh, G = fitCD( f, Z, p0, Cgeom=80.0 ) 
        where
        p0 = [fpeak0, Clow0, Chigh0] # initial guesses
    """
    lsq = fmin_plain( sumresidualsCD, p0, args=(Z, f), ftol=1e-8, xtol=1e-8, maxfun=2000, disp=False )
    fpeak, Clow, Chigh = abs(lsq)
        
    G = 2*pi*fpeak*Cgeom**2/sqrt(Clow*Chigh)*(Clow/Chigh - 1) * 1e-12   # see Eq. 2 of Ashoori's PRB
    return fpeak, Clow, Chigh, G
    
### here are routines for data correction:
###
def correct_CpGdata_forCbp( CpG ):
    CpG[:,1] -= Cbp*1e-12 # "parallel" capacitance
    return CpG

def convert_CpG_to_Zcomplex(CpG):
    f = CpG[:,0] # frequency
    C = CpG[:,1] # "parallel" capacitance
    G = CpG[:,2] # conductance
    B = 2*pi*f*C # susceptance
    
    Z = 1/sqrt(G**2 + B**2) # amplitude of impedance vector
    theta = -arctan(abs(B)/G) # see pp. 6-4 of the Agilent 4284A LCR manual

    # check angle on theta since arctan will only return a value between -pi/2 and pi/2:
    for i in range(Z.shape[0]):
        if G[i] < 0:
            theta[i] = theta[i] - pi
        elif B[i] < 0:
            theta[i] = theta[i] - 2*pi
        
    Zcomplex = Z*(cos(theta) + 1j*sin(theta))
    return Zcomplex  

def full_correction(devmeas,opencircuit,refmeas,refideal):#Torrents & Pallas-Areny, IEEE Trans Instr & Meas v51 pp310 (2002).
    f=devmeas[:,0]
    Imp_xm = convert_CpG_to_Zcomplex( devmeas )       # device impedance that we measured
    Imp_om = convert_CpG_to_Zcomplex( opencircuit )   # impedance of an open measurement in the setup
    Imp_sm = 0                                      # impedance of a short measurement in the setup (assumed zero here)
    Imp_ref = convert_CpG_to_Zcomplex( refideal )     # impedance of ideal reference
    Imp_rm = convert_CpG_to_Zcomplex( refmeas )       # impedance of ideal reference measured in the setup

    return f, (Imp_xm - Imp_sm) * Imp_ref/(Imp_rm - Imp_sm) * (Imp_om - Imp_rm)/(Imp_om - Imp_xm)

def open_correction(device,opencircuit):
    f=device[:,0]

    # in terms of admittances, multiply by 2*pi*f, as is commented out
    # but we can skip that since we just divide by it one step later,
    # and multiplication is associative.
    Bo = opencircuit[:,1]
    Bm = device[:,1]
    Bdut = Bm-Bo
    Cp = Bdut     

    Go = opencircuit[:,2]
    Gm = device[:,2]
    G = Gm-Go
    
    output = []
    output = np.array(np.zeros(3*G.shape[0])).reshape(G.shape[0],3) # all this just to add columns to a matrix? Jeez...
    output[:,0] = f
    output[:,1] = Cp
    output[:,2] = G
    return output

def loadpickle( pickle ):
    ifile = open( pickle,'rb' )
    d = cPickle.load( ifile )
    ifile.close()
    return sorted(d.iteritems(), key=lambda (k,v): int(k.split('_')[1])) # sort the dict into an iterable list


# Cython module for faster fitting of fluorescence lifetime spectra
# KJR Aug 2011

## this works but is only about 2x faster than pure python
## I guess because many of the expensive parts are either
## already in C or, like cspline1d_eval, not optimized with Cython.

import pylab
import numpy as np
cimport numpy as np
from cpython cimport bool
import cython
from scipy.signal import cspline1d, cspline1d_eval
from scipy.optimize import leastsq

cdef extern from "math.h":
    double exp(double)


cdef extern from "stdarg.h":
    ctypedef struct va_list:
        pass
    ctypedef struct nparray_type:
        pass
    void va_start(va_list, ...)
    void* va_arg(va_list, ...) 
    void va_end(va_list)

cdef extern from *:
    double double_type "double"

cdef class Function:
    cdef np.ndarray params, free_params, irf_generator, stddev
    cdef Py_ssize_t i, arg_len
    cdef double tpulse, irf_dt, irf_t0
    cdef bool has_key_t_ag, has_key_t_d3, has_key_a_fix
    cdef bool has_key_trise, has_key_tshift, has_key_l1
    cdef char* model
    def __init__( self,
            char* model,
            np.ndarray[np.float64_t, ndim=1] params,
            np.ndarray[np.float64_t, ndim=1] free_params,
            guess,
            double tpulse,
            double irf_dt,
            double irf_t0,
            np.ndarray[np.float64_t, ndim=1] irf_generator,
            np.ndarray[np.float64_t, ndim=1] stddev ):
        self.model = model
        self.params = params
        self.free_params = free_params
        self.arg_len = free_params.shape[0]
        self.has_key_t_ag = True if guess.has_key('t_ag') else False
        self.has_key_t_d3 = True if guess.has_key('t_d3') else False
        self.has_key_a_fix = True if guess.has_key('a_fix') else False
        self.has_key_tshift = True if guess.has_key('tshift') else False
        self.has_key_trise = True if guess.has_key('trise') else False
        self.has_key_l1 = True if guess.has_key('l1') else False
        self.tpulse = tpulse
        self.irf_generator = irf_generator
        self.irf_dt=irf_dt
        self.irf_t0=irf_t0
        self.stddev = stddev # used for weights during fitting
        
    cpdef double set_params(self, np.ndarray params) except *:
        self.params = params
    
    cpdef double set_free_params(self, np.ndarray[np.float64_t, ndim=1] free_params) except *:
        self.free_params = free_params
    
    cpdef print_params(self):
        print self.params
        
        

cdef class Convolved( Function ):
    def __init__( self, *args ):
        Function.__init__( self, *args )
        
    cdef np.ndarray[np.float64_t, ndim=1] multi_exponential( self, np.ndarray[np.float64_t, ndim=1] p, np.ndarray[np.float64_t, ndim=1] t ):
        """Typical multi-exponential fit. This can accomodate rising exponential (for saturation fitting)."""
        cdef Py_ssize_t i, arg_len, deduct
        cdef double arg, tshift, trise, t_ag, t_d3, scale, b, a, l
        cdef np.ndarray local_params, ideal, irf
        deduct = 1
        for i,arg in enumerate(p): self.params[ self.free_params[i] ] = arg
        local_params = self.params[:]
        
        if self.has_key_trise:
            trise = local_params[-deduct]
            deduct += 1
            
        if self.has_key_tshift:
            tshift = local_params[-deduct]
            deduct += 1
        else:
            tshift = 0.0

        if self.has_key_t_ag:
            t_ag = abs(local_params[-deduct])
        elif self.has_key_t_d3:
            t_d3 = abs(local_params[-deduct])
        elif self.has_key_a_fix:
            scale = local_params[-deduct]
        else:
            b = local_params[-deduct]
            
        ideal = np.zeros(len(t), dtype=np.float)
        for l,a in zip(local_params[:-2:2],local_params[1:-2:2]):
            if self.has_key_t_ag: l = 1.0/(1.0/l + 1.0/t_ag)
            if self.has_key_t_d3: l *= t_d3
            if self.has_key_a_fix: a *= scale
            ideal += abs(a)*np.exp(-t/abs(l))/(1.0-np.exp(-self.tpulse/abs(l)))
        if self.has_key_trise: ideal *= 1.0-np.exp(-t/abs(trise))
        irf = cspline1d_eval( self.irf_generator, t-tshift, dx=self.irf_dt, x0=self.irf_t0 )
        return np.real(np.fft.ifft( np.fft.fft(ideal)*np.fft.fft(irf) ))
        

    cdef np.ndarray[np.float64_t, ndim=1] stretched_exponential( self, np.ndarray[np.float64_t, ndim=1] p, np.ndarray[np.float64_t, ndim=1] t ):
        """Stretched single-exponential fit.
        See http://www.ncbi.nlm.nih.gov.ezp-prod1.hul.harvard.edu/pmc/articles/PMC1301608/pdf/11509343.pdf"""
        cdef Py_ssize_t i, arg_len, deduct
        cdef double arg, tshift, h, a, l, a1, h1, l1
        cdef np.ndarray local_params, ideal, irf
        cdef int N=10 # number previous pulses to include
        cdef int j
        deduct = 1
        for i,arg in enumerate(p): self.params[ self.free_params[i] ] = arg
        local_params = self.params[:]
        
        if self.has_key_tshift:
            tshift = local_params[-deduct]
            deduct += 1
        else:
            tshift = 0.0
        l, a, h = local_params[:3]
        a = abs(a)
        l = abs(l)
        if self.has_key_l1:
            l1, a1, h1 = local_params[3:6]
            a1 = abs(a1)
            l1 = abs(l1)
        ideal = np.zeros(len(t), dtype=np.float)
        for j from 0 <= j < N:
            ideal += a*np.exp(-((t+j*self.tpulse)/l)**(1.0/h)) # Kohlrausch function
            #ideal += a*np.exp(1-(1+(t+j*self.tpulse)/l)**(1.0/h))  # modified Kohlrausch function (see Berberan-Santos et al., 2005)
            if self.has_key_l1:
                ideal += a1*np.exp(-((t+j*self.tpulse)/l1)**(1.0/h1))
        irf = cspline1d_eval( self.irf_generator, t-tshift, dx=self.irf_dt, x0=self.irf_t0 )
        return np.real(np.fft.ifft( np.fft.fft(ideal)*np.fft.fft(irf) ))
        

    cpdef np.ndarray[np.float64_t, ndim=1] residuals(self,
            np.ndarray[np.float64_t, ndim=1] p,
            np.ndarray[np.float64_t, ndim=1] y,
            np.ndarray[np.float64_t, ndim=1] t):
        cdef np.ndarray res
        if self.model.decode() == u'multi_exp':
            res = (y-self.multi_exponential(p,np.array(t,dtype=np.float)))/self.stddev
        elif self.model.decode() == u'stretched_exp':
            res = (y-self.stretched_exponential(p,np.array(t,dtype=np.float)))/self.stddev
        
        return res

    def fit( self, np.ndarray[np.float64_t, ndim=1] t,
            np.ndarray[np.float64_t, ndim=1] data,
            np.ndarray[np.float64_t, ndim=1] initparams ):
        return leastsq( self.residuals, initparams, args=( data, t ), full_output=1 )
        # leastsq returns: (popt, pcov, infodict, errmsg, ier) = res
        

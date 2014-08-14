#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Just double-click a file of type .txt and this will try to plot its
    third column versus its first (like how our spectroscopy datasets work)"""
    
from pylab import *
import sys
from kasey_utils import despike

for arg in sys.argv[1:]:
    d = loadtxt(arg)
    x,y = despike(d[:,0],d[:,2],threshold=5.5) # use a lax threshold; better to miss a spike than cut data
    plot(x,y,label=arg.split('/')[-1])

#leg=legend(loc=6)
show()

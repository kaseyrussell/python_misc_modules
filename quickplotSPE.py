#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Just double-click a file of type .SPE and this will try to plot it
"""
    
import sys
import winspec
from pylab import figure, show

figure(1)

for arg in sys.argv[1:]:
    s = winspec.Spectrum( arg )
    s.plot()
    #ax = gca()
    #text(0.1,0.5, str(arg), transform = ax.transAxes)

show()


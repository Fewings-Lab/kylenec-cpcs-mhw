# -*- coding: utf-8 -*-
# filter.py
# Python 3.7
"""
@author: %(username)s
Created on %(date)s
Modified %(dates)s

Description:

"""

import numpy as np
# import xarray as xr
#import matplotlib.pyplot as plt

def lowpass(y,freq, win):
    
    # Get filter weights
    w = lanczos(freq,win)
    
    if sum(np.isnan(y)) == 0:
        ylow = np.convolve(y,w,mode='same')
    else:
        YLOW = np.convolve(y[np.isfinite(y)],w,mode='same')
        ylow = np.nan*np.ones(1,len(y))
        ylow[np.isfinite(y)] = YLOW
    
    # Throwing out ends of time series that show edge effects
    ylow[0:(win-1)/2] = np.nan*np.ones(1,win/2)
    ylow[-(win-1)/2:-1] = np.nan*np.ones(1,win/2)
    
    # High-passed is left behind after removing low-passed from raw
    yhigh = y - ylow
    
    return ylow, yhigh

def lanczos(freq,win):
    # freq is related by 2*pi/period
    # win is window length
    
    sig = 2*np.pi/(win)
    C = freq/np.pi
    t = np.arange(1,(win/2)+1)
    
    halfwin = C*np.sinc(freq*t/np.pi)*np.sinc(sig*t/np.pi)
    fullwin = np.concatenate((np.flip(halfwin),C,halfwin))
    w = fullwin/np.nansum(fullwin)    
    
    return w





#!/usr/bin/env python

# pyPL66.py
# 24 Aug 2020
# Kylene Cooley
# pl66tn.m but make it python
# will need to add appropriate credits to original filter and function

# import libraries
import numpy as np
import xarray as xr
import pandas as pd
import scipy.signal as sig
import matplotlib.pyplot as plt

def pypl66(x,dt=1,T=33):
    cutoff = T/dt
    fq = 1/cutoff
    nw = int(np.round((2*T)/dt))
    print("number of weights = "+str(nw))
    nw2 = 2*nw

    [npts,ncol] = x.shape
    if npts < ncol:
        x = x.T
        [npts,ncol] = x.shape
    xf = x

    # filter weights
    j = np.arange(1,nw+1)
    t = np.pi*j
    den = (fq**2)*(t**3)
    wts = (2*np.sin(2*fq*t)-np.sin(fq*t)-np.sin(3*fq*t))/den

    # make weights symmetric
    wts = np.concatenate((np.flip(wts),np.array([2*fq]),wts))
    #a = np.nansum(wts)
    wts /= np.sum(wts)

    # fold tapered time series
    cs = np.cos(t.T/nw2)
    jm = np.flip(j)
    
    for ic in range(ncol) :
        ind = np.invert(np.isnan(x))
        if ncol == 1:
            npts = len(x[ind])
        else:
            npts = len(x[ind,ic])
        if npts > nw2:
            if ncol == 1:
                xdt = sig.detrend(x[ind])
                trnd = x[ind]-xdt
                y = np.concatenate((cs[jm-1]*xdt[jm-1],xdt,cs[j-1]*xdt[npts-j]))
    
                # filter
                yf = sig.lfilter(wts, 1.0, y)
    
                # strip off extra points
                xf[ind] = yf[nw2:npts+nw2]
    
                # add back trend
                xf[ind] += trnd
            else:
            # detrend time series, then add back trend after filtering
                xdt = sig.detrend(x[ind,ic])
                trnd = x[ind,ic]-xdt
                y = [[cs[jm]*xdt[jm]],[xdt],[cs[j]*xdt[npts+1-j]]]
    
                # filter
                yf = sig.lfilter(wts, 1.0, y)
    
                # strip off extra points
                xf[ind,ic] = yf[nw2+1:npts+nw2]
    
                # add back trend
                xf[ind,ic] += trnd
        else:
            print("warning: time series is too short. column num = "+ str(ic+1))

    return xf

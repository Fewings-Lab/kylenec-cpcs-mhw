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

def pypl66(x,dt=1,T=33):
    cutoff = T/dt
    fq = 1/cutoff
    nw = np.round((2*T)/dt)
    print("number of weights = "+str(nw))
    nw2 = 2*nw

    [npts,ncol] = x.shape()
    if npts < ncol:
        x = x.T
        [npts,ncol] = x.shape()
    xf = x

    # filter weights
    j = np.arange(1,nw+1)
    t = np.pi*j
    den = (fq**2)*(t**3)
    wts = (2*np.sin(2*fq*t)-np.sin(3*fq*t))/den

    # make weights symmetric
    wts = [np.flip(wts),2*fq,wts]
    wts = wts/wts.sum()

    # fold tapered time series
    cs = cos(t.T/nw2)
    jm = np.flip(j)
    
    for ic in range(1, ncol+1) :
        jgd = x[:,ic].where(x != nan,drop=True)
        npts = jgd.size()
        if npts > nw2:
            # detrend time series, then add back trend after filtering
            xdt = sig.detrend(x[jgd,ic])
            trnd = x[jgd,ic]-xdt
            y = [[cs[jm]*xdt[jm]],[xdt],[cs[j]*xdt[npts+1-j]]]

            # filter
            yf = sig.lfilter(wts, 1.0, y)

            # strip off extra points
            xf[jgd,ic] = yf[nw2+1:npts+nw2]

            # add back trend
            xf[jgd,ic] += trnd
        else:
            print("warning: time series is too short")

    return xf

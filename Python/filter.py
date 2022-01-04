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
import xarray as xr
import pandas as pd
#import matplotlib.pyplot as plt

def lowpass(y,freq, win):
    
    # Get filter weights
    w = lanczos(freq,win)
    
    if np.sum(np.isnan(np.squeeze(np.array(y.data)))) == 0:
        ylow = np.convolve(y,w,mode='same')
    else:
        YLOW = np.convolve(y[np.isfinite(y)],w,mode='same')
        ylow = y.copy(data=np.nan*np.ones((len(y.data),)))
        ylow[np.isfinite(y)] = YLOW
    
    # Throwing out ends of time series that show edge effects
    ylow.data[0:(win-1)//2] = np.nan*np.ones((1,win//2))
    ylow.data[len(YLOW.data)-(win-1)//2:len(YLOW.data)] = np.nan*np.ones((1,win//2)) # This indexing with YLOW only works while the non-Nan data starts at the begining of the dataset. It was needed in this case because the end of the dataset that was excluded in the convolution is a length of Nan's longer than half the window.
    
    # High-passed is left behind after removing low-passed from raw
    yhigh = y - ylow
    
    return ylow, yhigh

def lanczos(freq,win):
    # freq is related by 2*pi/period
    # win is window length
    
    sig = 2*np.pi/(win)
    C = freq/np.pi
    t = np.arange(1,(win/2))
    
    halfwin = C*np.sinc(freq*t/np.pi)*np.sinc(sig*t/np.pi)
    fullwin = np.concatenate((np.flip(halfwin),[C],halfwin))
    w = fullwin/np.nansum(fullwin)    
    
    return w

def lowpass2(y,freq,win):
    
    # Get filter weights as a data array
    weights = xr.DataArray(lanczos(freq,win), dims=['window'])
    
    # create rolling object with centered window of length win along time dim
    # construct calls the last (outer most) dimension position 'window'
    # dot takes the inner product of the windowed and weights, equivalent to convolve, but any nan present makes result nan
    # Automatically any overlap of window off edges of array will insert nan's where there is overlap
    ylow = y.rolling(time=win, center=True, keep_attrs=True).construct('window').dot(weights) 
    
    # high-passed signal is the difference between original and low-passed
    yhigh = y - ylow
    
    return ylow, yhigh

def climatology(y, interp_time):
    
    time = pd.date_range(start="1979-01-01 00:00",end="2020-09-18 21:00",freq='D')
    # recordtime = y.time
    name0 = str('clim')+str(y.name)
    attrs0 = y.attrs
    latitude0 = y.latitude
    longitude0 = y.longitude
    expver0 = y.expver
    datLeap = []
    
    # groupby month or slice by month
    month = y.groupby("time.month")
    
    # loop for month
    for i in range(12):
        workdat = y[month.groups.get(i+1)] # get the data for each month
        # groupby day and take mean
        daymean = workdat.groupby("time.day").mean("time")
        # append to the 366 day record
        datLeap = np.append(datLeap,daymean)
    
    # datLeap = data0
    datNorm = np.delete(datLeap,59)
    data0 = []
    
    # loop to form a climatology for the length of the record 
    for year in np.arange(1979,2021):
        # If the selected year is a leap year then add a leap year set, otherwise use the 365-day year
        if np.isin(year,np.append(np.arange(2000,1977,-4),np.arange(2000,2021,4)))==1:
            data0 = np.hstack([data0,datLeap])
        else:
            data0 = np.hstack([data0,datNorm])    
    
    y0 = xr.DataArray(data=data0[0:len(time)],coords=[time],dims=["time"],attrs=attrs0,name=name0)
    
    # y0 = xr.DataArray(data=data0[0:len(time)-1],coords=[longitude0,latitude0,expver0,time],dims=["longitude","latitude","expver","time"],attrs=attrs0,name=name0) # This will be the line to create the climatology data array object with attributes and names related to the original dataset

    y0 = y0.resample(time=interp_time).interpolate()
    y0 = y0.loc["1979-01-01 00:00":"2020-09-17 21:00"]
    return y0


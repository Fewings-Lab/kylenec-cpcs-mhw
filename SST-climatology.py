# -*- coding: utf-8 -*-
# SST-climatology.py
# Python 3.7
"""
@author: %(username)s
Created on %(date)s
Modified %(dates)s

Description:

"""
# Import relevant libraries
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import xarray as xr
import filter as filt # Homemade filtering library

# Set path to the data files
data_path = Path("C:/Users/cooleyky/Desktop/kylenec-test/Data/")
data_files = data_path / "Constitucion-SST-*.nc"

# Open all SST files and combine into single dataset
SSTdataset = xr.open_mfdataset(str(data_files))

# Plot of unfiltered timeseries
sst = SSTdataset.sst.isel(expver=0,latitude=0,longitude=0) - 273.15
sst.attrs = SSTdataset.sst.attrs
sst.attrs["units"] = "$^\circ$C"
f1, ax1 = plt.subplots()
sst.plot(ax=ax1)
plt.title("ERA-5 SST at "+str(sst.latitude.data)+", "+str(sst.longitude.data))
plt.xlabel("Date")

# Plot of unfiltered timeseries in 2016
sst2016 = sst.loc["2015-12-31":"2017-01-01"]
f2, ax2 = plt.subplots()
sst2016.plot(ax=ax2)
plt.title("ERA-5 SST at "+str(sst.latitude.data)+", "+str(sst.longitude.data)+" in 2016")
plt.xlabel("Date")
currentaxis= plt.gca()
sstlabel = currentaxis.get_ylabel()

# Use a low-pass lanczos filter on timeseries
Tcutoff = 168 # This is the maximum period in hours (dt = 1 hr) we want to exclude from the low-passed signal, also to take care of possible effects from a change in AVHRR product used in ERA-5
# 168 hrs = 7 days (observations of events in CCS show periods of ~8 days)
freq = 2*np.pi/Tcutoff # Convert cutoff period to angular frequency
win = 363 # Window length of filter is odd to be symmetric
sstL, sstH = filt.lowpass2(sst,freq,win)

# Plotting filtered arrays
f3, ax3 = plt.subplots(nrows=2,sharex=True)

# Low-passed signal
sst.plot(ax=ax3[0],color=[0.80, 0.80, 0.80],linewidth=2,zorder=-1)
sstL.plot(ax=ax3[0],color="blue",linestyle='--')
ax3[0].set_title("Low-passed (T= "+str(Tcutoff)+" hr) ERA-5 SST at "+str(sst.latitude.data)+", "+str(sst.longitude.data))
#axes[0].xaxis.set_visible(False)
ax3[0].set_xlabel(' ')
ax3[0].set_ylabel(sstlabel)
plt.ylim((sst.min(dim="time"),sst.max(dim="time")))

# High-passed signal
#SST.plot(ax=axes[1],color=[0.80, 0.80, 0.80],linewidth=2,zorder=-1)
sstH.plot(ax=ax3[1],color="red")
ax3[1].set_title("High-passed ERA-5 SST at "+str(sst.latitude.data)+", "+str(sst.longitude.data))
plt.xlabel("Date")
ax3[1].set_ylabel(sstlabel)
plt.ylim((sst.min(dim="time")-sst.mean(dim="time"),sst.max(dim="time")-sst.mean(dim="time"))) #changes y-axis scale to be same as top plot

# Climatology that separates out Feb 29, homebrew function
clim = filt.climatology(sstL,"1H")
# Climatology that goes by yearday, built in function
sst0 = sstL.groupby("time.dayofyear").mean("time",skipna=True)

# Plot to compare these on the same axis
f4, ax4 = plt.subplots()
plt.plot(sst0.dayofyear,clim[0:8784:24].values)
sst0.plot(ax=ax4)
plt.title("Climatology ")
plt.legend({'Sorting by day of month','Sorting by day of year'})
ax4.set_ylabel(sstlabel)

# Obtain SST' using built in climatology method and plot
sstA = sstL.groupby("time.dayofyear")-sst0
sstA.drop_vars({"dayofyear","expver"})
mean = sstA.mean()
twosig = 2*sstA.std()
f5, ax5 = plt.subplots()
sstA.plot(ax=ax5,zorder=-1)
plt.axhline(mean,ls='dashed',color='k')
plt.axhline(mean+twosig,linestyle='dashed',color='r')
plt.axhline(mean-twosig,linestyle='dashed',color='r')
plt.title("SST' at "+str(sst.latitude.data)+", "+str(sst.longitude.data))
plt.xlabel("Date")
plt.ylabel(sstlabel)
plt.legend(["SST'",r"$\mu_{SST'}$",r"$\pm 2 \sigma_{SST'}$"])

# Histogram of time where abs(SST') exceeds 2std by month
centered = sstA - mean # anomalies detrended so that mean=0
isexceed = np.abs(centered) >= twosig # boolean for where these exceed 2std
anomalies = centered[isexceed] # picking out those values from the data
month = anomalies["time.month"] # array with month number for each point
f6, ax6 = plt.subplots() # open figure
month.plot.hist(bins=np.arange(0.5,13.5,1)) # histogram with month number bins
plt.title("SST' exceeding $\pm 2 \sigma$") # title 
plt.xticks(ticks=np.arange(1,13),labels=["J","F","M","A","M","J","J","A","S","O","N","D"]) #ticks at center of bins with first letter of month for label
plt.xlabel("Month") # x-label
plt.ylabel("Total hours") # y-label


# Histogram by season
season = anomalies["time.season"] # array with month number for each point
f7, ax7 = plt.subplots() # open figure
season.plot.hist() # histogram with month number bins
# season.plot.hist(bins=["DJF","MAM","JJA","SON"]) # histogram with month number bins
plt.title("SST' exceeding $\pm 2 \sigma$") # title 
plt.xlabel("Season") # x-label
plt.ylabel("Total hours") # y-label

# Histogram with high/low events separate
top = centered >= twosig # boolean for where these exceed 2std
bottom = centered <= -twosig # boolean for where these exceed 2std
anomtop = centered[top] # picking out those values from the data
anombot = centered[bottom] # picking out those values from the data
montht = anomtop["time.month"] # array with month number for each point
monthb = anombot["time.month"]
f8, ax8 = plt.subplots() # open figure
plt.hist([monthb, montht],bins=np.arange(0.5,13.5,1))
plt.title("SST' exceeding $\pm 2 \sigma$") # title 
plt.xticks(ticks=np.arange(1,13),labels=["J","F","M","A","M","J","J","A","S","O","N","D"]) #ticks at center of bins with first letter of month for label
plt.xlabel("Month") # x-label
plt.ylabel("Total hours") # y-label
plt.legend(["below $-2\sigma$","over $2\sigma$"])

# Histogram with abs(SST')<2sig high/low events separate
topW = (centered < twosig) & (centered > 0) # boolean where anom in [0,2sig)
bottomW = (centered > -twosig) & (centered < 0) # boolean where anom in [0,-2sig)
anomtopW = centered[topW] # picking out those values from the data
anombotW = centered[bottomW] # picking out those values from the data
monthtW = anomtopW["time.month"] # array with month number for each point
monthbW = anombotW["time.month"]
f9, ax9 = plt.subplots() # open figure
plt.hist([monthbW, monthtW],bins=np.arange(0.5,13.5,1))
# monthbW.plot.hist(bins=np.arange(0.5,13.5,1))
# monthtW.plot.hist(bins=np.arange(0.5,13.5,1)) # histogram with month number bins
plt.title("SST' between 0 and $\pm 2 \sigma$") # title 
plt.xticks(ticks=np.arange(1,13),labels=["J","F","M","A","M","J","J","A","S","O","N","D"]) #ticks at center of bins with first letter of month for label
plt.xlabel("Month") # x-label
plt.ylabel("Total hours") # y-label
plt.legend(["down to $-2\sigma$","up to $2\sigma$"])


# Take dSST'/dt 
dt1 = sstA.differentiate("time")
# plot SST' and time derivative in separate subplots
f10, [ax10, ax11] = plt.subplots(2,1,sharex=True)
sstA.plot(ax=ax10)
ax10.set_title('SST\'')
ax10.set_xlabel("")
ax10.set_ylabel(r"SST' [$^\circ$C]")
dt1.plot(ax=ax11)
ax11.set_title(r"$ \frac{\partial \mathrm{SST}' }{\partial t} $")
ax11.set_ylabel(r"Change in SST' [$^\circ$C/hr]")

# dSST'/dt using a first-order difference approximation
# Consider how rate of change is "rise/run" -- rise is the difference between successive points in time, run is the 1-hr between each point
dt2 = sstA[1:-1]-sstA[0:-2].values
dt2 = dt2[1:-1]
f11, ax12 = plt.subplots()
dt2.plot(ax=ax12)
ax12.set_title(r"$ \frac{\partial \mathrm{SST}' }{\partial t} $")
ax12.set_ylabel(r"Change in SST' [$^\circ$C/hr]")

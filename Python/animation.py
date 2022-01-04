# -*- coding: utf-8 -*-
# animation.py
# Python 3.7
"""
@author: %(username)s
Created on %(date)s
Modified %(dates)s

Description:

"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import xarray as xr
import filter as filt
import pandas as pd

# import warnings 
# warnings.filterwarnings("ignore") 

data_path = Path("C:/Users/cooleyky/Desktop/kylenec-test/")
data_files = data_path / "Data" / "ChileCoast*.nc"

ds = xr.open_mfdataset(str(data_files)) # opens multiple data files as one dataset

sst = ds.sst.isel(expver=0) - 273.15
sst.attrs = ds.sst.attrs # keeps attributes from original dataset
sst.attrs["units"] = "degC"

Tcutoff = 48 # This is the maximum period in hours (dt = 1 hr) we want to exclude from the low-passed signal
freq = 2*np.pi/Tcutoff # Convert cutoff period to angular frequency
win = 121 # Window length of filter is odd to be symmetric
sstL, sstH = filt.lowpass2(sst,freq,win)
sstL.attrs = sst.attrs
sstL = sstL.isel(time=slice((win//2),(-247-win//2)))

# def plot(temp:xr.DataArray):
#     import matplotlib.pyplot as plt
#     f,ax = plt.subplots(1,1)
#     temp.squeeze().plot()
#     fname = str(temp.time.squeeze().values).replace(":", "-")+".png"
#     f.savefig(data_path/"animation"/fname)
#     plt.close(f)
#     return temp.time

# plot(sstL.isel(time=54100))

### Using groupby to sort by date and do climatology

# First trimming down the domain and removing nan's from mean
sstLNW = sstL.sel(longitude=slice(-90,-70)) #cut out some land mass
sstLNW = sstLNW.fillna(value=-3.5) # fill in nan's over land

# Take climatology over NW ocean section with land mask
clim = sstLNW.groupby("time.dayofyear").mean("time",skipna=True)
sstAnomNW = sstLNW.groupby("time.dayofyear") - sstLNW.groupby("time.dayofyear").mean("time",skipna=True)

# mapped = (sstAnom.chunk({"time":1,"latitude":-1,"longitude":-1}).map_blocks(plot))

# mapped.compute()

# sstAnom.chunk({"time":1,"latitude":-1,"longitude":-1}).map_blocks(plot).compute()

# sstAnom = sstAnom.chunk({"time":-1,"latitude":-1,"longitude":-1}) # undoes chunks made by groupby in line 45
sstAnomNW = sstAnomNW.drop_vars({"dayofyear","expver"})
# sstAnom.chunk({"time":1,"latitude":-1,"longitude":-1}).map_blocks(plot)
time = pd.date_range(start="2016-02-01",end="2016-03-01",freq="6H")

# sstAnomNW = sstAnom.sel(longitude=slice(-90,-70)) #cut out land mass, hopefully speeds up some
# sstAnomNW = sstAnomNW.fillna(value=-3.5)
# sstAnomNW.chunk([61124,101,81])
# np.seterr(divide='ignore', invalid='ignore')
for i in time:
    f,ax = plt.subplots(1,1)
    sstAnomNW.sel(time=i).plot.pcolormesh(vmax=3.5) # sets colormap -3.5:3.5 with "center" for everyplot
    ax.set_aspect('equal')
    plt.title(sstAnomNW.sel(time=i).time.dt.strftime("%b %d, %Y").values)
    fname = str(sstAnomNW.sel(time=i).time.squeeze().values).replace(":", "-")+".png"
    f.savefig(data_path/"animation3"/fname)
    plt.close(f)

### Using climatology function I wrote, make climatology by lat and lon
# sstClim = sstL.copy()

# for lat in range (101):
#     for lon in range(121):
#         clim = filt.climatology(sstL.isel(latitude=lat,longitude=lon),"6H")
#         sstClim[:,lat,lon] = clim
        
# Then take anomaly by time step

import imageio
images = []
pngpath = data_path/"animation3"

for filename in pngpath.iterdir():
    images.append(imageio.imread(filename))
imageio.mimsave(data_path/"animation3"/"2016-01-SST-anomaly.gif", images)


        
# -*- coding: utf-8 -*-
# SST-climatology.py
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

data_path = Path("C:/Users/cooleyky/Desktop/kylenec-test/Data/")
data_files = data_path / "Constitucion-SST-*.nc"

SSTdataset = xr.open_mfdataset(str(data_files))

SST = SSTdataset.sst.isel(expver=0) - 273.15
SST.attrs = SSTdataset.sst.attrs
SST.attrs["units"] = "deg C"
plt.figure()
SST.plot()
plt.title("ERA-5 SST at "+str(SST.latitude.data[0])+"W, "+str(SST.longitude.data[0])+"S")
plt.xlabel("Date")
plt.show()

SST2016 = SST.loc["2015-12-31":"2017-01-01"]
plt.figure()
SST2016.plot()
plt.title("ERA-5 SST at "+str(SST.latitude.data[0])+"W, "+str(SST.longitude.data[0])+"S in 2016")
plt.xlabel("Date")
plt.show()
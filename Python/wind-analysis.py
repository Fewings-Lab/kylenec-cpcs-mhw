#!/usr/bin/env python

# wind-analysis.py
# 04 Sept 2020
# Kylene Cooley

import warnings 
warnings.filterwarnings("ignore") 

# import libraries
import datetime as dt
import xarray as xr
import fsspec
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from pyPL66 import pypl66
# import physoce.tseries
# make datasets display nicely
xr.set_options(display_style="html")

v10 = xr.open_dataset("V10_2010_2020.nc")

ldv = v10.isel(longitude=[5],latitude=[5],expver=[0])
ldv = ldv.to_array()
ldv = ldv[:,:,0,0,0]
# ldv.size

ldvLo = ldv.copy(data=(pypl66(ldv.data)).T)
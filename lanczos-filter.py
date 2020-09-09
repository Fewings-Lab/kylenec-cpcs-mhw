#!/usr/bin/env python

# lanczos-filter.py
# 08 Sept 2020

"""

Description: From SOI example on Iris 1.2 website https://scitools.org.uk/iris/docs/v1.2/examples/graphics/SOI_filtering.html

"""
import numpy as np
import matplotlib.pyplot as plt
# UDUNITS2_XML_PATH="C:/Users/cooleyky/AppData/Local/Continuum/miniconda3/Library/share/udunits/udunits2.xml"
import iris
import iris.plot as iplt


def low_pass_weights(window, cutoff):
    """Calculate weights for a low pass Lanczos filter.

    Args:

    window: int
        The length of the filter window.

    cutoff: float
        The cutoff frequency in inverse time steps.

    """
    order = ((window - 1) // 2 ) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1., n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
    w[n-1:0:-1] = firstfactor * sigma
    w[n+1:-1] = firstfactor * sigma
    return w[1:-1]


def main():

    # load the hourly 10-meter northward wind velocity time-series
    fname = 'V10_2010_2020.nc'
    nwv = iris.load_cube(fname)

    # window length for filters
    window = 12
    # This gives us a filter with a 123 weights

    # construct 33-hour low-pass filter to remove approx diurnal and higher freq signals
    # for the 10-m northward wind velocity which is hourly
    wgts33 = low_pass_weights(window, 1. / 33.)  #
    # wgts84 = low_pass_weights(window, 1. / 84.) # formerly for 7-year (84-month) low pass filters

    # apply the filters using the rolling_window method with the weights
    # keyword argument
    nwvL = nwv.rolling_window('time',
                               iris.analysis.MEAN,
                               len(wgts33),
                               weights=wgts33)
    # nwv84 =  nwv.rolling_window('time',
    #                             iris.analysis.MEAN,
    #                             len(wgts84),
    #                             weights=wgts84)

    nwvH = nwv - nwvL # subtracting the low-passed component from the raw signal to give us the high-passed part
    
    # plot the V10 time series and both filtered versions
    fig = plt.figure(figsize=(9, 4))
    iplt.plot(nwv, coords=['time'], color='0.7', linewidth=1., linestyle='-',
              alpha=1., label='no filter')
    iplt.plot(nwvL, coords=['time'], color='b', linewidth=2., linestyle='-',
              alpha=.7, label='Low-pass filter')
    iplt.plot(nwvH, coords=['time'], color='r', linewidth=2., linestyle='-',
              alpha=.7, label='High-pass filter')
    #plt.ylim([-4, 4])
    plt.title('10-meter Northward Wind Velocity')
    plt.xlabel('Time')
    plt.ylabel('V_10 [m/s]')
    plt.legend(fontsize=10)
    plt.show()


if __name__ == '__main__':
    main()
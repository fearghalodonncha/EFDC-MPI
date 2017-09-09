"""
Created on Wed Sep 6 13:49:11 2017

@author: ernestoa@ie.ibm.com

Synopsis

    TODO contour [-h,--help] [-v,--verbose] [--version]

Description

    TODO This describes how to use this script. This docstring
    will be printed by the script if there is an error or
    if the user requests help (-h or --help).

Examples

    TODO: Show some examples of how to use this script.

Exit Status

    TODO: List exit codes

Author

    TODO: Ernesto Arandia <ernestoa@ie.ibm.com>

License

    TODO: License information.

Version

    $Id$

"""


import argparse
import os
import traceback
import time
import sys
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import warnings


def getdataset(ncfilepath, datetime):
    ncfilename = ncfilepath + 'efdcout_' + datetime + '.nc'
    ncdataset = Dataset(ncfilename)
    return ncdataset


def getvariablevalues(data, variable, zlayer):
    var0 = data.variables[variable][:]
    var1 = var0[0]
    var2 = var1[zlayer - 1, :, :]
    if variable == 'Temperature':
        var2 = var2 - 273
    var2[var2 < -900] = np.nan
    return var2


def trimdata(data, lowerpercentile=1, upperpercentile=99):
    pc1 = np.nanpercentile(np.concatenate(data), lowerpercentile)
    pc2 = np.nanpercentile(np.concatenate(data), upperpercentile)
    data[data < pc1] = np.nan
    data[data > pc2] = np.nan
    return data


def getcoordinates(ncdataset):
    x = ncdataset.variables['X'][:]
    y = ncdataset.variables['Y'][:]
    return x, y


def makecontourplot(data, x1, y1, figfilename, trimflag=False):
    if trimflag:
        data = trimdata(data)
    plt.clf()
    plt.imshow(data, origin='lower', extent=[x1.min(), x1.max(), y1.min(), y1.max()])
    plt.contourf(x1, y1, data)
    plt.colorbar()
    plt.savefig(figfilename)


def main():

    # parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfiledir')
    parser.add_argument('datetime')
    parser.add_argument('zlayer')
    parser.add_argument('variable')
    parser.add_argument('outputfiguredir')
    args = parser.parse_args()

    # read the netcdf data
    dataset = getdataset(args.inputfiledir, args.datetime)

    # get the desired variable values
    var = getvariablevalues(dataset, args.variable, zlayer=int(args.zlayer))

    # get the XY coordinates
    x = getcoordinates(dataset)[0]
    y = getcoordinates(dataset)[1]

    # create contour plot
    figfilename = args.outputfiguredir + args.variable + '_contour_' + args.zlayer + '_' + args.datetime + '.pdf'
    makecontourplot(var, x, y, figfilename, True)


if __name__ == '__main__':
    start_time = time.time()
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    main()



import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

def read_drifter(filename):
    with open(filename) as f:
        lines = f.readlines()
    NPD = float(lines[3].split()[0])   ## NPD, number of particles specified on line 4
    times_list = lines[4::2]
    drifter_list = lines[5::2]
    times_np = np.zeros([len(times_list)])
    drift_x = np.zeros([len(times_list), int(NPD)])
    drift_y = np.zeros([len(times_list), int(NPD)])
    drift_z = np.zeros([len(times_list), int(NPD)])


    for t in range(0, len(times_list)):
        times_np[t] = float(times_list[t].split()[0])
        for d in range(0, int(NPD)):
            if t == 0:
                step = 3
                Lall = 1
            else:
                step = 3
                Lall = 1
            drift_x[t,d] = float(drifter_list[t].split()[1 - Lall + (d*step)])
            drift_y[t,d] = float(drifter_list[t].split()[2 - Lall + (d*step)])
            drift_z[t,d] = float(drifter_list[t].split()[3 - Lall + (d*step)])

    drift_x[drift_x == 0] = np.nan
    drift_y[drift_y == 0] = np.nan
    ## Convert drifter coordinates to km
    drift_x = drift_x/1000.
    drift_y = drift_y/1000.
    drift_z = drift_z # vertical stays as metre
    
    return times_np, drift_x, drift_y, drift_z

def getcoordinates(ncdataset, zeta = 'h'):
    """
    The ADCP data report coordinates in terms 
    of 'latitude', 'longitude', 'time', & 'depth'
    """
    x = ncdataset.variables['X'][:]
    y = ncdataset.variables['Y'][:]
    time = ncdataset.variables['Time'][:]
    depth = ncdataset.variables[zeta][:,:]
    return x, y, time, depth
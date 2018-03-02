import csv
import os
import numpy as np
import glob
from math import sin, cos, asin, sqrt, degrees, radians


def haversine(angle_radians):
    return sin(angle_radians / 2.0) ** 2  # haversine = sin^2(ang/2)


def inverse_haversine(h):
    return 2 * asin(sqrt(h))  # radians


def distance_between_points(lat1, lon1, lat2, lon2):
    # all args are in degrees
    # WARNING: loss of absolute precision when points are near-antipodal
    lat1 = radians(lat1)
    lat2 = radians(lat2)
    dlat = lat2 - lat1
    dlon = radians(lon2 - lon1)
    h = haversine(dlat) + cos(lat1) * cos(lat2) * haversine(dlon)
    return RADIUS * inverse_haversine(h)


def indices(a, func):
    return [i for (i, val) in enumerate(a) if func(val)]

def reconcile_efdc_state():
    # list all efdc_files and reconcile into 1
    path_to_write = './'
    efdc_files = glob.glob(path_to_write + 'EFDC_temp*.csv')

    recon = np.empty((0, 4))
    for file in efdc_files:
        a = np.loadtxt(file)
        recon = np.append(recon, a, axis=0)
        print(file)

    filewrite = path_to_write + 'EFDC_state.csv'
    np.savetxt(filewrite, recon, fmt='%d, %d, %10.5f, %10.5f')


def main():
    Earth_radius_km = 6371.0
    RADIUS = Earth_radius_km
    fname = 'CODGRID.INP'
    infile = fname
    reader = csv.reader(file(infile, 'r'), delimiter=',')
    values_read = set()
    for line in reader:
        lat5 = float(line[5])
        lon5 = float(line[4])
        I = int(float(line[0]))
        J = int(float(line[1]))
        a = lon5, lat5, I, J
        values_read.add(a)

    dd = np.array(list(values_read))


    # These define spatial extents for interpolation
    westbnd = 0
    eastbnd = 99999
    southbnd = 0
    northbnd = 99999

    efdc_delx = 0.5  # km
    codar_delx = 2  # km
    interp_rng = int(codar_delx / efdc_delx)
    strt = 1 - (interp_rng / 2)
    endt = 1 + (interp_rng / 2)

    westbnd = westbnd - strt  # set bounds on the domain read in from the EFDC model
    eastbnd = eastbnd - endt  # model based on the extents of the model grid and confidence
    southbnd = southbnd - strt  # in the Codar measurements around the footprint edges
    northbnd = northbnd - endt

    inf2 = 'codtemp.csv'
    reader2 = csv.reader(file(inf2, 'r'), delimiter=',')

    coords = list()
    ud = list()
    vd = list()
    ud_err = list()
    vd_err = list()
    firstline = True
    for line in reader2:
        if firstline:
            firstline = False
            continue
        lon = float(line[0])
        lat = float(line[1])
        u = float(line[2])
        v = float(line[3])
        uerr = float(line[4])
        verr = float(line[5])
        setmin = 5
        for j in dd:
            dist = distance_between_points(j[1], j[0], lat, lon)
            if (dist < setmin):
                setmin = dist
                map_i = j[0]
                map_j = j[1]
                x = j[2]
                y = j[3]
        if (setmin < 0.5 and uerr != 999 and verr != 999):
            if (x > westbnd and x < eastbnd):
                if (y > southbnd and y < northbnd):
                    for ii in range(strt, endt):
                        for jj in range(strt, endt):
                            a = x + ii, y + jj
                            coords.append(a)
                            ud.append(u)
                            vd.append(v)
                            ud_err.append(uerr)
                            vd_err.append(verr)
    du = {}
    dv = {}
    du_e = {}
    dv_e = {}
    for a, b in zip(coords, ud):
        du.setdefault(a, []).append(b)

    for a, b in zip(coords, vd):
        dv.setdefault(a, []).append(b)

    for a, b in zip(coords, ud_err):
        du_e.setdefault(a, []).append(b)

    for a, b in zip(coords, vd_err):
        dv_e.setdefault(a, []).append(b)

    avg_u = []
    avg_v = []
    avg_ue = []
    avg_ve = []
    lonlat = []

    for key in du:
        lonlat.append(key)
        avg_u.append(np.mean(du[key]))

    for key in dv:
        avg_v.append(np.mean(dv[key]))

    for key in du_e:
        avg_ue.append(np.mean(du_e[key]))

    for key in dv_e:
        avg_ve.append(np.mean(dv_e[key]))

    lonlat = np.array(lonlat)
    avg_u = np.array(avg_u)
    avg_v = np.array(avg_v)
    avg_ue = np.array(avg_ue)
    avg_ve = np.array(avg_ve)

    avg_u = avg_u.reshape(len(avg_u), 1)
    avg_v = avg_v.reshape(len(avg_v), 1)
    avg_ue = avg_ue.reshape(len(avg_ue), 1)
    avg_ve = avg_ve.reshape(len(avg_ve), 1)

    comb = np.concatenate((lonlat, avg_u, avg_v, avg_ue, avg_ve), axis=1)

    np.savetxt('Blue_files/CODAR.csv', comb, delimiter=',', fmt="%d,%d,%.3f,%.3f,%.3f,%.3f")


    # EFDC MPI model wrote individual files for each subdomain;
    # reconcile to a single file
    reconcile_efdc_state()

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
################################
# Program is part of TintPy    #
# Copyright (c) 2022, Lei Yuan #
# Author: Lei Yuan, 2022       #
################################

import argparse
import os
import sys
import numpy as np

EXAMPLE = """Example:
  # sar2geo (range, azimuth ==> longitude, latitude)
  python3 geocode.py lookup_fine dem_seg.par sar2geo 100 100

  # geo2sar (longitude, latitude ==> range, azimuth)
  python3 geocode.py lookup_fine dem_seg.par geo2sar 100.3 43.5
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(description="Perform coordinate transformation using lookup table (GAMMA format file).",
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)
    parser.add_argument('lookup', help='lookup table.')
    parser.add_argument('dem_seg_par', help='dem_seg.par.')
    parser.add_argument('flag', choices=['sar2geo', 'geo2sar', 'SAR2GEO', 'GEO2SAR'], help='flag (geo2sar, sar2geo).')
    parser.add_argument('xy', nargs=2, type=float, help='xy value in SAR or GEO coordinate.')

    inps = parser.parse_args()
    return inps


def read_gamma_par(par_file, keyword):
    """Extract value from par_file using keyword

    Args:
        par_file (str): GAMMA parameter file
        keyword (str): keyword like "reange_sample"
    """
    value = None
    with open(par_file, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            if line.count(keyword) == 1:
                value = line.split()[1].strip()

    return value


def read_lookup_table(file, lines):
    """Read lookup table file

    Args:
        file (str): lookup table
        lines (int): length of lookup table

    Returns:
        array: data
    """
    data = np.fromfile(file, dtype='complex64')
    # GAMMA output files are big-endian
    data.byteswap('True')
    data = data.reshape(lines, -1)

    return data


def geo2sar(lon_lat, dem_seg_par, lookup_table):
    """Transform longitude and latitude to SAR coordinates

    Args:
        lon_lat (list): longitude and latitude
        dem_seg_par (str): dem par
        lookup_table (str): lookup file

    Returns:
        tuple: (range, azimuth)
    """
    lines = int(read_gamma_par(dem_seg_par, 'nlines'))
    upper_left_lat = float(read_gamma_par(dem_seg_par, 'corner_lat'))
    upper_left_lon = float(read_gamma_par(dem_seg_par, 'corner_lon'))
    lat_step = float(read_gamma_par(dem_seg_par, 'post_lat'))
    lon_step = float(read_gamma_par(dem_seg_par, 'post_lon'))

    data = read_lookup_table(lookup_table, lines)

    range_data = np.real(data)
    azimuth_data = np.imag(data)

    lon, lat = lon_lat
    xx = int((lat - upper_left_lat) / lat_step)
    yy = int((lon - upper_left_lon) / lon_step)

    rng = int(range_data[xx, yy])
    azi = int(azimuth_data[xx, yy])

    return rng, azi


def sar2geo(rng_azi, dem_seg_par, lookup_table):
    """Transform range and azimuth to GEO coordinates

    Args:
        rng_azi (list): range and azimuth
        dem_seg_par (str): dem par
        lookup_table (str): lookup file

    Returns:
        tuple: (range, azimuth)
    """
    lines = int(read_gamma_par(dem_seg_par, 'nlines'))

    upper_left_lat = float(read_gamma_par(dem_seg_par, 'corner_lat'))
    upper_left_lon = float(read_gamma_par(dem_seg_par, 'corner_lon'))

    lat_step = float(read_gamma_par(dem_seg_par, 'post_lat'))
    lon_step = float(read_gamma_par(dem_seg_par, 'post_lon'))

    data = read_lookup_table(lookup_table, lines)

    cpx_input = complex(f"{rng_azi[0]}+{rng_azi[1]}j")

    distance = abs(cpx_input - data)

    idx = np.where(distance == distance.min())

    lat = upper_left_lat + idx[0] * lat_step
    lon = upper_left_lon + idx[1] * lon_step

    return lon[0], lat[0]


def main():
    # get inps
    inps = cmdline_parser()
    lookup = os.path.abspath(inps.lookup)
    dem_seg_par = os.path.abspath(inps.dem_seg_par)
    flag = inps.flag
    xy = inps.xy

    # check files
    if not os.path.isfile(lookup):
        sys.exit(f"{lookup} does not exist, please check it.")

    if not os.path.isfile(dem_seg_par):
        sys.exit(f"{dem_seg_par} does not exist, please check it.")

    if flag.upper() == 'SAR2GEO':
        x, y = xy
        if x <= 0:
            sys.exit("Value in SAR coordinate must be bigger than 0.")
        if y <= 0:
            sys.exit("Value in SAR coordinate must be bigger than 0.")
        lon, lat = sar2geo(xy, dem_seg_par, lookup)
        rng, azi = int(xy[0]), int(xy[1])
        print(f"SAR2GEO: {rng} {azi} ==> {lon} {lat}")

    if flag.upper() == 'GEO2SAR':
        rng, azi = geo2sar(xy, dem_seg_par, lookup)
        lon, lat = xy
        print(f"GEO2SAR: {lon} {lat} ==> {rng} {azi}")

if __name__ == "__main__":
    main()

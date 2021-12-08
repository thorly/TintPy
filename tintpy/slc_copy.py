#!/usr/bin/env python3
################################
# Program is part of TintPy    #
# Copyright (c) 2021, Lei Yuan #
# Author: Lei Yuan, 2021       #
################################

import argparse
import glob
import os
import re
import shutil
import sys
import zipfile
from xml.dom.minidom import parse

import numpy as np


def cmd_line_parser():
    parser = argparse.ArgumentParser(
        description='Copy a common segment from an existing set of SLCs',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLE)

    parser.add_argument('slc_dir', help='SLCs directory')
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('rlks', help='range looks', type=int)
    parser.add_argument('alks', help='azimuth looks', type=int)
    parser.add_argument(
        'method',
        help=
        'get ROI method, rdr for using radar coordinates， geo for using geographic coordinates',
        choices={'rdr', 'geo'})

    parser.add_argument('-r',
                        dest='roff',
                        help='offset to starting rangle sample',
                        type=int)
    parser.add_argument('-nr',
                        dest='nr',
                        help='number of range samples',
                        type=int)
    parser.add_argument('-l',
                        dest='loff',
                        help='offset to starting line',
                        type=int)
    parser.add_argument('-nl', dest='nl', help='number of lines', type=int)

    parser.add_argument('-d',
                        dest='dem_dir',
                        help='dem directory including *.dem and *.dem.par',
                        type=str)
    parser.add_argument('-k', dest='kml', help='ROI kml file', type=str)

    parser.add_argument(
        '-n',
        dest='num',
        help='number of slc used (default: -1, negative number for all)',
        type=int,
        default=-1)
    parser.add_argument('-e',
                        dest='extension',
                        help='file extension for SLCs (defaults: .rslc)',
                        default='.rslc')

    inps = parser.parse_args()

    return inps


EXAMPLE = """Example:
   python3 slc_copy.py /ly/rslc /ly/rslc_cut 8 2 rdr -r 1 -nr 1000 -l 1 -nl 1000
   python3 slc_copy.py /ly/rslc /ly/rslc_cut 8 2 rdr -r 1 -nr 1000 -l 1 -nl 1000 -n 1 -e slc

   python3 slc_copy.py /ly/rslc /ly/rslc_cut 8 2 geo -d dem -k area.kml
   python3 slc_copy.py /ly/rslc /ly/rslc_cut 8 2 geo -d dem -k area.kml -n 1 -e slc
"""


def read_lookup_table(file, lines):
    """Read lookup table file

    Args:
        file (str): lookup table
        lines (int): length of lookup table

    Returns:
        tuple: range and azimuth data
    """
    data = np.fromfile(file, dtype='complex64')
    # GAMMA output files are big-endian
    data.byteswap('True')
    data = data.reshape(lines, -1)

    range_data = np.real(data)
    azimuth_data = np.imag(data)

    return range_data, azimuth_data


def kml2polygon_dict(kml_file):
    """Parse kml or kmz file, get polygon

    Args:
        kml_file (str): kml or kmz file

    Returns:
        dict: polygon
    """
    flag = False
    # unzip kmz to get kml
    if kml_file.endswith('.kmz'):
        flag = True
        dir_name = os.path.dirname(kml_file)
        with zipfile.ZipFile(kml_file, 'r') as f:
            files = f.namelist()
            for file in files:
                if file.endswith('.kml'):
                    f.extract(file, dir_name)
        doc_kml = os.path.join(dir_name, file)
        kml_file = doc_kml

    dom_tree = parse(kml_file)

    if os.path.isfile(kml_file) and flag:
        os.remove(kml_file)

    # parse kml
    root_node = dom_tree.documentElement
    placemarks = root_node.getElementsByTagName('Placemark')

    polygon_dict = {}
    j = 0

    for placemark in placemarks:
        name = placemark.getElementsByTagName('name')[0].childNodes[0].data
        ploygon = placemark.getElementsByTagName('Polygon')[0]
        outerBoundaryIs = ploygon.getElementsByTagName('outerBoundaryIs')[0]
        LinearRing = outerBoundaryIs.getElementsByTagName('LinearRing')[0]
        coordinates = LinearRing.getElementsByTagName('coordinates')[0].childNodes[0].data
        lon_lat = [i.split(',')[0:2] for i in coordinates.strip().split()]
        polygon_dict[name + '_' + str(j)] = np.asarray(lon_lat, dtype='float32')
        j += 1

    return polygon_dict


def gen_grid(dem_seg_par):
    """Generate longitude and latitude grid

    Args:
        dem_seg_par (str): dem_seg.par file

    Returns:
        tuple: longitude and latitude grid
    """
    width = int(read_gamma_par(dem_seg_par, 'width'))
    length = int(read_gamma_par(dem_seg_par, 'nlines'))
    upper_left_lat = float(read_gamma_par(dem_seg_par, 'corner_lat'))
    upper_left_lon = float(read_gamma_par(dem_seg_par, 'corner_lon'))
    lat_step = float(read_gamma_par(dem_seg_par, 'post_lat'))
    lon_step = float(read_gamma_par(dem_seg_par, 'post_lon'))

    lon = np.linspace(upper_left_lon, upper_left_lon + lon_step * width, width)
    lat = np.linspace(upper_left_lat, upper_left_lat + lat_step * length, length)
    lons, lats = np.meshgrid(lon, lat)

    return lons, lats


def find_nearest_point(lons, lats, lon_step, lat_step, point):
    """Find nearest point

    Args:
        lons (array): longitude grid
        lats (array): latitude grid
        lon_step (float): longitude step
        lat_step (float): latitude step
        point (list): [lon, lat]

    Returns:
        tuple: (row, col)
    """
    r, c = None, None
    row, col = lons.shape
    for i in range(row):
        for j in range(col):
            if abs(lons[i, j] - point[0]) <= abs(lon_step / 2) and abs(lats[i, j] - point[1]) <= abs(lat_step / 2):
                r = i
                c = j
                return r, c

    return r, c


def get_rdr_from_kml(lookup_table, dem_seg_par, kml, rlks, alks):
    """Get radar cooridnates

    Args:
        lookup_table (str): lookup table file
        dem_seg_par (str): dem_seg.par file
        kml (str): kml file
        rlks (int): range looks
        alks (int): azimuth looks

    Returns:
        [type]: [description]
    """
    kml_area = kml2polygon_dict(kml)
    lines = int(read_gamma_par(dem_seg_par, 'nlines'))
    lat_step = float(read_gamma_par(dem_seg_par, 'post_lat'))
    lon_step = float(read_gamma_par(dem_seg_par, 'post_lon'))
    range_data, azimuth_data = read_lookup_table(lookup_table, lines)
    lons, lats = gen_grid(dem_seg_par)

    roff, nr, loff, nl = None, None, None, None

    for _, polygon in kml_area.items():
        if len(polygon) == 5:
            roff1 = []
            aoff1 = []
            for i in range(5):
                r, c = find_nearest_point(lons, lats, lon_step, lat_step, polygon[i])
                roff1.append(range_data[r, c] * rlks)
                aoff1.append(azimuth_data[r, c] * alks)

            roff_min, roff_max = int(min(roff1)), int(max(roff1))
            aoff_min, aoff_max = int(min(aoff1)), int(max(aoff1))
            roff = roff_min
            nr = roff_max - roff_min
            loff = aoff_min
            nl = aoff_max - aoff_min
        else:
            sys.exit('The kml file must be a four-point polygon')

    return roff, nr, loff, nl


def gc_map(slc, slc_par, dem, dem_par, rlks, alks, out_dir):
    """Calculate lookup table and DEM related products for terrain-corrected geocoding

    Args:
        slc (str): slc
        slc_par (str): slc par
        dem (str): dem
        dem_par (str): dem par
        rlks (int): range looks
        alks (int): azimuth looks
        out_dir (str): output directory

    Returns:
        tuple: lookup and dem_seg.par
    """
    os.chdir(out_dir)

    date = os.path.basename(slc)[0:8]

    call_str = f"multi_look {slc} {slc_par} {date}.mli {date}.mli.par {rlks} {alks}"
    os.system(call_str)

    dem_seg_par = os.path.join(out_dir, 'dem_seg.par')
    lookup_table = os.path.join(out_dir, 'lookup_table')

    call_str = f"gc_map {date}.mli.par  - {dem_par} {dem} {dem_seg_par} dem_seg {lookup_table} 1 1 sim_sar u v inc psi pix ls_map 8 1"
    os.system(call_str)

    return lookup_table, dem_seg_par


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


def slc_copy(slc, slc_par, roff, nr, loff, nl, out_slc, out_slc_par):
    """Copy a common segment from an existing set of SLCs

    Args:
        slc (str): slc
        slc_par (str): slc par
        roff (int): offset to starting rangle sample
        nr (int): number of range samples
        loff (int): offset to starting line
        nl (int): number of lines
        out_slc (str): output slc
        out_slc_par (str): output slc par
    """
    call_str = f"SLC_copy {slc} {slc_par} {out_slc} {out_slc_par} - - {roff} {nr} {loff} {nl} - -"
    os.system(call_str)


def slc2bmp(slc, slc_par, rlks, alks, bmp):
    """Generate 8-bit raster graphics image of intensity of complex (SLC) data

    Args:
        slc (str): slc file
        slc_par (str): slc par file
        rlks (int): range looks
        alks (int): azimuth looks
        bmp (str): output bmp
    """
    width = read_gamma_par(slc_par, 'range_samples')
    if width:
        call_str = f'rasSLC {slc} {width} 1 0 {rlks} {alks} 1. .35 1 0 0 {bmp}'
        os.system(call_str)


def main():
    # get inputs
    inps = cmd_line_parser()
    slc_dir = os.path.abspath(inps.slc_dir)
    out_dir = os.path.abspath(inps.out_dir)
    rlks = inps.rlks
    alks = inps.alks
    method = inps.method

    # for rdr method
    roff = inps.roff
    nr = inps.nr
    loff = inps.loff
    nl = inps.nl

    # for geo method
    dem_dir = inps.dem_dir
    kml = inps.kml

    num = inps.num
    extension = inps.extension

    # check slc_dir and out_dir
    if not os.path.isdir(slc_dir):
        sys.exit('{} does not exist.'.format(slc_dir))

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # get dates
    dates = sorted([i for i in os.listdir(slc_dir) if re.findall(r'^\d{8}$', i)])

    # check dates
    if len(dates) < 1:
        sys.exit('No slc in {}'.format(slc_dir))

    # check num
    if num < 0 or num > len(dates):
        length = len(dates)
    else:
        length = num

    # check extension
    if not extension.startswith('.'):
        extension = '.' + extension

    if method == 'rdr':
        if roff and nr and loff and nl:
            pass
        else:
            sys.exit('roff nr loff and nl are required parameters for rdr method')

    if method == 'geo':
        if dem_dir and kml:
            dem_dir = os.path.abspath(dem_dir)
            kml = os.path.abspath(kml)

            if not os.path.isdir(dem_dir):
                sys.exit("{} does not exist.".format(dem_dir))
            else:
                dems = glob.glob(dem_dir + '/*.dem')
                dem_pars = [i + '.par' for i in dems]
                for i, j in zip(dems, dem_pars):
                    if os.path.isfile(i) and os.path.isfile(j):
                        dem = dems[0]
                        dem_par = dem + '.par'
                        break
                    else:
                        sys.exit(
                            f'Cannot find *.dem and *.dem.par in {dem_dir}.')

            if not os.path.isfile(kml):
                sys.exit("{} does not exist.".format(kml))

            slc = os.path.join(slc_dir, dates[0], dates[0] + extension)
            slc_par = slc + '.par'

            tmp_dir = os.path.join(slc_dir, 'tmp_for_gc_map')
            if not os.path.isdir(tmp_dir):
                os.mkdir(tmp_dir)

            lookup_table, dem_seg_par = gc_map(slc, slc_par, dem, dem_par, rlks, alks, tmp_dir)
            roff, nr, loff, nl = get_rdr_from_kml(lookup_table, dem_seg_par, kml, rlks, alks)

            shutil.rmtree(tmp_dir)
        else:
            sys.exit(
                'dem_dir and kml file are required parameters for geo method')

    for date in dates[0:length]:
        slc = os.path.join(slc_dir, date, date + extension)
        slc_par = slc + '.par'

        out_slc_dir = os.path.join(out_dir, date)
        out_slc = os.path.join(out_slc_dir, date + extension)
        out_slc_par = out_slc + '.par'

        if os.path.isfile(slc) and os.path.isfile(slc_par):
            if not os.path.isdir(out_slc_dir):
                os.mkdir(out_slc_dir)

            slc_copy(slc, slc_par, roff, nr, loff, nl, out_slc, out_slc_par)

            bmp = out_slc + '.bmp'
            slc2bmp(out_slc, out_slc_par, rlks, alks, bmp)

    print('\nAll done, enjoy it.\n')


if __name__ == "__main__":
    main()

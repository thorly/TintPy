#!/usr/bin/env python3
################################
# Program is part of TintPy    #
# Copyright (c) 2022, Lei Yuan #
# Author: Lei Yuan, 2022       #
################################

import argparse
import os
import sys
import zipfile

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from osgeo import gdal


def read_gdal_file(file, band=1):
    """Read GDAL supported file using GDAL

    Args:
        file (str): GDAL supported file
        band (int, optional): band number. Defaults to 1.

    Returns:
        tuple: data, lon_lat, size
    """
    ds = gdal.Open(file, gdal.GA_ReadOnly)
    data = ds.GetRasterBand(band).ReadAsArray()

    trans = ds.GetGeoTransform()
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    lon_lat = [
        trans[0], trans[0] + xsize * trans[1], trans[3] + ysize * trans[5], trans[3]
    ]
    ds = None
    xstep = trans[1]
    ystep = trans[5]
    size = [xsize, ysize, xstep, ystep]

    return data, lon_lat, size


def draw_img(data, lon_lat, out_file, cmap='rainbow'):
    """Draw image from array

    Args:
        data (array): numpy array
        lon_lat (list): lon_lat of data
        out_file (str): output file name
        cmap (str, optional): colormap. Defaults to 'rainbow'.
    """
    plt.figure()
    ax = plt.gca()
    im = ax.imshow(data, extent=lon_lat, cmap=cmap)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    plt.colorbar(im, cax=cax)
    plt.savefig(out_file, bbox_inches='tight', dpi=300)


def mosaic_tifs(tifs, out_tif):
    """Mosaic tifs

    Args:
        tifs (list): tifs
        out_tif (str): mosaiced tif
    """
    if os.path.isfile(out_tif):
        os.remove(out_tif)
    try:
        gdal.Warp(out_tif, tifs)
    except Exception as e:
        print(e)


def write_gamma_par(par_file, lon_lat, size, data_format):
    """Write GAMMA format parameter

    Args:
        par_file (str): output par file
        lon_lat (list): longitude and latitude
        size (list): xsize ysize xstep ystep
        data_format (str): data format
    """
    upper_left_lon, upper_left_lat = lon_lat[0], lon_lat[3]
    xstep, ystep = size[2], size[3]
    width, length = size[0], size[1]

    f = open(par_file, 'w')

    title = os.path.basename(par_file)[:-4]
    proj = 'EQA'

    f.write("Gamma DIFF&GEO DEM/MAP parameter file\n")
    f.write("title:\t%s\n" % title)
    # Projection should be checked
    f.write("DEM_projection:     %s\n" % proj)
    # INTEGER*2 (short) OR REAL*4 (float32) should be modified
    f.write("data_format:        %s\n" % data_format)
    f.write("DEM_hgt_offset:          0.00000\n")
    f.write("DEM_scale:               1.00000\n")
    f.write("width:                %s\n" % width)
    f.write("nlines:               %s\n" % length)
    f.write("corner_lat:   %s  decimal degrees\n" % upper_left_lat)
    f.write("corner_lon:   %s  decimal degrees\n" % upper_left_lon)
    f.write("post_lat:   %s  decimal degrees\n" % -abs(ystep))
    f.write("post_lon:   %s  decimal degrees\n" % abs(xstep))
    f.write("\n")
    f.write("ellipsoid_name: WGS 84\n")
    f.write("ellipsoid_ra:        6378137.000   m\n")
    f.write("ellipsoid_reciprocal_flattening:  298.2572236\n")
    f.write("\n")
    f.write("datum_name: WGS 1984\n")
    f.write("datum_shift_dx:              0.000   m\n")
    f.write("datum_shift_dy:              0.000   m\n")
    f.write("datum_shift_dz:              0.000   m\n")
    f.write("datum_scale_m:         0.00000e+00\n")
    f.write("datum_rotation_alpha:  0.00000e+00   arc-sec\n")
    f.write("datum_rotation_beta:   0.00000e+00   arc-sec\n")
    f.write("datum_rotation_gamma:  0.00000e+00   arc-sec\n")
    f.write("datum_country_list Global Definition, WGS84, World\n")
    f.close()


def make_dem(tif, out_file):
    """Make GAMMA dem

    Args:
        tif (str): tif file
        out_file (str): output dem file name

    Returns:
        tuple: dem_data, lon_lat, size
    """
    dem_data, lon_lat, size = read_gdal_file(tif)
    ori_dem_data = dem_data.copy()

    if dem_data.dtype == 'float32':
        data_format = 'REAL*4'
    else:
        data_format = 'INTEGER*2'

    dem_data.byteswap(True)

    dem_data.tofile(out_file)
    write_gamma_par(out_file + '.par', lon_lat, size, data_format)

    dem_data = None

    return ori_dem_data, lon_lat, size


def uncompress_zips(zip_files, out_dir):
    """Uncompress zips file

    Args:
        zip_files (list): zip files
        out_dir (str): output directory

    Returns:
        list: uncompressed tifs
    """
    tifs = []
    for zip_file in zip_files:
        with zipfile.ZipFile(zip_file, mode='r') as f:
            files = f.namelist()
            for file in files:
                if file.endswith('hgt') or file.endswith('tif'):
                    f.extract(file, out_dir)
                    tifs.append(os.path.join(out_dir, file))
    return tifs


EXAMPLE = '''Example:
  python3 make_gamma_dem.py ly.dem 1.hgt 2.hgt
  python3 make_gamma_dem.py ly.dem 1.tif 2.tif
  python3 make_gamma_dem.py ly.dem srtm_57_08.zip
  python3 make_gamma_dem.py ly.dem N20E100.SRTMGL1.hgt.zip N21E100.SRTMGL1.hgt.zip
'''


def cmdline_parser():
    parser = argparse.ArgumentParser(description='Make GAMMA DEM for interferometry.',
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)
    parser.add_argument('out_file', type=str, help='output dem file name')
    parser.add_argument('tifs', type=str, nargs='+', help='tif(hgt) or zips files for making DEM')
    inps = parser.parse_args()

    return inps


def main():
    # get inputs
    inps = cmdline_parser()
    tifs = [os.path.abspath(t) for t in inps.tifs]
    out_file = os.path.abspath(inps.out_file)

    # check tif files
    not_exist = []
    for tif in tifs:
        if not os.path.isfile(tif):
            not_exist.append(tif)
    if not_exist:
        for i in not_exist:
            print(f"{i} does not exist, please check it.")
        sys.exit()

    # check output file name
    dir_name = os.path.dirname(out_file)
    if not os.path.isdir(dir_name):
        sys.exit(f"{dir_name} does not exist, please check it.")

    if not out_file.endswith('.dem'):
        out_file += '.dem'

    # uncompress zips
    flag = False
    if tifs[0].endswith('zip'):
        tifs = uncompress_zips(tifs, os.path.dirname(tifs[0]))
        flag = True

    # check tifs
    if len(tifs) == 0:
        sys.exit('No tif(hgt) files.')

    tif_for_dem = None
    if len(tifs) > 1:
        print('Mosaic data.')
        mosaiced_tif = out_file.replace('.dem', '.tif')
        mosaic_tifs(tifs, mosaiced_tif)

        tif_for_dem = mosaiced_tif

    if len(tifs) == 1:
        tif_for_dem = tifs[0]

    print(f"Make dem.")
    data, lon_lat, _ = make_dem(tif_for_dem, out_file)

    print('Draw image.')
    draw_img(data, lon_lat, out_file.replace('.dem', '.png'))

    data = None

    # delete uncompressed tifs
    if flag:
        for tif in tifs:
            os.remove(tif)

    print('All done, enjoy it!\n')


if __name__ == "__main__":
    main()

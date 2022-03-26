#!/usr/bin/env python3
################################################################
# Run two-pass D-InSAR using GAMMA                             #
# APS correction using Python (support .ztd and .ztd.tif file) #
# Copyright (c) 2022, Lei Yuan                                 #
################################################################

import argparse
import glob
import os
import re
import sys

import numpy as np
import pyresample
from osgeo import gdal

dinsar_script = """#!/bin/bash
m_rslc=m_rslc_flag
m_par=$m_rslc.par

s_rslc=s_rslc_flag
s_par=$s_rslc.par

dem=dem_flag
dem_par=dem_par_flag

rlks=rlks_flag
alks=alks_flag

cc_thres=cc_thres_flag
roff=roff_flag
loff=loff_flag

m_date=m_date_flag
s_date=s_date_flag
pair=$m_date\_$s_date


echo -ne "$pair\\n 0 0\\n 32 32\\n 64 64\\n 7.0\\n 0\\n\\n" > off_par.in
create_offset $m_par $s_par $pair.off 1 1 1 < off_par.in
rm -f off_par.in

init_offset_orbit $m_par $s_par $pair.off
init_offset $m_rslc $s_rslc $m_par $s_par $pair.off 1 1

SLC_intf $m_rslc $s_rslc $m_par $s_par $pair.off $pair.int $rlks $alks - - 1 1

width_rdc=$(awk '$1 == "interferogram_width:" {print $2}' $pair.off)
line_rdc=$(awk '$1 == "interferogram_azimuth_lines:" {print $2}' $pair.off)

multi_look $m_rslc $m_par $m_date.mli $m_date.mli.par $rlks $alks
multi_look $s_rslc $s_par $s_date.mli $s_date.mli.par $rlks $alks

raspwr $m_date.mli $width_rdc 1 0 1 1 1. .35 1 $m_date.mli.bmp
raspwr $s_date.mli $width_rdc 1 0 1 1 1. .35 1 $s_date.mli.bmp

rasmph_pwr $pair.int $m_date.mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.int_pwr.bmp

base_init $m_par $s_par $pair.off $pair.int $pair.base 0 1024 1024
base_perp $pair.base $m_par $pair.off > $pair.base.perp

cc_wave $pair.int $m_date.mli $s_date.mli $pair.cor $width_rdc - - 3

post=$(awk '$1 == "post_lon:" {print $2}' $dem_par)
mk_geo $m_date.mli $m_date.mli.par $dem $dem_par dem_seg dem_seg.par . $m_date $post 0 2
mk_geo $m_date.mli $m_date.mli.par $dem $dem_par dem_seg dem_seg.par . $m_date $post 2 2
mk_geo $m_date.mli $m_date.mli.par $dem $dem_par dem_seg dem_seg.par . $m_date $post 3 2
mk_geo $m_date.mli $m_date.mli.par $dem $dem_par dem_seg dem_seg.par . $m_date $post 4 2

mv ${m_date}_1.map_to_rdc $m_date.lookup_fine

phase_sim $m_par $pair.off $pair.base ${m_date}_dem.rdc $pair.sim_unw 0 0 - -

sub_phase $pair.int $pair.sim_unw $m_date.diff_par $pair.diff 1 0
rasmph $pair.diff $width_rdc 1 0 1 1 1. 0.35 1 $pair.diff.bmp
rasmph_pwr $pair.diff $m_date.mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.diff_pwr.bmp

adf $pair.diff $pair.adf.diff1 $pair.adf.cc1 $width_rdc 0.3 128
adf $pair.adf.diff1 $pair.adf.diff2 $pair.adf.cc2 $width_rdc 0.3 64
adf $pair.adf.diff2 $pair.adf.diff $pair.adf.cc $width_rdc 0.3

rm -rf $pair.adf.diff1 $pair.adf.diff2 $pair.adf.cc1 $pair.adf.cc2

rasmph $pair.adf.diff $width_rdc 1 0 1 1 1. 0.35 1 $pair.adf.diff.bmp
rasmph_pwr $pair.adf.diff $m_date.mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.adf.diff_pwr.bmp

# rascc_mask $pair.adf.cc $m_date.mli $width_rdc 1 1 0 1 1 $cc_thres 0. .1 .9 1. .35 1 $pair.adf.cc_mask.bmp
rascc_mask $pair.cor $m_date.mli $width_rdc 1 1 0 1 1 $cc_thres 0. .1 .9 1. .35 1 $pair.adf.cc_mask.bmp

mcf $pair.adf.diff $pair.cor $pair.adf.cc_mask.bmp $pair.adf.unw $width_rdc 1 0 0 - - 1 1 - $roff $loff 0

rasrmg $pair.adf.unw $m_date.mli $width_rdc 1 1 0 1 1 .6 1. .35 .0 1 $pair.adf.unw_pwr.bmp $pair.adf.cc 1 .2
rasrmg $pair.adf.unw - $width_rdc 1 1 0 1 1 .5 1. .35 .0 1 $pair.adf.unw.bmp $pair.adf.cc 1 .2

quad_fit $pair.adf.unw $m_date.diff_par 32 32 $pair.adf.cc_mask.bmp $pair.plot 3
quad_sub $pair.adf.unw $m_date.diff_par $pair.adf.unw.sub 0 0

rasrmg $pair.adf.unw.sub -  $width_rdc 1 1 0 1 1 .5 1. .35 .0 1 $pair.adf.unw.sub.bmp $pair.adf.cc 1 .2
rasrmg $pair.adf.unw.sub $m_date.mli $width_rdc 1 1 0 1 1 .6 1. .35 .0 1 $pair.adf.unw.sub_pwr.bmp $pair.adf.cc 1 .2

width_geo=$(awk '$1 == "width:" {print $2}' dem_seg.par)
line_geo=$(awk '$1 == "nlines:" {print $2}' dem_seg.par)

geocode_back $pair.adf.cc $width_rdc $m_date.lookup_fine $pair.adf.cc.geo $width_geo $line_geo 1 0
geocode_back $pair.cor $width_rdc $m_date.lookup_fine $pair.cor.geo $width_geo $line_geo 1 0
geocode_back $pair.adf.diff $width_rdc $m_date.lookup_fine $pair.adf.diff.geo $width_geo $line_geo 1 1
geocode_back $pair.adf.unw $width_rdc $m_date.lookup_fine $pair.adf.unw.geo $width_geo $line_geo 1 0
geocode_back $pair.adf.unw.sub $width_rdc $m_date.lookup_fine $pair.adf.unw.sub.geo $width_geo $line_geo 1 0
geocode_back $m_date.mli $width_rdc $m_date.lookup_fine $m_date.mli.geo $width_geo $line_geo 1 0

rascc $pair.adf.cc.geo - $width_geo 1 1 0 1 1 .1 .9 1. .35 1 $pair.adf.cc.geo.bmp
rascc $pair.cor.geo - $width_geo 1 1 0 1 1 .1 .9 1. .35 1 $pair.cor.geo.bmp

rasmph $pair.adf.diff.geo $width_geo 1 0 1 1 1. 0.35 1 $pair.adf.diff.geo.bmp
rasmph_pwr $pair.adf.diff.geo $m_date.mli.geo $width_geo 1 1 0 1 1 1. 0.35 1 $pair.adf.diff.geo_pwr.bmp

rasrmg $pair.adf.unw.geo - $width_geo 1 1 0 1 1 .18 1. .35 .0 1 $pair.adf.unw.geo.bmp $pair.adf.cc.geo 1 .2
rasrmg $pair.adf.unw.geo $m_date.mli.geo $width_geo 1 1 0 1 1 .18 1. .35 .0 1 $pair.adf.unw.geo_pwr.bmp $pair.adf.cc.geo 1 .2

rasrmg $pair.adf.unw.sub.geo - $width_geo 1 1 0 1 1 .18 1. .35 .0 1 $pair.adf.unw.sub.geo.bmp $pair.adf.cc.geo 1 .2
rasrmg $pair.adf.unw.sub.geo $m_date.mli.geo $width_geo 1 1 0 1 1 .18 1. .35 .0 1 $pair.adf.unw.sub.geo_pwr.bmp $pair.adf.cc.geo 1 .2
"""


aps_correction_script="""#!/bin/bash
m_date=m_date_flag
s_date=s_date_flag

cc_thres=cc_thres_flag
roff=roff_flag
loff=loff_flag

pair=$m_date\_$s_date

width_rdc=$(awk '$1 == "interferogram_width:" {print $2}' $pair.off)
line_rdc=$(awk '$1 == "interferogram_azimuth_lines:" {print $2}' $pair.off)
width_geo=$(awk '$1 == "width:" {print $2}' dem_seg.par)
line_geo=$(awk '$1 == "nlines:" {print $2}' dem_seg.par)

geocode $m_date.lookup_fine $pair.gacos $width_geo $pair.gacos.rdc $width_rdc $line_rdc 1 0
raspwr $pair.gacos.rdc $width_rdc 1 0 1 1 1. .35 1 $pair.gacos.rdc.bmp

sub_phase $pair.diff $pair.gacos.rdc $m_date.diff_par $pair.diff.gacos 1 0

rasmph_pwr $pair.diff.gacos $m_date.mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.diff.gacos_pwr.bmp

adf $pair.diff.gacos $pair.adf.diff.gacos1 $pair.adf.cc.gacos1 $width_rdc 0.3 128
adf $pair.adf.diff.gacos1 $pair.adf.diff.gacos2 $pair.adf.cc.gacos2 $width_rdc 0.3 64
adf $pair.adf.diff.gacos2 $pair.adf.diff.gacos $pair.adf.cc.gacos $width_rdc 0.3

rm -rf $pair.adf.diff.gacos1 $pair.adf.diff.gacos2 $pair.adf.cc.gacos1 $pair.adf.cc.gacos2

rasmph $pair.adf.diff.gacos $width_rdc 1 0 1 1 1. 0.35 1 $pair.adf.diff.gacos.bmp
rasmph_pwr $pair.adf.diff.gacos $m_date.mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.adf.diff.gacos_pwr.bmp

# rascc_mask $pair.adf.cc.gacos $m_date.mli $width_rdc 1 1 0 1 1 $cc_thres 0. .1 .9 1. .35 1 $pair.adf.cc_mask.bmp
rascc_mask $pair.cor $m_date.mli $width_rdc 1 1 0 1 1 $cc_thres 0. .1 .9 1. .35 1 $pair.adf.cc_mask.bmp

mcf $pair.adf.diff.gacos $pair.adf.cc.gacos $pair.adf.cc_mask.bmp $pair.adf.unw.gacos $width_rdc 1 0 0 - - 1 1 - $roff $loff 0

rasrmg $pair.adf.unw.gacos $m_date.mli $width_rdc 1 1 0 1 1 .6 1. .35 .0 1 $pair.adf.unw.gacos_pwr.bmp $pair.adf.cc.gacos 1 .2
rasrmg $pair.adf.unw.gacos  - $width_rdc 1 1 0 1 1 .5 1. .35 .0 1 $pair.adf.unw.gacos.bmp $pair.adf.cc.gacos 1 .2

quad_fit $pair.adf.unw.gacos $m_date.diff_par 32 32 $pair.adf.cc_mask.bmp $pair.plot 3
quad_sub $pair.adf.unw.gacos $m_date.diff_par $pair.adf.unw.sub.gacos 0 0

rasrmg $pair.adf.unw.sub.gacos $m_date.mli $width_rdc 1 1 0 1 1 .6 1. .35 .0 1 $pair.adf.unw.sub.gacos_pwr.bmp $pair.adf.cc.gacos 1 .2
rasrmg $pair.adf.unw.sub.gacos -  $width_rdc 1 1 0 1 1 .5 1. .35 .0 1 $pair.adf.unw.sub.gacos.bmp $pair.adf.cc.gacos 1 .2

geocode_back $pair.adf.diff.gacos $width_rdc $m_date.lookup_fine $pair.adf.diff.gacos.geo $width_geo $line_geo 1 1
geocode_back $pair.adf.unw.gacos $width_rdc $m_date.lookup_fine $pair.adf.unw.gacos.geo $width_geo $line_geo 1 0
geocode_back $pair.adf.unw.sub.gacos $width_rdc $m_date.lookup_fine $pair.adf.unw.sub.gacos.geo $width_geo $line_geo 1 0


rasmph $pair.adf.diff.gacos.geo $width_geo 1 0 1 1 1. 0.35 1 $pair.adf.diff.gacos.geo.bmp
rasmph_pwr $pair.adf.diff.gacos.geo $m_date.mli.geo $width_geo 1 1 0 1 1 1. 0.35 1 $pair.adf.diff.gacos.geo_pwr.bmp

rasrmg $pair.adf.unw.gacos.geo - $width_geo 1 1 0 1 1 .18 1. .35 .0 1 $pair.adf.unw.gacos.geo.bmp $pair.adf.cc.geo 1 .2
rasrmg $pair.adf.unw.gacos.geo $m_date.mli.geo $width_geo 1 1 0 1 1 .18 1. .35 .0 1 $pair.adf.unw.gacos.geo_pwr.bmp $pair.adf.cc.geo 1 .2

rasrmg $pair.adf.unw.sub.gacos.geo - $width_geo 1 1 0 1 1 .18 1. .35 .0 1 $pair.adf.unw.sub.gacos.geo.bmp $pair.adf.cc.geo 1 .2
rasrmg $pair.adf.unw.sub.gacos.geo $m_date.mli.geo $width_geo 1 1 0 1 1 .18 1. .35 .0 1 $pair.adf.unw.sub.gacos.geo_pwr.bmp $pair.adf.cc.geo 1 .2
"""

USAGE = """Example:
  # reference point of unwrapping is (0, 0) and no mask
  python3 dinsar1.py /ly/rslc 20210110 20210122 /ly/dem DInSAR 8 2
  # set reference point of unwrapping and mask
  python3 dinsar1.py /ly/rslc 20210110 20210122 /ly/dem DInSAR 8 2 -r 100 -l 100 -c 0.3
  # for aps correction (Sentinel-1 data)
  python3 dinsar1.py /ly/rslc 20210110 20210122 /ly/dem DInSAR 8 2 -g /ly/gacos -w 0.05546
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(
        description='Run two-pass D-InSAR using GAMMA.\nAPS correction using Python (support .ztd and .ztd.tif file)',
        formatter_class=argparse.RawTextHelpFormatter, epilog=USAGE)

    parser.add_argument('rslc_dir', help='RSLCs directory')
    parser.add_argument('m_date', help='date for master RSLC')
    parser.add_argument('s_date', help='date for slave RSLC')
    parser.add_argument('dem_dir', help='DEM directory contains *.dem and *.dem.par')
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('rlks', help='range looks', type=int)
    parser.add_argument('alks', help='azimuth looks', type=int)
    parser.add_argument('-e', dest='rslc_extension', type=str, default='.rslc', help='file extension for RSLC (defaults: .rslc)')
    parser.add_argument('-g', dest='gacos_dir', help='directory contains GACOS files for aps correction')
    parser.add_argument('-w', dest='wavelength', type=float, help='Microwave length (Sentinel-1: 0.05546576, ALOS: 0.23830879)')
    parser.add_argument('-r', dest='roff', help='phase reference range offset to unwrap (defaults: 0)', type=float, default=0)
    parser.add_argument('-l', dest='loff', help='phase reference azimuth offset to unwrap (defaults: 0)', type=float, default=0)
    parser.add_argument('-c', dest='cc_thres', type=float, default=0, help='threshold of correlation for creating the unwrapping mask (0.0 --> 1.0) (defaults: 0)')

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


def read_dem_par(dem_seg_par):
    """Extract value from dem_seg_par file

    Args:
        dem_seg_par (str): dem_seg.par file

    Returns:
        list: values
    """
    par_key = [
        'width', 'nlines', 'corner_lat', 'corner_lon', 'post_lat', 'post_lon'
    ]
    par_value = []
    with open(dem_seg_par, 'r') as f:
        for line in f.readlines():
            for i in par_key:
                if line.strip().startswith(i):
                    value = line.strip().split()[1]
                    if i in ['width', 'nlines']:
                        par_value.append(int(value))
                    else:
                        par_value.append(float(value))

    return par_value


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


def read_binary_file(file, rsc_file):
    """Read binary file
    Args:
        file (str): binary file
        rsc_file (str): parameter for binary file

    Returns:
        tuple: data, lon_lat, size
    """
    data = np.fromfile(file, dtype='float32')
    width = int(read_gamma_par(rsc_file, 'WIDTH'))
    length = int(read_gamma_par(rsc_file, 'FILE_LENGTH'))
    data = data.reshape((length, width))

    x_first = float(read_gamma_par(rsc_file, 'X_FIRST'))
    y_first = float(read_gamma_par(rsc_file, 'Y_FIRST'))
    x_step = float(read_gamma_par(rsc_file, 'X_STEP'))
    y_step = float(read_gamma_par(rsc_file, 'Y_STEP'))

    lon_lat = [x_first, x_first + width * x_step, y_first + length * y_step, y_first]
    size = [width, length, x_step, y_step]

    return data, lon_lat, size


def write_gamma(data, out_file, file_type):
    """Write GAMMA format file

    Args:
        data (array): array
        out_file (str): output file name
        file_type (str): file type
    """
    data = data.astype(file_type)
    data.byteswap('True')
    data.reshape(-1, 1)
    data.tofile(out_file)


def process_gacos(m_gacos, s_gacos, gacos_type, dem_seg_par, wavelength, incidence, out_file):
    """Process gacos aps data

    Args:
        m_gacos (str): gacos file for master date
        s_gacos (str): gacos file for slave date
        gacos_type (str): tif for geotiff, bin for binary file
        dem_seg_par (str): dem_seg.par file
        wavelength (float): microwave length (m)
        incidence (float): incidence angle (rad)
        out_file (str): output file
    """
    values = read_dem_par(dem_seg_par)
    width, length = values[0], values[1]
    upper_left_lon, upper_left_lat = values[3], values[2]
    lat_step, lon_step = values[4], values[5]

    lon = np.linspace(upper_left_lon, upper_left_lon + lon_step * width, width)
    lat = np.linspace(upper_left_lat, upper_left_lat + lat_step * length, length)

    lons, lats = np.meshgrid(lon, lat)

    if gacos_type == 'tif':
        m_ztd, lon_lat, size = read_gdal_file(m_gacos)
        s_ztd, _, _ = read_gdal_file(s_gacos)
    if gacos_type == 'bin':
        m_ztd, lon_lat, size = read_binary_file(m_gacos, m_gacos + '.rsc')
        s_ztd, _, _ = read_binary_file(s_gacos, s_gacos + '.rsc')
    width, length = size[0], size[1]

    lon = np.linspace(lon_lat[0], lon_lat[1], width)
    lat = np.linspace(lon_lat[3], lon_lat[2], length)

    lons_gacos, lats_gacos = np.meshgrid(lon, lat)

    ori_grid = pyresample.geometry.SwathDefinition(lons=lons_gacos, lats=lats_gacos)
    dst_grid = pyresample.geometry.SwathDefinition(lons=lons, lats=lats)

    diff_ztd = m_ztd - s_ztd

    radius = abs(lat_step) * np.pi / 180.0 * 6378122.65 * 3

    diff_ztd_interp = pyresample.kd_tree.resample_gauss(ori_grid, diff_ztd, dst_grid, radius_of_influence=radius, sigmas=25000)

    range2phase = -4 * np.pi / wavelength * np.cos(incidence)

    diff_ztd_interp_los = diff_ztd_interp * range2phase

    write_gamma(diff_ztd_interp_los, out_file, 'float32')


def get_gacos_extension(gacos_dir):
    """Get gacos file extension

    Args:
        gacos_dir (str): gacos directory

    Returns:
        str: file extension
    """
    if glob.glob(os.path.join(gacos_dir, '*.ztd')):
        return '.ztd'
    elif glob.glob(os.path.join(gacos_dir, '*.ztd.tif')):
        return '.ztd.tif'
    else:
        return None


def check_gacos(dates, gacos_dir, file_extension):
    """Check gacos files

    Args:
        dates (list): RSLC dates
        gacos_dir (str): gacos directory
        file_extension (str): file extension for gacos
    """
    not_exist = []
    if not os.path.isdir(gacos_dir):
        sys.exit('{} does not exist.'.format(gacos_dir))
    for date in dates:
        file = os.path.join(gacos_dir, date + file_extension)
        if not os.path.isfile(file):
            not_exist.append(file)
    if not_exist:
        for i in not_exist:
            print('{} does not exist.'.format(i))
        sys.exit()


if __name__ == "__main__":
    # get inputs
    inps = cmdline_parser()
    rslc_dir = os.path.abspath(inps.rslc_dir)
    m_date = inps.m_date
    s_date = inps.s_date
    dem_dir = os.path.abspath(inps.dem_dir)
    out_dir = os.path.abspath(inps.out_dir)
    rlks = inps.rlks
    alks = inps.alks

    rslc_extension = inps.rslc_extension

    gacos_dir = inps.gacos_dir
    wavelength = inps.wavelength

    roff = inps.roff
    loff = inps.loff
    cc_thres = inps.cc_thres

    # check rslc_dir
    if not os.path.isdir(rslc_dir):
        sys.exit('{} does not exist.'.format(rslc_dir))
    dates = sorted([i for i in os.listdir(rslc_dir) if re.match(r'^\d{8}$', i)])
    if not dates:
        sys.exit('Cannot find RSLCs in {}'.format(rslc_dir))

    # check RSLC extension
    if not rslc_extension.startswith('.'):
        rslc_extension = '.' + rslc_extension

    # check m_date rslc and s_date rslc
    m_rslc = os.path.join(rslc_dir, m_date, m_date + rslc_extension)
    s_rslc = os.path.join(rslc_dir, s_date, s_date + rslc_extension)
    if not os.path.isfile(m_rslc):
        sys.exit('Cannot find RSLC for {}'.format(m_date))
    if s_date not in dates:
        sys.exit('Cannot find RSLC for {}'.format(s_date))

    # check out_dir
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # check dem_dir
    if not os.path.isdir(dem_dir):
        sys.exit("{} does not exist.".format(dem_dir))
    else:
        dems = glob.glob(dem_dir + '/*.dem')
        dem_pars = [i + '.par' for i in dems]
        for i, j in zip(dems, dem_pars):
            if os.path.isfile(i) and os.path.isfile(j):
                dem = dems[0]
                dem_par = dem + '.par'
            else:
                sys.exit(f'Cannot find *.dem and *.dem.par in {dem_dir}.')

    # check gacos_dir
    if gacos_dir and wavelength:
        gacos_dir = os.path.abspath(gacos_dir)
        if not os.path.isdir(gacos_dir):
            sys.exit('{} does not exist.'.format(gacos_dir))

        gacos_extension = get_gacos_extension(gacos_dir)
        if gacos_extension is None:
            sys.exit(f'Cannot find *.ztd or *.ztd.tif files in {gacos_dir}')
        check_gacos([m_date, s_date], gacos_dir, gacos_extension)

    # DInSAR
    # replace value
    dinsar_script_out = dinsar_script.replace('m_date_flag', m_date)
    dinsar_script_out = dinsar_script_out.replace('s_date_flag', s_date)
    dinsar_script_out = dinsar_script_out.replace('m_rslc_flag', m_rslc)
    dinsar_script_out = dinsar_script_out.replace('s_rslc_flag', s_rslc)
    dinsar_script_out = dinsar_script_out.replace('dem_flag', dem)
    dinsar_script_out = dinsar_script_out.replace('dem_par_flag', dem_par)
    dinsar_script_out = dinsar_script_out.replace('rlks_flag', str(rlks))
    dinsar_script_out = dinsar_script_out.replace('alks_flag', str(alks))
    dinsar_script_out = dinsar_script_out.replace('cc_thres_flag', str(cc_thres))
    dinsar_script_out = dinsar_script_out.replace('roff_flag', str(roff))
    dinsar_script_out = dinsar_script_out.replace('loff_flag', str(loff))

    # write script
    os.chdir(out_dir)
    run_script = m_date + '_' + s_date + '_DInSAR.sh'
    with open(run_script, 'w+') as f:
        f.write(dinsar_script_out)

    # run script
    call_str = 'bash ' + run_script
    os.system(call_str)

    # aps correction
    if gacos_dir and wavelength:
        gacos_dir = os.path.abspath(gacos_dir)
        gacos_extension = get_gacos_extension(gacos_dir)

        os.chdir(out_dir)

        m_gacos = os.path.join(gacos_dir, m_date + gacos_extension)
        s_gacos = os.path.join(gacos_dir, s_date + gacos_extension)
        incidence = read_gamma_par(m_rslc + '.par', 'incidence_angle')
        incidence = np.deg2rad(float(incidence))
        out_file = m_date + '_' + s_date + '.gacos'

        if gacos_extension == '.ztd':
            gacos_type = 'bin'
        if gacos_extension == '.ztd.tif':
            gacos_type = 'tif'
        process_gacos(m_gacos, s_gacos, gacos_type, 'dem_seg.par', wavelength, incidence, out_file)

        # replace value
        aps_correction_script_out = aps_correction_script.replace('m_date_flag', m_date)
        aps_correction_script_out = aps_correction_script_out.replace('s_date_flag', s_date)
        aps_correction_script_out = aps_correction_script_out.replace('cc_thres_flag', str(cc_thres))
        aps_correction_script_out = aps_correction_script_out.replace('roff_flag', str(roff))
        aps_correction_script_out = aps_correction_script_out.replace('loff_flag', str(loff))

        # write script
        run_script = m_date + '_' + s_date + '_aps_correction.sh'
        with open(run_script, 'w+') as f:
            f.write(aps_correction_script_out)

        # run script
        call_str = 'bash ' + run_script
        os.system(call_str)

    print('\nAll done, enjoy it!\n')

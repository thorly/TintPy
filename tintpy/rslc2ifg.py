#!/usr/bin/env python3
###############################################################################
# Generate differential interferogram and unwrap them from RSLCs using GAMMA  #
# Copyright (c) 2022, Lei Yuan                                                #
###############################################################################

import argparse
import glob
import os
import re
import sys

import numpy as np
import pyresample
from osgeo import gdal

dinsar_script = """m_rslc=$rslc_dir/$m_date/$m_date.rslc
s_rslc=$rslc_dir/$s_date/$s_date.rslc
m_par=$m_rslc.par
s_par=$s_rslc.par

pair=$m_date_$s_date
off_par=$pair.off

echo -ne "$pair\\n 0 0\\n 32 32\\n 64 64\\n 7.0\\n 0\\n\\n" > off_par.in
create_offset $m_par $s_par $off_par 1 1 1 < off_par.in
rm -f off_par.in

init_offset_orbit $m_par $s_par $off_par
init_offset $m_rslc $s_rslc $m_par $s_par $off_par 1 1

SLC_intf $m_rslc $s_rslc $m_par $s_par $off_par $pair.int $rlks $alks - - 1 1

width_rdc=$(awk '$1 == "interferogram_width:" {print $2}' $off_par)
line_rdc=$(awk '$1 == "interferogram_azimuth_lines:" {print $2}' $off_par)

rasmph_pwr $pair.int $m_mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.int_pwr.bmp


base_init $m_par $s_par $off_par $pair.int $pair.base 0 1024 1024
base_perp $pair.base $m_par $off_par > $pair.base.perp


cc_wave $pair.int $m_par $s_par $pair.cor $width_rdc - - 3

gc_map $m_par $off_par $dem_par $dem $pair.dem_seg.par $pair.dem_seg $pair.lookup 1 1 $pair.sim_sar - - - - - - 8 1

width_geo=$(awk '$1 == "width:" {print $2}' $pair.dem_seg.par)
line_geo=$(awk '$1 == "nlines:" {print $2}' $pair.dem_seg.par)
col_post=$(awk '$1 == "post_lat:" {print $2}' $pair.dem_seg.par)
row_post=$(awk '$1 == "post_lon:" {print $2}' $pair.dem_seg.par)

rasshd $pair.dem_seg $width_geo $col_post $row_post 1 0 1 1 45. 135. 1 $pair.dem_seg.bmp

geocode $pair.lookup $pair.sim_sar $width_geo $pair.sim_sar_rdc $width_rdc $line_rdc 1 0
raspwr $pair.sim_sar_rdc $width_rdc 1 0 1 1 1. .35 1 $pair.sim_sar_rdc.bmp

echo -ne "$pair\\n 0 0\\n 64 64\\n 256 256\\n 7.0\\n" > diff_par.in
create_diff_par $off_par - $pair.diff_par 0 < diff_par.in
rm -f diff_par.in

init_offsetm $m_mli $pair.sim_sar_rdc $pair.diff_par 1 1 - - 0 0 7
offset_pwrm $m_mli $pair.sim_sar_rdc $pair.diff_par $pair.offs $pair.snr 256 256 $pair.offsets 2 100 100 7.0 2
offset_fitm $pair.offs $pair.snr $pair.diff_par $pair.coffs $pair.coffsets 8.0 6
offset_pwrm $m_mli $pair.sim_sar_rdc $pair.diff_par $pair.offs $pair.snr 64 64 $pair.offsets 2 300 300 9.0 2
offset_fitm $pair.offs $pair.snr $pair.diff_par $pair.coffs $pair.coffsets 10.0 6

gc_map_fine $pair.lookup $width_geo $pair.diff_par $pair.lookup_fine 0
geocode $pair.lookup_fine $pair.dem_seg $width_geo $pair.rdc_hgt $width_rdc $line_rdc 1 0
rashgt $pair.rdc_hgt $m_mli $width_rdc 1 1 0 1 1 20.0 1. .35 1 $pair.rdc_hgt_pwr.bmp


phase_sim $m_par $off_par $pair.base $pair.rdc_hgt $pair.sim_unw 0 0 - -

sub_phase $pair.int $pair.sim_unw $pair.diff_par $pair.diff 1 0
rasmph $pair.diff $width_rdc 1 0 1 1 1. 0.35 1 $pair.diff.bmp
rasmph_pwr $pair.diff $m_mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.diff_pwr.bmp

adf $pair.diff $pair.adf.diff1 $pair.adf.cc1 $width_rdc 0.3 128
adf $pair.adf.diff1 $pair.adf.diff2 $pair.adf.cc2 $width_rdc 0.3 64
adf $pair.adf.diff2 $pair.adf.diff $pair.adf.cc $width_rdc 0.3

rm -rf $pair.adf.diff1 $pair.adf.diff2 $pair.adf.cc1 $pair.adf.cc2

rasmph $pair.adf.diff $width_rdc 1 0 1 1 1. 0.35 1 $pair.adf.diff.bmp
rasmph_pwr $pair.adf.diff $m_mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.adf.diff_pwr.bmp


# rascc_mask $pair.adf.cc $m_mli $width_rdc 1 1 0 1 1 $cc_thres 0. .1 .9 1. .35 1 $pair.adf.cc_mask.bmp
rascc_mask $pair.cor $m_mli $width_rdc 1 1 0 1 1 $cc_thres 0. .1 .9 1. .35 1 $pair.adf.cc_mask.bmp

mcf $pair.adf.diff $pair.cor $pair.adf.cc_mask.bmp $pair.adf.unw $width_rdc 1 0 0 - - 1 1 - $roff $loff 0

rasrmg $pair.adf.unw $m_mli $width_rdc 1 1 0 1 1 .6 1. .35 .0 1 $pair.adf.unw_pwr.bmp $pair.adf.cc 1 .2
rasrmg $pair.adf.unw - $width_rdc 1 1 0 1 1 .5 1. .35 .0 1 $pair.adf.unw.bmp $pair.adf.cc 1 .2


quad_fit $pair.adf.unw $pair.diff_par 32 32 $pair.adf.cc_mask.bmp $pair.plot 3
quad_sub $pair.adf.unw $pair.diff_par $pair.adf.unw.sub 0 0

rasrmg $pair.adf.unw.sub $m_mli $width_rdc 1 1 0 1 1 .6 1. .35 .0 1 $pair.adf.unw.sub_pwr.bmp $pair.adf.cc 1 .2
rasrmg $pair.adf.unw.sub -  $width_rdc 1 1 0 1 1 .5 1. .35 .0 1 $pair.adf.unw.sub.bmp $pair.adf.cc 1 .2


# geocode_back $pair.adf.cc $width_rdc $pair.lookup_fine $pair.adf.cc.geo $width_geo $line_geo 1 0
# geocode_back $pair.adf.diff $width_rdc $pair.lookup_fine $pair.adf.diff.geo $width_geo $line_geo 1 1
# geocode_back $pair.adf.unw $width_rdc $pair.lookup_fine $pair.adf.unw.geo $width_geo $line_geo 1 0
# geocode_back $m_mli $width_rdc $pair.lookup_fine $pair.mli.geo $width_geo $line_geo 1 0

# rascc $pair.adf.cc.geo - $width_geo 1 1 0 1 1 .1 .9 1. .35 1 $pair.adf.cc.geo.bmp
# rascc $pair.adf.cc.geo $pair.mli.geo $width_geo 1 1 0 1 1 .1 .9 1. .35 1 $pair.adf.cc.geo_pwr.bmp
# rasmph $pair.adf.diff.geo $width_geo 1 0 1 1 1. 0.35 1 $pair.adf.diff.geo.bmp
# rasmph_pwr $pair.adf.diff.geo $pair.mli.geo $width_geo 1 1 0 1 1 1. 0.35 1 $pair.adf.diff.geo_pwr.bmp
# rasrmg $pair.adf.unw -  $width_rdc 1 1 0 1 1 .18 1. .35 .0 1 $pair.diff.unw.bmp $pair.adf.cc 1 .2
"""

aps_correction_script="""off_par=$pair.off
width_rdc=$(awk '$1 == "interferogram_width:" {print $2}' $off_par)
line_rdc=$(awk '$1 == "interferogram_azimuth_lines:" {print $2}' $off_par)
width_geo=$(awk '$1 == "width:" {print $2}' $pair.dem_seg.par)

geocode $pair.lookup_fine $pair.gacos $width_geo $pair.gacos.rdc $width_rdc $line_rdc 1 0
raspwr $pair.gacos.rdc $width_rdc 1 0 1 1 1. .35 1 $pair.gacos.rdc.bmp

sub_phase $pair.diff $pair.gacos.rdc $pair.diff_par $pair.diff.gacos 1 0

rasmph_pwr $pair.diff.gacos $m_mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.diff.gacos_pwr.bmp

adf $pair.diff.gacos $pair.adf.diff.gacos1 $pair.adf.cc.gacos1 $width_rdc 0.3 128
adf $pair.adf.diff.gacos1 $pair.adf.diff.gacos2 $pair.adf.cc.gacos2 $width_rdc 0.3 64
adf $pair.adf.diff.gacos2 $pair.adf.diff.gacos $pair.adf.cc.gacos $width_rdc 0.3

rm -rf $pair.adf.diff.gacos1 $pair.adf.diff.gacos2 $pair.adf.cc.gacos1 $pair.adf.cc.gacos2

rasmph $pair.adf.diff.gacos $width_rdc 1 0 1 1 1. 0.35 1 $pair.adf.diff.gacos.bmp
rasmph_pwr $pair.adf.diff.gacos $m_mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.adf.diff.gacos_pwr.bmp

rascc_mask $pair.adf.cc.gacos $m_mli $width_rdc 1 1 0 1 1 $cc_thres 0. .1 .9 1. .35 1 $pair.adf.cc.gacos_mask.bmp

mcf $pair.adf.diff.gacos $pair.adf.cc.gacos $pair.adf.cc.gacos_mask.bmp $pair.adf.unw.gacos $width_rdc 1 0 0 - - 1 1 - $roff $loff 0

rasrmg $pair.adf.unw.gacos $m_mli $width_rdc 1 1 0 1 1 .6 1. .35 .0 1 $pair.adf.unw.gacos_pwr.bmp $pair.adf.cc.gacos 1 .2
rasrmg $pair.adf.unw.gacos  - $width_rdc 1 1 0 1 1 .5 1. .35 .0 1 $pair.adf.unw.gacos.bmp $pair.adf.cc.gacos 1 .2

quad_fit $pair.adf.unw.gacos $pair.diff_par 32 32 $pair.adf.cc.gacos_mask.bmp $pair.plot 3
quad_sub $pair.adf.unw.gacos $pair.diff_par $pair.adf.unw.sub.gacos 0 0

rasrmg $pair.adf.unw.sub.gacos $m_mli $width_rdc 1 1 0 1 1 .6 1. .35 .0 1 $pair.adf.unw.sub.gacos_pwr.bmp $pair.adf.cc.gacos 1 .2
rasrmg $pair.adf.unw.sub.gacos -  $width_rdc 1 1 0 1 1 .5 1. .35 .0 1 $pair.adf.unw.sub.gacos.bmp $pair.adf.cc.gacos 1 .2
"""

USAGE = """Example:
  # reference point of unwrapping is (0, 0) and no mask
  python3 rslc2ifg.py /ly/rslc /ly/SBAS /ly/dem 20211229 8 2 200 60
  # set reference point of unwrapping and mask
  python3 rslc2ifg.py /ly/rslc /ly/SBAS /ly/dem 20211229 8 2 200 60 -r 100 -l 100 -c 0.3
  # for aps correction (Sentinel-1 data)
  python3 rslc2ifg.py /ly/rslc /ly/SBAS /ly/dem 20211229 8 2 200 60 -g /ly/gacos -w 0.05546
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(
        description='Generate differential interferogram and unwrap them from RSLCs using GAMMA.',
        formatter_class=argparse.RawTextHelpFormatter, epilog=USAGE)

    parser.add_argument('rslc_dir', help='RSLCs directory')
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('dem_dir', help='DEM directory contains *.dem and *.dem.par')
    parser.add_argument('ref_rslc', help='date of reference RSLC for calculating baseline')
    parser.add_argument('rlks', help='range looks', type=int)
    parser.add_argument('alks', help='azimuth looks', type=int)
    parser.add_argument('max_sb', type=float, help='maximum spatial baseline')
    parser.add_argument('max_tb', type=float, help='maximum temporal baseline')
    parser.add_argument('-e', dest='rslc_extension', type=str, default='.rslc', help='file extension for RSLCs (defaults: .rslc)')
    parser.add_argument('-g', dest='gacos_dir', help='directory contains GACOS files for aps correction')
    parser.add_argument('-w', dest='wavelength', type=float, help='Microwave length (Sentinel-1: 0.05546576, ALOS: 0.23830879)')
    parser.add_argument('-r', dest='roff', help='phase reference range offset to unwrap (defaults: 0)', type=float, default=0)
    parser.add_argument('-l', dest='loff', help='phase reference azimuth offset to unwrap (defaults: 0)', type=float, default=0)
    parser.add_argument('-c', dest='cc_thres', type=float, default=0, help='threshold of correlation for creating the unwrapping mask (0.0 --> 1.0) (defaults: 0)')

    inps = parser.parse_args()
    return inps


def mk_tab(slc_dir, slc_tab, slc_extension):
    """Generate SLC_tab for processing

    Args:
        slc_dir (str): slc directory
        slc_tab (str): tab file
        slc_extension (str): slc extension
    """
    dates = sorted(
        [i for i in os.listdir(slc_dir) if re.findall(r'^\d{8}$', i)])
    with open(slc_tab, 'w+', encoding='utf-8') as f:
        for date in dates:
            slc = os.path.join(slc_dir, date, date + slc_extension)
            slc_par = slc + '.par'
            f.write(slc + '    ' + slc_par + '\n')


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


def mli_all(slc_tab, out_dir, rlks, alks):
    """Calculate MLI images for a stack of SLCs

    Args:
        slc_tab (str): slc tab file including slc slc_par
        out_dir (str): output directory
        rlks (int): range looks
        alks (int): azimuth looks
    """
    os.chdir(out_dir)
    with open(slc_tab, 'r') as f:
        for line in f.readlines():
            if line.strip():
                slc = line.strip().split()[0]
                slc_par = line.strip().split()[1]

                date = os.path.basename(slc)[0:8]

                mli = date + '.rmli'
                mli_par = mli + '.par'

                call_str = f"multi_look {slc} {slc_par} {mli} {mli_par} {rlks} {alks}"
                os.system(call_str)

                width = read_gamma_par(mli_par, 'range_samples')

                call_str = f"raspwr {mli} {width} 1 0 1 1 1. .35 1"
                os.system(call_str)


def base_calc(slc_tab, slc_par, max_sb, max_tb, out_dir):
    """Generate baseline output file with perpendicular baselines and delta_T values

    Args:
        slc_tab (str): slc tab file including slc and slc_par
        slc_par (str): reference slc par
        max_sb (float): max spatial baseline
        max_tb (float): max time baseline
        out_dir (str): output directory

    Returns:
        str: baseline file
    """
    os.chdir(out_dir)
    bperp_file = os.path.join(out_dir, 'bperp_file')
    call_str = f'base_calc {slc_tab} {slc_par} {bperp_file} itab 1 1 0 {max_sb} 0 {max_tb}'
    os.system(call_str)

    return bperp_file


def get_pairs(bperp_file):
    """Get pairs from baseline file

    Args:
        bperp_file (str): baseline file

    Returns:
        list: pairs
    """
    pairs = []
    with open(bperp_file, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            if line:
                split_list = line.strip().split()
                date1 = split_list[1]
                date2 = split_list[2]
                if int(date1) > int(date2):
                    pairs.append(date2 + '_' + date1)
                else:
                    pairs.append(date1 + '_' + date2)

    return pairs


def select_pairs_sbas(slc_tab, slc_par, max_sb, max_tb, out_dir):
    """Select pairs using sbas method

    Args:
        slc_tab (str): slc tab file including slc and slc_par
        slc_par (str): reference slc par
        max_sb (float): max spatial baseline
        max_tb (float): max time baseline
        out_dir (str): output directory

    Returns:
        list: pairs
    """
    bperp_file = base_calc(slc_tab, slc_par, max_sb, max_tb, out_dir)

    pairs = get_pairs(bperp_file)

    return pairs


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

    diff_ztd_interp = pyresample.kd_tree.resample_gauss(ori_grid, diff_ztd, dst_grid, radius_of_influence=500000, sigmas=25000)

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


def main():
    # get inputs
    inps = cmdline_parser()
    rslc_dir = os.path.abspath(inps.rslc_dir)
    out_dir = os.path.abspath(inps.out_dir)
    dem_dir = os.path.abspath(inps.dem_dir)
    ref_rslc = inps.ref_rslc
    rlks = inps.rlks
    alks = inps.alks
    max_sb = inps.max_sb
    max_tb = inps.max_tb
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

    # check ref_slc
    if re.findall(r'^\d{8}$', ref_rslc):
        if not ref_rslc in dates:
            sys.exit('No RSLC for {}.'.format(ref_rslc))
    else:
        sys.exit('Error date for ref_rslc.')

    # check RSLC extension
    if not rslc_extension.startswith('.'):
        rslc_extension = '.' + rslc_extension
    for date in dates:
        rslc_path = os.path.join(rslc_dir, date, date + rslc_extension)
        if not os.path.isfile(rslc_path):
            sys.exit('Error file extension for RSLC.')

    # multi-look
    mli_dir = os.path.join(out_dir, 'mli')
    if not os.path.isdir(mli_dir):
        os.mkdir(mli_dir)

    rslc_tab = os.path.join(mli_dir, 'rslc_tab')
    mk_tab(rslc_dir, rslc_tab, rslc_extension)
    mli_all(rslc_tab, mli_dir, rlks, alks)

    # select pairs
    base_dir = os.path.join(out_dir, 'base_calc')
    if not os.path.isdir(base_dir):
        os.mkdir(base_dir)

    m_rslc_par = os.path.join(rslc_dir, ref_rslc, ref_rslc + rslc_extension + '.par')
    pairs = select_pairs_sbas(rslc_tab, m_rslc_par, max_sb, max_tb, base_dir)

    # diff and unwrap
    diff_dir = os.path.join(out_dir, 'diff')
    if not os.path.isdir(diff_dir):
        os.mkdir(diff_dir)
    os.chdir(diff_dir)

    for pair in pairs:
        m_date = pair[0:8]
        s_date = pair[9:17]

        # write D-InSAR script
        str_m_date = f"m_date={m_date}\n"
        str_s_date = f"s_date={s_date}\n"
        str_rslc_dir = f"rslc_dir={rslc_dir}\n"
        str_dem = f"dem={dem}\n"
        str_dem_par = f"dem_par={dem_par}\n"
        str_m_mli = f"m_mli={os.path.join(mli_dir, m_date + '.rmli')}\n"
        str_rlks = f"rlks={rlks}\n"
        str_alks = f"alks={alks}\n"
        str_cc_thres = f"cc_thres={cc_thres}\n"
        str_roff = f"roff={roff}\n"
        str_loff = f"loff={loff}\n"
        out_script = '#!/bin/bash\n' + str_m_date + str_s_date
        out_script += str_rslc_dir + str_dem + str_m_mli
        out_script += str_dem_par + str_rlks + str_alks
        out_script += str_cc_thres + str_roff + str_loff + dinsar_script
        out_script_path = os.path.join(diff_dir, pair + '_DInSAR.sh')
        with open(out_script_path, 'w+') as f:
            f.write(out_script)

        # run D-InSAR script
        call_str = 'bash ' + pair + '_DInSAR.sh'
        os.system(call_str)

    # aps correction
    if gacos_dir and wavelength:
        gacos_dir = os.path.abspath(gacos_dir)
        gacos_extension = get_gacos_extension(gacos_dir)
        if gacos_extension is None:
            sys.exit(f'Cannot find *.ztd or *.ztd.tif files in {gacos_dir}')
        check_gacos(dates, gacos_dir, gacos_extension)
        for pair in pairs:
            m_date = pair[0:8]
            s_date = pair[9:17]
            m_gacos = os.path.join(gacos_dir, m_date + gacos_extension)
            s_gacos = os.path.join(gacos_dir, s_date + gacos_extension)
            dem_seg_par = os.path.join(diff_dir, pair + '.dem_seg.par')
            incidence = read_gamma_par(m_rslc_par, 'incidence_angle')
            incidence = np.deg2rad(float(incidence))
            out_file = pair + '.gacos'
            if gacos_extension == '.ztd':
                process_gacos(m_gacos, s_gacos, 'bin', dem_seg_par, wavelength, incidence, out_file)
            if gacos_extension == '.ztd.tif':
                process_gacos(m_gacos, s_gacos, 'tif', dem_seg_par, wavelength, incidence, out_file)

            # write aps correction script
            str_pair = f"pair={pair}\n"
            str_m_mli = f"m_mli={os.path.join(mli_dir, m_date + '.rmli')}\n"
            str_cc_thres = f"cc_thres={cc_thres}\n"
            str_roff = f"roff={roff}\n"
            str_loff = f"loff={loff}\n"
            out_script = '#!/bin/bash\n' + str_pair + str_m_mli + str_cc_thres
            out_script += str_roff + str_loff + aps_correction_script
            out_script_path = os.path.join(diff_dir, pair + '_aps_correction.sh')
            with open(out_script_path, 'w+') as f:
                f.write(out_script)

            # run aps correction script
            call_str = 'bash ' + pair + '_aps_correction.sh'
            os.system(call_str)


    print('\nAll done, enjoy it!\n')

if __name__ == "__main__":
    main()

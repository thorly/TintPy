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
import sys

import numpy as np
import pyresample
from osgeo import gdal
import datetime

EXAMPLE = """Example:
  python3 rslc2ifg.py /ly/rslc /ly/stacking /ly/dem 20211229 20 5 -s 200 -t 60
  python3 rslc2ifg.py /ly/rslc /ly/stacking /ly/dem 20211229 20 5 -m sequential -n 4
  python3 rslc2ifg.py /ly/rslc /ly/stacking /ly/dem 20211229 20 5 -s 200 -t 60 -g /ly/gacos -w 0.05546576
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(description= 'Generate differential interferogram images and unwrap them from RSLCs using GAMMA.',
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)

    parser.add_argument('rslc_dir', help='RSLCs directory')
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('dem_dir', help='DEM directory contains *.dem and *.dem.par')
    parser.add_argument('ref_slc', help='reference RSLC for making rdc dem')
    parser.add_argument('rlks', help='range looks', type=int)
    parser.add_argument('alks', help='azimuth looks', type=int)
    parser.add_argument('-s', dest='max_sb', type=float, help='maximum spatial baseline')
    parser.add_argument('-t', dest='max_tb', type=float, help='maximum temporal baseline')
    parser.add_argument('-m', dest='method', default='sbas', choices={'sbas', 'sequential'},
        help='network selection method:\n' +
        'sbas - select based on the threshold values of the spatio-temporal baselines\n'
        +
        'sequential - select based on the sequential of the SAR acquisitions\n')
    parser.add_argument('-n', dest='con_num', type=int, default=4,
        help='Number of the neibour-connected SAR images at one side for sequential method (defaults: 4)')
    parser.add_argument('-e', dest='slc_extension', type=str, default='.rslc', help='file extension for RSLCs (defaults: .rslc)')
    parser.add_argument('-g', dest='gacos_dir', help='directory contains GACOS files (tif) for aps correction')
    parser.add_argument('-w', dest='wavelength', type=float, help='Microwave length (Sentinel-1: 0.05546576, ALOS: 0.23830879)')
    parser.add_argument('-r', dest='roff', help='phase reference range offset to unwrap (defaults: 0)', type=float, default=0)
    parser.add_argument('-l', dest='loff', help='phase reference azimuth offset to unwrap (defaults: 0)', type=float, default=0)
    parser.add_argument('-c', dest='cc_thres', type=float, default=0, help='threshold for correlation for creating the unwrapping mask (0.0 --> 1.0) (defaults: 0)')

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


def select_pairs_sequential(dates, con_num):
    """Select pairs using sequential method

    Args:
        dates (list): slc dates
        con_num (int): connection number

    Returns:
        list: pairs
    """
    pairs = []
    length = len(dates)
    for i in range(length):
        if i < length - con_num:
            for j in range(con_num):
                pairs.append(f"{dates[i]}_{dates[i + j + 1]}")
        else:
            for k in range(length - i - 1):
                pairs.append(f"{dates[i]}_{dates[i + k + 1]}")

    return pairs


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


def read_gamma(file, lines, file_type):
    """Read GAMMA format file

    Args:
        file (str): GAMMA format file
        lines (int): line of file
        file_type (str): file type

    Returns:
        array: data
    """
    # check file
    if not os.path.isfile(file):
        sys.exit('{} does not exist.'.format(file))
    data = np.fromfile(file, dtype=file_type)
    # GAMMA output files are big-endian
    data.byteswap('True')
    data = data.reshape(lines, -1)

    return data


def interp_gacos(gacos_file, dem_seg_par, wavelength, incidence, out_file):
    """Read gacos file, interpolate it and save it

    Args:
        gacos_file (str): gacos file
        dem_seg_par (str): dem_seg.par file
        wavelength (float): wave length
        incidence (float): incidence angle
        out_file (str): output file
    """
    values = read_dem_par(dem_seg_par)
    width, length = values[0], values[1]
    upper_left_lon, upper_left_lat = values[3], values[2]
    lat_step, lon_step = values[4], values[5]

    lon = np.linspace(upper_left_lon, upper_left_lon + lon_step * width, width)
    lat = np.linspace(upper_left_lat, upper_left_lat + lat_step * length, length)

    lons, lats = np.meshgrid(lon, lat)

    ztd, lon_lat, size = read_gdal_file(gacos_file)
    width, length = size[0], size[1]

    lon = np.linspace(lon_lat[0], lon_lat[1], width)
    lat = np.linspace(lon_lat[3], lon_lat[2], length)

    lons_gacos, lats_gacos = np.meshgrid(lon, lat)

    ori_grid = pyresample.geometry.SwathDefinition(lons=lons_gacos, lats=lats_gacos)
    dst_grid = pyresample.geometry.SwathDefinition(lons=lons, lats=lats)

    ztd_phase_interp = pyresample.kd_tree.resample_gauss( ori_grid, ztd, dst_grid, radius_of_influence=500000, sigmas=25000)

    range2phase = -4 * np.pi / wavelength * np.cos(incidence)

    los = ztd_phase_interp * range2phase

    write_gamma(los, out_file, 'float32')


def geocode(file, width_geo, lookup_table, width, length, out_file):
    """Forward geocoding transformation using a lookup table

    Args:
        file (str): file
        width_geo (int): geo width
        lookup_table (str): lookup table
        width (int): rdr width
        length (int): rdr length
        out_file (str): output file
    """
    call_str = f"geocode {lookup_table} {file} {width_geo} {out_file} {width} {length} 1 0"
    os.system(call_str)


def sub_gacos_from_int(int_file, length, m_gacos_file, s_gacos_file, int_correct_file):
    """Subtract gacos phase from interferogram file

    Args:
        int_file (str): int file
        length (int): rdr length
        m_gacos_file (str): gacos file1
        s_gacos_file (str): gacos file2
        int_correct_file (str): output int file
    """
    m_gacos = read_gamma(m_gacos_file, length, 'float32')
    s_gacos = read_gamma(s_gacos_file, length, 'float32')
    diff_gacos = m_gacos - s_gacos
    ztd_cpx = np.cos(diff_gacos) + np.sin(diff_gacos) * 1j

    diff_int = read_gamma(int_file, length, 'complex64')
    diff_correct = diff_int * np.conj(ztd_cpx)
    write_gamma(diff_correct, int_correct_file, 'complex64')


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


def make_rdc_dem(mli, mli_par, dem, dem_par, out_dir):
    """make radar coordinate dem

    Args:
        mli (str): multi-looked slc
        mli_par (str): multi-looked slc par
        dem (str): dem
        dem_par (str): dem par
        out_dir (str): output directory

    Returns:
        str: rdc dem file
    """
    os.chdir(out_dir)

    date = os.path.basename(mli)[0:8]

    call_str = f"gc_map {mli_par} - {dem_par} {dem} dem_seg.par dem_seg lookup_table 1 1 sim_sar u v inc psi pix ls_map 8 1"
    os.system(call_str)

    call_str = f"pixel_area {mli_par} dem_seg.par dem_seg lookup_table ls_map inc pix_sigma0 pix_gamma0"
    os.system(call_str)

    width_mli = read_gamma_par(mli_par, 'range_samples')

    call_str = f"raspwr pix_gamma0 {width_mli} - - - - - - - pix_gamma0.bmp"
    os.system(call_str)

    call_str = f"create_diff_par {mli_par} - {date}.diff_par 1 0"
    os.system(call_str)

    call_str = f"offset_pwrm pix_sigma0 {mli} {date}.diff_par offs snr 64 64 offsets 2 100 100 5.0"
    os.system(call_str)

    call_str = f"offset_fitm offs snr {date}.diff_par coffs coffsets 5.0 1"
    os.system(call_str)

    width_utm_dem = read_gamma_par('dem_seg.par', 'width')

    call_str = f"gc_map_fine lookup_table {width_utm_dem} {date}.diff_par lookup_table_fine 1"
    os.system(call_str)

    length_mli = read_gamma_par(mli_par, 'azimuth_lines')
    width_mli = read_gamma_par(mli_par, 'range_samples')

    call_str = f"geocode lookup_table_fine dem_seg {width_utm_dem} {date}.dem {width_mli} {length_mli} 2 0"
    os.system(call_str)

    call_str = f"rashgt {date}.dem {mli} {width_mli} - - - - - 50 - - - {date}.dem.bmp"
    os.system(call_str)

    mli_name = os.path.basename(mli)

    call_str = f"geocode_back {mli} {width_mli} lookup_table_fine {mli_name}.geo {width_utm_dem} - 2 0"
    os.system(call_str)

    call_str = f"raspwr {mli_name}.geo {width_utm_dem} 1 0 1 1 1. .35 1 {mli_name}.geo.bmp"
    os.system(call_str)

    rdc_dem = os.path.join(out_dir, f"{date}.dem")

    return rdc_dem


def del_file(file):
    """Delete file

    Args:
        file (str): file
    """
    if os.path.isfile(file):
        os.remove(file)


def calc_time_delta(pair):
    """Calculate time delta

    Args:
        pair (str): date1-date2

    Returns:
        int: days
    """
    dates = re.findall(r'\d{8}', pair)
    date1 = datetime.datetime.strptime(dates[0], '%Y%m%d')
    date2 = datetime.datetime.strptime(dates[1], '%Y%m%d')
    time_delta = (date2 - date1).days

    return time_delta


def comb_pic(pic_path, out_pic):
    """Combine all pictures to one

    Args:
        pic_path (str): picture path
        out_pic (str): output picture
    """
    cmd_str = f"montage -label %f -geometry +5+7 -tile +6 -resize 300x300 {pic_path} {out_pic}"
    os.system(cmd_str)


def main():
    # get inputs
    inps = cmdline_parser()
    rslc_dir = os.path.abspath(inps.rslc_dir)
    out_dir = os.path.abspath(inps.out_dir)
    dem_dir = os.path.abspath(inps.dem_dir)
    ref_slc = inps.ref_slc
    rlks = inps.rlks
    alks = inps.alks
    max_sb = inps.max_sb
    max_tb = inps.max_tb
    method = inps.method
    con_num = inps.con_num
    slc_extension = inps.slc_extension
    gacos_dir = inps.gacos_dir
    wavelength = inps.wavelength
    roff = inps.roff
    loff = inps.loff
    cc_thres = inps.cc_thres

    # check rslc_dir
    if not os.path.isdir(rslc_dir):
        sys.exit('{} does not exist.'.format(rslc_dir))
    dates = sorted(
        [i for i in os.listdir(rslc_dir) if re.match(r'^\d{8}$', i)])
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
    if re.findall(r'^\d{8}$', ref_slc):
        if not ref_slc in dates:
            sys.exit('No slc for {}.'.format(ref_slc))
    else:
        sys.exit('Error date for ref_slc.')

    # check extension
    if not slc_extension.startswith('.'):
        slc_extension = '.' + slc_extension

    m_rslc_par = os.path.join(rslc_dir, ref_slc, ref_slc + slc_extension + '.par')

    # check gacos
    if gacos_dir and wavelength:
        gacos_dir = os.path.abspath(gacos_dir)
        not_exist = []
        if not os.path.isdir(gacos_dir):
            sys.exit('{} does not exist'.format(gacos_dir))
        for date in dates:
            file = os.path.join(gacos_dir, date + '.ztd.tif')
            if not os.path.isfile(file):
                not_exist.append(file)
        if not_exist:
            for i in not_exist:
                print('{} does not exist'.format(i))
            sys.exit()

    # multi-look
    mli_dir = os.path.join(out_dir, 'mli')
    if not os.path.isdir(mli_dir):
        os.mkdir(mli_dir)

    slc_tab = os.path.join(mli_dir, 'slc_tab')
    mk_tab(rslc_dir, slc_tab, slc_extension)
    mli_all(slc_tab, mli_dir, rlks, alks)

    sm_mli = os.path.join(mli_dir, ref_slc + '.rmli')
    sm_mli_par = sm_mli + '.par'

    # select pairs
    if method == 'sbas':
        if max_sb and max_tb:
            base_dir = os.path.join(out_dir, 'base_calc')
            if not os.path.isdir(base_dir):
                os.mkdir(base_dir)

            pairs = select_pairs_sbas(slc_tab, m_rslc_par, max_sb, max_tb, base_dir)
        else:
            sys.exit('Spatio-temporal baselines are required for sbas method')
    if method == 'sequential':
        if con_num:
            pairs = select_pairs_sequential(dates, con_num)
        else:
            sys.exit('Connection number is required sequential method')

    geo_dir = os.path.join(out_dir, 'geo')
    if not os.path.isdir(geo_dir):
        os.mkdir(geo_dir)

    rdc_dem = make_rdc_dem(sm_mli, sm_mli_par, dem, dem_par, geo_dir)

    width_mli = read_gamma_par(sm_mli_par, 'range_samples')
    width_mli = int(width_mli)

    length_mli = read_gamma_par(sm_mli_par, 'azimuth_lines')
    length_mli = int(length_mli)

    dem_seg_par = os.path.join(geo_dir, 'dem_seg.par')
    lookup_fine = os.path.join(geo_dir, 'lookup_table_fine')
    width_geo = read_gamma_par(dem_seg_par, 'width')
    width_geo = int(width_geo)

    # process gacos
    if gacos_dir and wavelength:
        incidence = read_gamma_par(sm_mli_par, 'incidence_angle')
        incidence = np.deg2rad(float(incidence))

        os.chdir(gacos_dir)
        for date in dates:
            gacos_file = date + '.ztd.tif'
            out_file = gacos_file + '.phase'
            interp_gacos(gacos_file, dem_seg_par, wavelength, incidence, out_file)
            out_file2 = out_file + '.rdc'
            geocode(out_file, width_geo, lookup_fine, width_mli, length_mli, out_file2)

    diff_dir = os.path.join(out_dir, 'diff')
    if not os.path.isdir(diff_dir):
        os.mkdir(diff_dir)
    os.chdir(diff_dir)

    # diff and unwrap
    for pair in pairs:
        m_date = pair[0:8]
        s_date = pair[9:17]

        m_mli = os.path.join(mli_dir, m_date + '.rmli')
        s_mli = os.path.join(mli_dir, s_date + '.rmli')

        m_rslc = os.path.join(rslc_dir, m_date, m_date + slc_extension)
        m_rslc_par = m_rslc + '.par'

        s_rslc = os.path.join(rslc_dir, s_date, s_date + slc_extension)
        s_rslc_par = s_rslc + '.par'

        call_str = f'echo "{pair}\\n\\n\\n\\n\\n\\n\\n" > off_par.in'
        os.system(call_str)
        call_str = f"create_offset {m_rslc_par} {s_rslc_par} {pair}.off 1 {rlks} {alks} < off_par.in"
        os.system(call_str)

        call_str = f"base_orbit {m_rslc_par} {s_rslc_par} {pair}.base"
        os.system(call_str)

        call_str = f'echo " {pair}\\n\\n\\n\\n\\n" > diff_par.in'
        os.system(call_str)
        call_str = f"create_diff_par {pair}.off - {pair}.diff_par 0 < diff_par.in"
        os.system(call_str)

        call_str = f"phase_sim {m_rslc_par} {pair}.off {pair}.base {rdc_dem} {pair}.sim_unw 0 0 - {calc_time_delta(pair)} 1"
        os.system(call_str)

        call_str = f"SLC_diff_intf {m_rslc} {s_rslc} {m_rslc_par} {s_rslc_par} {pair}.off {pair}.sim_unw {pair}.diff {rlks} {alks} 1 1"
        os.system(call_str)

        call_str = f"rasmph_pwr {pair}.diff {sm_mli} {width_mli} 1 1 0 1 1 1. .35 1 {pair}.diff.bmp"
        os.system(call_str)

        call_str = f"cc_wave {pair}.diff {m_mli} {s_mli} {pair}.cc {width_mli} 5 5 1"
        os.system(call_str)

        call_str = f"rascc {pair}.cc {sm_mli} {width_mli} 1 1 0  1 1 0.1 0.9 1 0.35"
        os.system(call_str)

        call_str = f"adf {pair}.diff {pair}.adf.diff {pair}.adf.cc {width_mli} 0.75 32 5 4 0 0 .2"
        os.system(call_str)

        call_str = f"rasmph_pwr {pair}.adf.diff {sm_mli} {width_mli} 1 1 0 1 1 0.7 0.35"
        os.system(call_str)

        call_str = f"rascc {pair}.adf.cc {sm_mli} {width_mli} 1 1 0 1 1 0.1 0.9 0.7 0.35"
        os.system(call_str)

        call_str = f"rascc_mask {pair}.adf.cc {sm_mli} {width_mli} 1 1 0 1 1 {cc_thres} 0 0.1 0.9 1.0 0.35 1 {pair}.adf.cc_mask.bmp"
        os.system(call_str)

        call_str = f"mcf {pair}.adf.diff {pair}.adf.cc {pair}.adf.cc_mask.bmp {pair}.adf.unw {width_mli} 1 0 0 - - 1 1 - {roff} {loff} 0"
        os.system(call_str)

        call_str = f"rasrmg {pair}.adf.unw {sm_mli} {width_mli} 1 1 0 1 1 .5 1. .35 0.0 1"
        os.system(call_str)

        call_str = f"quad_fit {pair}.adf.unw {pair}.diff_par 32 32 {pair}.adf.cc_mask.bmp {pair}.plot 3"
        os.system(call_str)

        call_str = f"quad_sub {pair}.adf.unw {pair}.diff_par {pair}.adf.unw.sub 0 0"
        os.system(call_str)

        call_str = f"rasrmg {pair}.adf.unw.sub {sm_mli} {width_mli} 1 1 0 1 1 .6 1. .35 .0 1 {pair}.adf.unw.sub.bmp {pair}.adf.cc 1 .2"
        os.system(call_str)

        if gacos_dir and wavelength:
            m_gacos = os.path.join(gacos_dir, m_date + '.ztd.tif.phase.rdc')
            s_gacos = os.path.join(gacos_dir, s_date + '.ztd.tif.phase.rdc')

            int_file = f"{pair}.diff"
            int_correct_file = f"{pair}.diff.gacos"

            sub_gacos_from_int(int_file, length_mli, m_gacos, s_gacos, int_correct_file)

            call_str = f"adf {pair}.diff.gacos {pair}.adf.diff.gacos {pair}.adf.cc.gacos {width_mli} 0.75 32 5 4 0 0 .2"
            os.system(call_str)

            call_str = f"rasmph_pwr {pair}.adf.diff.gacos {sm_mli} {width_mli} 1 1 0 1 1 0.7 0.35"
            os.system(call_str)

            call_str = f"rascc {pair}.adf.cc.gacos {sm_mli} {width_mli} 1 1 0 1 1 0.1 0.9 0.7 0.35"
            os.system(call_str)

            call_str = f"rascc_mask {pair}.adf.cc.gacos {sm_mli} {width_mli} 1 1 0 1 1 {cc_thres} 0 0.1 0.9 1.0 0.35 1 {pair}.adf.cc_mask.gacos.bmp"
            os.system(call_str)

            call_str = f"mcf {pair}.adf.diff.gacos {pair}.adf.cc.gacos {pair}.adf.cc_mask.gacos.bmp {pair}.adf.unw.gacos {width_mli} 1 0 0 - - 1 1 - {roff} {loff} 0"
            os.system(call_str)

            call_str = f"rasrmg {pair}.adf.unw.gacos {sm_mli} {width_mli} 1 1 0 1 1 .5 1. .35 0.0 1"
            os.system(call_str)

            call_str = f"quad_fit {pair}.adf.unw.gacos {pair}.diff_par 32 32 {pair}.adf.cc_mask.gacos.bmp {pair}.plot 3"
            os.system(call_str)

            call_str = f"quad_sub {pair}.adf.unw.gacos {pair}.diff_par {pair}.adf.unw.gacos.sub 0 0"
            os.system(call_str)

            call_str = f"rasrmg {pair}.adf.unw.gacos.sub {sm_mli} {width_mli} 1 1 0 1 1 .6 1. .35 .0 1 {pair}.adf.unw.gacos.sub.bmp {pair}.adf.cc.gacos 1 .2"
            os.system(call_str)

        del_file('diff_par.in')
        del_file('off_par.in')
        del_file(f'{pair}.plot')

    comb_pic(diff_dir + '/*.adf.cc.bmp', out_dir + '/adf.cc.jpg')
    comb_pic(diff_dir + '/*.adf.diff.bmp', out_dir + '/adf.diff.jpg')
    comb_pic(diff_dir + '/*.adf.unw.bmp', out_dir + '/adf.unw.jpg')
    comb_pic(diff_dir + '/*.adf.unw.sub.bmp', out_dir + '/adf.unw.sub.jpg')

    if gacos_dir and wavelength:
        comb_pic(diff_dir + '/*.adf.cc.gacos.bmp', out_dir + '/adf.cc.gacos.jpg')
        comb_pic(diff_dir + '/*.adf.diff.gacos.bmp', out_dir + '/adf.diff.gacos.jpg')
        comb_pic(diff_dir + '/*.adf.unw.gacos.bmp', out_dir + '/adf.unw.gacos.jpg')
        comb_pic(diff_dir + '/*.adf.unw.gacos.sub.bmp', out_dir + '/adf.unw.gacos.sub.jpg')


if __name__ == "__main__":
    main()

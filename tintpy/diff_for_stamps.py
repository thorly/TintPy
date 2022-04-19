#!/usr/bin/env python3
#######################################################################
# Generate differential interferogram for StaMPS PS/SBAS using GAMMA  #
# APS correction using Python (support .ztd and .ztd.tif file)        #
# Copyright (c) 2022, Lei Yuan                                        #
#######################################################################

import argparse
from email.mime.nonmultipart import MIMENonMultipart
import glob
import os
import re
import sys

import numpy as np
import pyresample
from osgeo import gdal


USAGE = """Example:
  # for PS
  python3 diff_for_stamps.py /ly/rslc /ly/PS /ly/dem 20221229 1 1 9999 9999 -f 0
  # for SBAS
  python3 diff_for_stamps.py /ly/rslc /ly/SBAS /ly/dem 20221229 1 1 200 60 -f 1
  # for aps correction (Sentinel-1 data)
  python3 diff_for_stamps.py /ly/rslc /ly/SBAS /ly/dem 20221229 1 1 200 60 -f 1 -g /ly/gacos -w 0.05546
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(
        description='Generate differential interferogram for StaMPS PS/SBAS using GAMMA.\n' + 
        'APS correction using Python (support .ztd and .ztd.tif file)',
        formatter_class=argparse.RawTextHelpFormatter, epilog=USAGE)

    parser.add_argument('rslc_dir', help='RSLCs directory')
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('dem_dir', help='DEM directory contains *.dem and *.dem.par')
    parser.add_argument('ref_rslc', help='date of reference RSLC for calculating baseline')
    parser.add_argument('rlks', help='range looks', type=int)
    parser.add_argument('alks', help='azimuth looks', type=int)
    parser.add_argument('max_sb', type=float, help='maximum spatial baseline')
    parser.add_argument('max_tb', type=float, help='maximum temporal baseline')
    parser.add_argument('-f', dest='flag', type=int, default=1, help='flag(0 for PS, 1 for SBAS, defaults: 1)')
    parser.add_argument('-e', dest='rslc_extension', type=str, default='.rslc', help='file extension for RSLCs (defaults: .rslc)')
    parser.add_argument('-g', dest='gacos_dir', help='directory contains GACOS files for aps correction')
    parser.add_argument('-w', dest='wavelength', type=float, help='Microwave length (Sentinel-1: 0.05546576, ALOS: 0.23830879)')

    inps = parser.parse_args()
    return inps


def mk_tab(slc_dir, slc_tab, slc_extension):
    """Generate SLC_tab for processing

    Args:
        slc_dir (str): slc directory
        slc_tab (str): tab file
        slc_extension (str): slc extension
    """
    dates = sorted([i for i in os.listdir(slc_dir) if re.findall(r'^\d{8}$', i)])
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


def mk_mli_all(slc_tab, out_mli_dir, rlks, alks):
    """Calculate MLI images for a stack of SLCs

    Args:
        slc_tab (str): slc tab file including slc slc_par
        out_dir (str): output directory
        rlks (int): range looks
        alks (int): azimuth looks
    """
    cmd_str = f"mk_mli_all {slc_tab} {out_mli_dir} {rlks} {alks}"
    os.system(cmd_str)


def make_rdc_dem(mli_dir, dates, dem, dem_par, out_dir):
    """make radar coordinate dem

    Args:
        mli_dir (str): mli directory
        dates (list): dates for makeing rdc dem
        dem (str): dem
        dem_par (str): dem par
        out_dir (str): output directory
    """
    os.chdir(out_dir)

    for date in dates:
        mli = glob.glob(os.path.join(mli_dir, date + '*mli'))[0]
        mli_par = mli + '.par'

        dem_seg_par = date + '.dem_seg.par'
        dem_seg = date + '.dem_seg'
        lookup = date + '.lookup'
        inc = date + '.inc'
        ls_map = date + '.ls_map'

        call_str = f"gc_map {mli_par} - {dem_par} {dem} {dem_seg_par} {dem_seg} {lookup} 2 2 - - - {inc} - - {ls_map} 8 1"
        os.system(call_str)

        pix_sigma0 = date + '.pix_sigma0'
        pix_gamma0 = date + 'pix_gamma0'
        call_str = f"pixel_area {mli_par} {dem_seg_par} {dem_seg} {lookup} {ls_map} {inc} {pix_sigma0} {pix_gamma0}"
        os.system(call_str)

        width_mli = read_gamma_par(mli_par, 'range_samples')

        call_str = f"raspwr {pix_gamma0} {width_mli} - - - - - - - {pix_gamma0}.bmp"
        os.system(call_str)

        call_str = f"create_diff_par {mli_par} - {date}.diff_par 1 0"
        os.system(call_str)

        offs = date + '.offs'
        snr = date + '.snr'
        offsets = date + '.offsets'
        call_str = f"offset_pwrm {pix_sigma0} {mli} {date}.diff_par {offs} {snr} 64 64 {offsets} 2 100 100 5.0"
        os.system(call_str)

        coffsets = date + '.coffsets'
        call_str = f"offset_fitm {offs} {snr} {date}.diff_par coffs {coffsets} 5.0 1"
        os.system(call_str)

        width_utm_dem = read_gamma_par(f'{dem_seg_par}', 'width')
        lookup_fine = date + '.lookup_fine'

        call_str = f"gc_map_fine {lookup} {width_utm_dem} {date}.diff_par {lookup_fine} 1"
        os.system(call_str)

        length_mli = read_gamma_par(mli_par, 'azimuth_lines')
        width_mli = read_gamma_par(mli_par, 'range_samples')
        rdc_dem = date + '_dem.rdc'

        call_str = f"geocode {lookup_fine} {dem_seg} {width_utm_dem} {rdc_dem} {width_mli} {length_mli} 1 0"
        os.system(call_str)

        call_str = f"rashgt {rdc_dem} {mli} {width_mli} - - - - - 50 - - - {rdc_dem}.bmp"
        os.system(call_str)


def base_calc(slc_tab, slc_par, max_sb, max_tb, flag, out_dir):
    """Generate baseline output file with perpendicular baselines and delta_T values

    Args:
        slc_tab (str): slc tab file including slc and slc_par
        slc_par (str): reference slc par
        max_sb (float): max spatial baseline
        max_tb (float): max time baseline
        flag (int): flag for itab of base_calc
        out_dir (str): output directory

    Returns:
        str: baseline file
    """
    os.chdir(out_dir)
    bperp_file = os.path.join(out_dir, 'bperp_file')
    call_str = f'base_calc {slc_tab} {slc_par} {bperp_file} itab {flag} 1 0 {max_sb} 0 {max_tb}'
    os.system(call_str)

    return bperp_file


def select_pairs_sbas(slc_tab, slc_par, max_sb, max_tb, flag, out_dir):
    """Select pairs using sbas method

    Args:
        slc_tab (str): slc tab file including slc and slc_par
        slc_par (str): reference slc par
        max_sb (float): max spatial baseline
        max_tb (float): max time baseline
        flag (int): flag for itab of base_calc
        out_dir (str): output directory

    Returns:
        list: pairs
    """
    bperp_file = base_calc(slc_tab, slc_par, max_sb, max_tb, flag, out_dir)

    pairs = []
    with open(bperp_file, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            if line:
                split_list = line.strip().split()
                date1 = split_list[1]
                date2 = split_list[2]
                if flag == 0:
                    pairs.append(date1 + '_' + date2)
                if flag == 1:
                    if int(date1) > int(date2):
                        pairs.append(date2 + '_' + date1)
                    else:
                        pairs.append(date1 + '_' + date2)

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
    out_dir = os.path.abspath(inps.out_dir)
    dem_dir = os.path.abspath(inps.dem_dir)
    ref_rslc = inps.ref_rslc
    rlks = inps.rlks
    alks = inps.alks
    max_sb = inps.max_sb
    max_tb = inps.max_tb

    flag = inps.flag

    rslc_extension = inps.rslc_extension

    gacos_dir = inps.gacos_dir
    wavelength = inps.wavelength


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

    # check flag
    if flag not in [0, 1]:
        sys.exit("flag must be 0 or 1 (0 for PS, 1 for SBAS).")

    # check gacos_dir and wavelength
    if gacos_dir and wavelength is None:
        sys.exit('wavelength(-w) is required for aps correction.')
    if gacos_dir is None and wavelength:
        sys.exit('gacos_dir(-g) is required for aps correction.')
    if gacos_dir and wavelength:
        gacos_dir = os.path.abspath(gacos_dir)
        if not os.path.isdir(gacos_dir):
            sys.exit('{} does not exist.'.format(gacos_dir))

        gacos_extension = get_gacos_extension(gacos_dir)
        if gacos_extension is None:
            sys.exit(f'Cannot find *.ztd or *.ztd.tif files in {gacos_dir}')
        check_gacos(dates, gacos_dir, gacos_extension)

    mli_dir = os.path.join(out_dir, 'mli')
    if not os.path.isdir(mli_dir):
        os.mkdir(mli_dir)

    rslc_tab = os.path.join(mli_dir, 'rslc_tab')
    mk_tab(rslc_dir, rslc_tab, rslc_extension)

    # select pairs
    base_dir = os.path.join(out_dir, 'base_calc')
    if not os.path.isdir(base_dir):
        os.mkdir(base_dir)

    ref_rslc_par = os.path.join(rslc_dir, ref_rslc, ref_rslc + rslc_extension + '.par')
    pairs = select_pairs_sbas(rslc_tab, ref_rslc_par, max_sb, max_tb, flag, base_dir)

    # multi-look
    mk_mli_all(rslc_tab, mli_dir, rlks, alks)

    # mk_rdc_dem
    geo_dir = os.path.join(out_dir, 'geo')
    if not os.path.isdir(geo_dir):
        os.mkdir(geo_dir)

    # for PS
    if flag == 0:
        make_rdc_dem(mli_dir, [ref_rslc], dem, dem_par, geo_dir)
        pairs = [i for i in pairs if i[0:8] != i[9:17]]
    # for SBAS
    if flag == 1:
        make_rdc_dem(mli_dir, dates[0:-1], dem, dem_par, geo_dir)

    # diff
    diff_dir = os.path.join(out_dir, 'diff')
    if not os.path.isdir(diff_dir):
        os.mkdir(diff_dir)

    os.chdir(diff_dir)

    for pair in pairs:
        m_date = pair[0:8]
        s_date = pair[9:17]

        m_mli = glob.glob(os.path.join(mli_dir, m_date + '*mli'))[0]
        s_mli = glob.glob(os.path.join(mli_dir, s_date + '*mli'))[0]
        dem_rdc = os.path.join(geo_dir, m_date + '_dem.rdc')
        diff_par = os.path.join(geo_dir, m_date + '.diff_par')

        m_rslc = os.path.join(rslc_dir, m_date, m_date + rslc_extension)
        m_rslc_par = m_rslc + '.par'
        s_rslc = os.path.join(rslc_dir, s_date, s_date + rslc_extension)
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

        call_str = (f"phase_sim_orb {m_rslc_par} {s_rslc_par} {pair}.off {dem_rdc} {pair}.sim_unw {m_rslc_par} - - 1 1")
        os.system(call_str)

        call_str = f"SLC_diff_intf {m_rslc} {s_rslc} {m_rslc_par} {s_rslc_par} {pair}.off {pair}.sim_unw {pair}.diff {rlks} {alks} 1 1"
        os.system(call_str)

        width_mli = read_gamma_par(m_mli + '.par', 'range_samples')
        call_str = f"rasmph_pwr {pair}.diff {m_mli} {width_mli} 1 1 0 1 1 1. .35 1 {pair}.diff_pwr.bmp"
        os.system(call_str)

        call_str = f"cc_wave {pair}.diff {m_mli} {s_mli} {pair}.cc {width_mli} 5 5 1"
        os.system(call_str)

        call_str = f"rascc {pair}.cc {m_mli} {width_mli} 1 1 0  1 1 0.1 0.9 1 0.35"
        os.system(call_str)

    # aps correction
    if gacos_dir and wavelength:
        gacos_dir = os.path.abspath(gacos_dir)
        gacos_extension = get_gacos_extension(gacos_dir)

        os.chdir(diff_dir)

        for pair in pairs:
            m_date = pair[0:8]
            s_date = pair[9:17]
            m_gacos = os.path.join(gacos_dir, m_date + gacos_extension)
            s_gacos = os.path.join(gacos_dir, s_date + gacos_extension)
            dem_seg_par = os.path.join(geo_dir, m_date + '.dem_seg.par')
            incidence = read_gamma_par(m_rslc_par, 'incidence_angle')
            incidence = np.deg2rad(float(incidence))

            out_file = pair + '.gacos'
            if gacos_extension == '.ztd':
                process_gacos(m_gacos, s_gacos, 'bin', dem_seg_par, wavelength, incidence, out_file)
            if gacos_extension == '.ztd.tif':
                process_gacos(m_gacos, s_gacos, 'tif', dem_seg_par, wavelength, incidence, out_file)

            lookup = os.path.join(geo_dir, m_date + '.lookup_fine')
            width_geo = read_gamma_par(dem_seg_par, "width:")
            off_par = pair + '.off'
            width_rdc = read_gamma_par(off_par, "interferogram_width:")
            line_rdc = read_gamma_par(off_par, "interferogram_azimuth_lines:")

            cmd_str = f"geocode {lookup} {out_file} {width_geo} {out_file}.rdc {width_rdc} {line_rdc} 1 0"
            os.system(cmd_str)
            
            cmd_str = f"raspwr {out_file}.rdc $width_rdc 1 0 1 1 1. .35 1 {out_file}.rdc.bmp"
            os.system(cmd_str)

            cmd_str = f"sub_phase {pair}.diff {out_file}.rdc {pair}.diff_par {pair}.diff.gacos 1 0"
            os.system(cmd_str)

            m_mli = glob.glob(os.path.join(mli_dir, m_date + '*mli'))[0]
            cmd_str = f"rasmph_pwr {pair}.diff.gacos {m_mli} {width_rdc} 1 1 0 1 1 1. 0.35 1 {pair}.diff.gacos_pwr.bmp"
            os.system(cmd_str)

    print('\nAll done, enjoy it!\n')

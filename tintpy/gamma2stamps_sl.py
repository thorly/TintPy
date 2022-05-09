#!/usr/bin/env python3
################################
# Program is part of TintPy    #
# Copyright (c) 2021, Lei Yuan #
# Author: Lei Yuan, 2022       #
################################

# Expected directory structure:
# PS:
#   rslc/*.rslc
#   rslc/*.rslc.par
#   diff0/*.diff
#   diff0/*.base
#   geo/*dem.rdc
#   geo/*diff_par
#   geo/YYYYMMDD.lon (master)
#   geo/YYYYMMDD.lat (master)
#   dem/*_seg.par
#
# SB:
#   SMALL_BASELINES/YYYYMMDD_YYYYMMDD/YYYYMMDD.rslc (master)
#   SMALL_BASELINES/YYYYMMDD_YYYYMMDD/YYYYMMDD.rslc (slave)
#   SMALL_BASELINES/YYYYMMDD_YYYYMMDD/*.rslc.par
#   SMALL_BASELINES/YYYYMMDD_YYYYMMDD/*.diff
#   SMALL_BASELINES/YYYYMMDD_YYYYMMDD/*.base
#   geo/*dem.rdc
#   geo/*diff_par
#   geo/YYYYMMDD.lon (master)
#   geo/YYYYMMDD.lat (master)
#   dem/*_seg.par

import argparse
import glob
import os
import re
import shutil
import sys

import numpy as np


EXAMPLE = """Example:
  python3 gamma2stamps_sl.py /ly/rslc /ly/diff /ly/geo 20221229 PS /ly/PS
  python3 gamma2stamps_sl.py /ly/rslc /ly/diff /ly/geo 20221229 SBAS /ly/SBAS
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(
        description='Prepare single-looked files processed by GAMMA for StaMPS PS/SBAS processing.',
        epilog=EXAMPLE, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('rslc_dir', help='directory including multilooked RSLCs')
    parser.add_argument('diff_dir', help='directory including diff files')
    parser.add_argument('geo_dir', help='directory including lookup file')
    parser.add_argument('supermaster', help='supermaster rslc date')
    parser.add_argument('method', help='PS or SBAS method', choices=['ps', 'PS', 'sbas', 'SBAS'])
    parser.add_argument('output_dir', help='parent directory of INSAR_supermaster')

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
                value = line.strip().split()[1]

    return value


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


def mk_dir(dir):
    """Make directory

    Args:
        dir (str): directory
    """
    if not os.path.isdir(dir):
        os.mkdir(dir)


def is_dir(dir):
    """Check if it is directory

    Args:
        dir (str): directory
    """
    if not os.path.isdir(dir):
        sys.exit('{} does not exist.'.format(dir))


def gen_lon_lat(dem_seg_par, lookup, diff_par, supermaster, out_dir):
    """Generate lon and lat file for StaMPS processing

    Args:
        dem_seg_par (str): dem par file
        lookup (str): lookup table
        diff_par (str): diff_par file
        supermaster (str): master date
        out_dir (srt): output directory
    """
    width = read_gamma_par(diff_par, 'range_samp_1')
    length = read_gamma_par(diff_par, 'az_samp_1')

    # Extract geocoding information to generate the lon and lat matrices
    lat_max = float(read_gamma_par(dem_seg_par, 'corner_lat'))
    lon_min = float(read_gamma_par(dem_seg_par, 'corner_lon'))
    lat_step = float(read_gamma_par(dem_seg_par, 'post_lat'))
    lon_step = float(read_gamma_par(dem_seg_par, 'post_lon'))
    length_dem = int(read_gamma_par(dem_seg_par, 'nlines'))
    width_dem = int(read_gamma_par(dem_seg_par, 'width'))

    lat_min = lat_max + lat_step * (length_dem - 1)
    lon_max = lon_min + lon_step * (width_dem - 1)

    lon = np.linspace(lon_min, lon_max, width_dem)
    lat = np.linspace(lat_max, lat_min, length_dem)

    lons, lats = np.meshgrid(lon, lat)

    lon_geo = os.path.join(out_dir, 'lon_geo')
    lat_geo = os.path.join(out_dir, 'lat_geo')

    write_gamma(lons, lon_geo, 'float32')
    write_gamma(lats, lat_geo, 'float32')

    lon_file = os.path.join(out_dir, supermaster + '.lon')
    lat_file = os.path.join(out_dir, supermaster + '.lat')

    cmd_str = f'geocode {lookup} {lon_geo} {width_dem} {lon_file} {width} {length} 2 0'
    os.system(cmd_str)

    cmd_str = f'geocode {lookup} {lat_geo} {width_dem} {lat_file} {width} {length} 2 0'
    os.system(cmd_str)

    if os.path.isfile(lon_geo):
        os.remove(lon_geo)

    if os.path.isfile(lat_geo):
        os.remove(lat_geo)


def prep_files_for_ps(in_rslc_dir, in_diff_dir, in_geo_dir, supermaster, output_dir):
    """Prepare files for StaMPS PS processing

    Args:
        in_rslc_dir (str): rslc directory
        in_diff_dir (str): diff directory
        in_geo_dir (str): geo directory
        supermaster (str): supermaster rslc date
        output_dir (str): output directory
    """
    # create directories
    insar_dir = os.path.join(output_dir, 'INSAR_' + supermaster)
    mk_dir(insar_dir)

    rslc_dir = os.path.join(insar_dir, 'rslc')
    mk_dir(rslc_dir)

    geo_dir = os.path.join(insar_dir, 'geo')
    mk_dir(geo_dir)

    diff_dir = os.path.join(insar_dir, 'diff0')
    mk_dir(diff_dir)

    dem_dir = os.path.join(insar_dir, 'dem')
    mk_dir(dem_dir)

    print('prepare geo directory')

    # prep geo/*dem.rdc
    dem_rdc = os.path.join(in_geo_dir, supermaster + '_dem.rdc')
    dem_rdc_dst = os.path.join(geo_dir, supermaster + '.dem.rdc')
    shutil.copy(dem_rdc, dem_rdc_dst)

    # prep geo/*diff_par
    diff_par = os.path.join(in_geo_dir, supermaster + '.diff_par')
    diff_par_dst = os.path.join(geo_dir, supermaster + '.diff_par')
    shutil.copy(diff_par, diff_par_dst)

    # prep geo/YYYYMMDD.lon (master) geo/YYYYMMDD.lat (master)
    dem_seg_par = os.path.join(in_geo_dir, supermaster + '.dem_seg.par')
    lookup = os.path.join(in_geo_dir, supermaster + '.lookup_fine')
    gen_lon_lat(dem_seg_par, lookup, diff_par, supermaster, geo_dir)

    print('done')

    print('prepare dem directory')

    # prep dem/*_seg.par
    dem_seg_par_dst = os.path.join(dem_dir, supermaster + '_seg.par')
    shutil.copy(dem_seg_par, dem_seg_par_dst)

    print('done')

    print('prepare diff0 directory')
    # get ifg_pairs
    ifg_pairs = [i[0:17] for i in os.listdir(in_diff_dir) if re.match(r'\d{8}_\d{8}', i)]
    ifg_pairs = list(set(ifg_pairs))

    dates = []
    for ifg in ifg_pairs:
        ifg_out = ifg.replace('_', '-')

        # prep .diff
        diff = os.path.join(in_diff_dir, ifg + '.diff')
        diff_dst = os.path.join(diff_dir, ifg_out + '.diff')
        shutil.copy(diff, diff_dst)

        # prep .base
        baseline = os.path.join(in_diff_dir, ifg + '.base')
        baseline_dst = os.path.join(diff_dir, ifg_out + '.base')
        shutil.copy(baseline, baseline_dst)

        dates.append(ifg[0:8])
        dates.append(ifg[9:17])

    print('done')

    print('prepare rslc directory')

    # prep .rslc and .rslc.par
    dates = sorted(list(set(dates)))
    for date in dates:
        rslc = glob.glob(os.path.join(in_rslc_dir, date, date + '*slc'))[0]
        rslc_dst = os.path.join(rslc_dir, date + '.rslc')
        shutil.copy(rslc, rslc_dst)

        rslc_par = rslc + '.par'
        rslc_par_dst = rslc_dst + '.par'
        shutil.copy(rslc_par, rslc_par_dst)
    
    print('done')

    print('\nAll done, enjoy it, you can run mt_prep_gamma in terminal!\n')


def prep_files_for_sbas(rslc_dir, diff_dir, in_geo_dir, supermaster, output_dir):
    """Prepare files for StaMPS SBAS processing

    Args:
        rslc_dir (str): rslc directory
        diff_dir (str): diff directory
        in_geo_dir (str): geo directory
        supermaster (str): supermaster rslc date
        output_dir (str): output directory
    """
    # create directories
    insar_dir = os.path.join(output_dir, 'INSAR_' + supermaster)
    mk_dir(insar_dir)

    sb_dir = os.path.join(insar_dir, 'SMALL_BASELINES')
    mk_dir(sb_dir)

    geo_dir = os.path.join(insar_dir, 'geo')
    mk_dir(geo_dir)

    dem_dir = os.path.join(insar_dir, 'dem')
    mk_dir(dem_dir)

    print('prepare geo directory')

    # prep geo/*dem.rdc
    dem_rdc = os.path.join(in_geo_dir, supermaster + '_dem.rdc')
    dem_rdc_dst = os.path.join(geo_dir, supermaster + '.dem.rdc')
    shutil.copy(dem_rdc, dem_rdc_dst)

    # prep geo/*diff_par
    diff_par = os.path.join(in_geo_dir, supermaster + '.diff_par')
    diff_par_dst = os.path.join(geo_dir, supermaster + '.diff_par')
    shutil.copy(diff_par, diff_par_dst)

    # prep geo/YYYYMMDD.lon (master) geo/YYYYMMDD.lat (master)
    dem_seg_par = os.path.join(in_geo_dir, supermaster + '.dem_seg.par')
    lookup = os.path.join(in_geo_dir, supermaster + '.lookup_fine')
    gen_lon_lat(dem_seg_par, lookup, diff_par, supermaster, geo_dir)

    print('done')

    print('prepare dem directory')

    # dem/*_seg.par
    dem_seg_par_dst = os.path.join(dem_dir, supermaster + '_seg.par')
    shutil.copy(dem_seg_par, dem_seg_par_dst)

    print('done')

    print('prepare SMALL_BASELINES directory')
    # get ifg_pairs
    ifg_pairs = [i[0:17] for i in os.listdir(diff_dir) if re.match(r'\d{8}_\d{8}', i)]
    ifg_pairs = list(set(ifg_pairs))

    for ifg in ifg_pairs:
        ifg_out = ifg.replace('_', '-')

        # mkdir
        ifg_out_dir = os.path.join(sb_dir, ifg)
        mk_dir(ifg_out_dir)

        # prep .diff
        diff = os.path.join(diff_dir, ifg + '.diff')
        diff_dst = os.path.join(ifg_out_dir, ifg_out + '.diff')
        shutil.copy(diff, diff_dst)

        # prep .base
        baseline = os.path.join(diff_dir, ifg + '.base')
        baseline_dst = os.path.join(ifg_out_dir, ifg_out + '.base')
        shutil.copy(baseline, baseline_dst)

        # prep YYYYMMDD.rslc and .rslc.par
        m_rslc = glob.glob(os.path.join(rslc_dir, ifg[0:8], ifg[0:8] + '*slc'))[0]
        m_rslc_dst = os.path.join(ifg_out_dir, ifg[0:8] + '.rslc')
        shutil.copy(m_rslc, m_rslc_dst)

        m_rslc_par = m_rslc + '.par'
        m_rslc_par_dst = m_rslc_dst + '.par'
        shutil.copy(m_rslc_par, m_rslc_par_dst)

        s_rslc = glob.glob(os.path.join(rslc_dir, ifg[9:17], ifg[9:17] + '*slc'))[0]
        s_rslc_dst = os.path.join(ifg_out_dir, ifg[9:17] + '.rslc')
        shutil.copy(s_rslc, s_rslc_dst)

    print('done')

    print('\nAll done, enjoy it, you can run mt_prep_gamma in terminal!\n')



def main():
    # get inputs
    inps = cmdline_parser()
    rslc_dir = os.path.abspath(inps.rslc_dir)
    diff_dir = os.path.abspath(inps.diff_dir)
    geo_dir = os.path.abspath(inps.geo_dir)
    supermaster = inps.supermaster
    method = inps.method.lower()
    output_dir = os.path.abspath(inps.output_dir)

    mk_dir(output_dir)

    # run
    if method == 'ps':
        func = prep_files_for_ps
    if method == 'sbas':
        func = prep_files_for_sbas

    func(rslc_dir, diff_dir, geo_dir, supermaster, output_dir)


if __name__ == "__main__":
    main()

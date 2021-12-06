#!/usr/bin/env python3
################################
# Program is part of TintPy    #
# Copyright (c) 2021, Lei Yuan #
# Author: Lei Yuan, 2021       #
################################

# Expected directory structure:
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

EXAMPLE = """Example:
  python3 sbas_mli_ad.py /ly/mli /ly/diff /ly/geo 20211229 /ly/ad
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(
        description='Prepare necessary files(multilooked) processed by GAMMA for StaMPS SBAS processing.',
        epilog=EXAMPLE,formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('mli_dir', help='directory including multilooked RSLCs')
    parser.add_argument('diff_dir', help='directory including diff files')
    parser.add_argument('geo_dir', help='directory including lookup file')
    parser.add_argument('supermaster', help='supermaster for coregistration')
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
                value = line.split()[1].strip()

    return value


def gen_lon_lat(dem_par, lookup, diff_par, supermaster, out_dir):
    """Generate lon and lat file for StaMPS processing

    Args:
        dem_par (str): dem par file
        lookup (str): lookup table
        diff_par (str): diff_par file
        supermaster (str): master date
        out_dir (srt): output directory
    """
    # width of interferogram
    width = read_gamma_par(diff_par, 'range_samp_1')
    # length of inerferogram
    length = read_gamma_par(diff_par, 'az_samp_1')

    # Extract geocoding information to generate the lon and lat matrices
    lat = float(read_gamma_par(dem_par, 'corner_lat').split()[0])
    lon = float(read_gamma_par(dem_par, 'corner_lon').split()[0])
    lat_step = float(read_gamma_par(dem_par, 'post_lat').split()[0])
    lon_step = float(read_gamma_par(dem_par, 'post_lon').split()[0])
    length_dem = int(read_gamma_par(dem_par, 'nlines'))
    width_dem = int(read_gamma_par(dem_par, 'width'))
    lat1 = lat + lat_step * (length_dem - 1)
    lat1 = round(float(lat1), 6)
    lon1 = lon + lon_step * (width_dem - 1)
    lon1 = round(float(lon1), 6)

    # change path of running
    os.chdir(out_dir)

    # longitude file
    cmd = f"gmt grdmath -R{lon}/{lon1}/{lat1}/{lat} -I{width_dem}+/{length_dem}+ X = geo.grd"
    os.system(cmd)

    # take lons
    cmd = "gmt grd2xyz geo.grd -ZTLf > geo.raw"
    os.system(cmd)

    # set lons to 4-byte floats
    cmd = "swap_bytes geo.raw geolon.raw 4"
    os.system(cmd)

    # geocode
    cmd = f"geocode {lookup} geolon.raw {width_dem} {supermaster}.lon {width} {length} 2 0"
    os.system(cmd)

    # latitude file
    cmd = f"gmt grdmath -R{lon}/{lon1}/{lat1}/{lat} -I{width_dem}+/{length_dem}+ Y = geo.grd"
    os.system(cmd)

    # take lats
    cmd = "gmt grd2xyz geo.grd -ZTLf > geo.raw"
    os.system(cmd)

    # set lats to 4-byte floats
    cmd = "swap_bytes geo.raw geolat.raw 4"
    os.system(cmd)

    # geocode
    cmd = f"geocode {lookup} geolat.raw {width_dem} {supermaster}.lat {width} {length} 2 0"
    os.system(cmd)

    # cleaning
    cmd = "rm -rf geo.raw geolon.raw geolat.raw geo.grd gmt.history"
    os.system(cmd)


def prep_sb_dir(mli_dir, diff_dir, supermaster, sb_dir):
    """Prepare smallbaseline directory

    Args:
        mli_dir (str): mli directory
        diff_dir (str): diff directory
        supermaster (str): mater date
        sb_dir (str): output directory
    """
    # get content of mli.par and modify it
    mli_par = os.path.join(mli_dir, supermaster + '.rmli.par')
    with open(mli_par, 'r') as f:
        mli_par_content = f.read()
        mli_par_content = mli_par_content.replace('FLOAT', 'FCOMPLEX')

    width = read_gamma_par(mli_par, 'range_samples')

    # get ifg_pairs
    ifg_pairs = [i[0:17] for i in os.listdir(diff_dir) if re.match(r'\d{8}_\d{8}', i)]
    ifg_pairs = list(set(ifg_pairs))

    for ifg in ifg_pairs:
        ifg_out = ifg.replace('_', '-')

        # mkdir
        ifg_out_dir = os.path.join(sb_dir, ifg)
        if not os.path.isdir(ifg_out_dir):
            os.mkdir(ifg_out_dir)

        # prep .diff
        diff = os.path.join(diff_dir, ifg + '.diff')
        diff_dst = os.path.join(ifg_out_dir, ifg_out + '.diff')
        shutil.copy(diff, diff_dst)

        # prep .base
        baseline = os.path.join(diff_dir, ifg + '.base')
        baseline_dst = os.path.join(ifg_out_dir, ifg_out + '.base')
        shutil.copy(baseline, baseline_dst)

        # prep YYYYMMDD.rslc and .rslc.par
        # real_to_cpx
        mli1 = os.path.join(mli_dir, ifg[0:8] + '.rmli')
        mli1_dst = os.path.join(ifg_out_dir, ifg[0:8] + '.rslc')

        call_str = f"real_to_cpx {mli1} - {mli1_dst} {width} 0"
        os.system(call_str)

        mli2 = os.path.join(mli_dir, ifg[9:17] + '.rmli')
        mli2_dst = os.path.join(ifg_out_dir, ifg[9:17] + '.rslc')

        call_str = f"real_to_cpx {mli2} - {mli2_dst} {width} 0"
        os.system(call_str)

        # write .rslc.par
        mli_par_path = os.path.join(ifg_out_dir, supermaster + '.rslc.par')
        with open(mli_par_path, 'w+') as f:
            f.write(mli_par_content)


def prep_geo_dir(in_geo_dir, diff_dir, supermaster, out_geo_dir):
    """Prepare geo directory

    Args:
        in_geo_dir (str): geo directory
        diff_dir (str): diff directory
        supermaster (str): master date
        out_geo_dir (str): output directory
    """
    # prep geo/*dem.rdc
    dem_rdc = os.path.join(in_geo_dir, supermaster + '.dem')
    dem_rdc_dst = os.path.join(out_geo_dir, supermaster + '.dem.rdc')
    shutil.copy(dem_rdc, dem_rdc_dst)

    # prep geo/*diff_par
    diff_par = glob.glob(os.path.join(in_geo_dir, '*.diff_par'))[0]
    diff_par_dst = os.path.join(out_geo_dir, supermaster + '.diff_par')
    shutil.copy(diff_par, diff_par_dst)

    # prep geo/YYYYMMDD.lon (master) geo/YYYYMMDD.lat (master)
    dem_par = os.path.join(in_geo_dir, 'dem_seg.par')
    lookup = os.path.join(in_geo_dir, 'lookup_table_fine')
    gen_lon_lat(dem_par, lookup, diff_par, supermaster, out_geo_dir)


def prep_dem_dir(in_geo_dir, supermaster, dem_dir):
    """Prepare dem directory

    Args:
        in_geo_dir (str): geo directory
        supermaster (str): master date
        dem_dir (str): output directory
    """
    # dem/*_seg.par
    seg_par = os.path.join(in_geo_dir, 'dem_seg.par')
    seg_par_dst = os.path.join(dem_dir, supermaster + '_seg.par')
    shutil.copy(seg_par, seg_par_dst)


def is_dir(dir):
    """Check if it is directory

    Args:
        dir (str): directory
    """
    if not os.path.isdir(dir):
        sys.exit('{} does not exist.'.format(dir))


def mk_dir(dir):
    """Make directory

    Args:
        dir (str): directory
    """
    if not os.path.isdir(dir):
        os.mkdir(dir)


def main():
    # get inputs
    inps = cmdline_parser()
    mli_dir = os.path.abspath(inps.mli_dir)
    diff_dir = os.path.abspath(inps.diff_dir)
    in_geo_dir = os.path.abspath(inps.geo_dir)
    supermaster = inps.supermaster
    output_dir = os.path.abspath(inps.output_dir)

    # check dirs:
    is_dir(mli_dir)
    is_dir(diff_dir)
    is_dir(in_geo_dir)
    mk_dir(output_dir)

    # mkdir
    insar_dir = os.path.join(output_dir, 'INSAR_' + supermaster)
    mk_dir(insar_dir)

    # SMALL_BASELINES
    sb_dir = os.path.join(insar_dir, 'SMALL_BASELINES')
    mk_dir(sb_dir)

    prep_sb_dir(mli_dir, diff_dir, supermaster, sb_dir)

    # geo
    out_geo_dir = os.path.join(insar_dir, 'geo')
    mk_dir(out_geo_dir)

    prep_geo_dir(in_geo_dir, diff_dir, supermaster, out_geo_dir)

    # dem
    dem_dir = os.path.join(insar_dir, 'dem')
    mk_dir(dem_dir)

    prep_dem_dir(in_geo_dir, supermaster, dem_dir)

    print('\nAll done, enjoy it, you can run mt_prep_gamma in terminal!\n')


if __name__ == "__main__":
    main()

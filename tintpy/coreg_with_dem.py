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

EXAMPLE = """Example:
  python3 coreg_with_dem.py /ly/slc /ly/dem /ly/mli /ly/rslc 
  python3 coreg_with_dem.py /ly/slc /ly/dem /ly/mli /ly/rslc --rlks 2 --alks 10 --ref_slc 20201111
"""


def cmd_line_parser():
    parser = argparse.ArgumentParser(
        description='Coregister all of SLCs to a reference SLC with DEM.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLE)

    parser.add_argument('slc_dir', help='directory path of SLCs')
    parser.add_argument('dem_dir', help='directory of *.dem and *.dem.par')
    parser.add_argument('mli_dir',
                        help='output directory of multi-looked SLCs')
    parser.add_argument('rslc_dir', help='output directory of RSLCs')
    parser.add_argument('--rlks',
                        help='range looks (defaults: 8)',
                        default=8,
                        type=int)
    parser.add_argument('--alks',
                        help='azimuth looks (defaults: 2)',
                        default=2,
                        type=int)
    parser.add_argument('--ref_slc',
                        help='reference SLC (default: the first slc)',
                        default='0')
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
                value = line.split(':')[1].strip()

    return value


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


def mk_tab(slc_dir, slc_tab):
    """Generate SLC_tab for processing

    Args:
        slc_dir (str): slc directory
        slc_tab (str): tab file
    """
    dates = sorted(
        [i for i in os.listdir(slc_dir) if re.findall(r'^\d{8}$', i)])
    with open(slc_tab, 'w+') as f:
        for date in dates:
            slc = os.path.join(slc_dir, date, date + '.slc')
            slc_par = slc + '.par'
            f.write(slc + '    ' + slc_par + '\n')


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

    call_str = f"geocode_back {mli} {width_mli} lookup_table_fine {mli}.geo {width_utm_dem} - 2 0"
    os.system(call_str)

    call_str = f"raspwr {mli}.geo {width_utm_dem} 1 0 1 1 1. .35 1 {mli}.geo.bmp"
    os.system(call_str)

    rdc_dem = os.path.join(out_dir, f"{date}.dem")

    return rdc_dem


def print_coreg_quality(quality_files):
    """print coregistration quality

    Args:
        quality_files (str list): quality files
    """
    print('\n\n')
    print('-' * 36)
    print('|{:^34}|'.format('coregistration quality report'))
    print('-' * 36)
    print('|{:^12}|{:^10}|{:^10}|'.format('date', 'range', 'azimuth'))
    print('-' * 36)
    for file in quality_files:
        date = file.split('_')[1]
        with open(file, 'r', encoding='utf-8') as f:
            for line in f.readlines()[::-1]:
                if line.strip().startswith('final model fit std. dev'):
                    r_std = line.strip().split()[-3]
                    a_std = line.strip().split()[-1]
                    print('|{:^12}|{:^10}|{:^10}|'.format(date, r_std, a_std))
                    print('-' * 36)
                    break


def main():
    # get inputs
    inps = cmd_line_parser()
    slc_dir = os.path.abspath(inps.slc_dir)
    mli_dir = os.path.abspath(inps.mli_dir)
    rslc_dir = os.path.abspath(inps.rslc_dir)
    dem_dir = os.path.abspath(inps.dem_dir)
    rlks = inps.rlks
    alks = inps.alks
    ref_slc = inps.ref_slc

    # check slc_dir
    if not os.path.isdir(slc_dir):
        sys.exit("{} does not exist.".format(slc_dir))

    # check mli_dir
    if not os.path.isdir(mli_dir):
        os.mkdir(mli_dir)

    # check rslc_dir
    if not os.path.isdir(rslc_dir):
        os.mkdir(rslc_dir)

    # check dem
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
                sys.exit(f'Cannot find *.dem and *.dem.par in {dem_dir}.')

    # get all date
    dates = sorted(
        [i for i in os.listdir(slc_dir) if re.findall(r'^\d{8}$', i)])
    if len(dates) < 2:
        sys.exit('No enough SLCs.')

    # check ref_slc
    if ref_slc == '0':
        ref_slc = dates[0]
    else:
        if re.findall(r'^\d{8}$', ref_slc):
            if not ref_slc in dates:
                sys.exit('No slc for {}.'.format(ref_slc))
        else:
            sys.exit('Error date for ref_slc.')

    m_date = ref_slc
    # write slc_tab
    slc_tab = os.path.join(mli_dir, 'slc_tab')
    mk_tab(slc_dir, slc_tab)

    # multi-look slc
    call_str = f"mk_mli_all {slc_tab} {mli_dir} {rlks} {alks}"
    os.system(call_str)

    m_mli = os.path.join(mli_dir, f"{m_date}.mli")
    m_mli_par = m_mli + '.par'

    rdc_dem = make_rdc_dem(m_mli, m_mli_par, dem, dem_par, mli_dir)

    # Resample a set of SLCs to a common reference SLC using a lookup table and DEM
    m_slc = os.path.join(slc_dir, m_date, m_date + '.slc')
    m_slc_par = m_slc + '.par'
    rslc_tab = os.path.join(rslc_dir, 'rslc_tab')

    os.chdir(rslc_dir)
    for flag in range(5):
        call_str = f"SLC_resamp_lt_all {slc_tab} {m_slc} {m_slc_par} {m_mli_par} {rdc_dem} {mli_dir} {rslc_dir} {rslc_tab} {flag}"
        os.system(call_str)

    # move date.rslc date.rslc.par to date directory
    for date in dates:
        date_dir = os.path.join(rslc_dir, date)
        if not os.path.isdir(date_dir):
            os.mkdir(date_dir)

        rslc = date + '.rslc'
        rslc_par = rslc + '.par'
        if os.path.isfile(rslc) and os.path.isfile(rslc_par):
            shutil.move(rslc, date_dir)
            shutil.move(rslc_par, date_dir)

    # generate bmp for rslc
    rslc_files = glob.glob(rslc_dir + '/*/*.rslc')
    for rslc in rslc_files:
        bmp = rslc + '.bmp'
        slc2bmp(rslc, rslc + '.par', rlks, alks, bmp)

    # check quality
    log_files = sorted(glob.glob('*resamp_lt.log'))
    print_coreg_quality(log_files)

    # clean rslc_dir
    for file in os.listdir(rslc_dir):
        if os.path.isfile(file):
            if file not in log_files:
                os.remove(file)

    print('\nAll done, enjoy it.\n')


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
################################
# Program is part of TintPy    #
# Copyright (c) 2021, Lei Yuan #
# Author: Lei Yuan, 2021       #
################################

import argparse
import datetime
import glob
import os
import re
import shutil
import sys
import zipfile

EXAMPLE = """Example:
  python3 s1_zip2slc.py /ly/zip_dir /ly/orbits /ly/slc 1
  python3 s1_zip2slc.py /ly/zip_dir /ly/orbits /ly/slc 1 2 3 --pol vv --rlks 20 --alks 5
"""


def cmdLineParse():
    parser = argparse.ArgumentParser(
        description=
        'Generate SLC from Sentinel-1 raw data with orbit correction using GAMMA.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLE)

    parser.add_argument('s1_zip_dir', help='Sentinel-1 raw directory')
    parser.add_argument('orbit_dir', help='precise orbit directory')
    parser.add_argument('slc_dir', help='directory for saving slc')
    parser.add_argument('sub_swath',
                        help='sub_swath number',
                        type=int,
                        choices=[1, 2, 3],
                        nargs='+')
    parser.add_argument('--pol',
                        help='polarization(defaults: vv)',
                        nargs='+',
                        choices=['vv', 'vh'],
                        default=['vv'])
    parser.add_argument('--rlks',
                        help='range looks(default: 8)',
                        type=int,
                        default=8)
    parser.add_argument('--alks',
                        help='azimuth looks(default: 2)',
                        type=int,
                        default=2)
    inps = parser.parse_args()

    return inps


def get_date_from_zip(zip_file):
    """Get date from S1 zip file

    Args:
        zip_file (str): S1 zip file

    Returns:
        str: S1 date
    """
    file_name = os.path.basename(zip_file)
    date = re.findall(r'\d{8}', file_name)[0]

    return date


def get_dates_from_zips(zip_files):
    """Get dates from S1 zip files

    Args:
        zip_files (str list): S1 zip files

    Returns:
        str list: sorted S1 dates
    """
    dates = []
    for file in zip_files:
        dates.append(get_date_from_zip(file))
    dates = sorted(list(set(dates)))

    return dates


def get_date_from_par(slc_par):
    """Get date from slc par

    Args:
        slc_par (str): slc par file

    Returns:
        str: date
    """
    with open(slc_par, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            if line.strip().startswith('date'):
                tmp = line.strip().split(':')[1].strip().split()
                y, m, d = tmp[0], tmp[1], tmp[2]
                if len(m) == 1:
                    m = '0' + m
                if len(d) == 1:
                    d = '0' + d
                date = y + m + d

    return date


def get_sensor_from_par(slc_par):
    """Get sensor name from slc par

    Args:
        slc_par (str): slc par file

    Returns:
        str: sensor name, S1A or S1B
    """
    with open(slc_par, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            if line.strip().startswith('sensor'):
                if 'S1A' in line:
                    sensor = 'S1A'
                else:
                    sensor = 'S1B'

    return sensor


def date_math(date, day):
    """date operation (add or subtract day)

    Args:
        date (str): date
        day (int): day

    Returns:
        str: operated date
    """
    y = int(date[0:4])
    m = int(date[4:6])
    d = int(date[6:8])
    in_date = datetime.datetime(y, m, d)
    out_date = in_date + datetime.timedelta(days=day)
    out_date = out_date.strftime('%Y%m%d')

    return out_date


def get_satellite_name(zip_file):
    """Get satellite name from S1 zip file

    Args:
        zip_file (str): S1 zip file

    Returns:
        str: satellite name
    """
    file_name = os.path.basename(zip_file)
    if file_name.startswith('S1A_IW_SLC'):
        name = 'S1A'
    else:
        name = 'S1B'

    return name


def get_orbit_file(s1_date, sensor, orbit_dir):
    """Get orbit file using S1 date and sensor name

    Args:
        s1_date (str): S1 date
        sensor (str): satellite name, S1A or S1B
        orbit_dir (str): directory including orbit files

    Returns:
        str: orbit file including path
    """
    orbit_file = None
    orbit_date = date_math(s1_date, -1)
    for orb in os.listdir(orbit_dir):
        if orb.endswith('EOF') and orb.startswith(sensor):
            dates = re.findall(r'\d{8}', orb)
            if orbit_date == dates[-2]:
                orbit_file = os.path.join(orbit_dir, orb)

    return orbit_file


def lookup_from_list(str_list, pol, sub_swath):
    """lookup list element using pol and sub_swath

    Args:
        str_list (list): list
        pol (str): polarization
        sub_swath (str): sub_swath, iw1 iw2 iw3

    Returns:
        [str]: result
    """
    res = None
    for i in str_list:
        if pol in i and sub_swath in i:
            res = i

    return res


def s1_import_slc_from_zip(zip_file, pol, sub_swath, out_slc_dir):
    """Generate S1 TOPS SLC from zip file

    Args:
        zip_file (str): zip file
        pol (str list): polarization from ['vv', 'vh']
        sub_swath ([int list]): sub_swath number from [1, 2, 3]
        out_slc_dir (str): directory for saving SLC
    """
    # unzip file
    geo_tiffs = []
    annotation_xmls = []
    calibration_xmls = []
    noise_xmls = []
    all_files = []
    with zipfile.ZipFile(zip_file, mode='r') as f:
        files = f.namelist()
        for file in files:
            for p in pol:
                for i in sub_swath:
                    iw = 'iw' + str(i)
                    if file.endswith('tiff') and p in file and iw in file:
                        geo_tiffs.append(file)
                        all_files.append(file)
                    if file.endswith('xml') and p in file and iw in file:
                        if 'noise' in file:
                            noise_xmls.append(file)
                            all_files.append(file)
                        elif 'calibration' in file:
                            calibration_xmls.append(file)
                            all_files.append(file)
                        else:
                            annotation_xmls.append(file)
                            all_files.append(file)
        for file in all_files:
            f.extract(file, out_slc_dir)
    # generate slc
    date = get_date_from_zip(zip_file)
    for p in pol:
        for i in sub_swath:
            iw = 'iw' + str(i)
            slc = os.path.join(out_slc_dir, f'{date}.{iw}.{p}.slc')
            slc_par = slc + '.par'
            tops_par = slc + '.tops_par'

            geo_tiff = lookup_from_list(geo_tiffs, p, iw)
            annotation_xml = lookup_from_list(annotation_xmls, p, iw)
            calibration_xml = lookup_from_list(calibration_xmls, p, iw)
            noise_xml = lookup_from_list(noise_xmls, p, iw)

            geo_tiff = os.path.join(out_slc_dir, geo_tiff)
            annotation_xml = os.path.join(out_slc_dir, annotation_xml)
            calibration_xml = os.path.join(out_slc_dir, calibration_xml)

            if noise_xml:
                noise_xml = os.path.join(out_slc_dir, noise_xml)
                cmd_str = f'par_S1_SLC {geo_tiff} {annotation_xml} {calibration_xml} {noise_xml} {slc_par} {slc} {tops_par}'
            else:
                cmd_str = f'par_S1_SLC {geo_tiff} {annotation_xml} {calibration_xml} {slc_par} {slc} {tops_par}'
            os.system(cmd_str)
    # delete files
    shutil.rmtree(glob.glob(os.path.join(out_slc_dir, 'S1*SAFE'))[0])


def orbit_correct(slc_par, orbit_dir):
    """Extract Sentinel-1 OPOD state vectors and copy into the ISP image parameter file

    Args:
        slc_par (str): slc par
        orbit_dir (str): dectory including orbits

    Returns:
        str: orbit correct result
    """
    s1_date = get_date_from_par(slc_par)
    sensor = get_sensor_from_par(slc_par)
    orbit_file = get_orbit_file(s1_date, sensor, orbit_dir)
    if orbit_file:
        call_str = f'S1_OPOD_vec.vec {slc_par} {orbit_file}'
        os.system(call_str)
        return s1_date + ' YES'
    else:
        return s1_date + ' No'


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


def run_all_steps(zip_file, date_slc_dir, orb_dir, pol, sub_swath, rlks, alks):
    """run all step(unzip file, generate slc, orbit correct, generate bmp)

    Args:
        zip_file (str): zip file
        date_slc_dir (str): directory saving slc
        orb_dir (str): orbit directory
        pol (str list): polarization
        sub_swath (int list): sub_swath number
        rlks (int): range looks
        alks (int): azimuth looks

    Returns:
        str list: orbit correct result
    """
    orb_correct_res = []
    # unzip file and generate slc
    s1_import_slc_from_zip(zip_file, pol, sub_swath, date_slc_dir)
    # orbit correct and slc2bmp
    slcs = glob.glob(os.path.join(date_slc_dir, '*iw*slc'))
    for slc in slcs:
        slc_par = slc + '.par'
        res = orbit_correct(slc_par, orb_dir)
        orb_correct_res.append(res)
        bmp = slc + '.bmp'
        slc2bmp(slc, slc_par, rlks, alks, bmp)

    return orb_correct_res


def print_orb_correct_report(orb_correct_res):
    """print orbit correct report

    Args:
        orb_correct_res (str list): orbit correct report result
    """
    print('-' * 31)
    title = 'orbit correction report'
    print('|' + title.center(29) + '|')
    print('-' * 31)
    for res in orb_correct_res:
        date, flag = res.split()
        print('|' + date.center(14) + '|' + flag.center(14) + '|')
        print('-' * 31)


def main():
    # get inputs
    inps = cmdLineParse()
    zip_dir = os.path.abspath(inps.s1_zip_dir)
    orb_dir = os.path.abspath(inps.orbit_dir)
    slc_dir = os.path.abspath(inps.slc_dir)
    pol = inps.pol
    sub_swath = inps.sub_swath
    rlks = inps.rlks
    alks = inps.alks
    # check input
    if not os.path.isdir(zip_dir):
        sys.exit(f'{zip_dir} does not exist.')
    if not os.path.isdir(orb_dir):
        sys.exit(f'{orb_dir} does not exist.')
    if not os.path.isdir(slc_dir):
        os.mkdir(slc_dir)
    # get all raw files
    zip_files = glob.glob(os.path.join(zip_dir, 'S1*_IW_SLC*zip'))
    if len(zip_files) == 0:
        sys.exit(f'Cannot find zip files in {zip_dir}.')
    # get all dates
    dates = get_dates_from_zips(zip_files)
    orb_correct_res = []
    for date in dates:
        same_date_zips = glob.glob(
            os.path.join(zip_dir, f'S1*_IW_SLC*{date}*zip'))
        same_date_zips = sorted(same_date_zips, key=lambda i: i[26:32])
        # one date --> one zip
        if len(same_date_zips) == 1:
            zip_file = same_date_zips[0]
            # directory saving slc
            date_slc_dir = os.path.join(slc_dir, date)
            if not os.path.isdir(date_slc_dir):
                os.mkdir(date_slc_dir)
            # run all step
            res = run_all_steps(zip_file, date_slc_dir, orb_dir, pol,
                                sub_swath, rlks, alks)
            orb_correct_res += res
        # one date --> multi zips
        else:
            for i in range(len(same_date_zips)):
                zip_file = same_date_zips[i]
                # directory saving slc
                date_slc_dir = os.path.join(slc_dir, f'{date}-{i+1}')
                if not os.path.isdir(date_slc_dir):
                    os.mkdir(date_slc_dir)
                # run all step
                res = run_all_steps(zip_file, date_slc_dir, orb_dir, pol,
                                    sub_swath, rlks, alks)
                orb_correct_res += res

    print_orb_correct_report(orb_correct_res)

    print('\nAll done, enjoy it!\n')


if __name__ == '__main__':
    main()

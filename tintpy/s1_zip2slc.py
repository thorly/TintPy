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
  python3 s1_zip2slc.py /ly/zips /ly/orbits /ly/slc 1 -r 20 -a 5
  python3 s1_zip2slc.py /ly/zips /ly/orbits /ly/slc 1 2 3 -r 20 -a 5
"""


def cmdLineParse():
    parser = argparse.ArgumentParser(
        description='Generate SLC from Sentinel-1 raw data with orbit correction using GAMMA.',
        formatter_class=argparse.RawTextHelpFormatter,epilog=EXAMPLE)

    parser.add_argument('s1_zip_dir', help='Sentinel-1 raw directory')
    parser.add_argument('orbit_dir', help='precise orbit directory')
    parser.add_argument('slc_dir', help='directory for saving slc')
    parser.add_argument('sub_swath', help='sub_swath number', type=int, choices=[1, 2, 3], nargs='+')
    parser.add_argument('-r', dest='rlks', help='range looks', type=int, required=True)
    parser.add_argument('-a', dest='alks', help='azimuth looks', type=int, required=True)
    parser.add_argument('-p', dest='pol', help='polarization (defaults: vv)', nargs='+', choices=['vv', 'vh'], default=['vv'])
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


def gen_number_table(tops_par_files, out_file):
    """Generate number table for each iw

    Args:
        tops_par_files (list): tops par files
        out_file (str): output file
    """
    num = None
    time = []

    f_out = open(out_file, 'w+', encoding='utf-8')

    for tops_par_file in tops_par_files:
        name = os.path.basename(tops_par_file)
        iw = name.split('.')[1]

        with open(tops_par_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            for line in lines:
                if 'number_of_bursts' in line:
                    num = line.strip().split()[1]
                if 'burst_asc_node_' in line:
                    time.append(line.strip().split()[1])

        if num and time:
            f_out.write(f"{iw}_number_of_bursts: {num}\n")
            f_out.write(f"{iw}_first_burst:      {time[0]}\n")
            f_out.write(f"{iw}_last_burst:       {time[-1]}\n")
        num = None
        time = []

    f_out.close()


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
    date = get_date_from_zip(zip_file)
    print(f"unzip necessary files for date: {date}")
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
                cmd_str = f'par_S1_SLC {geo_tiff} {annotation_xml} {calibration_xml} - {slc_par} {slc} {tops_par}'
            os.system(cmd_str)

    tops_pars = glob.glob(os.path.join(out_slc_dir, f'*{pol[0]}.slc.tops_par'))
    out_file = os.path.join(out_slc_dir, f'{date}.burst_number_table')
    gen_number_table(tops_pars, out_file)
    # delete files
    shutil.rmtree(glob.glob(os.path.join(out_slc_dir, 'S1*SAFE'))[0])


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
    """
    # unzip file and generate slc
    s1_import_slc_from_zip(zip_file, pol, sub_swath, date_slc_dir)
    # orbit correct and slc2bmp
    slcs = glob.glob(os.path.join(date_slc_dir, '*iw*slc'))
    os.chdir(date_slc_dir)
    for slc in slcs:
        slc_par = slc + '.par'
        # orbit correct
        call_str = f'OPOD_vec {slc_par} {orb_dir}'
        os.system(call_str)
        # slc2bmp
        slc2bmp(slc, slc_par, rlks, alks, slc + '.bmp')


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

    for date in dates:
        same_date_zips = glob.glob(os.path.join(zip_dir, f'S1*_IW_SLC*{date}*zip'))
        same_date_zips = sorted(same_date_zips, key=lambda i: i[26:32])

        # one date --> one zip
        if len(same_date_zips) == 1:
            zip_file = same_date_zips[0]
            # directory saving slc
            date_slc_dir = os.path.join(slc_dir, date)
            if not os.path.isdir(date_slc_dir):
                os.mkdir(date_slc_dir)
            # run all step
            run_all_steps(zip_file, date_slc_dir, orb_dir, pol, sub_swath, rlks, alks)

        # one date --> multi zips
        else:
            for i in range(len(same_date_zips)):
                zip_file = same_date_zips[i]
                # directory saving slc
                date_slc_dir = os.path.join(slc_dir, f'{date}-{i+1}')
                if not os.path.isdir(date_slc_dir):
                    os.mkdir(date_slc_dir)
                # run all step
                run_all_steps(zip_file, date_slc_dir, orb_dir, pol, sub_swath, rlks, alks)

    print('\nAll done, enjoy it!\n')


if __name__ == '__main__':
    main()

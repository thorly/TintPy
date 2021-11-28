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

EXAMPLE = """Example:
  python3 alos_raw2slc.py /ly/ALOS_PALSAR_raw /ly/slc_fbd palsar_ant_20061024.dat 2 10
  python3 alos_raw2slc.py /ly/ALOS_PALSAR_raw /ly/slc_fbs palsar_ant_20061024.dat 4 6
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(description='Convert ALOS PALSAR raw data (level 1.0) to SLC data',
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)
    parser.add_argument('data_dir', help='directory including ALOS level 1.0 data (**uncompress data first**)')
    parser.add_argument('output_dir', help='directory for saving ALOS SLC')
    parser.add_argument('antenna_file', help='ALOS PALSAR JAXA antenna pattern file')
    parser.add_argument('rlks', help='range looks', type=int)
    parser.add_argument('alks', help='azimuth looks', type=int)

    inps = parser.parse_args()
    return inps


def get_raw_names(data_dir):
    names = []
    files = os.listdir(data_dir)
    for file in files:
        if re.search(r'ALPSRP.*L1\.0', file):
            if os.path.isdir(os.path.join(data_dir, file)):
                names.append(file)
    return names


def get_date(workreport):
    """get date for ALOS PALSAR"""
    date = None
    with open(workreport, 'r') as f:
        for line in f.readlines():
            if line.startswith('Img_SceneCenterDateTime'):
                res = re.findall(r'\d{8}', line)
                if res:
                    date = res[0]
    return date


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


def main():
    inps = cmdline_parser()
    data_dir = inps.data_dir
    output_dir = inps.output_dir
    antenna_file = inps.antenna_file
    rlks = inps.rlks
    alks = inps.alks

    # check
    data_dir = os.path.abspath(data_dir)
    output_dir = os.path.abspath(output_dir)
    antenna_file = os.path.abspath(antenna_file)
    if not os.path.isdir(data_dir):
        sys.exit('Cannot find this directory: {}'.format(data_dir))
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    if not os.path.isfile(antenna_file):
        sys.exit('Cannot find file: {}'.format(antenna_file))

    # run
    for name in get_raw_names(data_dir):
        workreport_file = os.path.join(data_dir, name, 'workreport')
        date = get_date(workreport_file)

        # create date directory
        slc_dir = os.path.join(output_dir, date)
        if not os.path.isdir(slc_dir):
            os.mkdir(slc_dir)

        os.chdir(slc_dir)

        LED_file = glob.glob(os.path.join(data_dir, name, 'LED*'))[0]
        IMG_file = glob.glob(os.path.join(data_dir, name, 'IMG-HH*'))[0]

        cmd_str = f"PALSAR_proc {LED_file} palsar.par p{date}.slc.par {IMG_file} {date}.raw 0 0"
        os.system(cmd_str)

        log_file = os.path.join(slc_dir, date + '.log')

        os.system(f"echo {cmd_str} > {log_file}")

        cmd_str = f"PALSAR_antpat palsar.par p{date}.slc.par {antenna_file} PALSAR_antpat_MSP.dat"
        os.system(cmd_str)
        os.system(f"echo {cmd_str} >> {log_file}")

        cmd_str = f"dop_mlcc palsar.par p{date}.slc.par {date}.raw dop_mlcc.dat -"
        os.system(cmd_str)
        os.system(f"echo {cmd_str} >> {log_file}")

        cmd_str = f"azsp_IQ palsar.par p{date}.slc.par {date}.raw {date}.azsp 0 - - 1 - 0"
        os.system(cmd_str)
        os.system(f"echo {cmd_str} >> {log_file}")

        cmd_str = f"doppler palsar.par p{date}.slc.par {date}.raw {date}.dop - - 0 - - 0"
        os.system(cmd_str)
        os.system(f"echo {cmd_str} >> {log_file}")

        cmd_str = f"rspec_IQ palsar.par p{date}.slc.par {date}.raw {date}.rspec - - - - 0"
        os.system(cmd_str)
        os.system(f"echo {cmd_str} >> {log_file}")

        cmd_str = f"pre_rc palsar.par p{date}.slc.par {date}.raw {date}.rc - - - - - - - - - - 1 - - -"
        os.system(cmd_str)
        os.system(f"echo {cmd_str} >> {log_file}")

        for i in range(6):
            cmd_str = f"autof palsar.par p{date}.slc.par {date}.rc {date}.autof 5.0"
            os.system(cmd_str)
            os.system(f"echo {cmd_str} >> {log_file}")

        cmd_str = f"az_proc palsar.par p{date}.slc.par {date}.rc {date}.slc 16384 0 -49.8 0 2.12"
        os.system(cmd_str)
        os.system(f"echo {cmd_str} >> {log_file}")

        cmd_str = f"par_MSP palsar.par p{date}.slc.par {date}.slc.par"
        os.system(cmd_str)
        os.system(f"echo {cmd_str} >> {log_file}")

        slc2bmp(f"{date}.slc", f"{date}.slc.par", rlks, alks, f"{date}.slc.bmp")

        # clean directory
        save_files = []
        save_files.append(date + '.slc')
        save_files.append(date + '.slc.par')
        save_files.append(date + '.slc.bmp')
        save_files.append(date + '.log')

        for file in os.listdir(slc_dir):
            if file not in save_files:
                os.remove(file)

        print('\nAll done, enjoy it!\n')

if __name__ == "__main__":
    main()

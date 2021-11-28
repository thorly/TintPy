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
  python3 alos2slc.py /ly/ALOS_PALSAR /ly/slc 1 2 10
  python3 alos2slc.py /ly/ALOS_PALSAR2 /ly/slc 2 8 16
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(
        description='Reformat EORC PALSAR + PALSAR2 level 1.1 CEOS format SLC data and generate the ISP parameter file',
        formatter_class=argparse.RawTextHelpFormatter,epilog=EXAMPLE)
    parser.add_argument('data_dir', help='directory including ALOS level 1.1 data (**uncompress data first**)')
    parser.add_argument('output_dir', help='directory saving SLC')
    parser.add_argument('flag', choices=[1, 2], help='flag for data type (1 for PALSAR, 2 for PALSAR2)', type=int)
    parser.add_argument('rlks', help='range looks', type=int)
    parser.add_argument('alks', help='azimuth looks', type=int)

    inps = parser.parse_args()
    return inps


def get_slc_dir_names(data_dir, flag):
    names = []
    files = os.listdir(data_dir)
    # ALOS PALSAR
    if flag == 1:
        for file in files:
            if re.search(r'ALPSRP.*', file):
                if os.path.isdir(os.path.join(data_dir, file)):
                    names.append(file)
    # ALOS PALSAR2
    if flag == 2:
        for file in files:
            if re.search(r'\d{10}_\d{6}_ALOS\d{10}-\d{6}', file):
                if os.path.isdir(os.path.join(data_dir, file)):
                    names.append(file)
    return names


def get_date(summary_file):
    """get date for ALOS PALSAR"""
    date = None
    with open(summary_file, 'r') as f:
        for line in f.readlines():
            if line.startswith('Lbi_ObservationDate'):
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


def reformat_alos(data_dir, output_dir, flag, rlks, alks):
    slc_dir_names = get_slc_dir_names(data_dir, flag)
    # reformat ALOS data
    if slc_dir_names:
        for name in slc_dir_names:
            slc_dir = os.path.join(data_dir, name)
            if flag == 1:
                img = glob.glob(os.path.join(slc_dir, 'IMG-HH-ALPSRP*'))
                led = glob.glob(os.path.join(slc_dir, 'LED-ALPSRP*'))

                summary_file = os.path.join(slc_dir, 'summary.txt')
                date = get_date(summary_file)

                if date:
                    date_dir = os.path.join(output_dir, date)
                    if not os.path.isdir(date_dir):
                        os.mkdir(date_dir)

                    if os.path.isfile(img[0]) and os.path.isfile(led[0]):
                        out_slc = os.path.join(date_dir, date + '.slc')
                        out_slc_par = out_slc + '.par'

                        cmd_str = f"par_EORC_PALSAR {led[0]} {out_slc_par} {img[0]} {out_slc}"
                        os.system(cmd_str)

                        slc2bmp(out_slc, out_slc_par, rlks, alks,
                                out_slc + '.bmp')
                else:
                    sys.exit('Cannot find date in summary.txt')
            if flag == 2:
                img = glob.glob(os.path.join(slc_dir, 'IMG-HH-ALOS*'))
                led = glob.glob(os.path.join(slc_dir, 'LED-ALOS*'))

                date = '20' + slc_dir[-6:]
                date_dir = os.path.join(output_dir, date)

                if not os.path.isdir(date_dir):
                    os.mkdir(date_dir)

                if os.path.isfile(img[0]) and os.path.isfile(led[0]):
                    out_slc = os.path.join(date_dir, date + '.slc')
                    out_slc_par = out_slc + '.par'

                    cmd_str = f"par_EORC_PALSAR {led[0]} {out_slc_par} {img[0]} {out_slc}"
                    os.system(cmd_str)

                    slc2bmp(out_slc, out_slc_par, rlks, alks, out_slc + '.bmp')
    else:
        sys.exit('Cannot find any slc in {}'.format(output_dir))


def main():
    inps = cmdline_parser()
    data_dir = inps.data_dir
    output_dir = inps.output_dir
    flag = inps.flag
    rlks = inps.rlks
    alks = inps.alks
    # check
    data_dir = os.path.abspath(data_dir)
    output_dir = os.path.abspath(output_dir)
    if not os.path.isdir(data_dir):
        sys.exit('cannot find this directory: {}'.format(data_dir))
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    # run
    reformat_alos(data_dir, output_dir, flag, rlks, alks)

    print('\nAll done, enjoy it!\n')


if __name__ == "__main__":
    main()

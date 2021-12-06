#!/usr/bin/env python3
################################
# Program is part of TintPy    #
# Copyright (c) 2021, Lei Yuan #
# Author: Lei Yuan, 2021       #
################################

import argparse
import os
import re
import sys

EXAMPLE = """Example:
  python3 deramp_s1_b.py /ly/slc /ly/slc_deramp 2 -r 8 -a 2
  python3 deramp_s1_b.py /ly/slc /ly/slc_deramp 1 2 -r 8 -a 2
"""


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


def write_tab(iw_slc_path, out_file):
    """Write path of slc slc_par and tops_par to tab file

    Args:
        iw_slc_path (str): slc file
        out_file (str): output file
    """
    with open(out_file, 'w+') as f:
        f.write(
            f"{iw_slc_path} {iw_slc_path + '.par'} {iw_slc_path + '.tops_par'}\n"
        )


def cmd_line_parser():
    parser = argparse.ArgumentParser(
        description=
        'Calculate and subtract S1 TOPS Doppler phase from burst SLC data before coregistration.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLE)

    parser.add_argument('slc_dir', help='directory path of SLCs')
    parser.add_argument('out_dir', help='directory path of deramped SLCs')
    parser.add_argument('iw_num',
                        choices=['1', '2', '3'],
                        nargs='+',
                        type=str,
                        help='IW num for deramp')
    parser.add_argument('-r',
                        dest='rlks',
                        help='range looks',
                        required=True,
                        type=int)
    parser.add_argument('-a',
                        dest='alks',
                        help='azimuth looks',
                        required=True,
                        type=int)
    inps = parser.parse_args()

    return inps


def main():
    inps = cmd_line_parser()
    slc_dir = inps.slc_dir
    out_dir = inps.out_dir
    iw_num = inps.iw_num
    rlks = inps.rlks
    alks = inps.alks

    # check slc_dir
    slc_dir = os.path.abspath(slc_dir)
    if not os.path.isdir(slc_dir):
        sys.exit('Cannot find directory {}'.format(slc_dir))

    # check out_dir
    out_dir = os.path.abspath(out_dir)
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # get all dates
    dates = sorted(
        [i for i in os.listdir(slc_dir) if re.findall(r'^\d{8}$', i)])

    for date in dates:
        slc_date_dir = os.path.join(slc_dir, date)

        out_date_dir = os.path.join(out_dir, date)
        if not os.path.isdir(out_date_dir):
            os.mkdir(out_date_dir)

        tab = os.path.join(out_date_dir, date + '_tab')
        tab_deramp = os.path.join(out_date_dir, date + '_tab.deramp')

        if len(iw_num) == 1:
            iw_slc = os.path.join(slc_date_dir, f"{date}.iw{iw_num[0] * 2}.slc")

            iw_slc_deramp = os.path.join(out_date_dir, f"{date}.iw{iw_num[0] * 2}.slc")

            write_tab(iw_slc, tab)

            write_tab(iw_slc_deramp, tab_deramp)

        if len(iw_num) == 2:
            iw_slc1 = os.path.join(slc_date_dir, f"{date}.iw{iw_num[0] * 2}.slc")
            iw_slc2 = os.path.join(slc_date_dir, f"{date}.iw{iw_num[1] * 2}.slc")

            iw_slc1_deramp = os.path.join(out_date_dir, f"{date}.iw{iw_num[0] * 2}.slc")
            iw_slc2_deramp = os.path.join(out_date_dir, f"{date}.iw{iw_num[1] * 2}.slc")

            write_tab(iw_slc1, tab)
            write_tab(iw_slc2, tab)

            write_tab(iw_slc1_deramp, tab_deramp)
            write_tab(iw_slc2_deramp, tab_deramp)

        if len(iw_num) == 3:
            iw_slc1 = os.path.join(slc_date_dir, f"{date}.iw{iw_num[0] * 2}.slc")
            iw_slc2 = os.path.join(slc_date_dir, f"{date}.iw{iw_num[1] * 2}.slc")
            iw_slc3 = os.path.join(slc_date_dir, f"{date}.iw{iw_num[2] * 2}.slc")

            iw_slc1_deramp = os.path.join(out_date_dir, f"{date}.iw{iw_num[0] * 2}.slc")
            iw_slc2_deramp = os.path.join(out_date_dir, f"{date}.iw{iw_num[1] * 2}.slc")
            iw_slc3_deramp = os.path.join(out_date_dir, f"{date}.iw{iw_num[2] * 2}.slc")

            write_tab(iw_slc1, tab)
            write_tab(iw_slc2, tab)
            write_tab(iw_slc3, tab)

            write_tab(iw_slc1_deramp, tab_deramp)
            write_tab(iw_slc2_deramp, tab_deramp)
            write_tab(iw_slc3_deramp, tab_deramp)

        os.chdir(out_date_dir)

        call_str = f"SLC_deramp_S1_TOPS {tab} {tab_deramp} 0 1"
        os.system(call_str)

        slc_deramp = os.path.join(out_date_dir, date + ".slc")
        slc_deramp_par = slc_deramp + '.par'

        call_str = f"SLC_mosaic_S1_TOPS {tab_deramp} {slc_deramp} {slc_deramp_par} {rlks} {alks} 1"
        os.system(call_str)

        slc2bmp(slc_deramp, slc_deramp_par, rlks, alks, slc_deramp + '.bmp')

    print('\nAll done, enjoy it!\n')


if __name__ == "__main__":
    main()

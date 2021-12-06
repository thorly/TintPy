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
  python3 deramp_s1_a.py /ly/rslc 2 20181229
  python3 deramp_s1_a.py /ly/rslc 1 2 20181229 -r 8 -a 2
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
        'Calculate and subtract S1 TOPS Doppler phase from burst SLC data after coregistration.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLE)

    parser.add_argument('rslc_dir', help='directory path of RSLCs')
    parser.add_argument('iw_num',
                        choices=['1', '2', '3'],
                        nargs='+',
                        type=str,
                        help='IW num for deramp')
    parser.add_argument('ref_slc', help='reference SLC')
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
    rslc_dir = inps.rslc_dir
    iw_num = inps.iw_num
    iw_num = iw_num.split()
    rlks = inps.rlks
    alks = inps.alks
    ref_slc = inps.ref_slc

    # check rslc_dir
    rslc_dir = os.path.abspath(rslc_dir)
    if not os.path.isdir(rslc_dir):
        sys.exit('Cannot find directory {}'.format(rslc_dir))

    # get all date
    all_date = sorted(
        [i for i in os.listdir(rslc_dir) if re.findall(r'^\d{8}$', i)])
    if len(all_date) < 2:
        sys.exit('No enough RSLCs.')

    # check ref_slc
    if re.findall(r'^\d{8}$', ref_slc):
        if not ref_slc in all_date:
            sys.exit('No slc for {}.'.format(ref_slc))
    else:
        sys.exit('Error date for ref_slc.')

    s_dates = all_date
    s_dates.remove(ref_slc)

    dates = [ref_slc] + s_dates

    for date in dates:
        rslc_date_dir = os.path.join(rslc_dir, date)
        tab = os.path.join(rslc_date_dir, date + '_tab')

        if len(iw_num) == 1:
            iw_slc = os.path.join(rslc_date_dir, f"{date}.iw{iw_num[0] * 2}.slc")

            write_tab(iw_slc, tab)

        if len(iw_num) == 2:
            iw_slc1 = os.path.join(rslc_date_dir, f"{date}.iw{iw_num[0] * 2}.slc")
            iw_slc2 = os.path.join(rslc_date_dir, f"{date}.iw{iw_num[1] * 2}.slc")

            write_tab(iw_slc1, tab)
            write_tab(iw_slc2, tab)

        if len(iw_num) == 3:
            iw_slc1 = os.path.join(rslc_date_dir, f"{date}.iw{iw_num[0] * 2}.slc")
            iw_slc2 = os.path.join(rslc_date_dir, f"{date}.iw{iw_num[1] * 2}.slc")
            iw_slc3 = os.path.join(rslc_date_dir, f"{date}.iw{iw_num[2] * 2}.slc")

            write_tab(iw_slc1, tab)
            write_tab(iw_slc2, tab)
            write_tab(iw_slc3, tab)

        if date == ref_slc:
            call_str = f"S1_deramp_TOPS_reference {tab}"
            os.system(call_str)
        else:
            reference_tab = os.path.join(rslc_dir, ref_slc, ref_slc + '_tab')
            call_str = f"S1_deramp_TOPS_slave {tab} {date} {reference_tab} {rlks} {alks} 0"
            os.system(call_str)

        rslc_deramp = os.path.join(rslc_date_dir, date + ".rslc.deramp")
        rslc_deramp_par = rslc_deramp + '.par'

        if date == ref_slc:
            call_str = f"SLC_mosaic_S1_TOPS {date}_tab.deramp {rslc_deramp} {rslc_deramp_par} {rlks} {alks} 1"
            os.system(call_str)

        slc2bmp(rslc_deramp, rslc_deramp_par, rlks, alks, rslc_deramp + '.bmp')

    print('\nAll done, enjoy it!\n')


if __name__ == "__main__":
    main()

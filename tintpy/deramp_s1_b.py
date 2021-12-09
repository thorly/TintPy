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


def write_tab(slc_date_dir, sub_swath, pol, out_file):
    """Write path of slc slc_par and tops_par to tab file

    Args:
        slc_date_dir (str): slc directory
        sub_swath (list): sub_swath number
        pol (str): polarization
        out_file (str): output file
    """
    date = os.path.basename(slc_date_dir)[0:8]
    with open(out_file, 'w+', encoding='utf-8') as f:
            for i in sub_swath:
                iw_rslc = os.path.join(slc_date_dir, f'{date}.iw{i}.{pol}.slc')
                iw_rslc_par = iw_rslc + '.par'
                iw_rslc_tops_par = iw_rslc + '.tops_par'
                f.write(f'{iw_rslc} {iw_rslc_par} {iw_rslc_tops_par}\n')


def cmd_line_parser():
    parser = argparse.ArgumentParser(
        description='Calculate and subtract S1 TOPS Doppler phase from burst SLC data before coregistration.',
        formatter_class=argparse.RawTextHelpFormatter,epilog=EXAMPLE)

    parser.add_argument('slc_dir', help='directory path of SLCs')
    parser.add_argument('out_dir', help='directory path of deramped SLCs')
    parser.add_argument('sub_swath',
                        type=int,
                        nargs='+',
                        choices=[1, 2, 3],
                        help='sub_swath number for deramp')
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
    parser.add_argument('-p',
                        dest='pol',
                        help='polarization(defaults: vv)',
                        choices=['vv', 'vh'],
                        default='vv')
    inps = parser.parse_args()

    return inps


def main():
    inps = cmd_line_parser()
    slc_dir = os.path.abspath(inps.slc_dir)
    out_dir = os.path.abspath(inps.out_dir)
    sub_swath = inps.sub_swath
    rlks = inps.rlks
    alks = inps.alks
    pol = inps.pol

    # check slc_dir
    if not os.path.isdir(slc_dir):
        sys.exit('{} does not exist.'.format(slc_dir))

    # check out_dir
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # get all dates
    dates = sorted([i for i in os.listdir(slc_dir) if re.findall(r'^\d{8}$', i)])

    if len(dates) < 1:
        sys.exit('No SLC.')

    for date in dates:
        slc_date_dir = os.path.join(slc_dir, date)

        out_date_dir = os.path.join(out_dir, date)
        if not os.path.isdir(out_date_dir):
            os.mkdir(out_date_dir)

        tab = os.path.join(out_date_dir, date + '_tab')
        tab_deramp = os.path.join(out_date_dir, date + '_tab.deramp')

        write_tab(slc_date_dir, sub_swath, pol, tab)
        write_tab(out_date_dir, sub_swath, pol, tab_deramp)

        os.chdir(out_date_dir)

        call_str = f"SLC_deramp_S1_TOPS {tab} {tab_deramp} 0 1"
        os.system(call_str)

        slc_deramp = os.path.join(out_date_dir, date + ".slc")
        slc_deramp_par = slc_deramp + '.par'

        call_str = f"SLC_mosaic_S1_TOPS {tab_deramp} {slc_deramp} {slc_deramp_par} {rlks} {alks} 1"
        os.system(call_str)

        os.remove(tab)
        os.remove(tab_deramp)

        slc2bmp(slc_deramp, slc_deramp_par, rlks, alks, slc_deramp + '.bmp')

    print('\nAll done, enjoy it!\n')


if __name__ == "__main__":
    main()

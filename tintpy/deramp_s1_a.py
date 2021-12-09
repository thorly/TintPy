#!/usr/bin/env python3
################################
# Program is part of TintPy    #
# Copyright (c) 2021, Lei Yuan #
# Author: Lei Yuan, 2021       #
################################

import argparse
import os
import re
import shutil
import sys

EXAMPLE = """Example:
  python3 deramp_s1_a.py rslc rslc_deramp 20211229 3 -r 32 -a 8
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


def cmd_line_parser():
    parser = argparse.ArgumentParser(
        description=
        'Calculate and subtract S1 TOPS Doppler phase from burst SLC data after coregistration.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLE)

    parser.add_argument('rslc_dir', help='RSLCs directory')
    parser.add_argument('out_dir', help='output directory')
    parser.add_argument('ref_slc', help='reference SLC')
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
    rslc_dir = os.path.abspath(inps.rslc_dir)
    out_dir = os.path.abspath(inps.out_dir)
    sub_swath = inps.sub_swath
    rlks = inps.rlks
    alks = inps.alks
    ref_slc = inps.ref_slc
    pol = inps.pol

    # check rslc_dir
    if not os.path.isdir(rslc_dir):
        sys.exit('{} does not exist.'.format(rslc_dir))

    # check out_dir
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # get all date
    all_date = sorted([i for i in os.listdir(rslc_dir) if re.findall(r'^\d{8}$', i)])

    if len(all_date) < 2:
        sys.exit('No enough RSLCs.')

    # check ref_slc
    if re.findall(r'^\d{8}$', ref_slc):
        if not ref_slc in all_date:
            sys.exit('No slc for {}.'.format(ref_slc))
    else:
        sys.exit('Error date for ref_slc.')

    s_dates = all_date.copy()
    s_dates.remove(ref_slc)

    dates = [ref_slc] + s_dates

    for date in dates:
        rslc_date_dir = os.path.join(rslc_dir, date)

        tab = os.path.join(rslc_date_dir, date + '_tab')

        with open(tab, 'w+', encoding='utf-8') as f:
            for i in sub_swath:
                iw_rslc = os.path.join(rslc_date_dir, f'{date}.iw{i}.{pol}.rslc')
                iw_rslc_par = iw_rslc + '.par'
                iw_rslc_tops_par = iw_rslc + '.tops_par'
                f.write(f'{iw_rslc} {iw_rslc_par} {iw_rslc_tops_par}\n')

        os.chdir(rslc_date_dir)

        rslc_deramp = os.path.join(rslc_date_dir, date + ".rslc.deramp")
        rslc_deramp_par = rslc_deramp + '.par'

        if date == ref_slc:
            call_str = f"S1_deramp_TOPS_reference {tab}"
            os.system(call_str)

            call_str = f"SLC_mosaic_S1_TOPS {date}_tab.deramp {rslc_deramp} {rslc_deramp_par} {rlks} {alks} 1"
            os.system(call_str)
        else:
            reference_tab = os.path.join(rslc_dir, ref_slc, ref_slc + '_tab')
            call_str = f"S1_deramp_TOPS_slave {tab} {date} {reference_tab} {rlks} {alks} 1"
            os.system(call_str)

        slc2bmp(rslc_deramp, rslc_deramp_par, rlks, alks, rslc_deramp + '.bmp')

        out_date_dir = os.path.join(out_dir, date)
        if not os.path.isdir(out_date_dir):
            os.mkdir(out_date_dir)

        shutil.move(rslc_deramp, out_date_dir)
        shutil.move(rslc_deramp_par, out_date_dir)
        shutil.move(rslc_deramp + '.bmp', out_date_dir)

    for date in dates:
        tab = os.path.join(rslc_dir, date, date + '_tab')
        tab_deramp = tab + '.deramp'
        if date == ref_slc:
            with open(tab_deramp, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                for line in lines:
                    if line.strip():
                        for i in line.strip().split():
                            if os.path.isfile(i):
                                os.remove(i)
                            if os.path.isfile(i + '.dph'):
                                os.remove(i + '.dph')
        os.remove(tab)
        os.remove(tab_deramp)

    print('\nAll done, enjoy it!\n')


if __name__ == "__main__":
    main()

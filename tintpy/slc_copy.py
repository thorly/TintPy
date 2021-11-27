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


def cmd_line_parser():
    parser = argparse.ArgumentParser(description='Copy a common segment from an existing set of SLCs', formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)

    parser.add_argument('slc_dir', help='SLCs directory')
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('roff', help='offset to starting rangle sample', type=int)
    parser.add_argument('nr', help='number of range samples', type=int)
    parser.add_argument('loff', help='offset to starting line', type=int)
    parser.add_argument('nl', help='number of lines', type=int)
    parser.add_argument('--rlks', help='range looks', type=int, default=8)
    parser.add_argument('--alks', help='azimuth looks', type=int, default=2)
    parser.add_argument('--num', help='number of slc used (default: -1, negative number for all)', type=int, default=-1)
    parser.add_argument('--extension', help='file extension for SLCs (defaults: .rslc)', default='.rslc')

    inps = parser.parse_args()

    return inps


EXAMPLE = """Example:
   python3 slc_copy.py ./rslc ./rslc_cut 1 1000 1 1000
   python3 slc_copy.py ./rslc ./rslc_cut 1 1000 1 1000 --rlks 8 --rlks 2 --num 1
   python3 slc_copy.py ./rslc ./rslc_cut 1 1000 1 1000 --rlks 20 --alks 5 --extension slc
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


def slc_copy(slc, slc_par, roff, nr, loff, nl, out_slc, out_slc_par):
    """Copy a common segment from an existing set of SLCs

    Args:
        slc (str): slc
        slc_par (str): slc par
        roff (int): offset to starting rangle sample
        nr (int): number of range samples
        loff (int): offset to starting line
        nl (int): number of lines
        out_slc (str): output slc
        out_slc_par (str): output slc par
    """
    call_str = f"SLC_copy {slc} {slc_par} {out_slc} {out_slc_par} - - {roff} {nr} {loff} {nl} - -"
    os.system(call_str)


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
    # get inputs
    inps = cmd_line_parser()
    slc_dir = os.path.abspath(inps.slc_dir)
    out_dir = os.path.abspath(inps.out_dir)
    roff = inps.roff
    nr = inps.nr
    loff = inps.loff
    nl = inps.nl
    rlks = inps.rlks
    alks = inps.alks
    num = inps.num
    extension = inps.extension

    # check slc_dir and out_dir
    if not os.path.isdir(slc_dir):
        sys.exit('{} does not exist.'.format(slc_dir))

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # get dates
    dates = sorted(
        [i for i in os.listdir(slc_dir) if re.findall(r'^\d{8}$', i)])

    # check dates
    if len(dates) < 1:
        sys.exit('No slc in {}'.format(slc_dir))

    # check num
    if num < 0 or num > len(dates):
        length = len(dates)
    else:
        length = num

    # check extension
    if not extension.startswith('.'):
        extension = '.' + extension

    for date in dates[0:length]:
        slc = os.path.join(slc_dir, date, date + extension)
        slc_par = slc + '.par'

        out_slc_dir = os.path.join(out_dir, date)
        out_slc = os.path.join(out_slc_dir, date + extension)
        out_slc_par = out_slc + '.par'

        if os.path.isfile(slc) and os.path.isfile(slc_par):
            if not os.path.isdir(out_slc_dir):
                os.mkdir(out_slc_dir)

            slc_copy(slc, slc_par, roff, nr, loff, nl, out_slc, out_slc_par)

            bmp = out_slc + '.bmp'
            slc2bmp(slc, slc_par, rlks, alks, bmp)

    print('\nAll done, enjoy it.\n')


if __name__ == "__main__":
    main()

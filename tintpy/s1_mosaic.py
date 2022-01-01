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

EXAMPLE = '''Example:
  python3 s1_mosaic.py /ly/slc /ly/slc -s 1 -r 20 -a 5
  python3 s1_mosaic.py /ly/slc /ly/slc -s 1 2 3 -r 20 -a 5 -p vv -n 1
'''


def cmdLineParse():
    parser = argparse.ArgumentParser(description= 'Calculate SLC mosaic of Sentinel-1 TOPS burst SLC data.',
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)

    parser.add_argument('slc_dir', help='slc directory')
    parser.add_argument('out_slc_dir', help='directory saving processed slc')
    parser.add_argument('-s', dest='sub_swath', help='sub_swath number', type=int, nargs='+', choices=[1, 2, 3], required=True)
    parser.add_argument('-r', dest='rlks', help='range looks', type=int, required=True)
    parser.add_argument('-a', dest='alks', help='azimuth looks', type=int, required=True)
    parser.add_argument('-p', dest='pol', help='polarization(defaults: vv)', choices=['vv', 'vh'], default='vv')
    parser.add_argument('-n', dest='num', help='number of slc used (default: -1, negative number for all)', type=int, default=-1)
    parser.add_argument('-e', dest='slc_extension', type=str, default='slc', help='file extension for SLCs (defaults: slc)')

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


def slc_mosaic(date_slc_dir, sub_swaths, pol, rlks, alks, extension, out_dir):
    """Calculate SLC mosaic of Sentinel-1 TOPS burst SLC data

    Args:
        date_slc_dir (str): slc directory
        sub_swaths (list): sub_swath number
        pol (str): polarization
        rlks (int): range looks
        alks (int): azimuth looks
        extension (int): slc extension
        out_dir (str): output directory

    Returns:
        tuple: mosaiced slc and par
    """
    date = os.path.basename(date_slc_dir)

    slc_tab = os.path.join(date_slc_dir, 'slc_tab')
    with open(slc_tab, 'w+', encoding='utf-8') as f:
        for i in sub_swaths:
            slc = os.path.join(date_slc_dir, f'{date}.iw{i}.{pol}.{extension}')
            slc_par = slc + '.par'
            tops_par = slc + '.tops_par'
            f.write(f'{slc} {slc_par} {tops_par}\n')

    slc_out = os.path.join(out_dir, date + '.' + extension)
    slc_par_out = slc_out + '.par'

    call_str = f"SLC_mosaic_S1_TOPS {slc_tab} {slc_out} {slc_par_out} {rlks} {alks}"
    os.system(call_str)

    os.remove(slc_tab)

    return slc_out, slc_par_out


def main():
    # get inputs
    inps = cmdLineParse()
    slc_dir = os.path.abspath(inps.slc_dir)
    out_slc_dir = os.path.abspath(inps.out_slc_dir)
    sub_swath = inps.sub_swath
    pol = inps.pol
    rlks = inps.rlks
    alks = inps.alks
    extension = inps.slc_extension
    slc_num = inps.num

    # check inputs
    if not os.path.isdir(slc_dir):
        sys.exit('{} does not exist.'.format(slc_dir))

    if not os.path.isdir(out_slc_dir):
        os.mkdir(out_slc_dir)

    if len(sub_swath) == 2:
        if sorted(sub_swath) == [1, 2] or sorted(sub_swath) == [2, 3]:
            pass
        else:
            sys.exit("Error sub_swath for two sub_swaths, it must be 1 2 or 2 3.")
    
    if extension.startswith('.'):
        extension = extension[1:]

    # get all dates
    dates = sorted([i for i in os.listdir(slc_dir) if re.findall(r'^\d{8}$', i)])

    if dates:
        if slc_num < 0 or slc_num > len(dates):
            length = len(dates)
        else:
            length = slc_num

        for date in dates[0:length]:
            date_slc_dir = os.path.join(slc_dir, date)
            out_date_slc_dir = os.path.join(out_slc_dir, date)

            if not os.path.isdir(out_date_slc_dir):
                os.mkdir(out_date_slc_dir)

            slc, slc_par = slc_mosaic(date_slc_dir, sub_swath, pol, rlks, alks, extension, out_date_slc_dir)
            slc2bmp(slc, slc_par, rlks, alks, slc + '.bmp')

        print('\nAll done, enjoy it!\n')
    else:
        print('\nCannot find any data in {}.\n'.format(slc_dir))


if __name__ == "__main__":
    main()

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
  # one sub_swath (sub_swath1 start_burst: 1 end_burst: 3)
  python3 s1_copy_mosaic.py /ly/slc /ly/slc_extract -s 1 -b 1 3 -r 20 -a 5

  # multi sub_swaths (sub_swath1 start_burst: 1 end_burst: 3, sub_swath2 start_burst: 2 end_burst: 4)
  python3 s1_copy_mosaic.py /ly/slc /ly/slc_extract -s 1 2 -b 1 3 2 4 -r 20 -a 5 -p vv
'''


def cmdLineParse():
    parser = argparse.ArgumentParser(description= 'Copy required bursts and mosaic them for Sentinel-1 TOPS SLC.',
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)

    parser.add_argument('slc_dir', help='slc directory')
    parser.add_argument('out_slc_dir', help='directory saving processed slc')
    parser.add_argument('-s', dest='sub_swath', help='sub_swath number', type=int, nargs='+', choices=[1, 2, 3], required=True)
    parser.add_argument('-b', dest='burst_num', help='burst number', type=int, nargs='+', required=True)
    parser.add_argument('-r', dest='rlks', help='range looks', type=int, required=True)
    parser.add_argument('-a', dest='alks', help='azimuth looks', type=int, required=True)
    parser.add_argument('-p', dest='pol', help='polarization(defaults: vv)', choices=['vv', 'vh'], default='vv')
    parser.add_argument('-n', dest='num', help='number of slc used (default: -1, negative number for all)', type=int, default=-1)

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


def write_tab(slc, slc_par, tops_par, tab_file):
    """Write path of slc slc_par and tops_par to tab file

    Args:
        slc (str): slc file
        slc_par (str): slc par file
        tops_par (str): tops file
        tab_file (str): output file
    """
    with open(tab_file, 'w+', encoding='utf-8') as f:
        f.write(f'{slc} {slc_par} {tops_par}\n')


def merge_tab(tab_files, out_tab):
    """Merge multi tabs to one

    Args:
        tab_files (list): tab files
        out_tab (str): out tab
    """
    with open(out_tab, 'w+', encoding='utf-8') as f_out:
        for file in tab_files:
            with open(file, 'r', encoding='utf-8') as f_in:
                f_out.write(f_in.readline().strip() + '\n')


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


def slc_copy(date_slc_dir, sub_swath, pol, start_burst, end_burst, rlks, alks, out_dir):
    """Copy multiple bursts from a Sentinel-1 TOPS mode SLC

    Args:
        date_slc_dir (str): slc directory
        sub_swath (int): sub_swath number
        pol (str): polarization
        start_burst (int): start burst number
        end_burst (int): last burst number
        rlks (int): range looks
        alks (int): azimuth looks
        out_dir (str): output slc directory

    Returns:
        str: slc_tab_out
    """
    date = os.path.basename(date_slc_dir)
    slc = os.path.join(date_slc_dir, f'{date}.iw{sub_swath}.{pol}.slc')
    slc_par = slc + '.par'
    tops_par = slc + '.tops_par'
    # output files
    slc_out = os.path.join(out_dir, os.path.basename(slc))
    slc_par_out = slc_out + '.par'
    tops_par_out = slc_out + '.tops_par'

    slc_tab_in = os.path.join(out_dir, f'{date}.iw{sub_swath}.{pol}.slc_tab_in')
    write_tab(slc, slc_par, tops_par, slc_tab_in)

    slc_tab_out = os.path.join(out_dir, f'{date}.iw{sub_swath}.{pol}.slc_tab_out')
    write_tab(slc_out, slc_par_out, tops_par_out, slc_tab_out)

    call_str = f"SLC_copy_S1_TOPS {slc_tab_in} {slc_tab_out} 1 {start_burst} 1 {end_burst}"
    os.system(call_str)

    os.remove(slc_tab_in)

    bmp = slc_out + '.bmp'
    slc2bmp(slc_out, slc_par_out, rlks, alks, bmp)

    return slc_tab_out


def slc_mosaic(slc_tab, rlks, alks, out_dir):
    """Calculate SLC mosaic of Sentinel-1 TOPS burst SLC data

    Args:
        slc_tab (str): slc tab file
        rlks (int): range looks
        alks (int): azimuth looks
        out_dir (str): output directory
    """
    # get date
    with open(slc_tab, 'r') as f:
        tmp = f.readline().strip().split()[0]
        date = os.path.basename(tmp)[0:8]

    slc_out = os.path.join(out_dir, date + '.slc')
    slc_par_out = os.path.join(out_dir, date + '.slc.par')

    call_str = f"SLC_mosaic_S1_TOPS {slc_tab} {slc_out} {slc_par_out} {rlks} {alks}"
    os.system(call_str)

    os.remove(slc_tab)

    bmp = slc_out + '.bmp'
    slc2bmp(slc_out, slc_par_out, rlks, alks, bmp)


def copy_mosaic_one(date_slc_dir, sub_swath, pol, bursts, rlks, alks, out_dir):
    """Copy and mosaic one sub_swath S1 TOPS SLC

    Args:
        date_slc_dir (str): slc directory
        sub_swath (int list): sub_swath number
        pol (str): polarization
        bursts (int list): burst number
        rlks (int): range looks
        alks (int): azimuth looks
        out_dir (str): output directory
    """
    slc_tab = slc_copy(date_slc_dir, sub_swath[0], pol, bursts[0], bursts[1], rlks, alks, out_dir)
    slc_mosaic(slc_tab, rlks, alks, out_dir)


def copy_mosaic_two(date_slc_dir, sub_swath, pol, bursts, rlks, alks, out_dir):
    """Copy and mosaic two sub_swaths S1 TOPS SLC

    Args:
        date_slc_dir (str): slc directory
        sub_swath (int list): sub_swath number
        pol (str): polarization
        bursts (int list): burst number
        rlks (int): range looks
        alks (int): azimuth looks
        out_dir (str): output directory
    """
    slc_tab1 = slc_copy(date_slc_dir, sub_swath[0], pol, bursts[0], bursts[1], rlks, alks, out_dir)
    slc_tab2 = slc_copy(date_slc_dir, sub_swath[1], pol, bursts[2], bursts[3], rlks, alks, out_dir)

    slc_tab = os.path.join(out_dir, 'slc_tab')
    merge_tab([slc_tab1, slc_tab2], slc_tab)
    os.remove(slc_tab1)
    os.remove(slc_tab2)

    slc_mosaic(slc_tab, rlks, alks, out_dir)


def copy_mosaic_three(date_slc_dir, sub_swath, pol, bursts, rlks, alks,
                      out_dir):
    """Copy and mosaic three sub_swaths S1 TOPS SLC

    Args:
        date_slc_dir (str): slc directory
        sub_swath (int list): sub_swath number
        pol (str): polarization
        bursts (int list): burst number
        rlks (int): range looks
        alks (int): azimuth looks
        out_dir (str): output directory
    """
    slc_tab1 = slc_copy(date_slc_dir, sub_swath[0], pol, bursts[0], bursts[1], rlks, alks, out_dir)
    slc_tab2 = slc_copy(date_slc_dir, sub_swath[1], pol, bursts[2], bursts[3], rlks, alks, out_dir)
    slc_tab3 = slc_copy(date_slc_dir, sub_swath[2], pol, bursts[4], bursts[5], rlks, alks, out_dir)

    slc_tab = os.path.join(out_dir, 'slc_tab')
    merge_tab([slc_tab1, slc_tab2, slc_tab3], slc_tab)
    os.remove(slc_tab1)
    os.remove(slc_tab2)
    os.remove(slc_tab3)

    slc_mosaic(slc_tab, rlks, alks, out_dir)


def main():
    # get inputs
    inps = cmdLineParse()
    slc_dir = os.path.abspath(inps.slc_dir)
    out_slc_dir = os.path.abspath(inps.out_slc_dir)
    sub_swath = inps.sub_swath
    burst_num = inps.burst_num
    pol = inps.pol
    rlks = inps.rlks
    alks = inps.alks
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
    if len(sub_swath) == 3:
        if sorted(sub_swath) != [1, 2, 3]:
            sys.exit("Error sub_swath for three sub_swaths, must be 1 2 3.")

    if 2 * len(sub_swath) != len(burst_num):
        sys.exit('Error sub_swath or butst_num.')

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

            if len(sub_swath) == 1:
                exec_func = copy_mosaic_one
            elif len(sub_swath) == 2:
                exec_func = copy_mosaic_two
            else:
                exec_func = copy_mosaic_three

            exec_func(date_slc_dir, sub_swath, pol, burst_num, rlks, alks,out_date_slc_dir)

        print('\nAll done, enjoy it!\n')
    else:
        print('\nCannot find any data in {}.\n'.format(slc_dir))


if __name__ == "__main__":
    main()

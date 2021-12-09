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
  python3 s1_copy_bursts_auto.py slc slc_extract ref.number_table -r 32 -a 8 -p vv -n 1
'''


def cmdLineParse():
    parser = argparse.ArgumentParser(description= 'Copy multiple bursts automatically from Sentinel-1 TOPS mode SLC..',
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)

    parser.add_argument('slc_dir', help='slc directory')
    parser.add_argument('out_slc_dir', help='directory saving processed slc')
    parser.add_argument('ref_table', help='reference burst number table')
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
        tuple: slc and par
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
    os.remove(slc_tab_out)

    return slc_out, slc_par_out


def parse_table(table):
    """Parse burst_number_table file

    Args:
        table (str): table file

    Returns:
        dict: iw_number_time
    """
    iw_number_time = {}
    with open(table, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        for line in lines:
            for i in [1, 2, 3]:
                if f'iw{i}_number_of_bursts' in line:
                    index = lines.index(line)
                    num = line.strip().split()[1]
                    first = lines[index + 1].strip().split()[-1]
                    last = lines[index + 2].strip().split()[-1]
                    iw_number_time[f'iw{i}'] = num + ' ' + first + ' ' + last

    return iw_number_time


def get_iw_burst(ref_table, table):
    """Get iw and burst

    Args:
        ref_table (str): reference table
        table (str): table

    Returns:
        dict: iw_burst
    """
    iw_burst = {}
    iw_number_time1 = parse_table(ref_table)
    iw_number_time2 = parse_table(table)

    for i in iw_number_time1.keys():
        time1 = iw_number_time1.get(i)
        time2 = iw_number_time2.get(i, None)
        if time2:
            tmp1 = time1.split()
            num1, first1, _ = int(tmp1[0]), float(tmp1[1]), float(tmp1[2])
            tmp2 = time2.split()
            num2, first2, _ = int(tmp2[0]), float(tmp2[1]), float(tmp2[2])
            for j in range(num2):
                if abs(first1 - first2 - j) < 0.1:
                    if j + num1 <= num2:
                        iw_burst[i] = [j+1, j+num1]

    return iw_burst


def main():
    # get inputs
    inps = cmdLineParse()
    slc_dir = os.path.abspath(inps.slc_dir)
    out_slc_dir = os.path.abspath(inps.out_slc_dir)
    ref_table = os.path.abspath(inps.ref_table)
    pol = inps.pol
    rlks = inps.rlks
    alks = inps.alks
    slc_num = inps.num

    # check inputs
    if not os.path.isdir(slc_dir):
        sys.exit('{} does not exist.'.format(slc_dir))

    if not os.path.isdir(out_slc_dir):
        os.mkdir(out_slc_dir)

    if not os.path.isfile(ref_table):
        sys.exit('{} does not exist.'.format(ref_table))

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

            table = os.path.join(date_slc_dir, f'{date}.burst_number_table')
            iw_burst = get_iw_burst(ref_table, table)

            sub_swath = [int(i[2:]) for i in iw_burst.keys()]
            burst_num = []
            for i in iw_burst.values():
                burst_num += i

            for i in range(len(sub_swath)):
                slc, slc_par = slc_copy(date_slc_dir, sub_swath[i], pol, burst_num[2*i], burst_num[2*i+1], rlks, alks, out_date_slc_dir)
                slc2bmp(slc, slc_par, rlks, alks, slc + '.bmp')

        print('\nAll done, enjoy it!\n')
    else:
        print('\nCannot find any data in {}.\n'.format(slc_dir))


if __name__ == "__main__":
    main()

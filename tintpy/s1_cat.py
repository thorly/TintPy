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
import shutil
import sys

EXAMPLE = """Example:
  [Note: This script only concatenates adjacent SLC processed by s1_zip2slc.py]
  python3 s1_cat.py slc slc_cat 1 --rlks 8 --rlks 2
  python3 s1_cat.py slc slc_cat 1 2 --rlks 8 --rlks 2
  python3 s1_cat.py slc slc_cat 1 2 3 --rlks 20 --alks 5 --pol vv
"""


def cmdline_parse():
    parser = argparse.ArgumentParser(description='Concatenate adjacent Sentinel-1 TOPS SLC images(support two or three images)',
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)
    parser.add_argument('slc_dir', help='slc directory')
    parser.add_argument('save_dir', help='directory saving concatenated slc')
    parser.add_argument('sub_swath', help='sub_swath number', type=int, nargs='+', choices=[1, 2, 3])
    parser.add_argument('--rlks', help='range looks', type=int, required=True)
    parser.add_argument('--alks', help='azimuth looks', type=int, required=True)
    parser.add_argument('--pol', help='polarization(defaults: vv)', nargs='+', choices=['vv', 'vh'], default=['vv'])
    inps = parser.parse_args()

    return inps


def get_state_vector(par_file):
    """Get all state_vectors from par file

    Args:
        par_file (str): par file

    Returns:
        str list: state vector
    """
    state_vector = []
    temp_list = []
    with open(par_file) as f:
        for line in f.readlines():
            if line.startswith('state_vector_position_') or line.startswith('state_vector_velocity_'):
                temp_list.append(line.strip())
    for s in temp_list:
        line = re.split(r'\s+', s)
        state_vector.append(line[1:4])

    return state_vector


def get_nonrepeated_state_vector(state_vector_1, state_vector_2):
    """Compare two state_vectors, return non-repeated state vector and first state vector in state_vector_2

    Args:
        state_vector_1 (str list): state vector1
        state_vector_2 (str list): state vector2

    Returns:
        str list: non-repeated state vector
    """
    nonrepeated_state_vector = []
    index = 1
    for sv in state_vector_2:
        if sv not in state_vector_1:
            nonrepeated_state_vector.append(sv)
    for i in range(len(state_vector_2)):
        if state_vector_2[i] == nonrepeated_state_vector[0]:
            index = i

    return nonrepeated_state_vector, index


def get_time_of_first_state_vector(par_file):
    """get time_of_first_state_vector from .par file

    Args:
        par_file (str): par file

    Returns:
        str: time of first state vector
    """
    time_of_first_state_vector = ''
    with open(par_file, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            if line.startswith('time_of_first_state_vector'):
                time_of_first_state_vector = re.split(r'\s+', line)[1]

    return time_of_first_state_vector


def get_content_brfore(par_file):
    """Get content before the line of number_of_state_vectors

    Args:
        par_file (str): par file

    Returns:
        str: content before the line of number_of_state_vectors
    """
    content = ''
    with open(par_file, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            if not line.startswith('number_of_state_vectors'):
                content += line
            else:
                break

    return content


def get_new_par(par_file1, par_file2, new_par_file2):
    """Compare two par file, write non-repeated .par file

    Args:
        par_file1 (str): par file1
        par_file2 (str): par file2
        out_par_file (str): output new par file
    """
    state_vector_1 = get_state_vector(par_file1)
    state_vector_2 = get_state_vector(par_file2)
    nonrepeated_state_vector, index = get_nonrepeated_state_vector(state_vector_1, state_vector_2)
    number_of_state_vectors = int(len(nonrepeated_state_vector) / 2)

    # calculate the new time_of_first_state_vector
    time_of_first_state_vector = get_time_of_first_state_vector(par_file2)
    time_of_first_state_vector = float(time_of_first_state_vector) + index / 2 * 10
    time_of_first_state_vector = str(time_of_first_state_vector) + '00000'

    # get content before 'number_of_state_vectors'
    content = get_content_brfore(par_file2)

    # write new par file
    with open(new_par_file2, 'w+', encoding='utf-8') as f:
        f.write(content)
        f.write('number_of_state_vectors:' + str(number_of_state_vectors).rjust(21, ' ') + '\n')
        f.write('time_of_first_state_vector:' + time_of_first_state_vector.rjust(18, ' ') + '   s' + '\n')
        f.write('state_vector_interval:              10.000000   s' + '\n')
        nr = nonrepeated_state_vector
        for i in range(number_of_state_vectors):
            f.write('state_vector_position_' + str(i + 1) + ':' + nr[i * 2][0].rjust(15, ' ') + nr[i * 2][1].rjust(16, ' ') + nr[i * 2][2].rjust(16, ' ') + '   m   m   m' + '\n')
            f.write('state_vector_velocity_' + str(i + 1) + ':' + nr[i * 2 + 1][0].rjust(15, ' ') + nr[i * 2 + 1][1].rjust(16, ' ') + nr[i * 2 + 1][2].rjust(16, ' ') + '   m/s m/s m/s' + '\n')


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


def get_time_and_direction(par_file):
    """get time and orbit direction form .par file

    Args:
        par_file (str): par file

    Returns:
        tuple: time and direction
    """
    with open(par_file, 'r') as f:
        for l in f.readlines():
            if 'start_time' in l:
                start_time = l.strip().split()[1]
            if 'heading' in l:
                heading = l.strip().split()[1]
    heading = float(heading)
    start_time = float(start_time)
    if heading > -180 and heading < -90:
        direction = 'DES'
    else:
        direction = 'ASC'

    return start_time, direction


def get_date_for_cat(slc_dir, date_num):
    """get the date which has date_num images

    Args:
        slc_dir (str): slc directory
        date_num (int): date number

    Returns:
        list: date
    """
    dates = []
    for i in os.listdir(slc_dir):
        if re.search(r'^\d{8}-\d{1}$', i):
            dates.append(i[0:8])
    dates = set(j for j in dates if dates.count(j) == date_num)

    return sorted(list(dates))


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


def cat_two_slc(date_slc_dirs, sub_swath, pol, rlks, alks, out_dir):
    """Concatenate two SLCs

    Args:
        date_slc_dirs (str list): slc directorys
        sub_swath (int list): sub_swath number
        pol (str list): polariaztion
        rlks (int): range looks
        alks (int): azimuth looks
        out_dir (str): output directory
    """
    date = os.path.basename(date_slc_dirs[0])[0:8]
    for p in pol:
        for s in sub_swath:
            slc1 = os.path.join(date_slc_dirs[0], f'{date}.iw{s}.{p}.slc')
            slc_par1 = slc1 + '.par'
            tops_par1 = slc1 + '.tops_par'

            slc2 = os.path.join(date_slc_dirs[1], f'{date}.iw{s}.{p}.slc')
            slc_par2 = slc2 + '.par'
            tops_par2 = slc2 + '.tops_par'

            slc_out = os.path.join(out_dir, f'{date}.iw{s}.{p}.slc')
            slc_par_out = slc_out + '.par'
            tops_par_out = slc_out + '.tops_par'

            # backup par file
            par1_backup = slc_par1 + '.backup'
            par2_backup = slc_par2 + '.backup'
            if not os.path.isfile(par1_backup):
                shutil.copy(slc_par1, par1_backup)
            if not os.path.isfile(par2_backup):
                shutil.copy(slc_par2, par2_backup)

            # get image start time and direction
            start_time1, direction1 = get_time_and_direction(slc_par1)
            start_time2, direction2 = get_time_and_direction(slc_par2)

            # exchange values for the following cases
            des = (direction1 == 'DES') and (start_time1 > start_time2)
            asc = (direction1 == 'ASC') and (start_time1 > start_time2)
            if asc or des:
                slc1, slc2 = slc2, slc1
                slc_par1, slc_par2 = slc_par2, slc_par1
                tops_par1, tops_par2 = tops_par2, tops_par1

            # delete repeated state vector and write new par file
            get_new_par(slc_par1, slc_par2, slc_par2)

            slc_tab_in1 = os.path.join(out_dir, 'slc_tab_in1')
            slc_tab_in2 = os.path.join(out_dir, 'slc_tab_in2')
            slc_tab_out = os.path.join(out_dir, 'slc_tab_out')
            write_tab(slc1, slc_par1, tops_par1, slc_tab_in1)
            write_tab(slc2, slc_par2, tops_par2, slc_tab_in2)
            write_tab(slc_out, slc_par_out, tops_par_out, slc_tab_out)

            # concatenate adjacent Sentinel-1 TOPS SLC images
            call_str = f'SLC_cat_S1_TOPS.TOPS {slc_tab_in1} {slc_tab_in2} {slc_tab_out}'
            os.system(call_str)

            bmp = slc_out + '.bmp'
            slc2bmp(slc_out, slc_par_out, rlks, alks, bmp)

            os.remove(slc_tab_in1)
            os.remove(slc_tab_in2)
            os.remove(slc_tab_out)


def cat_three_slc(date_slc_dirs, sub_swath, pol, rlks, alks, out_dir):
    """Concatenate three SLCs

    Args:
        date_slc_dirs (str list): slc directorys
        sub_swath (int list): sub_swath number
        pol (str list): polariaztion
        rlks (int): range looks
        alks (int): azimuth looks
        out_dir (str): output directory
    """
    date = os.path.basename(date_slc_dirs[0])[0:8]

    # get start time
    slc_par1 = os.path.join(date_slc_dirs[0], f'{date}.iw{sub_swath[0]}.{pol[0]}.slc.par')
    slc_par2 = os.path.join(date_slc_dirs[1], f'{date}.iw{sub_swath[0]}.{pol[0]}.slc.par')
    slc_par3 = os.path.join(date_slc_dirs[2], f'{date}.iw{sub_swath[0]}.{pol[0]}.slc.par')
    start_time1, _ = get_time_and_direction(slc_par1)
    start_time2, _ = get_time_and_direction(slc_par2)
    start_time3, _ = get_time_and_direction(slc_par3)

    middle_time = sorted([start_time1, start_time2, start_time3])[1]

    # create tmp_dir
    tmp_dir = date_slc_dirs[0] + '_tmp'
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    # first case: slc1 or slc2 is middle slc
    if middle_time in [start_time1, start_time2]:
        cat_two_slc(date_slc_dirs[0], date_slc_dirs[1], sub_swath, pol, rlks, alks, out_dir)
        cat_two_slc(date_slc_dirs[2], tmp_dir, sub_swath, pol, rlks, alks, out_dir)
    # second case: slc3 is middle slc
    else:
        cat_two_slc(date_slc_dirs[0], date_slc_dirs[2], sub_swath, pol, rlks, alks, out_dir)
        cat_two_slc(date_slc_dirs[1], tmp_dir, sub_swath, pol, rlks, alks, out_dir)

    # delete tmp_dir
    shutil.rmtree(tmp_dir)


def main():
    # get inputs
    inps = cmdline_parse()
    slc_dir = os.path.abspath(inps.slc_dir)
    save_dir = os.path.abspath(inps.save_dir)
    sub_swath = inps.sub_swath
    pol = inps.pol
    rlks = inps.rlks
    alks = inps.alks

    # check inputs
    if not os.path.isdir(slc_dir):
        sys.exit('{} does not exist.'.format(slc_dir))

    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    # cat slc
    for num in [2, 3]:
        dates = get_date_for_cat(slc_dir, num)

        for date in dates:
            cat_dir = os.path.join(save_dir, date)
            if not os.path.isdir(cat_dir):
                os.mkdir(cat_dir)

            slc_dirs = glob.glob(os.path.join(slc_dir, date + '-[0-9]'))

            if num == 2:
                cat_func = cat_two_slc
            if num == 3:
                cat_func = cat_three_slc

            cat_func(slc_dirs, sub_swath, pol, rlks, alks, cat_dir)

    print('\nAll done, enjoy it!\n')


if __name__ == "__main__":
    main()

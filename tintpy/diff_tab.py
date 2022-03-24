#!/usr/bin/env python3
################################
# Program is part of TintPy    #
# Copyright (c) 2021, Lei Yuan #
# Author: Lei Yuan, 2021       #
################################

import argparse
import datetime
import os
import re
import sys

import numpy as np

EXAMPLE = """Example:
  python3 diff_tab.py /ly/diff 0.5 diff_tab
  python3 diff_tab.py /ly/diff 0.5 diff_tab -u adf.unw.gacos -c adf.cc.gacos -d diff.par
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(description='Write diff_tab for Phase-Stacking using GAMMA.',
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)

    parser.add_argument('diff_dir', help='diff directory including *.unw and *.cc.')
    parser.add_argument('geo_dir', help='geo directory including *.diff_par.')
    parser.add_argument('coh',help='coherence threshold (smaller than this will not be used) for Phase-Stacking.', type=float)
    parser.add_argument('out_tab', help='diff_tab file for Phase-Stacking.')
    parser.add_argument('-u', dest='unw_extension', help='filename extension of unwrapped file (defaults: adf.unw).', default='adf.unw')
    parser.add_argument('-c', dest='cc_extension', help='filename extension of coherence file (defaults: adf.cc).', default='adf.cc')
    parser.add_argument('-d', dest='diff_par_extension', help='filename extension of diff_par file (defaults: diff_par).', default='diff_par')
    inps = parser.parse_args()

    return inps


def calc_time_delta(pair):
    """Calculate time delta

    Args:
        pair (str): date1-date2

    Returns:
        int: days
    """
    dates = re.findall(r'\d{8}', pair)
    date1 = datetime.datetime.strptime(dates[0], '%Y%m%d')
    date2 = datetime.datetime.strptime(dates[1], '%Y%m%d')
    time_delta = (date2 - date1).days

    return time_delta


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
                value = line.split()[1].strip()

    return value


def read_gamma(file, lines, file_type):
    """Read GAMMA format file

    Args:
        file (str): GAMMA format file
        lines (int): line of file
        file_type (str): file type

    Returns:
        array: data
    """
    # check file
    if not os.path.isfile(file):
        sys.exit('{} does not exist.'.format(file))
    data = np.fromfile(file, dtype=file_type)
    # GAMMA output files are big-endian
    data.byteswap('True')
    data = data.reshape(lines, -1)

    return data


def check_extension(extension):
    """Check file extension

    Args:
        extension (str): extension

    Returns:
        str: extension
    """
    if not extension.startswith('.'):
        extension = '.' + extension

    return extension


def get_pairs(diff_dir):
    """Get all pairs

    Args:
        diff_dir (str): directory cluding pairs

    Returns:
        list: pairs
    """
    pairs = [i[0:17] for i in os.listdir(diff_dir) if re.match(r'\d{8}_\d{8}', i)]
    pairs = sorted(list(set(pairs)))

    return pairs


def get_diff_par(diff_dir, geo_dir, diff_par_extension):
    """Get diff_par file

    Args:
        diff_dir (str): directory contains unw cc
        geo_dir (str): directory contains diff_par
        diff_par_extension (str): diff_par file extension

    Returns:
        list: diff_par file path
    """
    pairs = get_pairs(diff_dir)
    diff_pars = []

    for pair in pairs:
        diff_par = os.path.join(geo_dir, pair[0:8] + '.diff_par')
        diff_pars.append(diff_par)

    return diff_pars


def get_unw_cc(diff_dir, unw_extension, cc_extension):
    """Get unw and cc file

    Args:
        diff_dir (str): directory contains unw and cc
        unw_extension (str): unw extension
        cc_extension (cc): cc extension

    Returns:
        tuple: unws and ccs
    """
    pairs = get_pairs(diff_dir)

    unws = []
    ccs = []

    for pair in pairs:
        unw = os.path.join(diff_dir, pair + unw_extension)
        cc = os.path.join(diff_dir, pair + cc_extension)
        if os.path.isfile(unw) and os.path.isfile(cc):
            unws.append(unw)
            ccs.append(cc)

    return unws, ccs


def get_mean_coh(ccs, diff_pars):
    """Calculate mean coherence

    Args:
        ccs (list): cc files
        diff_pars (list): diff_par file

    Returns:
        array: mean coherence array
    """
    mean_coh_array = np.zeros(len(ccs))

    for cc, diff_par in zip(ccs, diff_pars):
        length = int(read_gamma_par(diff_par, 'az_samp_1'))
        mean_coh = np.mean(read_gamma(cc, length, 'float32'))
        mean_coh_array[ccs.index(cc)] = mean_coh

    np.save('mean_coh_array', mean_coh_array)

    return mean_coh_array


def count_coh(start, end, coh_array):
    """Count coherence rate in [start end]

    Args:
        start (float): start point
        end (float): end point
        coh_array (array): mean coherence array
    """
    num = 0
    for i in coh_array:
        if i > start and i <= end:
            num += 1
    rate = num / coh_array.shape[0]
    print(f"{start}~{end}: {round(rate, 2)}")


def statistic_coh(mean_coh_array, step=0.05):
    """Statistic coherence

    Args:
        mean_coh_array (array): mean coherence array
        step (float, optional): step. Defaults to 0.05.
    """
    for i in np.arange(0, 1, step):
        start = round(i, 2)
        end = round(start + step, 2)
        count_coh(start, end, mean_coh_array)


def write_diff_tab(unws, mean_coh_array, coh_thres, tab_file):
    """Write diff_tab file

    Args:
        unws (list): unw files
        mean_coh_array (array): mean coherence array
        coh_thres (float): coherence threshold
        tab_file (str): output tab file
    """
    print('\nWriting data to {}'.format(tab_file))
    unw_used = 0
    with open(tab_file, 'w+') as f:
        for unw in unws:
            unw_name = os.path.basename(unw)
            time_delta = calc_time_delta(unw_name[0:17])
            coh = mean_coh_array[unws.index(unw)]
            if coh >= coh_thres:
                f.write(f"{unw}    {time_delta}\n")
                unw_used += 1

    rate = round(unw_used / len(unws), 2) * 100
    print('unwrapped file used: {}/{} [{}%]'.format(unw_used, len(unws), rate))


def main():
    # get inputs
    inps = cmdline_parser()
    diff_dir = os.path.abspath(inps.diff_dir)
    geo_dir = os.path.abspath(inps.geo_dir)
    coh_thres = inps.coh
    tab_file = inps.out_tab
    unw_extension = inps.unw_extension
    cc_extension = inps.cc_extension
    diff_par_extension = inps.diff_par_extension

    # check diff_dir
    if not os.path.isdir(diff_dir):
        sys.exit('{} does not exist'.format(diff_dir))

    # check diff_dir
    if not os.path.isdir(geo_dir):
        sys.exit('{} does not exist'.format(geo_dir))

    # check extension
    unw_extension = check_extension(unw_extension)
    cc_extension = check_extension(cc_extension)
    diff_par_extension = check_extension(diff_par_extension)

    # get pairs
    pairs = get_pairs(diff_dir)

    if len(pairs) < 1:
        sys.exit('No pair in {}.'.format(diff_dir))

    # get diff_par
    diff_pars = get_diff_par(diff_dir, geo_dir, diff_par_extension)

    # get unws and ccs:
    unws, ccs = get_unw_cc(diff_dir, unw_extension, cc_extension)

    if len(unws) < 1:
        sys.exit('No unw file in {}.'.format(diff_dir))

    # statistic coherence
    print('Statisticing coherence:')
    if os.path.isfile('mean_coh_array.npy'):
        mean_coh_array = np.load('mean_coh_array.npy')
    else:
        mean_coh_array = get_mean_coh(ccs, diff_pars)

    statistic_coh(mean_coh_array)

    # write diff_tab
    write_diff_tab(unws, mean_coh_array, coh_thres, tab_file)

    print('\nAll done, enjot it!\n')


if __name__ == "__main__":
    main()

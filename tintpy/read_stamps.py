#!/usr/bin/env python3
################################
# Program is part of TintPy    #
# Copyright (c) 2021, Lei Yuan #
# Author: Lei Yuan, 2021       #
################################

import argparse
import datetime
import os
import sys

import numpy as np
import scipy.io as scio

EXAMPLE = """Example:
  python3 read_stamps.py /ly/ps_plot_v-do.mat /ly/ps_plot_ts_v-do.mat
  python3 read_stamps.py /ly/ps_plot_v-do /ly/ps_plot_ts_v-do -v vel_hy.txt -t ts_hy.txt
"""


def cmdline_parser():
    parser = argparse.ArgumentParser("Read results processed by StaMPS SBAS or PS.",
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)
    parser.add_argument('vel_mat', help="velocity mat [use ps_plot('V-do', -1) to get it].")
    parser.add_argument('ts_mat', help="time-series mat [use ps_plot('V-do', 'ts') to get it].")
    parser.add_argument('-v', dest='vel_out', default='vel.txt', help="output velocity file.")
    parser.add_argument('-t', dest='ts_out', default='ts.txt', help="output time-series file.")

    inps = parser.parse_args()

    return inps


def get_dates(days):
    """Get date using days

    Args:
        days (list): days

    Returns:
        list: date
    """
    dates = []
    for day in days:
        date = datetime.datetime.fromordinal(int(day)) + datetime.timedelta(
            days=int(day) % 1) - datetime.timedelta(days=366)
        date_str = date.strftime('%Y%m%d')
        dates.append(int(date_str))

    return dates


def save_vel_ts(ts_mat, vel_mat, out_ts_file, out_vel_file):
    """Read vel and ts mat and save them

    Args:
        ts_mat (str): ts file
        vel_mat (str): vel file
        out_ts_file (str): output ts file
        out_vel_file (str): output vel file
    """
    # load mat
    data_ts = scio.loadmat(ts_mat)
    data_vel = scio.loadmat(vel_mat)

    # get ts vel date lon lat
    vel = data_vel['ph_disp']
    lonlat = data_ts['lonlat']
    ts = data_ts['ph_mm']
    days = data_ts['day']

    # day --> date
    dates = get_dates(days)
    master_day = data_ts['master_day']
    master_date = get_dates(master_day)
    all_dates = sorted(dates + master_date)

    # get master index
    master_index = all_dates.index(master_date[0])
    master_disp = np.zeros((ts.shape[0], 1))

    # add master_disp to ts
    if master_index == 0:
        ts = np.hstack((master_index, ts))
    elif master_index == len(all_dates) - 1:
        ts = np.hstack((ts, master_index))
    else:
        ts_before = ts[:, 0:master_index].reshape((-1, master_index))
        ts_after = ts[:, master_index:].reshape((ts.shape[0], -1))
        ts = np.hstack((ts_before, master_disp, ts_after))

    ts = ts - ts[:, 0].reshape((-1, 1))

    # save vel file
    num = np.arange(0, vel.shape[0]).reshape((-1, 1))
    vel_out = np.hstack((num, lonlat, vel))
    print('Writing data to {}'.format(out_vel_file))
    np.savetxt(out_vel_file, vel_out, fmt='%4f')
    print('done.')

    # save ts file
    all_dates = np.asarray(all_dates, dtype='int64').reshape((1, -1))
    dates_out = np.hstack((np.array([-1, -1, -1, -1]).reshape((1, -1)), all_dates))
    tmp = np.hstack((vel_out, ts))
    ts_out = np.vstack((dates_out, tmp))
    print('Writing data to {}'.format(out_ts_file))
    np.savetxt(out_ts_file, ts_out, fmt='%4f')
    print('done.')

def main():
    # get inputs
    inps = cmdline_parser()
    vel_mat = os.path.abspath(inps.vel_mat)
    ts_mat = os.path.abspath(inps.ts_mat)
    vel_out = inps.vel_out
    ts_out = inps.ts_out
    # check vel_mat and ts_mat
    if not vel_mat.endswith('.mat'):
        vel_mat += '.mat'

    if not os.path.isfile(vel_mat):
        sys.exit('Cannot find {}, please check it!'.format(vel_mat))

    if not ts_mat.endswith('.mat'):
        ts_mat += '.mat'

    if not os.path.isfile(ts_mat):
        sys.exit('Cannot find {}, please check it!'.format(ts_mat))

    # save files
    save_vel_ts(ts_mat, vel_mat, ts_out, vel_out)

    print('\nAll done, enjoy it!\n')


if __name__ == "__main__":
    main()

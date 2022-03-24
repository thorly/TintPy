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


def cmd_line_parser():
    parser = argparse.ArgumentParser(
        description='Prepare files for MintPy time series processing.',
        formatter_class=argparse.RawTextHelpFormatter,epilog=EXAMPLE)
    parser.add_argument('diff_dir', help='diff directory')
    parser.add_argument('mli_dir', help='mli directory')
    parser.add_argument('geo_dir', help='geo directory')
    parser.add_argument('mintpy_dir', help='output directory')
    parser.add_argument('-m', dest='mli_par_extension', default='rmli.par', help='mli parameter file extension (defaults: rmli.par)')
    parser.add_argument('-u', dest='unw_extension', default='adf.unw', help='unw file extension (defaults: adf.unw)')
    parser.add_argument('-c', dest='cc_extension', default='adf.cc', help='cc file extension (defaults: adf.cc)')
    inps = parser.parse_args()

    return inps


EXAMPLE = """Example:
  python3 gamma2mintpy2d.py /ly/diff /ly/mli /ly/geo /ly/MintPy
  python3 gamma2mintpy2d.py /ly/diff /ly/mli /ly/geo /ly/MintPy -u adf.unw.gacos -c adf.cc.gacos
"""


def is_dir(dir):
    """Check directory

    Args:
        dir (str): directory
    """
    if not os.path.isdir(dir):
        sys.exit('{} does not exist'.format(dir))

def mk_dir(dir):
    """Make directory

    Args:
        dir (str): directory
    """
    if not os.path.isdir(dir):
        os.mkdir(dir)

def get_rlks(mli_par):
    """Get rlks from par file

    Args:
        mli_par (str): par file

    Returns:
        str: rlks
    """
    rlks = None
    with open(mli_par, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            if line.count('range_looks') == 1:
                rlks = line.split(':')[1].strip()

    return rlks

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


def main():
    # get inputs
    inps = cmd_line_parser()
    diff_dir = os.path.abspath(inps.diff_dir)
    mli_dir = os.path.abspath(inps.mli_dir)
    geo_dir = os.path.abspath(inps.geo_dir)
    mintpy_dir = os.path.abspath(inps.mintpy_dir)
    mli_par_extension = inps.mli_par_extension
    unw_extension = inps.unw_extension
    cc_extension = inps.cc_extension

    # check
    is_dir(diff_dir)
    is_dir(mli_dir)
    is_dir(geo_dir)
    mk_dir(mintpy_dir)

    mli_par_extension = check_extension(mli_par_extension)
    unw_extension = check_extension(unw_extension)
    cc_extension = check_extension(cc_extension)

    # directory for mintpy
    geometry_dir = os.path.join(mintpy_dir, 'geometry')
    interferograms_dir = os.path.join(mintpy_dir, 'interferograms')
    mk_dir(geometry_dir)
    mk_dir(interferograms_dir)

    # get pairs
    pairs = [i[0:17] for i in os.listdir(diff_dir) if re.match(r'\d{8}_\d{8}', i)]
    pairs = sorted(list(set(pairs)))

    if len(pairs) < 1:
        sys.exit('No pair in {}.'.format(diff_dir))

    # get rlks
    date = pairs[0][0:8]
    mli_par = os.path.join(mli_dir, date + mli_par_extension)
    rlks = get_rlks(mli_par)

    print('Copying files into {}'.format(geometry_dir))

    UTM_TO_RDC_path = os.path.join(geo_dir, date + '.lookup_fine')
    dst_UTM_TO_RDC = 'sim_' + date + '_' + rlks + 'rlks.UTM_TO_RDC'
    dst_UTM_TO_RDC_path = os.path.join(geometry_dir, dst_UTM_TO_RDC)

    diff_par_path = os.path.join(geo_dir, date + '.diff_par')
    dst_diff_par = 'sim_' + date + '_' + rlks + 'rlks.diff_par'
    dst_diff_par_path = os.path.join(geometry_dir, dst_diff_par)

    rdc_dem_path = os.path.join(geo_dir, date + '_dem.rdc')
    dst_rdc_dem = 'sim_' + date + '_' + rlks + 'rlks.rdc.dem'
    dst_rdc_dem_path = os.path.join(geometry_dir, dst_rdc_dem)

    utm_dem_path = os.path.join(geo_dir, date + '.dem_seg')
    dst_utm_dem = 'sim_' + date + '_' + rlks + 'rlks.utm.dem'
    dst_utm_dem_path = os.path.join(geometry_dir, dst_utm_dem)

    utm_dem_par_path = utm_dem_path + '.par'
    dst_utm_dem_par = dst_utm_dem + '.par'
    dst_utm_dem_par_path = os.path.join(geometry_dir, dst_utm_dem_par)

    shutil.copy(UTM_TO_RDC_path, dst_UTM_TO_RDC_path)
    shutil.copy(diff_par_path, dst_diff_par_path)
    shutil.copy(rdc_dem_path, dst_rdc_dem_path)
    shutil.copy(utm_dem_path, dst_utm_dem_path)
    shutil.copy(utm_dem_par_path, dst_utm_dem_par_path)

    print('Copying files into {}'.format(interferograms_dir))

    for pair in pairs:
        m_date = pair[0:8]
        s_date = pair[-8:]

        off = pair + '.off'
        off_path = os.path.join(diff_dir, off)
        e1 = os.path.isfile(off_path)

        m_amp_par = m_date + mli_par_extension
        m_amp_par_path = os.path.join(mli_dir, m_amp_par)
        e2 = os.path.isfile(m_amp_par_path)

        s_amp_par = s_date + mli_par_extension
        s_amp_par_path = os.path.join(mli_dir, s_amp_par)
        e3 = os.path.isfile(s_amp_par_path)

        cor = pair + cc_extension
        cor_path = os.path.join(diff_dir, cor)
        e3 = os.path.isfile(cor_path)

        unw = pair + unw_extension
        unw_path = os.path.join(diff_dir, unw)
        e4 = os.path.isfile(unw_path)

        if e1 and e2 and e3 and e4:
            dst_ifg_dir = os.path.join(interferograms_dir, pair)
            mk_dir(dst_ifg_dir)

            dst_off = pair + '_' + rlks + 'rlks.off'
            dst_off_path = os.path.join(dst_ifg_dir, dst_off)
            shutil.copy(off_path, dst_off_path)

            dst_m_amp_par = m_date + '_' + rlks + 'rlks.amp.par'
            dst_m_amp_par_path = os.path.join(dst_ifg_dir, dst_m_amp_par)
            shutil.copy(m_amp_par_path, dst_m_amp_par_path)

            dst_s_amp_par = s_date + '_' + rlks + 'rlks.amp.par'
            dst_s_amp_par_path = os.path.join(dst_ifg_dir, dst_s_amp_par)
            shutil.copy(s_amp_par_path, dst_s_amp_par_path)

            dst_cor = 'filt_' + pair + '_' + rlks + 'rlks.cor'
            dst_cor_path = os.path.join(dst_ifg_dir, dst_cor)
            shutil.copy(cor_path, dst_cor_path)

            dst_unw = 'diff_' + pair + '_' + rlks + 'rlks.unw'
            dst_unw_path = os.path.join(dst_ifg_dir, dst_unw)
            shutil.copy(unw_path, dst_unw_path)

    # copy notebook to mintpy_dir
    try:
        TINTPY_HOME = os.environ['TINTPY_HOME']
        notebook1 = os.path.join(TINTPY_HOME, 'tintpy', 'notebook', 'MintPy.ipynb')
        notebook2 = os.path.join(TINTPY_HOME, 'tintpy', 'notebook', 'read_mintpy.ipynb')
        shutil.copy(notebook1, mintpy_dir)
        shutil.copy(notebook2, mintpy_dir)
    except Exception as e:
        print(e)

    print('\nAll done, enjoy it!\n')


if __name__ == "__main__":
    main()

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
import sys

EXAMPLE = """Example:
  python3 MAI.py /ly/rslc /ly/stacking /ly/dem 20211229 200 60 20 5
  python3 MAI.py /ly/rslc /ly/stacking /ly/dem 20211229 200 60 20 5 -e rslc.deramp
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(description= 'MAI processing from RSLCs using GAMMA.',
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)

    parser.add_argument('rslc_dir', help='RSLCs directory')
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('dem_dir', help='DEM directory contains *.dem and *.dem.par')
    parser.add_argument('ref_slc', help='reference RSLC for making lookup table')
    parser.add_argument('max_sb', type=float, help='maximum spatial baseline')
    parser.add_argument('max_tb', type=float, help='maximum temporal baseline')
    parser.add_argument('rlks', help='range looks', type=int)
    parser.add_argument('alks', help='azimuth looks', type=int)
    parser.add_argument('-e', dest='slc_extension', type=str, default='.rslc', help='file extension for RSLCs (defaults: .rslc)')

    inps = parser.parse_args()
    return inps


def mk_tab(slc_dir, slc_tab, slc_extension):
    """Generate SLC_tab for processing

    Args:
        slc_dir (str): slc directory
        slc_tab (str): tab file
    """
    dates = sorted(
        [i for i in os.listdir(slc_dir) if re.findall(r'^\d{8}$', i)])
    with open(slc_tab, 'w+', encoding='utf-8') as f:
        for date in dates:
            slc = os.path.join(slc_dir, date, date + slc_extension)
            slc_par = slc + '.par'
            f.write(slc + '    ' + slc_par + '\n')


def mli_all(slc_tab, out_dir, rlks, alks):
    """Calculate MLI images for a stack of SLCs

    Args:
        slc_tab (str): slc tab file including slc slc_par
        out_dir (str): output directory
        rlks (int): range looks
        alks (int): azimuth looks
    """
    os.chdir(out_dir)
    with open(slc_tab, 'r') as f:
        for line in f.readlines():
            if line.strip():
                slc = line.strip().split()[0]
                slc_par = line.strip().split()[1]

                date = os.path.basename(slc)[0:8]

                mli = date + '.rmli'
                mli_par = mli + '.par'

                call_str = f"multi_look {slc} {slc_par} {mli} {mli_par} {rlks} {alks}"
                os.system(call_str)

                width = read_gamma_par(mli_par, 'range_samples')

                call_str = f"raspwr {mli} {width} 1 0 1 1 1. .35 1"
                os.system(call_str)


def base_calc(slc_tab, slc_par, max_sb, max_tb, out_dir):
    """Generate baseline output file with perpendicular baselines and delta_T values

    Args:
        slc_tab (str): slc tab file including slc and slc_par
        slc_par (str): reference slc par
        max_sb (float): max spatial baseline
        max_tb (float): max time baseline
        out_dir (str): output directory

    Returns:
        str: baseline file
    """
    os.chdir(out_dir)
    bperp_file = os.path.join(out_dir, 'bperp_file')
    call_str = f'base_calc {slc_tab} {slc_par} {bperp_file} itab 1 1 0 {max_sb} 0 {max_tb}'
    os.system(call_str)

    return bperp_file


def get_pairs(bperp_file):
    """Get pairs from baseline file

    Args:
        bperp_file (str): baseline file

    Returns:
        list: pairs
    """
    pairs = []
    with open(bperp_file, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            if line:
                split_list = line.strip().split()
                date1 = split_list[1]
                date2 = split_list[2]
                if int(date1) > int(date2):
                    pairs.append(date2 + '_' + date1)
                else:
                    pairs.append(date1 + '_' + date2)

    return pairs


def select_pairs_sbas(slc_tab, slc_par, max_sb, max_tb, out_dir):
    """Select pairs using sbas method

    Args:
        slc_tab (str): slc tab file including slc and slc_par
        slc_par (str): reference slc par
        max_sb (float): max spatial baseline
        max_tb (float): max time baseline
        out_dir (str): output directory

    Returns:
        list: pairs
    """
    bperp_file = base_calc(slc_tab, slc_par, max_sb, max_tb, out_dir)

    pairs = get_pairs(bperp_file)

    return pairs


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


def make_rdc_dem(mli, mli_par, dem, dem_par, out_dir):
    """make radar coordinate dem

    Args:
        mli (str): multi-looked slc
        mli_par (str): multi-looked slc par
        dem (str): dem
        dem_par (str): dem par
        out_dir (str): output directory

    Returns:
        str: rdc dem file
    """
    os.chdir(out_dir)

    date = os.path.basename(mli)[0:8]

    call_str = f"gc_map {mli_par} - {dem_par} {dem} dem_seg.par dem_seg lookup_table 1 1 sim_sar u v inc psi pix ls_map 8 1"
    os.system(call_str)

    call_str = f"pixel_area {mli_par} dem_seg.par dem_seg lookup_table ls_map inc pix_sigma0 pix_gamma0"
    os.system(call_str)

    width_mli = read_gamma_par(mli_par, 'range_samples')

    call_str = f"raspwr pix_gamma0 {width_mli} - - - - - - - pix_gamma0.bmp"
    os.system(call_str)

    call_str = f"create_diff_par {mli_par} - {date}.diff_par 1 0"
    os.system(call_str)

    call_str = f"offset_pwrm pix_sigma0 {mli} {date}.diff_par offs snr 64 64 offsets 2 100 100 5.0"
    os.system(call_str)

    call_str = f"offset_fitm offs snr {date}.diff_par coffs coffsets 5.0 1"
    os.system(call_str)

    width_utm_dem = read_gamma_par('dem_seg.par', 'width')

    call_str = f"gc_map_fine lookup_table {width_utm_dem} {date}.diff_par lookup_table_fine 1"
    os.system(call_str)

    length_mli = read_gamma_par(mli_par, 'azimuth_lines')
    width_mli = read_gamma_par(mli_par, 'range_samples')

    call_str = f"geocode lookup_table_fine dem_seg {width_utm_dem} {date}.dem {width_mli} {length_mli} 2 0"
    os.system(call_str)

    call_str = f"rashgt {date}.dem {mli} {width_mli} - - - - - 50 - - - {date}.dem.bmp"
    os.system(call_str)

    mli_name = os.path.basename(mli)

    call_str = f"geocode_back {mli} {width_mli} lookup_table_fine {mli_name}.geo {width_utm_dem} - 2 0"
    os.system(call_str)

    call_str = f"raspwr {mli_name}.geo {width_utm_dem} 1 0 1 1 1. .35 1 {mli_name}.geo.bmp"
    os.system(call_str)

    rdc_dem = os.path.join(out_dir, f"{date}.dem")

    return rdc_dem


def del_file(file):
    """Delete file

    Args:
        file (str): file
    """
    if os.path.isfile(file):
        os.remove(file)


def comb_pic(pic_path, out_pic):
    """Combine all pictures to one

    Args:
        pic_path (str): picture path
        out_pic (str): output picture
    """
    cmd_str = f"montage -label %f -geometry +5+7 -tile +6 -resize 300x300 {pic_path} {out_pic}"
    os.system(cmd_str)


def main():
    # get inputs
    inps = cmdline_parser()
    rslc_dir = os.path.abspath(inps.rslc_dir)
    out_dir = os.path.abspath(inps.out_dir)
    dem_dir = os.path.abspath(inps.dem_dir)
    ref_slc = inps.ref_slc
    max_sb = inps.max_sb
    max_tb = inps.max_tb
    rlks = inps.rlks
    alks = inps.alks
    slc_extension = inps.slc_extension

    # check rslc_dir
    if not os.path.isdir(rslc_dir):
        sys.exit('{} does not exist.'.format(rslc_dir))
    dates = sorted(
        [i for i in os.listdir(rslc_dir) if re.match(r'^\d{8}$', i)])
    if not dates:
        sys.exit('Cannot find RSLCs in {}'.format(rslc_dir))

    # check out_dir
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # check dem_dir
    if not os.path.isdir(dem_dir):
        sys.exit("{} does not exist.".format(dem_dir))
    else:
        dems = glob.glob(dem_dir + '/*.dem')
        dem_pars = [i + '.par' for i in dems]
        for i, j in zip(dems, dem_pars):
            if os.path.isfile(i) and os.path.isfile(j):
                dem = dems[0]
                dem_par = dem + '.par'
            else:
                sys.exit(f'Cannot find *.dem and *.dem.par in {dem_dir}.')

    # check ref_slc
    if ref_slc not in dates:
        sys.exit("RSLC for {} does not exist".format(ref_slc))

    # check extension
    if not slc_extension.startswith('.'):
        slc_extension = '.' + slc_extension

    m_rslc = os.path.join(rslc_dir, ref_slc, ref_slc + slc_extension)
    m_rslc_par = m_rslc + '.par'

    # select pairs
    base_dir = os.path.join(out_dir, 'base_calc')
    if not os.path.isdir(base_dir):
        os.mkdir(base_dir)

    slc_tab = os.path.join(base_dir, 'slc_tab')
    mk_tab(rslc_dir, slc_tab, slc_extension)

    pairs = select_pairs_sbas(slc_tab, m_rslc_par, max_sb, max_tb, base_dir)

    # multi-look
    mli_dir = os.path.join(out_dir, 'mli')
    if not os.path.isdir(mli_dir):
        os.mkdir(mli_dir)

    mli_all(slc_tab, mli_dir, rlks, alks)

    # making lookup table
    geo_dir = os.path.join(out_dir, 'geo')

    if not os.path.isdir(geo_dir):
        os.mkdir(geo_dir)

    sm_mli = os.path.join(mli_dir, ref_slc + '.rmli')
    sm_mli_par = sm_mli + '.par'

    _ = make_rdc_dem(sm_mli, sm_mli_par, dem, dem_par, geo_dir)

    diff_dir = os.path.join(out_dir, 'diff')

    if not os.path.isdir(diff_dir):
        os.mkdir(diff_dir)

    os.chdir(diff_dir)

    # diff and unwrap
    for pair in pairs:
        m_date = pair[0:8]
        s_date = pair[9:17]

        m_rslc = os.path.join(rslc_dir, m_date, m_date + slc_extension)
        m_rslc_par = m_rslc + '.par'

        s_rslc = os.path.join(rslc_dir, s_date, s_date + slc_extension)
        s_rslc_par = s_rslc + '.par'

        call_str = f"SBI_INT_F_F {m_rslc} {m_rslc_par} {s_rslc} {s_rslc_par} {pair}.sbi.int {pair}.off {m_date}.mli {m_date}.mli.par 0.5 {rlks} {alks} 0 1"
        os.system(call_str)

        width_mli = read_gamma_par(f"{m_date}.mli.par", 'range_samples')

        call_str = f"rasmph_pwr {pair}.sbi.int {m_date}.mli {width_mli} 1 1 0 1 1 1. 0.35 1 {pair}.sbi.int.bmp"
        os.system(call_str)

        call_str = f"adf {pair}.sbi.int {pair}.sbi.int.sm1 {pair}.sm.cc1 {width_mli} 0.3 64"
        os.system(call_str)
        call_str = f"adf {pair}.sbi.int.sm1 {pair}.sbi.int.sm2 {pair}.sm.cc2 {width_mli} 0.4 32"
        os.system(call_str)
        call_str = f"adf {pair}.sbi.int.sm2 {pair}.sbi.int.sm {pair}.sm.cc {width_mli} 0.5 16"
        os.system(call_str)

        del_file(f"{pair}.sbi.int.sm1")
        del_file(f"{pair}.sbi.int.sm2")
        del_file(f"{pair}.sm.cc1")
        del_file(f"{pair}.sm.cc2")

        call_str = f"rasmph_pwr {pair}.sbi.int.sm {m_date}.mli {width_mli} 1 1 0 1 1 1. 0.35 1 {pair}.sbi.int.sm.bmp"
        os.system(call_str)

        call_str = f"cpx_to_real {pair}.sbi.int.sm {pair}.sbi.int.sm.phase {width_mli} 4"
        os.system(call_str)

        call_str = f"sbi_offset {pair}.sbi.int.sm.phase {m_rslc_par}f {s_rslc_par}b {pair}.off {pair}.azi.off"
        os.system(call_str)

        call_str = f"rasrmg {pair}.azi.off {m_date}.mli {width_mli} 1 1 0 1 1 .6 1. .35 .0 1 {pair}.azi.off.bmp {pair}.sm.cc 1 .2"
        os.system(call_str)

        call_str = f"create_diff_par {m_date}.mli.par - {pair}.diff_par 1 0"
        os.system(call_str)

        call_str = f"quad_fit {pair}.azi.off {pair}.diff_par 32 32 - - 3"
        os.system(call_str)

        call_str = f"quad_sub {pair}.azi.off {pair}.diff_par {pair}.azi.off.sub 0 0"
        os.system(call_str)

        call_str = f"rasrmg {pair}.azi.off.sub {m_date}.mli {width_mli} 1 1 0 1 1 .6 1. .35 .0 1 {pair}.azi.off.sub.bmp {pair}.sm.cc 1 .2"
        os.system(call_str)

    comb_pic(diff_dir + '/*.sbi.int.bmp', out_dir + '/sbi.int.jpg')
    comb_pic(diff_dir + '/*.sbi.int.sm.bmp', out_dir + '/sbi.int.sm.jpg')
    comb_pic(diff_dir + '/*.azi.off.bmp', out_dir + '/azi.off.jpg')
    comb_pic(diff_dir + '/*.azi.off.sub.bmp', out_dir + '/azi.off.sub.jpg')


if __name__ == "__main__":
    main()

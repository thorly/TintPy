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
  **You have to follow the steps in order**

  # Step 1 (mk_diff_2d --> mk_adf_2d)
  python3 rslc2ifg_2d.py /ly/rslc /ly/dem /ly/diff_2d 10 2 20211229 1 -s 200 -t 60

  # Step 2 (set the starting point for unwrapping [disras *.adf.diff.bmp] --> mk_unw_2d)
  python3 rslc2ifg_2d.py /ly/rslc /ly/dem /ly/diff_2d 10 2 20211229 2 -r 100 -l 200

  # Step 3 (mk_base_2d --> mk_diff_2d --> mk_adf_2d --> mk_unw_2d)
  python3 rslc2ifg_2d.py /ly/rslc /ly/dem /ly/diff_2d 10 2 20211229 3 -r 100 -l 200
"""


def cmd_line_parser():
    parser = argparse.ArgumentParser(description='Calculate 2D diff, adf and unw using GAMMA.',
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)
    parser.add_argument('rslc_dir', help='directory path of RSLCs')
    parser.add_argument('dem_dir', help='directory path of xxx.dem and xxx.dem.par')
    parser.add_argument('out_dir', help='output directory for saving results')
    parser.add_argument('rlks', help='range looks', type=int)
    parser.add_argument('alks', help='azimuth looks', type=int)
    parser.add_argument('ref_slc', help='reference SLC for making rdc dem', type=str)
    parser.add_argument('step', help='step number', type=int, choices=[1, 2, 3])
    parser.add_argument('-s', dest='max_sb', help='maximum spatial baseline', type=float)
    parser.add_argument('-t', dest='max_tb', help='maximum temporal baseline', type=float)
    parser.add_argument('-r', dest='roff', help='offset to starting range of section to unwrap', type=float)
    parser.add_argument('-l', dest='loff', help='offset to starting line of section to unwrap', type=float)
    parser.add_argument('-e', dest='slc_extension', type=str, default='.rslc', help='file extension for RSLCs (defaults: .rslc)')
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
                value = line.split()[1].strip()

    return value


def mk_tab(slc_dir, slc_tab, slc_extension):
    """Generate SLC_tab for processing

    Args:
        slc_dir (str): slc directory
        slc_tab (str): tab file
        slc_extension (str): slc extension
    """
    dates = sorted([i for i in os.listdir(slc_dir) if re.findall(r'^\d{8}$', i)])
    with open(slc_tab, 'w+') as f:
        for date in dates:
            slc = os.path.join(slc_dir, date, date + slc_extension)
            slc_par = slc + '.par'
            f.write(slc + '    ' + slc_par + '\n')


def comb_pic(pic_path, out_pic):
    """Combine all pictures to one

    Args:
        pic_path (str): picture path
        out_pic (str): output picture
    """
    cmd_str = f"montage -label %f -geometry +5+7 -tile +6 -resize 300x300 {pic_path} {out_pic}"
    os.system(cmd_str)


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

    call_str = f"geocode_back {mli} {width_mli} lookup_table_fine {mli}.geo {width_utm_dem} - 2 0"
    os.system(call_str)

    call_str = f"raspwr {mli}.geo {width_utm_dem} 1 0 1 1 1. .35 1 {mli}.geo.bmp"
    os.system(call_str)

    rdc_dem = os.path.join(out_dir, f"{date}.dem")

    return rdc_dem


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


def main():
    # get inputs
    inps = cmd_line_parser()
    rslc_dir = os.path.abspath(inps.rslc_dir)
    dem_dir = os.path.abspath(inps.dem_dir)
    out_dir = os.path.abspath(inps.out_dir)
    rlks = inps.rlks
    alks = inps.alks
    ref_slc = inps.ref_slc
    max_sb = inps.max_sb
    max_tb = inps.max_tb
    roff = inps.roff
    loff = inps.loff
    step = inps.step
    extension = inps.slc_extension

    # check rslc_dir
    if not os.path.isdir(rslc_dir):
        sys.exit("{} does not exist.".format(rslc_dir))

    # check dem
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

    # check out_dir
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # get all date
    all_date = sorted([i for i in os.listdir(rslc_dir) if re.findall(r'^\d{8}$', i)])
    if len(all_date) < 2:
        sys.exit('No enough RSLCs.')

    # check ref_slc
    if re.findall(r'^\d{8}$', ref_slc):
        if ref_slc not in all_date:
            sys.exit('No rslc for {}.'.format(ref_slc))
    else:
        sys.exit('Error date for ref_slc.')

    # check extension
    if not extension.startswith('.'):
        extension = '.' + extension

    rmli_dir = os.path.join(out_dir, 'mli')
    rslc_tab = os.path.join(rmli_dir, 'rslc_tab')

    m_rmli = os.path.join(rmli_dir, f"{ref_slc}.rmli")
    m_rmli_par = m_rmli + '.par'

    geo_dir = os.path.join(out_dir, 'geo')
    rdc_dem = os.path.join(geo_dir, f"{ref_slc}.dem")

    base_dir = os.path.join(out_dir, 'base_calc')
    diff1_dir = os.path.join(out_dir, 'diff')
    diff2_dir = os.path.join(out_dir, 'diff_refine')

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Step 1 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if step == 1:
        if max_sb and max_tb:
            # multilook RSLCs
            if not os.path.isdir(rmli_dir):
                os.mkdir(rmli_dir)

            # write rslc_tab
            mk_tab(rslc_dir, rslc_tab, extension)

            mli_all(rslc_tab, rmli_dir, rlks, alks)

            # make rdc dem
            if not os.path.isdir(geo_dir):
                os.mkdir(geo_dir)

            rdc_dem = make_rdc_dem(m_rmli, m_rmli_par, dem, dem_par, geo_dir)

            # base_calc
            if not os.path.isdir(base_dir):
                os.mkdir(base_dir)

            m_rslc_par = os.path.join(rslc_dir, ref_slc, ref_slc + extension + '.par')

            os.chdir(base_dir)
            cmd_str = f"base_calc {rslc_tab} {m_rslc_par} bperp_file itab 1 1 0 {max_sb} 0 {max_tb}"
            os.system(cmd_str)

            # mk_diff_2d
            if not os.path.isdir(diff1_dir):
                os.mkdir(diff1_dir)

            cmd_str = f"mk_diff_2d {rslc_tab} itab 0 {rdc_dem} - {m_rmli} {rmli_dir} {diff1_dir} {rlks} {alks} 5 1 1 0"
            os.system(cmd_str)

            # mk_adf_2d
            cmd_str = f"mk_adf_2d {rslc_tab} itab {m_rmli} {diff1_dir} 5 0.75 32 4"
            os.system(cmd_str)

            # com_pic
            comb_pic(diff1_dir + '/*.adf.cc.bmp', out_dir + '/adf.cc.jpg')
            comb_pic(diff1_dir + '/*.adf.diff.bmp', out_dir + '/adf.diff.jpg')

            print(
                '\nStep 1 is done, you can run step 2 now.\npython3 diff_2d.py <rslc_dir> <dem_dir> <out_dir> <rlks> <alks> <ref_slc> 2 -r <roff> -l <loff>\n'
            )
        else:
            sys.exit('max_sb max_tb are required parameters for step 1.')

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Step 2 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if step == 2:
        if roff and loff:
            os.chdir(base_dir)

            # mk_unw_2d
            cmd_str = f"mk_unw_2d {rslc_tab} itab {m_rmli} {diff1_dir} 0 0 3 1 1 1 {roff} {loff} 1"
            os.system(cmd_str)

            # com_pic
            comb_pic(diff1_dir + '/*.adf.unw.bmp', out_dir + '/adf.unw.jpg')

            print(
                '\nStep 2 is done, you can run step 3 now.\npython3 diff_2d.py <rslc_dir> <dem_dir> <out_dir> <rlks> <alks> <ref_slc> 3 -r <roff> -l <loff>\n'
            )
        else:
            sys.exit('roff loff are required parameters for step 2.')

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Step 3 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if step == 3:
        if roff and loff:
            if not os.path.isdir(diff2_dir):
                os.mkdir(diff2_dir)

            os.chdir(base_dir)

            # mk_base_2d
            cmd_str = f"mk_base_2d {rslc_tab} itab {rdc_dem} {diff1_dir} pbase - 32 32 7 1"
            os.system(cmd_str)

            # copy base file to diff2_dir
            cmd_str = f"cp {os.path.join(diff1_dir, '*.base')} {diff2_dir}"
            os.system(cmd_str)

            # mk_diff_2d
            cmd_str = f"mk_diff_2d {rslc_tab} itab 0 {rdc_dem} - {m_rmli} {rmli_dir} {diff2_dir} {rlks} {alks} 5 1 1 0"
            os.system(cmd_str)

            # mk_adf_2d
            cmd_str = f"mk_adf_2d {rslc_tab} itab {m_rmli} {diff2_dir} 5 0.75 32 4"
            os.system(cmd_str)

            # mk_unw_2d
            cmd_str = f"mk_unw_2d {rslc_tab} itab {m_rmli} {diff2_dir} 0 0 3 1 1 1 {roff} {loff} 1 - 0 0 - - -d diff_tab"
            os.system(cmd_str)

            # com_pic
            comb_pic(diff2_dir + '/*.adf.cc.bmp', out_dir + '/adf.cc_refine.jpg')
            comb_pic(diff2_dir + '/*.adf.diff.bmp', out_dir + '/adf.diff_refine.jpg')
            comb_pic(diff2_dir + '/*.adf.unw.bmp', out_dir + '/adf.unw_refine.jpg')

            print('\nAll done, enjoy it!\n')
        else:
            sys.exit(
                'roff loff are required parameters for step 3.'
            )


if __name__ == "__main__":
    main()

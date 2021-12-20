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
  python3 s1_coreg.py /ly/slc /ly/rslc /ly/dem 20211229 2 -r 8 -a 2
  python3 s1_coreg.py /ly/slc /ly/rslc /ly/dem 20211229 1 2 -r 8 -a 2
  python3 s1_coreg.py /ly/slc /ly/rslc /ly/dem 20211229 1 2 3 -r 8 -a 2 -f t
"""


def cmd_line_parser():
    parser = argparse.ArgumentParser(description='Coregistrate all of the Sentinel-1 TOPS SLCs to a reference SLC.',
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)

    parser.add_argument('slc_dir', help='SLCs directory')
    parser.add_argument('rslc_dir', help='RSLCs directory')
    parser.add_argument('dem_dir', help='dem directory including *.dem and *.dem.par')
    parser.add_argument('ref_slc', help='reference SLC', type=str)
    parser.add_argument('sub_swath', type=int, nargs='+', choices=[1, 2, 3], help='sub_swath number for coregistration')
    parser.add_argument('-r', dest='rlks', help='range looks', type=int, required=True)
    parser.add_argument('-a', dest='alks', help='azimuth looks', type=int, required=True)
    parser.add_argument('-p', dest='pol', help='polarization(defaults: vv)', choices=['vv', 'vh'], default='vv')
    parser.add_argument('-f', dest='deramp_flag', choices=['t', 'T', 'f', 'F'], help='flag for deramp later (t for YES, f for No, defaults: f)', default='f')
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


def copy_file(src_file, dst_file):
    """Copy file

    Args:
        src_file (str): src file
        dst_file (str): dst file
    """
    if not os.path.isfile(dst_file):
        shutil.copy(src_file, dst_file)


def make_rdc_dem(slc, slc_par, dem, dem_par, rlks, alks, out_dir):
    """make radar coordinate dem

    Args:
        slc (str): slc
        slc_par (str): slc par
        dem (str): dem
        dem_par (str): dem par
        rlks (int): range looks
        alks (int): azimuth looks
        out_dir (str): output directory

    Returns:
        str: rdc dem file
    """
    os.chdir(out_dir)

    date = os.path.basename(slc)[0:8]

    mli = f"{date}.mli"
    mli_par = mli + '.par'

    call_str = f"multi_look {slc} {slc_par} {mli} {mli_par} {rlks} {alks}"
    os.system(call_str)

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

    mli_name = os.path.basename(mli)

    call_str = f"geocode_back {mli} {width_mli} lookup_table_fine {mli_name}.geo {width_utm_dem} - 2 0"
    os.system(call_str)

    call_str = f"raspwr {mli_name}.geo {width_utm_dem} 1 0 1 1 1. .35 1 {mli_name}.geo.bmp"
    os.system(call_str)

    length_mli = read_gamma_par(mli_par, 'azimuth_lines')

    call_str = f"geocode lookup_table_fine dem_seg {width_utm_dem} {date}.dem {width_mli} {length_mli} 2 0"
    os.system(call_str)

    call_str = f"rashgt {date}.dem {mli} {width_mli} 1 1 0 1 1 160.0 1. .35 1 {date}.dem.bmp"
    os.system(call_str)

    rdc_dem = os.path.join(out_dir, f"{date}.dem")

    return rdc_dem


def write_tab(slcs, tab_file):
    """Write path of slc slc_par and tops_par to tab file

    Args:
        slcs (str list): slc files
        tab_file (str): output file
    """
    with open(tab_file, 'w+', encoding='utf-8') as f:
        for slc in slcs:
            f.write(f'{slc} {slc}.par {slc}.tops_par\n')


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


def slc_mosaic(date_slc_dir, sub_swaths, pol, rlks, alks, out_dir):
    """Calculate SLC mosaic of Sentinel-1 TOPS burst SLC data

    Args:
        date_slc_dir (str): slc directory
        sub_swaths (list): sub_swath number
        pol (str): polarization
        rlks (int): range looks
        alks (int): azimuth looks
        out_dir (str): output directory

    Returns:
        tuple: mosaiced slc and par
    """
    date = os.path.basename(date_slc_dir)

    slc_tab = os.path.join(date_slc_dir, 'slc_tab')
    with open(slc_tab, 'w+', encoding='utf-8') as f:
        for i in sub_swaths:
            slc = os.path.join(date_slc_dir, f'{date}.iw{i}.{pol}.slc')
            slc_par = slc + '.par'
            tops_par = slc + '.tops_par'
            f.write(f'{slc} {slc_par} {tops_par}\n')

    slc_out = os.path.join(out_dir, date + '.slc')
    slc_par_out = os.path.join(out_dir, date + '.slc.par')

    call_str = f"SLC_mosaic_S1_TOPS {slc_tab} {slc_out} {slc_par_out} {rlks} {alks}"
    os.system(call_str)

    os.remove(slc_tab)

    return slc_out, slc_par_out


def print_coreg_quality(quality_files):
    """print coregistration quality

    Args:
        quality_files (str list): quality files
    """
    print('-' * 37)
    print('|{:^35}|'.format('coregistration quality report'))
    print('-' * 37)
    print('|{:^16}|{:^10}|       |'.format('date', 'daz10000'))
    print('-' * 37)

    for file in sorted(quality_files):
        with open(file, 'r', encoding='utf-8') as f:
            for i in f.readlines()[::-1]:
                if i.startswith('azimuth_pixel_offset'):
                    daz = i.split()[1]
                    date = file[-22:-14]
                    daz10000 = float(daz) * 10000
                    daz10000 = round(daz10000, 2)
                    if daz10000 > 5:
                        print('|{:^16}|{:^10}|  >  5 |'.format(date, daz10000))
                        print('-' * 37)
                    elif daz10000 < -5:
                        print('|{:^16}|{:^10}|  < -5 |'.format(date, daz10000))
                        print('-' * 37)
                    else:
                        print('|{:^16}|{:^10}|       |'.format(date, daz10000))
                        print('-' * 37)
                    break


def main():
    # get inputs
    inps = cmd_line_parser()
    slc_dir = os.path.abspath(inps.slc_dir)
    rslc_dir = os.path.abspath(inps.rslc_dir)
    dem_dir = os.path.abspath(inps.dem_dir)
    sub_swath = inps.sub_swath
    pol = inps.pol
    rlks = inps.rlks
    alks = inps.alks
    ref_slc = inps.ref_slc
    flag = inps.deramp_flag.lower()

    # check slc_dir
    if not os.path.isdir(slc_dir):
        sys.exit("{} does not exist.".format(slc_dir))

    # check rslc_dir
    if not os.path.isdir(rslc_dir):
        os.mkdir(rslc_dir)

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
                break
            else:
                sys.exit(f'Cannot find *.dem and *.dem.par in {dem_dir}.')

    # get all date
    dates = sorted(
        [i for i in os.listdir(slc_dir) if re.findall(r'^\d{8}$', i)])
    if len(dates) < 2:
        sys.exit('No enough SLCs.')

    # check ref_slc
    if re.findall(r'^\d{8}$', ref_slc):
        if not ref_slc in dates:
            sys.exit('No slc for {}.'.format(ref_slc))
    else:
        sys.exit('Error date for ref_slc.')

    # copy reference iw slc to rslc_dir for deramp
    m_date = ref_slc
    m_slc_dir = os.path.join(slc_dir, m_date)
    m_rslc_dir = os.path.join(rslc_dir, m_date)

    if not os.path.isdir(m_rslc_dir):
        os.mkdir(m_rslc_dir)

    if flag == 't':
        print('Copy reference iw slc to {}'.format(m_rslc_dir))
        for i in sub_swath:
            iw_slc = os.path.join(m_slc_dir, f'{m_date}.iw{i}.{pol}.slc')
            iw_slc_par = iw_slc + '.par'
            iw_slc_tops_par = iw_slc + '.tops_par'

            iw_rslc = os.path.join(m_rslc_dir, f'{m_date}.iw{i}.{pol}.rslc')
            iw_rslc_par = iw_rslc + '.par'
            iw_rslc_tops_par = iw_rslc + '.tops_par'

            copy_file(iw_slc, iw_rslc)
            copy_file(iw_slc_par, iw_rslc_par)
            copy_file(iw_slc_tops_par, iw_rslc_tops_par)

    # get slave dates
    s_dates = dates.copy()
    s_dates.remove(m_date)

    # make rdr dem
    geo_dir = os.path.join(rslc_dir, 'geo')

    if os.path.isdir(geo_dir):
        shutil.rmtree(geo_dir)

    os.mkdir(geo_dir)

    m_slc, m_slc_par = slc_mosaic(m_slc_dir, sub_swath, pol, rlks, alks, geo_dir)

    rdc_dem = make_rdc_dem(m_slc, m_slc_par, dem, dem_par, rlks, alks, geo_dir)

    copy_ref_flag = True

    for s_date in s_dates:
        s_slc_dir = os.path.join(slc_dir, s_date)

        s_rslc_dir = os.path.join(rslc_dir, s_date)

        if not os.path.isdir(s_rslc_dir):
            os.mkdir(s_rslc_dir)

        os.chdir(s_rslc_dir)

        s_iw_slcs = []
        m_iw_slcs = []
        s_iw_rslcs = []

        for i in sub_swath:
            s_iw_slc = os.path.join(s_slc_dir, f"{s_date}.iw{i}.{pol}.slc")
            m_iw_slc = os.path.join(m_slc_dir, f"{m_date}.iw{i}.{pol}.slc")
            s_iw_rslc = os.path.join(s_rslc_dir, f"{s_date}.iw{i}.{pol}.rslc")

            s_iw_slcs.append(s_iw_slc)
            m_iw_slcs.append(m_iw_slc)
            s_iw_rslcs.append(s_iw_rslc)

        write_tab(s_iw_slcs, 'slc_tab_s')
        write_tab(m_iw_slcs, 'slc_tab_m')
        write_tab(s_iw_rslcs, 'rslc_tab')

        call_str = f"S1_coreg_TOPS slc_tab_m {m_date} slc_tab_s {s_date} rslc_tab {rdc_dem} {rlks} {alks} - - 0.7 0.001 0.7 1"
        os.system(call_str)

        if copy_ref_flag:
            dst_rslc = os.path.join(m_rslc_dir, m_date + '.rslc')
            dst_rslc_par = dst_rslc + '.par'
            if not os.path.isfile(dst_rslc):
                shutil.move(m_date + '.rslc', m_rslc_dir)
            if not os.path.isfile(dst_rslc_par):
                shutil.move(m_date + '.rslc.par', m_rslc_dir)
            copy_ref_flag = False

        # clean s_rslc dir
        save_files = []
        save_files.append(s_date + '.rslc')
        save_files.append(s_date + '.rslc.par')
        save_files.append(m_date + '_' + s_date + '.coreg_quality')
        save_files.append(m_date + '_' + s_date + '.diff.bmp')

        # save iw rslc for deramp
        if flag == 't':
            for i in sub_swath:
                iw_rslc = f'{s_date}.iw{i}.{pol}.rslc'
                iw_rslc_par = iw_rslc + '.par'
                iw_rslc_tops_par = iw_rslc + '.tops_par'
                save_files.append(iw_rslc)
                save_files.append(iw_rslc_par)
                save_files.append(iw_rslc_tops_par)

        for file in os.listdir(s_rslc_dir):
            if file not in save_files:
                os.remove(file)

    # del geo_dir
    if os.path.isdir(geo_dir):
        shutil.rmtree(geo_dir)

    # generate bmp for rslc
    rslc_files = glob.glob(rslc_dir + '/*/*.rslc')
    for rslc in rslc_files:
        slc2bmp(rslc, rslc + '.par', rlks, alks, rslc + '.bmp')

    # check coreg_quality
    quality_files = glob.glob(rslc_dir + '/*/*.coreg_quality')
    print_coreg_quality(quality_files)

    print('\nAll done, enjoy it!\n')


if __name__ == "__main__":
    main()

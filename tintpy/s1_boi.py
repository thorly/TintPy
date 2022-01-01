#!/usr/bin/env python3
################################
# Program is part of TintPy    #
# Copyright (c) 2022, Lei Yuan #
# Author: Lei Yuan, 2022       #
################################

import argparse
import glob
import os
import re
import sys

import numpy as np

EXAMPLE = """Example:
  python3 s1_boi.py /ly/rslc /ly/BOI /ly/dem 20211229 1 -r 20 -a 5 -s 200 -t 60
  python3 s1_boi.py /ly/rslc /ly/BOI /ly/dem 20211229 1 2 -r 20 -a 5 -s 200 -t 60
  python3 s1_boi.py /ly/rslc /ly/BOI /ly/dem 20211229 1 2 3 -r 20 -a 5 -s 200 -t 60
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(
        description=
        'Sentinel-1 burst overlap double-difference interferometry from RSLCs using GAMMA.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLE)

    parser.add_argument('rslc_dir',
                        help='RSLCs directory (iw*.**.rslc and *.rslc)')
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('dem_dir',
                        help='DEM directory contains *.dem and *.dem.par')
    parser.add_argument('ref_slc', help='reference RSLC for making rdc dem')
    parser.add_argument('sub_swath',
                        type=int,
                        nargs='+',
                        choices=[1, 2, 3],
                        help='sub_swath number')
    parser.add_argument('-r',
                        dest='rlks',
                        help='range looks',
                        type=int,
                        required=True)
    parser.add_argument('-a',
                        dest='alks',
                        help='azimuth looks',
                        type=int,
                        required=True)
    parser.add_argument('-s',
                        dest='max_sb',
                        type=float,
                        help='maximum spatial baseline',
                        required=True)
    parser.add_argument('-t',
                        dest='max_tb',
                        type=float,
                        help='maximum temporal baseline',
                        required=True)
    parser.add_argument('-e',
                        dest='slc_extension',
                        type=str,
                        default='rslc',
                        help='file extension for RSLCs (defaults: rslc)')
    parser.add_argument('-p',
                        dest='pol',
                        help='polarization(defaults: vv)',
                        choices=['vv', 'vh'],
                        default='vv')

    inps = parser.parse_args()
    return inps


def mk_tab(slc_dir, slc_tab, slc_extension):
    """Generate SLC_tab for processing

    Args:
        slc_dir (str): slc directory
        slc_tab (str): tab file
        slc_extension (str): slc extension
    """
    dates = sorted(
        [i for i in os.listdir(slc_dir) if re.findall(r'^\d{8}$', i)])
    with open(slc_tab, 'w+', encoding='utf-8') as f:
        for date in dates:
            slc = os.path.join(slc_dir, date, date + '.' + slc_extension)
            slc_par = slc + '.par'
            f.write(slc + '    ' + slc_par + '\n')


def write_tab(rslc_dir, date, sub_swath, pol, extension, tab_file):
    """Write path of slc slc_par and tops_par to tab file

    Args:
        rslc_dir (str): rslc directory
        date (str): date
        sub_swath (int list): sub_swath number
        pol (str): polarization
        extension (str): file extension for RSLC
        tab_file (str): output file
    """
    with open(tab_file, 'w+', encoding='utf-8') as f:
        for i in sub_swath:
            iw_slc = os.path.join(rslc_dir, date,
                                  f"{date}.iw{i}.{pol}.{extension}")
            iw_slc_par = iw_slc + '.par'
            iw_slc_tops_par = iw_slc + '.tops_par'
            f.write(f"{iw_slc} {iw_slc_par} {iw_slc_tops_par}\n")


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


def write_gamma(data, out_file, file_type):
    """Write GAMMA format file

    Args:
        data (array): array
        out_file (str): output file name
        file_type (str): file type
    """
    data = data.astype(file_type)
    data.byteswap('True')
    data.reshape(-1, 1)
    data.tofile(out_file)


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


def make_lookup(slc, slc_par, dem, dem_par, rlks, alks, out_dir):
    """make lookup table

    Args:
        slc (str): slc
        slc_par (str): slc par
        dem (str): dem
        dem_par (str): dem par
        rlks (int): range looks
        alks (int): azimuth looks
        out_dir (str): output directory

    Returns:
        str: multi-looked par
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

    return os.path.join(out_dir, mli_par)


def del_file(file):
    """Delete file

    Args:
        file (str): file
    """
    if os.path.isfile(file):
        os.remove(file)


def mk_dir(dir):
    """Make directory

    Args:
        dir (str): directory
    """
    if not os.path.isdir(dir):
        os.mkdir(dir)


def get_dates(rslc_dir):
    """Get dates using rslc_dir

    Args:
        rslc_dir (str): RSLC directory

    Returns:
        list: sorted dates
    """
    dates = sorted([i for i in os.listdir(rslc_dir) if re.findall(r'^\d{8}$', i)])

    return dates


def get_overlap_poly(rslc_dir, date, sub_swath, pol, extension, rlks, alks, flag, out_dir):
    """Get polygons for overlap regions in Sentinel-1 TOPS data

    Args:
        rslc_dir (str): RSLC directory
        date (str): date
        sub_swath (list): sub_swath number
        pol (str): polariaztion
        extension (str): rslc extension
        rlks (int): range looks
        alks (int): azimuth looks
        flat (int): flag for azimuth(1) or range overlap(2)
        out_dir (str): output directory

    Returns:
        array: overlap polygons
    """
    tab_file = os.path.join(out_dir, f"{date}_tab")
    write_tab(rslc_dir, date, sub_swath, pol, extension, tab_file)

    poly_file = os.path.join(out_dir, f"{date}.azi.poly")
    if flag == 2:
        poly_file = os.path.join(out_dir, f"{date}.range.poly")

    cmd_str = f"S1_poly_overlap3 {tab_file} {rlks} {alks} {poly_file} {flag}"
    os.system(cmd_str)

    poly = np.loadtxt(poly_file, np.int16)

    return poly


def extract_burst_overlap(rslc_dir, sub_swath, pol, extension, rlks, alks, out_bo_dir):
    """Extract Sentinel-1 azimuth overlap region

    Args:
        rslc_dir (str): RSLC directory
        sub_swath (list): sub_swath number
        pol (array): overlap polygons
        extension (str): rslc extension
        rlks (int): range looks
        alks (int): azimuth looks
        out_bo_dir (str): output directory

    Returns:
        list: burst numbers for per sub_swath
    """
    bursts_numbers = []
    dates = get_dates(rslc_dir)
    for date in dates:
        rslc_date_dir = os.path.join(rslc_dir, date)

        out_date_dir = os.path.join(out_bo_dir, date)
        mk_dir(out_date_dir)

        for i in sub_swath:
            iw_slc = os.path.join(rslc_date_dir, f"{date}.iw{i}.{pol}.{extension}")
            iw_slc_par = iw_slc + '.par'
            iw_slc_tops_par = iw_slc + '.tops_par'

            # lines offset between start of burst1 and start of burst2
            azimuth_line_time = float(read_gamma_par(iw_slc_par, 'azimuth_line_time'))
            burst_start_time_1 = float(read_gamma_par(iw_slc_tops_par, 'burst_start_time_1'))
            burst_start_time_2 = float(read_gamma_par(iw_slc_tops_par, 'burst_start_time_2'))

            lines_offset = int(0.5 + (burst_start_time_2 - burst_start_time_1) / azimuth_line_time)

            bursts_number = int(read_gamma_par(iw_slc_tops_par, 'number_of_bursts'))
            bursts_numbers.append(bursts_number)

            lines_per_burst = int(read_gamma_par(iw_slc_tops_par, 'lines_per_burst'))
            lines_overlap = lines_per_burst - lines_offset

            range_samples = int(read_gamma_par(iw_slc_par, 'range_samples'))

            for j in range(1, bursts_number):
                # extract SLC sections for overlap region i (i=1 --> overlap between bursts 1 and 2)
                for k in range(1, 3):
                    starting_line = lines_offset + (j - 1) * lines_per_burst
                    if k == 2:
                        starting_line = j * lines_per_burst

                    slc_copy = os.path.join(
                        out_date_dir,
                        f"{date}.iw{i}.{pol}.{extension}.{j}.{k}")

                    cmd_str = f"SLC_copy {iw_slc} {iw_slc_par} {slc_copy} {slc_copy}.par - 1. 0 {range_samples} {starting_line} {lines_overlap}"
                    os.system(cmd_str)

                    cmd_str = f"rasSLC {slc_copy} {range_samples} 1 0 {rlks} {alks} 1. .35"
                    os.system(cmd_str)

    return bursts_numbers[0:len(sub_swath)]


def burst_overlap_ddi(bo_dir, pairs, sub_swath, pol, extension, bursts_numbers, rlks, alks, diff_dir):
    """Sentinel-1 burst overlap double-difference interferometry

    Args:
        bo_dir (str): burst overlap directory
        pairs (list): ifg pairs
        sub_swath (list): sub_swath number
        pol (str): polariaztion
        extension (str): rslc extension
        bursts_numbers (list): bursts numbers for per sub_swath
        rlks (int): range looks
        alks (int): azimuth looks
        diff_dir (str): directory for saving results
    """
    os.chdir(diff_dir)

    for pair in pairs:
        m_date = pair[0:8]
        s_date = pair[9:17]
        m_date_dir = os.path.join(bo_dir, m_date)
        s_date_dir = os.path.join(bo_dir, s_date)

        for i in sub_swath:
            index = sub_swath.index(i)

            for j in range(1, bursts_numbers[index]):
                for k in range(1, 3):
                    m_slc = os.path.join(m_date_dir, f"{m_date}.iw{i}.{pol}.{extension}.{j}.{k}")
                    s_slc = os.path.join(s_date_dir, f"{s_date}.iw{i}.{pol}.{extension}.{j}.{k}")

                    off = f"{m_date}_{s_date}.iw{i}.{j}.off{k}"
                    cmd_str = f"create_offset {m_slc}.par {m_slc}.par {off} 1 1 1 0"
                    os.system(cmd_str)

                    inf = f"{m_date}_{s_date}.iw{i}.{j}.int{k}"
                    # cmd_str = f"SLC_intf {m_slc} {s_slc} {m_slc}.par {s_slc}.par {off} {inf} 1 1 0 - 0 0"
                    cmd_str = f"SLC_intf {m_slc} {s_slc} {m_slc}.par {s_slc}.par {off} {inf} {rlks} {alks} 0 - 0 0"
                    os.system(cmd_str)

                diff_par = f"{m_date}_{s_date}.iw{i}.{j}.diff_par"
                off1, off2 = off[:-1] + '1', off[:-1] + '2'
                cmd_str = f"create_diff_par {off1} {off2} {diff_par} 0 0"
                os.system(cmd_str)

                range_samples = read_gamma_par(diff_par, 'range_samp_1')
                inf1, inf2 = inf[:-1] + '1', inf[:-1] + '2'
                cmd_str = f"cpx_to_real {inf2} tmp {range_samples} 4"
                os.system(cmd_str)

                diff = f"{m_date}_{s_date}.iw{i}.{j}.diff"
                cmd_str = f"sub_phase {inf1} tmp {diff_par} {diff} 1 0"
                os.system(cmd_str)

                del_file('tmp')

                diff20 = diff + '20'
                off20 = off[:-1] + '20'
                # cmd_str = f"multi_cpx {diff} {off1} {diff20} {off20} {rlks} {alks}"
                cmd_str = f"multi_cpx {diff} {off1} {diff20} {off20} 1 1"
                os.system(cmd_str)

                range_samples20 = int(read_gamma_par(off20, 'interferogram_width'))
                # azimuth_lines20 = int(read_gamma_par(off20, 'interferogram_azimuth_lines'))
                cmd_str = f"cc_wave {diff20} - - {diff20}.cc {range_samples20} 5 5 0"
                os.system(cmd_str)

                # adf filtering of double differential interferogram
                cmd_str = f"adf {diff20} {diff20}.adf1 {diff20}.adf.cc1 {range_samples20} 0.6 32 5 4"
                os.system(cmd_str)

                cmd_str = f"adf {diff20}.adf1 {diff20}.adf2 {diff20}.adf.cc2 {range_samples20} 0.4 16 5 2"
                os.system(cmd_str)

                cmd_str = f"adf {diff20}.adf2 {diff20}.adf {diff20}.adf.cc {range_samples20} 0.2 8 5 1"
                os.system(cmd_str)

                del_file(f"{diff20}.adf1")
                del_file(f"{diff20}.adf2")
                del_file(f"{diff20}.adf.cc1")
                del_file(f"{diff20}.adf.cc2")

                cmd_str = f"rascc {diff20}.adf.cc - {range_samples20} 1 1 0 1 1 - - 1. .35 1 {diff20}.adf.cc.bmp"
                os.system(cmd_str)

                cc_thresh = 0.6
                cmd_str = f"rascc_mask {diff20}.adf.cc - {range_samples20} 1 1 0 1 1 {cc_thresh} - 0.0 1.0 1. .35 1 {diff20}.adf.cc_mask.bmp"
                os.system(cmd_str)

                cmd_str = f"rasmph {diff20}.adf {range_samples20} 1 0 1 1 1. .35 1 {diff20}.adf.bmp"
                os.system(cmd_str)

                ref_x = 0
                ref_y = 0
                cmd_str = f"mcf {diff20}.adf {diff20}.cc {diff20}.adf.cc_mask.bmp {diff20}.phase {range_samples20} 0 0 0 - - 1 1 512 {ref_x} {ref_y}"
                os.system(cmd_str)

                cmd_str = f"rasrmg {diff20}.phase - {range_samples20} 1 1 0 1 1 .5 1. .35 0.0 1 {diff20}.phase.bmp {diff20}.adf.cc 1 .2"
                os.system(cmd_str)


def merge_data(diff_dir, pairs, sub_swath, bursts_numbers, azi_poly, width, lines):
    """Merge overlap region into one file for phase and cc

    Args:
        diff_dir (str): diff directory
        pairs (list): ifg pairs
        sub_swath (list): sub_swath nubmer
        bursts_numbers (list): bursts numbers for per sub_swath
        azi_poly (array): azimuth overlap polygons
        width (int): width for merged file
        lines (int): lines for merged file
    """
    merged_phase = np.ones((lines, width), dtype='float32') * 99999
    merged_cc = merged_phase.copy()

    os.chdir(diff_dir)

    for pair in pairs:
        for i in sub_swath:
            index = sub_swath.index(i)
            bursts_number = bursts_numbers[index]

            for j in range(1, bursts_number):
                off20 = f"{pair}.iw{i}.{j}.off20"
                line_sm = int(read_gamma_par(off20, 'interferogram_azimuth_lines'))
                width_sm = int(read_gamma_par(off20, 'interferogram_width'))

                line_start = azi_poly[4 * j - 4 * (index + 1) + 4 * sum(bursts_numbers[0:index]), 1] - 1
                range_start = azi_poly[4 * j - 4 * (index + 1) + 4 * sum(bursts_numbers[0:index]), 0] - 1

                phase_file = f"{pair}.iw{i}.{j}.diff20.phase"
                phase = read_gamma(phase_file, line_sm, 'float32')
                merged_phase[line_start:line_start + line_sm, range_start:range_start + width_sm] = phase

                cc_file = f"{pair}.iw{i}.{j}.diff20.adf.cc"
                cc = read_gamma(cc_file, line_sm, 'float32')
                merged_cc[line_start:line_start + line_sm, range_start:range_start + width_sm] = cc

        merged_phase_file = f"{pair}.boi.disp"
        write_gamma(merged_phase, merged_phase_file, 'float32')

        merged_cc_file = f"{pair}.boi.cc"
        write_gamma(merged_cc, merged_cc_file, 'float32')


def main():
    # get inputs
    inps = cmdline_parser()
    rslc_dir = os.path.abspath(inps.rslc_dir)
    out_dir = os.path.abspath(inps.out_dir)
    dem_dir = os.path.abspath(inps.dem_dir)
    ref_slc = inps.ref_slc
    sub_swath = inps.sub_swath
    rlks = inps.rlks
    alks = inps.alks
    max_sb = inps.max_sb
    max_tb = inps.max_tb
    slc_extension = inps.slc_extension
    pol = inps.pol

    # check rslc_dir
    if not os.path.isdir(rslc_dir):
        sys.exit('{} does not exist.'.format(rslc_dir))
    dates = get_dates(rslc_dir)
    if not dates:
        sys.exit('Cannot find RSLCs in {}'.format(rslc_dir))

    # check out_dir
    mk_dir(out_dir)

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
    if re.findall(r'^\d{8}$', ref_slc):
        if not ref_slc in dates:
            sys.exit('No slc for {}.'.format(ref_slc))
    else:
        sys.exit('Error date for ref_slc.')

    # check extension
    if slc_extension.startswith('.'):
        slc_extension = slc_extension[1:]

    # check sub_swath
    error_date = {}
    for i in sub_swath:
        error_date[i] = []
        for date in dates:
            iw_slc = os.path.join(rslc_dir, date, f'{date}.iw{i}.{pol}.{slc_extension}')
            iw_slc_par = iw_slc + '.par'
            iw_slc_tops_par = iw_slc + '.tops_par'

            e1 = os.path.isfile(iw_slc)
            e2 = os.path.isfile(iw_slc_par)
            e3 = os.path.isfile(iw_slc_tops_par)

            if not (e1 and e2 and e3):
                error_date[i].append(date)

    if error_date[list(error_date.keys())[0]]:
        for key in error_date.keys():
            for date in error_date[key]:
                print(f'No slc or slc_par or tops_par for {date} sub_swath {key}')
        sys.exit('\nPlease check it.')

    # make lookup table
    geo_dir = os.path.join(out_dir, 'geo')
    mk_dir(geo_dir)

    m_rslc = os.path.join(rslc_dir, ref_slc, f"{ref_slc}.{slc_extension}")
    m_rslc_par = m_rslc + '.par'
    mli_par = make_lookup(m_rslc, m_rslc_par, dem, dem_par, rlks, alks, geo_dir)

    # select pairs
    base_dir = os.path.join(out_dir, 'base_calc')
    mk_dir(base_dir)

    slc_tab = os.path.join(base_dir, 'slc_tab')
    mk_tab(rslc_dir, slc_tab, slc_extension)

    pairs = select_pairs_sbas(slc_tab, m_rslc_par, max_sb, max_tb, base_dir)

    # extract overlap area
    bo_dir = os.path.join(out_dir, 'burst_overlap')
    mk_dir(bo_dir)

    bursts_numbers = extract_burst_overlap(rslc_dir, sub_swath, pol, slc_extension, rlks, alks, bo_dir)

    # burst overlap double-difference interferometry
    diff_dir = os.path.join(out_dir, 'diff')
    mk_dir(diff_dir)

    burst_overlap_ddi(bo_dir, pairs, sub_swath, pol, slc_extension, bursts_numbers, rlks, alks, diff_dir)

    # merge data
    azi_poly = get_overlap_poly(rslc_dir, ref_slc, sub_swath, pol, slc_extension, rlks, alks, 1, diff_dir)
    width = int(read_gamma_par(mli_par, 'range_samples'))
    lines = int(read_gamma_par(mli_par, 'azimuth_lines'))

    merge_data(diff_dir, pairs, sub_swath, bursts_numbers, azi_poly, width, lines)

    print("\nAll done, enjoy it!\n")

if __name__ == "__main__":
    main()

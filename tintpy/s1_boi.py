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

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

from subprocess import Popen, PIPE

EXAMPLE = """Example:
  python3 s1_boi.py /ly/rslc /ly/BOI /ly/dem 20211229 1 -r 20 -a 5 -s 200 -t 60
  python3 s1_boi.py /ly/rslc /ly/BOI /ly/dem 20211229 1 2 -r 20 -a 5 -s 200 -t 60
  python3 s1_boi.py /ly/rslc /ly/BOI /ly/dem 20211229 1 2 3 -r 20 -a 5 -s 200 -t 60
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(
        description='Sentinel-1 burst overlap double-difference interferometry from RSLCs using GAMMA.',
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)

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
    dates = sorted([i for i in os.listdir(slc_dir) if re.findall(r'^\d{8}$', i)])
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
            iw_slc = os.path.join(rslc_dir, date, f"{date}.iw{i}.{pol}.{extension}")
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


def make_lookup(mli, mli_par, dem, dem_par, out_dir):
    """make lookup table

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

                    slc_copy = os.path.join(out_date_dir, f"{date}.iw{i}.{pol}.{extension}.{j}.{k}")

                    cmd_str = f"SLC_copy {iw_slc} {iw_slc_par} {slc_copy} {slc_copy}.par - 1. 1 {range_samples} {starting_line} {lines_overlap}"
                    os.system(cmd_str)

                    cmd_str = f"rasSLC {slc_copy} {range_samples} 1 0 {rlks} {alks} 1. .35"
                    os.system(cmd_str)

    return bursts_numbers[0:len(sub_swath)]


def burst_overlap_ddi(bo_dir, pairs, sub_swath, pol, extension, bursts_numbers, rlks, alks, ddi_dir):
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
        ddi_dir (str): directory for saving results
    """
    os.chdir(ddi_dir)

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


def get_unit(rslc_dir, date, sub_swath, extension, pol):
    """Get unit of transforming BOI displacement

    Args:
        rslc_dir (str): RSLC directory
        date (str): date
        sub_swath (list): sub_swath number
        extension (str): slc extension
        pol (str): polarization

    Returns:
        list: unit
    """
    unit = []

    # using GAMMA command `af_SLC` to get it
    Vs = 7180.5433

    slc_par = os.path.join(rslc_dir, date, f"{date}.{extension}.par")

    radar_frequency = float(read_gamma_par(slc_par, 'radar_frequency'))
    light_speed = 299792458
    wavelength = light_speed / radar_frequency

    prf = float(read_gamma_par(slc_par, 'prf'))
    Ts = 1 / prf

    azimuth_pixel_spacing = float(read_gamma_par(slc_par, 'azimuth_pixel_spacing'))

    for i in sub_swath:
        iw_slc = os.path.join(rslc_dir, date, f'{date}.iw{i}.{pol}.{extension}')
        iw_slc_par = iw_slc + '.par'
        iw_slc_tops_par = iw_slc + '.tops_par'

        burst_time1 = float(read_gamma_par(iw_slc_tops_par, 'burst_start_time_1'))
        burst_time2 = float(read_gamma_par(iw_slc_tops_par, 'burst_start_time_2'))
        Tc = burst_time2 - burst_time1

        az_steering_rate = abs(float(read_gamma_par(iw_slc_tops_par, 'az_steering_rate')))
        center_range = float(read_gamma_par(iw_slc_par, 'center_range_slc'))

        Ka = -2 * Vs * Vs / wavelength / center_range

        Ks = 2 * np.deg2rad(az_steering_rate) * Vs / wavelength

        Kt = Ka * Ks / (Ka - Ks)

        fovl = abs(Kt) * Tc

        u = 1 / (2 * np.pi * fovl * Ts / azimuth_pixel_spacing)

        unit.append(u)

    return unit


def draw_img(data, out_file, cmap='jet'):
    """Draw image from array

    Args:
        data (array): numpy array
        out_file (str): output file name
        cmap (str, optional): colormap. Defaults to 'jet'.
    """
    plt.figure()
    ax = plt.gca()
    ax.axis('off')

    data[data==0.0] = np.nan
    im = ax.imshow(data, cmap=cmap)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    plt.colorbar(im, cax=cax)
    plt.savefig(out_file, bbox_inches='tight', dpi=200)


def merge_data(diff_dir, pairs, sub_swath, bursts_numbers, azi_poly, width, lines, unit, out_dir):
    """Merge overlap region into one file for phase and cc

    Args:
        diff_dir (str): diff directory
        pairs (list): ifg pairs
        sub_swath (list): sub_swath nubmer
        bursts_numbers (list): bursts numbers for per sub_swath
        azi_poly (array): azimuth overlap polygons
        width (int): width for merged file
        lines (int): lines for merged file
        unit (list): unit of transforming BOI displacement
        out_dir (str): directory of saving merged resluts
    """
    merged_phase = np.zeros((lines, width), dtype='float32')
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
                merged_phase[line_start:line_start + line_sm, range_start:range_start + width_sm] = phase * unit[index]

                cc_file = f"{pair}.iw{i}.{j}.diff20.adf.cc"
                cc = read_gamma(cc_file, line_sm, 'float32')
                merged_cc[line_start:line_start + line_sm, range_start:range_start + width_sm] = cc

        merged_phase_file = os.path.join(out_dir, f"{pair}.boi.disp")
        write_gamma(merged_phase, merged_phase_file, 'float32')
        draw_img(merged_phase, merged_phase_file + '.png')
        merged_phase = None

        merged_cc_file = os.path.join(out_dir, f"{pair}.boi.cc")
        write_gamma(merged_cc, merged_cc_file, 'float32')
        draw_img(merged_cc, merged_cc_file + '.png', cmap='gray_r')
        merged_cc = None

def create_off(rslc_dir, pairs, rlks, alks, extension, out_dir):
    """Create off par for MintPy

    Args:
        rslc_dir (str): RSLC directory
        pairs (list): ifg pairs
        rlks (int): range looks
        alks (int): azimuth looks
        extension (str): slc extension
        out_dir (str): output directory
    """
    os.chdir(out_dir)

    for pair in pairs:
        m_date = pair[0:8]
        s_date = pair[9:17]

        m_rslc = os.path.join(rslc_dir, m_date, m_date + '.' + extension)
        m_rslc_par = m_rslc + '.par'

        s_rslc = os.path.join(rslc_dir, s_date, s_date + '.' + extension)
        s_rslc_par = s_rslc + '.par'

        call_str = f'echo "{pair}\\n\\n\\n\\n\\n\\n\\n" > off_par.in'
        os.system(call_str)

        call_str = f"create_offset {m_rslc_par} {s_rslc_par} {pair}.off 1 {rlks} {alks} < off_par.in"
        os.system(call_str)

        del_file('off_par.in')


def check_dem_size(slc_par, dem_par):
    """Check whether the dem completely covers the research area

    Args:
        slc_par (str): SLC parameter file
        dem_par (str): dem parameter file
    """
    # calc lon and lat of dem
    width = int(read_gamma_par(dem_par, 'width'))
    nlines = int(read_gamma_par(dem_par, 'nlines'))
    corner_lat = float(read_gamma_par(dem_par, 'corner_lat'))
    corner_lon = float(read_gamma_par(dem_par, 'corner_lon'))
    post_lat = float(read_gamma_par(dem_par, 'post_lat'))
    post_lon = float(read_gamma_par(dem_par, 'post_lon'))

    north = corner_lat
    south = north + post_lat * nlines
    west = corner_lon
    east = west + post_lon * width

    print(f"longitude and latitude of DEM: {west:>9}{east:>9}{south:>9}{north:>9}")

    # calc lon and lat of slc
    out = Popen(f"SLC_corners {slc_par}", shell=True, stdout=PIPE)
    out_info = out.stdout.readlines()
    for line in out_info:
        line = str(line, 'utf-8')
        if "upper left corner" in line:
            sp = line.split()
            min_lon, max_lat = float(sp[-1]), float(sp[-2])
        if "lower right corner" in line:
            sp = line.split()
            max_lon, min_lat = float(sp[-1]), float(sp[-2])

    print(f"longitude and latitude of SLC: {min_lon:>9}{max_lon:>9}{min_lat:>9}{max_lat:>9}")

    # check
    e1 = (min_lon >= west)
    e2 = (max_lon <= east)
    e3 = (min_lat >= south)
    e4 = (max_lat <= north)

    if not (e1 and e2 and e3 and e4):
        print("SLC out of DEM coverage, please remake DEM!")
        sys.exit()


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

    # check dem
    sm_rslc_par = os.path.join(rslc_dir, ref_slc, f"{ref_slc}.{slc_extension}.par")
    check_dem_size(sm_rslc_par, dem_par)

    # multi-look
    mli_dir = os.path.join(out_dir, 'mli')
    mk_dir(mli_dir)

    slc_tab = os.path.join(mli_dir, 'slc_tab')
    mk_tab(rslc_dir, slc_tab, slc_extension)

    mli_all(slc_tab, mli_dir, rlks, alks)

    # make lookup table
    geo_dir = os.path.join(out_dir, 'geo')
    mk_dir(geo_dir)

    sm_mli = os.path.join(mli_dir, ref_slc + '.rmli')
    sm_mli_par = sm_mli + '.par'

    make_lookup(sm_mli, sm_mli_par, dem, dem_par, geo_dir)

    # select pairs
    base_dir = os.path.join(out_dir, 'base_calc')
    mk_dir(base_dir)

    pairs = select_pairs_sbas(slc_tab, sm_rslc_par, max_sb, max_tb, base_dir)

    # extract overlap area
    bo_dir = os.path.join(out_dir, 'burst_overlap')
    mk_dir(bo_dir)

    bursts_numbers = extract_burst_overlap(rslc_dir, sub_swath, pol, slc_extension, rlks, alks, bo_dir)

    # burst overlap double-difference interferometry
    ddi_dir = os.path.join(out_dir, 'ddi')
    mk_dir(ddi_dir)

    burst_overlap_ddi(bo_dir, pairs, sub_swath, pol, slc_extension, bursts_numbers, rlks, alks, ddi_dir)

    # merge data
    merge_dir = os.path.join(out_dir, 'merge')
    mk_dir(merge_dir)

    azi_poly = get_overlap_poly(rslc_dir, ref_slc, sub_swath, pol, slc_extension, rlks, alks, 1, merge_dir)
    width = int(read_gamma_par(sm_mli_par, 'range_samples'))
    lines = int(read_gamma_par(sm_mli_par, 'azimuth_lines'))
    unit = get_unit(rslc_dir, ref_slc, sub_swath, slc_extension, pol)

    merge_data(ddi_dir, pairs, sub_swath, bursts_numbers, azi_poly, width, lines, unit, merge_dir)

    # create off file
    create_off(rslc_dir, pairs, rlks, alks, slc_extension, merge_dir)

    print("\nAll done, enjoy it!\n")

if __name__ == "__main__":
    main()

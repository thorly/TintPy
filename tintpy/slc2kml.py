#!/usr/bin/env python3
################################
# Program is part of TintPy    #
# Copyright (c) 2021, Lei Yuan #
# Author: Lei Yuan, 2021       #
################################

import argparse
import glob
import os
import sys
import zipfile

EXAMPLE = """Example:
  python3 slc2kml.py /ly/slc/20201229 /ly/dem /ly/kml 28 7
  python3 slc2kml.py /ly/slc/20201229 /ly/dem /ly/kml 28 7 -e .slc
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(description="Make kml with dem and SLC using GAMMA.",
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)
    parser.add_argument('slc_dir', help='slc directory contains slc and slc_par.', type=str)
    parser.add_argument('dem_dir', help='dem directory contains *.dem and *.dem.par.', type=str)
    parser.add_argument('out_dir', help='directory for saving kml.', type=str)
    parser.add_argument('rlks', help='range looks.', type=int)
    parser.add_argument('alks', help='azimuth looks.', type=int)
    parser.add_argument('-e', dest='extension', help='file extension for slc (defaults: .rslc)', default='.rslc')

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


def main():
    # get inps
    inps = cmdline_parser()
    slc_dir = inps.slc_dir
    slc_dir = os.path.abspath(slc_dir)
    dem_dir = inps.dem_dir
    dem_dir = os.path.abspath(dem_dir)
    out_dir = inps.out_dir
    out_dir = os.path.abspath(out_dir)
    rlks = inps.rlks
    alks = inps.alks
    extension = inps.extension

    # check extension
    if not extension.startswith('.'):
        extension = '.' + extension

    # check slc_dir
    if not os.path.isdir(slc_dir):
        sys.exit(f"{slc_dir} does not exist.")

    slcs = glob.glob(slc_dir + '/[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]' + extension)
    slc_pars = [i + '.par' for i in slcs]
    for i, j in zip(slcs, slc_pars):
        if os.path.isfile(i) and os.path.isfile(j):
            pass
        else:
            sys.exit('Cannot find slc or slc_par in {}'.format(slc_dir))

    slc, slc_par = slcs[0], slc_pars[0]

    # get date
    slc_name = os.path.basename(slc)
    slc_date = slc_name[0:8]

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
                break
            else:
                sys.exit(f'Cannot find *.dem and *.dem.par in {dem_dir}.')

    # check out_dir
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    os.chdir(out_dir)

    # geocode slc
    m_mli = os.path.join(out_dir, f"{slc_date}.mli")
    m_mli_par = m_mli + '.par'

    call_str = f"multi_look {slc} {slc_par} {m_mli} {m_mli_par} {rlks} {alks}"
    os.system(call_str)

    utm_dem = os.path.join(out_dir, 'dem_seg')
    utm_dem_par = utm_dem + '.par'
    utm2rdc = os.path.join(out_dir, 'lookup_table')
    sim_sar_utm = os.path.join(out_dir, f"{slc_date}.sim_sar")
    u = os.path.join(out_dir, 'u')
    v = os.path.join(out_dir, 'v')
    inc = os.path.join(out_dir, 'inc')
    psi = os.path.join(out_dir, 'psi')
    pix = os.path.join(out_dir, 'pix')
    ls_map = os.path.join(out_dir, 'ls_map')

    call_str = f"gc_map {m_mli_par} - {dem_par} {dem} {utm_dem_par} {utm_dem} {utm2rdc} 1 1 {sim_sar_utm} {u} {v} {inc} {psi} {pix} {ls_map} 8 1"
    os.system(call_str)

    width_utm_dem = read_gamma_par(utm_dem_par, 'width')
    width_mli = read_gamma_par(m_mli_par, 'range_samples')
    geo_m_mli = os.path.join(out_dir, f"geo_{slc_date}.mli")

    call_str = f"geocode_back {m_mli} {width_mli} {utm2rdc} {geo_m_mli} {width_utm_dem} - 2 0"
    os.system(call_str)

    geo_m_mli_bmp = geo_m_mli + '.bmp'
    call_str = f"raspwr {geo_m_mli} {width_utm_dem} 1 0 1 1 1. .35 1 {geo_m_mli_bmp}"
    os.system(call_str)

    mli_kml = os.path.basename(geo_m_mli + '.kml')
    geo_m_mli_bmp = os.path.basename(geo_m_mli_bmp)

    call_str = f"kml_map {geo_m_mli_bmp} {utm_dem_par} {mli_kml}"
    os.system(call_str)

    # unzip bmp and kml
    kmz_name = os.path.basename(geo_m_mli) + '.kmz'
    with zipfile.ZipFile(kmz_name, 'w') as f:
        f.write(mli_kml)
        f.write(geo_m_mli_bmp)

    # delete files
    files = os.listdir(out_dir)
    for file in files:
        if file not in [kmz_name, mli_kml, geo_m_mli_bmp]:
            os.remove(file)


if __name__ == "__main__":
    main()

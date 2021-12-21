#!/usr/bin/env python3
################################
# Program is part of TintPy    #
# Copyright (c) 2021, Lei Yuan #
# Author: Lei Yuan, 2021       #
################################

import argparse
import os
import sys
import zipfile

EXAMPLE = """Example:
  python3 ph2kmz.py ph_rate lookup_fine dem_seg.par cc mli diff_par phase_stacking_res 6
  python3 ph2kmz.py ph_rate lookup_fine dem_seg.par cc mli diff_par phase_stacking_res 5 -t 0.3
  python3 ph2kmz.py ph_rate lookup_fine dem_seg.par cc mli diff_par phase_stacking_res 3 4 5 -t 0.3
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(
        description='Display wrapped phase derived by Stacking in Google Earth.',
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)

    parser.add_argument('ph_rate', type=str, help='radar coordinate phase rate.')
    parser.add_argument('lookup', type=str, help='lookup_fine.')
    parser.add_argument('dem_seg_par', type=str, help='dem_seg.par.')
    parser.add_argument('cc', type=str, help='radar coordinate coherence file.')
    parser.add_argument('pwr', type=str, help='radar coordinate pwr file.')
    parser.add_argument('diff_par', type=str, help='*.diff_par file.')
    parser.add_argument('out_dir', type=str, help='output directory.')
    parser.add_argument('cycles', type=float, nargs='+', help='data value per color cycle (6 means -3 ~ 3).')
    parser.add_argument('-t', dest='threshold', type=float, help='display coherence threshold (default: 0.0, display all).', default=0.0)

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


def geocode(infile, lookup_file, outfile, width_rdr, width_geo, lines_geo):
    """Geocoding of image data using lookup table values

    Args:
        infile (str): file
        lookup_file (str): lookup table
        outfile (str): geocoded file
        width_rdr (int): rdr width
        width_geo (int): geo width
        lines_geo (int): geo length
    """
    call_str = f"geocode_back {infile} {width_rdr} {lookup_file} {outfile} {width_geo} {lines_geo} 1 0"
    os.system(call_str)


def check_file(file):
    """Check if the file exists 

    Args:
        file (str): file
    """
    if not os.path.isfile(file):
        sys.exit("{} does not exist.".format(file))


def main():
    # get inps
    inps = cmdline_parser()
    ph_rate = os.path.abspath(inps.ph_rate)
    lookup = os.path.abspath(inps.lookup)
    dem_seg_par = os.path.abspath(inps.dem_seg_par)
    cc = os.path.abspath(inps.cc)
    pwr = os.path.abspath(inps.pwr)
    diff_par = os.path.abspath(inps.diff_par)
    out_dir = os.path.abspath(inps.out_dir)
    cycles = inps.cycles
    threshold = inps.threshold

    # check file
    check_file(ph_rate)
    check_file(lookup)
    check_file(cc)
    check_file(pwr)
    check_file(dem_seg_par)
    check_file(diff_par)

    # check out_dir
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # read gamma paramaters
    width_rdr = read_gamma_par(diff_par, 'range_samp_1')
    width = read_gamma_par(dem_seg_par, 'width')
    nlines = read_gamma_par(dem_seg_par, 'nlines')

    # geocoding
    geo_ph_rate = ph_rate + '.geo'
    geocode(ph_rate, lookup, geo_ph_rate, width_rdr, width, nlines)

    geo_pwr = pwr + '.geo'
    geocode(pwr, lookup, geo_pwr, width_rdr, width, nlines)

    geo_cc = cc + '.geo'
    geocode(cc, lookup, geo_cc, width_rdr, width, nlines)

    for cycle in cycles:
        # generate raster graphics image of phase + intensity data
        out_bmp = f"{geo_ph_rate}_{cycle}.bmp"
        call_str = f"rasdt_pwr24 {geo_ph_rate} {geo_pwr} {width} 1 1 0 1 1 {cycle} 1. .35 1 {out_bmp} {geo_cc} 1 {threshold}"
        os.system(call_str)

        # create kml XML file with link to image
        geo_ph_rate_name = os.path.basename(geo_ph_rate)
        kml_name = f"{geo_ph_rate_name}_{cycle}.kml"
        out_bmp_name = os.path.basename(out_bmp)

        os.chdir(out_dir)
        call_str = f"kml_map {out_bmp_name} {dem_seg_par} {kml_name}"
        os.system(call_str)

        # unzip bmp and kml
        kmz_name = f"{geo_ph_rate_name}_{cycle}.kmz"
        with zipfile.ZipFile(kmz_name, 'w') as f:
            f.write(kml_name)
            os.remove(kml_name)
            f.write(out_bmp_name)
            os.remove(out_bmp_name)

    print('\nAll done, enjoy it!\n')


if __name__ == "__main__":
    main()

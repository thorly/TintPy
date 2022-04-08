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

from subprocess import Popen, PIPE

EXAMPLE = """Example:
  python3 s1_coreg2.py /ly/slc /ly/rslc /ly/dem 20211229 2 -r 8 -a 2
  python3 s1_coreg2.py /ly/slc /ly/rslc /ly/dem 20211229 1 2 -r 8 -a 2
  python3 s1_coreg2.py /ly/slc /ly/rslc /ly/dem 20211229 1 2 3 -r 8 -a 2 -f t
"""


def cmd_line_parser():
    parser = argparse.ArgumentParser(description="Coregistrate all of the Sentinel-1 TOPS SLCs to a reference SLC.",
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)

    parser.add_argument("slc_dir", help="SLCs directory")
    parser.add_argument("rslc_dir", help="RSLCs directory")
    parser.add_argument("dem_dir", help="dem directory including *.dem and *.dem.par")
    parser.add_argument("ref_slc", help="reference SLC", type=str)
    parser.add_argument("sub_swath", type=int, nargs="+", choices=[1, 2, 3], help="sub_swath number for coregistration")
    parser.add_argument("-r", dest="rlks", help="range looks", type=int, required=True)
    parser.add_argument("-a", dest="alks", help="azimuth looks", type=int, required=True)
    parser.add_argument("-p", dest="pol", help="polarization(defaults: vv)", choices=["vv", "vh"], default="vv")
    parser.add_argument("-f", dest="flag", choices=["t", "T", "f", "F"], help="flag for saving iw*.rslc (t for YES, f for No, defaults: f)", default="f")
    inps = parser.parse_args()

    return inps


def run(call_str):
    """Run command

    Args:
        call_str (str): command
    """
    os.system(call_str)


def read_gamma_par(par_file, keyword):
    """Extract value from par_file using keyword

    Args:
        par_file (str): GAMMA parameter file
        keyword (str): keyword like "reange_sample"
    """
    value = None
    with open(par_file, "r", encoding="utf-8") as f:
        for line in f.readlines():
            if line.count(keyword) == 1:
                value = line.split()[1].strip()

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
    mli_par = mli + ".par"

    run(f"multi_look {slc} {slc_par} {mli} {mli_par} {rlks} {alks}")
    
    run(f"gc_map {mli_par} - {dem_par} {dem} dem_seg.par dem_seg lookup_table 1 1 sim_sar u v inc psi pix ls_map 8 1")
    
    run(f"pixel_area {mli_par} dem_seg.par dem_seg lookup_table ls_map inc pix_sigma0 pix_gamma0")
    

    width_mli = read_gamma_par(mli_par, "range_samples")

    run(f"raspwr pix_gamma0 {width_mli} - - - - - - - pix_gamma0.bmp")
    
    run(f"create_diff_par {mli_par} - {date}.diff_par 1 0")

    run(f"offset_pwrm pix_sigma0 {mli} {date}.diff_par offs snr 64 64 offsets 2 100 100 5.0")
    
    run(f"offset_fitm offs snr {date}.diff_par coffs coffsets 5.0 1")
    
    width_utm_dem = read_gamma_par("dem_seg.par", "width")

    run(f"gc_map_fine lookup_table {width_utm_dem} {date}.diff_par lookup_table_fine 1")
    
    mli_name = os.path.basename(mli)

    run(f"geocode_back {mli} {width_mli} lookup_table_fine {mli_name}.geo {width_utm_dem} - 2 0")
    
    run(f"raspwr {mli_name}.geo {width_utm_dem} 1 0 1 1 1. .35 1 {mli_name}.geo.bmp")
    
    length_mli = read_gamma_par(mli_par, "azimuth_lines")

    run(f"geocode lookup_table_fine dem_seg {width_utm_dem} {date}.dem {width_mli} {length_mli} 2 0")
    
    run(f"rashgt {date}.dem {mli} {width_mli} 1 1 0 1 1 160.0 1. .35 1 {date}.dem.bmp")
    
    rdc_dem = os.path.join(out_dir, f"{date}.dem")

    return rdc_dem


def write_tab(slcs, tab_file):
    """Write path of slc slc_par and tops_par to tab file

    Args:
        slcs (str list): slc files
        tab_file (str): output file
    """
    with open(tab_file, "w+", encoding="utf-8") as f:
        for slc in slcs:
            f.write(f"{slc} {slc}.par {slc}.tops_par\n")


def slc2bmp(slc, slc_par, rlks, alks, bmp):
    """Generate 8-bit raster graphics image of intensity of complex (SLC) data

    Args:
        slc (str): slc file
        slc_par (str): slc par file
        rlks (int): range looks
        alks (int): azimuth looks
        bmp (str): output bmp
    """
    width = read_gamma_par(slc_par, "range_samples")
    if width:
        run(f"rasSLC {slc} {width} 1 0 {rlks} {alks} 1. .35 1 0 0 {bmp}")


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

    slc_tab = os.path.join(date_slc_dir, "slc_tab")
    with open(slc_tab, "w+", encoding="utf-8") as f:
        for i in sub_swaths:
            slc = os.path.join(date_slc_dir, f"{date}.iw{i}.{pol}.slc")
            slc_par = slc + ".par"
            tops_par = slc + ".tops_par"
            f.write(f"{slc} {slc_par} {tops_par}\n")

    slc_out = os.path.join(out_dir, date + ".slc")
    slc_par_out = os.path.join(out_dir, date + ".slc.par")

    run(f"SLC_mosaic_S1_TOPS {slc_tab} {slc_out} {slc_par_out} {rlks} {alks}")

    os.remove(slc_tab)

    return slc_out, slc_par_out


def check_dem_size(slc_par, dem_par):
    """Check whether the dem completely covers the research area

    Args:
        slc_par (str): SLC parameter file
        dem_par (str): dem parameter file
    """
    # calc lon and lat of dem
    width = int(read_gamma_par(dem_par, "width"))
    nlines = int(read_gamma_par(dem_par, "nlines"))
    corner_lat = float(read_gamma_par(dem_par, "corner_lat"))
    corner_lon = float(read_gamma_par(dem_par, "corner_lon"))
    post_lat = float(read_gamma_par(dem_par, "post_lat"))
    post_lon = float(read_gamma_par(dem_par, "post_lon"))

    north = corner_lat
    south = north + post_lat * nlines
    west = corner_lon
    east = west + post_lon * width

    print(f"longitude and latitude of DEM: {west:>9}{east:>9}{south:>9}{north:>9}")

    # calc lon and lat of slc
    out = Popen(f"SLC_corners {slc_par}", shell=True, stdout=PIPE)
    out_info = out.stdout.readlines()
    for line in out_info:
        line = str(line, "utf-8")
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
                    date = os.path.basename(file)[0:8]
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
    flag = inps.flag.lower()

    # check slc_dir
    if not os.path.isdir(slc_dir):
        sys.exit("{} does not exist.".format(slc_dir))

    # check dem
    if not os.path.isdir(dem_dir):
        sys.exit("{} does not exist.".format(dem_dir))
    else:
        dems = glob.glob(dem_dir + "/*.dem")
        dem_pars = [i + ".par" for i in dems]
        for i, j in zip(dems, dem_pars):
            if os.path.isfile(i) and os.path.isfile(j):
                dem = dems[0]
                dem_par = dem + ".par"
                break
            else:
                sys.exit(f"Cannot find *.dem and *.dem.par in {dem_dir}.")

    # get all date
    dates = sorted([i for i in os.listdir(slc_dir) if re.findall(r"^\d{8}$", i)])
    if len(dates) < 2:
        sys.exit("No enough SLCs.")

    # check ref_slc
    if re.findall(r"^\d{8}$", ref_slc):
        if not ref_slc in dates:
            sys.exit("No slc for {}.".format(ref_slc))
    else:
        sys.exit("Error date for ref_slc.")

    # check sub_swath
    error_date = {}
    for i in sub_swath:
        error_date[i] = []
        for date in dates:
            iw_slc = os.path.join(slc_dir, date,  f"{date}.iw{i}.{pol}.slc")
            iw_slc_par = iw_slc + ".par"
            iw_slc_tops_par = iw_slc + ".tops_par"

            e1 = os.path.isfile(iw_slc)
            e2 = os.path.isfile(iw_slc_par)
            e3 = os.path.isfile(iw_slc_tops_par)

            if not (e1 and e2 and e3):
                error_date[i].append(date)

    if error_date[list(error_date.keys())[0]]:
        for key in error_date.keys():
            for date in error_date[key]:
                print(f"No slc or slc_par or tops_par for {date} sub_swath {key}")
        sys.exit("\nPlease check it.")

    # check rslc_dir
    if not os.path.isdir(rslc_dir):
        os.mkdir(rslc_dir)

    # get slave dates
    m_date = ref_slc
    s_dates = dates.copy()
    s_dates.remove(m_date)

    # mosaic iw slc for m_date
    geo_dir = os.path.join(rslc_dir, "geo")

    if os.path.isdir(geo_dir):
        shutil.rmtree(geo_dir)

    os.mkdir(geo_dir)

    m_slc_dir = os.path.join(slc_dir, m_date)

    m_slc, m_slc_par = slc_mosaic(m_slc_dir, sub_swath, pol, rlks, alks, geo_dir)

    # check dem
    check_dem_size(m_slc_par, dem_par)

    # copy reference iw slc to rslc_dir
    m_rslc_dir = os.path.join(rslc_dir, m_date)

    if not os.path.isdir(m_rslc_dir):
        os.mkdir(m_rslc_dir)

    # copy m_date slc to m_rslc_dir
    m_rslc = os.path.join(m_rslc_dir, m_date + ".rslc")
    m_rslc_par = m_rslc + ".par"
    if not os.path.isfile(m_rslc):
        shutil.copy(m_slc, m_rslc)
    if not os.path.isfile(m_rslc_par):
        shutil.copy(m_slc_par, m_rslc_par)

    if flag == "t":
        print("Copy reference iw slc to {}".format(m_rslc_dir))
        for i in sub_swath:
            iw_slc = os.path.join(m_slc_dir, f"{m_date}.iw{i}.{pol}.slc")
            iw_slc_par = iw_slc + ".par"
            iw_slc_tops_par = iw_slc + ".tops_par"

            iw_rslc = os.path.join(m_rslc_dir, f"{m_date}.iw{i}.{pol}.rslc")
            iw_rslc_par = iw_rslc + ".par"
            iw_rslc_tops_par = iw_rslc + ".tops_par"

            copy_file(iw_slc, iw_rslc)
            copy_file(iw_slc_par, iw_rslc_par)
            copy_file(iw_slc_tops_par, iw_rslc_tops_par)

    # make rdc dem
    rdc_dem = make_rdc_dem(m_slc, m_slc_par, dem, dem_par, rlks, alks, geo_dir)

    for s_date in s_dates:
        s_slc_dir = os.path.join(slc_dir, s_date)
        s_rslc_dir = os.path.join(rslc_dir, s_date)

        if not os.path.isdir(s_rslc_dir):
            os.mkdir(s_rslc_dir)

        os.chdir(s_rslc_dir)

        # write tab files
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

        write_tab(s_iw_slcs, "slc_tab_s")
        write_tab(m_iw_slcs, "slc_tab_m")
        write_tab(s_iw_rslcs, "rslc_tab")

        # mosaic slc
        s_slc, s_slc_par = slc_mosaic(s_slc_dir, sub_swath, pol, rlks, alks, s_rslc_dir)

        s_mli = f"{s_date}.mli"
        s_mli_par = f"{s_mli}.par"
        run(f"multi_look {s_slc} {s_slc_par} {s_mli} {s_mli_par} {rlks} {alks}")

        m_mli = os.path.join(geo_dir, f"{m_date}.mli")
        m_mli_par = f"{m_mli}.par"

        lt = f"{s_date}.lt"
        run(f"rdc_trans {m_mli_par} {rdc_dem} {s_mli_par} {lt}")

        s_rslc = s_date + ".rslc"
        s_rslc_par = s_rslc + ".par"
        run(f"SLC_interp_lt_S1_TOPS slc_tab_s {s_slc_par} slc_tab_m {m_slc_par} {lt} {m_mli_par} {s_mli_par} - rslc_tab {s_rslc} {s_rslc_par}")

        pair = f"{m_date}-{s_date}"
        off_file = pair + ".off"
        run(f"create_offset {m_slc_par} {s_slc_par} {off_file} 1 {rlks} {alks} 0")

        offs_file = pair + ".offs"
        snr_file = pair + ".snr"
        run(f"offset_pwr {m_slc} {s_rslc} {m_slc_par} {s_rslc_par} {off_file} {offs_file} {snr_file} 256 64 - 1 100 100 15.0 4 0 0")

        run(f"offset_fit {offs_file} {snr_file} {off_file}  - - 0.5 1 0")

        run(f"SLC_interp_lt_S1_TOPS slc_tab_s {s_slc_par} slc_tab_m {m_slc_par} {lt} {m_mli_par} {s_mli_par} {off_file} rslc_tab {s_rslc} {s_rslc_par}")

        off1_file = f"{pair}.off1"
        run(f"create_offset {m_slc_par} {s_slc_par} {off1_file} 1 {rlks} {alks} 0")

        offs1_file = f"{pair}.offs1"
        snr1_file = f"{pair}.snr1"
        coreg_file = f"{pair}.coreg"
        run(f"offset_pwr {m_slc} {s_rslc} {m_slc_par} {s_rslc_par} {off1_file} {offs1_file} {snr1_file} 256 64 {coreg_file} 1 200 200 16.0 4 0 0")

        run(f"offset_fit {offs1_file} {snr1_file} {off1_file}  - - 0.5 1 0")

        off_total_file = f"{pair}.off.total"
        run(f"offset_add {off_file} {off1_file} {off_total_file}")

        run(f"SLC_interp_lt_S1_TOPS slc_tab_s {s_slc_par} slc_tab_m {m_slc_par} {lt} {m_mli_par} {s_mli_par} {off_total_file} rslc_tab {s_rslc} {s_rslc_par}")

        off_corrected1 = f"{pair}.off.corrected1"
        run(f"S1_coreg_overlap slc_tab_m rslc_tab {pair} {off_total_file} {off_corrected1} > {s_date}.coreg_quality")
        run(f"SLC_interp_lt_S1_TOPS slc_tab_s {s_slc_par} slc_tab_m {m_slc_par} {lt} {m_mli_par} {s_mli_par} {off_corrected1} rslc_tab {s_rslc} {s_rslc_par}")

        off_corrected2 = f"{pair}.off.corrected2"
        run(f"S1_coreg_overlap slc_tab_m rslc_tab {pair} {off_corrected1} {off_corrected2} >> {s_date}.coreg_quality")
        run(f"SLC_interp_lt_S1_TOPS slc_tab_s {s_slc_par} slc_tab_m {m_slc_par} {lt} {m_mli_par} {s_mli_par} {off_corrected2} rslc_tab {s_rslc} {s_rslc_par}")

        off_corrected3 = f"{pair}.off.corrected3"
        run(f"S1_coreg_overlap slc_tab_m rslc_tab {pair} {off_corrected2} {off_corrected3} >> {s_date}.coreg_quality")
        run(f"SLC_interp_lt_S1_TOPS slc_tab_s {s_slc_par} slc_tab_m {m_slc_par} {lt} {m_mli_par} {s_mli_par} {off_corrected3} rslc_tab {s_rslc} {s_rslc_par}")

        off_corrected4 = f"{pair}.off.corrected4"
        run(f"S1_coreg_overlap slc_tab_m rslc_tab {pair} {off_corrected3} {off_corrected4} >> {s_date}.coreg_quality")
        run(f"SLC_interp_lt_S1_TOPS slc_tab_s {s_slc_par} slc_tab_m {m_slc_par} {lt} {m_mli_par} {s_mli_par} {off_corrected4} rslc_tab {s_rslc} {s_rslc_par}")

        run(f"phase_sim_orb {m_slc_par} {s_slc_par} {off_corrected4} {rdc_dem} {pair}.sim_unw {m_slc_par} - - 1 1")
        run(f"SLC_diff_intf {m_slc} {s_rslc} {m_slc_par} {s_rslc_par} {off_corrected4} {pair}.sim_unw {pair}.diff {rlks} {alks} 0 0 0.25 1 1")

        width_mli = read_gamma_par(m_mli_par, "range_samples")
        run(f"rasmph {pair}.diff {width_mli} 1 0 1 1 1. .35 1 {pair}.diff.bmp")
        run(f"rasmph_pwr {pair}.diff {m_mli} {width_mli} 1 1 0 1 1 1. .35 1 {pair}.diff.pwr.bmp")

        # clean s_rslc dir
        save_files = []
        save_files.append(s_date + ".rslc")
        save_files.append(s_date + ".rslc.par")
        save_files.append(pair + ".diff.bmp")
        save_files.append(pair + ".diff.pwr.bmp")
        save_files.append(s_date + '.coreg_quality')

        # save iw rslc
        if flag == "t":
            for i in sub_swath:
                iw_rslc = f"{s_date}.iw{i}.{pol}.rslc"
                iw_rslc_par = iw_rslc + ".par"
                iw_rslc_tops_par = iw_rslc + ".tops_par"
                save_files.append(iw_rslc)
                save_files.append(iw_rslc_par)
                save_files.append(iw_rslc_tops_par)

        # delete files
        for file in os.listdir(s_rslc_dir):
            if file not in save_files:
                os.remove(file)

    # del geo_dir
    if os.path.isdir(geo_dir):
        shutil.rmtree(geo_dir)

    # generate bmp for rslc
    rslc_files = glob.glob(rslc_dir + "/*/*.rslc")
    for rslc in rslc_files:
        slc2bmp(rslc, rslc + ".par", rlks, alks, rslc + ".bmp")

    # check coreg_quality
    quality_files = glob.glob(rslc_dir + '/*/*.coreg_quality')
    print_coreg_quality(quality_files)

    print('\ndaz10000 is a relative indicator which cannot guarantee the coregistration is ok, you should also check the diff files!')

    print("\nAll done, enjoy it!\n")


if __name__ == "__main__":
    main()

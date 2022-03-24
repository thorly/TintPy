#!/usr/bin/env python3
###############################################################################
# Generate differential interferogram and unwrap them from RSLCs using GAMMA  #
# APS correction using MATLAB (Only support .ztd file)                        #
# Copyright (c) 2022, Lei Yuan                                                #
###############################################################################

import argparse
import glob
import os
import re
import sys

dinsar_script = """#!/bin/bash
m_date=m_date_flag
s_date=s_date_flag
rslc_dir=rslc_dir_flag
m_rslc=$rslc_dir/$m_date/$m_date.rslc
m_par=$m_rslc.par
s_rslc=$rslc_dir/$s_date/$s_date.rslc
s_par=$s_rslc.par

dem=dem_flag
dem_par=dem_par_flag

m_mli=m_mli_flag
dem_rdc=dem_rdc_flag
diff_par=diff_par_flag

rlks=rlks_flag
alks=alks_flag

cc_thres=cc_thres_flag
roff=roff_flag
loff=loff_flag

pair=$m_date\_$s_date
off_par=$pair.off

echo -ne "$pair\\n 0 0\\n 32 32\\n 64 64\\n 7.0\\n 0\\n\\n" > off_par.in
create_offset $m_par $s_par $off_par 1 1 1 < off_par.in
rm -f off_par.in

init_offset_orbit $m_par $s_par $off_par
init_offset $m_rslc $s_rslc $m_par $s_par $off_par 1 1

SLC_intf $m_rslc $s_rslc $m_par $s_par $off_par $pair.int $rlks $alks - - 1 1

width_rdc=$(awk '$1 == "interferogram_width:" {print $2}' $off_par)
line_rdc=$(awk '$1 == "interferogram_azimuth_lines:" {print $2}' $off_par)

rasmph_pwr $pair.int $m_mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.int_pwr.bmp

base_init $m_par $s_par $off_par $pair.int $pair.base 0 1024 1024
base_perp $pair.base $m_par $off_par > $pair.base.perp

cc_wave $pair.int $m_par $s_par $pair.cor $width_rdc - - 3

phase_sim $m_par $off_par $pair.base $dem_rdc $pair.sim_unw 0 0 - -

sub_phase $pair.int $pair.sim_unw $diff_par $pair.diff 1 0
rasmph $pair.diff $width_rdc 1 0 1 1 1. 0.35 1 $pair.diff.bmp
rasmph_pwr $pair.diff $m_mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.diff_pwr.bmp

adf $pair.diff $pair.adf.diff1 $pair.adf.cc1 $width_rdc 0.3 128
adf $pair.adf.diff1 $pair.adf.diff2 $pair.adf.cc2 $width_rdc 0.3 64
adf $pair.adf.diff2 $pair.adf.diff $pair.adf.cc $width_rdc 0.3

rm -rf $pair.adf.diff1 $pair.adf.diff2 $pair.adf.cc1 $pair.adf.cc2

rasmph $pair.adf.diff $width_rdc 1 0 1 1 1. 0.35 1 $pair.adf.diff.bmp
rasmph_pwr $pair.adf.diff $m_mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.adf.diff_pwr.bmp


# rascc_mask $pair.adf.cc $m_mli $width_rdc 1 1 0 1 1 $cc_thres 0. .1 .9 1. .35 1 $pair.adf.cc_mask.bmp
rascc_mask $pair.cor $m_mli $width_rdc 1 1 0 1 1 $cc_thres 0. .1 .9 1. .35 1 $pair.adf.cc_mask.bmp

mcf $pair.adf.diff $pair.cor $pair.adf.cc_mask.bmp $pair.adf.unw $width_rdc 1 0 0 - - 1 1 - $roff $loff 0

rasrmg $pair.adf.unw $m_mli $width_rdc 1 1 0 1 1 .6 1. .35 .0 1 $pair.adf.unw_pwr.bmp $pair.adf.cc 1 .2
rasrmg $pair.adf.unw - $width_rdc 1 1 0 1 1 .5 1. .35 .0 1 $pair.adf.unw.bmp $pair.adf.cc 1 .2


quad_fit $pair.adf.unw $diff_par 32 32 $pair.adf.cc_mask.bmp $pair.plot 3
quad_sub $pair.adf.unw $diff_par $pair.adf.unw.sub 0 0

rasrmg $pair.adf.unw.sub $m_mli $width_rdc 1 1 0 1 1 .6 1. .35 .0 1 $pair.adf.unw.sub_pwr.bmp $pair.adf.cc 1 .2
rasrmg $pair.adf.unw.sub -  $width_rdc 1 1 0 1 1 .5 1. .35 .0 1 $pair.adf.unw.sub.bmp $pair.adf.cc 1 .2
"""

matlab_code="""function [] = process_gacos(m_gacos, s_gacos, dem_seg_par, wavelength, incidence, out_file)
%PROCESS_GACOS
% read dem_seg.par
fid=fopen(dem_seg_par,'r');
parms = textscan(fid,'%s','headerlines', 6);
fclose(fid);
line=str2double(parms{1,1}{4,1});
sample=str2double(parms{1,1}{2,1});
left_lon=str2double(parms{1,1}{10,1});
left_lat=str2double(parms{1,1}{6,1});
cha=str2double(parms{1,1}{18,1});
% generate grid
A=left_lon:cha:left_lon+(sample-1)*cha;
B=left_lat:-cha:left_lat-(line-1)*cha;
[lon_dem,lat_dem]=meshgrid(A,B);
% read rsc
rscfile=[m_gacos,'.rsc'];
fid=fopen(rscfile,'r');
parms=textscan(fid, '%s');
fclose(fid);
wid=str2double(parms{1,1}{2,1});
len=str2double(parms{1,1}{4,1});
left_lon=str2double(parms{1,1}{14,1});
left_lat=str2double(parms{1,1}{16,1});
xstep=str2double(parms{1,1}{18,1});
cha=xstep;
A=left_lon:cha:left_lon+(wid-1)*cha;
B=left_lat:-cha:left_lat-(len-1)*cha;
[lon_aps,lat_aps]=meshgrid(A,B);
% read gacos
if(max(max(lon_dem))>max(max(lon_aps)) || max(max(lat_dem))>max(max(lat_aps)) || min(min(lon_dem))<min(min(lat_aps)) || min(min(lat_dem))<min(min(lat_aps)))
    disp('Error: the gacos image is smaller than the dem!');
    return;
else
    % read m_gacos
    fid=fopen(m_gacos,'rb');
    [m_ztd, ~]=fread(fid,[wid, len],'float');
    fclose(fid);
    % read s_gacos
    fid=fopen(s_gacos,'rb');
    [s_ztd, ~]=fread(fid,[wid, len],'float');
    fclose(fid);
    % diff and convert
    diff_ztd = m_ztd - s_ztd;
    diff_ztd=rot90(flip(diff_ztd),3);
    diff_ztd_phs=-4*pi*diff_ztd./wavelength.*cosd(incidence);
    % interpolate
    inter_diff_ztd_phs=interp2(lon_aps,lat_aps,diff_ztd_phs,lon_dem,lat_dem,'cubic');
    if exist(out_file, 'file')
        delete(out_file)
    end
    fid = fopen(out_file,'a','b');
    inter_diff_ztd_phs = inter_diff_ztd_phs.';
    fwrite(fid,inter_diff_ztd_phs,'float32');
    fclose(fid);
    clear lon_aps;clear lat_aps;clear lon_dem;clear lat_dem;
    clear m_ztd;clear s_ztd;clear diff_ztd;clear diff_ztd_phs;clear inter_diff_ztd_phs;
end
end
"""

aps_correction_script="""#!/bin/bash
m_date=m_date_flag
s_date=s_date_flag

dem_seg_par=dem_seg_par_flag
lookup=lookup_flag
m_mli=m_mli_flag

cc_thres=cc_thres_flag
roff=roff_flag
loff=loff_flag
wavelength=wavelength_flag

off_par=off_par_flag
gacos_dir=gacos_dir_flag
diff_par=diff_par_flag

pair=$m_date\_$s_date

inc_angle=$(awk '$1 == "incidence_angle:" {print $2}' $m_mli.par)
m_gacos=$gacos_dir/$m_date.ztd
s_gacos=$gacos_dir/$s_date.ztd
matlab -nodesktop -nosplash -r "process_gacos('$m_gacos','$s_gacos','$dem_seg_par',$wavelength,$inc_angle,'$pair.gacos');quit;"

width_rdc=$(awk '$1 == "interferogram_width:" {print $2}' $off_par)
line_rdc=$(awk '$1 == "interferogram_azimuth_lines:" {print $2}' $off_par)
width_geo=$(awk '$1 == "width:" {print $2}' $dem_seg_par)

geocode $lookup $pair.gacos $width_geo $pair.gacos.rdc $width_rdc $line_rdc 1 0
raspwr $pair.gacos.rdc $width_rdc 1 0 1 1 1. .35 1 $pair.gacos.rdc.bmp

sub_phase $pair.diff $pair.gacos.rdc $diff_par $pair.diff.gacos 1 0

rasmph_pwr $pair.diff.gacos $m_mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.diff.gacos_pwr.bmp

adf $pair.diff.gacos $pair.adf.diff.gacos1 $pair.adf.cc.gacos1 $width_rdc 0.3 128
adf $pair.adf.diff.gacos1 $pair.adf.diff.gacos2 $pair.adf.cc.gacos2 $width_rdc 0.3 64
adf $pair.adf.diff.gacos2 $pair.adf.diff.gacos $pair.adf.cc.gacos $width_rdc 0.3

rm -rf $pair.adf.diff.gacos1 $pair.adf.diff.gacos2 $pair.adf.cc.gacos1 $pair.adf.cc.gacos2

rasmph $pair.adf.diff.gacos $width_rdc 1 0 1 1 1. 0.35 1 $pair.adf.diff.gacos.bmp
rasmph_pwr $pair.adf.diff.gacos $m_mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.adf.diff.gacos_pwr.bmp

rascc_mask $pair.adf.cc.gacos $m_mli $width_rdc 1 1 0 1 1 $cc_thres 0. .1 .9 1. .35 1 $pair.adf.cc.gacos_mask.bmp

mcf $pair.adf.diff.gacos $pair.adf.cc.gacos $pair.adf.cc.gacos_mask.bmp $pair.adf.unw.gacos $width_rdc 1 0 0 - - 1 1 - $roff $loff 0

rasrmg $pair.adf.unw.gacos $m_mli $width_rdc 1 1 0 1 1 .6 1. .35 .0 1 $pair.adf.unw.gacos_pwr.bmp $pair.adf.cc.gacos 1 .2
rasrmg $pair.adf.unw.gacos  - $width_rdc 1 1 0 1 1 .5 1. .35 .0 1 $pair.adf.unw.gacos.bmp $pair.adf.cc.gacos 1 .2

quad_fit $pair.adf.unw.gacos $diff_par 32 32 $pair.adf.cc.gacos_mask.bmp $pair.plot 3
quad_sub $pair.adf.unw.gacos $diff_par $pair.adf.unw.sub.gacos 0 0

rasrmg $pair.adf.unw.sub.gacos $m_mli $width_rdc 1 1 0 1 1 .6 1. .35 .0 1 $pair.adf.unw.sub.gacos_pwr.bmp $pair.adf.cc.gacos 1 .2
rasrmg $pair.adf.unw.sub.gacos -  $width_rdc 1 1 0 1 1 .5 1. .35 .0 1 $pair.adf.unw.sub.gacos.bmp $pair.adf.cc.gacos 1 .2
"""

USAGE = """Example:
  # reference point of unwrapping is (0, 0) and no mask
  python3 diff_by_baseline2.py /ly/rslc /ly/SBAS /ly/dem 20211229 8 2 200 60
  # set reference point of unwrapping and mask
  python3 diff_by_baseline2.py /ly/rslc /ly/SBAS /ly/dem 20211229 8 2 200 60 -r 100 -l 100 -c 0.3
  # for aps correction (Sentinel-1 data)
  python3 diff_by_baseline2.py /ly/rslc /ly/SBAS /ly/dem 20211229 8 2 200 60 -g /ly/gacos -w 0.05546
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(
        description='Generate differential interferogram and unwrap them from RSLCs using GAMMA.\n' + 
        'APS correction using MATLAB (Only support .ztd file)',
        formatter_class=argparse.RawTextHelpFormatter, epilog=USAGE)

    parser.add_argument('rslc_dir', help='RSLCs directory')
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('dem_dir', help='DEM directory contains *.dem and *.dem.par')
    parser.add_argument('ref_rslc', help='date of reference RSLC for calculating baseline')
    parser.add_argument('rlks', help='range looks', type=int)
    parser.add_argument('alks', help='azimuth looks', type=int)
    parser.add_argument('max_sb', type=float, help='maximum spatial baseline')
    parser.add_argument('max_tb', type=float, help='maximum temporal baseline')
    parser.add_argument('-e', dest='rslc_extension', type=str, default='.rslc', help='file extension for RSLCs (defaults: .rslc)')
    parser.add_argument('-g', dest='gacos_dir', help='directory contains GACOS files for aps correction')
    parser.add_argument('-w', dest='wavelength', type=float, help='Microwave length (Sentinel-1: 0.05546576, ALOS: 0.23830879)')
    parser.add_argument('-r', dest='roff', help='phase reference range offset to unwrap (defaults: 0)', type=float, default=0)
    parser.add_argument('-l', dest='loff', help='phase reference azimuth offset to unwrap (defaults: 0)', type=float, default=0)
    parser.add_argument('-c', dest='cc_thres', type=float, default=0, help='threshold of correlation for creating the unwrapping mask (0.0 --> 1.0) (defaults: 0)')

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
            slc = os.path.join(slc_dir, date, date + slc_extension)
            slc_par = slc + '.par'
            f.write(slc + '    ' + slc_par + '\n')


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


def mk_mli_all(slc_tab, out_mli_dir, rlks, alks):
    """Calculate MLI images for a stack of SLCs

    Args:
        slc_tab (str): slc tab file including slc slc_par
        out_dir (str): output directory
        rlks (int): range looks
        alks (int): azimuth looks
    """
    cmd_str = f"mk_mli_all {slc_tab} {out_mli_dir} {rlks} {alks}"
    os.system(cmd_str)


def mk_geo_all(mli_dir, dates, dem, dem_par, out_geo_dir):
    """Terrain geocoding of SAR images with lookup table refinement and resample DEM to SAR Range-Doppler Coordinates (RDC)

    Args:
        mli_dir (str): mli directory
        dates (list): dates for makeing geo
        dem (str): dem file
        dem_par (str): dem parameter file
        out_geo_dir (str): output directory
    """
    os.chdir(out_geo_dir)
    post = read_gamma_par(dem_par, 'post_lon')
    for date in dates:
        mli = glob.glob(os.path.join(mli_dir, date + '*mli'))[0]
        mli_par = mli + '.par'
        dem_seg = date + '.dem_seg'
        dem_seg_par = date + '.dem_seg.par'
        for i in range(5):
            cmd_str = f"mk_geo {mli} {mli_par} {dem} {dem_par} {dem_seg} {dem_seg_par} {out_geo_dir} {date} {post} {i} 2"
            os.system(cmd_str)


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


def check_gacos(dates, gacos_dir, file_extension):
    """Check gacos files

    Args:
        dates (list): RSLC dates
        gacos_dir (str): gacos directory
        file_extension (str): file extension for gacos
    """
    not_exist = []
    if not os.path.isdir(gacos_dir):
        sys.exit('{} does not exist.'.format(gacos_dir))
    for date in dates:
        file = os.path.join(gacos_dir, date + file_extension)
        if not os.path.isfile(file):
            not_exist.append(file)
    if not_exist:
        for i in not_exist:
            print('{} does not exist.'.format(i))
        sys.exit()


if __name__ == "__main__":
    # get inputs
    inps = cmdline_parser()
    rslc_dir = os.path.abspath(inps.rslc_dir)
    out_dir = os.path.abspath(inps.out_dir)
    dem_dir = os.path.abspath(inps.dem_dir)
    ref_rslc = inps.ref_rslc
    rlks = inps.rlks
    alks = inps.alks
    max_sb = inps.max_sb
    max_tb = inps.max_tb

    rslc_extension = inps.rslc_extension

    gacos_dir = inps.gacos_dir
    wavelength = inps.wavelength

    roff = inps.roff
    loff = inps.loff
    cc_thres = inps.cc_thres

    # check rslc_dir
    if not os.path.isdir(rslc_dir):
        sys.exit('{} does not exist.'.format(rslc_dir))
    dates = sorted([i for i in os.listdir(rslc_dir) if re.match(r'^\d{8}$', i)])
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
    if re.findall(r'^\d{8}$', ref_rslc):
        if not ref_rslc in dates:
            sys.exit('No RSLC for {}.'.format(ref_rslc))
    else:
        sys.exit('Error date for ref_rslc.')

    # check RSLC extension
    if not rslc_extension.startswith('.'):
        rslc_extension = '.' + rslc_extension
    for date in dates:
        rslc_path = os.path.join(rslc_dir, date, date + rslc_extension)
        if not os.path.isfile(rslc_path):
            sys.exit('Error file extension for RSLC.')

    # check gacos_dir
    if gacos_dir and wavelength:
        gacos_dir = os.path.abspath(gacos_dir)
        if not os.path.isdir(gacos_dir):
            sys.exit('{} does not exist.'.format(gacos_dir))

        check_gacos(dates, gacos_dir, '.ztd')

    # multi-look
    mli_dir = os.path.join(out_dir, 'mli')
    if not os.path.isdir(mli_dir):
        os.mkdir(mli_dir)

    rslc_tab = os.path.join(mli_dir, 'rslc_tab')
    mk_tab(rslc_dir, rslc_tab, rslc_extension)
    mk_mli_all(rslc_tab, mli_dir, rlks, alks)

    # mk_geo
    geo_dir = os.path.join(out_dir, 'geo')
    if not os.path.isdir(geo_dir):
        os.mkdir(geo_dir)

    mk_geo_all(mli_dir, dates[0:-1], dem, dem_par, geo_dir)

    # select pairs
    base_dir = os.path.join(out_dir, 'base_calc')
    if not os.path.isdir(base_dir):
        os.mkdir(base_dir)

    m_rslc_par = os.path.join(rslc_dir, ref_rslc, ref_rslc + rslc_extension + '.par')
    pairs = select_pairs_sbas(rslc_tab, m_rslc_par, max_sb, max_tb, base_dir)

    # diff and unwrap
    diff_dir = os.path.join(out_dir, 'diff')
    if not os.path.isdir(diff_dir):
        os.mkdir(diff_dir)

    os.chdir(diff_dir)

    for pair in pairs:
        m_date = pair[0:8]
        s_date = pair[9:17]

        m_mli = glob.glob(os.path.join(mli_dir, m_date + '*mli'))[0]
        dem_rdc = os.path.join(geo_dir, m_date + '_dem.rdc')
        diff_par = os.path.join(geo_dir, m_date + '.diff_par')

        # replace value
        dinsar_script = dinsar_script.replace('m_date_flag', m_date)
        dinsar_script = dinsar_script.replace('s_date_flag', s_date)
        dinsar_script = dinsar_script.replace('rslc_dir_flag', rslc_dir)
        dinsar_script = dinsar_script.replace('dem_flag', dem)
        dinsar_script = dinsar_script.replace('dem_par_flag', dem_par)
        dinsar_script = dinsar_script.replace('m_mli_flag', m_mli)
        dinsar_script = dinsar_script.replace('dem_rdc_flag', dem_rdc)
        dinsar_script = dinsar_script.replace('diff_par_flag', diff_par)
        dinsar_script = dinsar_script.replace('rlks_flag', str(rlks))
        dinsar_script = dinsar_script.replace('alks_flag', str(alks))
        dinsar_script = dinsar_script.replace('cc_thres_flag', str(cc_thres))
        dinsar_script = dinsar_script.replace('roff_flag', str(roff))
        dinsar_script = dinsar_script.replace('loff_flag', str(loff))

        # write  script
        run_script = pair + '_DInSAR.sh'
        with open(run_script, 'w+') as f:
            f.write(dinsar_script)

        # run script
        call_str = 'bash ' + run_script
        os.system(call_str)

    # aps correction
    if gacos_dir and wavelength:
        gacos_dir = os.path.abspath(gacos_dir)

        os.chdir(diff_dir)

        for pair in pairs:
            m_date = pair[0:8]
            s_date = pair[9:17]
            # replace value
            dem_seg_par = os.path.join(geo_dir, m_date + '.dem_seg.par')
            lookup = os.path.join(geo_dir, m_date + '_1.map_to_rdc')
            m_mli = glob.glob(os.path.join(mli_dir, m_date + '*mli'))[0]
            off_par = os.path.join(diff_dir, pair + '.off')
            diff_par = os.path.join(geo_dir, m_date + '.diff_par')

            aps_correction_script = aps_correction_script.replace('m_date_flag', m_date)
            aps_correction_script = aps_correction_script.replace('s_date_flag', s_date)
            aps_correction_script = aps_correction_script.replace('dem_seg_par_flag', dem_seg_par)
            aps_correction_script = aps_correction_script.replace('lookup_flag', lookup)
            aps_correction_script = aps_correction_script.replace('m_mli_flag', m_mli)
            aps_correction_script = aps_correction_script.replace('cc_thres_flag', str(cc_thres))
            aps_correction_script = aps_correction_script.replace('roff_flag', str(roff))
            aps_correction_script = aps_correction_script.replace('loff_flag', str(loff))
            aps_correction_script = aps_correction_script.replace('off_par_flag', off_par)
            aps_correction_script = aps_correction_script.replace('wavelength_flag', str(wavelength))
            aps_correction_script = aps_correction_script.replace('gacos_dir_flag', gacos_dir)
            aps_correction_script = aps_correction_script.replace('diff_par_flag', diff_par)

            # write script
            with open('process_gacos.m', 'w+') as f:
                f.write(matlab_code)

            run_script = pair + '_aps_correction.sh'
            with open(run_script, 'w+') as f:
                f.write(aps_correction_script)

            # run script
            call_str = 'bash ' + run_script
            os.system(call_str)


    print('\nAll done, enjoy it!\n')

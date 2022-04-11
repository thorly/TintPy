#!/usr/bin/env python3
########################################################
# Run two-pass D-InSAR using GAMMA                     #
# APS correction using MATLAB (Only support .ztd file) #
# Copyright (c) 2022, Lei Yuan                         #
########################################################

import argparse
import glob
import os
import re
import sys

dinsar_script = """#!/bin/bash
m_rslc=m_rslc_flag
m_par=$m_rslc.par

s_rslc=s_rslc_flag
s_par=$s_rslc.par

dem=dem_flag
dem_par=dem_par_flag

rlks=rlks_flag
alks=alks_flag

cc_thres=cc_thres_flag
roff=roff_flag
loff=loff_flag

m_date=m_date_flag
s_date=s_date_flag
pair=$m_date\_$s_date


echo -ne "$pair\\n 0 0\\n 32 32\\n 64 64\\n 7.0\\n 0\\n\\n" > off_par.in
create_offset $m_par $s_par $pair.off 1 1 1 < off_par.in
rm -f off_par.in

init_offset_orbit $m_par $s_par $pair.off
init_offset $m_rslc $s_rslc $m_par $s_par $pair.off 1 1

SLC_intf $m_rslc $s_rslc $m_par $s_par $pair.off $pair.int $rlks $alks - - 1 1

width_rdc=$(awk '$1 == "interferogram_width:" {print $2}' $pair.off)
line_rdc=$(awk '$1 == "interferogram_azimuth_lines:" {print $2}' $pair.off)

multi_look $m_rslc $m_par $m_date.mli $m_date.mli.par $rlks $alks
multi_look $s_rslc $s_par $s_date.mli $s_date.mli.par $rlks $alks

raspwr $m_date.mli $width_rdc 1 0 1 1 1. .35 1 $m_date.mli.bmp
raspwr $s_date.mli $width_rdc 1 0 1 1 1. .35 1 $s_date.mli.bmp

rasmph_pwr $pair.int $m_date.mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.int_pwr.bmp

base_init $m_par $s_par $pair.off $pair.int $pair.base 0 1024 1024
base_perp $pair.base $m_par $pair.off > $pair.base.perp

cc_wave $pair.int $m_date.mli $s_date.mli $pair.cor $width_rdc - - 3

post=$(awk '$1 == "post_lon:" {print $2}' $dem_par)
mk_geo $m_date.mli $m_date.mli.par $dem $dem_par dem_seg dem_seg.par . $m_date $post 0 2
mk_geo $m_date.mli $m_date.mli.par $dem $dem_par dem_seg dem_seg.par . $m_date $post 2 2
mk_geo $m_date.mli $m_date.mli.par $dem $dem_par dem_seg dem_seg.par . $m_date $post 3 2
mk_geo $m_date.mli $m_date.mli.par $dem $dem_par dem_seg dem_seg.par . $m_date $post 4 2

mv ${m_date}_1.map_to_rdc $m_date.lookup_fine

phase_sim $m_par $pair.off $pair.base ${m_date}_dem.rdc $pair.sim_unw 0 0 - -

sub_phase $pair.int $pair.sim_unw $m_date.diff_par $pair.diff 1 0
rasmph $pair.diff $width_rdc 1 0 1 1 1. 0.35 1 $pair.diff.bmp
rasmph_pwr $pair.diff $m_date.mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.diff_pwr.bmp

adf $pair.diff $pair.adf.diff1 $pair.adf.cc1 $width_rdc 0.3 128
adf $pair.adf.diff1 $pair.adf.diff2 $pair.adf.cc2 $width_rdc 0.3 64
adf $pair.adf.diff2 $pair.adf.diff $pair.adf.cc $width_rdc 0.3

rm -rf $pair.adf.diff1 $pair.adf.diff2 $pair.adf.cc1 $pair.adf.cc2

rasmph $pair.adf.diff $width_rdc 1 0 1 1 1. 0.35 1 $pair.adf.diff.bmp
rasmph_pwr $pair.adf.diff $m_date.mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.adf.diff_pwr.bmp

# rascc_mask $pair.adf.cc $m_date.mli $width_rdc 1 1 0 1 1 $cc_thres 0. .1 .9 1. .35 1 $pair.adf.cc_mask.bmp
rascc_mask $pair.cor $m_date.mli $width_rdc 1 1 0 1 1 $cc_thres 0. .1 .9 1. .35 1 $pair.adf.cc_mask.bmp

mcf $pair.adf.diff $pair.cor $pair.adf.cc_mask.bmp $pair.adf.unw $width_rdc 1 0 0 - - 1 1 - $roff $loff 0

rasrmg $pair.adf.unw $m_date.mli $width_rdc 1 1 0 1 1 .6 1. .35 .0 1 $pair.adf.unw_pwr.bmp $pair.adf.cc 1 .2
rasrmg $pair.adf.unw - $width_rdc 1 1 0 1 1 .5 1. .35 .0 1 $pair.adf.unw.bmp $pair.adf.cc 1 .2

quad_fit $pair.adf.unw $m_date.diff_par 32 32 $pair.adf.cc_mask.bmp $pair.plot 3
quad_sub $pair.adf.unw $m_date.diff_par $pair.adf.unw.sub 0 0

rasrmg $pair.adf.unw.sub -  $width_rdc 1 1 0 1 1 .5 1. .35 .0 1 $pair.adf.unw.sub.bmp $pair.adf.cc 1 .2
rasrmg $pair.adf.unw.sub $m_date.mli $width_rdc 1 1 0 1 1 .6 1. .35 .0 1 $pair.adf.unw.sub_pwr.bmp $pair.adf.cc 1 .2

width_geo=$(awk '$1 == "width:" {print $2}' dem_seg.par)
line_geo=$(awk '$1 == "nlines:" {print $2}' dem_seg.par)

geocode_back $pair.adf.cc $width_rdc $m_date.lookup_fine $pair.adf.cc.geo $width_geo $line_geo 1 0
geocode_back $pair.cor $width_rdc $m_date.lookup_fine $pair.cor.geo $width_geo $line_geo 1 0
geocode_back $pair.adf.diff $width_rdc $m_date.lookup_fine $pair.adf.diff.geo $width_geo $line_geo 1 1
geocode_back $pair.adf.unw $width_rdc $m_date.lookup_fine $pair.adf.unw.geo $width_geo $line_geo 1 0
geocode_back $pair.adf.unw.sub $width_rdc $m_date.lookup_fine $pair.adf.unw.sub.geo $width_geo $line_geo 1 0
geocode_back $m_date.mli $width_rdc $m_date.lookup_fine $m_date.mli.geo $width_geo $line_geo 1 0

rascc $pair.adf.cc.geo - $width_geo 1 1 0 1 1 .1 .9 1. .35 1 $pair.adf.cc.geo.bmp
rascc $pair.cor.geo - $width_geo 1 1 0 1 1 .1 .9 1. .35 1 $pair.cor.geo.bmp

rasmph $pair.adf.diff.geo $width_geo 1 0 1 1 1. 0.35 1 $pair.adf.diff.geo.bmp
rasmph_pwr $pair.adf.diff.geo $m_date.mli.geo $width_geo 1 1 0 1 1 1. 0.35 1 $pair.adf.diff.geo_pwr.bmp

rasrmg $pair.adf.unw.geo - $width_geo 1 1 0 1 1 .18 1. .35 .0 1 $pair.adf.unw.geo.bmp $pair.adf.cc.geo 1 .2
rasrmg $pair.adf.unw.geo $m_date.mli.geo $width_geo 1 1 0 1 1 .18 1. .35 .0 1 $pair.adf.unw.geo_pwr.bmp $pair.adf.cc.geo 1 .2

rasrmg $pair.adf.unw.sub.geo - $width_geo 1 1 0 1 1 .18 1. .35 .0 1 $pair.adf.unw.sub.geo.bmp $pair.adf.cc.geo 1 .2
rasrmg $pair.adf.unw.sub.geo $m_date.mli.geo $width_geo 1 1 0 1 1 .18 1. .35 .0 1 $pair.adf.unw.sub.geo_pwr.bmp $pair.adf.cc.geo 1 .2
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

cc_thres=cc_thres_flag
roff=roff_flag
loff=loff_flag
wavelength=wavelength_flag

pair=$m_date\_$s_date

m_gacos=m_gacos_flag
s_gacos=s_gacos_flag

inc_angle=$(awk '$1 == "incidence_angle:" {print $2}' $m_date.mli.par)

matlab -nodesktop -nosplash -r "process_gacos('$m_gacos','$s_gacos','dem_seg.par',$wavelength,$inc_angle,'$pair.gacos');quit;"

width_rdc=$(awk '$1 == "interferogram_width:" {print $2}' $pair.off)
line_rdc=$(awk '$1 == "interferogram_azimuth_lines:" {print $2}' $pair.off)
width_geo=$(awk '$1 == "width:" {print $2}' dem_seg.par)
line_geo=$(awk '$1 == "nlines:" {print $2}' dem_seg.par)

geocode $m_date.lookup_fine $pair.gacos $width_geo $pair.gacos.rdc $width_rdc $line_rdc 1 0
raspwr $pair.gacos.rdc $width_rdc 1 0 1 1 1. .35 1 $pair.gacos.rdc.bmp

sub_phase $pair.diff $pair.gacos.rdc $m_date.diff_par $pair.diff.gacos 1 0

rasmph_pwr $pair.diff.gacos $m_date.mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.diff.gacos_pwr.bmp

adf $pair.diff.gacos $pair.adf.diff.gacos1 $pair.adf.cc.gacos1 $width_rdc 0.3 128
adf $pair.adf.diff.gacos1 $pair.adf.diff.gacos2 $pair.adf.cc.gacos2 $width_rdc 0.3 64
adf $pair.adf.diff.gacos2 $pair.adf.diff.gacos $pair.adf.cc.gacos $width_rdc 0.3

rm -rf $pair.adf.diff.gacos1 $pair.adf.diff.gacos2 $pair.adf.cc.gacos1 $pair.adf.cc.gacos2

rasmph $pair.adf.diff.gacos $width_rdc 1 0 1 1 1. 0.35 1 $pair.adf.diff.gacos.bmp
rasmph_pwr $pair.adf.diff.gacos $m_date.mli $width_rdc 1 1 0 1 1 1. 0.35 1 $pair.adf.diff.gacos_pwr.bmp

# rascc_mask $pair.adf.cc.gacos $m_date.mli $width_rdc 1 1 0 1 1 $cc_thres 0. .1 .9 1. .35 1 $pair.adf.cc_mask.bmp
rascc_mask $pair.cor $m_date.mli $width_rdc 1 1 0 1 1 $cc_thres 0. .1 .9 1. .35 1 $pair.adf.cc_mask.bmp

mcf $pair.adf.diff.gacos $pair.adf.cc.gacos $pair.adf.cc_mask.bmp $pair.adf.unw.gacos $width_rdc 1 0 0 - - 1 1 - $roff $loff 0

rasrmg $pair.adf.unw.gacos $m_date.mli $width_rdc 1 1 0 1 1 .6 1. .35 .0 1 $pair.adf.unw.gacos_pwr.bmp $pair.adf.cc.gacos 1 .2
rasrmg $pair.adf.unw.gacos  - $width_rdc 1 1 0 1 1 .5 1. .35 .0 1 $pair.adf.unw.gacos.bmp $pair.adf.cc.gacos 1 .2

quad_fit $pair.adf.unw.gacos $m_date.diff_par 32 32 $pair.adf.cc_mask.bmp $pair.plot 3
quad_sub $pair.adf.unw.gacos $m_date.diff_par $pair.adf.unw.sub.gacos 0 0

rasrmg $pair.adf.unw.sub.gacos $m_date.mli $width_rdc 1 1 0 1 1 .6 1. .35 .0 1 $pair.adf.unw.sub.gacos_pwr.bmp $pair.adf.cc.gacos 1 .2
rasrmg $pair.adf.unw.sub.gacos -  $width_rdc 1 1 0 1 1 .5 1. .35 .0 1 $pair.adf.unw.sub.gacos.bmp $pair.adf.cc.gacos 1 .2

geocode_back $pair.adf.diff.gacos $width_rdc $m_date.lookup_fine $pair.adf.diff.gacos.geo $width_geo $line_geo 1 1
geocode_back $pair.adf.unw.gacos $width_rdc $m_date.lookup_fine $pair.adf.unw.gacos.geo $width_geo $line_geo 1 0
geocode_back $pair.adf.unw.sub.gacos $width_rdc $m_date.lookup_fine $pair.adf.unw.sub.gacos.geo $width_geo $line_geo 1 0


rasmph $pair.adf.diff.gacos.geo $width_geo 1 0 1 1 1. 0.35 1 $pair.adf.diff.gacos.geo.bmp
rasmph_pwr $pair.adf.diff.gacos.geo $m_date.mli.geo $width_geo 1 1 0 1 1 1. 0.35 1 $pair.adf.diff.gacos.geo_pwr.bmp

rasrmg $pair.adf.unw.gacos.geo - $width_geo 1 1 0 1 1 .18 1. .35 .0 1 $pair.adf.unw.gacos.geo.bmp $pair.adf.cc.geo 1 .2
rasrmg $pair.adf.unw.gacos.geo $m_date.mli.geo $width_geo 1 1 0 1 1 .18 1. .35 .0 1 $pair.adf.unw.gacos.geo_pwr.bmp $pair.adf.cc.geo 1 .2

rasrmg $pair.adf.unw.sub.gacos.geo - $width_geo 1 1 0 1 1 .18 1. .35 .0 1 $pair.adf.unw.sub.gacos.geo.bmp $pair.adf.cc.geo 1 .2
rasrmg $pair.adf.unw.sub.gacos.geo $m_date.mli.geo $width_geo 1 1 0 1 1 .18 1. .35 .0 1 $pair.adf.unw.sub.gacos.geo_pwr.bmp $pair.adf.cc.geo 1 .2
"""

USAGE = """Example:
  # reference point of unwrapping is (0, 0) and no mask
  python3 dinsar2.py /ly/rslc 20210110 20210122 /ly/dem DInSAR 8 2
  # set reference point of unwrapping and mask
  python3 dinsar2.py /ly/rslc 20210110 20210122 /ly/dem DInSAR 8 2 -r 100 -l 100 -c 0.3
  # for aps correction (Sentinel-1 data)
  python3 dinsar2.py /ly/rslc 20210110 20210122 /ly/dem DInSAR 8 2 -g /ly/gacos -w 0.05546
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(
        description='Run two-pass D-InSAR using GAMMA.\nAPS correction using MATLAB (Only support .ztd file)',
        formatter_class=argparse.RawTextHelpFormatter, epilog=USAGE)

    parser.add_argument('rslc_dir', help='RSLCs directory')
    parser.add_argument('m_date', help='date for master RSLC')
    parser.add_argument('s_date', help='date for slave RSLC')
    parser.add_argument('dem_dir', help='DEM directory contains *.dem and *.dem.par')
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('rlks', help='range looks', type=int)
    parser.add_argument('alks', help='azimuth looks', type=int)
    parser.add_argument('-e', dest='rslc_extension', type=str, default='.rslc', help='file extension for RSLC (defaults: .rslc)')
    parser.add_argument('-g', dest='gacos_dir', help='directory contains GACOS files for aps correction')
    parser.add_argument('-w', dest='wavelength', type=float, help='Microwave length (Sentinel-1: 0.05546576, ALOS: 0.23830879)')
    parser.add_argument('-r', dest='roff', help='phase reference range offset to unwrap (defaults: 0)', type=float, default=0)
    parser.add_argument('-l', dest='loff', help='phase reference azimuth offset to unwrap (defaults: 0)', type=float, default=0)
    parser.add_argument('-c', dest='cc_thres', type=float, default=0, help='threshold of correlation for creating the unwrapping mask (0.0 --> 1.0) (defaults: 0)')

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
    m_date = inps.m_date
    s_date = inps.s_date
    dem_dir = os.path.abspath(inps.dem_dir)
    out_dir = os.path.abspath(inps.out_dir)
    rlks = inps.rlks
    alks = inps.alks

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

    # check RSLC extension
    if not rslc_extension.startswith('.'):
        rslc_extension = '.' + rslc_extension

    # check m_date rslc and s_date rslc
    m_rslc = os.path.join(rslc_dir, m_date, m_date + rslc_extension)
    s_rslc = os.path.join(rslc_dir, s_date, s_date + rslc_extension)
    if not os.path.isfile(m_rslc):
        sys.exit('Cannot find RSLC for {}'.format(m_date))
    if s_date not in dates:
        sys.exit('Cannot find RSLC for {}'.format(s_date))

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

    # check gacos_dir and wavelength
    if gacos_dir and wavelength is None:
        sys.exit('wavelength(-w) is required for aps correction.')
    if gacos_dir is None and wavelength:
        sys.exit('gacos_dir(-g) is required for aps correction.')
    if gacos_dir and wavelength:
        gacos_dir = os.path.abspath(gacos_dir)
        if not os.path.isdir(gacos_dir):
            sys.exit('{} does not exist.'.format(gacos_dir))

        check_gacos([m_date, s_date], gacos_dir, '.ztd')

    # DInSAR
    # replace value
    dinsar_script_out = dinsar_script.replace('m_date_flag', m_date)
    dinsar_script_out = dinsar_script_out.replace('s_date_flag', s_date)
    dinsar_script_out = dinsar_script_out.replace('m_rslc_flag', m_rslc)
    dinsar_script_out = dinsar_script_out.replace('s_rslc_flag', s_rslc)
    dinsar_script_out = dinsar_script_out.replace('dem_flag', dem)
    dinsar_script_out = dinsar_script_out.replace('dem_par_flag', dem_par)
    dinsar_script_out = dinsar_script_out.replace('rlks_flag', str(rlks))
    dinsar_script_out = dinsar_script_out.replace('alks_flag', str(alks))
    dinsar_script_out = dinsar_script_out.replace('cc_thres_flag', str(cc_thres))
    dinsar_script_out = dinsar_script_out.replace('roff_flag', str(roff))
    dinsar_script_out = dinsar_script_out.replace('loff_flag', str(loff))

    # write script
    os.chdir(out_dir)
    run_script = m_date + '_' + s_date + '_DInSAR.sh'
    with open(run_script, 'w+') as f:
        f.write(dinsar_script_out)

    # run script
    call_str = 'bash ' + run_script
    os.system(call_str)

    # aps correction
    if gacos_dir and wavelength:
        gacos_dir = os.path.abspath(gacos_dir)

        os.chdir(out_dir)

        # write matlab script
        with open('process_gacos.m', 'w+') as f:
                f.write(matlab_code)

        m_gacos = os.path.join(gacos_dir, m_date + '.ztd')
        s_gacos = os.path.join(gacos_dir, s_date + '.ztd')

        # replace value
        aps_correction_script_out = aps_correction_script.replace('m_date_flag', m_date)
        aps_correction_script_out = aps_correction_script_out.replace('s_date_flag', s_date)
        aps_correction_script_out = aps_correction_script_out.replace('m_gacos_flag', m_gacos)
        aps_correction_script_out = aps_correction_script_out.replace('s_gacos_flag', s_gacos)
        aps_correction_script_out = aps_correction_script_out.replace('wavelength_flag', str(wavelength))
        aps_correction_script_out = aps_correction_script_out.replace('cc_thres_flag', str(cc_thres))
        aps_correction_script_out = aps_correction_script_out.replace('roff_flag', str(roff))
        aps_correction_script_out = aps_correction_script_out.replace('loff_flag', str(loff))

        # write script
        run_script = m_date + '_' + s_date + '_aps_correction.sh'
        with open(run_script, 'w+') as f:
            f.write(aps_correction_script_out)

        # run script
        call_str = 'bash ' + run_script
        os.system(call_str)

    print('\nAll done, enjoy it!\n')

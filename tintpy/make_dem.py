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

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from osgeo import gdal


def read_gdal_file(file, band=1):
    """Read GDAL supported file using GDAL

    Args:
        file (str): GDAL supported file
        band (int, optional): band number. Defaults to 1.

    Returns:
        tuple: data, lon_lat, size
    """
    ds = gdal.Open(file, gdal.GA_ReadOnly)
    data = ds.GetRasterBand(band).ReadAsArray()

    trans = ds.GetGeoTransform()
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    lon_lat = [
        trans[0], trans[0] + xsize * trans[1], trans[3] + ysize * trans[5], trans[3]
    ]
    ds = None
    xstep = trans[1]
    ystep = trans[5]
    size = [xsize, ysize, xstep, ystep]

    return data, lon_lat, size


def draw_img(data, lon_lat, out_file, cmap='rainbow'):
    """Draw image from array

    Args:
        data (array): numpy array
        lon_lat (list): lon_lat of data
        out_file (str): output file name
        cmap (str, optional): colormap. Defaults to 'rainbow'.
    """
    plt.figure()
    ax = plt.gca()
    im = ax.imshow(data, extent=lon_lat, cmap=cmap)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    plt.colorbar(im, cax=cax)
    plt.savefig(out_file, bbox_inches='tight', dpi=300)


def mosaic_tifs(tifs, out_tif):
    """Mosaic tifs

    Args:
        tifs (list): tifs
        out_tif (str): mosaiced tif
    """
    if os.path.isfile(out_tif):
        os.remove(out_tif)
    try:
        gdal.Warp(out_tif, tifs)
    except Exception as e:
        print(e)


def write_gamma_par(par_file, lon_lat, size, data_format):
    """Write GAMMA format parameter

    Args:
        par_file (str): output par file
        lon_lat (list): longitude and latitude
        size (list): xsize ysize xstep ystep
        data_format (str): data format
    """
    upper_left_lon, upper_left_lat = lon_lat[0], lon_lat[3]
    xstep, ystep = size[2], size[3]
    width, length = size[0], size[1]

    f = open(par_file, 'w')

    title = os.path.basename(par_file)[:-4]
    proj = 'EQA'

    f.write("Gamma DIFF&GEO DEM/MAP parameter file\n")
    f.write("title:\t%s\n" % title)
    # Projection should be checked
    f.write("DEM_projection:     %s\n" % proj)
    # INTEGER*2 (short) OR REAL*4 (float32) should be modified
    f.write("data_format:        %s\n" % data_format)
    f.write("DEM_hgt_offset:          0.00000\n")
    f.write("DEM_scale:               1.00000\n")
    f.write("width:                %s\n" % width)
    f.write("nlines:               %s\n" % length)
    f.write("corner_lat:   %s  decimal degrees\n" % upper_left_lat)
    f.write("corner_lon:   %s  decimal degrees\n" % upper_left_lon)
    f.write("post_lat:   %s  decimal degrees\n" % -abs(ystep))
    f.write("post_lon:   %s  decimal degrees\n" % abs(xstep))
    f.write("\n")
    f.write("ellipsoid_name: WGS 84\n")
    f.write("ellipsoid_ra:        6378137.000   m\n")
    f.write("ellipsoid_reciprocal_flattening:  298.2572236\n")
    f.write("\n")
    f.write("datum_name: WGS 1984\n")
    f.write("datum_shift_dx:              0.000   m\n")
    f.write("datum_shift_dy:              0.000   m\n")
    f.write("datum_shift_dz:              0.000   m\n")
    f.write("datum_scale_m:         0.00000e+00\n")
    f.write("datum_rotation_alpha:  0.00000e+00   arc-sec\n")
    f.write("datum_rotation_beta:   0.00000e+00   arc-sec\n")
    f.write("datum_rotation_gamma:  0.00000e+00   arc-sec\n")
    f.write("datum_country_list Global Definition, WGS84, World\n")
    f.write("\n")
    f.close()


def write_saracape_par(hdr_file, sml_file, lon_lat, size):
    """Write Sarscape format parameter

    Args:
        hdr_file (str): output hdr file
        sml_file (str): output sml file
        lon_lat (list): longitude and latitude
        size (list): xsize ysize xstep ystep
    """
    hdr = open(hdr_file, 'w')

    hdr.write("ENVI\n")
    hdr.write("description = {\n")
    hdr.write("ANCILLARY INFO = DEM.\n")
    hdr.write("File generated with SARscape  5.2.1 }\n")
    hdr.write("\n")
    hdr.write(f"samples                   = {size[0]}\n")
    hdr.write(f"lines                     = {size[1]}\n")
    hdr.write("bands                     = 1\n")
    hdr.write("headeroffset              = 0\n")
    hdr.write("file type                 = ENVI Standard\n")
    hdr.write("data type                 = 2\n")
    hdr.write("sensor type               = Unknown\n")
    hdr.write("interleave                = bsq\n")
    hdr.write("byte order                = 0\n")
    hdr.write(
        "map info = {Geographic Lat/Lon, 1, 1, %s, %s, %s, %s, WGS-84, \n" %
        (lon_lat[0], lon_lat[3], abs(size[2]), abs(size[3])))
    hdr.write("units=Degrees}\n")
    hdr.write("x start                   = 1\n")
    hdr.write("y start                   = 1\n")
    hdr.close()

    sml = open(sml_file, 'w')

    sml.write('<?xml version="1.0" ?>\n')
    sml.write(
        '<HEADER_INFO xmlns="http://www.sarmap.ch/xml/SARscapeHeaderSchema"\n')
    sml.write('    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\n')
    sml.write(
        '    xsi:schemaLocation="http://www.sarmap.ch/xml/SARscapeHeaderSchema\n'
    )
    sml.write(
        '    http://www.sarmap.ch/xml/SARscapeHeaderSchema/SARscapeHeaderSchema_version_1.0.xsd">\n'
    )
    sml.write('<RegistrationCoordinates>\n')
    sml.write(f'    <LatNorthing>{lon_lat[3]}</LatNorthing>\n')
    sml.write(f'    <LonEasting>{lon_lat[0]}</LonEasting>\n')
    sml.write(
        f'    <PixelSpacingLatNorth>{-abs(size[3])}</PixelSpacingLatNorth>\n')
    sml.write(
        f'    <PixelSpacingLonEast>{abs(size[2])}</PixelSpacingLonEast>\n')
    sml.write('    <Units>DEGREES</Units>\n')
    sml.write('</RegistrationCoordinates>\n')
    sml.write('<RasterInfo>\n')
    sml.write('    <HeaderOffset>0</HeaderOffset>\n')
    sml.write('    <RowPrefix>0</RowPrefix>\n')
    sml.write('    <RowSuffix>0</RowSuffix>\n')
    sml.write('    <FooterLen>0</FooterLen>\n')
    sml.write('    <CellType>SHORT</CellType>\n')
    sml.write('    <DataUnits>DEM</DataUnits>\n')
    sml.write('    <NullCellValue>-32768</NullCellValue>\n')
    sml.write(f'    <NrOfPixelsPerLine>{size[1]}</NrOfPixelsPerLine>\n')
    sml.write(f'    <NrOfLinesPerImage>{size[0]}</NrOfLinesPerImage>\n')
    sml.write('    <GeocodedImage>OK</GeocodedImage>\n')
    sml.write('    <BytesOrder>LSBF</BytesOrder>\n')
    sml.write('    <OtherInfo>\n')
    sml.write(
        '        <MatrixString NumberOfRows = "1" NumberOfColumns = "2">\n')
    sml.write('            <MatrixRowString ID = "0">\n')
    sml.write('            <ValueString ID = "0">SOFTWARE</ValueString>\n')
    sml.write(
        '            <ValueString ID = "1">SARscape ENVI  5.1.0 Sep  8 2014  W64</ValueString>\n'
    )
    sml.write('            </MatrixRowString>\n')
    sml.write('        </MatrixString>\n')
    sml.write('    </OtherInfo>\n')
    sml.write('</RasterInfo>\n')
    sml.write('<CartographicSystem>\n')
    sml.write('    <State>GEO-GLOBAL</State>\n')
    sml.write('    <Hemisphere></Hemisphere>\n')
    sml.write('    <Projection>GEO</Projection>\n')
    sml.write('    <Zone></Zone>\n')
    sml.write('    <Ellipsoid>WGS84</Ellipsoid>\n')
    sml.write('    <DatumShift></DatumShift>\n')
    sml.write('</CartographicSystem>\n')
    sml.write('<DEMCoordinates>\n')
    sml.write(
        f'    <EastingCoordinateBegin>{lon_lat[0]}</EastingCoordinateBegin>\n')
    sml.write(
        f'    <EastingCoordinateEnd>{lon_lat[1]}</EastingCoordinateEnd>\n')
    sml.write(f'    <EastingGridSize>{abs(size[3])}</EastingGridSize>\n')
    sml.write(
        f'    <NorthingCoordinateBegin>{lon_lat[2]}</NorthingCoordinateBegin>\n'
    )
    sml.write(
        f'    <NorthingCoordinateEnd>{lon_lat[3]}</NorthingCoordinateEnd>\n')
    sml.write(f'    <NorthingGridSize>{abs(size[2])}</NorthingGridSize>\n')
    sml.write('</DEMCoordinates>\n')
    sml.write('</HEADER_INFO>\n')
    sml.close()


def make_dem(processor, tif, out_file):
    """Make GAMMA and Sarscape dem

    Args:
        processor (str): InSAR processor
        tif (str): tif file
        out_file (str): output dem file name

    Returns:
        tuple: dem_data, lon_lat, size
    """
    dem_data, lon_lat, size = read_gdal_file(tif)
    ori_dem_data = dem_data.copy()

    if processor == 'GAMMA':
        if dem_data.dtype == 'float32':
            data_format = 'REAL*4'
        else:
            data_format = 'INTEGER*2'

        dem_data.byteswap(True)

        if not out_file.endswith('.dem'):
            out_file += '.dem'

        dem_data.tofile(out_file)
        write_gamma_par(out_file + '.par', lon_lat, size, data_format)

    if processor == 'SARSCAPE':
        if not out_file.endswith('_dem'):
            out_file += '_dem'

        dem_data.tofile(out_file)
        write_saracape_par(out_file + '.hdr', out_file + '.sml', lon_lat, size)

    dem_data = None
    
    return ori_dem_data, lon_lat, size


def uncompress_zips(zip_files, out_dir):
    """Uncompress zips file

    Args:
        zip_files (list): zip files
        out_dir (str): output directory

    Returns:
        list: uncompressed tifs
    """
    tifs = []
    for zip_file in zip_files:
        with zipfile.ZipFile(zip_file, mode='r') as f:
            files = f.namelist()
            for file in files:
                if file.endswith('hgt') or file.endswith('tif'):
                    f.extract(file, out_dir)
                    tifs.append(os.path.join(out_dir, file))
    return tifs


EXAMPLE = '''Example:
  python3 make_dem.py sarscape ly 1.tif
  python3 make_dem.py sarscape ly 1.tif 2.tif
  python3 make_dem.py gamma ly 1.hgt 2.hgt
  python3 make_dem.py gamma ly srtm_57_08.zip
  python3 make_dem.py gamma ly N20E100.SRTMGL1.hgt.zip N21E100.SRTMGL1.hgt.zip
'''


def cmdline_parser():
    parser = argparse.ArgumentParser(description='Make DEM used in interferometry both for GAMMA and SARSCAPE processor.',
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)
    parser.add_argument('processor', type=str, choices=['GAMMA', 'gamma', 'SARSCAPE', 'sarscape'],
        help='interferometry processor [ gamma or sarscape ]')
    parser.add_argument('out_file', type=str, help='output dem file name')
    parser.add_argument('tifs', type=str, nargs='+', help='tif(hgt) or zips files for making DEM')
    inps = parser.parse_args()

    return inps


def main():
    # get inputs
    inps = cmdline_parser()
    tifs = [os.path.abspath(t) for t in inps.tifs]
    processor = inps.processor.upper()
    out_file = os.path.abspath(inps.out_file)

    # check tif files
    not_exist = []
    for tif in tifs:
        if not os.path.isfile(tif):
            not_exist.append(tif)
    if not_exist:
        for i in not_exist:
            print(f"{i} does not exist, please check it.")
        sys.exit()

    # check output file nname
    dir_name = os.path.dirname(out_file)
    if not os.path.isdir(dir_name):
        sys.exit(f"{dir_name} does not exist, please check it.")

    # uncompress zips
    flag = False
    if tifs[0].endswith('zip'):
        tifs = uncompress_zips(tifs, os.path.dirname(tifs[0]))
        flag = True

    # make dem with multi tifs
    if len(tifs) > 1:
        print('Mosaic data.')
        mosaiced_tif = out_file + '.tif'
        mosaic_tifs(tifs, mosaiced_tif)

        if os.path.isfile(mosaiced_tif):
            print(f"Make {processor} dem.")
            data, lon_lat, _ = make_dem(processor, mosaiced_tif, out_file)

            print('Plot data.')
            draw_img(data, lon_lat, out_file + '.png')

            data = None
            print('All done, enjoy it!\n')
    # make dem with one tif
    elif len(tifs) == 1:
        print(f"Make {processor} dem.")
        data, lon_lat, _ = make_dem(processor, tifs[0], out_file)

        print('Plot data.')
        draw_img(data, lon_lat, out_file + '.png')

        data = None
        print('All done, enjoy it!\n')
    else:
        sys.exit('No tif(hgt) files.')

    # delete uncompressed tifs
    if flag:
        for tif in tifs:
            os.remove(tif)


if __name__ == "__main__":
    main()

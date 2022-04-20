#!/usr/bin/env python3
################################
# Program is part of TintPy    #
# Copyright (c) 2022, Lei Yuan #
# Author: Lei Yuan, 2022       #
################################

import argparse
import os
import sys
import zipfile

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from lxml import etree
from pykml.factory import KML_ElementMaker as KML

EXAMPLE = """Example:
  # geo files
  python3 ph2kmz.py geo ph_rate dem_seg.par mli wrapped_res 6 -d 1200
  # rdc files
  python3 ph2kmz.py rdc ph_rate dem_seg.par mli wrapped_res 6 -l lookup_fine -p diff.par -d 1200
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(
        description='Display wrapped phase with intensity in Google Earth.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLE)

    parser.add_argument('flag',
                        type=str,
                        choices=['geo', 'GEO', 'rdc', 'RDC'],
                        help='file flag (geo or rdc).')
    parser.add_argument('phase', type=str, help='phase file.')
    parser.add_argument('dem_seg_par', type=str, help='dem_seg.par.')
    parser.add_argument('mli', type=str, help='intensity file.')
    parser.add_argument('out_dir', type=str, help='output directory.')
    parser.add_argument('cycles',
                        type=float,
                        nargs='+',
                        help='data value per color cycle (6 means -3 ~ 3).')
    parser.add_argument('-l', dest='lookup', type=str, help='lookup file.')
    parser.add_argument('-p', dest='diff_par', type=str, help='diff_par file.')
    parser.add_argument('-d',
                        dest='dpi',
                        type=int,
                        default=600,
                        help='image dpi (defaults:600)')

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


def geocode_back(infile, lookup_file, outfile, width_rdr, width_geo,
                 lines_geo):
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


def draw_img_with_pwr(img_data, cmap, pwr_data=None, vlim=None, out_file=None, withcmap=False, dpi=600):
    """Draw image with intensity data

    Args:
        img_data (array): data for drawing image
        cmap (str): colorbar name
        pwr_data (array, optional): intensity data. Defaults to None.
        vlim (list, optional): value limitation. Defaults to None.
        out_file (str, optional): output image name. Defaults to None.
        withcmap (bool, optional): whether to draw colorbar. Defaults to False.
        dpi (int, optional): output image dpi. Defaults to 600.
    """
    if vlim:
        minv, maxv = vlim
    else:
        minv, maxv = np.nanmin(img_data), np.nanmax(img_data)

    img_data[img_data == 0.0] = np.nan

    # scale data to the range 0 -> 1
    dt = (img_data - minv) / (maxv - minv)

    # create RGB display of the data in the range 0 -> 1
    cmap = mpl.cm.get_cmap(cmap)
    img_rgba = cmap(dt, bytes=True)
    r, g, b, a = np.rollaxis(img_rgba, axis=-1)

    if pwr_data is not None:
        exp, sc = 0.35, 1.0
        pwr = np.abs(pwr_data)
        # all points that are non-zero
        idx = np.nonzero(pwr)
        # scale these data values and evaluate the mean value
        pwr = np.power(pwr, exp)
        ave = pwr[idx].mean()
        # determine default clip value
        cv = 2.0 * ave / sc
        v = (pwr.clip(0., cv)) / cv
        v[v == 0.0] = np.nan
        # multiply by v to get scaled RGB
        r = (r.astype(np.float32) * v).astype(np.uint8)
        b = (b.astype(np.float32) * v).astype(np.uint8)
        g = (g.astype(np.float32) * v).astype(np.uint8)

    img_rgba = np.dstack([r, g, b, a])

    fig, ax = plt.subplots()
    ax.axis('off')

    if pwr_data is None:
        ax.imshow(img_rgba)
    else:
        ax.imshow(v, zorder=0, cmap='gray')
        ax.imshow(img_rgba, zorder=1)

    if withcmap:
        norm = mpl.colors.Normalize(vmin=minv, vmax=maxv)
        fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                     ax=ax,
                     orientation='vertical',
                     pad=0.02,
                     ticks=[minv, 0, maxv])

    if out_file:
        fig.savefig(out_file, pad_inches=0.0, bbox_inches='tight', dpi=dpi, transparent=True)


def draw_colorbar(vlim, cmap, label, out_file, figsize=(0.18, 3.6)):
    """Draw colorbar

    Args:
        vlim (list): value limitation
        cmap (str): colorbar name
        label (str): label for colorbar
        out_file (str): output colorbar name
        figsize (tuple, optional): colorbar size. Defaults to (0.18, 3.6).
    """
    fig, ax = plt.subplots(figsize=figsize)

    cmap = plt.get_cmap(cmap)
    minv, maxv = vlim
    norm = mpl.colors.Normalize(vmin=minv, vmax=maxv)
    fig.patch.set_facecolor('white')
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                 cax=ax,
                 orientation='vertical',
                 label=label,
                 ticks=[minv, 0, maxv])
    fig.savefig(out_file, bbox_inches='tight', dpi=300)


def make_kmz(data, lon_lat, out_kmz, img_dpi=600, pwr_data=None, cmap='jet', vlim=None, label=None):
    """Make kmz

    Args:
        data (array): data for making kmz
        lon_lat (list): longitude and latitude
        out_kmz (str): output kmz name
        img_dpi (int, optional): image dpi. Defaults to 600.
        pwr_data (array, optional): intensity data. Defaults to None.
        cmap (str, optional): colorbar name. Defaults to 'jet'.
        vlim (list, optional): value limitation. Defaults to None.
        label (str, optional): label for colorbar. Defaults to None.
    """
    out_kmz = os.path.abspath(out_kmz)
    if out_kmz.endswith('.kmz'):
        out_kmz = out_kmz[:-4]

    out_dir = os.path.dirname(out_kmz)

    os.chdir(out_dir)

    kmz = os.path.basename(out_kmz)
    cbar = 'cbar.png'
    img = 'img.png'

    if vlim is None:
        vlim = [np.nanmin(data), np.nanmax(data)]

    print('drawing colorbar')
    draw_colorbar(vlim, cmap, label, cbar)

    if pwr_data is not None:
        print('drawing image with pwr data')
    else:
        print('drawing image')
    draw_img_with_pwr(data, cmap, pwr_data=pwr_data, vlim=vlim, out_file=img, withcmap=False, dpi=img_dpi)

    south, north, west, east = lon_lat

    doc = KML.kml(KML.Folder(KML.name(kmz)))
    img_displayed = KML.GroundOverlay(
        KML.name('img'), KML.Icon(KML.href(img)),
        KML.LatLonBox(KML.north(str(north)), KML.south(str(south)),
                      KML.east(str(east)), KML.west(str(west))))
    doc.Folder.append(img_displayed)

    legend = KML.ScreenOverlay(
        KML.name('colorbar'), KML.Icon(KML.href(cbar), KML.viewBoundScale(0.75)),
        KML.overlayXY(
            x="0.0",
            y="1",
            xunits="fraction",
            yunits="fraction",
        ), KML.screenXY(
            x="0.0",
            y="1",
            xunits="fraction",
            yunits="fraction",
        ),
        KML.rotationXY(
            x="0.",
            y="1.",
            xunits="fraction",
            yunits="fraction",
        ), KML.size(
            x="0",
            y="250",
            xunits="pixel",
            yunits="pixel",
        ), KML.visibility(1), KML.open(0))

    doc.Folder.append(legend)

    kml_str = etree.tostring(doc, pretty_print=True)

    kml = kmz + '.kml'
    with open(kml, 'wb') as f:
        f.write(kml_str)

    kmz = kmz + '.kmz'
    print('writing data to {}'.format(kmz))

    with zipfile.ZipFile(kmz, 'w') as f:
        f.write(kml)
        os.remove(kml)
        f.write(img)
        os.remove(img)
        f.write(cbar)
        os.remove(cbar)


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


def get_lon_lat(dem_par):
    """Get longitude and latitude from dem_seg.par

    Args:
        dem_par (str): dem_seg.par

    Returns:
        list: [S, N, W, E]
    """
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

    return [south, north, west, east]


def wrap_data(data, wrap_range=[-np.pi, np.pi]):
    """Wrap data

    Args:
        data (array): data for wrapping
        wrap_range (list, optional): wrap range. Defaults to [-np.pi, np.pi].

    Returns:
        array: wrapped data
    """
    w0, w1 = wrap_range
    wrapped_data = w0 + np.mod(data - w0, w1 - w0)

    return wrapped_data


# register rmg cmap
cdict = {
    'red': [(0.0, 1.0, 1.0), (0.16665, 1.0, 1.0), (0.5, 0.0, 0.0), (0.83333, 1.0, 1.0), (1.0, 1.0, 1.0)],
    'green': [(0.0, 0.5, 0.5), (0.16665, 1.0, 1.0), (0.5, 1.0, 1.0), (0.83333, 0.0, 0.0), (1.0, 0.5, 0.5)],
    'blue': [(0.0, 0.5, 0.5), (0.16665, 0.0, 0.0), (0.5, 1.0, 1.0), (0.83333, 1.0, 1.0), (1.0, 0.5, 0.5)]
}
cm_rmg = mpl.colors.LinearSegmentedColormap('rmg', cdict)
plt.cm.register_cmap(cmap=cm_rmg)


def main():
    # get inps
    inps = cmdline_parser()
    flag = inps.flag.lower()
    phase = os.path.abspath(inps.phase)
    dem_seg_par = os.path.abspath(inps.dem_seg_par)
    mli = os.path.abspath(inps.mli)
    out_dir = os.path.abspath(inps.out_dir)
    cycles = inps.cycles
    dpi = inps.dpi

    # check file
    check_file(phase)
    check_file(dem_seg_par)
    check_file(mli)
    if flag == 'rdc':
        lookup = os.path.abspath(inps.lookup)
        diff_par = os.path.abspath(inps.diff_par)
        check_file(lookup)
        check_file(diff_par)

    # check out_dir
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    os.chdir(out_dir)

    if flag == 'rdc':
        width_rdc = read_gamma_par(diff_par, 'range_samp_1')
        width = read_gamma_par(dem_seg_par, 'width')
        nlines = read_gamma_par(dem_seg_par, 'nlines')

        phase_geo = os.path.basename(phase) + '.geo'
        mli_geo = os.path.basename(mli) + '.geo'

        if not os.path.isfile(phase_geo):
            geocode_back(phase, lookup, phase_geo, width_rdc, width, nlines)
        if not os.path.isfile(mli_geo):
            geocode_back(mli, lookup, mli_geo, width_rdc, width, nlines)

        phase = phase_geo
        mli = mli_geo

    nlines = int(read_gamma_par(dem_seg_par, 'nlines'))

    phase_data = read_gamma(phase, nlines, 'float32')
    pwr_data = read_gamma(mli, nlines, 'float32')

    lon_lat = get_lon_lat(dem_seg_par)

    for cycle in cycles:
        wrapped_phase = wrap_data(phase_data, [-cycle / 2, cycle / 2])
        cycle = str(cycle)
        if cycle.endswith('.0'):
            cycle = cycle[0:-2]
        kmz = 'wrap_' + cycle + '.kmz'
        make_kmz(wrapped_phase, lon_lat, kmz, img_dpi=dpi, pwr_data=pwr_data, cmap='rmg')

    print('\nAll done, enjoy it!\n')


if __name__ == "__main__":
    main()

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
  python3 gamma2kmz.py geo ph_rate.geo dem_seg.par res -r 6 -m mli.geo
  python3 gamma2kmz.py geo ph_rate.geo dem_seg.par res -r 6 10 20 -m mli.geo
  # rdc files
  python3 gamma2kmz.py rdc ph_rate dem_seg.par res -r 6 -m mli -l lookup_fine -p diff.par
  python3 gamma2kmz.py rdc adf.diff dem_seg.par res -m mli -l lookup_fine -p diff.par
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(
        description='Convert data supported by GAMMA to kmz for display.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLE)

    parser.add_argument('flag',
                        type=str,
                        choices=['geo', 'GEO', 'rdc', 'RDC'],
                        help='file flag (geo or rdc).')
    parser.add_argument('file', type=str, help='file for making kmz.')
    parser.add_argument('dem_seg_par', type=str, help='dem_seg.par.')
    parser.add_argument('out_dir', type=str, help='output directory.')
    parser.add_argument('-r',
                        dest='wrap_range',
                        type=float,
                        nargs='+',
                        help='wrap range (6 means -3 ~ 3).')
    parser.add_argument('-m', dest='mli', type=str, help='intensity file.')
    parser.add_argument('-l', dest='lookup', type=str, help='lookup file.')
    parser.add_argument('-p', dest='diff_par', type=str, help='diff_par file.')
    parser.add_argument('-c',
                        dest='cmap',
                        default='rmg',
                        type=str,
                        help='colorbar (defaults: rmg).')
    parser.add_argument('-t',
                        dest='type',
                        default='float32',
                        type=str,
                        choices=['float32', 'complex64'],
                        help='file type (defaults: float32).')
    parser.add_argument('-d',
                        dest='dpi',
                        type=int,
                        default=1200,
                        help='image dpi (defaults: 1200)')

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


def geocode_back(infile, lookup_file, outfile, width_rdr, width_geo, lines_geo,
                 file_type):
    """Geocoding of image data using lookup table values

    Args:
        infile (str): file
        lookup_file (str): lookup table
        outfile (str): geocoded file
        width_rdr (int): rdr width
        width_geo (int): geo width
        lines_geo (int): geo length
        file_type (int): input/output file type (0 for float, 1 for fcomplex)
    """
    call_str = f"geocode_back {infile} {width_rdr} {lookup_file} {outfile} {width_geo} {lines_geo} 1 {file_type}"
    os.system(call_str)


def check_file(file):
    """Check if the file exists 

    Args:
        file (str): file
    """
    if not os.path.isfile(file):
        sys.exit("{} does not exist.".format(file))


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


def scale_pwr(pwr_data):
    """Scale pwr data

    Args:
        pwr_data (array): intensity data

    Returns:
        array: scaled pwr
    """
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

    return v


def draw_img_with_pwr(data,
                      cmap,
                      pwr=None,
                      vlim=None,
                      out_file=None,
                      withcmap=False,
                      dpi=600):
    """Draw image with pwr data

    Args:
        data (array): data for drawing image
        cmap (str): colorbar name
        pwr (array, optional): scaled pwr data. Defaults to None.
        vlim (list, optional): value limitation. Defaults to None.
        out_file (str, optional): output image name. Defaults to None.
        withcmap (bool, optional): whether to draw colorbar. Defaults to False.
        dpi (int, optional): output image dpi. Defaults to 600.
    """
    if vlim:
        minv, maxv = vlim
    else:
        minv, maxv = np.nanmin(data), np.nanmax(data)

    data[data == 0.0] = np.nan

    # scale data to the range 0 -> 1
    dt = (data - minv) / (maxv - minv)

    # create RGB display of the data in the range 0 -> 1
    cmap = mpl.cm.get_cmap(cmap)
    img_rgba = cmap(dt, bytes=True)
    r, g, b, a = np.rollaxis(img_rgba, axis=-1)

    if pwr is not None:
        # multiply by v to get scaled RGB
        r = (r.astype(np.float32) * pwr).astype(np.uint8)
        b = (b.astype(np.float32) * pwr).astype(np.uint8)
        g = (g.astype(np.float32) * pwr).astype(np.uint8)

    img_rgba = np.dstack([r, g, b, a])

    fig, ax = plt.subplots()
    ax.axis('off')

    if pwr is None:
        ax.imshow(img_rgba)
    else:
        ax.imshow(pwr, zorder=0, cmap='gray')
        ax.imshow(img_rgba, zorder=1)

    if withcmap:
        norm = mpl.colors.Normalize(vmin=minv, vmax=maxv)
        fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                     ax=ax,
                     orientation='vertical',
                     pad=0.02,
                     ticks=[minv, 0, maxv])

    if out_file:
        print(f"drawing {out_file}")
        fig.savefig(out_file,
                    pad_inches=0.0,
                    bbox_inches='tight',
                    dpi=dpi,
                    transparent=True)
        print("done")


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
    print(f"drawing {out_file}")
    fig.savefig(out_file, bbox_inches='tight', dpi=300)
    print("done")


def make_kmz(img, cbar, lon_lat, kmz_name, out_dir):
    """Make kmz

    Args:
        img (str): image file
        cbar (str): colorbar file
        lon_lat (list): [S N W E]
        kmz_name (str): output kmz file name
        out_dir (str): output directory
    """
    if kmz_name.endswith('.kmz'):
        kmz_name = kmz_name[:-4]

    s, n, w, e = lon_lat

    doc = KML.kml(KML.Folder(KML.name(kmz_name)))
    img_displayed = KML.GroundOverlay(
        KML.name('image'), KML.Icon(KML.href(img)),
        KML.LatLonBox(KML.north(str(n)), KML.south(str(s)), KML.east(str(e)),
                      KML.west(str(w))))
    doc.Folder.append(img_displayed)

    legend = KML.ScreenOverlay(
        KML.name('colorbar'), KML.Icon(KML.href(cbar),
                                       KML.viewBoundScale(0.75)),
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

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    os.chdir(out_dir)

    kml = kmz_name + '.kml'
    with open(kml, 'wb') as f:
        f.write(kml_str)

    kmz = kmz_name + '.kmz'
    print('writing data to {}'.format(kmz))

    with zipfile.ZipFile(kmz, 'w') as f:
        f.write(kml)
        os.remove(kml)
        f.write(img)
        f.write(cbar)
    print('done')


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
    'red': [(0.0, 1.0, 1.0), (0.16665, 1.0, 1.0), (0.5, 0.0, 0.0),
            (0.83333, 1.0, 1.0), (1.0, 1.0, 1.0)],
    'green': [(0.0, 0.5, 0.5), (0.16665, 1.0, 1.0), (0.5, 1.0, 1.0),
              (0.83333, 0.0, 0.0), (1.0, 0.5, 0.5)],
    'blue': [(0.0, 0.5, 0.5), (0.16665, 0.0, 0.0), (0.5, 1.0, 1.0),
             (0.83333, 1.0, 1.0), (1.0, 0.5, 0.5)]
}
cm_rmg = mpl.colors.LinearSegmentedColormap('rmg', cdict)
plt.cm.register_cmap(cmap=cm_rmg)


def main():
    # get inps
    inps = cmdline_parser()
    flag = inps.flag.lower()
    main_file = os.path.abspath(inps.file)
    dem_seg_par = os.path.abspath(inps.dem_seg_par)
    out_dir = os.path.abspath(inps.out_dir)

    wrap_range = inps.wrap_range
    mli = inps.mli
    lookup = inps.lookup
    diff_par = inps.diff_par
    cmap = inps.cmap
    dpi = inps.dpi
    file_type = inps.type

    # check file
    check_file(main_file)
    check_file(dem_seg_par)
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    if flag == 'rdc':
        if not lookup:
            sys.exit("lookup(-l) is required for rdc method.")
        if not diff_par:
            sys.exit("diff_par(-p) is required for rdc method.")
        lookup = os.path.abspath(lookup)
        diff_par = os.path.abspath(diff_par)
        check_file(lookup)
        check_file(diff_par)

        width_rdc = read_gamma_par(diff_par, 'range_samp_1')
        width = read_gamma_par(dem_seg_par, 'width')
        nlines = read_gamma_par(dem_seg_par, 'nlines')

        main_file_geo = main_file + '.geo'

        if not os.path.isfile(main_file_geo):
            if file_type == 'float32':
                f_type = 0
            else:
                f_type = 1
            geocode_back(main_file, lookup, main_file_geo, width_rdc, width,
                         nlines, f_type)

        main_file = main_file_geo

    nlines = int(read_gamma_par(dem_seg_par, 'nlines'))
    data = read_gamma(main_file, nlines, file_type)
    if file_type == 'complex64':
        data = np.angle(data)

    if mli:
        mli = os.path.abspath(mli)
        check_file(mli)
        if flag == 'rdc':
            mli_geo = mli + '.geo'
            if not os.path.isfile(mli_geo):
                geocode_back(mli, lookup, mli_geo, width_rdc, width, nlines, 0)
            mli = mli_geo

        pwr = read_gamma(mli, nlines, 'float32')
        scaled_pwr = scale_pwr(pwr)

    lon_lat = get_lon_lat(dem_seg_par)

    os.chdir(out_dir)

    if wrap_range:
        for r in wrap_range:
            wrapped_data = wrap_data(data, [-r / 2, r / 2])
            r = str(r)
            if r.endswith('.0'):
                r = r[0:-2]
            img_name = f"image_{r}.png"
            cbar_name = f"colorbar_{r}.png"
            vlim = [np.nanmin(wrapped_data), np.nanmax(wrapped_data)]
            draw_colorbar(vlim, cmap=cmap, label=None, out_file=cbar_name)
            if mli:
                draw_img_with_pwr(wrapped_data,
                                  cmap=cmap,
                                  pwr=scaled_pwr,
                                  vlim=vlim,
                                  out_file=img_name,
                                  dpi=dpi)
            else:
                draw_img_with_pwr(wrapped_data,
                                  cmap=cmap,
                                  pwr=None,
                                  vlim=vlim,
                                  out_file=img_name,
                                  dpi=dpi)
            kmz_name = f"{os.path.basename(main_file)}_wrap_{r}.kmz"
            make_kmz(img_name, cbar_name, lon_lat, kmz_name, out_dir)
    else:
        img_name = "image.png"
        cbar_name = "colorbar.png"
        vlim = [np.nanmin(data), np.nanmax(data)]
        draw_colorbar(vlim, cmap=cmap, label=None, out_file=cbar_name)
        if mli:
            draw_img_with_pwr(data,
                              cmap=cmap,
                              pwr=scaled_pwr,
                              vlim=vlim,
                              out_file=img_name,
                              dpi=dpi)
        else:
            draw_img_with_pwr(data,
                              cmap=cmap,
                              pwr=None,
                              vlim=vlim,
                              out_file=img_name,
                              dpi=dpi)
        kmz_name = f"{os.path.basename(main_file)}.kmz"
        make_kmz(img_name, cbar_name, lon_lat, kmz_name, out_dir)

    print('\nAll done, enjoy it!\n')


if __name__ == "__main__":
    main()

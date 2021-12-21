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
import base64

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from lxml import etree
from pykml.factory import KML_ElementMaker as KML

EXAMPLE = r"""Example:
  python3 make_km_velz.py vels.txt vels.kmz
  python3 make_kmz_vel.py vels.txt vels.kmz -v -100 100 -c jet -s 0.6 -n f 
  # data format (the first column is not important)
  num1 lon1 lat1 vel1
  num2 lon2 lat2 vel2
  num3 lon3 lat3 vel3
  ...
"""

DOT_STR = b"""iVBORw0KGgoAAAANSUhEUgAAAEAAAABACAYAAACqaXHeAAAABGdBTUEAALGPC
/xhBQAAAAlwSFlzAAAPYQAAD2EBqD+naQAABNhJREFUeNrtmmlSGlEUhRvnGVxByA7YgZ0VJO6A
Kkv9KVlBcAEW4FTiUNEdsAM7KwjZgVlBcJ4l53T1JZeXlsR002h8r+rWa+CHnu+dN93bKeeVt5Q
FYAFYABaABWABWAAWgAVgAVgAFoAF0KU2NzeXQecickEvrRn8Dw1EPQhvb2+v8V8AgHCKLSDeDw
wMOP39/X709fU5qdSvP91sNp2Hhwfn/v7ej7u7uwN8vQ8Q3osEEAgvQWhuaGjIGRwcdAhAIFC8Q
BDx7APxftzc3Di3t7cefvsIEPUXASCwehHCloaHhx2KZ08AAuExAOIAiodwPwjh+vqaUQGEwrMG
EIg/hMjcyMiIwxAI2gUyDSheAOjRZ0/xFK4AOJeXl3X89i7u9SEVk3gubocQmRkbG3NGR0d98Rq
CuCAMgDH/2wBcXV1RvPR1/BYrhFRMI/8VQrMUziCETi6g/c0pIADE+iI+GP1WnJ+fxwohDgAceX
d8fNwXLg6QIADtApkGBCA7gJ4Ceu4TgDjg4uKiFYBQq1arsz0HAPGfIKZI0RMTE04YhDAXaABsd
IBe/QWAHnklXqKws7NT6RkAiM/S+hCbEfHsKV5AaBcQhN4N9BQw7S+jr2zf6s/OzqRv4Lu3UadC
FACfISY/OTnpCxcIOvR6IFPhTy3M+mrUWxAYp6enRUyF5cQBBAvfD4okAIY8ixPM6UDb/22jGyh
cRl0Lh2j/mT1dgOdILvhXAEuwcVnEiwtMN8iU0PP9KRAolCBkxAWA9EHkt7e3D5IGcIhRdTUAE4
YE5/u/Nq4JpviQqG1sbMwmDaCpBU9NTbU9ixto/aiNDpBRPzk5aRMffG5UKpXpxABA/Aws7VGoC
H+sf8q8f6xxa1Ri23oJfH67u7t7lBQAXm1rAsAEIUEXxNVM0SYAhIt14EtSAD5h+yum0+mWcP0s
ALj4xdW4GIaJPz4+ls/u1tZWcgCwnxf1yAsAHdz742qyDuhRV+L5bAEkCWAGU8B7zVNgBoug12k
RlC0xwUUwi0Xwe2LngPn5+WanbVAiyiFIH4ZEsD4LKPGNUqmU3DlAToIQ6IYdhHTEeRAyQ4Gora
6uJn4SXIK4Mvd6U7wcgeU5jqOwvgCFRB5H4cTvAmkIa4RdhvStUCLKZSjkCqzvBY2VlZXpKA6Ll
A+AuLx5+ZEwEyNPORbz+GtmgSQRogGgL66trS33CsAbbId1iM2YwqUX8TpB2qkxO8SMEMUzMWLk
AducwFwAIlutVo97mhOEyKJOgnTKCzIdxrygLo9JVthMiOqcoOQFjCmRX19fP4i6yEbOCmNLPIR
oV9veTI/rpKjOCptpcakICQCdFzScUCuXy73PCsuCyJMhhOcoWgMwawNmecwsjYWUxH7LDLMugN
6Nav24K0NvAghZXRPolBIPqw3qaaALIypB6ovHsTcW8bEB0E6A8JwefV0gNUvkui6gU+OPuMADi
A9xio8VgLTFxcUShBd0RYi9lMd1aUxXhnRxRAMInFDEaW/Z6ULr1vsBM7wyI1xdGNULoOmAsPI4
AHgAUNjc3PzmdKl19Q0R7BC8OucBIW/avwOABqIGAPsQ/sXpckvsJamFhQVeo10AyGHxyxgAjoL
wkhDdEwDPtVkAFoAFYAFYABaABWABWAAWgAXwOttPxi1EjO7EVwYAAAAASUVORK5CYII="""


def cmdline_parser():
    parser = argparse.ArgumentParser(description='Display velocity derived by InSAR in Google Earth.',
        formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)
    parser.add_argument('vel_file', help='formatted velocity file for making KMZ file')
    parser.add_argument('out_file', help='output KMZ file')
    parser.add_argument('-v', dest='vlim', nargs=2, metavar=('MIN', 'MAX'), type=float,
        default=(-60, 60), help='velocity limits (defaults: -60 60)')
    parser.add_argument('-c', dest='colormap', default='jet', help='colormap (defaults: jet)')
    parser.add_argument('-s', dest='scale', default=0.8, type=float, help='scale of point for display (defaults: 0.8)')
    parser.add_argument('-n', dest='num_flag', default='t', choices=['t', 'T', 'f', 'F'],
        help='first column of data is point number [t] or not [f] (defaults: t)')

    inps = parser.parse_args()
    return inps


def plot_colorbar(out_file, vlim, cmap, figsize=(0.36, 3.6), nbins=6):
    """Draw colorbar for making KMZ

    Args:
        out_file (str): output colorbar file
        vlim (list): value limit
        cmap (str): colormap
        figsize (tuple, optional): figure size. Defaults to (0.36, 3.6).
        nbins (int, optional): nbins. Defaults to 7.
    """
    fig, cax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(left=0.1, right=0.5)

    norm = mpl.colors.Normalize(vmin=vlim[0], vmax=vlim[1])
    cmap = plt.get_cmap(cmap)
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')

    label = 'mean LOS velocity [mm/yr]'
    cbar.set_label(label, fontsize=6)
    cbar.locator = mpl.ticker.MaxNLocator(nbins=nbins)
    cbar.update_ticks()
    cbar.ax.tick_params(which='both', labelsize=6)
    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(1)

    fig.savefig(out_file, bbox_inches='tight', facecolor=fig.get_facecolor(), dpi=300)


def load_vel(vel_file, num_flag):
    """Load velocity data

    Args:
        vel_file (str): vel file
        num_flag (str): number flag (t or f)

    Returns:
        tuple: velocity data
    """
    content = np.loadtxt(vel_file, np.float64)
    # add points number
    if num_flag == 'f':
        number = np.arange(-1, content.shape[0] - 1)
        content = np.hstack((number.reshape(-1, 1), content))

    return content


def get_description_string(lon, lat, vel):
    """Description information of each data point

    Args:
        lon (float): longitude
        lat (float): latitude
        vel (float): velocity

    Returns:
        str: description string
    """
    des_str = '<html xmlns:fo="http://www.w3.org/1999/XSL/Format" xmlns:msxsl="urn:schemas-microsoft-com:xslt">'
    des_str += '<head>'
    des_str += '<META http-equiv="Content-Type" content="text/html">'
    des_str += '<meta http-equiv="content-type" content="text/html; charset=UTF-8">'
    des_str += '<style>td{padding:5px;}</style>'
    des_str += '</head>'
    des_str += '<body style="margin:0px 0px 0px 0px;overflow:auto;background:#FFFFFF;">'
    des_str += '<table style="font-family:Arial,Verdana,Times;font-size:12px;text-align:left;width:100%;border-spacing:1px;padding:3px 3px 3px 3px;">'

    des_str += f'<tr bgcolor="#F5F5F5"><td>Longitude</td><td>{lon}</td></tr>'
    des_str += f'<tr bgcolor="#F5F5F5"><td>Latitude</td><td>{lat}</td></tr>'
    des_str += f'<tr bgcolor="#F5F5F5"><td>Mean LOS velocity</td><td>{vel}</td></tr>'
    des_str += '</table>'
    des_str += '</body>'
    des_str += '</html>'

    return des_str


def get_hex_color(v, colormap, norm):
    """Get color name in hex format

    Args:
        v (float): number of interest
        colormap (instance): matplotlib.colors.Colormap instance
        norm (instance): matplotlib.colors.Normalize instance

    Returns:
        str: color name in hex format
    """
    # get rgba color components for point velocity
    rgba = colormap(norm(v))
    c_hex = mpl.colors.to_hex([rgba[3], rgba[2], rgba[1], rgba[0]], keep_alpha=True)[1:]

    return c_hex


def make_kmz(vel_file, scale, num_flag, vlim, cmap, dot_str, out_file):
    """Make kmz

    Args:
        vel_file (str): timeseries file
        scale (float): point scale
        num_flag (str): number flag, t or f
        vlim (list): value limits
        cmap (str): colormap name
        dot_str (str): dot_str for draw dot
        out_file (str): output kmz file
    """
    print('Writing data to {}'.format(out_file))

    # get lons, lats, vels, ts, dates
    try:
        data = load_vel(vel_file, num_flag)
        nums, lons = data[:, 0], data[:, 1]
        lats, vels = data[:, 2], data[:, 3]
    except Exception:
        sys.exit('Error format of {}'.format(vel_file))

    # create document
    doc_name = os.path.basename(out_file)[:-4]
    doc = KML.kml(KML.Document(KML.Folder(KML.name(doc_name))))

    # set and normalize colormap to defined vlim
    colormap = mpl.cm.get_cmap(cmap)
    norm = mpl.colors.Normalize(vmin=vlim[0], vmax=vlim[1])

    # create placemark
    for j in range(lons.shape[0]):
        num, lon, lat, vel = nums[j], lons[j], lats[j], vels[j]
        description_str = get_description_string(lon, lat, vel)
        description = KML.description(description_str)
        name = KML.name(str(int(num)))

        color = get_hex_color(vel, colormap, norm)
        style = KML.Style(
            KML.IconStyle(KML.color(color), KML.scale(str(scale)),
                          KML.Icon(KML.href("shaded_dot.png"))),
            KML.LabelStyle(KML.color('00000000'), KML.scale('0.0')))
        point = KML.Point(KML.altitudeMode('clampToGround'),
                          KML.coordinates(f"{lon},{lat}"))
        placemark = KML.Placemark(name, style, description, point)
        doc.Document.Folder.append(placemark)

    # create legend
    legend = KML.ScreenOverlay(
        KML.name('colorbar'),
        KML.Icon(KML.href('colorbar.png')),
        KML.overlayXY(
            x="0.0",
            y="1",
            xunits="fraction",
            yunits="fraction",
        ),
        KML.screenXY(
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
        ),
        KML.size(
            x="0.08",
            y="0.48",
            xunits="fraction",
            yunits="fraction",
        ),
    )
    doc.Document.Folder.append(legend)

    kml_str = etree.tostring(doc, pretty_print=True)

    # write kml file
    dir_name = os.path.dirname(out_file)
    os.chdir(dir_name)
    kml_file = doc_name + '.kml'
    with open(kml_file, 'wb') as f:
        f.write(kml_str)

    # draw colorbar and dot
    cbar_file = 'colorbar.png'
    plot_colorbar(cbar_file, vlim, colormap)

    dot_file = 'shaded_dot.png'
    with open(dot_file, 'wb') as f:
        f.write(base64.b64decode(dot_str))

    # unzip kml, dot, colorbar, dygraph_file
    with zipfile.ZipFile(out_file, 'w') as f:
        f.write(kml_file)
        os.remove(kml_file)
        f.write(cbar_file)
        os.remove(cbar_file)
        f.write(dot_file)
        os.remove(dot_file)

    print('\nAll done, enjoy it!\n')


if __name__ == "__main__":
    inps = cmdline_parser()
    vel_file = os.path.abspath(inps.vel_file)
    out_file = os.path.abspath(inps.out_file)
    vlim = inps.vlim
    cmap = inps.colormap
    scale = inps.scale
    flag = inps.num_flag

    # check vel_file
    if not os.path.isfile(vel_file):
        sys.exit("{} does not exist".format(vel_file))

    # check out_file
    if not out_file.endswith('.kmz'):
        out_file += '.kmz'
    dir_name = os.path.dirname(out_file)
    if not os.path.isdir(dir_name):
        sys.exit("{} does not exist".format(dir_name))

    # check vlim
    if vlim[0] > vlim[1]:
        vlim = vlim[::-1]

    # check scale
    if scale <= 0:
        sys.exit('scale must bigger than 0')

    make_kmz(vel_file, scale, flag, vlim, cmap, DOT_STR, out_file)

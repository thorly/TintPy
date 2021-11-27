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
from multiprocessing import Process
from xml.dom.minidom import parse

import numpy as np

EXAMPLE = """Example:
  python3 clip_vel_ts.py vel.txt v cut.kml
  python3 clip_vel_ts.py vel.txt v cut_multi.kml -n f -a m -o clip_data
  python3 clip_vel_ts.py ts.txt t cut.kmz
  python3 clip_vel_ts.py ts.txt t cut_multi.kmz -n f -a m -o clip_data
"""


def cmdline_parser():
    parser = argparse.ArgumentParser(
        description="Clip velocity or timeseries data using kml or kmz",
        formatter_class=argparse.RawTextHelpFormatter,epilog=EXAMPLE)
    parser.add_argument('data_file', help='original velocity or timeseries data')
    parser.add_argument('data_flag', choices=['t', 'T', 'v', 'V'], help='data flag, t for timeseries, v for velocity')
    parser.add_argument('kml', help='kml or kmz file for clipping data')
    parser.add_argument('-n', dest='num_flag', choices=['t', 'T', 'f', 'F'],
        default='t', help='first column of data is point number [t] or not [f] (defaults: t)'
    )
    parser.add_argument('-a', dest='area_flag', default='s', choices=['s', 'S', 'm', 'M'],
        help='kml file contain single area [s] or multi area [m] (defaults: s)'
    )

    parser.add_argument('-o', dest='out_dir', default='.', help='output directory (defaults: .)')

    inps = parser.parse_args()
    return inps


def ray_intersect_segment(point, s_point, e_point):
    """Determine whether the ray intersects the line segment

    Args:
        point (list): point [lon, lat]
        s_point (list): start point [lon, lat]
        e_point (list): end point

    Returns:
        bool: True or False
    """
    # parallel and coincident with the ray，s_point coincides with s_point
    if s_point[1] == e_point[1]:
        return False
    # line segment is above the ray
    if s_point[1] > point[1] and e_point[1] > point[1]:
        return False
    # line segment under the ray
    if s_point[1] < point[1] and e_point[1] < point[1]:
        return False
    # point coincides with s_point
    if s_point[1] == point[1] and e_point[1] > point[1]:
        return False
    # point coincides with e_point
    if e_point[1] == point[1] and s_point[1] > point[1]:
        return False
    # line segment is to the left of the ray
    if s_point[0] < point[0] and e_point[0] < point[0]:
        return False
    # find the intersection
    xseg = e_point[0] - (e_point[0] - s_point[0]) * (e_point[1] - point[1]) / (e_point[1] - s_point[1])
    # intersection is to the left of point
    if xseg < point[0]:
        return False
    # exclude the above cases
    return True


def point_within_polygon(point, polygon):
    """Determine whether the point is within the polyline

    Args:
        point (list): point [lon, lat]
        polygon (list): polygon [[lon, lat], [lon, lat]...]

    Returns:
        bool: True or False
    """
    # number of intersection
    num = 0
    for i in range(len(polygon) - 1):
        if ray_intersect_segment(point, polygon[i], polygon[i + 1]):
            num += 1

    if num % 2 == 1:
        return True
    else:
        return False


def kml2polygon_dict(kml_file):
    """Parse kml or kmz file, get polygon

    Args:
        kml_file (str): kml or kmz file

    Returns:
        dict: polygon
    """
    # unzip kmz to get kml
    if kml_file.endswith('.kmz'):
        dir_name = os.path.dirname(kml_file)
        with zipfile.ZipFile(kml_file, 'r') as f:
            files = f.namelist()
            for file in files:
                if file.endswith('.kml'):
                    f.extract(file, dir_name)
        doc_kml = os.path.join(dir_name, file)
        kml_file = doc_kml

    dom_tree = parse(kml_file)

    if os.path.isfile(doc_kml):
        os.remove(doc_kml)

    # parse kml
    root_node = dom_tree.documentElement
    placemarks = root_node.getElementsByTagName('Placemark')

    polygon_dict = {}
    j = 0

    for placemark in placemarks:
        name = placemark.getElementsByTagName('name')[0].childNodes[0].data
        ploygon = placemark.getElementsByTagName('Polygon')[0]
        outerBoundaryIs = ploygon.getElementsByTagName('outerBoundaryIs')[0]
        LinearRing = outerBoundaryIs.getElementsByTagName('LinearRing')[0]
        coordinates = LinearRing.getElementsByTagName('coordinates')[0].childNodes[0].data
        lon_lat = [i.split(',')[0:2] for i in coordinates.strip().split()]
        polygon_dict[name + '_' + str(j)] = np.asarray(lon_lat, dtype='float32')
        j += 1

    return polygon_dict


def clip_ts_with_lonlat(data, lon_lat):
    """Clip timeseries using lon lat

    Args:
        data (array): timeseries data (first column is number)
        lon_lat (tuple): min max lon_lat

    Returns:
        array: clipped result
    """
    lon_min, lon_max, lat_min, lat_max = lon_lat
    lon_data = data[1:, 1]
    lat_data = data[1:, 2]

    lon_index = ((lon_data > lon_min) == (lon_data < lon_max))
    lat_index = ((lat_data > lat_min) == (lat_data < lat_max))
    index = (lon_index & lat_index)

    first_row = data[0, :]
    left_row = data[1:, :]
    clipped_data = np.vstack((first_row, left_row[index, :]))

    return clipped_data


def clip_vel_with_lonlat(data, lon_lat):
    """Clip velocity using lon lat

    Args:
        data (array): velocity data (first column is number)
        lon_lat (tuple): min max lon_lat

    Returns:
        array: clipped result
    """
    lon_min, lon_max, lat_min, lat_max = lon_lat
    lon_data = data[:, 1]
    lat_data = data[:, 2]

    lon_index = ((lon_data > lon_min) == (lon_data < lon_max))
    lat_index = ((lat_data > lat_min) == (lat_data < lat_max))
    index = (lon_index & lat_index)

    clipped_data = data[index, :]

    return clipped_data


def clip_data_with_lonlat(data, lon_lat, data_flag, num_flag):
    """Clip timeseries or velocity using lon lat

    Args:
        data (array): timeseries or velocity data
        lon_lat (tuple): min max lon_lat
        data_flag (str): t for timeseries or v for velocity
        num_flag (str): t for first column is number f for not

    Returns:
        array: clipped result (first column is number)
    """
    if num_flag == 't':
        if data_flag == 't':
            return clip_ts_with_lonlat(data, lon_lat)
        else:
            return clip_vel_with_lonlat(data, lon_lat)
    else:
        if data_flag == 't':
            # add number to first column
            number = np.arange(-1, data.shape[0] - 1).reshape((-1, 1))
            data = np.hstack((number, data))

            return clip_ts_with_lonlat(data, lon_lat)
        else:
            # add number to first column
            number = np.arange(0, data.shape[0]).reshape((-1, 1))
            data = np.hstack((number, data))

            return clip_vel_with_lonlat(data, lon_lat)


def get_extremum(polygon):
    """Get min max lon_lat from polygon

    Args:
        polygon (array): polygon

    Returns:
        tuple: min max lon_lat
    """
    lon_min = np.min(polygon[:, 0])
    lon_max = np.max(polygon[:, 0])
    lat_min = np.min(polygon[:, 1])
    lat_max = np.max(polygon[:, 1])

    return lon_min, lon_max, lat_min, lat_max


def clip_single_vel_with_kml(kml_file, vel_file, num_flag, out_dir):
    """Clip velocity data with kml which contains single area

    Args:
        kml_file (str): kml or kmz file
        vel_file (str): velocity file
        num_flag (str): t for first column is number f for not
        out_dir (str): output directory
    """
    polygon_dict = kml2polygon_dict(kml_file)
    print('Loading data')
    vel = np.loadtxt(vel_file)

    for _, polygon in polygon_dict.items():
        lon_lat = get_extremum(polygon)
        first_clipped_data = clip_data_with_lonlat(vel, lon_lat, 'v', num_flag)

        out_data = np.arange(vel.shape[1])

        for line in first_clipped_data:
            if point_within_polygon(line[1:3], polygon):
                if num_flag == 't':
                    out_data = np.vstack((out_data, line))
                else:
                    out_data = np.vstack((out_data, line[1:]))

        if out_data.size > vel.shape[1]:
            vel_file_name = os.path.splitext(os.path.basename(vel_file))[0]
            out_file = os.path.join(out_dir, vel_file_name + '_clip.txt')
            print(f'Writing data to {out_file}')
            np.savetxt(out_file, out_data[1:, :], fmt='%4f')
            print('\nAll done, enjoy it!\n')


def clip_single_ts_with_kml(kml_file, ts_file, num_flag, out_dir):
    """Clip timeseries data with kml which contains single area

    Args:
        kml_file (str): kml or kmz file
        ts_file (str): timeseries file
        num_flag (str): t for first column is number f for not
        out_dir (str): output directory
    """
    polygon_dict = kml2polygon_dict(kml_file)
    print('Loading data')
    data = np.loadtxt(ts_file)

    for _, polygon in polygon_dict.items():
        lon_lat = get_extremum(polygon)
        first_clipped_data = clip_data_with_lonlat(data, lon_lat, 't', num_flag)

        out_data = data[0, :]

        for line in first_clipped_data[1:, :]:
            if point_within_polygon(line[1:3], polygon):
                if num_flag == 't':
                    out_data = np.vstack((out_data, line))
                else:
                    out_data = np.vstack((out_data, line[1:]))

        if out_data.size > first_clipped_data[1:, :].shape[1]:
            ts_file_name = os.path.splitext(os.path.basename(ts_file))[0]
            out_file = os.path.join(out_dir, ts_file_name + '_clip.txt')
            print(f'Writing data to {out_file}')
            np.savetxt(out_file, out_data, fmt='%4f')
            print('\nAll done, enjoy it!\n')


def clip_multi_vel_with_kml(kml_file, vel_file, num_flag, out_dir):
    """Clip velocity data with kml which contains multi area

    Args:
        kml_file (str): kml or kmz file
        vel_file (str): velocity file
        num_flag (str): t for first column is number f for not
        out_dir (str): output directory
    """
    print('Loading data')
    vel = np.loadtxt(vel_file)
    polygon_dict = kml2polygon_dict(kml_file)

    num = len(polygon_dict)
    i = 0
    for name, polygon in polygon_dict.items():
        i += 1
        print(f'\rProcessing: {i}/{num}', end=" ", flush=True)

        out_data = np.arange(vel.shape[1])

        lon_lat = get_extremum(polygon)
        first_clipped_data = clip_data_with_lonlat(vel, lon_lat, 'v', num_flag)

        for line in first_clipped_data:
            if point_within_polygon(line[1:3], polygon):
                if num_flag == 't':
                    out_data = np.vstack((out_data, line))
                else:
                    out_data = np.vstack((out_data, line[1:]))

        if out_data.size > vel.shape[1]:
            vel_file_name = os.path.splitext(os.path.basename(vel_file))[0]
            out_file = os.path.join(out_dir, f'{vel_file_name}_{name}.txt')
            np.savetxt(out_file, out_data[1:, :], fmt='%4f')

    print(f'\rProcessed: {i}/{num}, enjoy it.', end=" ", flush=True)
    print('\n')


def clip_multi_ts_with_kml(kml_file, ts_file, num_flag, out_dir):
    """Clip timeseries data with kml which contains multi area

    Args:
        kml_file (str): kml or kmz file
        ts_file (str): timeseries file
        num_flag (str): t for first column is number f for not
        out_dir (str): output directory
    """
    print('Loading data')
    data = np.loadtxt(ts_file)
    polygon_dict = kml2polygon_dict(kml_file)

    num = len(polygon_dict)
    i = 0
    for name, polygon in polygon_dict.items():
        i += 1
        print(f'\rProcessing: {i}/{num}', end="", flush=True)

        out_data = data[0, :]

        lon_lat = get_extremum(polygon)
        first_clipped_data = clip_data_with_lonlat(data, lon_lat, 't', num_flag)

        for line in first_clipped_data[1:, :]:
            if point_within_polygon(line[1:3], polygon):
                if num_flag == 't':
                    out_data = np.vstack((out_data, line))
                else:
                    out_data = np.vstack((out_data, line[1:]))

        if out_data.size > first_clipped_data[1:, :].shape[1]:
            ts_file_name = os.path.splitext(os.path.basename(ts_file))[0]
            out_file = os.path.join(out_dir, f'{ts_file_name}_{name}.txt')
            np.savetxt(out_file, out_data, fmt='%4f')

    print(f'\rProcessed: {i}/{num}, enjoy it.', end=" ", flush=True)
    print('\n')


def main():
    inps = cmdline_parser()
    data_file = os.path.abspath(inps.data_file)
    kml = os.path.abspath(inps.kml)
    data_flag = inps.data_flag.lower()
    num_flag = inps.num_flag.lower()
    area_flag = inps.area_flag.lower()
    out_dir = os.path.abspath(inps.out_dir)

    # check data_file
    if not os.path.isfile(data_file):
        sys.exit('{} does not exist'.format(data_file))
    # check kml
    if not os.path.isfile(kml):
        sys.exit('{} does not exist'.format(kml))
    # check out_dir
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # cut timeseries data
    if data_flag == 't':
        if area_flag == 's':
            p = Process(target=clip_single_ts_with_kml, args=(kml, data_file, num_flag, out_dir))
        else:
            p = Process(target=clip_multi_ts_with_kml, args=(kml, data_file, num_flag, out_dir))
    # cut velocity data
    else:
        if area_flag == 's':
            p = Process(target=clip_single_vel_with_kml, args=(kml, data_file, num_flag, out_dir))
        else:
            p = Process(target=clip_multi_vel_with_kml, args=(kml, data_file, num_flag, out_dir))

    p.start()
    p.join()


if __name__ == "__main__":
    main()

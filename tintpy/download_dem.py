#!/usr/bin/env python3
################################
# Program is part of TintPy    #
# Copyright (c) 2021, Lei Yuan #
# Author: Lei Yuan, 2021       #
################################

import argparse
import math
import os
import sys
import urllib.request


def format_num(num, length, flag):
    """Get ALOS dem lon lat

    Args:
        num (int): number
        length (int): length of result
        flag (str): lon or lat

    Returns:
        str: result
    """
    num_str = str(abs(num))
    zero_num = length - len(num_str)
    if zero_num > 0:
        formated_num = '0' * zero_num + num_str
    else:
        formated_num = num_str

    if num < 0:
        if flag == 'lat':
            res = 'S' + formated_num
        if flag == 'lon':
            res = 'W' + formated_num
    else:
        if flag == 'lat':
            res = 'N' + formated_num
        if flag == 'lon':
            res = 'E' + formated_num

    return res


def alos_dem(s, n, w, e):
    """Get ALOS dem download urls and lon_lats

    Args:
        s (int): south latitude
        n (int): north latitude
        w (int): west longitude
        e (int): east longitude

    Returns:
        tuple: urls and lon_lats
    """
    HEADER = "https://www.eorc.jaxa.jp/ALOS/aw3d30/data/release_v1903/"
    urls = []
    lon_lats = []

    # calculate latitude and longitude, it must be a multiple of 5
    lon_min = int(w // 5 * 5)
    lon_max = math.ceil(e / 5) * 5
    lat_min = int(s // 5 * 5)
    lat_max = math.ceil(n / 5) * 5

    # traverse to get download urls
    for i in range(lat_min, lat_max, 5):
        for j in range(lon_min, lon_max, 5):
            i_5, j_5 = i + 5, j + 5

            i_0 = format_num(i, 3, 'lat')
            j_0 = format_num(j, 3, 'lon')
            i_5 = format_num(i_5, 3, 'lat')
            j_5 = format_num(j_5, 3, 'lon')

            lon_lat = "({}° ~ {}° , {}° ~ {}°)".format(j_0, j_5, i_0, i_5)
            lon_lats.append(lon_lat)

            name = "{}{}_{}{}.tar.gz".format(i_0, j_0, i_5, j_5)
            url = HEADER + name
            urls.append(url)

    return lon_lats, urls


def srtm_dem9055(s, n, w, e):
    """Get SRTM 90m dem (5*5 degrees) download urls and lon_lats

    Args:
        s (int): south latitude
        n (int): north latitude
        w (int): west longitude
        e (int): east longitude

    Returns:
        tuple: urls and lon_lats
    """
    HEADER = "http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/srtm_"

    urls = []
    lon_lats = []

    # latitude and longitude of starting point for calculating dem number
    lon_min = -180
    lat_max = 60

    # calculate dem number
    num_min_lon = (w - lon_min) / 5 + 1
    num_max_lon = (e - lon_min) / 5 + 1
    num_min_lat = (lat_max - n) / 5 + 1
    num_max_lat = (lat_max - s) / 5 + 1

    if num_min_lon > int(num_min_lon):
        num_min_lon = int(num_min_lon)

    if num_max_lon > int(num_max_lon):
        num_max_lon = int(num_max_lon + 1)

    if num_min_lat > int(num_min_lat):
        num_min_lat = int(num_min_lat)

    if num_max_lat > int(num_max_lat):
        num_max_lat = int(num_max_lat + 1)

    # traverse to get download urls
    for i in range(int(num_min_lon), int(num_max_lon)):
        for j in range(int(num_min_lat), int(num_max_lat)):
            # calculate the longitude and latitude range of dem
            w = i * 5 - 180
            e = w - 5
            s = 60 - j * 5
            n = s + 5

            # add 0
            w = "0" + str(w) if len(str(w)) == 1 else str(w)
            e = "0" + str(e) if len(str(e)) == 1 else str(e)
            s = "0" + str(s) if len(str(s)) == 1 else str(s)
            n = "0" + str(n) if len(str(n)) == 1 else str(n)

            lon_lat = "({}° ~ {}° , {}° ~ {}°)".format(e, w, s, n)
            lon_lats.append(lon_lat)

            num_lon = str(i)
            num_lat = str(j)
            if len(num_lon) == 1:
                num_lon = "0" + num_lon
            if len(num_lat) == 1:
                num_lat = "0" + num_lat

            url = "{}{}_{}.zip".format(HEADER, num_lon, num_lat)
            urls.append(url)

    return lon_lats, urls


def srtm_dem11(s, n, w, e, flag):
    """Get SRTM dem (1*1 degrees) download urls and lon_lats

    Args:
        s (int): south latitude
        n (int): north latitude
        w (int): west longitude
        e (int): east longitude
        flag (str): srtm30 or srtm90

    Returns:
        tuple: urls and lon_lats
    """
    number = 1
    if flag.upper() == 'SRTM90':
        number = 3
    HEADER = f"https://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL{number}.003/2000.02.11/"

    urls = []
    lon_lats = []

    lon_min = math.floor(w)
    lon_max = math.ceil(e)
    lat_min = math.floor(s)
    lat_max = math.ceil(n)

    def add_zero(i, flag):
        s = str(abs(i))
        if flag == 'SN':
            if len(s) == 1:
                return '0' + s
            else:
                return s
        if flag == 'WE':
            if len(s) == 1:
                return '00' + s
            elif len(s) == 2:
                return '0' + s
            else:
                return s

    for i in range(lat_min, lat_max):
        for j in range(lon_min, lon_max):
            lon_lat = "({}° ~ {}° , {}° ~ {}°)".format(j, j + 1, i, i + 1)
            lon_lats.append(lon_lat)

            lon = add_zero(i, 'SN')
            lat = add_zero(j, 'WE')

            if i >= 0:
                if j >= 0:
                    name = "N{}E{}.SRTMGL{}.hgt.zip".format(lon, lat, number)
                else:
                    name = "N{}W{}.SRTMGL{}.hgt.zip".format(lon, lat, number)
            else:
                if j >= 0:
                    name = "S{}E{}.SRTMGL{}.hgt.zip".format(lon, lat, number)
                else:
                    name = "S{}W{}.SRTMGL{}.hgt.zip".format(lon, lat, number)

            urls.append(HEADER + name)

    return lon_lats, urls


def get_urls(flag, bound):
    """Get dem download urls

    Args:
        flag (str): dem flag
        bound (int list): lon and lat

    Returns:
        str list: dem urls
    """
    lon_lats = []
    urls = []

    w, e, s, n = bound[0], bound[1], bound[2], bound[3]

    if w >= e or s >= n:
        sys.exit('Error bound, please reset it (W E S N)!')

    if 'SRTM' in flag:
        if s < -60 or n > 60:
            sys.exit(
                'Error bound, latitude of SRTM dem must be in [-60, 60]!')

        if flag == 'SRTM9011':
            lon_lats, urls = srtm_dem11(s, n, w, e, flag)
        elif flag == 'SRTM3011':
            lon_lats, urls = srtm_dem11(s, n, w, e, flag)
        else:
            lon_lats, urls = srtm_dem9055(s, n, w, e)
    else:
        lon_lats, urls = alos_dem(s, n, w, e)

    print("\nDownloading urls of all DEMs:")
    for i, j in zip(lon_lats, urls):
        if lon_lats.index(i) == len(lon_lats) - 1:
            print("{}: {}\n".format(i, j))
        else:
            print("{}: {}".format(i, j))

    return urls


def download_progress(blocknum, blocksize, totalsize):
    """Print download progress

    Args:
        blocknum (int): download block number
        blocksize (int): block size
        totalsize (int): file size
    """
    percent = 100.0 * blocknum * blocksize / totalsize
    if percent > 100:
        percent = 100
        print("\rDownloaded: " + "#" * int(percent / 2) + " %.2f%%" % percent, end="\n", flush=True)
    else:
        print("\rDownloading: " + "#" * int(percent / 2) + " %.2f%%" % percent, end=" ", flush=True)


def download_dem(url, out_dir):
    """download dem

    Args:
        url (str): dem url
        out_dir (str): output directory
    """
    file = os.path.join(out_dir, url.split('/')[-1])
    try:
        urllib.request.urlretrieve(url, file, download_progress)
    except Exception as e:
        print(e)


EXAMPLE = '''Example:
  # only get urls
  python3 download_dem.py srtm3011 100 101 20 24

  # get urls and download them (for srtm3011 and srtm9011, you must install wget)
  python3 download_dem.py srtm9055 100 101 20 24 -o /ly/dem
'''


def cmdline_parser():
    parser = argparse.ArgumentParser(description='Download SRTM DEM [30m 90m] and ALOS DSM [30m].',
        formatter_class=argparse.RawTextHelpFormatter,epilog=EXAMPLE)
    parser.add_argument('flag',type=str,
        choices=['alos', 'ALOS', 'srtm9011', 'srtm3011', 'SRTM9011', 'SRTM3011', 'srtm9055', 'SRTM9055'], help='DEM type')
    parser.add_argument('bound', type=float, nargs=4,help='DEM bound (W E S N)')
    parser.add_argument('-o',dest='out_dir', type=str, help='directory for saving DEM')

    inps = parser.parse_args()

    return inps


def main():
    # get inputs
    inps = cmdline_parser()
    flag = inps.flag.upper()
    bound = inps.bound
    out_dir = inps.out_dir

    # get urls of DEM
    urls = get_urls(flag, bound)

    # download DEM
    if out_dir:
        out_dir = os.path.abspath(out_dir)
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        if flag in ['ALOS', 'SRTM9055']:
            for url in urls:
                print(f"\nStart to download {url.split('/')[-1]}")
                download_dem(url, out_dir)

        if flag in ['SRTM3011', 'SRTM9011']:
            os.chdir(out_dir)
            for url in urls:
                print(f"\nStart to download {url.split('/')[-1]}")
                cmd_str = f"wget {url} --user=leiyuan --password=PVmg2NeCSLatf3v"
                os.system(cmd_str)


if __name__ == "__main__":
    main()

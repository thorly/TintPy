#!/usr/bin/env python3
################################
# Program is part of TintPy    #
# Copyright (c) 2021, Lei Yuan #
# Author: Lei Yuan, 2021       #
################################

# Modified from https://github.com/isce-framework/isce2


import requests
import re
import os
import argparse
import datetime
from html.parser import HTMLParser
import sys
import glob

server = 'https://scihub.copernicus.eu/gnss/'

orbitMap = [('precise', 'AUX_POEORB'), ('restituted', 'AUX_RESORB')]

datefmt = "%Y%m%dT%H%M%S"
queryfmt = "%Y-%m-%d"
queryfmt2 = "%Y/%m/%d/"

# Generic credentials to query and download orbit files
credentials = ('gnssguest', 'gnssguest')

EXAMPLE = """Example:
  python3 fetchOrbit.py /ly/sentinel.txt /ly/orbits
"""


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(
        'Download Sentinel-1 A/B precise orbits',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLE)
    parser.add_argument(
        'input',
        type=str,
        help='Path of text or directory for getting images names')
    parser.add_argument('outdir',
                        type=str,
                        default='.',
                        help='Path to output directory')

    return parser.parse_args()


def FileToTimeStamp(safename):
    '''
    Return timestamp from SAFE name.
    '''
    safename = os.path.basename(safename)
    fields = safename.split('_')
    sstamp = [
    ]  # sstamp for getting SAFE file start time, not needed for orbit file timestamps

    try:
        tstamp = datetime.datetime.strptime(fields[-4], datefmt)
        sstamp = datetime.datetime.strptime(fields[-5], datefmt)
    except:
        p = re.compile(r'(?<=_)\d{8}')
        dt2 = p.search(safename).group()
        tstamp = datetime.datetime.strptime(dt2, '%Y%m%d')

    satName = fields[0]

    return tstamp, satName, sstamp


class MyHTMLParser(HTMLParser):
    def __init__(self, url):
        HTMLParser.__init__(self)
        self.fileList = []
        self._url = url

    def handle_starttag(self, tag, attrs):
        for name, val in attrs:
            if name == 'href':
                if val.startswith("https://scihub.copernicus.eu/gnss/odata") and val.endswith(")/"):
                    pass
                else:
                    downloadLink = val.strip()
                    downloadLink = downloadLink.split("/Products('Quicklook')")
                    downloadLink = downloadLink[0] + downloadLink[-1]
                    self._url = downloadLink

    def handle_data(self, data):
        if data.startswith("S1") and data.endswith(".EOF"):
            self.fileList.append((self._url, data.strip()))


def download_file(url, outdir='.', session=None):
    '''
    Download file to specified directory.
    '''

    if session is None:
        session = requests.session()

    path = outdir
    print('Downloading URL: ', url)
    request = session.get(url, stream=True, verify=True, auth=credentials)

    try:
        val = request.raise_for_status()
        success = True
    except:
        success = False

    if success:
        with open(path, 'wb') as f:
            for chunk in request.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
                    f.flush()

    return success


def fileToRange(fname):
    '''
    Derive datetime range from orbit file name.
    '''

    fields = os.path.basename(fname).split('_')
    start = datetime.datetime.strptime(fields[-2][1:16], datefmt)
    stop = datetime.datetime.strptime(fields[-1][:15], datefmt)
    mission = fields[0]

    return (start, stop, mission)


def get_sentinel1_from_zip(zip_dir):
    """Get Sentinel-1 names from directory

    Args:
        zip_dir (str): directory including Sentinel-1 zips

    Returns:
        list: Sentinel-1 names
    """
    image_name = []
    date_mission = []
    files = os.listdir(zip_dir)
    for file in files:
        if re.search(r'S1\w{65}\.zip', file):
            dm = re.findall(r"\d{8}", file)[0] + file[0:3]
            if dm not in date_mission:
                date_mission.append(dm)
                image_name.append(file)

    return sorted(image_name)


def get_sentinel1_from_text(txt_file):
    """Get Sentinel-1 names from file

    Args:
        txt_file (str): file including Sentinel-1 names

    Returns:
        list: Sentinel-1 names
    """
    image_name = []
    date_mission = []
    with open(txt_file, encoding='utf-8') as file:
        content = file.read()
        names = re.findall(r"S1\w{65}", content)
        if names:
            for name in names:
                dm = re.findall(r"\d{8}", name)[0] + name[0:3]
                if dm not in date_mission:
                    date_mission.append(dm)
                    image_name.append(name)

    return sorted(image_name)


def get_sentinel1(path):
    """Get Sentinel-1 names from file or directory

    Args:
        path (str): file or directory

    Returns:
        list: Sentinel-1 names
    """
    if os.path.isdir(path):
        return get_sentinel1_from_zip(path)
    else:
        return get_sentinel1_from_text(path)


def check_orbit(orb_file):
    """Check if the file is complete

    Args:
        orb_file (str): orbit file

    Returns:
        bool: True for complete, False for not
    """
    with open(orb_file, 'r', encoding='utf-8') as f:
        content = f.read()
        if '</Earth_Explorer_File>' in content and '<Earth_Explorer_File>' in content:
            return True
        else:
            return False


def date_math(date, day):
    """date operation (add or subtract day)

    Args:
        date (str): date
        day (int): day

    Returns:
        str: operated date
    """
    y = int(date[0:4])
    m = int(date[4:6])
    d = int(date[6:8])
    in_date = datetime.datetime(y, m, d)
    out_date = in_date + datetime.timedelta(days=day)
    out_date = out_date.strftime('%Y%m%d')

    return out_date


def get_left(orb_dir, sentinel_names):
    """Get no-orbit Sentinel-1 names

    Args:
        orb_dir (str): orbit directory
        sentinel_names (list): Sentinel-1 names

    Returns:
        list: no-orbit Sentinel-1 names
    """
    no_orb = sentinel_names.copy()
    orbits = glob.glob(os.path.join(orb_dir, 'S1*.EOF'))

    for name in sentinel_names:
        sensor1 = name[:3]
        date1 = re.findall(r'\d{8}', name)[0]
        dm1 = sensor1 + date1

        for orb in orbits:
            orb_name = os.path.basename(orb)
            sensor2 = orb_name[:3]
            date2 = re.findall(r'\d{8}', orb_name)[-2]
            date2 = date_math(date2, 1)
            dm2 = sensor2 + date2

            if dm1 == dm2 and check_orbit(orb):
                no_orb.remove(name)

    return no_orb


if __name__ == '__main__':
    '''
    Main driver.
    '''

    inps = cmdLineParse()
    input = os.path.abspath(inps.input)
    outdir = os.path.abspath(inps.outdir)

    if not os.path.exists(input):
        sys.exit('{} does not exist'.format(input))

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    sentinel_names = get_sentinel1(input)

    sentinel_names = get_left(outdir, sentinel_names)

    if not sentinel_names:
        print('No orbit need to download.')

    for name in sentinel_names:
        fileTS, satName, fileTSStart = FileToTimeStamp(name)
        print('Task: ' + str(sentinel_names.index(name) + 1) + '/' + str(len(sentinel_names)))
        print('Reference time: ', fileTS)
        print('Satellite name: ', satName)
        match = None
        session = requests.Session()

        for spec in orbitMap:
            oType = spec[0]
            delta = datetime.timedelta(days=1)
            timebef = (fileTS - delta).strftime(queryfmt)
            timeaft = (fileTS + delta).strftime(queryfmt)
            url = server + 'search?q=( beginPosition:[{0}T00:00:00.000Z TO {1}T23:59:59.999Z] AND endPosition:[{0}T00:00:00.000Z TO {1}T23:59:59.999Z] ) AND ( (platformname:Sentinel-1 AND filename:{2}_* AND producttype:{3}))&start=0&rows=100'.format(
                timebef, timeaft, satName, spec[1])

            success = False
            match = None

            try:
                r = session.get(url, verify=True, auth=credentials)
                r.raise_for_status()
                parser = MyHTMLParser(url)
                parser.feed(r.text)
                for resulturl, result in parser.fileList:
                    tbef, taft, mission = fileToRange(os.path.basename(result))
                    if (tbef <= fileTSStart) and (taft >= fileTS):
                        matchFileName = result
                        match = resulturl

                if match is not None:
                    success = True
            except:
                pass

            if success:
                break

        if match is not None:
            output = os.path.join(outdir, matchFileName)
            res = download_file(match, output, session)
            if res is False:
                print('Failed to download URL: ', match)
        else:
            print('Failed to find {1} orbits for tref {0}'.format(fileTS, satName))

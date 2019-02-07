#!/usr/bin/env python

from zipfile import ZipFile
from lxml import html
import os
import concurrent.futures
import time
from tqdm import tqdm
from gnssr.targets import *

def main():
    print('Searching oil platform')

    # Search all files
    kmz_files = []
    for root, subdirs, files in os.walk(os.path.join(os.environ['TDS_ROOT'], 'raw/L1B_Catalogue')):
        for file in files:
            if '.kmz' in file:
                filename = os.path.join(root,file)
                kmz_files.append(filename)

    # Single execution
    #for file in kmz_files:
    #    print(file)
    #    catalog_search(file)

    # Parallel execution
    catalog_path = os.path.join(os.environ['TDS_ROOT'],'search_target/catalog_search_output.txt')
    open(catalog_path, 'wt').close() # Clear file
    progress_bar = tqdm(total=len(kmz_files), unit_scale=False, unit='Files')
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
        for results in executor.map(catalog_search, kmz_files, chunksize=1):
            progress_bar.update()
            with open(catalog_path, 'at') as f:
                for line in results:
                    f.write(line);
                    print(line)
    print('\n%s: Download finished\n' % time.ctime())

def catalog_search(filename):
    target = targets['atlantis_pq']
    # 0.5 deg error approx 55 km error
    search_error = 0.3

    results=[]
    lat = 0
    lon = 0
    saved_file=False;
    kmz = ZipFile(filename, 'r')
    for kml_name in kmz.namelist():
        if 'doc.kml' in kml_name:
            continue
        kml = kmz.open(kml_name, 'r').read()
        doc = html.fromstring(kml)
        for pm in doc.cssselect('Document Placemark'):
            tmp = pm.cssselect('track')
            if len(tmp):
                # Track Placemark
                tmp = tmp[0]  # always one element by definition
                skip = 20 
                skip_count = 0
                for desc in tmp.iterdescendants():
                    # Skip some placemarks to speed up
                    if skip_count < skip:
                        skip_count = skip_count + 1
                        continue
                    skip_count = 0
                    content = desc.text_content()
                    if desc.tag == 'coord':
                        lon = float(content.split()[0])
                        lat = float(content.split()[1])
                        if (abs((lat%360) - (target.lat%360)) <= search_error) and (abs((lon%360) - (target.lon%360)) <= search_error):
                            if not saved_file:
                                results.append('\nFile: ' + filename + '\n')
                                saved_file=True
                            results.append(content + '\n')
    return results

if __name__ == '__main__':
    main()

#!/usr/bin/env python

from zipfile import ZipFile
from lxml import html
import os
import concurrent.futures

def process_file(filename):
    results=[]
    try:
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
                    for desc in tmp.iterdescendants():
                        content = desc.text_content()
                        if desc.tag == 'coord':
                            lon = float(content.split()[0])
                            lat = float(content.split()[1])
                            search_lat = 29
                            search_lon = -87
                            search_error = 2
                            if (abs((lat%360) - (search_lat%360)) <= search_error and abs((lon%360) - (search_lon%360)) <= search_error):
                                if not saved_file:
                                    results.append('\nFile: ' + filename + '\n')
                                    saved_file=True
                                results.append(content + '\n')
    except:
        pass
    return results

def main():
    # Search all files
    kmz_files = []
    for root, subdirs, files in os.walk('../raw/L1B_Catalogue'):
        for file in files:
            if '.kmz' in file:
                filename = os.path.join(root,file)
                kmz_files.append(filename)

    # Parallel execution
    f=open('kmz_search_output.txt', 'wt')
    count = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
        for results in executor.map(process_file, kmz_files, chunksize=1):
            print("{} / {}".format(count,len(kmz_files)))
            count = count + 1
            for line in results:
                f.write(line);
            pass
    f.close()
    print('Done')

if __name__ == '__main__':
    main()

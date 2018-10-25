#!/usr/bin/env python
import ftplib
from tqdm import tqdm
import time
import os

def download_file(remote_file):
    root_dir = os.path.join(os.environ['TDS_ROOT'],'raw/L1B/')
    local_file = os.path.join(root_dir, remote_file.replace('/Data/L1B/','').replace('/','-'))
    os.makedirs(os.path.dirname(local_file), exist_ok=True)
    try:
        if not os.path.exists(local_file) or os.path.getsize(local_file) != ftp.size(remote_file):
                    with open(local_file, 'wb') as file:
                        print('\n%s: Begin to download: %s\n' % (time.ctime(), remote_file))
                        progress_bar = tqdm(total=ftp.size(remote_file), unit_scale=True, unit='bit')
                        def update_download(data):
                            file.write(data)
                            progress_bar.update(len(data))
                        ftp.retrbinary('RETR ' + remote_file, update_download)
                        print('\n%s: Download finished\n' % time.ctime())
    except ftplib.error_perm:
        return

ip = 'ftp.merrbys.co.uk'
login = 'jblasco'
password = '@jYhr=M4rF'
with ftplib.FTP(ip, login, password) as ftp:
    search_list=[]
    with open(os.path.join(os.environ['TDS_ROOT'],'lat_lon_search/catalog_search_output.txt')) as kmz_search_file:
        for line in kmz_search_file:
            if 'File' in line:
                search_list.append(line);

    for count, line in enumerate(search_list):
        print('\nFile %s of %s\n' % (count,len(search_list)))
        file_string = line.split(' ')[1].strip()
        ddm_remote_file = os.path.join('/Data/L1B/', os.path.basename(file_string).replace('.','/').replace('kmz','DDMs.nc'))
        metadata_remote_file = os.path.join('/Data/L1B/', os.path.basename(file_string).replace('.','/').replace('kmz','metadata.nc'))
        download_file(ddm_remote_file)
        download_file(metadata_remote_file)
        if count == 5:
            break

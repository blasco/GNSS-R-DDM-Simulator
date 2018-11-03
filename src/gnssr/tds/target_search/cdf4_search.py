#!/usr/bin/env python

from netCDF4 import Dataset
import os
from gnssr.math_lib.utils import *
from gnssr.targets import *

file_root_name = 'raw/L1B/2017-11-13-H18'
target = targets['petronius']
# 0.5 deg error approx 55 km error
search_error = 0.25

def search_lat_lon(filename, search_lat, search_lon):
    metagrp = Dataset(filename, "r", format="NETCDF4")
    print(filename)
    search_output = filename + '\n' 
    for group in metagrp.groups:
        print(group + '/' + str(len(metagrp.groups)))
        #for index, datenum in enumerate(metagrp.groups[group].variables['IntegrationMidPointTime']):
        step = 1
        for index in range(0,len(metagrp.groups[group].variables['IntegrationMidPointTime']),step):
            datenum = metagrp.groups[group].variables['IntegrationMidPointTime'][index]
            lat = metagrp.groups[group].variables['SpecularPointLat'][index]
            lon = metagrp.groups[group].variables['SpecularPointLon'][index]
            # 0.5 deg error approx 55 km error
            if (abs((lat%360) - (search_lat%360)) <= search_error and abs((lon%360) - (search_lon%360)) <= search_error):
                string = 'G: ' + group + ' I: ' + str(index) + ' - ' + \
                        str(datenum) + ' - ' + str(datenum_to_pytime(float(datenum))) + ' - Lat: ' + \
                        str(lat) + ' Lon: ' + str(lon) + '\n'
                print(string)
                search_output = search_output + string
    return search_output

with open(os.path.join(os.environ['TDS_ROOT'],'lat_lon_search/cdf4_output.txt'),'wt') as file:
    filename = os.path.join(os.environ['TDS_ROOT'], file_root_name+'-metadata.nc')
    file.write(search_lat_lon(filename, search_lat, search_lon))

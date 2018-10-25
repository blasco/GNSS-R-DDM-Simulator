#!/usr/bin/env python

from netCDF4 import Dataset
from datetime import *
from dateutil import tz
import os

def datenum_to_pytime(matlab_datenum):
    python_datetime = datetime.fromordinal(int(matlab_datenum)) + timedelta(days=matlab_datenum%1) - timedelta(days = 366)
    return python_datetime

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
            search_error = 0.5
            if (abs((lat%360) - (search_lat%360)) <= search_error and abs((lon%360) - (search_lon%360)) <= search_error):
                string = 'G: ' + group + ' I: ' + str(index) + ' - ' + \
                        str(datenum) + ' - ' + str(datenum_to_pytime(float(datenum))) + ' - Lat: ' + \
                        str(lat) + ' Lon: ' + str(lon) + '\n'
                print(string)
                search_output = search_output + string
    return search_output

with open(os.path.join(os.environ['TDS_ROOT'],'lat_lon_search/cdf4_output.txt'),'wt') as file:
    #filename = os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2017-11-11-H18-metadata.nc')

    # Di Simone
    #filename = os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2017-11-07-H06-metadata.nc')
    filename = os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2017-11-07-H06-metadata.nc')
    search_lat = 46.75009
    search_lon = -48.78161

    # Petronius Oil Platform
    #filename = os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2017-11-26-H06-metadata.nc')
    #search_lat = 29.10499958
    #search_lon = -87.93832958

    # Atlantis PQ
    #filename = os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2017-11-08-H18-metadata.nc')
    #search_lat = 27.195278 
    #search_lon = -90.026944

    # Devils's Tower
    #filename = os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2018-07-30-H18-metadata.nc')
    #search_lat = 28.19013
    #search_lon = -88.49552

    # Songa Mercur
    #filename = os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2017-11-26-H18-metadata.nc')
    #search_lat = 8.375064
    #search_lon = 108.706184

    file.write(search_lat_lon(filename, search_lat, search_lon))

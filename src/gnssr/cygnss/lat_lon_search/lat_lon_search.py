#!/usr/bin/env python

from netCDF4 import Dataset
from datetime import *
from dateutil import tz
import os

def search_lat_lon(filename, search_lat, search_lon):
    root = Dataset(filename, "r", format="NETCDF4")
    print(filename)
    search_output = filename + '\n' 
    time_coverage_start = root.time_coverage_start
    step = 5
    for index in range(0,len(root.variables['ddm_timestamp_utc']),step):
        print(str(index) + '/' + str(len(root.variables['ddm_timestamp_utc'])))
        ddm_timestamp_utc = root.variables['ddm_timestamp_utc'][index]
        #print(str(root.variables['sp_lat'][index][0]) + ' : ' +  str(root.variables['sp_lon'][index][0]) )
        for channel in range(4):
            lat = root.variables['sp_lat'][index][channel]
            lon = root.variables['sp_lon'][index][channel]
            search_error = 1
            if (abs((lat%360) - (search_lat%360)) <= search_error) and (abs((lon%360) - (search_lon%360)) <= search_error):
                string = 'C: ' + str(channel) + ' I: ' + str(index) + ' - ' + str(time_coverage_start) + ' - ' + \
                        'Lat: ' + str(lat) + ' Lon: ' + str(lon) + '\n'
                print(string)
                search_output = search_output + string
    return search_output

with open('lat_lon_search_output.txt','wt') as file:
    filename = os.path.join(os.environ['CYGNSS_ROOT'], 'raw/cyg08.ddmi.s20180813-000000-e20180813-235959.l1.power-brcs.a20.d20.nc')
    # Petronius oil platform 29.10499958 -87.93832958
    #search_lat = 29.10499958
    #search_lon = -87.93832958

    # Devils's Tower
    search_lat = 28.125290
    search_lon = -88.442520

    # Songa Mercur
    #search_lat = 8.375064
    #search_lon = 108.706184

    file.write(search_lat_lon(filename, search_lat, search_lon))

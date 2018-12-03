#!/usr/bin/env python

from netCDF4 import Dataset
import os
from gnssr.utils import *
from gnssr.targets import *

# Example usage
def main():
    file_root_name = 'raw/L1B/2017-11-13-H18'
    target = targets['petronius']
    # 0.5 deg error approx 55 km error
    search_error = 0.5
    with open(os.path.join(os.environ['TDS_ROOT'],'search_target/cdf4_output.txt'),'wt') as file:
        file.write(cdf4_search(filen_root_name, target, search_error))

def cdf4_search(file_root_name, target, search_error):
    filename = os.path.join(os.environ['TDS_ROOT'], file_root_name+'-metadata.nc')
    metagrp = Dataset(filename, "r", format="NETCDF4")
    print(filename)
    search_output = filename + '\n' 
    for group in metagrp.groups:
        print(group + '/' + str(len(metagrp.groups)))
        #for index, datenum in enumerate(metagrp.groups[group].variables['IntegrationMidPointTime']):
        step = 3 
        for index in range(0,len(metagrp.groups[group].variables['IntegrationMidPointTime']),step):
            datenum = metagrp.groups[group].variables['IntegrationMidPointTime'][index]
            lat = metagrp.groups[group].variables['SpecularPointLat'][index]
            lon = metagrp.groups[group].variables['SpecularPointLon'][index]
            # 0.5 deg error approx 55 km error
            if (abs((lat%360) - (target.lat%360)) <= search_error and abs((lon%360) - (target.lon%360)) <= search_error):
                string = 'G: ' + group + ' I: ' + str(index) + ' - ' + \
                        str(datenum) + ' - ' + str(datenum_to_pytime(float(datenum))) + ' - Lat: ' + \
                        str(lat) + ' Lon: ' + str(lon) + '\n'
                print(string)
                search_output = search_output + string
    return search_output

if __name__ == '__main__':
    main()

#!/usr/bin/env python

import os
import numpy as np

from gnssr.tds.tds_data import *
from gnssr.utils import *
from gnssr.targets import *

def main():
    file_root_name = 'raw/L1B/2017-11-13-H18'
    target = targets['petronius']
    # 0.5 deg error approx 55 km error
    search_error = 10*np.pi/180
    with open(os.path.join(os.environ['TDS_ROOT'],'search_target/cdf4_geometry_output.txt'),'wt') as file:
        file.write(cdf4_search(file_root_name, target, search_error))

def cdf4_search(file_root_name, target, search_error):
    filename = os.path.join(os.environ['TDS_ROOT'], file_root_name+'-metadata.nc')
    print(filename)
    search_output = filename + '\n' 
    tds = tds_data(file_root_name)
    group_len = len(tds.metagrp.groups)
    for group in tds.metagrp.groups:
        step = 10
        index_len = len(tds.metagrp.groups[group].variables['IntegrationMidPointTime'])
        for index in range(0,index_len,step):
            print(str(group) + '/' + str(group_len) + ' : ' + str(index) + '/' + str(index_len))
            tds.set_group_index(group, index)
            datenum = tds.metagrp.groups[group].variables['IntegrationMidPointTime'][index]
            # 0.5 deg error approx 55 km error
            if (angle_between(tds.r_r, tds.r_t) <= search_error):
                string = 'G: ' + group + ' I: ' + str(index) + ' - ' + \
                        str(datenum) + ' - ' + str(datenum_to_pytime(float(datenum))) + ' - Lat: ' + \
                        str(tds.lat_sp_tds) + ' Lon: ' + str(tds.lon_sp_tds) + \
                        'Angle (deg): ' + str(angle_between(tds.r_r, tds.r_t)*180/np.pi) + '\n'
                print(string)
                search_output = search_output + string
    return search_output

if __name__ == '__main__':
    main()

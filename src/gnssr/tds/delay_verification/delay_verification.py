#!/usr/bin/env python

"""

Verified Delay Scale. Same results as in following publication: 

[1]P. Addabbo, G. Giangregorio, C. Galdi, and M. di Bisceglie, “Simulation of 
TechDemoSat-1 Delay-Doppler Maps for GPS Ocean Reflectometry,” IEEE Journal of 
Selected Topics in Applied Earth Observations and Remote Sensing, vol. 10, no. 
9, pp. 4256–4268, Sep. 2017.  

"""

from gnssr.tds.tds_data import *
from gnssr.utils import *

def main(): 
    # Verification 1
    #file_root_name = 'raw/L1B/2015-02-11-H18'
    #group = '000139'
    #index = 839

    # Verification 2
    file_root_name = 'raw/L1B/2015-02-11-H18'
    group = '000139'
    index = 592

    tds = tds_data(file_root_name)
    tds.set_group_index(group, index)
    tds.plot_ddm()

    datenum = tds.metagrp.groups[tds.group].variables['IntegrationMidPointTime'][tds.index]
    lat = tds.metagrp.groups[tds.group].variables['SpecularPointLat'][tds.index]
    lon = tds.metagrp.groups[tds.group].variables['SpecularPointLon'][tds.index]

    string = 'G: ' + group + ' I: ' + str(index) + ' - ' + \
            str(datenum_to_pytime(float(datenum))) + \
            ' - Lat: ' + str(lat) + ' Lon: ' + str(lon) + '\n'
    print(string)

    plt.show(block=False)
    plt.show()

if __name__ == '__main__':
    main()

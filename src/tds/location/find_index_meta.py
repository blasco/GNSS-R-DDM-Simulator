from netCDF4 import Dataset
import numpy as np

def find_index_meta(rootgrp, metagrp, group, index):
    datenum = rootgrp.groups[group].variables['IntegrationMidPointTime'][index]
    datenums_meta = np.array(metagrp.variables['IntegrationMidPointTime'])
    for i, datenum_meta in enumerate(datenums_meta):
        if datenum_meta == datenum:
            #print("Foud: {0} | {1}".format(i, datenum_meta))
            return i
            break

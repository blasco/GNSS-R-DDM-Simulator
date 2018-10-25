#!/usr/bin/env python

from netCDF4 import Dataset
from tqdm import tqdm
import concurrent.futures
import os
import numpy as np

rootgrp = Dataset(os.path.join(os.environ['TDS_ROOT'],"raw/2015-04-01-H00-ddm.nc"), "r", format="NETCDF4")
metagrp = Dataset(os.path.join(os.environ['TDS_ROOT'],"raw/2015-04-01-H00-metadata.nc"), "r", format="NETCDF4")
group = '000095'
index = 525;

datenum = rootgrp.groups[group].variables['IntegrationMidPointTime'][index]
print(datenum)

datenums_meta = np.array(metagrp.variables['IntegrationMidPointTime'])
print(datenums_meta)

for i, datenum_meta in enumerate(datenums_meta):
    if datenum_meta == datenum:
        print("Foud: {0} | {1}".format(i, datenum_meta))
        break

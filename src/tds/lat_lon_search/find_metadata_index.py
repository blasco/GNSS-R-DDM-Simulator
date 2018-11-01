#!/usr/bin/env python

from netCDF4 import Dataset
from tqdm import tqdm
import concurrent.futures
import os
import numpy as np

rootgrp = Dataset(root_path, "r", format="NETCDF4")
metagrp = Dataset(meta_path, "r", format="NETCDF4")

def find_metadata_index(group, index):
    datenum = rootgrp.groups[group].variables['IntegrationMidPointTime'][index]
    #print(datenum)

    datenums_meta = np.array(metagrp.variables['IntegrationMidPointTime'])
    #print(datenums_meta)

    for i, datenum_meta in enumerate(datenums_meta):
        if datenum_meta == datenum:
            #print("Foud: {0} | {1}".format(i, datenum_meta))
            return i
            break

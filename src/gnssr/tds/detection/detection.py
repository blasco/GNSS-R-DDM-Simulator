#!/usr/bin/env python

from netCDF4 import Dataset
import os
import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
import matplotlib.patches as mpatches
from datetime import *

# Petronius Oil Platform
#file_dir = os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2017-11-26-H06-')
#rootgrp = Dataset(file_dir+"DDMs.nc", "r", format="NETCDF4")
#metagrp = Dataset(file_dir+"metadata.nc", "r", format="NETCDF4")
#group = '000047'
#start_0 = 590

# Di Simone Oil Platform
#rootgrp = Dataset("../raw/2015-04-01-H00-ddm.nc", "r", format="NETCDF4")
#metagrp = Dataset("../raw/2015-04-01-H00-metadata.nc", "r", format="NETCDF4")
#group = '000095'
#start_0 = 525

file_root_name = 'raw/L1B/2016-09-03-H00'
rootgrp = Dataset(os.path.join(os.environ['TDS_ROOT'], file_root_name+'-DDMs.nc') , "r", format="NETCDF4")
metagrp = Dataset(os.path.join(os.environ['TDS_ROOT'], file_root_name+'-metadata.nc') , "r", format="NETCDF4")
group = '000006'
index_start = 295
index_end = index_start + 10

# Di Simone Oil Platform
#filename = os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2017-11-11-H18-DDMs.nc')
#rootgrp = Dataset(filename, "r", format="NETCDF4")
#group = '000028'
#start_0 = 670

# Di Simone Oil Platform 
#File: raw/L1B_Catalogue/2017-11/01/H06/2017-11.01.H06.kmz
#file_dir = os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2017-11-07-H06-')
#rootgrp = Dataset(file_dir+'DDMs.nc', 'r', format='NETCDF4')
#metagrp = Dataset(file_dir+'metadata.nc', 'r', format='NETCDF4')
#group = '000021'
#start_0 = 426 + 8

# Atlantis PQ Oil Platform
#file_dir = os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2017-11-08-H18-')
#rootgrp = Dataset(file_dir+"DDMs.nc", "r", format="NETCDF4")
#metagrp = Dataset(file_dir+"metadata.nc", "r", format="NETCDF4")
#group = '000059'
#start_0 = 580;

# Songa Mercur Oil Platform
#file_dir = os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2017-11-26-H18-')
#rootgrp = Dataset(file_dir+"DDMs.nc", "r", format="NETCDF4")
#metagrp = Dataset(file_dir+"metadata.nc", "r", format="NETCDF4")
#group = '000041'
#start_0 = 625;

# Devil's Tower Oil Platform
#search_lat = 28.19013
#search_lon = -88.49552
#file_dir = os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2018-07-30-H18-')
#rootgrp = Dataset(file_dir+"DDMs.nc", "r", format="NETCDF4")
#metagrp = Dataset(file_dir+"metadata.nc", "r", format="NETCDF4")
#group = '000076'
#start_0 = 75 - 30

min_col = 60
max_col = 85

def datenum_to_pytime(matlab_datenum):
    python_datetime = datetime.fromordinal(int(matlab_datenum)) + timedelta(days=matlab_datenum%1) - timedelta(days = 366)
    return python_datetime

def normalize(mat):
    return (mat - np.min(mat))/(np.max(mat)-np.min(mat))

def cut_noise_region(ddm,ddm_ref,new_value=0):
    '''
    For each pixel in ddm, replaces the pixel with new_value if the pixel 
    with same indeces in ddm_ref is less than 1.4*cut_threshold.
    cut_threshold is computed as the mean of the noise region of the ddm_ref 
    (from column 0 to 40)
    '''
    noise_mean = np.min(ddm_ref)
    for row_i, row in enumerate(ddm_ref):
        for col_i, val in enumerate(row):
            if(col_i >= 0 and col_i <= 40):
                noise_mean = 0.5*noise_mean + 0.5*ddm_ref[row_i][col_i]

    ddm_cut = np.zeros(ddm.shape) + new_value
    cut_threshold = noise_mean*3
    for row_i, row in enumerate(ddm_cut):
        for col_i, val in enumerate(row):
            if(col_i >= min_col and col_i <= max_col):
                if(ddm_ref[row_i][col_i]>=cut_threshold):
                    ddm_cut[row_i][col_i] = ddm[row_i][col_i]
    return ddm_cut

def next_max(ddm,col):
    next_max = 0
    for row_i, row in enumerate(ddm):
        for col_i, val in enumerate(row):
            if col_i >= col:
                if ddm[row_i][col_i] > next_max:
                    next_max = ddm[row_i][col_i] 
    return next_max

for start in range(index_start,index_end):
    # original
    ddm_original = np.array(rootgrp.groups[group].variables['DDM'][start].data)

    # sum low pass filter
    #n = 20
    #offset = 0
    #tau=0.08
    #ddm_sum = np.array(rootgrp.groups[group].variables['DDM'][start-offset].data)
    #for i in range(n):
    #    ddm_i = np.array(rootgrp.groups[group].variables['DDM'][start-offset-i].data)
    #    for row_i, row in enumerate(ddm_sum):
    #        for col_i, val in enumerate(row):
    #            val_i = ddm_i[row_i][col_i]
    #            val = ddm_sum[row_i][col_i]
    #            ddm_sum[row_i][col_i] = val + tau*(val_i-val)

    # sum min filter
    n = 0
    offset = 0
    ddm_sum_0 = np.array(rootgrp.groups[group].variables['DDM'][start-offset].data)
    for i in range(n):
        ddm_i = np.array(rootgrp.groups[group].variables['DDM'][start-offset-i].data)
        for row_i, row in enumerate(ddm_sum_0):
            for col_i, val in enumerate(row):
                val_i = ddm_i[row_i][col_i]
                val = ddm_sum_0[row_i][col_i]
                min,max = np.sort([val_i,val])
                ddm_sum_0[row_i][col_i] = min

    # sum min low pass filter
    tau=0.2
    ddm_sum = ddm_sum_0
    '''
    for i in range(n):
        ddm_i = np.array(rootgrp.groups[group].variables['DDM'][start-i].data)
        for row_i, row in enumerate(ddm_sum):
            for col_i, val in enumerate(row):
                val_i = ddm_i[row_i][col_i]
                val = ddm_sum[row_i][col_i]
                min,max = np.sort([val_i,val])
                ddm_sum[row_i][col_i] = val + tau*(min - val)
                '''
    ddm_sum_cut = cut_noise_region(ddm_sum, ddm_sum)
    ddm_original_cut = cut_noise_region(ddm_original, ddm_sum)

    # substracted
    ddm_sub = ddm_original_cut - 0.85*ddm_sum_cut
    if (np.min(ddm_sub) < 0):
        ddm_sub = ddm_sub - np.min(ddm_sub)
    cut_value_fig = np.min(ddm_sub)
    ddm_sub_fig = cut_noise_region(ddm_sub,ddm_sum,cut_value_fig)
    ddm_sub = cut_noise_region(ddm_sub,ddm_sum)

    # detections
    ddm_detections = ddm_sub*0
    threshold = np.max(ddm_sum)*0.4
    print(np.max(ddm_sub))
    print(threshold)
    print(start)
    for row_i, row in enumerate(ddm_sub):
        for col_i, val in enumerate(row):
            if(col_i >= min_col and col_i <= max_col):
                if(ddm_sub[row_i][col_i] >= threshold):
                    ddm_detections[row_i][col_i] = 1

    fig_original = plt.figure(figsize=(10, 4))
    im_original = plt.imshow(ddm_original, cmap='viridis', extent=(0,127,-10,9))
    plt.show(block=False)

    #fig_sum = plt.figure(figsize=(10, 4))
    #im_sum = plt.imshow(ddm_sum, cmap='viridis', extent=(0,127,-10,9))    
    #plt.show(block=False)

    #fig_original_cut = plt.figure(figsize=(10, 4))
    #im_original_cut = plt.imshow(ddm_original_cut, cmap='viridis', extent=(0,127,-10,9))
    #plt.show(block=False)

    #fig_sum_cut = plt.figure(figsize=(10, 4))
    #im_sum_cut = plt.imshow(ddm_sum_cut, cmap='viridis', extent=(0,127,-10,9))    
    #plt.show(block=False)

    #fig_sub = plt.figure(figsize=(10, 4))
    #im_sub = plt.imshow(ddm_sub, cmap='viridis', extent=(0,127,-10,9))
    #plt.show(block=False)

    #fig_detections = plt.figure(figsize=(10, 4))
    #im_detections = plt.imshow(ddm_detections, cmap='viridis', extent=(0,127,-10,9))
    #plt.show(block=False)

    #all_labels = label(ddm_detections)
    #fig_labels,ax_labels = plt.subplots(1,figsize=(10, 4))
    #ax_labels.imshow(ddm_original, cmap='viridis')
    #plt.show(block=False)

    #for region in regionprops(all_labels):
    #    minr, minc, maxr, maxc = region.bbox
    #    l = 1 
    #    rect = mpatches.Rectangle((minc-l, minr-l), maxc - minc + 2*l-1, maxr - minr + 2*l - 1, fill=False, edgecolor='red', linewidth=2)
    #    ax_labels.add_patch(rect)

    datenum = metagrp.groups[group].variables['IntegrationMidPointTime'][start]
    lat = metagrp.groups[group].variables['SpecularPointLat'][start]
    lon = metagrp.groups[group].variables['SpecularPointLon'][start]
    string = 'G: ' + group + ' I: ' + str(start) + ' - ' + \
            str(datenum) + ' - ' + str(datenum_to_pytime(float(datenum))) + ' - Lat: ' + \
            str(lat) + ' Lon: ' + str(lon) + '\n'
    t = plt.text(5, 5, string, {'color': 'w', 'fontsize': 12})

plt.show()

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
rootgrp = Dataset(os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2015-04-01-H00-DDMs.nc') , "r", format="NETCDF4")
metagrp = Dataset(os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2015-04-01-H00-metadata.nc') , "r", format="NETCDF4")
group = '000095'
index = 525

# The region of interest lies between delay column 60 and delay column 80
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
    with same indeces in ddm_ref is less than 0.3
    cut_threshold is computed as the mean of the noise region of the ddm_ref 
    (from column 0 to 40)
    '''
    ddm_cut = np.zeros(ddm.shape) + new_value
    cut_threshold = 0.15
    for row_i, row in enumerate(ddm_cut):
        for col_i, val in enumerate(row):
            if(col_i >= min_col and col_i <= max_col):
                if(ddm_ref[row_i][col_i]>=cut_threshold):
                    ddm_cut[row_i][col_i] = ddm[row_i][col_i]
    return ddm_cut


# 1. Original
ddm_original = normalize(np.array(rootgrp.groups[group].variables['DDM'][index].data))

# 2. Sea clutter estimation. 

# As target appear as bright spots, the initial estimation is based of the 
# composition of the minimum values for each pixel for the last n measurements
n = 10
offset = 0
sea_clutter_0 = normalize(np.array(rootgrp.groups[group].variables['DDM'][index-offset].data))
for i in range(n):
    ddm_i = normalize(np.array(rootgrp.groups[group].variables['DDM'][index-offset-i].data))
    sea_clutter_0 = normalize(sea_clutter_0)
    for row_i, row in enumerate(sea_clutter_0):
        for col_i, val in enumerate(row):
            val_i = ddm_i[row_i][col_i]
            val = sea_clutter_0[row_i][col_i]
            min,max = np.sort([val_i,val])
            sea_clutter_0[row_i][col_i] = 0.5*min +0.5*max
sea_clutter_0 = normalize(sea_clutter_0)

# Using the sea_clutter_0 as initial estimation a low pass filter that gives 
# more weight to the lower values is applied 
n = 200
tau = 0.08
sea_clutter = np.array(sea_clutter_0)
for i in range(n):
    ddm_i = normalize(np.array(rootgrp.groups[group].variables['DDM'][index-offset-i].data))
    sea_clutter = normalize(sea_clutter)
    for row_i, row in enumerate(sea_clutter):
        for col_i, val in enumerate(row):
            val_i = ddm_i[row_i][col_i]
            val = sea_clutter[row_i][col_i]
            #min,max = np.sort([val_i,val])
            #new_val = 0.7*min + 0.3*max;
            sea_clutter[row_i][col_i] = val + tau*(val_i - val)
sea_clutter = normalize(sea_clutter)

# Only the region of the wake is relevant for detection, so we cut the 
# irrelevant region
sea_clutter_cut = cut_noise_region(sea_clutter, sea_clutter)
ddm_original_cut = cut_noise_region(ddm_original, sea_clutter)

# 3. Sea clutter substracted Difference Map
ddm_diff = ddm_original_cut - 0.85*sea_clutter_cut
if (np.min(ddm_diff) < 0):
    ddm_diff = ddm_diff - np.min(ddm_diff)
cut_value_fig = np.min(ddm_diff)
ddm_diff_fig = cut_noise_region(ddm_diff,sea_clutter,cut_value_fig)
ddm_diff = cut_noise_region(ddm_diff,sea_clutter)

# 4. Over threshold detection
ddm_detections = np.zeros(ddm_diff.shape)
threshold = 0.38
print(np.max(ddm_diff))
print(threshold)
print(index)
for row_i, row in enumerate(ddm_diff):
    for col_i, val in enumerate(row):
        if(col_i >= min_col and col_i <= max_col):
            if(ddm_diff[row_i][col_i] >= threshold):
                ddm_detections[row_i][col_i] = 1

# Plotting
fig_original = plt.figure(figsize=(10, 4))
im_original = plt.imshow(ddm_original, cmap='viridis', extent=(0,127,-10,9))
plt.show(block=False)

fig_sea_clutter = plt.figure(figsize=(10, 4))
im_sea_clutter = plt.imshow(sea_clutter, cmap='viridis', extent=(0,127,-10,9))    
plt.show(block=False)

fig_sub = plt.figure(figsize=(10, 4))
im_sub = plt.imshow(ddm_diff, cmap='viridis', extent=(0,127,-10,9))
plt.show(block=False)

fig_detections = plt.figure(figsize=(10, 4))
im_detections = plt.imshow(ddm_detections, cmap='viridis', extent=(0,127,-10,9))
plt.show(block=False)

all_labels = label(ddm_detections)
fig_labels,ax_labels = plt.subplots(1,figsize=(10, 4))
ax_labels.imshow(ddm_original, cmap='viridis')
plt.show(block=False)

for region in regionprops(all_labels):
    minr, minc, maxr, maxc = region.bbox
    l = 1 
    rect = mpatches.Rectangle((minc-l, minr-l), maxc - minc + 2*l-1, maxr - minr + 2*l - 1, fill=False, edgecolor='red', linewidth=2)
    ax_labels.add_patch(rect)

datenum = metagrp.groups[group].variables['IntegrationMidPointTime'][index]
lat = metagrp.groups[group].variables['SpecularPointLat'][index]
lon = metagrp.groups[group].variables['SpecularPointLon'][index]
string = 'G: ' + group + ' I: ' + str(index) + ' - ' + \
        str(datenum) + ' - ' + str(datenum_to_pytime(float(datenum))) + ' - Lat: ' + \
        str(lat) + ' Lon: ' + str(lon) + '\n'
t = plt.text(5, 5, string, {'color': 'w', 'fontsize': 12})

plt.show()

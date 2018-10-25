#!/usr/bin/env python

from netCDF4 import Dataset
import os
import numpy
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
import matplotlib.patches as mpatches

#filename = os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2017-11-11-H18-DDMs.nc')
#rootgrp = Dataset(filename, "r", format="NETCDF4")
#group = '000028'
#start = 662

# Di Simone Oil Platform 
#File: raw/L1B_Catalogue/2017-11/01/H06/2017-11.01.H06.kmz
file_dir = os.path.join(os.environ['TDS_ROOT'], 'raw/L1B/2017-11-07-H06-')
rootgrp = Dataset(file_dir+'DDMs.nc', 'r', format='NETCDF4')
metagrp = Dataset(file_dir+'metadata.nc', 'r', format='NETCDF4')
group = '000021'
start = 426 + 11

#rootgrp = Dataset("../raw/2015-04-01-H00-ddm.nc", "r", format="NETCDF4")
#group = '000095'
#start = 525

min_col = 60
max_col = 80

def normalize(mat):
    return (mat - numpy.min(mat))/(numpy.max(mat)-numpy.min(mat))

def cut_noise_region(ddm,ddm_ref):
    noise_ref = ddm_ref[0][0]
    for row_i, row in enumerate(ddm_ref):
        for col_i, val in enumerate(row):
            if(col_i >= 0 and col_i <= 30):
                noise_ref = 0.5*noise_ref + 0.5*ddm_ref[row_i][col_i]

    ddm_out = ddm*0
    cut_threshold = noise_ref + 0.15
    for row_i, row in enumerate(ddm_out):
        for col_i, val in enumerate(row):
            if(col_i >= min_col and col_i <= max_col):
                if(ddm_ref[row_i][col_i]>=cut_threshold):
                    ddm_out[row_i][col_i] = ddm[row_i][col_i]
    return ddm_out

# normalized wiht the cut!
# original
ddm_original = normalize(numpy.array(rootgrp.groups[group].variables['DDM'][start].data))

# sum min filter
n = 30
offset = 0
ddm_sum_0 = normalize(numpy.array(rootgrp.groups[group].variables['DDM'][start-offset].data))
for i in range(n):
    ddm_i = normalize(numpy.array(rootgrp.groups[group].variables['DDM'][start-offset-i].data))
    for row_i, row in enumerate(ddm_sum_0):
        for col_i, val in enumerate(row):
            val_i = ddm_i[row_i][col_i]
            val = ddm_sum_0[row_i][col_i]
            min,max = numpy.sort([val_i,val])
            ddm_sum_0[row_i][col_i] = 0.95*min + 0.05*max

# sum min low pass filter
tau=0.2
ddm_sum = ddm_sum_0
for i in range(n):
    ddm_i = normalize(numpy.array(rootgrp.groups[group].variables['DDM'][start-i].data))
    for row_i, row in enumerate(ddm_sum):
        for col_i, val in enumerate(row):
            val_i = ddm_i[row_i][col_i]
            val = ddm_sum[row_i][col_i]
            min,max = numpy.sort([val_i,val])
            ddm_sum[row_i][col_i] = val + tau*(0.95*min + 0.05*max - val)

ddm_sum_cut = normalize(cut_noise_region(ddm_sum, ddm_sum))
ddm_original_cut = normalize(cut_noise_region(ddm_original, ddm_sum))

# substracted
ddm_sub = ddm_original_cut - ddm_sum_cut
for row_i, row in enumerate(ddm_sub):
    for col_i, val in enumerate(row):
        if ddm_sub[row_i][col_i] < 0:
            ddm_sub[row_i][col_i] = 0

# cut
ddm_cut = cut_noise_region(ddm_sub,ddm_sum)

# detections
ddm_detections = ddm_cut*1
print(numpy.max(ddm_cut))
threshold = 0.22
for row_i, row in enumerate(ddm_cut):
    for col_i, val in enumerate(row):
        if(col_i >= min_col and col_i <= max_col):
            if(ddm_detections[row_i][col_i] >= threshold):
                ddm_detections[row_i][col_i] = 1
            else:
                ddm_detections[row_i][col_i] = 0
        else:
            ddm_detections[row_i][col_i] = 0

fig_original = plt.figure(figsize=(10, 4))
im_original = plt.imshow(ddm_original, cmap='viridis', extent=(0,127,-10,9))
plt.show(block=False)

fig_sum = plt.figure(figsize=(10, 4))
im_sum = plt.imshow(ddm_sum, cmap='viridis', extent=(0,127,-10,9))    
plt.show(block=False)

fig_original_cut = plt.figure(figsize=(10, 4))
im_original_cut = plt.imshow(ddm_original_cut, cmap='viridis', extent=(0,127,-10,9))
plt.show(block=False)

fig_sum_cut = plt.figure(figsize=(10, 4))
im_sum_cut = plt.imshow(ddm_sum_cut, cmap='viridis', extent=(0,127,-10,9))    
plt.show(block=False)

fig_sub = plt.figure(figsize=(10, 4))
im_sub = plt.imshow(ddm_sub, cmap='viridis', extent=(0,127,-10,9))
plt.show(block=False)

fig_cut = plt.figure(figsize=(10, 4))
im_cut = plt.imshow(ddm_cut, cmap='viridis', extent=(0,127,-10,9))
plt.show(block=False)

fig_detections = plt.figure(figsize=(10, 4))
im_detections = plt.imshow(ddm_detections, cmap='viridis', extent=(0,127,-10,9))
plt.show(block=False)

all_labels = label(ddm_detections)
fig_labels,ax_labels = plt.subplots(1,figsize=(10, 4))
ax_labels.imshow(ddm_original, cmap='viridis')

plt.show(block=False)

for region in regionprops(all_labels):
    # draw rectangle around segmented coins
    minr, minc, maxr, maxc = region.bbox
    l = 1 
    rect = mpatches.Rectangle((minc-l, minr-l), maxc - minc + 2*l-1, maxr - minr + 2*l - 1, fill=False, edgecolor='red', linewidth=2)
    ax_labels.add_patch(rect)

plt.show()

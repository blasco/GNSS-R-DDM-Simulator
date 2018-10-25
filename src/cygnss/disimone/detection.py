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
#start_0 = 670

rootgrp = Dataset("../raw/2015-04-01-H00-ddm.nc", "r", format="NETCDF4")
group = '000095'
start_0 = 524

min_col = 60
max_col = 80

def normalize(mat):
    return mat.astype(float)

def cut_noise_region(ddm,ddm_ref):
    ddm_out = ddm*0
    cut_threshold = noise_ref*1.075
    for row_i, row in enumerate(ddm_out):
        for col_i, val in enumerate(row):
            if(col_i >= min_col and col_i <= max_col):
                if(ddm_ref[row_i][col_i]>=cut_threshold):
                    ddm_out[row_i][col_i] = ddm[row_i][col_i]
    return ddm_out

def next_max(ddm,col):
    next_max = 0
    for row_i, row in enumerate(ddm):
        for col_i, val in enumerate(row):
            if col_i >= col:
                if ddm[row_i][col_i] > next_max:
                    next_max = ddm[row_i][col_i] 
    return next_max


for start in range(start_0,start_0+10):
    # normalized wiht the cut!
    # original
    ddm_original = normalize(numpy.array(rootgrp.groups[group].variables['DDM'][start].data))

    # sum
    n = 200
    offset = 0
    tau=0.05
    ddm_sum = normalize(numpy.array(rootgrp.groups[group].variables['DDM'][start-offset].data))
    for i in range(n):
        ddm_i = normalize(numpy.array(rootgrp.groups[group].variables['DDM'][start-offset-i].data))
        for row_i, row in enumerate(ddm_sum):
            for col_i, val in enumerate(row):
                val_i = ddm_i[row_i][col_i]
                val = ddm_sum[row_i][col_i]
                ddm_sum[row_i][col_i] = val + tau*(val_i-val)

    noise_ref = numpy.min(ddm_sum)
    for row_i, row in enumerate(ddm_sum):
        for col_i, val in enumerate(row):
            if(col_i >= 0 and col_i <= 40):
                noise_ref = 0.5*noise_ref + 0.5*ddm_sum[row_i][col_i]

    ddm_sum_cut = (cut_noise_region(ddm_sum, ddm_sum))
    ddm_original_cut = normalize(cut_noise_region(ddm_original, ddm_sum))

    threshold_ref = noise_ref
    for row_i, row in enumerate(ddm_original):
        for col_i, val in enumerate(row):
            if(col_i >= min_col and col_i <= max_col):
                threshold_ref = 0.5*threshold_ref + 0.5*ddm_sum[row_i][col_i]

    # substracted
    ddm_sub = ddm_original_cut - ddm_sum_cut
    if (numpy.min(ddm_sub) < 0):
        ddm_sub = ddm_sub - numpy.min(ddm_sub)
    ddm_sub = cut_noise_region(ddm_sub,ddm_sum)

    # detections
    beta = 1.35 
    ddm_detections = ddm_sub*0
    prev_col = 0
    threshold = 0
    print(start)
    for row_i, row in enumerate(ddm_sub):
        for col_i, val in enumerate(row):
            if(col_i >= min_col and col_i <= max_col):
                threshold = 0.15*threshold_ref*(next_max(ddm_sum,col_i)/threshold_ref)**1.8
                if(ddm_sub[row_i][col_i] >= threshold):
                    ddm_detections[row_i][col_i] = 1

    #fig_original = plt.figure(figsize=(10, 4))
    #im_original = plt.imshow(ddm_original, cmap='viridis', extent=(0,127,-10,9))
    #plt.show(block=False)

    #fig_sum = plt.figure(figsize=(10, 4))
    #im_sum = plt.imshow(ddm_sum, cmap='viridis', extent=(0,127,-10,9))    
    #plt.show(block=False)

    #fig_original_cut = plt.figure(figsize=(10, 4))
    #im_original_cut = plt.imshow(ddm_original_cut, cmap='viridis', extent=(0,127,-10,9))
    #plt.show(block=False)

    #fig_sum_cut = plt.figure(figsize=(10, 4))
    #im_sum_cut = plt.imshow(ddm_sum_cut, cmap='viridis', extent=(0,127,-10,9))    
    #plt.show(block=False)

    fig_sub = plt.figure(figsize=(10, 4))
    im_sub = plt.imshow(ddm_sub, cmap='viridis', extent=(0,127,-10,9))
    plt.show(block=False)

    #fig_detections = plt.figure(figsize=(10, 4))
    #im_detections = plt.imshow(ddm_detections, cmap='viridis', extent=(0,127,-10,9))
    #plt.show(block=False)

    all_labels = label(ddm_detections)
    fig_labels,ax_labels = plt.subplots(1,figsize=(10, 4))
    ax_labels.imshow(ddm_original, cmap='viridis')
    plt.show(block=False)

    for region in regionprops(all_labels):
        minr, minc, maxr, maxc = region.bbox
        l = 1 
        rect = mpatches.Rectangle((minc-l, minr-l), maxc - minc + 2*l-1, maxr - minr + 2*l - 1, fill=False, edgecolor='red', linewidth=2)
        ax_labels.add_patch(rect)

plt.show()

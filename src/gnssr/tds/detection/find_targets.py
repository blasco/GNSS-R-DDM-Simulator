#!/usr/bin/env python
from netCDF4 import Dataset
import os
import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
import matplotlib.patches as mpatches
from datetime import *

from gnssr.tds.tds_data import *
from gnssr.utils import *
from gnssr.targets import *
from gnssr.tds.search_database.cdf4_search import *

def datenum_to_pytime(matlab_datenum):
    python_datetime = datetime.fromordinal(int(matlab_datenum)) + timedelta(days=matlab_datenum%1) - timedelta(days = 366)
    return python_datetime

class target_processor:

    def __init__(self):
        self.min_col = 60
        #self.max_col = 85
        self.max_col = 95
        self.ddm_list = []

    def cut_noise_region(self, ddm,ddm_ref,new_value=0):
        """
        For each pixel in ddm, replaces the pixel with new_value if the pixel 
        with same indeces in ddm_ref is less than 0.3
        cut_threshold is computed as the mean of the noise region of the ddm_ref 
        (from column 0 to 40)
        """
        ddm_cut = np.zeros(ddm.shape) + new_value
        cut_threshold = 0.15
        for row_i, row in enumerate(ddm_cut):
            for col_i, val in enumerate(row):
                if(col_i >= self.min_col and col_i <= self.max_col):
                    if(ddm_ref[row_i][col_i]>=cut_threshold):
                        ddm_cut[row_i][col_i] = ddm[row_i][col_i]
        return ddm_cut

    def process_ddm(self, ddm_raw):
        self.ddm_original = normalize(ddm_raw) 

        n = 30 
        if len(self.ddm_list) < n :
            self.ddm_list.insert(0, self.ddm_original)
            return

        # 1. Sea clutter estimation. 
        # As target appear as bright spots, the initial estimation is based of the 
        # composition of the minimum values for each pixel for the last n measurements
        n = 30 
        self.sea_clutter_0 = self.ddm_list[0]
        for i in range(1, n):
            ddm_i = normalize(self.ddm_list[i])
            self.sea_clutter_0 = normalize(self.sea_clutter_0)
            for row_i, row in enumerate(self.sea_clutter_0):
                for col_i, val in enumerate(row):
                    val_i = ddm_i[row_i][col_i]
                    val = self.sea_clutter_0[row_i][col_i]
                    min,max = np.sort([val_i,val])
                    self.sea_clutter_0[row_i][col_i] = 0.5*min +0.5*max
        self.sea_clutter_0 = normalize(self.sea_clutter_0)

        # Using the self.sea_clutter_0 as initial estimation a low pass filter that gives 
        # more weight to the lower values is applied 
        n = 30 
        tau = 0.08
        self.sea_clutter = np.array(self.sea_clutter_0)
        for i in range(1, n):
            ddm_i = normalize(self.ddm_list[i])
            self.sea_clutter = normalize(self.sea_clutter)
            for row_i, row in enumerate(self.sea_clutter):
                for col_i, val in enumerate(row):
                    val_i = ddm_i[row_i][col_i]
                    val = self.sea_clutter[row_i][col_i]
                    #min,max = np.sort([val_i,val])
                    #new_val = 0.7*min + 0.3*max;
                    self.sea_clutter[row_i][col_i] = val + tau*(val_i - val)
        self.sea_clutter = normalize(self.sea_clutter)

        # Only the region of the wake is relevant for detection, so we cut the 
        # irrelevant region
        self.sea_clutter_cut = self.cut_noise_region(self.sea_clutter, self.sea_clutter)
        self.ddm_original_cut = self.cut_noise_region(self.ddm_original, self.sea_clutter)

        # 2. Sea clutter substracted Difference Map
        ddm_diff = self.ddm_original_cut - 0.85*self.sea_clutter_cut
        if (np.min(ddm_diff) < 0):
            ddm_diff = ddm_diff - np.min(ddm_diff)
        cut_value_fig = np.min(ddm_diff)
        ddm_diff_fig = self.cut_noise_region(ddm_diff,self.sea_clutter,cut_value_fig)
        ddm_diff = self.cut_noise_region(ddm_diff,self.sea_clutter)

        # 3. Over threshold detection
        ddm_detections = np.zeros(ddm_diff.shape)
        #threshold = 0.38
        threshold = 0.5
        print(np.max(ddm_diff))
        print(threshold)
        for row_i, row in enumerate(ddm_diff):
            for col_i, val in enumerate(row):
                if(col_i >= self.min_col and col_i <= self.max_col):
                    if(ddm_diff[row_i][col_i] >= threshold):
                        ddm_detections[row_i][col_i] = 1

        # Plotting
        """
        fig_original = plt.figure(figsize=(10, 4))
        im_original = plt.imshow(self.ddm_original, cmap='viridis', extent=(0,127,-10,9))
        plt.show(block=False)

        fig_self.sea_clutter = plt.figure(figsize=(10, 4))
        im_self.sea_clutter = plt.imshow(self.sea_clutter, cmap='viridis', extent=(0,127,-10,9))    
        plt.show(block=False)

        fig_sub = plt.figure(figsize=(10, 4))
        im_sub = plt.imshow(ddm_diff, cmap='viridis', extent=(0,127,-10,9))
        plt.show(block=False)

        fig_detections = plt.figure(figsize=(10, 4))
        im_detections = plt.imshow(ddm_detections, cmap='viridis', extent=(0,127,-10,9))
        plt.show(block=False)
        """

        """
        datenum = metagrp.groups[group].variables['IntegrationMidPointTime'][index]
        lat = metagrp.groups[group].variables['SpecularPointLat'][index]
        lon = metagrp.groups[group].variables['SpecularPointLon'][index]
        string = 'G: ' + group + ' I: ' + str(index) + ' - ' + \
                str(datenum) + ' - ' + str(datenum_to_pytime(float(datenum))) + ' - Lat: ' + \
                str(lat) + ' Lon: ' + str(lon) + '\n'
        t = plt.text(5, 5, string, {'color': 'w', 'fontsize': 12})
        """

        """
        # Plot with pause
        plt.draw()
        plt.pause(1) 
        input("<Hit Enter To Continue>")
        plt.close(fig_original)
        plt.close(fig_self.sea_clutter)
        plt.close(fig_sub)
        plt.close(fig_detections)
        plt.close(fig_labels)
        """

    def plot_targets():
        self.fig_labels, self.ax_labels = plt.subplots(1,figsize=(10, 4))
        [p.remove() for p in reversed(self.ax_labels.patches)] # Clear previous patches
        selfall_labels = label(ddm_detections)
        self.ax_labels.imshow(self.ddm_original, cmap='viridis')
        target_flag = False
        for region in regionprops(selfall_labels):
            target_flag = True
            minr, minc, maxr, maxc = region.bbox
            l = 1 
            rect = mpatches.Rectangle((minc-l, minr-l), maxc - minc + 2*l-1, maxr - minr + 2*l - 1, fill=False, edgecolor='red', linewidth=2)
            self.ax_labels.add_patch(rect)

        plt.draw()
        plt.pause(0.0001)


def main():
    # Di Simone Oil Platform Data
    file_root_name = 'raw/L1B/2015-04-01-H00'
    target = targets['hibernia']
    group = '000095'
    index = 525
    tds = tds_data(file_root_name)
    p = target_processor();
    for i in range(index - 200, index + 10):
        ddm = tds.rootgrp.groups[group].variables['DDM'][i].data
        p.process_ddm(ddm)

if __name__ == '__main__':
    main()

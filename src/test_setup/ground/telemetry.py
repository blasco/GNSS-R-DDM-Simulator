#!/usr/bin/env python

from netCDF4 import Dataset
import os
import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
import matplotlib.patches as mpatches
from datetime import *
import pickle

from gnssr.tds.tds_data import *
from gnssr.utils import *
from gnssr.targets import *
from gnssr.tds.search_target.cdf4_search import *

def datenum_to_pytime(matlab_datenum):
    python_datetime = datetime.fromordinal(int(matlab_datenum)) + timedelta(days=matlab_datenum%1) - timedelta(days = 366)
    return python_datetime

def normalize(mat):
    return (mat - np.min(mat))/(np.max(mat)-np.min(mat))

class processor:

    def __init__(self):
        self.fig_labels, self.ax_labels = plt.subplots(1,figsize=(10, 4))

    def load_telemetry(self, telemetry):
        self.metadata = telemetry['metadata']
        self.ddm_original = telemetry['ddm_original']
        self.all_labels = telemetry['all_labels']

    def plot_telemetry(self):
        [p.remove() for p in reversed(self.ax_labels.patches)] # Clear previous patches
        self.ax_labels.imshow(self.ddm_original, cmap='viridis')
        target_flag = False
        for region in regionprops(self.all_labels):
            target_flag = True
            minr, minc, maxr, maxc = region.bbox
            l = 1 
            rect = mpatches.Rectangle((minc-l, minr-l), maxc - minc + 2*l-1, maxr - minr + 2*l - 1, fill=False, edgecolor='red', linewidth=2)
            self.ax_labels.add_patch(rect)
        if target_flag:
            print('TM: ' + self.metadata)
        plt.draw()
        plt.pause(0.0001)

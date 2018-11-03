#!/usr/bin/env python

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
from datetime import *
from dateutil import tz

earth_a = 6378137 # meters
earth_b = 6356752.314245 # meters

def datenum_to_pytime(matlab_datenum):
    python_datetime = datetime.fromordinal(int(matlab_datenum)) + timedelta(days=matlab_datenum%1) - timedelta(days=366)
    return python_datetime.strftime('%Y-%m-%d %H:%M:%S')

class tds_problem:

    def __init__(self, file_root_name):
        rootgrp_path = os.path.join(os.environ['TDS_ROOT'], file_root_name+'-DDMs.nc') 
        self.rootgrp = Dataset(rootgrp_path, "r", format="NETCDF4")
        metagrp_path = os.path.join(os.environ['TDS_ROOT'], file_root_name+'-metadata.nc') 
        self.metagrp = Dataset(metagrp_path, "r", format="NETCDF4")

    def find_index_meta(self):
        datenum = self.rootgrp.groups[self.group].variables['IntegrationMidPointTime'][self.index]
        datenums_meta = np.array(self.metagrp.variables['IntegrationMidPointTime'])
        for i, datenum_meta in enumerate(datenums_meta):
            if datenum_meta == datenum:
                #print("Foud: {0} | {1}".format(i, datenum_meta))
                return i
                break

    def set_group_index(self, group, index):
        self.group = group
        self.index = index
        self.index_meta = self.find_index_meta()

        r_tx = self.metagrp.groups[self.group].variables['TransmitterPositionX'][self.index]
        r_ty = self.metagrp.groups[self.group].variables['TransmitterPositionY'][self.index]
        r_tz = self.metagrp.groups[self.group].variables['TransmitterPositionZ'][self.index]
        self.r_t = np.array([r_tx, r_ty, r_tz])

        v_tx = self.metagrp.groups[self.group].variables['TransmitterVelocityX'][self.index]
        v_ty = self.metagrp.groups[self.group].variables['TransmitterVelocityY'][self.index]
        v_tz = self.metagrp.groups[self.group].variables['TransmitterVelocityZ'][self.index]
        self.v_t = np.array([v_tx, v_ty, v_tz])

        r_rx = self.metagrp.variables['ReceiverPositionX'][self.index_meta]
        r_ry = self.metagrp.variables['ReceiverPositionY'][self.index_meta]
        r_rz = self.metagrp.variables['ReceiverPositionZ'][self.index_meta]
        self.r_r = np.array([r_rx, r_ry, r_rz])

        v_rx = self.metagrp.variables['ReceiverVelocityX'][self.index_meta]
        v_ry = self.metagrp.variables['ReceiverVelocityY'][self.index_meta]
        v_rz = self.metagrp.variables['ReceiverVelocityZ'][self.index_meta]
        self.v_r = np.array([v_rx, v_ry, v_rz])

        r_spx = self.metagrp.groups[self.group].variables['SpecularPointPositionX'][self.index]
        r_spy = self.metagrp.groups[self.group].variables['SpecularPointPositionY'][self.index]
        r_spz = self.metagrp.groups[self.group].variables['SpecularPointPositionZ'][self.index]
        self.r_sp_tds = np.array([r_spx, r_spy, r_spz])
        self.r_sp = self.find_sp()

        self.lat_sp_tds = self.metagrp.groups[group].variables['SpecularPointLat'][index].data
        self.lon_sp_tds = self.metagrp.groups[group].variables['SpecularPointLon'][index].data

    def find_sp(self):
        r_sp_estimate = self.r_sp_tds
        r_center = np.array(r_sp_estimate)
        for it in range(3):
            r_increment = 100e3/10**(it)
            r_step = 1e3/10**(it)
            x = np.arange(r_center[0]-r_increment, r_center[0]+r_increment, r_step)
            y = np.arange(r_center[1]-r_increment, r_center[1]+r_increment, r_step)
            xx, yy = np.meshgrid(x, y)
            z = earth_b*(1-(xx/earth_a)**2-(yy/earth_a)**2)**(1/2)
            min_path = np.inf
            for x_i, x_val in enumerate(x):
                for y_i, y_val in enumerate(y):
                    z_val = z[y_i,x_i]
                    r = np.array([x_val, y_val, z_val])
                    path = np.linalg.norm(r-self.r_t) + np.linalg.norm(self.r_r-r)
                    if path < min_path:
                        min_path = path
                        r_sp = r
            self.lat_sp = np.arcsin(abs(r_sp[2]/np.linalg.norm(r_sp)))*180/np.pi
            self.lon_sp = np.arctan2(r_sp[1], r_sp[0])*180/np.pi

        return r_sp

    def calculate_delay_increment_seconds(self, delay_pixel):
        delay_in_samples = delay_pixel * self.metagrp.groups[self.group].CodeDelaySpacingSamplesBetweenPixels
        delay_in_seconds = delay_in_samples / self.metagrp.groups[self.group].SamplingFrequency
        delay_of_tracking_point_in_seconds = (self.metagrp.groups[self.group].TrackingOffsetDelayNs \
                - self.metagrp.groups[self.group].variables['SpecularPathRangeOffset'][self.index])*1e-9
        delay_increment_seconds = delay_in_seconds - delay_of_tracking_point_in_seconds
        return delay_increment_seconds

    def calculate_delay_increment_chips(self, delay_pixel):
        chips_per_second = 1.023e6
        return calculate_delay_increment_seconds(delay_pixel)*chips_per_second

    def calculate_doppler_increment(self, doppler_pixel):
        return doppler_pixel*self.metagrp.groups[self.group].DopplerResolution \
                - self.metagrp.groups[self.group].TrackingOffsetDopplerHz

    def show_ddm(self):
        datenum = self.rootgrp.groups[self.group].variables['IntegrationMidPointTime'][self.index]
        ddm = self.rootgrp.groups[self.group].variables['DDM'][self.index].data
        string = str(datenum_to_pytime(float(datenum))) \
            + ' Lat: ' + "{0:.2f}".format(self.lat_sp_tds) \
            + ' Lon: ' + "{0:.2f}".format(self.lon_sp_tds)
        delay_0 = self.calculate_delay_increment_seconds(0)
        delay_1 = self.calculate_delay_increment_seconds(127)
        doppler_0 = self.calculate_doppler_increment(-10)
        doppler_1 = self.calculate_doppler_increment(9)
        print(delay_0,delay_1,doppler_0,doppler_1)
        fig, ax = plt.subplots(1,figsize=(10, 4))
        im = ax.imshow(ddm, cmap='viridis', 
                extent=(delay_0, delay_1, doppler_0, doppler_1), 
                aspect=(20/128)/np.abs(doppler_0/delay_0)
                )
        plt.show(block=False)


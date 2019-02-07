#!/usr/bin/env python

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
from gnssr.utils import *

class tds_data:

    def __init__(self, file_root_name, group, index):
        rootgrp_path = os.path.join(os.environ['TDS_ROOT'], file_root_name+'-DDMs.nc') 
        self.rootgrp = Dataset(rootgrp_path, "r", format="NETCDF4")
        metagrp_path = os.path.join(os.environ['TDS_ROOT'], file_root_name+'-metadata.nc') 
        self.metagrp = Dataset(metagrp_path, "r", format="NETCDF4")

        self.set_group_index(group, index)

    def find_index_meta(self):
        """
        The rootgrp index and the metagrp index do not match, the matching 
        key is the IntegrationMidPointTime. This function finds the 
        corresponding meta_index
        """
        datenum = self.rootgrp.groups[self.group].variables['IntegrationMidPointTime'][self.index]
        datenums_meta = np.array(self.metagrp.variables['IntegrationMidPointTime'])
        for i, datenum_meta in enumerate(datenums_meta):
            if datenum_meta == datenum:
                #print("Foud: {0} | {1}".format(i, datenum_meta))
                return i
                break

    def set_group_index(self, group, index):
        """
        Configures the problem associated to the specified group an index.
        """
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
        # inertial velocity
        self.v_t_i = self.v_t + np.cross(np.array([0,0,earth_inertial_angular_speed]), self.r_t)

        r_rx = self.metagrp.variables['ReceiverPositionX'][self.index_meta]
        r_ry = self.metagrp.variables['ReceiverPositionY'][self.index_meta]
        r_rz = self.metagrp.variables['ReceiverPositionZ'][self.index_meta]
        self.r_r = np.array([r_rx, r_ry, r_rz])

        v_rx = self.metagrp.variables['ReceiverVelocityX'][self.index_meta]
        v_ry = self.metagrp.variables['ReceiverVelocityY'][self.index_meta]
        v_rz = self.metagrp.variables['ReceiverVelocityZ'][self.index_meta]
        self.v_r = np.array([v_rx, v_ry, v_rz])
        # inertial velocity
        self.v_r_i = self.v_r + np.cross(np.array([0,0,earth_inertial_angular_speed]), self.r_r)

        r_spx = self.metagrp.groups[self.group].variables['SpecularPointPositionX'][self.index]
        r_spy = self.metagrp.groups[self.group].variables['SpecularPointPositionY'][self.index]
        r_spz = self.metagrp.groups[self.group].variables['SpecularPointPositionZ'][self.index]
        self.r_sp_tds = np.array([r_spx, r_spy, r_spz])

        self.lat_sp_tds = self.metagrp.groups[group].variables['SpecularPointLat'][index].data
        self.lon_sp_tds = self.metagrp.groups[group].variables['SpecularPointLon'][index].data

        self.sp_incidence_tds = self.metagrp.groups[group].variables['SPIncidenceAngle'][index].data

    def find_sp(self):
        """
        Computes a more precise specular point location using the WGS-84 based on the TDS 
        estimation.

        Returns:
            r_sp (numpy.ndarray size(3): 
                Location of the specular point in ECEF coordinates.
            lat_sp (float):
                Latitude of the specular point.
            lon_sp (float)
                Longitude of the specular point.
        """
        r_sp_estimate = self.r_sp_tds
        for it in range(4):
            r_center = np.array(r_sp_estimate)
            r_increment = 100e3/(10**(it))
            r_step = 1e3/(10**(it))
            x = np.arange(r_center[0]-r_increment, r_center[0]+r_increment, r_step)
            y = np.arange(r_center[1]-r_increment, r_center[1]+r_increment, r_step)
            xx, yy = np.meshgrid(x, y)
            z = earth_semiminor_axis*(1-(xx/earth_semimajor_axis)**2-(yy/earth_semimajor_axis)**2)**(1/2)
            min_path = np.inf
            for x_i, x_val in enumerate(x):
                for y_i, y_val in enumerate(y):
                    z_val = z[y_i,x_i]
                    r = np.array([x_val, y_val, z_val])
                    path = np.linalg.norm(r-self.r_t) + np.linalg.norm(self.r_r-r)
                    if path < min_path:
                        min_path = path
                        r_sp_estimate = r
        r_sp = r_sp_estimate
        lat_sp = np.arcsin(abs(r_sp[2]/np.linalg.norm(r_sp)))*180/np.pi
        lon_sp = np.arctan2(r_sp[1], r_sp[0])*180/np.pi
        return r_sp, lat_sp, lon_sp

    def calculate_delay_increment_seconds(self, delay_pixel):
        """
        Calculates the delay increment in seconds of the specified delay pixel.
        """
        code_delay_spacing_samples_between_pixels = self.metagrp.groups[self.group].CodeDelaySpacingSamplesBetweenPixels
        sampling_frequency = self.metagrp.groups[self.group].SamplingFrequency
        tracking_offset_delay_ns = self.metagrp.groups[self.group].TrackingOffsetDelayNs
        specular_path_range_offset = self.metagrp.groups[self.group].variables['SpecularPathRangeOffset'][self.index]
        self.time_delay_resolution = code_delay_spacing_samples_between_pixels/sampling_frequency

        delay_in_samples = delay_pixel * code_delay_spacing_samples_between_pixels
        delay_in_seconds = delay_in_samples / sampling_frequency
        delay_of_tracking_point_in_seconds = (tracking_offset_delay_ns - specular_path_range_offset)*1e-9
        delay_increment_seconds = delay_in_seconds - delay_of_tracking_point_in_seconds
        return delay_increment_seconds

    def calculate_delay_increment_chips(self, delay_pixel):
        """
        Calculates the delay increment in chips of the specified delay pixel.
        """
        # From https://en.wikipedia.org/wiki/GPS_signals
        #The C/A code, for civilian use, transmits data at 1.023 million chips per 
        #second, whereas the P code, for U.S. military use, transmits at 10.23 million 
        #chips per second.  
        return self.calculate_delay_increment_seconds(delay_pixel)*gps_ca_chips_per_second

    def calculate_doppler_increment(self, doppler_pixel):
        """
        Calculates the doppler increment in HZ of the specified pixel.
        """
        self.doppler_resolution = self.metagrp.groups[self.group].DopplerResolution
        tracking_offset_doppler_hz = self.metagrp.groups[self.group].TrackingOffsetDopplerHz
        return doppler_pixel*self.doppler_resolution - tracking_offset_doppler_hz

    def plot_ddm(self):
        """
        Plots delay Doppler map.
        """
        datenum = self.rootgrp.groups[self.group].variables['IntegrationMidPointTime'][self.index]
        ddm = self.rootgrp.groups[self.group].variables['DDM'][self.index].data

        string = str(datenum_to_pytime(float(datenum))) \
            + ' Lat: ' + "{0:.2f}".format(self.lat_sp_tds) \
            + ' Lon: ' + "{0:.2f}".format(self.lon_sp_tds)

        number_of_delay_pixels = self.metagrp.groups[self.group].NumberOfDelayPixels
        number_of_doppler_pixels = self.metagrp.groups[self.group].NumberOfDopplerPixels

        delay_start = self.calculate_delay_increment_chips(0)
        delay_end = self.calculate_delay_increment_chips(number_of_delay_pixels-1)
        doppler_start = self.calculate_doppler_increment(-np.floor(number_of_doppler_pixels/2))
        doppler_end = self.calculate_doppler_increment(np.floor(number_of_doppler_pixels/2 - 0.5))

        fig, ax = plt.subplots(1,figsize=(10, 4))
        ax.set_ylabel('Hz')
        ax.set_xlabel('C/A chips')
        im = ax.imshow(ddm, cmap='viridis', 
                extent=(delay_start, delay_end, doppler_start, doppler_end), 
                aspect=(number_of_doppler_pixels/number_of_delay_pixels)/np.abs(doppler_start/delay_start)
                )
        t = plt.text(0.01, 0.80, string, {'color': 'w', 'fontsize': 12}, transform=ax.transAxes)
        plt.show(block=False)

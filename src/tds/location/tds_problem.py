#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import os

from find_index_meta import * 

class tds_problem:
    def __init__(self, file_root_name):
        rootgrp_path = os.path.join(os.environ['TDS_ROOT'], file_root_name+'-DDMs.nc') 
        self.rootgrp = Dataset(rootgrp_path, "r", format="NETCDF4")
        metagrp_path = os.path.join(os.environ['TDS_ROOT'], file_root_name+'-metadata.nc') 
        self.metagrp = Dataset(metagrp_path, "r", format="NETCDF4")

    def set_group_index(self, group, index):
        self.group = group
        self.index = index
        self.index_meta = find_index_meta(self.rootgrp, self.metagrp, self.group, self.index)

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

        ellip_norm = lambda x,y,z: np.array([2*x/a**2,2*y/a**2,2*z/b**2])

    def find_sp(self):
        a = 6378137 # meters
        b = 6356752.314245 # meters
        r_sp_estimate = self.r_sp_tds
        r_center = np.array(r_sp_estimate)
        for it in range(3):
            r_inc = 100e3/10**(it)
            r_step = 1e3/10**(it)
            x = np.arange(r_center[0]-r_inc, r_center[0]+r_inc, r_step)
            y = np.arange(r_center[1]-r_inc, r_center[1]+r_inc, r_step)
            xx, yy = np.meshgrid(x, y)
            z = b*(1-(xx/a)**2-(yy/a)**2)**(1/2)
            min = np.inf
            for x_i, x_val in enumerate(x):
                for y_i, y_val in enumerate(y):
                    z_val = z[y_i,x_i]
                    r = np.array([x_val, y_val, z_val])
                    d = np.linalg.norm(r-self.r_t) + np.linalg.norm(self.r_r-r)
                    if d < min:
                        min = d
                        r_min = r
            r_sol = r_min
            self.lat_sp = np.arcsin(abs(r_sol[2]/np.linalg.norm(r_sol)))*180/np.pi
            self.lon_sp = np.arctan2(r_sol[1],r_sol[0])*180/np.pi
            r_center = r_sol
        return r_center

#!/usr/bin/env python

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import griddata
import xml.etree.ElementTree as ET

class tds_antenna_gain:

    def __init__(self):
        antenna_path = os.path.join(os.environ['TDS_ROOT'], 'antenna/gain_map.xml') 
        tree = ET.parse(antenna_path)
        root = tree.getroot()

        self.azimuth_pixels = 0
        for iter in root.iter('AzimuthPixels'):
            self.azimuth_pixels = int(iter.text) + 1

        self.elevation_pixels = 0
        for iter in root.iter('ElevationPixels'):
            self.elevation_pixels = int(iter.text)

        azimuths = np.linspace(-180, 180, self.azimuth_pixels) 
        elevations = np.linspace(90, 0, self.elevation_pixels) 
        azimuth_mesh, elevation_mesh = np.meshgrid(azimuths, elevations)
        self.points = np.array([azimuth_mesh.flatten(), elevation_mesh.flatten()]).T

        self.gain_map = np.zeros((self.elevation_pixels, self.azimuth_pixels))

        for index, element in enumerate(root.iter('float')):
            elevation_index = int(index/(self.azimuth_pixels-1))
            azimuth_index = index%(self.azimuth_pixels-1)
            self.gain_map[elevation_index][azimuth_index] = np.exp(float(element.text)/20)
        # The data contains the values for azimuth from -180 to 179.
        # The last column (azimuth=180) is aproximiated with the 
        # previous column (azimuth=179)
        self.gain_map[:,-1] = self.gain_map[:,-2] 

    def gain(self, azimuth, elevation):
        #np.place(azimuth, azimuth>179, 179)
        return griddata(self.points, self.gain_map.flatten(), (azimuth, elevation))

def surface_gain(x,y,tds_antenna):
    h_r = 680e3
    elevation = np.arctan2(h_r,np.sqrt(x**2+y**2))*180/np.pi
    azimuth = np.arctan2(-y,-x)*180/np.pi
    return tds_antenna.gain(azimuth, elevation)

def main():
    antenna = tds_antenna_gain()
    extent = 1500e3
    extent_x0 =  -extent
    extent_x1 =  extent
    extent_y0 =  -extent
    extent_y1 =  extent
    linsapce_delta = 500
    X, Y = np.meshgrid(
            np.linspace(extent_x0, extent_x1, linsapce_delta), 
            np.linspace(extent_y0, extent_y1, linsapce_delta)
            )
    Z_gain = surface_gain(X,Y,antenna)

    fig_surface, ax_surface = plt.subplots(1,figsize=(10, 4))
    plt.xlabel('[km]')
    plt.ylabel('[km]')
    ax_surface.set_title('Gain Map')
    cs_surface = ax_surface.contourf(X, Y, Z_gain, 55, cmap='jet')
    ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/1000))
    ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000))
    ax_surface.xaxis.set_major_formatter(ticks_x)
    ax_surface.yaxis.set_major_formatter(ticks_y)
    ticks_z = ticker.FuncFormatter(lambda x, pos: '{:.2f}'.format(20*np.log10(x)))
    cbar = fig_surface.colorbar(cs_surface, label='Gain [dBi]')
    cbar.formatter= ticks_z
    cbar.update_ticks()

    A, E = np.meshgrid(
            np.linspace(-180, 180, 100),
            np.linspace(90, 0, 100) 
            )

    Z_gain_elev_azi = antenna.gain(A, E)

    fig_gain, ax_gain = plt.subplots(1,figsize=(10, 4))
    plt.xlabel('azimuth [deg]')
    plt.ylabel('elevation [deg]')
    ax_gain.set_title('Gain Map')
    ax_gain.contourf(A, E, Z_gain_elev_azi, 25, cmap='jet')

    plt.show()

if __name__ == '__main__':
    main()

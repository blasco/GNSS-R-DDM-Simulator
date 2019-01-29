#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import xml.etree.ElementTree as ET
tree = ET.parse('antenna_gain_map.xml')
root = tree.getroot()

h_r = 680e3

def surface_gain(x,y):
    elevation = np.arctan2(h_r,np.sqrt(x**2+y**2))*180/np.pi
    azimuth = np.arctan2(-y,-x)*180/np.pi
    return receiver_antenna_gain(elevation, azimuth)

def receiver_antenna_gain(elevation, azimuth):
    elevations = np.linspace(90, 0, elevation_pixels) 
    azimuths = np.linspace(-180, 180 - 1, azimuth_pixels) 
    elevation_mesh, azimuth_mesh = np.meshgrid(elevations, azimuths)
    points = np.array([azimuth_mesh.flatten(), elevation_mesh.flatten()]).T
    return griddata(points, gain_map.flatten(), (elevation,azimuth))

azimuth_pixels = 0
for iter in root.iter('AzimuthPixels'):
    azimuth_pixels = int(iter.text)

elevation_pixels = 0
for iter in root.iter('ElevationPixels'):
    elevation_pixels = int(iter.text)

gain_map = np.zeros((elevation_pixels, azimuth_pixels))

for index, element in enumerate(root.iter('float')):
    elevation_index = int(index/azimuth_pixels)
    azimuth_index = index%azimuth_pixels
    gain_map[elevation_index][azimuth_index] = float(element.text)

extent = 100e3
extent_x0 =  -extent
extent_x1 =  extent
extent_y0 =  -extent
extent_y1 =  extent
linsapce_delta = 500
X, Y = np.meshgrid(
        np.linspace(extent_x0, extent_x1, linsapce_delta), 
        np.linspace(extent_y0, extent_y1, linsapce_delta)
        )
#surface_gain(10,10)
Z_gain = surface_gain(X,Y)

fig_gain, ax_gain = plt.subplots(1,figsize=(10, 4))
plt.xlabel('[km]')
plt.ylabel('[km]')
ax_gain.set_title('Gain Map')
cs_time = ax_gain.contourf(X, Y, Z_gain, 25, cmap='viridis')

E, A = np.meshgrid(
        np.linspace(90, 0, 100), 
        np.linspace(-180, 180, 100)
        )

Z_gain_elev_azi = receiver_antenna_gain(E,A)

fig_gain, ax_gain = plt.subplots(1,figsize=(10, 4))
plt.xlabel('azimuth [deg]')
plt.ylabel('elevation [deg]')
ax_gain.set_title('Gain Map')
cs_time = ax_gain.contourf(E, A, Z_gain_elev_azi, 25, cmap='viridis')

plt.show()

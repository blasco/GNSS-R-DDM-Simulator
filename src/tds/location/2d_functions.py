#!/usr/bin/env python

import numpy as np
import os
from shapely import geometry
import matplotlib.pyplot as plt
from netCDF4 import Dataset

from find_contour_intersections import *
from tds_problem import *

def unit_vector(r):
    return r/np.linalg.norm(r)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2' """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.dot(v1_u, v2_u))

earth_a = 6378137 # meters
earth_b = 6356752.314245 # meters

tds = tds_problem('raw/L1B/2018-07-29-H18')
group = '000050'
index = 380
tds.set_group_index(group, index)

ellip_norm = lambda r: np.array([2*r[0]/earth_a**2,2*r[1]/earth_a**2,2*r[2]/earth_b**2])
n_z = unit_vector(ellip_norm(tds.r_sp))
n_x = unit_vector(np.cross(n_z, tds.r_sp-tds.r_t))
n_y = np.cross(n_z, n_x)

h = np.dot((tds.r_r - tds.r_sp), n_z)
elev = angle_between(n_y, tds.r_t-tds.r_sp)

def time_eq(x, y):
    c = 299792458 # m/s
    return (1/c)*((x**2 +(y+h/np.tan(elev))**2 + h**2)**(1/2) - h/np.sin(elev) - y*np.cos(elev))

def doppler_eq(x, y):
    # GPS L1 center frequency
    c = 299792458 # m/s
    f_c = 1575.42e6 # Hz 
    v_ty = np.dot(tds.v_t, n_y)
    v_tz = np.dot(tds.v_t, n_z)
    v_rx = np.dot(tds.v_r, n_x)
    v_ry = np.dot(tds.v_r, n_y)
    v_rz = np.dot(tds.v_r, n_z)
    return (-f_c/c)*(-v_ty*np.cos(elev) - v_tz*np.sin(elev) + (v_rx*x + v_ry*(y+h/np.tan(elev)) - v_rz*h)/(x**2 + (y+h/np.tan(elev))**2 + h**2)**(1/2))

def doppler_inc_eq(x, y):
    return doppler_eq(x,y) - doppler_eq(0,0)

extent_x0 =  -50e3
extent_x1 =  50e3
extent_y0 =  -50e3
extent_y1 =  50e3
linsapce_delta = 500
X, Y = np.meshgrid(
        np.linspace(extent_x0, extent_x1, linsapce_delta), 
        np.linspace(extent_y0, extent_y1, linsapce_delta)
        )
Z_time = time_eq(X,Y)
Z_doppler = doppler_inc_eq(X,Y)

# Iso-Delay Contour
fig_time, ax_time = plt.subplots(1,figsize=(10, 4))
ax_time.set_title('Eq Time')
cs_time = ax_time.contourf(X, Y, Z_time, 25, cmap='viridis')
fig_time.colorbar(cs_time)
ax_time.grid(c='k', ls='-', alpha=0.3)
plt.show(block=False)

# Iso-Doppler Contour
fig_doppler, ax_doppler = plt.subplots(1,figsize=(10, 4))
ax_doppler.set_title('Eq doppler')
cs_doppler = ax_doppler.contourf(X, Y, Z_doppler, 25, cmap='viridis')
fig_doppler.colorbar(cs_doppler)
ax_doppler.grid(c='k', ls='-', alpha=0.3)
plt.show(block=False)

# Iso-Delay and Iso-Doppler Contour
fig_contour, ax_contour = plt.subplots(1,figsize=(10, 4))
c_time = ax_contour.contour(X, Y, Z_time,
        cmap='viridis'
        )
fig_contour.colorbar(c_time)
c_doppler = ax_contour.contour(X, Y, Z_doppler, 
        cmap='plasma'
        )
fig_contour.colorbar(c_doppler)

# Iso-Delay and Iso-Doppler Lines
time_delay = 2500
c_time = ax_time.contour(X, Y, Z_time, [time_delay], colors='r', linestyles='dashdot')
fig, ax = plt.subplots(1,figsize=(10, 4))
contour_time = ax.contour(X, Y, Z_time, [time_delay],
        colors='k', 
        linestyles='solid',
        extent=(extent_x0, extent_x1, extent_y0, extent_y1)
        )
plt.show(block=False)

doppler_delay = -1.3750
c_doppler = ax_doppler.contour(X, Y, Z_doppler, [doppler_delay], colors='r', linestyles='dashdot')
contour_doppler = ax.contour(X, Y, Z_doppler, [doppler_delay],
        colors='b', 
        linestyles='dashdot',
        extent=(extent_x0, extent_x1, extent_y0, extent_y1)
        )
plt.show(block=False)

# Intersections
intersections = find_intersection(contour_time, contour_doppler)
try:
    for i in intersections:
        ax.scatter(i.x, i.y)
except TypeError as te:
    print ('No intersections')
plt.show(block=False)

plt.show()

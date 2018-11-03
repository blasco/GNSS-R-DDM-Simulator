#!/usr/bin/env python

import numpy as np
import os
from shapely import geometry
import matplotlib.pyplot as plt
from netCDF4 import Dataset

from find_contour_intersections import *
from tds_problem import *
from gnssr.utils import *

file_root_path = 'raw/L1B/2017-11-13-H18'
group = '000066'
index = 375

search_lat_deg = 29.10793
search_lon_deg = -87.94369

time_delay = 5e-6
doppler_delay = 1300

tds = tds_problem(file_root_path)
tds.set_group_index(group, index)
tds.show_ddm()

n_z = unit_vector(ellip_norm(tds.r_sp))
n_x = unit_vector(np.cross(n_z, tds.r_sp-tds.r_t))
n_y = unit_vector(np.cross(n_z, n_x))

h = np.dot((tds.r_r - tds.r_sp), n_z)
h_0 = np.dot((tds.r_t - tds.r_sp), n_z)
elev = angle_between(n_y, tds.r_t-tds.r_sp)
print("elve: {0} elev: {1}".format(elev,  angle_between(-n_y, tds.r_r-tds.r_sp)))
print("h: {0} h_0: {1}".format(h, h_0))

def z_sp(x, y):
    R = np.linalg.norm(tds.r_sp)
    return (R**2 -x**2 -y**2 )**(1/2) - R

def time_eq(x, y):
    c = 299792458 # m/s
    return (1/c)*( \
           (x**2 + (y-h_0/np.tan(elev))**2 + (z_sp(x,y)-h_0)**2)**(1/2) + \
           (x**2 + (-h/np.tan(elev) -y)**2 + (h - z_sp(x,y))**2)**(1/2) \
           )
    #return (1/c)*((x**2 +(y+h/np.tan(elev))**2 + h**2)**(1/2) - h/np.sin(elev) - y*np.cos(elev))

def time_inc_eq(x, y):
    return time_eq(x,y) - time_eq(0,0)

def doppler_eq(x, y):
    # GPS L1 center frequency
    c = 299792458 # m/s
    f_c = 1575.42e6 # Hz 
    v_tx = np.dot(tds.v_t, n_x)
    v_ty = np.dot(tds.v_t, n_y)
    v_tz = np.dot(tds.v_t, n_z)
    v_rx = np.dot(tds.v_r, n_x)
    v_ry = np.dot(tds.v_r, n_y)
    v_rz = np.dot(tds.v_r, n_z)
    return (f_c/c)*( \
            (v_tx*(x)  + v_ty*(y-h_0/np.tan(elev)) + v_tz*(z_sp(x,y)-h_0))  / (x**2 + (y-h_0/np.tan(elev))**2 + (z_sp(x,y)-h_0)**2)**(1/2) \
           -(v_rx*(-x) + v_ry*(-h/np.tan(elev)-y)  + v_rz*(h-z_sp(x,y))   ) / (x**2 + (-h/np.tan(elev) -y)**2 + (h - z_sp(x,y))**2)**(1/2) \
            )
    #return (-f_c/c)*(-v_ty*np.cos(elev) - v_tz*np.sin(elev) + (v_rx*x + v_ry*(y+h/np.tan(elev)) - v_rz*h)/(x**2 + (y+h/np.tan(elev))**2 + h**2)**(1/2))

def doppler_inc_eq(x, y):
    return doppler_eq(x,y) - doppler_eq(0,0)

extent = 80e3
extent_x0 =  -extent
extent_x1 =  extent
extent_y0 =  -extent
extent_y1 =  extent
linsapce_delta = 500
X, Y = np.meshgrid(
        np.linspace(extent_x0, extent_x1, linsapce_delta), 
        np.linspace(extent_y0, extent_y1, linsapce_delta)
        )
Z_time = time_inc_eq(X,Y)
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
c_time = ax_contour.contour(X, Y, Z_time, 25, cmap='winter')
fig_contour.colorbar(c_time)
c_doppler = ax_contour.contour(X, Y, Z_doppler, 25, cmap='autumn')
fig_contour.colorbar(c_doppler)

# Iso-Delay and Iso-Doppler Lines
c_time = ax_time.contour(X, Y, Z_time, [time_delay], colors='r', linestyles='dashdot')
fig, ax = plt.subplots(1,figsize=(10, 4))
contour_time = ax.contour(X, Y, Z_time, [time_delay],
        colors='k', 
        linestyles='solid',
        extent=(extent_x0, extent_x1, extent_y0, extent_y1)
        )
plt.show(block=False)

c_doppler = ax_doppler.contour(X, Y, Z_doppler, [doppler_delay], colors='r', linestyles='dashdot')
contour_doppler = ax.contour(X, Y, Z_doppler, [doppler_delay],
        colors='b', 
        linestyles='dashdot',
        extent=(extent_x0, extent_x1, extent_y0, extent_y1)
        )
plt.show(block=False)

# Searched point
search_lat = search_lat_deg*np.pi/180
search_lon = search_lon_deg*np.pi/180

r_search_x = np.linalg.norm(tds.r_sp)*np.cos(search_lat)*np.cos(search_lon)
r_search_y = np.linalg.norm(tds.r_sp)*np.cos(search_lat)*np.sin(search_lon)
r_search_z = np.linalg.norm(tds.r_sp)*np.sin(search_lat)
r_search =  np.array([r_search_x, r_search_y, r_search_z])
x_search = np.dot(n_x,(r_search - tds.r_sp))
y_search = np.dot(n_y,(r_search - tds.r_sp))

ax.scatter(x_search, y_search)
lon = np.arctan2(r_search[1],r_search[0])*180/np.pi
lat = np.arcsin(abs(r_search[2]/np.linalg.norm(r_search)))*180/np.pi
print("lat search: {0} lon search: {1}".format(lat, lon))

# Intersections
intersections = find_intersection(contour_time, contour_doppler)
try:
    for i in intersections:
        ax.scatter(i.x, i.y)
        r_sol = tds.r_sp + n_x*i.x + n_y*i.y
        lon = np.arctan2(r_sol[1],r_sol[0])*180/np.pi
        lat = np.arcsin(abs(r_sol[2]/np.linalg.norm(r_sol)))*180/np.pi
        print("lat: {0} lon: {1}".format(lat, lon))

except TypeError as te:
    print ('No intersections')
plt.show(block=False)

plt.show()

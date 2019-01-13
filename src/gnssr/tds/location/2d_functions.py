#!/usr/bin/env python

import numpy as np
import os
from shapely import geometry
import tkinter
import matplotlib.pyplot as plt
from netCDF4 import Dataset

from find_contour_intersections import *
from gnssr.tds.tds_data import *
from gnssr.utils import *
from gnssr.targets import *
from gnssr.tds.search_target.cdf4_search import *

# Di Simone
file_root_name = 'raw/L1B/2015-04-01-H00'
target = targets['hibernia']
group = '000095'
index = 525

# 0.5 deg error approx 55 km error
if index == 0:
    search_error = 0.7 
    cdf4_search(file_root_name, target, search_error)

search_lat_deg = target.lat
search_lon_deg = target.lon

time_delay = 3.910308326579947e-0620/1.023e6
doppler_delay = 1500

tds = tds_data(file_root_name)
tds.set_group_index(group, index)
tds.plot_ddm()

print("delay feature: {0}".format(tds.calculate_delay_increment_seconds(80)))

print("t vel: {}".format(np.linalg.norm(tds.v_t)))

r_sp, lat_sp, lon_sp = tds.find_sp();
n_z = unit_vector(ellip_norm(r_sp))
n_x = unit_vector(np.cross(n_z, r_sp-tds.r_t))
n_y = unit_vector(np.cross(n_z, n_x))

#TODO: imporve sp auto calculation
h = np.dot((tds.r_r - r_sp), n_z)
print("h: {}".format(h))
h_0 = np.dot((tds.r_t - r_sp), n_z)
print("h_0: {}".format(h_0))
elev = angle_between(n_y, tds.r_t-r_sp)
print("elev: {}".format(elev))

print("elev tds: {}".format(90 -tds.sp_incidence_tds))
print("elve: {0} elev: {1}".format(elev*180/np.pi,  angle_between(-n_y, tds.r_r-r_sp)*180/np.pi))
print("h: {0} h_0: {1}".format(h, h_0))

print("Transmitter: ")
print("r_tx: {0} \nr_ty: {1} \nr_tz: {2}".format(tds.r_t[0], tds.r_t[1], tds.r_t[2]))
print("v_tx: {0} \nv_ty: {1} \nv_tz: {2}".format(tds.v_t[0], tds.v_t[1], tds.v_t[2]))

print("Receiver: ")
print("r_rx: {0} \nr_ry: {1} \nr_rz: {2}".format(tds.r_r[0], tds.r_r[1], tds.r_r[2]))
print("v_rx: {0} \nv_ry: {1} \nv_rz: {2}".format(tds.v_r[0], tds.v_r[1], tds.v_r[2]))

def z_sp(x, y):
    R = np.linalg.norm(r_sp)
    print("r: {}".format(R))
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
    #f_c = 1575.42e6 # Hz 
    f_0 = 10.23e6 # Hz
    f_c = 154*f_0;
    v_tx = np.dot(tds.v_t, n_x)
    v_ty = np.dot(tds.v_t, n_y)
    v_tz = np.dot(tds.v_t, n_z)
    v_rx = np.dot(tds.v_r, n_x)
    v_ry = np.dot(tds.v_r, n_y)
    v_rz = np.dot(tds.v_r, n_z)
    print("v: {0}, {1}, {2}".format(v_rx,v_ry,v_rz))
    return (f_c/c)*( \
            (v_tx*(x)  + v_ty*(y-h_0/np.tan(elev)) + v_tz*(z_sp(x,y)-h_0))  / (x**2 + (y-h_0/np.tan(elev))**2 + (z_sp(x,y)-h_0)**2)**(1/2) \
           -(v_rx*(-x) + v_ry*(-h/np.tan(elev)-y)  + v_rz*(h-z_sp(x,y))   ) / (x**2 + (-h/np.tan(elev) -y)**2 + (h - z_sp(x,y))**2)**(1/2) \
            )
    #return (-f_c/c)*(-v_ty*np.cos(elev) - v_tz*np.sin(elev) + (v_rx*x + v_ry*(y+h/np.tan(elev)) - v_rz*h)/(x**2 + (y+h/np.tan(elev))**2 + h**2)**(1/2))

def doppler_inc_eq(x, y):
    return doppler_eq(x,y) - doppler_eq(0,0)

extent = 50e3
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
plt.xlabel('[m]')
plt.ylabel('[m]')
ax_time.set_title('Time Delay Contour')
cs_time = ax_time.contourf(X, Y, Z_time, 25, cmap='viridis')
fig_time.colorbar(cs_time, label='seconds')
ax_time.grid(c='k', ls='-', alpha=0.3)
plt.show(block=False)

# Iso-Doppler Contour
fig_doppler, ax_doppler = plt.subplots(1,figsize=(10, 4))
ax_doppler.set_title('Doppler Shift Contour')
plt.xlabel('[m]')
plt.ylabel('[m]')
cs_doppler = ax_doppler.contourf(X, Y, Z_doppler, 25, cmap='viridis')
fig_doppler.colorbar(cs_doppler, label='Hz')
ax_doppler.grid(c='k', ls='-', alpha=0.3)
plt.show(block=False)

# Iso-Delay and Iso-Doppler Lines
fig_lines, ax_lines = plt.subplots(1,figsize=(10, 4))
ax_lines.set_title('Iso-Delay and Iso-Doppler Lines')
plt.xlabel('[m]')
plt.ylabel('[m]')

print("res t: {}".format(tds.time_delay_resolution))
print("res h: {}".format(tds.doppler_resolution))

iso_delay_values = list(np.arange(0,tds.time_delay_resolution*50,tds.time_delay_resolution))
c_time = ax_lines.contour(X, Y, Z_time, iso_delay_values, cmap='winter')
fig_lines.colorbar(c_time, label='seconds')
iso_doppler_values = list(np.arange(-6000,6000,tds.doppler_resolution))
c_doppler = ax_lines.contour(X, Y, Z_doppler, iso_doppler_values, cmap='autumn')
fig_lines.colorbar(c_doppler, label='Hz')

# Iso-Delay and Iso-Doppler Line
c_time = ax_time.contour(X, Y, Z_time, [time_delay], colors='r', linestyles='dashdot')
fig_line, ax_line = plt.subplots(1,figsize=(10, 4))
plt.xlabel('[m]')
plt.ylabel('[m]')
contour_time = ax_line.contour(X, Y, Z_time, [time_delay],
        colors='k', 
        linestyles='solid',
        extent=(extent_x0, extent_x1, extent_y0, extent_y1)
        )
ax_lines.contour(X, Y, Z_time, [time_delay],
        colors='r', 
        linewidths = 2,
        linestyles='dashed',
        extent=(extent_x0, extent_x1, extent_y0, extent_y1)
        )
plt.show(block=False)

c_doppler = ax_doppler.contour(X, Y, Z_doppler, [doppler_delay], colors='r', linestyles='dashdot')
contour_doppler = ax_line.contour(X, Y, Z_doppler, [doppler_delay],
        colors='b', 
        linestyles='dashdot',
        extent=(extent_x0, extent_x1, extent_y0, extent_y1)
        )
ax_lines.contour(X, Y, Z_doppler, [doppler_delay],
        colors='r', 
        linewidths = 2,
        linestyles='dashed',
        extent=(extent_x0, extent_x1, extent_y0, extent_y1)
        )
plt.show(block=False)

# Searched point
search_lat = search_lat_deg*np.pi/180
search_lon = search_lon_deg*np.pi/180

r_search_x = np.linalg.norm(r_sp)*np.cos(search_lat)*np.cos(search_lon)
r_search_y = np.linalg.norm(r_sp)*np.cos(search_lat)*np.sin(search_lon)
r_search_z = np.linalg.norm(r_sp)*np.sin(search_lat)
r_search =  np.array([r_search_x, r_search_y, r_search_z])
x_search = np.dot(n_x,(r_search - r_sp))
y_search = np.dot(n_y,(r_search - r_sp))

x_sp_tds = np.dot(n_x, (tds.r_sp_tds - r_sp))
y_sp_tds = np.dot(n_y, (tds.r_sp_tds - r_sp))
ax_line.scatter(x_sp_tds, x_sp_tds, s=70, marker=(5, 2))
ax_lines.scatter(x_sp_tds, y_sp_tds, s=70, marker=(5, 2))

ax_line.scatter(x_search, y_search, s=70)
ax_lines.scatter(x_search, y_search, s=70)
lon = np.arctan2(r_search[1],r_search[0])*180/np.pi
lat = np.arcsin(abs(r_search[2]/np.linalg.norm(r_search)))*180/np.pi
print("lat search: {0} lon search: {1}".format(lat, lon))

# Intersections
intersections = find_intersection(contour_time, contour_doppler)
try:
    for i in intersections:
        ax_line.scatter(i.x, i.y, s=70)
        ax_lines.scatter(i.x, i.y, s=70)
        r_i = np.array([i.x, i.y])
        r_target = np.array([x_search, y_search])
        print("error: {0}".format(np.linalg.norm(r_target - r_i)/1e3))
        r_sol = r_sp + n_x*i.x + n_y*i.y
        lon = np.arctan2(r_sol[1],r_sol[0])*180/np.pi
        lat = np.arcsin(abs(r_sol[2]/np.linalg.norm(r_sol)))*180/np.pi
        print("lat: {0} lon: {1}".format(lat, lon))

except TypeError as te:
    print ('No intersections')
plt.show(block=False)

plt.show()

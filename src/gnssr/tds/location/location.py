#!/usr/bin/env python

import numpy as np
import os
from shapely import geometry
import tkinter
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from netCDF4 import Dataset

from gnssr.tds.location.find_contour_intersections import *
from gnssr.tds.tds_data import *
from gnssr.utils import *
from gnssr.targets import *
from gnssr.tds.search_database.cdf4_search import *

def main():

    #file_root_name = '2018-03-31-H06'
    #target = targets['hibernia']
    #group = '000021'
    #index = 300 

    file_root_name = '2015-12-04-H18'
    target = targets['devils_tower']
    group = '000066'
    index =  385

    #File: /home/woowapdabug/projects/thesis/python/src/tds/raw/L1B_Catalogue/2018-03/31/H06/2018-03.31.H06.kmz

    # 0.5 deg error approx 55 km error
    if index == 0:
        search_error = 0.7 
        cdf4_search(file_root_name, target, search_error)

    target_delay_increment = 14
    target_doppler_increment = -3000

    tds = tds_data(file_root_name, group, index)
    tds.set_group_index(group, index)
    tds.plot_ddm()
    plt.xlim([-5,15])

    print("t vel: {}".format(np.linalg.norm(tds.v_t)))

    r_sp, lat_sp, lon_sp = tds.find_sp();
    n_z = unit_vector(ellip_norm(r_sp))
    n_x = unit_vector(np.cross(n_z, r_sp-tds.r_t))
    n_y = unit_vector(np.cross(n_z, n_x))

    #TODO: imporve sp auto calculation
    h_r = np.dot((tds.r_r - r_sp), n_z)
    print("h_r: {}".format(h_r))
    h_t = np.dot((tds.r_t - r_sp), n_z)
    print("h_t: {}".format(h_t))
    elev = angle_between(n_y, tds.r_t-r_sp)
    print("elev: {}".format(elev))

    print("elev tds: {}".format(90 -tds.sp_incidence_tds))
    print("elve: {0} elev: {1}".format(elev*180/np.pi,  angle_between(-n_y, tds.r_r-r_sp)*180/np.pi))
    print("h_r: {0} h_t: {1}".format(h_r, h_t))

    def z_sp(x, y):
        R = np.linalg.norm(r_sp)
        print("r: {}".format(R))
        return (R**2 -x**2 -y**2 )**(1/2) - R

    light_speed = 299792458 # m/s
    def time_inc_eq(x, y):
        return (np.sqrt(x ** 2 + (y + h_r / np.tan(elev)) ** 2 + h_r ** 2) - h_r / np.sin(elev) - y * np.cos(elev)) / light_speed

    def doppler_eq(x, y):
        # GPS L1 center frequency
        f_0 = 10.23e6 # Hz
        f_carrier = 154*f_0; # 1575.42e6 Hz 
        v_tx = np.dot(tds.v_t, n_x)
        v_ty = np.dot(tds.v_t, n_y)
        v_tz = np.dot(tds.v_t, n_z)
        v_rx = np.dot(tds.v_r, n_x)
        v_ry = np.dot(tds.v_r, n_y)
        v_rz = np.dot(tds.v_r, n_z)
        # With the very far way transmitter approximation:
        f_D_0 =  -f_carrier / light_speed * (-v_ty * np.cos(elev) - v_tz * np.sin(elev) + (v_rx * x + v_ry * (y + h_r / np.tan(elev)) - v_rz * h_r) * (x ** 2 + (y + h_r / np.tan(elev)) ** 2 + h_r ** 2) ** (-0.1e1 / 0.2e1))
        return f_D_0

    def doppler_inc_eq(x, y):
        return doppler_eq(x,y) - doppler_eq(0,0)

    extent = 200e3
    extent_x0 =  -extent
    extent_x1 =  extent
    extent_y0 =  -extent
    extent_y1 =  extent
    linsapce_delta = 500
    X, Y = np.meshgrid(
            np.linspace(extent_x0, extent_x1, linsapce_delta), 
            np.linspace(extent_y0, extent_y1, linsapce_delta)
            )
    delay_chip =  1/1.023e6 # seconds
    Z_time_chip = time_inc_eq(X,Y)/delay_chip
    Z_doppler = doppler_inc_eq(X,Y)

    fig_surface, ax_surface = plt.subplots(1,figsize=(10, 4))
    plt.xlabel('[km]')
    plt.ylabel('[km]')
    ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000))
    ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/1000))
    ax_surface.xaxis.set_major_formatter(ticks_x)
    ax_surface.yaxis.set_major_formatter(ticks_y)

    # Iso-Delay Contour
    ax_surface.set_title('Time Delay Contour')
    contour_delay = ax_surface.contour(X, Y, Z_time_chip, np.arange(0,15,0.5), cmap='viridis')
    fig_surface.colorbar(contour_delay, label='C/A chips')
    ax_surface.grid(c='k', ls='-', alpha=0.3)
    plt.show(block=False)

    # Iso-Doppler Contour
    ax_surface.set_title('Doppler Shift Contour')
    plt.xlabel('[km]')
    plt.ylabel('[km]')
    contour_doppler = ax_surface.contour(X, Y, Z_doppler, np.arange(-5000,5000,500), cmap='viridis')
    fig_surface.colorbar(contour_doppler, label='Hz')
    ax_surface.grid(c='k', ls='-', alpha=0.3)
    #plt.xlim([-20e3,90e3])
    #plt.ylim([-20e3,90e3])

    target_iso_delay = ax_surface.contour(X, Y, Z_time_chip, [target_delay_increment-0.3],
            colors='red', 
            linewidths = 2.5,
            linestyles='dashed',
            extent=(extent_x0, extent_x1, extent_y0, extent_y1)
            )
    target_iso_delay = ax_surface.contour(X, Y, Z_time_chip, [target_delay_increment+0.3],
            colors='red', 
            linewidths = 2.5,
            linestyles='dashed',
            extent=(extent_x0, extent_x1, extent_y0, extent_y1)
            )
    target_iso_doppler = ax_surface.contour(X, Y, Z_doppler, [target_doppler_increment-250],
            colors='red', 
            linewidths = 2.5,
            linestyles='dashed',
            extent=(extent_x0, extent_x1, extent_y0, extent_y1)
            )
    target_iso_doppler = ax_surface.contour(X, Y, Z_doppler, [target_doppler_increment+250],
            colors='red', 
            linewidths = 2.5,
            linestyles='dashed',
            extent=(extent_x0, extent_x1, extent_y0, extent_y1)
            )

    # TDS specular point in local coordinates
    x_sp_tds = np.dot(n_x, (tds.r_sp_tds - r_sp))
    y_sp_tds = np.dot(n_y, (tds.r_sp_tds - r_sp))

    ax_surface.scatter(x_sp_tds, x_sp_tds, s=70, marker=(5, 2), zorder=4, color='green')

    # Target in ECEF coordinates
    target_lat = target.lat_deg*np.pi/180
    target_lon = target.lon_deg*np.pi/180
    r_target_x = np.linalg.norm(r_sp)*np.cos(target_lat)*np.cos(target_lon)
    r_target_y = np.linalg.norm(r_sp)*np.cos(target_lat)*np.sin(target_lon)
    r_target_z = np.linalg.norm(r_sp)*np.sin(target_lat)
    r_target =  np.array([r_target_x, r_target_y, r_target_z])

    # Target in local coordinates
    x_target = np.dot(n_x,(r_target - r_sp))
    y_target = np.dot(n_y,(r_target - r_sp))

    ax_surface.scatter(x_target, y_target, s=70, zorder=4, color='black')

    lon = np.arctan2(r_target[1],r_target[0])*180/np.pi
    lat = np.arcsin(abs(r_target[2]/np.linalg.norm(r_target)))*180/np.pi
    print("lat target: {0} lon target: {1}".format(lat, lon))

    ## Intersection of target iso-delay and iso-doppler intersections
    #intersections = find_contour_intersection(target_iso_delay, target_iso_doppler)
    #try:
    #    for i in intersections:
    #        ax_surface.scatter(i.x, i.y, s=70, zorder=4, color='orange')
    #        r_i = np.array([i.x, i.y])
    #        r_target = np.array([x_target, y_target])
    #        print("error: {0}".format(np.linalg.norm(r_target - r_i)/1e3))
    #        r_sol = r_sp + n_x*i.x + n_y*i.y
    #        lon = np.arctan2(r_sol[1],r_sol[0])*180/np.pi
    #        lat = np.arcsin(abs(r_sol[2]/np.linalg.norm(r_sol)))*180/np.pi
    #        print("lat: {0} lon: {1}".format(lat, lon))

    #except TypeError as te:
    #    print ('No intersections')

    plt.show(block=False)

    plt.show()

if __name__ == '__main__':
    main()

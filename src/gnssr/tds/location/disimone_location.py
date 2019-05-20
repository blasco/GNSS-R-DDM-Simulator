#!/usr/bin/env python
from gnssr.utils import *
from gnssr.simulator.isolines import *
import gnssr.simulator.jacobian.spherical as spherical

from gnssr.simulator.simulation_configuration import *
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
    # Di Simone
    file_root_name = '2015-04-01-H00'
    target = targets['hibernia']
    group = '000095'
    index = 525

    # Hibernia
    # Really good detection
    #file_root_name = 'raw/L1B/2017-03-12-H18'
    #target = targets['hibernia']
    #group = '000035'
    #index = 681

    # Devils Tower
    #file_root_name = 'raw/L1B/2015-12-04-H18'
    #target = targets['devils_tower']
    #group = '000066'
    #index =  376

    #file_root_name = 'raw/L1B/2016-11-14-H06'
    #target = targets['petronius']
    #group = '000057'
    #index = 0 

    # 0.5 deg error approx 55 km error
    if index == 0:
        search_error = 0.7 
        cdf4_search(file_root_name, target, search_error)

    target_delay_increment = 4
    target_doppler_increment = 2000

    tds = tds_data(file_root_name, group, index)
    tds.set_group_index(group, index)
    tds.plot_ddm()

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
        f_D_0 =  f_carrier / light_speed * (-v_ty * np.cos(elev) - v_tz * np.sin(elev) + (v_rx * x + v_ry * (y + h_r / np.tan(elev)) - v_rz * h_r) * (x ** 2 + (y + h_r / np.tan(elev)) ** 2 + h_r ** 2) ** (-0.1e1 / 0.2e1))
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

    target_iso_delay = ax_surface.contour(X, Y, Z_time_chip, [target_delay_increment],
            colors='red', 
            linewidths = 2.5,
            linestyles='dashed',
            extent=(extent_x0, extent_x1, extent_y0, extent_y1)
            )
    target_iso_doppler = ax_surface.contour(X, Y, Z_doppler, [target_doppler_increment],
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

    # Intersection of target iso-delay and iso-doppler intersections
    intersections = find_contour_intersection(target_iso_delay, target_iso_doppler)
    try:
        for i in intersections:
            ax_surface.scatter(i.x, i.y, s=70, zorder=4, color='orange')
            r_i = np.array([i.x, i.y])
            r_target = np.array([x_target, y_target])
            print("error: {0}".format(np.linalg.norm(r_target - r_i)/1e3))
            r_sol = r_sp + n_x*i.x + n_y*i.y
            lon = np.arctan2(r_sol[1],r_sol[0])*180/np.pi
            lat = np.arcsin(abs(r_sol[2]/np.linalg.norm(r_sol)))*180/np.pi
            print("lat: {0} lon: {1}".format(lat, lon))

    except TypeError as te:
        print ('No intersections')

    
    n_z = unit_vector(ellip_norm(tds.r_sp_tds))
    n_x = unit_vector(np.cross(n_z, tds.r_sp_tds-tds.r_t))
    n_y = unit_vector(np.cross(n_z, n_x))

    v_tx = np.dot(tds.v_t, n_x)
    v_ty = np.dot(tds.v_t, n_y)
    v_tz = np.dot(tds.v_t, n_z)

    v_rx = np.dot(tds.v_r, n_x)
    v_ry = np.dot(tds.v_r, n_y)
    v_rz = np.dot(tds.v_r, n_z)

    sim_config = simulation_configuration()
    sim_config.set_scenario_local_ref(
        h_r = np.dot((tds.r_r - tds.r_sp_tds), n_z),
        h_t = np.dot((tds.r_t - tds.r_sp_tds), n_z),
        elevation = angle_between(n_y, tds.r_t-tds.r_sp_tds),
        v_t = np.array([v_tx,v_ty,v_tz]),
        v_r = np.array([v_rx,v_ry,v_rz])
        )

    doppler_specular_point =  eq_doppler_absolute_shift(np.array([0,0,0]), sim_config)
    delay = np.array([4*delay_chip])
    f_doppler = np.array([doppler_specular_point + 1800])
    r_1, r_2 = spherical.delay_doppler_to_local_surface(delay, f_doppler, sim_config)
    print(r_1)
    print(r_2)
    ax_surface.scatter(r_1[0], r_1[1], s=70, zorder=4, color='blue')
    ax_surface.scatter(r_2[0], r_2[1], s=70, zorder=4, color='blue')

    plt.show(block=False)

    plt.show()

if __name__ == '__main__':
    main()

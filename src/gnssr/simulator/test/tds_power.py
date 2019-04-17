#!/usr/bin/env python

import numpy as np

import matplotlib.pyplot as plt

from gnssr.simulator.ddm import *
from gnssr.simulator.simulation_configuration import *

from gnssr.targets import *
from gnssr.tds.tds_data import *
from gnssr.tds.detection.find_targets import *
from gnssr.utils import *
from gnssr.simulator.waf import *

import cv2

def rescale(ddm_original, n_row_res, n_col_res):
    n_row, n_col = ddm_original.shape 
    assert n_row > n_row_res, "Cannot rescale to a biger size"
    assert n_col > n_col_res, "Cannot rescale to a biger size"
    assert n_col % n_col_res == 0, "low res should be a multiple"
    assert n_row % n_row_res == 0, "low res should be a multiple"
    n_row_res = int(n_row/int(n_row/n_row_res))
    n_col_res = int(n_col/int(n_col/n_col_res))
    ddm_res = np.zeros((n_row_res, n_col_res))
    for row_i, row in enumerate(ddm_original):
        for col_i, val in enumerate(row):
            row_i_res = int(row_i/(n_row/n_row_res))
            col_i_res = int(col_i/(n_col/n_col_res))
            ddm_res[row_i_res,col_i_res] += val
    return ddm_res

def main():

    sim_config = simulation_configuration()
    sim_config.jacobian_type = 'spherical'

    delay_chip = sim_config.delay_chip

    file_root_name = '2017-03-12-H18'
    target = targets['hibernia']
    group = '000035'
    index = 675

    tds = tds_data(file_root_name, group, index)

    mean_wind = tds.get_wind()
    p = target_processor_power();
    n = 30
    p.n = n
    for i in range(index - n, index + 2):
        tds.set_group_index(group, i)
        ddm_i = normalize(tds.rootgrp.groups[group].variables['DDM'][i].data)*tds.peak_power()
        p.process_ddm(ddm_i)
        wind = tds.get_wind() 
        if (wind != None):
            mean_wind += wind
            mean_wind /= 2
    ddm_tds = np.copy(p.sea_clutter)
    print("mean wind: {0}".format(mean_wind))

    sim_config.u_10 = 9.5 
    sim_config.phi_0 = -120*np.pi/180

    # Plot TDS DDM sea clutter
    tds.set_group_index(group, index)
    r_sp, lat_sp, lon_sp = tds.find_sp();
    datenum = tds.rootgrp.groups[tds.group].variables['IntegrationMidPointTime'][tds.index]

    string = str(datenum_to_pytime(float(datenum))) \
        + ' Lat: ' + "{0:.2f}".format(tds.lat_sp_tds) \
        + ' Lon: ' + "{0:.2f}".format(tds.lon_sp_tds)

    tds_number_of_delay_pixels = tds.metagrp.groups[tds.group].NumberOfDelayPixels
    tds_number_of_doppler_pixels = tds.metagrp.groups[tds.group].NumberOfDopplerPixels

    tds_delay_start = tds.calculate_delay_increment_chips(0)
    tds_delay_end = tds.calculate_delay_increment_chips(tds_number_of_delay_pixels-1)
    tds_delay_resolution = (tds_delay_end-tds_delay_start)/128

    tds_doppler_start = tds.calculate_doppler_increment(-np.floor(tds_number_of_doppler_pixels/2))
    tds_doppler_end = tds.calculate_doppler_increment(np.floor(tds_number_of_doppler_pixels/2 - 0.5))
    tds_doppler_resolution = 500

    fig_tds, ax_tds = plt.subplots(1,figsize=(10, 4))
    plt.title('TDS-1 Experimental Data')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    contour_tds = ax_tds.imshow(ddm_tds, cmap='jet', 
            extent=(tds_delay_start, tds_delay_end, tds_doppler_end, tds_doppler_start), 
            aspect=(tds_number_of_doppler_pixels/tds_number_of_delay_pixels)/np.abs(tds_doppler_start/tds_delay_start)
            )
    t = plt.text(0.01, 0.85, string, {'color': 'w', 'fontsize': 12}, transform=ax_tds.transAxes)
    cbar = fig_tds.colorbar(contour_tds, label='Power')

    # Load TDS Geometry in simulation configuration
    n_z = unit_vector(ellip_norm(tds.r_sp_tds))
    n_x = unit_vector(np.cross(n_z, tds.r_sp_tds-tds.r_t))
    n_y = unit_vector(np.cross(n_z, n_x))

    v_tx = np.dot(tds.v_t, n_x)
    v_ty = np.dot(tds.v_t, n_y)
    v_tz = np.dot(tds.v_t, n_z)

    v_rx = np.dot(tds.v_r, n_x)
    v_ry = np.dot(tds.v_r, n_y)
    v_rz = np.dot(tds.v_r, n_z)

    sim_config.set_scenario_local_ref(
        h_r = np.dot((tds.r_r - tds.r_sp_tds), n_z),
        h_t = np.dot((tds.r_t - tds.r_sp_tds), n_z),
        elevation = angle_between(n_y, tds.r_t-tds.r_sp_tds),
        v_t = np.array([v_tx,v_ty,v_tz]),
        v_r = np.array([v_rx,v_ry,v_rz])
        )

    # DDM
    sim_config.delay_increment_start = tds_delay_start*delay_chip
    sim_config.delay_increment_end = tds_delay_end*delay_chip
    sim_config.delay_resolution = (sim_config.delay_increment_end - sim_config.delay_increment_start)/tds_number_of_delay_pixels/3
    sim_config.doppler_increment_start = tds_doppler_start
    sim_config.doppler_increment_end = tds_doppler_end + tds_doppler_resolution
    sim_config.doppler_resolution = (sim_config.doppler_increment_end - sim_config.doppler_increment_start)/tds_number_of_doppler_pixels/3

    ddm_sim = (simulate_ddm(sim_config))

    fig_ddm, ax_ddm = plt.subplots(1,figsize=(10, 4))
    plt.title('DDM original simulation')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    contour_ddm = ax_ddm.imshow(ddm_sim, cmap='jet', 
            extent=(sim_config.delay_increment_start, sim_config.delay_increment_end, sim_config.doppler_increment_end, sim_config.doppler_increment_start), 
            aspect="auto"
            )
    cbar = fig_ddm.colorbar(contour_ddm, label='Normalized Power')

    # Image downscaling to desired resolution:
    fig_ddm_rescaled, ax_ddm_rescaled = plt.subplots(1,figsize=(10, 4))
    plt.title('Simulation Rescaled')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')

    ddm_rescaled = rescale(ddm_sim, tds_number_of_doppler_pixels, tds_number_of_delay_pixels) 

    # Noise
    waf_delay_increment_values = list(np.arange(
        -tds_delay_end*delay_chip, 
        tds_delay_end*delay_chip + tds_delay_resolution*delay_chip, 
        tds_delay_resolution*delay_chip
        ))
    waf_doppler_increment_values = list(np.arange(
        tds_doppler_start, 
        tds_doppler_end + tds_doppler_resolution, 
        tds_doppler_resolution
        ))
    waf_delay_grid, waf_doppler_grid = np.meshgrid(waf_delay_increment_values, waf_doppler_increment_values)
    waf_matrix = woodward_ambiguity_function(waf_delay_grid, waf_doppler_grid, sim_config)**2

    T_noise_receiver = 200
    k_b = 1.38e-23 # J/K
    y_noise = 1/sim_config.coherent_integration_time*k_b*T_noise_receiver

    p1 = target_processor_power();
    n = 30000 
    p1.n = n
    p1.tau = 0.08
    ddm_noise = np.zeros(ddm_rescaled.shape)
    for i in range(n+1):
        print("i: {0}".format(i))
        noise_i = y_noise*(np.random.rand(ddm_rescaled.shape[0], ddm_rescaled.shape[1])-0.5)
        ddm_noise_i = np.abs(signal.convolve2d(noise_i, waf_matrix, mode='same'))
        p1.process_ddm(np.abs(ddm_rescaled + ddm_noise_i))
    ddm_rescaled = p1.sea_clutter

    contour_res = ax_ddm_rescaled.imshow(ddm_rescaled, cmap='jet', 
            extent=(tds_delay_start, tds_delay_end, tds_doppler_end, tds_doppler_start), 
            aspect=(tds_number_of_doppler_pixels/tds_number_of_delay_pixels)/np.abs(tds_doppler_start/tds_delay_start)
            )
    cbar = fig_ddm_rescaled.colorbar(contour_res, label='Normalized Power', shrink=0.35)

    fig_diff, ax_diff = plt.subplots(1,figsize=(10, 4))
    plt.title('Difference')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')

    ddm_diff = np.copy(ddm_rescaled)
    for row_i, row in enumerate(ddm_diff):
        for col_i, val in enumerate(row):
            val_tds = ddm_tds[row_i,col_i]
            val = ddm_rescaled[row_i,col_i]
            ddm_diff[row_i,col_i] = (val-val_tds)/val_tds
    #np.place(ddm_diff, ddm_diff < 0, np.nan)

    im = ax_diff.imshow(ddm_diff, cmap='jet', 
            extent=(tds_delay_start, tds_delay_end, tds_doppler_end, tds_doppler_start), 
            aspect=(tds_number_of_doppler_pixels/tds_number_of_delay_pixels)/np.abs(tds_doppler_start/tds_delay_start)
            )
    cbar = fig_diff.colorbar(im, label='Normalized Power', shrink=0.35)

    plt.show()

if __name__ == '__main__':
    main()

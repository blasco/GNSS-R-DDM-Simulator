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
from gnssr.simulator.noise import *

import cv2

def main():

    sim_config = simulation_configuration()
    sim_config.jacobian_type = 'spherical'

    delay_chip = sim_config.delay_chip

    # Decent noise results
    #file_root_name = '2017-03-12-H18'
    #target = targets['hibernia']
    #group = '000035'
    #index = 675

    file_root_name = '2017-03-12-H18'
    target = targets['hibernia']
    group = '000035'
    index = 675-35

    tds = tds_data(file_root_name, group, index)

    '''
    p = target_processor_power();
    p.n = n
    for i in range(index - n, index + 2):
        tds.set_group_index(group, i)
        ddm_i = normalize(tds.rootgrp.groups[group].variables['DDM'][i].data)*tds.peak_power()
        p.process_ddm(ddm_i)
        wind = tds.get_wind() 
        if (wind != None):
            print("wind: {0}".format(wind))
            mean_wind += wind
            mean_wind /= 2
    ddm_tds = np.copy(p.sea_clutter)
    print("mean wind: {0}".format(mean_wind))
    '''

    n_tds = 20 
    n = n_tds 
    ddm_tds = np.zeros(tds.rootgrp.groups[group].variables['DDM'][index].data.shape)
    mean_wind = 0
    n_wind = 0
    noise_antenna_mean = 0
    noise_rx_mean = 0
    for i in range(n):
        print("i: {0}".format(i))

        tds.set_group_index(group, index + i)
        power_i, noise_i = tds.peak_power()
        ddm_i = normalize(tds.rootgrp.groups[group].variables['DDM'][index + i].data)*power_i
        print("noise power estimate: {0}".format(noise_i))
        ddm_tds += ddm_i

        wind = tds.get_wind() 

        #noise_antenna_mean_i = tds.metagrp.groups[group].variables['AntennaTemperature'][tds.index].data
        #if(np.isnan(noise_antenna_mean_i)):
        #    noise_antenna_mean_i = tds.metagrp.groups[group].variables['AntennaTemperatureExtRef'][tds.index].data
        #    if(np.isnan(noise_antenna_mean_i)):
        #        noise_antenna_mean_i = 246.85

        #noise_rx_mean_i = tds.metagrp.variables['RxTemperature'][tds.index_meta].data
        #if(np.isnan(noise_rx_mean_i)):
        #    noise_rx_mean_i = 232.09

        #noise_antenna_mean += noise_antenna_mean_i
        #noise_rx_mean += noise_rx_mean_i
        #print("noise antenna : {}".format(noise_antenna_mean_i))

        if (wind != None):
            print("wind: {0}".format(wind))
            mean_wind += wind
            n_wind += 1

    mean_wind /= n_wind
    ddm_tds /= n
    #noise_antenna_mean /= n 
    #noise_rx_mean /= n

    noise_antenna_mean = 246.85
    noise_rx_mean = 232.09

    print("noise_antenna_mean: {0}".format(noise_antenna_mean))
    print("noise_rx_mean: {0}".format(noise_rx_mean))
    print("mean wind: {0}".format(mean_wind))

    sim_config.u_10 = 11 
    sim_config.phi_0 = 90*np.pi/180 # max diff: 0.64
    #sim_config.phi_0 = 0*np.pi/180 # max diff:   0.62

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
            #aspect=(tds_number_of_doppler_pixels/tds_number_of_delay_pixels)/np.abs(tds_doppler_start/tds_delay_start)
            aspect='auto'
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


    h_r = np.dot((tds.r_r - tds.r_sp_tds), n_z)
    h_t = np.dot((tds.r_t - tds.r_sp_tds), n_z)
    elevation = angle_between(n_y, tds.r_t-tds.r_sp_tds)
    v_t = np.array([v_tx,v_ty,v_tz])
    v_r = np.array([v_rx,v_ry,v_rz])


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
    cbar = fig_ddm.colorbar(contour_ddm, label='Correlated Power [Watts]')

    # DDM Rescaled
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

    noise_temperature = noise_antenna_mean + noise_rx_mean
    ddm_rescaled = add_thermal_noise(tds_number_of_doppler_pixels, tds_number_of_delay_pixels, 1*1000, noise_temperature, ddm_rescaled, sim_config)

    contour_res = ax_ddm_rescaled.imshow(ddm_rescaled, cmap='jet', 
            extent=(tds_delay_start, tds_delay_end, tds_doppler_end, tds_doppler_start), 
            #aspect=(tds_number_of_doppler_pixels/tds_number_of_delay_pixels)/np.abs(tds_doppler_start/tds_delay_start)
            aspect='auto'
            )
    cbar = fig_ddm_rescaled.colorbar(contour_res, label='Correlated Power [Watts]' )

    # Difference

    fig_diff, ax_diff = plt.subplots(1,figsize=(10, 4))
    plt.title('Difference')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')

    ddm_diff = np.copy(ddm_rescaled)

    for row_i, row in enumerate(ddm_diff):
        for col_i, val in enumerate(row):
            col_i_shift = col_i
            row_i_shift = row_i
            col_shift = 0 
            row_shift = 1 
            if (col_i + col_shift) >= 0 and (col_i + col_shift) < 128:
                col_i_shift += col_shift
            if (row_i + row_shift) >= 0 and (row_i + row_shift) < 20:
                row_i_shift += row_shift
            val_tds = ddm_tds[row_i_shift,col_i_shift]
            val = ddm_rescaled[row_i,col_i]
            if row_i == 0 :
                ddm_diff[row_i,col_i] = 0
            elif col_i == 0:
                ddm_diff[row_i,col_i] = 0
            else:
                rel = np.abs((val-val_tds)/val_tds)
                ddm_diff[row_i,col_i] = rel
    ddm_diff[:,-1] = 0
    ddm_diff[0,0] = 0
    ddm_diff[0,1] = 1

    im = ax_diff.imshow(ddm_diff, cmap='jet', 
            extent=(tds_delay_start, tds_delay_end, tds_doppler_end, tds_doppler_start), 
            #aspect=(tds_number_of_doppler_pixels/tds_number_of_delay_pixels)/np.abs(tds_doppler_start/tds_delay_start)
            aspect='auto'
            )
    cbar = fig_diff.colorbar(im, label='Relative error' )

    # SNR

    fig_snr, ax_snr = plt.subplots(1,figsize=(10, 4))
    plt.title('SNR')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')

    T_noise_receiver = noise_antenna_mean + noise_rx_mean
    k_b = 1.38e-23 # J/K
    y_noise = 4/sim_config.coherent_integration_time*k_b*T_noise_receiver/np.sqrt(1000)

    ddm_rescaled_snr = rescale(ddm_sim, tds_number_of_doppler_pixels, tds_number_of_delay_pixels) 
    ddm_snr = np.copy(10*np.log10(np.abs(ddm_rescaled_snr)/y_noise))
    contour_snr = ax_snr.imshow(ddm_snr, cmap='jet', 
            extent=(tds_delay_start, tds_delay_end, tds_doppler_end, tds_doppler_start), 
            #aspect=(tds_number_of_doppler_pixels/tds_number_of_delay_pixels)/np.abs(tds_doppler_start/tds_delay_start)
            aspect='auto'
            )
    cbar = fig_snr.colorbar(contour_snr, label='dB' )

    print(file_root_name)
    print("h_r: {}".format(h_r))
    print("h_t: {}".format(h_t))
    print("elevation: {}".format(elevation*180/np.pi))
    print("v_r: {}".format(np.linalg.norm(v_r)))
    print("v_t: {}".format(np.linalg.norm(v_t)))
    print("max diff: {}".format(np.max(ddm_diff[1:,1:])))

    plt.show()

if __name__ == '__main__':
    main()

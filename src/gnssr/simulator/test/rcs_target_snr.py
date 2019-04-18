#!/usr/bin/env python

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import gnssr.simulator.rcs.target_rcs as target_rcs
import gnssr.simulator.rcs.sea_rcs as sea_rcs
from gnssr.simulator.isolines import *
from gnssr.simulator.simulation_configuration import *
from gnssr.simulator.ddm import *

import cv2

def main():

    sim_config = simulation_configuration()
    #sim_config.jacobian_type = 'spherical'

    sim_config.receiver_antenna_gain = lambda p1,p2: 12.589

    sim_config.set_scenario_local_ref(
            h_t = 20e6, # m
            h_r = 20e3, # meters
            elevation = 60.0*np.pi/180,
            v_t = np.array([-2684.911, 1183.799, -671.829]), # m/s
            v_r = np.array([20, 20, 20]) # m/s
            )

    #sim_config.rcs = sea_rcs.radar_cross_section
    sim_config.rcs = lambda p1,p2: target_rcs.radar_cross_section(p1, 0.0, p2)
    sim_config.u_10 = 8.0 # m/s

    #sim_config.delay_chip = 1/gps_ca_chips_per_second # seconds
    delay_chip = sim_config.delay_chip

    sim_config.doppler_increment_start = -100
    sim_config.doppler_increment_end = 100
    sim_config.delay_increment_start = -0.2*delay_chip
    sim_config.delay_increment_end = 3*delay_chip
    sim_config.coherent_integration_time = 2e-2 # sec

    delay_increment_start = sim_config.delay_increment_start 
    delay_increment_end = sim_config.delay_increment_end 
    delay_resolution = sim_config.delay_resolution

    doppler_increment_start = sim_config.doppler_increment_start
    doppler_increment_end = sim_config.doppler_increment_end
    doppler_resolution = sim_config.doppler_resolution

    doppler_specular_point = eq_doppler_absolute_shift(np.array([0,0,0]), sim_config)

    number_of_delay_pixels = 128 - 50
    number_of_doppler_pixels = 20 + 50
    rescaled_doppler_resolution = (doppler_increment_end - doppler_increment_start)/number_of_doppler_pixels
    rescaled_delay_resolution_chips = (delay_increment_end - delay_increment_start)/number_of_delay_pixels
    sim_config.delay_resolution = rescaled_delay_resolution_chips/4
    sim_config.doppler_resolution = rescaled_doppler_resolution/4

    print("doppler res: {0}".format(rescaled_doppler_resolution))
    print("delay res: {0}".format(rescaled_delay_resolution_chips))

    # Surface mesh
    x_0 = 0
    x_1 = 6e3 # meters
    n_x = 800

    y_0 = -1e3
    y_1 = 6e3 # meters
    n_y = 800

    x_grid, y_grid = np.meshgrid(
       np.linspace(x_0, x_1, n_x), 
       np.linspace(y_0, y_1, n_y)
       )

    r = np.array([x_grid, y_grid, 0])

    # Isolines and RCS
    z_grid_delay_chip = eq_delay_incremet(r, sim_config)/delay_chip

    doppler_specular_point = eq_doppler_absolute_shift(np.array([0,0,0]), sim_config)
    z_grid_doppler_increment = eq_doppler_absolute_shift(r, sim_config) - doppler_specular_point

    z_rcs = sim_config.rcs(r, sim_config)

    # Plot
    fig_rcs, ax_rcs = plt.subplots(1,figsize=(10, 4))

    #contour_delay_chip = ax_rcs.contour(
    #        x_grid, y_grid, z_grid_delay_chip, 
    #        np.arange(0, 2, 0.1), 
    #        cmap='winter', alpha = 0.3
    #        )
    #contour_doppler = ax_rcs.contour(
    #        x_grid, y_grid, z_grid_doppler_increment, 
    #        np.arange(-50, 50, 1), 
    #        cmap='jet', alpha = 0.3
    #        )
    contour_rcs = ax_rcs.contourf(x_grid, y_grid, z_rcs, 55, cmap='jet', alpha = 0.8)

    ax_rcs.set_title('RCS')
    plt.xlabel('[km]')
    plt.ylabel('[km]')
    #fig_rcs.colorbar(contour_delay_chip, label='C/A chips')
    #fig_rcs.colorbar(contour_doppler, label='Hz')
    fig_rcs.colorbar(contour_rcs, label='Gain')

    target_delay_increment = 0.54
    target_doppler_increment = 17.35

    target_delay = 1.26
    target_doppler = 22
    target_iso_delay = ax_rcs.contour(x_grid, y_grid, z_grid_delay_chip, 
            #[target_delay_increment-0.1],
            [target_delay - rescaled_delay_resolution_chips/delay_chip],
            colors='red', 
            linewidths = 2.5,
            linestyles='dashed',
            )
    target_iso_delay = ax_rcs.contour(x_grid, y_grid, z_grid_delay_chip, 
            #[target_delay_increment+0.1],
            [target_delay + rescaled_delay_resolution_chips/delay_chip],
            colors='red', 
            linewidths = 2.5,
            linestyles='dashed',
            )
    target_iso_delay = ax_rcs.contour(x_grid, y_grid, z_grid_doppler_increment, 
            #[target_doppler_increment-0.5],
            [target_doppler - rescaled_doppler_resolution],
            colors='red', 
            linewidths = 2.5,
            linestyles='dashed',
            )
    target_iso_delay = ax_rcs.contour(x_grid, y_grid, z_grid_doppler_increment, 
            #[target_doppler_increment+0.5],
            [target_doppler + rescaled_doppler_resolution],
            colors='red', 
            linewidths = 2.5,
            linestyles='dashed',
            )

    ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/1000))
    ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000))
    ax_rcs.xaxis.set_major_formatter(ticks_x)
    ax_rcs.yaxis.set_major_formatter(ticks_y)

    # DDM Noise
    T_noise_receiver = 225
    k_b = 1.38e-23 # J/K
    y_noise = 1/sim_config.coherent_integration_time*k_b*T_noise_receiver
    print("y_noise: {0}".format(y_noise))

    # DDM 
    ddm_sim = np.copy(simulate_ddm(sim_config))
    ddm_sim_snr  = np.copy(10*np.log10(np.abs(ddm_sim)/y_noise))
    sim_config_1 = sim_config
    sim_config_1.rcs = sea_rcs.radar_cross_section
    #sim_config.rcs = lambda p1,p2: target_rcs.radar_cross_section(p1, 0, p2)
    sim_config_1.u_10 = 8.0
    ddm_sim_1 = np.copy(simulate_ddm(sim_config_1))

    ddm_diff = np.copy(np.abs(ddm_sim - ddm_sim_1))
    ddm_diff_snr = np.copy(10*np.log10(np.abs(ddm_sim - ddm_sim_1)/y_noise))

    ddm_sim_res = rescale(ddm_sim, number_of_doppler_pixels, number_of_delay_pixels)
    ddm_sim_1_res = rescale(ddm_sim_1, number_of_doppler_pixels, number_of_delay_pixels)
    ddm_diff_res = np.abs(ddm_sim_res - ddm_sim_1_res)

    # Plotting
    fig_ddm, ax_ddm = plt.subplots(1,figsize=(10, 4))
    plt.title('DDM simulation')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    contour_sim = ax_ddm.imshow(ddm_sim, cmap='jet', 
            extent=(
                delay_increment_start/delay_chip, delay_increment_end/delay_chip, 
                doppler_increment_end, doppler_increment_start), 
            aspect="auto"
            )
    fig_ddm.colorbar(contour_sim, label='Correlated Power [W]')

    fig_ddm_snr, ax_ddm_snr = plt.subplots(1,figsize=(10, 4))
    plt.title('DDM simulation SNR')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    contour_sim_snr = ax_ddm_snr.imshow(ddm_sim_snr, cmap='jet', 
            extent=(
                delay_increment_start/delay_chip, delay_increment_end/delay_chip, 
                doppler_increment_end, doppler_increment_start), 
            aspect="auto"
            )
    fig_ddm_snr.colorbar(contour_sim_snr, label='Correlated Power [W]')

    fig_diff, ax_diff = plt.subplots(1,figsize=(10, 4))
    plt.title('DDM diff')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    contour_diff = ax_diff.imshow(ddm_diff, cmap='jet', 
            extent=(
                delay_increment_start/delay_chip, delay_increment_end/delay_chip, 
                doppler_increment_end, doppler_increment_start), 
            aspect="auto"
            )
    fig_diff.colorbar(contour_diff, label='SNR [dB]')

    fig_diff_snr, ax_diff_snr = plt.subplots(1,figsize=(10, 4))
    plt.title('DDM diff SNR')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    contour_diff_snr = ax_diff_snr.imshow(ddm_diff_snr, cmap='jet', 
            extent=(
                delay_increment_start/delay_chip, delay_increment_end/delay_chip, 
                doppler_increment_end, doppler_increment_start), 
            aspect="auto"
            )
    fig_diff_snr.colorbar(contour_diff_snr, label='SNR [dB]')

    fig_ddm_rescaled, ax_ddm_rescaled = plt.subplots(1,figsize=(10, 4))
    plt.title('DDM dif Rescaled')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')

    contour_diff_res = ax_ddm_rescaled.imshow(ddm_diff_res, cmap='jet', 
            extent=(
                delay_increment_start/delay_chip, delay_increment_end/delay_chip, 
                doppler_increment_end, doppler_increment_start), 
            aspect='auto'
            )
    fig_ddm_rescaled.colorbar(contour_diff_res, label='Correlated Power [W]')

    plt.show()

if __name__ == '__main__':
    main()

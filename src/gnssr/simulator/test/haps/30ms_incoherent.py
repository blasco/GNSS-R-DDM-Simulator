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

    sim_config.set_scenario_local_ref(
            h_t = 13.82e6, # m
            h_r = 20e3, # meters
            elevation = 70.0*np.pi/180,
            v_t = np.array([-2684.911, 1183.799, -671.829]), # m/s
            v_r = np.array([33, 0, 0]) # m/s
            )

    sim_config.receiver_antenna_gain = lambda p1,p2: 1
    sim_config.rcs = lambda p1,p2: target_rcs.radar_cross_section(p1, 0, p2)
    sim_config.target_x = 5e3;
    sim_config.target_y = 0.5e3;

    sim_config.fresnel = 0.65;
    sim_config.convolve_type = 'fft'

    u_10 = 10.0
    sim_config.u_10 = u_10

    sim_config.delay_chip /= 10
    delay_chip = sim_config.delay_chip

    number_of_delay_pixels = 128 - 50
    number_of_doppler_pixels = 20 + 50
    #number_of_delay_pixels = (128 - 50)*2
    #number_of_doppler_pixels = (20 + 50)*2

    sim_config.doppler_increment_start = -70
    sim_config.doppler_increment_end = 70
    sim_config.doppler_resolution = (sim_config.doppler_increment_end - sim_config.doppler_increment_start)/number_of_doppler_pixels/8
    sim_config.delay_increment_start = -1*delay_chip
    sim_config.delay_increment_end = 40*delay_chip
    #sim_config.delay_resolution = 0.01*delay_chip
    sim_config.delay_resolution = (sim_config.delay_increment_end - sim_config.delay_increment_start)/number_of_delay_pixels/4
    sim_config.coherent_integration_time = 20e-3 # sec

    delay_increment_start = sim_config.delay_increment_start 
    delay_increment_end = sim_config.delay_increment_end 
    delay_resolution = sim_config.delay_resolution

    doppler_increment_start = sim_config.doppler_increment_start
    doppler_increment_end = sim_config.doppler_increment_end
    doppler_resolution = sim_config.doppler_resolution

    doppler_specular_point = eq_doppler_absolute_shift(np.array([0,0,0]), sim_config)

    rescaled_doppler_resolution = (doppler_increment_end - doppler_increment_start)/number_of_doppler_pixels
    rescaled_delay_resolution_chips = (delay_increment_end - delay_increment_start)/delay_chip/number_of_delay_pixels

    print("doppler res: {0}".format(rescaled_doppler_resolution))
    print("delay res: {0}".format(rescaled_delay_resolution_chips))

    # Surface mesh
    x_0 = -8e3 
    x_1 = 8e3 # meters
    n_x = 800

    y_0 = -8e3
    y_1 = 8e3 # meters
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

    contour_delay_chip = ax_rcs.contour(
            x_grid, y_grid, z_grid_delay_chip, 
            np.arange(0, 10, delay_resolution/delay_chip), 
            cmap='jet', alpha = 0.3
            )
    contour_doppler = ax_rcs.contour(
            x_grid, y_grid, z_grid_doppler_increment, 
            np.arange(-70, 70, doppler_resolution), 
            cmap='jet', alpha = 0.3
            )

    contour_delay_chip = ax_rcs.contour(
            x_grid, y_grid, z_grid_delay_chip, 
            np.arange(0, 10, 1), 
            cmap='winter', alpha = 0.4,
            linewidths = 2.5,
            linestyles='dashed',
            )
    contour_doppler = ax_rcs.contour(
            x_grid, y_grid, z_grid_doppler_increment, 
            np.arange(-70, 70, 1), 
            cmap='jet', alpha = 0.4,
            linewidths = 2.5,
            linestyles='dashed',
            )

    contour_rcs = ax_rcs.contourf(x_grid, y_grid, z_rcs, 55, cmap='jet', alpha = 0.8)

    ax_rcs.set_title('RCS')
    plt.xlabel('[km]')
    plt.ylabel('[km]')
    #fig_rcs.colorbar(contour_delay_chip, label='C/A chips')
    #fig_rcs.colorbar(contour_doppler, label='Hz')
    fig_rcs.colorbar(contour_rcs, label='Radar Cross Section')

    target_delay = 2.1
    target_doppler = 29
    target_iso_delay = ax_rcs.contour(x_grid, y_grid, z_grid_delay_chip, 
            #[target_delay_increment-0.1],
            [target_delay - rescaled_delay_resolution_chips],
            colors='red', 
            linewidths = 2.5,
            linestyles='dashed',
            )
    target_iso_delay = ax_rcs.contour(x_grid, y_grid, z_grid_delay_chip, 
            #[target_delay_increment+0.1],
            [target_delay + rescaled_delay_resolution_chips],
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

    # DDM

    ddm_sim = np.copy(simulate_ddm(sim_config))
    sim_config.rcs = sea_rcs.radar_cross_section
    #sim_config.rcs = lambda p1,p2: target_rcs.radar_cross_section(p1, 0, p2)
    sim_config.u_10 = 10.8
    ddm_sim_1 = np.copy(simulate_ddm(sim_config))
    ddm_diff = np.abs(ddm_sim - ddm_sim_1)

    fig_diff, ax_diff = plt.subplots(1,figsize=(10, 4))
    plt.title('DDM diff simulation')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    contour_diff = ax_diff.imshow(ddm_diff, cmap='jet', 
            extent=(
                delay_increment_start/delay_chip, delay_increment_end/delay_chip, 
                doppler_increment_end, doppler_increment_start), 
            aspect="auto"
            )
    fig_diff.colorbar(contour_diff, label='Correlated Power [Watts]')

    fig_ddm, ax_ddm = plt.subplots(1,figsize=(10, 4))
    plt.title('DDM original simulation')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    contour_sim = ax_ddm.imshow(ddm_sim, cmap='jet', 
            extent=(
                delay_increment_start/delay_chip, delay_increment_end/delay_chip, 
                doppler_increment_end, doppler_increment_start), 
            aspect="auto"
            )
    fig_ddm.colorbar(contour_sim, label='Correlated Power [Watts]')

    fig_ddm_1, ax_ddm_1 = plt.subplots(1,figsize=(10, 4))
    plt.title('DDM original simulation 1')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    contour_sim_1 = ax_ddm_1.imshow(ddm_sim_1, cmap='jet', 
            extent=(
                delay_increment_start/delay_chip, delay_increment_end/delay_chip, 
                doppler_increment_end, doppler_increment_start), 
            aspect="auto"
            )
    fig_ddm_1.colorbar(contour_sim_1, label='Correlated Power [Watts]')

    # Image downscaling to desired resolution:
    # TODO: This is just an average of the pixels around the area
    # This is not valid, summation i srequired:
    # Di Simone > From a physical viewpoint, 
    # such an approach should call for summation instead of averaging
    # https://stackoverflow.com/questions/48121916/numpy-resize-rescale-image
    fig_ddm_rescaled_diff, ax_ddm_rescaled_diff = plt.subplots(1,figsize=(10, 4))
    plt.title('Wake to Sea Clutter DDM Difference')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    #ddm_rescaled = rescale(ddm_diff, number_of_doppler_pixels, number_of_delay_pixels)

    ddm_sim_res = rescale(ddm_sim, number_of_doppler_pixels, number_of_delay_pixels)
    ddm_sim_1_res = rescale(ddm_sim_1, number_of_doppler_pixels, number_of_delay_pixels)
    ddm_diff_res = np.abs(ddm_sim_res - ddm_sim_1_res)

    contour_diff = ax_ddm_rescaled_diff.imshow(ddm_diff_res, cmap='jet', 
            extent=(
                delay_increment_start/delay_chip, delay_increment_end/delay_chip, 
                doppler_increment_end, doppler_increment_start), 
            aspect='auto'
            )
    fig_ddm_rescaled_diff.colorbar(contour_diff, label='Correlated Power [Watts]')

    fig_ddm_rescaled, ax_ddm_rescaled = plt.subplots(1,figsize=(10, 4))
    plt.title('Ship Wake DDM Simulation')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')

    contour_rescaled = ax_ddm_rescaled.imshow(ddm_sim_res, cmap='jet', 
            extent=(
                delay_increment_start/delay_chip, delay_increment_end/delay_chip, 
                doppler_increment_end, doppler_increment_start), 
            aspect='auto'
            )
    fig_ddm_rescaled.colorbar(contour_rescaled, label='Correlated Power [Watts]')

    # DDM Noise
    T_noise_receiver = 232.09 + 246.85
    k_b = 1.38e-23 # J/K
    y_noise = 2/sim_config.coherent_integration_time*k_b*T_noise_receiver
    print("expected SNR: {}".format(y_noise))

    fig_snr, ax_snr = plt.subplots(1,figsize=(10, 4))
    plt.title('SNR')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')

    ddm_snr = 10*np.log10(np.abs(ddm_diff_res)/y_noise)
    np.place(ddm_snr, ddm_snr < -30, np.nan)
    contour_snr = ax_snr.imshow(ddm_snr, cmap='jet', 
            extent=(
                delay_increment_start/delay_chip, delay_increment_end/delay_chip, 
                doppler_increment_end, doppler_increment_start), 
            aspect='auto'
            )
    fig_snr.colorbar(contour_snr, label='Correlated Power [Watts]')

    plt.show()

if __name__ == '__main__':
    main()

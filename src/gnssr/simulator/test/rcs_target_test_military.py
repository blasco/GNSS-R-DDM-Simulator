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

def rescale(ddm_original, n_row_res, n_col_res):
    n_row, n_col = ddm_original.shape 
    assert n_row > n_row_res, "Cannot rescale to a biger size"
    assert n_row > n_row_res, "Cannot rescale to a biger size"
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

    sim_config.set_scenario_local_ref(
            h_t = 13.82e6, # m
            h_r = 18e3, # meters
            elevation = 60.0*np.pi/180,
            v_t = np.array([-2684.911, 1183.799, -671.829]), # m/s
            v_r = np.array([25, 25, 25]) # m/s
            )

    sim_config.rcs = target_rcs.radar_cross_section
    sim_config.u_10 = 5 # m/s

    sim_config.delay_chip = 1/10.23e6 # s
    delay_chip = sim_config.delay_chip

    sim_config.doppler_increment_start = -45
    sim_config.doppler_increment_end = 45 
    sim_config.doppler_resolution = 0.3 
    sim_config.delay_increment_start = -0.2*delay_chip
    sim_config.delay_increment_end = 10*delay_chip
    sim_config.delay_resolution = 0.05*delay_chip
    sim_config.coherent_integration_time = 1e-2 # sec

    delay_increment_start = sim_config.delay_increment_start 
    delay_increment_end = sim_config.delay_increment_end 
    delay_resolution = sim_config.delay_resolution

    doppler_increment_start = sim_config.doppler_increment_start
    doppler_increment_end = sim_config.doppler_increment_end
    doppler_resolution = sim_config.doppler_resolution

    doppler_specular_point = eq_doppler_absolute_shift(np.array([0,0,0]), sim_config)

    # Surface mesh
    x_0 =  0
    x_1 =  4e3 # meters
    n_x = 800

    y_0 =  0
    y_1 =  4e3 # meters
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
            np.arange(0, delay_increment_end/delay_chip, delay_resolution/delay_chip), 
            cmap='winter', alpha = 0.4
            )
    contour_doppler = ax_rcs.contour(
            x_grid, y_grid, z_grid_doppler_increment, 
            np.arange(doppler_increment_start, doppler_increment_end, doppler_resolution), 
            cmap='jet', alpha = 0.4
            )
    contour_rcs = ax_rcs.contourf(x_grid, y_grid, z_rcs, 55, cmap='jet', alpha = 0.8)

    ax_rcs.set_title('RCS')
    plt.xlabel('[km]')
    plt.ylabel('[km]')
    #fig_rcs.colorbar(contour_delay_chip, label='C/A chips')
    #fig_rcs.colorbar(contour_doppler, label='Hz')
    fig_rcs.colorbar(contour_rcs, label='Gain')


    target_delay_increment = 0.54
    target_doppler_increment = 17.35
    target_iso_delay = ax_rcs.contour(x_grid, y_grid, z_grid_delay_chip, 
            #[target_delay_increment-0.1],
            [4.5],
            colors='red', 
            linewidths = 2.5,
            linestyles='dashed',
            )
    target_iso_delay = ax_rcs.contour(x_grid, y_grid, z_grid_delay_chip, 
            #[target_delay_increment+0.1],
            [5],
            colors='red', 
            linewidths = 2.5,
            linestyles='dashed',
            )
    target_iso_delay = ax_rcs.contour(x_grid, y_grid, z_grid_doppler_increment, 
            #[target_doppler_increment-0.5],
            [14],
            colors='red', 
            linewidths = 2.5,
            linestyles='dashed',
            )
    target_iso_delay = ax_rcs.contour(x_grid, y_grid, z_grid_doppler_increment, 
            #[target_doppler_increment+0.5],
            [16],
            colors='red', 
            linewidths = 2.5,
            linestyles='dashed',
            )

    ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/1000))
    ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000))
    ax_rcs.xaxis.set_major_formatter(ticks_x)
    ax_rcs.yaxis.set_major_formatter(ticks_y)

    # DDM
    ddm_sim = simulate_ddm(sim_config) 
    sim_config.rcs = sea_rcs.radar_cross_section
    ddm_sim_sea = simulate_ddm(sim_config) 
    ddm_diff = np.abs(ddm_sim - ddm_sim_sea)

    fig_diff, ax_diff = plt.subplots(1,figsize=(10, 4))
    plt.title('DDM diff simulation')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    im = ax_diff.imshow(ddm_diff, cmap='jet', 
            extent=(
                delay_increment_start/delay_chip, delay_increment_end/delay_chip, 
                doppler_increment_end, doppler_increment_start), 
            aspect="auto"
            )

    fig_ddm, ax_ddm = plt.subplots(1,figsize=(10, 4))
    plt.title('DDM original simulation')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    im = ax_ddm.imshow(ddm_sim, cmap='jet', 
            extent=(
                delay_increment_start/delay_chip, delay_increment_end/delay_chip, 
                doppler_increment_end, doppler_increment_start), 
            aspect="auto"
            )

    # Image downscaling to desired resolution:
    # TODO: This is just an average of the pixels around the area
    # This is not valid, summation i srequired:
    # Di Simone > From a physical viewpoint, 
    # such an approach should call for summation instead of averaging
    # https://stackoverflow.com/questions/48121916/numpy-resize-rescale-image
    fig_ddm_rescaled, ax_ddm_rescaled = plt.subplots(1,figsize=(10, 4))
    plt.title('Rescaled Diff')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    number_of_delay_pixels = 128
    number_of_doppler_pixels = 20
    ddm_sim_res = rescale(ddm_sim, number_of_doppler_pixels, number_of_delay_pixels)
    ddm_sea_res = rescale(ddm_sim_sea, number_of_doppler_pixels, number_of_delay_pixels)
    ddm_diff_res = np.abs(ddm_sim_res - ddm_sea_res)

    im = ax_ddm_rescaled.imshow(ddm_diff_res, cmap='jet', 
            extent=(
                delay_increment_start/delay_chip, delay_increment_end/delay_chip, 
                doppler_increment_end, doppler_increment_start), 
            aspect='auto'
            )

    plt.show()

if __name__ == '__main__':
    main()

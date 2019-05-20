#!/usr/bin/env python

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from gnssr.simulator.rcs.sea_rcs import *
from gnssr.simulator.isolines import *
from gnssr.simulator.simulation_configuration import *
from gnssr.simulator.jacobian.planar import *

def main():

    sim_config = simulation_configuration()

    sim_config.set_scenario_local_ref(
            h_t = 13.82e6, # m
            h_r = 20e3, # meters
            elevation = 50.0*np.pi/180,
            v_t = np.array([-2684.911, 1183.799, -671.829]), # m/s
            v_r = np.array([20, 20, 20]) # m/s
            )

    sim_config.u_10 = 10
    sim_config.phi_0 = 0

    delay_chip = sim_config.delay_chip

    delay_increment_start = sim_config.delay_increment_start 
    delay_increment_end = sim_config.delay_increment_end 
    delay_resolution = sim_config.delay_resolution

    doppler_increment_start = sim_config.doppler_increment_start
    doppler_increment_end = sim_config.doppler_increment_end
    doppler_resolution = sim_config.doppler_resolution

    doppler_specular_point = eq_doppler_absolute_shift(np.array([0,0,0]), sim_config)

    # Surface mesh

    x_0 =  -40e3 # meters
    x_1 =  40e3 # meters
    n_x = 500

    y_0 =  -40e3 # meters
    y_1 =  40e3 # meters
    n_y = 500

    x_grid, y_grid = np.meshgrid(
       np.linspace(x_0, x_1, n_x), 
       np.linspace(y_0, y_1, n_y)
       )

    r = np.array([x_grid, y_grid, 0])

    # Isolines and Antenna gain
    z_grid_delay_chip = eq_delay_incremet(r, sim_config)/delay_chip

    doppler_specular_point = eq_doppler_absolute_shift(np.array([0,0,0]), sim_config)
    z_grid_doppler_increment = eq_doppler_absolute_shift(r, sim_config) - doppler_specular_point

    z_rcs = sim_config.rcs(r, sim_config)

    # Iso lines plot
    fig_isolines, ax_isolines = plt.subplots(1,figsize=(10, 4))

    #contour_rcs = ax_isolines.contourf(x_grid, y_grid, z_rcs, 55, cmap='jet', alpha = 0.8)

    contour_delay_chip = ax_isolines.contour(
            x_grid, y_grid, z_grid_delay_chip, 
            np.arange(0, 11, 0.14), 
            cmap='summer', alpha = 0.6
            )
    contour_doppler = ax_isolines.contour(
            x_grid, y_grid, z_grid_doppler_increment, 
            np.arange(-70, 80, 2), 
            cmap='winter', alpha = 0.8
            )
    #fig_isolines.colorbar(contour_rcs, label='Radar Cross Section')
    fig_isolines.colorbar(contour_delay_chip, label='Delay [C/A chips]')
    fig_isolines.colorbar(contour_doppler, label='Doppler [Hz]')

    '''
    test_delay = np.array([3*delay_chip])
    test_doppler = np.array([1000]) +  doppler_specular_point
    print('Finding intersection for d:{0}, f:{1}'.format(test_delay, test_delay))
    x_s_1 = x_delay_doppler_1(test_delay, test_doppler, sim_config)
    y_s_1 = y_delay_doppler_1(test_delay, test_doppler, sim_config)
    x_s_2 = x_delay_doppler_2(test_delay, test_doppler, sim_config)
    y_s_2 = y_delay_doppler_2(test_delay, test_doppler, sim_config)
    ax_isolines.scatter(x_s_1, y_s_1, s=70, marker=(5, 2), zorder=4)
    ax_isolines.scatter(x_s_2, y_s_2, s=70, marker=(5, 2), zorder=4)

    fig_isolines.colorbar(contour_delay_chip, label='C/A chips')
    fig_isolines.colorbar(contour_doppler, label='Hz')
    '''

    ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/1000))
    ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000))
    ax_isolines.xaxis.set_major_formatter(ticks_x)
    ax_isolines.yaxis.set_major_formatter(ticks_y)
    plt.xlabel('[km]')
    plt.ylabel('[km]')


    plt.show()

if __name__ == '__main__':
    main()

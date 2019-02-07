#!/usr/bin/env python

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from gnssr.simulator.rcs.sea_rcs import *
from gnssr.simulator.isolines import *
from gnssr.simulator.simulation_configuration import *

def main():

    sim_config = simulation_configuration()
    sim_config.u_10 = 15 # m/s

    delay_chip = sim_config.delay_chip

    delay_increment_start = sim_config.delay_increment_start 
    delay_increment_end = sim_config.delay_increment_end 
    delay_resolution = sim_config.delay_resolution

    doppler_increment_start = sim_config.doppler_increment_start
    doppler_increment_end = sim_config.doppler_increment_end
    doppler_resolution = sim_config.doppler_resolution

    doppler_specular_point = eq_doppler_absolute_shift(np.array([0,0,0]), sim_config)

    # Surface mesh
    x_0 =  -150e3 # meters
    x_1 =  150e3 # meters
    n_x = 500

    y_0 =  -150e3 # meters
    y_1 =  150e3 # meters
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

    z_rcs = radar_cross_section(r, sim_config)

    # Plot

    fig_antenna, ax_antenna = plt.subplots(1,figsize=(10, 4))

    contour_delay_chip = ax_antenna.contour(
            x_grid, y_grid, z_grid_delay_chip, 
            np.arange(0, delay_increment_end/delay_chip, 1), 
            cmap='winter', alpha = 0.8
            )
    contour_doppler = ax_antenna.contour(
            x_grid, y_grid, z_grid_doppler_increment, 
            np.arange(doppler_increment_start, doppler_increment_end, 500), 
            cmap='jet', alpha = 0.8
            )
    contour_antenna = ax_antenna.contourf(x_grid, y_grid, z_rcs, 55, cmap='jet', alpha = 0.5)

    ax_antenna.set_title('Antenna')
    plt.xlabel('[km]')
    plt.ylabel('[km]')
    fig_antenna.colorbar(contour_delay_chip, label='C/A chips')
    fig_antenna.colorbar(contour_doppler, label='Hz')
    fig_antenna.colorbar(contour_antenna, label='Gain')

    ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/1000))
    ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000))
    ax_antenna.xaxis.set_major_formatter(ticks_x)
    ax_antenna.yaxis.set_major_formatter(ticks_y)

    plt.show()

if __name__ == '__main__':
    main()

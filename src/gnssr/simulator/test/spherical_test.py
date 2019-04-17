#!/usr/bin/env python

import numpy as np
from gnssr.utils import *

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from gnssr.simulator.rcs.sea_rcs import *
from gnssr.simulator.isolines import *
from gnssr.simulator.simulation_configuration import *

import gnssr.simulator.jacobian.spherical as spherical
import gnssr.simulator.jacobian.planar as planar

def main():

    '''
    sim_config = simulation_configuration()

    h     = 540e3; 
    h0    = 20e6; 
    Vt    = np.array([-1e3,1e3,0]); 
    Vr    = np.array([500,1000,1e2]); 
    gamma = 70*np.pi/180;
    wavelength = light_speed/sim_config.f_carrier

    Re     = 6371e3;

    r      = np.array([100, -57, 0]);
    r[2]   = np.sqrt(Re**2 - r[0]**2 - r[1]**2) - Re;

    Rt = np.array([0, h0/np.tan(gamma), h0]); 
    Rr = np.array([0, -h/np.tan(gamma), h]); 

    ni = -(Rt-r)/np.linalg.norm(Rt-r); 
    ns = (Rr-r)/np.linalg.norm(Rr-r); 

    sim_config.set_scenario_local_ref(
            h_t = h0, # m
            h_r = h,  # m
            elevation = gamma, # radians
            v_t = Vt, # m/s
            v_r = Vr # m/s
            )

    delay = np.array(sim_config.delay_chip*1)
    doppler_specular_point = np.array(eq_doppler_absolute_shift(np.array([0,0,0]), sim_config))
    fd = doppler_specular_point

    spherical.compute_transformation(delay, fd, sim_config)

    '''

    sim_config = simulation_configuration()

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

    # Iso lines plot
    fig_isolines, ax_isolines = plt.subplots(1,figsize=(10, 4))

    contour_delay_chip = ax_isolines.contour(
            x_grid, y_grid, z_grid_delay_chip, 
            np.arange(0, delay_increment_end/delay_chip, 1), 
            cmap='winter', alpha = 0.6
            )
    contour_doppler = ax_isolines.contour(
            x_grid, y_grid, z_grid_doppler_increment, 
            np.arange(doppler_increment_start, doppler_increment_end, 500), 
            cmap='jet', alpha = 0.8
            )

    test_delay = np.array([1*delay_chip])
    test_doppler = np.array([500]) +  doppler_specular_point

    print('Finding intersection for d:{0}, f:{1}'.format(test_delay, test_delay))

    '''
    x_s_1 = planar.x_delay_doppler_1(test_delay, test_doppler, sim_config)
    y_s_1 = planar.y_delay_doppler_1(test_delay, test_doppler, sim_config)
    x_s_2 = planar.x_delay_doppler_2(test_delay, test_doppler, sim_config)
    y_s_2 = planar.y_delay_doppler_2(test_delay, test_doppler, sim_config)
    '''
    r_1, r_2 = spherical.delay_doppler_to_local_surface(test_delay, test_doppler, sim_config)
    x_s_1 = r_1[0]
    y_s_1 = r_1[1]
    x_s_2 = r_2[0]
    y_s_2 = r_2[1]
    ax_isolines.scatter(x_s_1, y_s_1, s=70, marker=(5, 2), zorder=4)
    ax_isolines.scatter(x_s_2, y_s_2, s=70, marker=(5, 2), zorder=4)

    fig_isolines.colorbar(contour_delay_chip, label='C/A chips')
    fig_isolines.colorbar(contour_doppler, label='Hz')
    ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/1000))
    ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000))
    ax_isolines.xaxis.set_major_formatter(ticks_x)
    ax_isolines.yaxis.set_major_formatter(ticks_y)
    plt.xlabel('[km]')
    plt.ylabel('[km]')


    plt.show()

if __name__ == '__main__':
    main()


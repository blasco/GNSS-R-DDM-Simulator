#!/usr/bin/env python

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import gnssr.simulator.rcs.target_rcs as target_rcs
import gnssr.simulator.rcs.sea_rcs as sea_rcs
from gnssr.simulator.isolines import *
from gnssr.simulator.simulation_configuration import *
from gnssr.simulator.ddm import *

from gnssr.utils import *
from gnssr.simulator.jacobian.planar import *

import concurrent.futures
import pickle

def main():

    # Surface mesh
    x_0 = -8e3 
    x_1 = 8e3 # meters
    n_x = 20

    y_0 = -8e3
    y_1 = 8e3 # meters
    n_y = 20

    x_grid, y_grid = np.meshgrid(
       np.linspace(x_0, x_1, n_x), 
       np.linspace(y_0, y_1, n_y)
       )

    z_grid = np.zeros(x_grid.shape)

    parallel_grid = []
    for x_pos, col in enumerate(x_grid):
        for y_pos, _ in enumerate(col):
            parallel_grid.append([
                x_pos,y_pos,
                x_grid[x_pos][y_pos],
                y_grid[x_pos][y_pos] 
                ])

    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
        for results in executor.map(compute_snr, parallel_grid, chunksize=5):
            x_pos, y_pos, z_val = results
            z_grid[x_pos][y_pos] = z_val

    import pdb; pdb.set_trace() # break
    with open('results_snr_map.pkl', 'wb') as f: 
        pickle.dump([x_grid, y_grid, z_grid], f)

    plt.show()

def compute_snr(parallel_input):
    print("parallel_input: {}".format(parallel_input))
    x_pos, y_pos, x_val, y_val = parallel_input
    sim_config = simulation_configuration()
    sim_config.set_scenario_local_ref(
            h_t = 13.82e6, # m
            h_r = 20e3, # meters
            elevation = 70.0*np.pi/180,
            v_t = np.array([-2684.911, 1183.799, -671.829]), # m/s
            v_r = np.array([20, 20, 20]) # m/s
            )

    #sim_config.jacobian_type = 'spherical'
    sim_config.convolve_type = 'fft'
    sim_config.receiver_antenna_gain = lambda p1,p2: 1
    sim_config.rcs = lambda p1,p2: target_rcs.radar_cross_section(p1, 0, p2)
    sim_config.u_10 = 10.00

    #sim_config.delay_chip = 1/gps_ca_chips_per_second # seconds
    delay_chip = sim_config.delay_chip

    number_of_delay_pixels = 128 - 50
    number_of_doppler_pixels = 20 + 50
    sim_config.doppler_increment_start = -70
    sim_config.doppler_increment_end = 70
    sim_config.doppler_resolution = (sim_config.doppler_increment_end - sim_config.doppler_increment_start)/number_of_doppler_pixels/8
    sim_config.delay_increment_start = -1*delay_chip
    sim_config.delay_increment_end = 10*delay_chip
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

    doppler_specular_point = eq_doppler_absolute_shift(np.array([0,0,0]), sim_config)

    sim_config.target_x = x_val
    sim_config.target_y = y_val

    # DDM
    sim_config.rcs = lambda p1,p2: target_rcs.radar_cross_section(p1, 0, p2)
    ddm_sim = np.copy(simulate_ddm(sim_config))

    sim_config.rcs = sea_rcs.radar_cross_section
    ddm_sim_1 = np.copy(simulate_ddm(sim_config))

    ddm_sim_res = rescale(ddm_sim, number_of_doppler_pixels, number_of_delay_pixels)
    ddm_sim_1_res = rescale(ddm_sim_1, number_of_doppler_pixels, number_of_delay_pixels)
    ddm_diff_res = np.abs(ddm_sim_res - ddm_sim_1_res)

    # DDM Noise
    T_noise_receiver = 232.09 + 246.85
    k_b = 1.38e-23 # J/K
    y_noise = 2/sim_config.coherent_integration_time*k_b*T_noise_receiver

    ddm_snr = np.copy(10*np.log10(np.abs(ddm_diff_res)/y_noise))
    np.place(ddm_snr, ddm_snr == np.nan, -1000)
    np.place(ddm_snr, ddm_snr == np.inf, -1000)
    np.place(ddm_snr, ddm_snr == -np.inf, -1000)

    max_indices = np.where(ddm_snr == np.amax(ddm_snr))

    snr_max = np.amax(ddm_snr)

    delay_target = delay_increment_start/delay_chip + (delay_increment_end/delay_chip - delay_increment_start/delay_chip)/number_of_delay_pixels*max_indices[1]
    doppler_target = doppler_increment_start + (doppler_increment_end - doppler_increment_start)/number_of_doppler_pixels*max_indices[0] + doppler_specular_point

    x_snr_1 = x_delay_doppler_1(delay_target, doppler_target, sim_config)
    y_snr_1 = y_delay_doppler_1(delay_target, doppler_target, sim_config)
    x_snr_2 = x_delay_doppler_2(delay_target, doppler_target, sim_config)
    y_snr_2 = y_delay_doppler_2(delay_target, doppler_target, sim_config)

    if (x_snr_1.size > 1 or x_snr_2.size > 1):
        snr_max = np.nan

    return x_pos, y_pos, snr_max

if __name__ == '__main__':
    main()

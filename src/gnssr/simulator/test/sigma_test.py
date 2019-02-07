#!/usr/bin/env python

import cv2
import numpy as np
import matplotlib.pyplot as plt
from gnssr.simulator.sigma import *
from gnssr.simulator.simulation_configuration import *
from gnssr.simulator.isolines import *
from gnssr.utils import gps_ca_chips_per_second

def main():

    sim_config = simulation_configuration()
    sim_config.doppler_resolution = 50 # Hz
    sim_config.delay_resolution = 0.2/gps_ca_chips_per_second # seconds

    delay_increment_start = sim_config.delay_increment_start 
    delay_increment_end = sim_config.delay_increment_end 
    delay_resolution = sim_config.delay_resolution

    doppler_increment_start = sim_config.doppler_increment_start
    doppler_increment_end = sim_config.doppler_increment_end
    doppler_resolution = sim_config.doppler_resolution

    delay_chip = sim_config.delay_chip

    # Delay Doppler grid
    delay_increment_values = list(np.arange(
        delay_increment_start, 
        delay_increment_end, 
        delay_resolution
        ))
    doppler_increment_values = list(np.arange(
        doppler_increment_start, 
        doppler_increment_end, 
        doppler_resolution
        ))

    doppler_specular_point = eq_doppler_absolute_shift(np.array([0,0,0]), sim_config)
    doppler_absolute_values = doppler_increment_values + doppler_specular_point

    delay_grid, doppler_grid = np.meshgrid(delay_increment_values, doppler_absolute_values)

    # Sigma
    sigma_matrix = sigma(delay_grid, doppler_grid, sim_config)

    # Plot
    fig_sigma, ax_sigma = plt.subplots(1,figsize=(10, 4))
    ax_sigma.set_title('Sigma')
    im = ax_sigma.imshow(sigma_matrix, cmap='viridis', 
            extent=(delay_increment_start/delay_chip, delay_increment_end/delay_chip, doppler_increment_start, doppler_increment_end),
            aspect="auto"
            )

    plt.show()

if __name__ == '__main__':
    main()

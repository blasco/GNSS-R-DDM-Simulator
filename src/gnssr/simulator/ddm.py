#!/usr/bin/env python

import numpy as np

from gnssr.simulator.waf import *
from gnssr.simulator.sigma import *
from gnssr.simulator.isolines import *

from scipy import signal

def simulate_ddm(sim_config):
    """
    Simulates the DDM based on the Z-V model for the system defined in the 
    problem_definition file .

        V. U. Zavorotny and A. G. Voronovich, “Scattering of GPS signals from the 
        ocean with wind remote sensing application,” IEEE Transactions on Geoscience and 
        Remote Sensing, vol. 38, no. 2, pp. 951–964, Mar. 2000.  

    Args:
        sim_config: Instance of simulation_configuration class.

    Returns:
        numpy.ndarray with size:
            len(list(np.arange(delay_increment_start, delay_increment_end, delay_resolution))),
            len(list(np.arange(doppler_increment_start, doppler_increment_end, doppler_resolution)))
    """
    delay_chip = sim_config.delay_chip

    delay_increment_start = sim_config.delay_increment_start 
    delay_increment_end = sim_config.delay_increment_end 
    delay_resolution = sim_config.delay_resolution

    doppler_increment_start = sim_config.doppler_increment_start
    doppler_increment_end = sim_config.doppler_increment_end
    doppler_resolution = sim_config.doppler_resolution

    delay_increments_sigma = list(np.arange(
        delay_increment_start, 
        delay_increment_end, 
        delay_resolution
        ))
    delay_increments_waf = list(np.arange(
        -delay_increment_end, 
        delay_increment_end, 
        delay_resolution
        ))
    doppler_increments = list(np.arange(
        doppler_increment_start, 
        doppler_increment_end, 
        doppler_resolution
        ))

    doppler_specular_point =  eq_doppler_absolute_shift(np.array([0,0,0]), sim_config)

    delay_grid_sigma, doppler_grid_sigma = np.meshgrid(delay_increments_sigma, doppler_increments + doppler_specular_point)
    delay_grid_waf, doppler_grid_waf = np.meshgrid(delay_increments_waf, doppler_increments)

    sigma_matrix = sigma(delay_grid_sigma, doppler_grid_sigma, sim_config)
    waf_matrix = woodward_ambiguity_function(delay_grid_waf, doppler_grid_waf, sim_config)**2
    p = signal.convolve2d(sigma_matrix, waf_matrix, mode='same')
    return p

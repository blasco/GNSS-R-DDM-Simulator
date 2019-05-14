#!/usr/bin/env python

import numpy as np

import matplotlib.pyplot as plt

from gnssr.utils import *
from gnssr.simulator.waf import *
from gnssr.simulator.simulation_configuration import *

def add_thermal_noise(n_rows, n_cols, n_incoherent, noise_temperature, ddm_power, sim_config):
    sim_config = simulation_configuration()

    delay_start = sim_config.delay_increment_start
    delay_end = sim_config.delay_increment_end
    delay_resolution = (delay_end - delay_start)/n_cols
    k_b = 1.38e-23 # J/K
    #T_r = 225.7 # K
    T_r = noise_temperature

    delay_values = list(np.arange(
                            delay_start, 
                            delay_end,
                            delay_resolution
                    ))

    covariance = np.zeros((len(delay_values),len(delay_values)))
    for i, col in enumerate(covariance):
        for j, val in enumerate(col):
            covariance[i,j] = 2/1e-3*2*k_b*T_r*waf_delay(delay_values[j] - delay_values[i], sim_config)

    fig_covar, ax_covar = plt.subplots(1,figsize=(10, 4))
    contour_ddm = ax_covar.imshow(covariance, cmap='jet', 
            aspect="auto"
            )

    ddm_noise = np.zeros((n_rows,len(delay_values)))
    for i in range(n_incoherent):
        print("i: {0}".format(i))
        ddm_noise_i = np.power(np.random.multivariate_normal(np.zeros(len(delay_values)),covariance, (n_rows)),2)
        ddm_noise += ddm_noise_i/np.sqrt(n_incoherent) + ddm_power

    ddm_noise /= n_incoherent

    fig_ddm_noise, ax_ddm_noise = plt.subplots(1,figsize=(10, 4))
    contour_ddm_noise = ax_ddm_noise.imshow(ddm_noise, cmap='jet', 
            aspect="auto"
            )
    cbar = fig_ddm_noise.colorbar(contour_ddm_noise, label='Power')

    return ddm_noise


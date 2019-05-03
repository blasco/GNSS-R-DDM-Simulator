#!/usr/bin/env python

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
from gnssr.utils import *

from gnssr.simulator.waf import *
from gnssr.simulator.simulation_configuration import *

def main():
    P_t = 26.61 # W
    G_t = 20 
    G_r = 25
    wavelength = 0.19  # m
    h_t = 13.82e6 # m
    h_r = 20e3  # m
    k_b = 1.38e-23 # J/K
    T_r = 225.7 # K
    T_i = 1e-3 # s
    sigma_target = 12

    # From Gleason page 76
    P_n = k_b*T_r/T_i

    # Processed power
    Y_n = 1/T_i*k_b*T_r
    print("Y_n = {0}".format(Y_n))

    std = np.sqrt(k_b*T_r/T_i)
    x = np.random.normal(0,std,1000)
    print("x: {0}".format(np.sqrt(np.mean(np.power(x,2)))))

    sim_config = simulation_configuration()

    delay_chip = 1/1.023*1e-6
    delay_values = list(np.arange(
                            -5*delay_chip, 
                            5*delay_chip,
                            0.1*delay_chip
                    ))

    x = np.zeros((len(delay_values),len(delay_values)))
    for i, col in enumerate(x):
        for j, val in enumerate(col):
            x[i,j] = np.power((2/1e-3*2*k_b*T_r*waf_delay(delay_values[j] - delay_values[i], sim_config)),2)

    fig_covar, ax_covar = plt.subplots(1,figsize=(10, 4))
    contour_ddm = ax_covar.imshow(x, cmap='jet', 
            aspect="auto"
            )

    '''
    ddm_noise = np.zeros((40, len(delay_values)))
    for i in range(40):
        import pdb; pdb.set_trace() # break
        ddm_noise[i,:] = np.random.normal(0,x)
        '''

    fig_ddm_noise, ax_ddm_noise = plt.subplots(1,figsize=(10, 4))
    ddm_noise = np.zeros((50,len(delay_values)))
    n = 1000
    for i in range(n):
        print("i: {0}".format(i))
        ddm_noise_i = np.abs(np.random.multivariate_normal(np.zeros(len(delay_values)),x, (50)))
        ddm_noise += ddm_noise_i

    ddm_noise /= n

    contour_ddm = ax_ddm_noise.imshow(ddm_noise, cmap='jet', 
            aspect="auto"
            )

    P_r = (P_t*G_t*wavelength**2*sigma_target*G_r) / (
            (4*np.pi)**3 * (h_t*h_r)**2 \
    )

    print("P_r = {0}".format(P_r))
    print("P_n = {0}".format(P_n))
    print("SNR = {0}".format(10*np.log10(P_r/P_n)))

    plt.show()

if __name__ == '__main__':
    main()

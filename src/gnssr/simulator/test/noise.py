#!/usr/bin/env python

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
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
    T_r = 2*260 # K
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
    '''
    for i, col in enumerate(x):
        for j, val in enumerate(col):
            x[i,j] = (2/1e-3*2*k_b*T_r*waf_delay(delay_values[j] - delay_values[i], sim_config))
    '''

    for i, col in enumerate(x):
        x[i,:] =  (2/1e-3*2*k_b*T_r*waf_delay(delay_values - delay_values[i], sim_config))

    fig_covar, ax_covar = plt.subplots(1,figsize=(10, 4))
    contour_covar = ax_covar.imshow(x, cmap='jet', 
            extent=[
                delay_values[0]/delay_chip, delay_values[-1]/delay_chip, 
                delay_values[0]/delay_chip, delay_values[-1]/delay_chip
                ],
            aspect="auto"
            )

    #ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(delax_values[pos]/delax_chip))
    #ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(delay_values[pos]/delay_chip))
    #ax_covar.xaxis.set_major_formatter(ticks_x)
    #ax_covar.yaxis.set_major_formatter(ticks_y)
    fig_covar.colorbar(contour_covar, label='[Watts]')
    plt.xlabel('C/A chips')
    plt.ylabel('C/A chips')


    '''
    ddm_noise = np.zeros((40, len(delay_values)))
    for i in range(40):
        import pdb; pdb.set_trace() # break
        ddm_noise[i,:] = np.random.normal(0,x)
        '''

    fig_ddm_noise, ax_ddm_noise = plt.subplots(1,figsize=(10, 4))
    ddm_noise = np.zeros((50,len(delay_values)))
    ddm_noise_integrated = np.copy(ddm_noise)
    n_int = 1
    for i in range(n_int):
        n = 1000
        for i in range(n):
            print("i: {0}".format(i))
            ddm_noise_i = np.power(np.random.multivariate_normal(np.zeros(len(delay_values)),x, (50)),2)
            ddm_noise += ddm_noise_i
        ddm_noise /= n
        ddm_noise -= np.mean(ddm_noise)
        ddm_noise = np.abs(ddm_noise)

        ddm_noise_integrated += ddm_noise
    ddm_noise_integrated /= n_int 

    print("mean cleaned noise expected: {}".format(2/1e-3*2*k_b*T_r/np.sqrt(n)))
    print("mean cleaned noise: {}".format(np.mean(ddm_noise_integrated)))


    contour_ddm = ax_ddm_noise.imshow(ddm_noise_integrated, cmap='jet', 
            aspect="auto"
            )

    fig_ddm_noise_i, ax_ddm_noise_i = plt.subplots(1,figsize=(10, 4))
    ddm_noise_i = np.power(np.random.multivariate_normal(np.zeros(len(delay_values)),x, (50)),2)
    contour_ddm = ax_ddm_noise_i.imshow(ddm_noise_integrated - 
            ddm_noise_i/np.sqrt(n), cmap='jet', 
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

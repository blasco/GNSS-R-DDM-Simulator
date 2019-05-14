#!/usr/bin/env python

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
from gnssr.utils import *

from gnssr.simulator.waf import *
from gnssr.simulator.simulation_configuration import *

from gnssr.simulator.isolines import *

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
    print("Y_n amplitude = {0}".format(np.sqrt(Y_n*12)))

    std = np.sqrt(k_b*T_r/T_i)
    x = np.random.normal(0,std,1000)
    print("x: {0}".format(np.sqrt(np.mean(np.power(x,2)))))

    sim_config = simulation_configuration()

    delay_chip = 1/1.023*1e-6
    delay_values = list(np.arange(
                            -15*delay_chip, 
                            15*delay_chip,
                            (15*delay_chip-(-15*delay_chip))/128
                    ))

    doppler_values = list(np.arange(
                            -5000*delay_chip, 
                            5000*delay_chip,
                            (5000*delay_chip-(-5000*delay_chip))/20
                    ))


    doppler_specular_point = eq_doppler_absolute_shift(np.array([0,0,0]), sim_config)

    '''
    ddm_noise = np.zeros((40, len(delay_values)))
    for i in range(40):
        import pdb; pdb.set_trace() # break
        ddm_noise[i,:] = np.random.normal(0,x)
        '''

    n_row = 20
    ddm_noise = np.zeros((n_row,len(delay_values)))
    n_incoh = 10
    x_0 = np.zeros((len(delay_values),len(delay_values)), dtype=np.complex_)
    for i, col in enumerate(x_0):
        tau_increment = delay_values - delay_values[i]
        a1 = 1/1e-3*1*k_b*T_r*waf_delay(tau_increment, sim_config)
        x_0[i,:] = a1


    for i_incoh in range(n_incoh):
        print("i: {0}".format(i_incoh))
        ddm_noise_i = np.zeros((n_row,len(delay_values)))
        for i_row in range(n_row):
            x = np.zeros((len(delay_values),len(delay_values)), dtype=np.complex_)
            for i, col in enumerate(x):
                tau_increment = delay_values - delay_values[i]
                a1 = x_0[i,:]
                a2 = np.exp(-2j*np.pi*(doppler_values[i_row] + doppler_specular_point)*tau_increment)
                a = np.multiply(a1,a2)
                x[i,:] = a
            '''
            fig_var, ax_var = plt.subplots(1,figsize=(10, 4))
            contour_ddm_noise = ax_var.imshow(np.absolute(x), cmap='jet', 
                    aspect="auto"
                    )
            plt.show()
            '''
            ddm_noise_i[i_row,:] = np.power(np.absolute(np.random.multivariate_normal(np.zeros(len(delay_values)),x)),2)
        ddm_noise += ddm_noise_i

    ddm_noise /= n_incoh

    fig_ddm_noise, ax_ddm_noise = plt.subplots(1,figsize=(10, 4))
    contour_ddm_noise = ax_ddm_noise.imshow(ddm_noise, cmap='jet', 
            aspect="auto"
            )
    cbar = fig_ddm_noise.colorbar(contour_ddm_noise, label='Power')

    P_r = (P_t*G_t*wavelength**2*sigma_target*G_r) / (
            (4*np.pi)**3 * (h_t*h_r)**2 \
    )

    print("P_r = {0}".format(P_r))
    print("P_n = {0}".format(P_n))
    print("SNR = {0}".format(10*np.log10(P_r/P_n)))

    plt.show()

if __name__ == '__main__':
    main()

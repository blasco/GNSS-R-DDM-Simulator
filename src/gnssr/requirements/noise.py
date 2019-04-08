#!/usr/bin/env python

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
from gnssr.utils import *

def main():
    P_t = 26.61 # W
    G_t = 20 
    G_r = 25
    wavelength = 0.19  # m
    h_t = 13.82e6 # m
    h_r = 12e3  # m
    k_b = 1.38e-23 # J/K
    T_r = 225.7 # K
    T_i = 1e-2 # s
    sigma_target = 12

    # From Gleason page 76
    P_n = k_b*T_r/T_i

    # Processed power
    Y_n = T_i*k_b*T_r
    print("Y_n = {0}".format(Y_n))

    P_r = (P_t*G_t*wavelength**2*sigma_target*G_r) / (
            (4*np.pi)**3 * (h_t*h_r)**2 \
    )

    print("P_r = {0}".format(P_r))
    print("P_n = {0}".format(P_n))
    print("SNR = {0}".format(10*np.log10(P_r/P_n)))

if __name__ == '__main__':
    main()

#!/usr/bin/env python

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
from gnssr.utils import *

def main():
    wavelength = 0.19
    w = 2*np.pi*1/7.2
    A = 2.5
    t_coh = wavelength/(4*w*A)

    print("t_coh: {0}".format(t_coh))

    chip = 1/gps_ca_chips_per_second/10 # seconds
    light_speed = 299792458.0 # m/s
    h_t = 13.82e6, # m
    h_r = 20e3, # meters
    theta = 40*np.pi/180
    rs = h_r/np.cos(theta)
    ts = h_t/np.cos(theta)
    a = (rs*ts*chip*light_speed/(rs+ts))**(1/2)
    print(a)
    t_coh_1 = 1.22*wavelength*(18000)/(2*a*30)

    print("t_coh_1: {0}".format(t_coh_1))
if __name__ == '__main__':
    main()

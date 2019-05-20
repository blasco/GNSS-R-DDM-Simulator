#!/usr/bin/env python

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import gnssr.simulator.rcs.target_rcs as target_rcs
import gnssr.simulator.rcs.sea_rcs as sea_rcs
from gnssr.simulator.isolines import *
from gnssr.simulator.simulation_configuration import *
from gnssr.simulator.ddm import *

import cv2

def main():

    u = np.arange(0,25,0.001)
    #t_c = 0.19/(4*2*3.14*(0.13*9.8/u)*0.27*u**2/9.8)
    t_c = 0.19/(0.28*np.pi*u)

    plt.plot(u, t_c)
    plt.yscale('log')
    plt.grid(True, which = "both")


    u = 5
    print(0.19/(4*2*3.14*(0.13*9.8/u)*0.27*u**2/9.8) )

    plt.figure()

    phi = np.arange(30,90,0.01)
    t_c_2 = 1.22*0.19*20e3/30/(2/np.sin(phi*np.pi/180)*np.sqrt(20e3*18e6*(1/1023)*1e-3*3e8/(20e3+18e6)))

    plt.plot(phi, t_c_2)
    plt.grid()
 
    plt.show()


if __name__ == '__main__':
    main()

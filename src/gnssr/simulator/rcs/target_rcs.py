#!/usr/bin/env python

import numpy as np
from gnssr.simulator.geometry.geometry import *
import gnssr.simulator.rcs.sea_rcs as sea_rcs

def radar_cross_section(r, t, sim_config):
    """
    Args:
        r (numpy.ndarray with size(3,)): Position on the local coordinates.
        sim_config: 
        sim_config: Instance of simulation_configuration class.
    Returns:
        numpy.ndarray with  size(1,).
    """

    #100x25m 
    v_x = 33 # m/s
    wake_x = 4e3 - v_x*t#m
    wake_y = 0.5e3 #m
    wake_x_size = 150 #m
    wake_y_size = wake_x_size*np.tan(20*np.pi/180)
    wake = np.logical_and.reduce((
            np.abs(r[0]-wake_x) < wake_x_size, 
            wake_y_size/wake_x_size*(r[0]-wake_x) - (r[1]-wake_y) >= 0, 
            wake_y_size/wake_x_size*(r[0]-wake_x) + (r[1]-wake_y) >= 0
            ))
    sea = sea_rcs.radar_cross_section(r, sim_config)
    #target = rcs_sea(r, sim_config, sim_config.u_10*2)
    np.place(sea, wake, 0)
    #np.place(target, np.logical_not(wake), 0)
    return sea #+ target

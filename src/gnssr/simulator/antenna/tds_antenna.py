#!/usr/bin/env python

import numpy as np

from gnssr.tds.antenna.gain import *

sim_antenna = tds_antenna_gain()
def receiver_antenna_gain(r, sim_config):
    """
    Args:
        r (numpy.ndarray with size(3,)): Position on the earth surface in ECEF coordinates.
        sim_config: Instance of simulation_configuration class.
    Returns:
        numpy.ndarray with size(1,).
    """
    r_antenna = r[0:2] - sim_config.r_r[0:2] 
    elevation = np.arctan2(sim_config.h_r,np.sqrt(r_antenna[0]**2+r_antenna[1]**2))*180/np.pi
    azimuth = np.arctan2(-r_antenna[1],-r_antenna[0])*180/np.pi
    return sim_antenna.gain(azimuth, elevation)

def transmitting_antenna_gain(r, sim_config):
    """
    Args:
        r (numpy.ndarray with size(3,)): Position on the earth surface in ECEF coordinates.
        sim_config: Instance of simulation_configuration class.
    Returns:
        numpy.ndarray with size(1,)
    """
    # Coulson, “The GPS at GTO Experiment.” 1996.
    gps_isotropic_antenna_gain = np.power(10, 14.4/10)
    return gps_isotropic_antenna_gain


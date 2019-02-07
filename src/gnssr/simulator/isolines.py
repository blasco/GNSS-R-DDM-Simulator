#!/usr/bin/env python

import numpy as np

from gnssr.utils import *

def eq_doppler_absolute_shift(r, sim_config):
    """
    Doppler shift as a contribution of the relative motion of transmitter and 
    receiver as well as the reflection point. 

    Implements equation 14:
        V. U. Zavorotny and A. G. Voronovich, “Scattering of GPS signals from 
        the ocean with wind remote sensing application,” IEEE Transactions on 
        Geoscience and Remote Sensing, vol. 38, no. 2, pp. 951–964, Mar. 2000.  

    Args:
        r (numpy.ndarray with size(3,)): Position on the earh surface in ECEF coordinates.
        sim_config: Instance of simulation_configuration class.

    Returns:
        numpy.ndarray with size(1,).
    """
    #f_surface = scattering_vector(r)*v_surface(r)/2*pi
    f_surface = 0

    f_carrier = sim_config.f_carrier
    v_t = sim_config.v_t
    v_r = sim_config.v_r
    h_r = sim_config.h_r
    elevation = sim_config.elevation

    # With the very far way transmitter approximation solved in Maple:
    f_D_0 =  f_carrier / light_speed * (-v_t[1] * np.cos(elevation) - v_t[2] * np.sin(elevation) + (v_r[0] * r[0] + v_r[1] * (r[1] + h_r / np.tan(elevation)) - v_r[2] * h_r) * (r[0] ** 2 + (r[1] + h_r / np.tan(elevation)) ** 2 + h_r ** 2) ** (-0.1e1 / 0.2e1))
    return f_D_0 + f_surface

def eq_doppler_absolute_shift_increment(r, sim_config):
    """
    Args:
        r (numpy.ndarray with size(3,)): Position on the earth surface in ECEF coordinates.
        sim_config: Instance of simulation_configuration class.

    Returns:
        numpy.ndarray with size(1,).
    """
    return eq_doppler_absolute_shift(r, sim_config) - eq_doppler_absolute_shift(np.array([0,0,0]))

def eq_delay_incremet(r, sim_config):
    """
    Args:
        r (numpy.ndarray with size(3,)): Position on the earth surface in ECEF coordinates.
        sim_config: Instance of simulation_configuration class.

    Returns:
        numpy.ndarray with  size(1,).
    """
    h_r = sim_config.h_r
    elevation = sim_config.elevation
    return (np.sqrt(r[0] ** 2 + (r[1] + h_r / np.tan(elevation)) ** 2 + h_r ** 2) - h_r / np.sin(elevation) - r[1] * np.cos(elevation)) / light_speed


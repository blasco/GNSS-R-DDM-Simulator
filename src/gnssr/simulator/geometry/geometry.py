#!/usr/bin/env python

import numpy as np

def scattering_vector(r, sim_config):
    """
    Implements equation 7:
        V. U. Zavorotny and A. G. Voronovich, “Scattering of GPS signals from 
        the ocean with wind remote sensing application,” IEEE Transactions on 
        Geoscience and Remote Sensing, vol. 38, no. 2, pp. 951–964, Mar. 2000.  
    Args:
        r (numpy.ndarray with size(3,)): Position on the earth surface in ECEF coordinates.
        sim_config: Instance of simulation_configuration class.
    Returns:
        numpy.ndarray with size(3,).
    """
    scattering_vector = (scattered_direction(r, sim_config) - incident_direction(r, sim_config))
    scattering_vector_norm = np.linalg.norm(scattering_vector)
    scattering_vector[0] /= scattering_vector_norm
    scattering_vector[1] /= scattering_vector_norm
    scattering_vector[2] /= scattering_vector_norm
    return scattering_vector

def scattered_direction(r, sim_config):
    r"""
    The scattered unit vector: 

    .. math::
        \frac{r_r - r}{|r_r - r|}

    Args:
        r (numpy.ndarray with size(3,)): Position on the earth surface in ECEF coordinates.
        sim_config: Instance of simulation_configuration class.
    Returns:
        numpy.ndarray with size(3,).
    """
    r_r = sim_config.r_r

    scattered_direction = (r_r - r)
    scattered_direction_norm = np.linalg.norm(r_r - r)
    scattered_direction[0] /= scattered_direction_norm
    scattered_direction[1] /= scattered_direction_norm
    scattered_direction[2] /= scattered_direction_norm
    return scattered_direction

def incident_direction(r, sim_config):
    r"""
    The incident unit vector: 
    
    .. math::
        \frac{r - r_t}{|r - r_t|}

    Args:
        r (numpy.ndarray with size(3,)): Position on the earth surface in ECEF coordinates.
        sim_config: Instance of simulation_configuration class.
    Returns:
        numpy.ndarray with size(3,).
    """
    r_t = sim_config.r_t

    incident_direction = (r - r_t)
    incident_direction_norm = np.linalg.norm(r - r_t)
    incident_direction[0] /= incident_direction_norm
    incident_direction[1] /= incident_direction_norm
    incident_direction[2] /= incident_direction_norm
    return  incident_direction


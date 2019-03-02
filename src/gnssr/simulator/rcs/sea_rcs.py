#!/usr/bin/env python

import numpy as np
from gnssr.simulator.geometry.geometry import *

def radar_cross_section(r, sim_config):
    """
    Args:
        r (numpy.ndarray with size(3,)): Position on the local coordinate.
        sim_config: 
        sim_config: Instance of simulation_configuration class.
    Returns:
        numpy.ndarray with  size(1,).
    """
    return rcs_sea(r, sim_config)

def rcs_sea(r, sim_config):
    '''
    Radar Cross Section of the sea surface.
    Implements equation 34:
        [1]V. U. Zavorotny and A. G. Voronovich, “Scattering of GPS signals from 
        the ocean with wind remote sensing application,” IEEE Transactions on 
        Geoscience and Remote Sensing, vol. 38, no. 2, pp. 951–964, Mar. 2000.  
    '''
    fresnel_coefficient = sim_config.fresnel_coefficient

    q = scattering_vector(r, sim_config)
    q_norm = np.linalg.norm(q)
    q_tangent = q[0:2]
    q_z = [2]
    ocean_surface_slope = -q_tangent/q_z

    return np.pi*(fresnel_coefficient**2)*((q_norm/q_z)**4) * \
            slope_probability_density_function(
                    ocean_surface_slope, 
                    sim_config.u_10, 
                    sim_config.phi_0
                    )

def slope_probability_density_function(x, u_10, phi_0):
    '''
    Implements equation 4:
        [1]Q. Yan and W. Huang, “GNSS-R Delay-Doppler Map Simulation Based on the 
        2004 Sumatra-Andaman Tsunami Event,” Journal of Sensors, vol. 2016, pp. 
        1–14, 2016.  
    '''
    local_to_wind_reference_frame = np.array([
        [np.cos(phi_0), np.sin(phi_0)],
        [-np.sin(phi_0), np.cos(phi_0)]
    ])

    s =  local_to_wind_reference_frame.dot(x)

    rms_c = np.sqrt(variance_crosswind(u_10))
    rms_u = np.sqrt(variance_upwind(u_10))

    covariance = np.array([
        [rms_u**2, 0],
        [0, rms_c**2]
    ])

    hermite_coeficients = np.zeros((5,5))
    hermite_coeficients[0,0] = 1
    hermite_coeficients[2,2] = 0.45*(0.01 - 0.0086*f_u_10(u_10))
    hermite_coeficients[3,0] = 0.45*(0.04 -0.033*f_u_10(u_10))
    hermite_coeficients[0,4] = 0.45*(0.4)
    hermite_coeficients[2,2] = 0.45*(0.12)
    hermite_coeficients[4,0] = 0.45*(0.23)

    result = 1/(2*np.pi*rms_u*rms_c) * \
            np.exp(
                -1/2*((s[0]/rms_u)**2+(s[1]/rms_c)**2) \
            ) * \
            np.polynomial.hermite.hermval2d(s[0]/rms_u, s[1]/rms_c, hermite_coeficients)

    np.place(result, result < 0, 0)
    return result

def f_u_10(u_10):
    return np.piecewise(u_10, 
        [
            u_10 <= 3.49,
            np.logical_and(u_10 > 3.49, u_10 <= 46),
            u_10 > 46
            
        ],
        [
            lambda x: x,
            lambda x: 6*np.log(x) - 4,
            lambda x: 0.411*x
        ])

def variance_upwind(u_10):
    ''' 
    Based on the 'clean surface mean square slope model' of Cox and Munk
    Implements Equation 4
        Q. Yan and W. Huang, “GNSS-R Delay-Doppler Map Simulation Based on the 
        2004 Sumatra-Andaman Tsunami Event,” Journal of Sensors, vol. 2016, pp. 
        1–14, 2016.  

    Args: 
        u_10:   Wind speed at 10 meters above sea surface

    Returns:
        upwind variance
    '''
    return (0.45*(3.16e-3*f_u_10(u_10)))

def variance_crosswind(u_10):
    ''' 
    Based on the 'clean surface mean square slope model' of Cox and Mux
    Implements Equation 4
        Q. Yan and W. Huang, “GNSS-R Delay-Doppler Map Simulation Based on the 
        2004 Sumatra-Andaman Tsunami Event,” Journal of Sensors, vol. 2016, pp. 
        1–14, 2016.  

    Args:
        u_10:   Wind speed at 10 meters above sea surface

    Returns:    
        crosswind variance
    '''
    return (0.45*(0.003 + 1.92e-3*f_u_10(u_10)))

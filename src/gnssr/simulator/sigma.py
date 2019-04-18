#!/usr/bin/env python

import numpy as np

import gnssr.simulator.jacobian.planar as planar
import gnssr.simulator.jacobian.spherical as spherical
from gnssr.utils import *

def sigma(delay, doppler, sim_config):
    """
    Accounts for the surface geometry, the antenna patterns, and the bistatic 
    radar cross section of each cell. Assigns a weighting factor to each 
    delay-Doppler cell of the scene. 

    Implements equation 10:
        J. F. Marchan-Hernandez, A. Camps, N. Rodriguez-Alvarez, E. Valencia, X.  
        Bosch-Lluis, and I. Ramos-Perez, “An Efficient Algorithm to the 
        Simulation of Delay–Doppler Maps of Reflected Global Navigation 
        Satellite System Signals,” IEEE Transactions on Geoscience and Remote 
        Sensing, vol. 47, no.  8, pp. 2733–2740, Aug. 2009.  

    Args:
        delay (numpy.ndarray with size(1,)): Delay increment.
        doppler (numpy.ndarray with size(1,)): Doppler increment
        sim_config: Instance of simulation_configuration class.

    Returns:
        numpy.ndarray with  size(1,).
    """
    
    if (sim_config.jacobian_type == 'planar'):
        r_t = sim_config.r_t
        r_r = sim_config.r_r
        v_t = sim_config.v_t
        v_r = sim_config.v_r
        elevation = sim_config.elevation
        f_carrier = sim_config.f_carrier
        h_r = sim_config.h_r
        transmitting_power = sim_config.transmitting_power
        coherent_integration_time = sim_config.coherent_integration_time

        x_1 = planar.x_delay_doppler_1(delay, doppler, sim_config).real
        y_1 = planar.y_delay_doppler_1(delay, doppler, sim_config).real
        r_1 = np.array([x_1,y_1,0])

        x_2 = planar.x_delay_doppler_2(delay, doppler, sim_config).real
        y_2 = planar.y_delay_doppler_2(delay, doppler, sim_config).real
        r_2 = np.array([x_2,y_2,0])

        radar_cross_section = sim_config.rcs
        wavelength = light_speed/sim_config.f_carrier

        p = transmitting_power*wavelength**2/(4*np.pi)**3 * ( \
                    radar_cross_section(r_1, sim_config)/( \
                        np.linalg.norm(r_1-r_t)**2* \
                        np.linalg.norm(r_r-r_1)**2 \
                    ) * \
                    planar.delay_doppler_jacobian_1(delay, doppler, sim_config) * \
                    sim_config.receiver_antenna_gain(r_1, sim_config) * \
                    sim_config.transmitting_antenna_gain(r_1, sim_config) \
                    + 
                    radar_cross_section(r_2, sim_config)/( \
                        np.linalg.norm(r_2-r_t)**2* \
                        np.linalg.norm(r_r-r_2)**2 \
                    ) * \
                    planar.delay_doppler_jacobian_2(delay, doppler, sim_config) * \
                    sim_config.receiver_antenna_gain(r_2, sim_config) * \
                    sim_config.transmitting_antenna_gain(r_2, sim_config) \
                )*sim_config.doppler_resolution*sim_config.delay_resolution

    if (sim_config.jacobian_type == 'spherical'):
        r_t = sim_config.r_t
        r_r = sim_config.r_r
        v_t = sim_config.v_t
        v_r = sim_config.v_r
        elevation = sim_config.elevation
        f_carrier = sim_config.f_carrier
        h_r = sim_config.h_r
        transmitting_power = sim_config.transmitting_power
        coherent_integration_time = sim_config.coherent_integration_time

        r_1, r_2 = spherical.delay_doppler_to_local_surface(delay, doppler, sim_config)

        r_1 = np.array([r_1[0].real,r_1[1].real,r_1[2].real,0])[0:3]
        r_2 = np.array([r_2[0].real,r_2[1].real,r_2[2].real,0])[0:3]

        j1, j2 = spherical.delay_doppler_jacobian(delay, doppler, sim_config)

        radar_cross_section = sim_config.rcs
        wavelength = light_speed/sim_config.f_carrier

        p = transmitting_power*wavelength**2/(4*np.pi)**3 * ( \
                    radar_cross_section(r_1, sim_config)/( \
                        np.linalg.norm(r_1-r_t)**2* \
                        np.linalg.norm(r_r-r_1)**2 \
                    ) * \
                    j1 * \
                    sim_config.receiver_antenna_gain(r_1, sim_config) * \
                    sim_config.transmitting_antenna_gain(r_1, sim_config) \
                    + 
                    radar_cross_section(r_2, sim_config)/( \
                        np.linalg.norm(r_2-r_t)**2* \
                        np.linalg.norm(r_r-r_2)**2 \
                    ) * \
                    j2 * \
                    sim_config.receiver_antenna_gain(r_2, sim_config) * \
                    sim_config.transmitting_antenna_gain(r_2, sim_config) \
                )*sim_config.doppler_resolution*sim_config.delay_resolution

    # Plot
    #fig_sigma, ax_sigma = plt.subplots(1,figsize=(10, 4))
    #ax_sigma.set_title('Sigma')
    #im = ax_sigma.imshow(p, cmap='viridis', 
    #        aspect="auto"
    #        )

    #noise = 0.01*np.max(p)*np.random.normal(0.02,0.05, p.shape)
    return p
     

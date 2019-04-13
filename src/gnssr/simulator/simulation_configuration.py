#!/usr/bin/env python

import numpy as np
from gnssr.utils import *
from gnssr.simulator.rcs.sea_rcs import *

class simulation_configuration:

    def __init__(self):
        f_0 = 10.23e6 # Hz
        self.f_carrier = 154*f_0
        """
        GPS L1 center frequency is defined in relation to a reference frequency 
        f_0 = 10.23e6, so that f_carrier = 154*f_0 = 1575.42e6 # Hz 
        Explained in section 'DESCRIPTION OF THE EMITTED GPS SIGNAL' in Zarotny 
        and Voronovich 2000.
        """

        #self.fresnel_coefficient = 0.65
        self.fresnel_coefficient = 1

        self.transmitting_power = np.power(10, 14.25/10)
        """From Coulson, “The GPS at GTO Experiment.” 1996."""

        self.coherent_integration_time = 1e-3 # seconds
        self.delay_chip =  1/gps_ca_chips_per_second # seconds

        self.u_10 = 3.0 # m/s 
        """Wind speed at 10 meters above sea surface [m/s]."""
        self.phi_0 = (180+125)*np.pi/180 # rad 
        """Angle between the up-down wind direction and the x-axis."""

        self.delay_increment_start = -15*self.delay_chip # seconds
        self.delay_increment_end = 15*self.delay_chip # seconds
        self.delay_resolution = 0.2*self.delay_chip # seconds

        self.doppler_increment_start = -4500 # Hz
        self.doppler_increment_end = 4500 # Hz
        self.doppler_resolution = 50.0 # Hz

        self.set_scenario_local_ref()

        self.rcs = radar_cross_section

    def set_scenario_local_ref(self,
            h_t = 13.82e6, # m
            h_r = 620e3,  # m
            elevation = 80*np.pi/180, # radians
            v_t = np.array([-2684.911, 1183.799, -671.829]), # m/s
            v_r = np.array([-3605, 3843.852, -850]) # m/s
            ):
        """
        Local Coordinate Frame as defined in Figure 2:
             J. F. Marchan-Hernandez, A. Camps, N. Rodriguez-Alvarez, E. Valencia, X. 
             Bosch-Lluis, and I. Ramos-Perez, “An Efficient Algorithm to the Simulation of 
             Delay–Doppler Maps of Reflected Global Navigation Satellite System Signals,” 
             IEEE Transactions on Geoscience and Remote Sensing, vol. 47, no. 8, pp. 
             2733–2740, Aug. 2009.  
        """
        self.h_t = h_t
        self.h_r = h_r
        self.elevation = elevation
        self.v_t = v_t
        self.v_r = v_r
        self.r_t = np.array([0,h_t/np.tan(elevation),h_t])
        self.r_r = np.array([0,-h_r/np.tan(elevation),h_r])

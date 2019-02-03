#!/usr/bin/env python

import numpy as np

h_t = 16836198.649 # meters
#h_r = 12487.690 # meters
h_r = 622487.690 # meters
elevation = 80*np.pi/180 # rad

# Coordinate Frame as defined in Figure 2
#      J. F. Marchan-Hernandez, A. Camps, N. Rodriguez-Alvarez, E. Valencia, X. 
#      Bosch-Lluis, and I. Ramos-Perez, “An Efficient Algorithm to the Simulation of 
#      Delay–Doppler Maps of Reflected Global Navigation Satellite System Signals,” 
#      IEEE Transactions on Geoscience and Remote Sensing, vol. 47, no. 8, pp. 
#      2733–2740, Aug. 2009.  
r_t = np.array([0,h_t/np.tan(elevation),h_t])
r_r = np.array([0,-h_r/np.tan(elevation),h_r])

# Velocity
v_t = np.array([-2684.911, 1183.799, -671.829]) # m/s
#v_r = np.array([-25,32.852, -3]) # m/s
v_r = np.array([-2605,3243.852, -300]) # m/s

light_speed = 299792458.0 # m/s

# GPS L1 center frequency is defined in relation to a reference frequency 
# f_0 = 10.23e6, so that f_carrier = 154*f_0 = 1575.42e6 # Hz 
# Explained in section 'DESCRIPTION OF THE EMITTED GPS SIGNAL' in Zarotny 
# and Voronovich 2000
f_0 = 10.23e6 # Hz
f_carrier = 154*f_0

fresnel_coefficient = 1 # TODO

# From: 
# Coulson, “The GPS at GTO Experiment.” 1996.
transmitting_power = np.exp(14.25/10)
gps_isotropic_antenna_gain = np.exp(14.4/10)

integration_time = 1e-3 # seconds
delay_chip =  1/1.023e6 # seconds
# Wind speed at 10 meters above sea surface
u_10 = 6.0 # m/s 
# Angle between the up-down wind direction and the x-axis
#phi_0 = (135)*np.pi/180 # rad 
phi_0 = (135)*np.pi/180 # rad 

delay_increment_start = -15*delay_chip # seconds
delay_increment_end = 30*delay_chip # seconds
delay_resolution = 0.2*delay_chip # seconds

doppler_increment_start = -8000 # Hz
doppler_increment_end = 8000 # Hz
#doppler_start = -2500 # Hz
#doppler_end = -1200 # Hz
doppler_resolution = 50.0 # Hz

#!/usr/bin/env python

import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from radar_cross_section_sea import *

h_t = 20000e3 # meters
h_r = 500e3 # meters
elevation = 60*np.pi/180 # rad

# Coordinate Frame as defined in Figure 2
#      J. F. Marchan-Hernandez, A. Camps, N. Rodriguez-Alvarez, E. Valencia, X. 
#      Bosch-Lluis, and I. Ramos-Perez, “An Efficient Algorithm to the Simulation of 
#      Delay–Doppler Maps of Reflected Global Navigation Satellite System Signals,” 
#      IEEE Transactions on Geoscience and Remote Sensing, vol. 47, no. 8, pp. 
#      2733–2740, Aug. 2009.  
r_t = np.array([0,h_t/np.tan(elevation),h_t])
r_r = np.array([0,-h_r/np.tan(elevation),h_r])

# Velocity
v_t = np.array([2121, 2121, 5]) # m/s
v_r = np.array([2210, 7299, 199]) # m/s

light_speed = 299792458 # m/s

# GPS L1 center frequency is defined in relation to a reference frequency 
# f_0 = 10.23e6, so that f_carrier = 154*f_0 = 1575.42e6 # Hz 
# Explained in section 'DESCRIPTION OF THE EMITTED GPS SIGNAL' in Zarotny 
# and Voronovich 2000
f_0 = 10.23e6 # Hz
f_carrier = 154*f_0;

# -----------------------
# Bistatic radar equation
# ----------------------

def woodward_ambiguity_function(delay_increment, frequency_increment): 
    ''' 
    The ambiguity function can be approximated by the product of a function 
    dependent on the delay and a function dependent on the frequency. 

    Implements Equation 3
        J. F. Marchan-Hernandez, A. Camps, N. Rodriguez-Alvarez, E. Valencia, X. 
        Bosch-Lluis, and I. Ramos-Perez, “An Efficient Algorithm to the Simulation 
        of Delay–Doppler Maps of Reflected Global Navigation Satellite System 
        Signals,” IEEE Transactions on Geoscience and Remote Sensing, vol. 47, no. 
        8, pp. 2733–2740, Aug. 2009.  
    ''' 

    return list(waf_delay(delay_increment))*waf_frequency(frequency_increment)

def waf_delay(delays):
    return np.where(np.abs(delays) <= delay_chip*(1+delay_chip/integration_time),
                    1 - np.abs(delays)/delay_chip,
                    -delay_chip/integration_time)

def waf_frequency(frequency_increment):
    return np.sin(np.pi*frequency_increment*integration_time) / \
                 (np.pi*frequency_increment*integration_time) * \
                 np.exp(-np.pi*1j*frequency_increment*integration_time)

integration_time = 1e-3 # seconds
delay_chip =  1/1.023e6 # seconds

def doppler_shift(r):
    ''' 
    Doppler shift as a contribution of the relative motion of transmitter and 
    receiver as well as the reflection point. 

    Implements Equation 14
        V. U. Zavorotny and A. G. Voronovich, “Scattering of GPS signals from 
        the ocean with wind remote sensing application,” IEEE Transactions on 
        Geoscience and Remote Sensing, vol. 38, no. 2, pp. 951–964, Mar. 2000.  
    '''
    wavelength = light_speed/f_carrier
    f_D_0 = (1/wavelength)*(
                np.inner(v_t, incident_vector(r)) \
               -np.inner(v_r, reflection_vector(r))
            )
    #f_surface = scattering_vector(r)*v_surface(r)/2*pi
    f_surface = 0
    return f_D_0 + f_surface

def doppler_increment(r):
    return doppler_shift(r) - doppler_shift(np.array([0,0,0]))

def scattering_vector(r):
    '''
    Implements Equation 7
        V. U. Zavorotny and A. G. Voronovich, “Scattering of GPS signals from 
        the ocean with wind remote sensing application,” IEEE Transactions on 
        Geoscience and Remote Sensing, vol. 38, no. 2, pp. 951–964, Mar. 2000.  
    '''
    K = 2*np.pi*f_carrier/light_speed
    return K*(reflection_vector(r,t) - incident_vector(r,t))

def reflection_vector(r):
    reflection_vector = (r_r - r)
    reflection_vector_norm = np.linalg.norm(r_r - r)
    reflection_vector[0] /= reflection_vector_norm
    reflection_vector[1] /= reflection_vector_norm
    reflection_vector[2] /= reflection_vector_norm
    return reflection_vector

#def incident_vector(r):
#    return (r - r_t)/np.linalg.norm(r - r_t)

def incident_vector(r):
    incident_vector = (r - r_t)
    incident_vector_norm = np.linalg.norm(r - r_t)
    incident_vector[0] /= incident_vector_norm
    incident_vector[1] /= incident_vector_norm
    incident_vector[2] /= incident_vector_norm
    return  incident_vector

def time_delay(r):
    path_r = np.linalg.norm(r-r_t) + np.linalg.norm(r_r-r)
    path_specular = np.linalg.norm(r_t) + np.linalg.norm(r_r)
    return (1/light_speed)*(path_r - path_specular)

def radar_cross_section(r):
    return rcs_sea(r)

# -------------------------------------
# Sea Surface Radar Cross Section Model
# -------------------------------------

def rcs_sea(r):
    '''
    Radar Cross Section of the sea surface.
    Implements Equation 34
        [1]V. U. Zavorotny and A. G. Voronovich, “Scattering of GPS signals from 
        the ocean with wind remote sensing application,” IEEE Transactions on 
        Geoscience and Remote Sensing, vol. 38, no. 2, pp. 951–964, Mar. 2000.  
    '''
    fresnel_coefficient = 1 # TODO

    q = scattering_vector(r)
    q_norm = np.linalg.norm(scattering_vector(r))
    q_tangent = np.array([q[0],q[1])
    q_z = [2]
    ocean_surface_slope = -q_tangent/q_z

    return np.pi*(fresnel_coefficient**2)*((q_norm/q_z)**4)* \
            slope_probability_density_function(ocean_surface_slope)

def slope_probability_density_function(x):
    '''
    Implements Equation 4
        [1]Q. Yan and W. Huang, “GNSS-R Delay-Doppler Map Simulation Based on the 
        2004 Sumatra-Andaman Tsunami Event,” Journal of Sensors, vol. 2016, pp. 
        1–14, 2016.  
    '''
    # phi_0 is the angle between the up-down wind direction and the x-axis
    phi_0 = 0
    wind_rotation = np.array([
        [np.cos(phi_0), -np.sin(phi_0)],
        [np.sin(phi_0),  np.cos(phi_0)]
        ])
    covariance = np.array([
        [variance_upwind(wind_speed_10m_above_sea), 0],
        [0, variance_crosswind(wind_speed_10m_above_sea)]
        ])
    W = (wind_rotation.dot(covariance)).dot(np.transpose(wind_rotation)

    return 1/(2*np.pi*(numpy.linalg.det(W)**(1/2))) \
            *np.exp( \
                -1/2*(np.transpose(x).dot(np.linalg.inv(W))).dot(x) \
            )

# --------------------

# Plotting Area

x_0 =  -100e3 # meters
x_1 =  100e3 # meters
n_x = 500

y_0 =  -100e3 # meters
y_1 =  100e3 # meters
n_y = 500

x_grid, y_grid = np.meshgrid(
   np.linspace(x_0, x_1, n_x), 
   np.linspace(y_0, y_1, n_y)
   )

r = [x_grid, y_grid, 0]
z_grid_delay = time_delay(r)/delay_chip
z_grid_doppler = doppler_increment(r)

delay_start = 0 # C/A chips
delay_increment = 0.5 # C/A chips
delay_end = 15 # C/A chips
iso_delay_values = list(np.arange(delay_start, delay_end, delay_increment))

doppler_start = -3000 # Hz
doppler_increment = 500 # Hz
doppler_end = 3000 # Hz
iso_doppler_values = list(np.arange(doppler_start, doppler_end, doppler_increment))

fig_lines, ax_lines = plt.subplots(1,figsize=(10, 4))
contour_delay = ax_lines.contour(x_grid, y_grid, z_grid_delay, iso_delay_values, cmap='winter')
fig_lines.colorbar(contour_delay, label='C/A chips', )

contour_doppler = ax_lines.contour(x_grid, y_grid, z_grid_doppler, iso_doppler_values, cmap='winter')
fig_lines.colorbar(contour_doppler, label='Hz', )

ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/1000))
ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000))
ax_lines.xaxis.set_major_formatter(ticks_x)
ax_lines.yaxis.set_major_formatter(ticks_y)
plt.xlabel('[km]')
plt.ylabel('[km]')

plt.show()

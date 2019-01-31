#!/usr/bin/env python

import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from gnssr.tds.antenna_gain.antenna_gain import *
from scipy import signal
from delay_doppler_jacobian import *
from problem_definition import *

# -----------------------
# Bistatic radar equation
# ----------------------

def waf_squared(delay, f_doppler): 
    ''' 
    The Woodward Ambiguity Function (waf) squared can be approximated by the 
    product of a function dependent on the delay and a function dependent on the 
    frequency. 

    Implements Equation 3.
        J. F. Marchan-Hernandez, A. Camps, N. Rodriguez-Alvarez, E. Valencia, X. 
        Bosch-Lluis, and I. Ramos-Perez, “An Efficient Algorithm to the Simulation 
        of Delay–Doppler Maps of Reflected Global Navigation Satellite System 
        Signals,” IEEE Transactions on Geoscience and Remote Sensing, vol. 47, no. 
        8, pp. 2733–2740, Aug. 2009.  
    ''' 
    delay_sp =  time_delay(np.array([0,0,0]))
    f_doppler_sp = doppler_shift(np.array([0,0,0]))
    return waf_delay(delay)**2 * np.abs(waf_frequency(f_doppler_sp - f_doppler))**2

#def waf_delay(delay):
#    ''' 
#    Voronovich implementation
#    '''
#    return np.where(np.abs(delay) <= delay_chip*(1+delay_chip/integration_time),
#                    1 - np.abs(delay)/delay_chip,
#                    -delay_chip/integration_time)

def waf_delay(delay):
    '''
    Defined in the text after Equation 1.
        J. F. Marchan-Hernandez, A. Camps, N. Rodriguez-Alvarez, E. Valencia, X.  
        Bosch-Lluis, and I. Ramos-Perez, “An Efficient Algorithm to the 
        Simulation of Delay–Doppler Maps of Reflected Global Navigation 
        Satellite System Signals,” IEEE Transactions on Geoscience and Remote 
        Sensing, vol. 47, no.  8, pp. 2733–2740, Aug. 2009.  
    '''
    t = np.where(np.abs(delay) <= delay_chip, 1 - np.abs(delay/delay_chip), 0)
    return t

#def waf_frequency(frequency_increment):
#    '''
#    Voronovich implementation
#    '''
#    return np.sin(np.pi*frequency_increment*integration_time) / \
#                 (np.pi*frequency_increment*integration_time) * \
#                 np.exp(-np.pi*1j*frequency_increment*integration_time)

def waf_frequency(f):
    '''
    Defined in the text after Equation 1.
        J. F. Marchan-Hernandez, A. Camps, N. Rodriguez-Alvarez, E. Valencia, X.  
        Bosch-Lluis, and I. Ramos-Perez, “An Efficient Algorithm to the 
        Simulation of Delay–Doppler Maps of Reflected Global Navigation 
        Satellite System Signals,” IEEE Transactions on Geoscience and Remote 
        Sensing, vol. 47, no.  8, pp. 2733–2740, Aug. 2009.  
    '''
    f *= integration_time
    sol =  np.where(np.abs(f) <= 0.2, 
            np.abs(1-(np.pi**2*f**2)/6+(np.pi**4*f**4)/120), # Taylor expansion around 0
            np.abs(np.sin(np.pi*f)/(np.pi*f))
            )
    return sol

def sigma(delay, f_doppler):
    '''
    Implements Equation 10.
        J. F. Marchan-Hernandez, A. Camps, N. Rodriguez-Alvarez, E. Valencia, X.  
        Bosch-Lluis, and I. Ramos-Perez, “An Efficient Algorithm to the 
        Simulation of Delay–Doppler Maps of Reflected Global Navigation 
        Satellite System Signals,” IEEE Transactions on Geoscience and Remote 
        Sensing, vol. 47, no.  8, pp. 2733–2740, Aug. 2009.  
    '''
    x_1 = x_delay_doppler_1(delay, f_doppler).real
    y_1 = y_delay_doppler_1(delay, f_doppler).real
    r_1 = np.array([x_1,y_1,0])

    x_2 = x_delay_doppler_2(delay, f_doppler).real
    y_2 = y_delay_doppler_2(delay, f_doppler).real
    r_2 = np.array([x_2,y_2,0])

    result = receiver_antenna_gain(r_1) + receiver_antenna_gain(r_2)

    #result =  transmitting_power*integration_time**2/(4*np.pi) * ( \
    #            radar_cross_section(r_1)/( \
    #                np.linalg.norm(r_1-r_t)**2* \
    #                np.linalg.norm(r_r-r_1)**2 \
    #            ) * \
    #            delay_doppler_jacobian_1(delay, f_doppler) * \
    #            receiver_antenna_gain(r_1) * \
    #            transmitting_antenna_gain(r_1) \
    #            + 
    #            radar_cross_section(r_2)/( \
    #                np.linalg.norm(r_2-r_t)**2* \
    #                np.linalg.norm(r_r-r_2)**2 \
    #            ) * \
    #            delay_doppler_jacobian_2(delay, f_doppler) * \
    #            receiver_antenna_gain(r_2) * \
    #            transmitting_antenna_gain(r_2) \
    #        )
    np.place(result, result < 0, 0)
    return result.real

antenna = tds_antenna()
def receiver_antenna_gain(r):
    #r_antenna = r[0:2] - r_r[0:2] 
    r_antenna = r[0:2]
    elevation = np.arctan2(h_r,np.sqrt(r_antenna[0]**2+r_antenna[1]**2))*180/np.pi
    azimuth = np.arctan2(-r_antenna[1],-r_antenna[0])*180/np.pi
    import pdb; pdb.set_trace() # break
    return antenna.gain(azimuth, elevation)

def transmitting_antenna_gain(r):
    #distance_to_sp = np.linalg.norm(r[0:2])
    #s = np.cos(distance_to_sp/200e3*np.pi/2)
    #if s < 0:
    #    s = 0
    return 1

def doppler_shift(r):
    ''' 
    Doppler shift as a contribution of the relative motion of transmitter and 
    receiver as well as the reflection point. 

    Implements Equation 14
        V. U. Zavorotny and A. G. Voronovich, “Scattering of GPS signals from 
        the ocean with wind remote sensing application,” IEEE Transactions on 
        Geoscience and Remote Sensing, vol. 38, no. 2, pp. 951–964, Mar. 2000.  
    '''
    #wavelength = light_speed/f_carrier
    #f_D_0 = (1/wavelength)*( \
    #           +np.inner(v_t, incident_direction(r)) \
    #           -np.inner(v_r, scattered_direction(r))
    #        )
    #f_surface = scattering_vector(r)*v_surface(r)/2*pi
    f_surface = 0
    # With the very far way transmitter approximation:
    f_D_0 =  f_carrier / light_speed * (-v_t[1] * np.cos(elevation) - v_t[2] * np.sin(elevation) + (v_r[0] * r[0] + v_r[1] * (r[1] + h_r / np.tan(elevation)) - v_r[2] * h_r) * (r[0] ** 2 + (r[1] + h_r / np.tan(elevation)) ** 2 + h_r ** 2) ** (-0.1e1 / 0.2e1))
    return f_D_0 + f_surface

def doppler_shift_increment(r):
    return doppler_shift(r) - doppler_shift(np.array([0,0,0]))

def scattering_vector(r):
    '''
    Implements Equation 7
        V. U. Zavorotny and A. G. Voronovich, “Scattering of GPS signals from 
        the ocean with wind remote sensing application,” IEEE Transactions on 
        Geoscience and Remote Sensing, vol. 38, no. 2, pp. 951–964, Mar. 2000.  
    '''
    scattering_vector = (scattered_direction(r) - incident_direction(r))
    scattering_vector_norm = np.linalg.norm(scattering_vector)
    scattering_vector[0] /= scattering_vector_norm
    scattering_vector[1] /= scattering_vector_norm
    scattering_vector[2] /= scattering_vector_norm
    return scattering_vector

def scattered_direction(r):
    scattered_direction = (r_r - r)
    scattered_direction_norm = np.linalg.norm(r_r - r)
    scattered_direction[0] /= scattered_direction_norm
    scattered_direction[1] /= scattered_direction_norm
    scattered_direction[2] /= scattered_direction_norm
    return scattered_direction

#def incident_direction(r):
#    return (r - r_t)/np.linalg.norm(r - r_t)

def incident_direction(r):
    incident_direction = (r - r_t)
    incident_direction_norm = np.linalg.norm(r - r_t)
    incident_direction[0] /= incident_direction_norm
    incident_direction[1] /= incident_direction_norm
    incident_direction[2] /= incident_direction_norm
    return  incident_direction

def time_delay(r):
    return (np.sqrt(r[0] ** 2 + (r[1] + h_r / np.tan(elevation)) ** 2 + h_r ** 2) - h_r / np.sin(elevation) - r[1] * np.cos(elevation)) / light_speed

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

    q = scattering_vector(r)
    q_norm = np.linalg.norm(scattering_vector(r))
    q_tangent = q[0:2]
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
    wind_rotation = np.array([
        [np.cos(phi_0), -np.sin(phi_0)],
        [np.sin(phi_0),  np.cos(phi_0)]
    ])
    covariance = np.array([
        [variance_upwind(u_10), 0],
        [0, variance_crosswind(u_10)]
    ])
    w = (wind_rotation.dot(covariance)).dot(np.transpose(wind_rotation))
    return 1/np.sqrt(np.linalg.det(2*np.pi*w)) * \
            np.exp(
                -1/2*(np.transpose(x).dot(np.linalg.inv(w))).dot(x)
            )

def variance_upwind(u_10):
    ''' 
    Based on the 'clean surface mean square slope model' of Cox and Mux
    Implements Equation 4
        Q. Yan and W. Huang, “GNSS-R Delay-Doppler Map Simulation Based on the 
        2004 Sumatra-Andaman Tsunami Event,” Journal of Sensors, vol. 2016, pp. 
        1–14, 2016.  

    Args: 
        u_10:   Wind speed at 10 meters above sea surface

    Returns:
        upwind variance
    '''
    f = np.piecewise(u_10, 
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
    return (0.45*(3.16e-3*f))

def variance_crosswind(wind_speed_10m_above_sea):
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
    f = np.piecewise(u_10, 
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
    return (0.45*(0.003 + 1.92e-3*f))

# --------------------

# Plotting Area

x_0 =  -200e3 # meters
x_1 =  200e3 # meters
n_x = 500

y_0 =  -200e3 # meters
y_1 =  200e3 # meters
n_y = 500

x_grid, y_grid = np.meshgrid(
   np.linspace(x_0, x_1, n_x), 
   np.linspace(y_0, y_1, n_y)
   )

r = [x_grid, y_grid, 0]
z_grid_delay = time_delay(r)/delay_chip
z_grid_doppler_increment = doppler_shift(r) - doppler_shift(np.array([0,0,0]))

delay_sp =  time_delay(np.array([0,0,0]))
f_doppler_sp = doppler_shift(np.array([0,0,0]))

delay_values = list(np.arange(delay_start, delay_end, delay_resolution))
doppler_values = list(np.arange(
                        f_doppler_sp + doppler_start, 
                        f_doppler_sp + doppler_end, 
                        doppler_resolution
                        ))
doppler_values_increment = doppler_values - doppler_shift(np.array([0,0,0]))

fig_lines, ax_lines = plt.subplots(1,figsize=(10, 4))
contour_delay = ax_lines.contour(x_grid, y_grid, z_grid_delay, [i/delay_chip for i in delay_values], cmap='winter')
fig_lines.colorbar(contour_delay, label='C/A chips', )

contour_doppler = ax_lines.contour(x_grid, y_grid, z_grid_doppler_increment, doppler_values_increment, cmap='winter')
fig_lines.colorbar(contour_doppler, label='Hz', )

ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/1000))
ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000))
ax_lines.xaxis.set_major_formatter(ticks_x)
ax_lines.yaxis.set_major_formatter(ticks_y)
plt.xlabel('[km]')
plt.ylabel('[km]')

#print('Finding intersection for d:{0}, f:{1}'.format(delay_values[4], doppler_values[2]))
#x_s_1 = x_delay_doppler_2(delay_values[4], doppler_values[int(len(doppler_values)/2)])
#y_s_1 = y_delay_doppler_2(delay_values[4], doppler_values[int(len(doppler_values)/2)])
#x_s_2 = x_delay_doppler_2(delay_values[8], doppler_values[int(len(doppler_values)/2)])
#y_s_2 = y_delay_doppler_2(delay_values[8], doppler_values[int(len(doppler_values)/2)])
#ax_lines.scatter(x_s_1, y_s_1, s=70, marker=(5, 2), zorder=4)
#ax_lines.scatter(x_s_2, y_s_2, s=70, marker=(5, 2), zorder=4)

delay_grid, doppler_grid = np.meshgrid(delay_values, doppler_values)

waf_matrix = np.zeros([len(delay_values), len(doppler_values)])
waf_matrix = waf_squared(delay_grid, doppler_grid)

sigma_matrix = np.zeros([len(delay_values), len(doppler_values)])
sigma_matrix = sigma(delay_grid, doppler_grid)

ddm = np.zeros([len(delay_values), len(doppler_values)])
ddm = signal.convolve2d(waf_matrix, sigma_matrix)

fig_waf, ax_waf = plt.subplots(1,figsize=(10, 4))
ax_waf.set_title('WAF')
im = ax_waf.imshow(waf_matrix, cmap='viridis', 
        extent=(delay_start/delay_chip, delay_end/delay_chip, doppler_start, doppler_end),
        aspect="auto"
    )

fig_sigma, ax_sigma = plt.subplots(1,figsize=(10, 4))
ax_sigma.set_title('Sigma')
im = ax_sigma.imshow(sigma_matrix, cmap='viridis', 
        extent=(delay_start/delay_chip, delay_end/delay_chip, doppler_start, doppler_end),
        aspect="auto"
    )

fig_ddm, ax_ddm = plt.subplots(1,figsize=(10, 4))
ax_ddm.set_title('DDM')
im = ax_ddm.imshow(ddm, cmap='viridis', 
        extent=(delay_start/delay_chip, delay_end/delay_chip, doppler_start, doppler_end),
        aspect="auto"
    )

fig_waf_delay, ax_waf_delay = plt.subplots(1,figsize=(10, 4))
waf_delay_values = np.zeros(len(delay_values))
waf_delay_values = waf_delay(np.array(delay_values))**2
ax_waf_delay.plot([i/delay_chip for i in delay_values], waf_delay_values)
ax_waf_delay.set_title('waf_delay')

fig_waf_frequency, ax_waf_frequency = plt.subplots(1,figsize=(10, 4))
waf_frequency_values = np.zeros(len(doppler_values))
waf_frequency_values = waf_frequency(f_doppler_sp - np.array(doppler_values))**2
ax_waf_frequency.plot(doppler_values_increment, waf_frequency_values)
ax_waf_frequency.set_title('waf_freq')

plt.show()

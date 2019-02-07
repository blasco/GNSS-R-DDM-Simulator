#!/usr/bin/env python

import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from gnssr.tds.antenna.gain import *
from scipy import signal
from delay_doppler_jacobian import *
from problem_definition import *
import cv2

# -----------------------
# Bistatic radar equation
# ----------------------

def waf_squared(delay, doppler): 
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
    return waf_delay(delay)**2 * np.abs(waf_frequency(doppler))**2

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

def sigma(delay, doppler):
    '''
    Implements Equation 10.
        J. F. Marchan-Hernandez, A. Camps, N. Rodriguez-Alvarez, E. Valencia, X.  
        Bosch-Lluis, and I. Ramos-Perez, “An Efficient Algorithm to the 
        Simulation of Delay–Doppler Maps of Reflected Global Navigation 
        Satellite System Signals,” IEEE Transactions on Geoscience and Remote 
        Sensing, vol. 47, no.  8, pp. 2733–2740, Aug. 2009.  
    '''
    x_1 = x_delay_doppler_1(delay, doppler).real
    y_1 = y_delay_doppler_1(delay, doppler).real
    r_1 = np.array([x_1,y_1,0])
    x_2 = x_delay_doppler_2(delay, doppler).real
    y_2 = y_delay_doppler_2(delay, doppler).real
    r_2 = np.array([x_2,y_2,0])

    result =  transmitting_power*integration_time**2/(4*np.pi) * ( \
                radar_cross_section(r_1)/( \
                    np.linalg.norm(r_1-r_t)**2* \
                    np.linalg.norm(r_r-r_1)**2 \
                ) * \
                delay_doppler_jacobian_1(delay, doppler) * \
                receiver_antenna_gain(r_1) * \
                transmitting_antenna_gain(r_1) \
                + 
                radar_cross_section(r_2)/( \
                    np.linalg.norm(r_2-r_t)**2* \
                    np.linalg.norm(r_r-r_2)**2 \
                ) * \
                delay_doppler_jacobian_2(delay, doppler) * \
                receiver_antenna_gain(r_2) * \
                transmitting_antenna_gain(r_2) \
            )
    return result.real

antenna = tds_antenna_gain()
def receiver_antenna_gain(r):
    r_antenna = r[0:2] - r_r[0:2] 
    #r_antenna = r[0:2]
    elevation = np.arctan2(h_r,np.sqrt(r_antenna[0]**2+r_antenna[1]**2))*180/np.pi
    azimuth = np.arctan2(-r_antenna[1],-r_antenna[0])*180/np.pi
    return antenna.gain(azimuth, elevation)

def transmitting_antenna_gain(r):
    return gps_isotropic_antenna_gain

def eq_doppler_absolute_shift(r):
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

def eq_doppler_absolute_shift_increment(r):
    return eq_doppler_absolute_shift(r) - eq_doppler_absolute_shift(np.array([0,0,0]))

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

def incident_direction(r):
    incident_direction = (r - r_t)
    incident_direction_norm = np.linalg.norm(r - r_t)
    incident_direction[0] /= incident_direction_norm
    incident_direction[1] /= incident_direction_norm
    incident_direction[2] /= incident_direction_norm
    return  incident_direction

def eq_delay_incremet(r):
    return (np.sqrt(r[0] ** 2 + (r[1] + h_r / np.tan(elevation)) ** 2 + h_r ** 2) - h_r / np.sin(elevation) - r[1] * np.cos(elevation)) / light_speed

def radar_cross_section(r):
    return rcs_sea(r)

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
    hermite_coeficients[1,2] = (0.01 - 0.0086*(u_10))
    hermite_coeficients[3,0] = (0.04 -0.033*(u_10))
    hermite_coeficients[0,4] = (0.4)
    hermite_coeficients[2,2] = (0.12)
    hermite_coeficients[4,0] = (0.23)

    result = 1/(2*np.pi*rms_u*rms_c) * \
            np.exp(
                -1/2*((s[0]/rms_u)**2+(s[1]/rms_c)**2) \
            ) * \
            np.polynomial.hermite.hermval2d(s[0]/rms_u, s[1]/rms_c, hermite_coeficients)

    #eta = s[0]/rms_u
    #xi = s[1]/rms_c
    #result = 1/(2*np.pi*rms_u*rms_c)*np.exp(-1/2*(xi**2+eta**2)) * \
    #        ( \
    #            1 - \
    #            1/2*c_21*(xi**2-1)*eta - \
    #            1/6*c_03*(eta**3-3*eta) + \
    #            1/24*c_40*(xi**4-6*xi**2+3) + \
    #            1/4*c_22*(xi**2-1)*(eta**2-1) + \
    #            1/24*c_04*(eta**4-6*eta**2+3) \
    #        )
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

# TODO: Fresnel model:
# Di Simone > It even contains the Klein-Swift model for computation of the dielectric constant of sea water:
# CYGNSS - Algorithm Theoretical Basis Document Level 2 Mean-Square Slope Retrieval
#This may be useful for computing the reflection coefficient R, whose expressions are reported here (eq. 12.23-12.26)
#F. Ulaby, R. Moore, and A. Fung, Microwave Remote Sensing, Active
#and Passive. Vol. II: Radar Remote Sensing and Surface Scattering and
#Emission Theory. Reading, MA: Addison-Wesley, 1982.
#to be mapped in circular polarization as explained in (eq. (9))
#Di Bisceglie, M., Di Martino, G., Di Simone, A., Galdi, C., Iodice, A., Riccio, D., & Ruella, G. (2018, July). Two-Scale Model for the Evaluation of Sea-Surface Scattering in GNSS-R Ship-Detection Applications. In IGARSS 2018-2018 IEEE International Geoscience and Remote Sensing Symposium (pp. 3181-3184). IEEE.

def fresnel_model(r):
    return 1

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

def main():

    # Surface mesh
    x_0 =  -150e3 # meters
    x_1 =  150e3 # meters
    n_x = 500
    y_0 =  -150e3 # meters
    y_1 =  150e3 # meters
    n_y = 500
    x_grid, y_grid = np.meshgrid(
       np.linspace(x_0, x_1, n_x), 
       np.linspace(y_0, y_1, n_y)
       )

    # Dealy-Doppler values
    delay_increment_values = list(np.arange(delay_increment_start, delay_increment_end, delay_resolution))
    doppler_increment_values = list(np.arange(doppler_increment_start, doppler_increment_end, doppler_resolution))
    waf_delay_increment_values = list(np.arange(-delay_increment_end, delay_increment_end, delay_resolution))
    waf_doppler_increment_values = list(np.arange(doppler_increment_start, doppler_increment_end, doppler_resolution))
    doppler_absolute_values = doppler_increment_values + eq_doppler_absolute_shift(np.array([0,0,0]))

    r = np.array([x_grid, y_grid, 0])

    # Iso-delay and iso-doppler maps
    z_grid_delay_chip = eq_delay_incremet(r)/delay_chip
    z_grid_doppler_increment = eq_doppler_absolute_shift(r) - eq_doppler_absolute_shift(np.array([0,0,0]))

    # Receiver Antenna Gain
    z_antenna = receiver_antenna_gain(r)

    # Radar cross section
    z_rcs = radar_cross_section(r)

    # Iso lines plot
    fig_isolines, ax_isolines = plt.subplots(1,figsize=(10, 4))

    contour_delay_chip = ax_isolines.contour(
            x_grid, y_grid, z_grid_delay_chip, 
            np.arange(0, delay_increment_end/delay_chip, 0.5), 
            cmap='winter', alpha = 0.6
            )
    contour_doppler = ax_isolines.contour(
            x_grid, y_grid, z_grid_doppler_increment, 
            np.arange(doppler_increment_start, doppler_increment_end, 250), 
            cmap='jet', alpha = 0.8
            )

    test_delay = np.array([12*delay_chip])
    test_doppler = np.array([1000]) +  eq_doppler_absolute_shift(np.array([0,0,0]))

    print('Finding intersection for d:{0}, f:{1}'.format(test_delay, test_delay))
    x_s_1 = x_delay_doppler_1(test_delay, test_doppler)
    y_s_1 = y_delay_doppler_1(test_delay, test_doppler)
    x_s_2 = x_delay_doppler_2(test_delay, test_doppler)
    y_s_2 = y_delay_doppler_2(test_delay, test_doppler)
    ax_isolines.scatter(x_s_1, y_s_1, s=70, marker=(5, 2), zorder=4)
    ax_isolines.scatter(x_s_2, y_s_2, s=70, marker=(5, 2), zorder=4)

    fig_isolines.colorbar(contour_delay_chip, label='C/A chips')
    fig_isolines.colorbar(contour_doppler, label='Hz')
    ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/1000))
    ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000))
    ax_isolines.xaxis.set_major_formatter(ticks_x)
    ax_isolines.yaxis.set_major_formatter(ticks_y)
    plt.xlabel('[km]')
    plt.ylabel('[km]')

    # Antenna pattern plot
    fig_antenna, ax_antenna = plt.subplots(1,figsize=(10, 4))
    ax_antenna.set_title('Antenna')

    contour_delay_chip = ax_antenna.contour(
            x_grid, y_grid, z_grid_delay_chip, 
            np.arange(0, delay_increment_end/delay_chip, 0.5), 
            cmap='winter', alpha = 0.3
            )
    contour_doppler = ax_antenna.contour(
            x_grid, y_grid, z_grid_doppler_increment, 
            np.arange(doppler_increment_start, doppler_increment_end, 250), 
            cmap='jet', alpha = 0.3
            )
    contour_antenna = ax_antenna.contourf(x_grid, y_grid, z_antenna, 55, cmap='jet', alpha = 1.0)

    test_delay = np.array([12*delay_chip])
    test_doppler = np.array([1000]) +  eq_doppler_absolute_shift(np.array([0,0,0]))
    print('Finding intersection for d:{0}, f:{1}'.format(test_delay, test_delay))

    fig_antenna.colorbar(contour_delay_chip, label='C/A chips')
    fig_antenna.colorbar(contour_doppler, label='Hz')
    fig_antenna.colorbar(contour_antenna, label='Gain')
    ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/1000))
    ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000))
    ax_antenna.xaxis.set_major_formatter(ticks_x)
    ax_antenna.yaxis.set_major_formatter(ticks_y)
    plt.xlabel('[km]')
    plt.ylabel('[km]')

    # RCS surface plot
    fig_rcs, ax_rcs = plt.subplots(1,figsize=(10, 4))
    ax_rcs.set_title('RCS')

    #contour_delay_chip = ax_rcs.contour(
    #        x_grid, y_grid, z_grid_delay_chip, 
    #        np.arange(0, delay_increment_end/delay_chip, 1), 
    #        cmap='jet', alpha = 0.5
    #        )
    #contour_doppler = ax_rcs.contour(
    #        x_grid, y_grid, z_grid_doppler_increment, 
    #        np.arange(doppler_increment_start, doppler_increment_end, 250), 
    #        cmap='jet', alpha = 0.5
    #        )
    contour_rcs = ax_rcs.contourf(x_grid, y_grid, z_rcs, 55, cmap='viridis', alpha = 1.0)

    test_delay = np.array([12*delay_chip])
    test_doppler = np.array([1000]) +  eq_doppler_absolute_shift(np.array([0,0,0]))
    print('Finding intersection for d:{0}, f:{1}'.format(test_delay, test_delay))

    #fig_rcs.colorbar(contour_delay_chip, label='C/A chips')
    #fig_rcs.colorbar(contour_doppler, label='Hz')
    fig_rcs.colorbar(contour_rcs, label='RCS')
    ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/1000))
    ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000))
    ax_rcs.xaxis.set_major_formatter(ticks_x)
    ax_rcs.yaxis.set_major_formatter(ticks_y)
    plt.xlabel('[km]')
    plt.ylabel('[km]')

    # Delay-Doppler mesh
    delay_grid, doppler_grid = np.meshgrid(delay_increment_values, doppler_absolute_values)
    waf_delay_grid, waf_doppler_grid = np.meshgrid(waf_delay_increment_values, waf_doppler_increment_values)

    # Retrived power
    waf_matrix = waf_squared(waf_delay_grid, waf_doppler_grid)
    sigma_matrix = sigma(delay_grid, doppler_grid)
    ddm = signal.convolve2d(sigma_matrix, waf_matrix, mode='same')

    fig_waf, ax_waf = plt.subplots(1,figsize=(10, 4))
    ax_waf.set_title('WAF')
    im = ax_waf.imshow(waf_matrix, cmap='viridis', 
            extent=(-delay_increment_end/delay_chip, delay_increment_end/delay_chip, doppler_increment_start, doppler_increment_end),
            aspect="auto"
            )

    fig_sigma, ax_sigma = plt.subplots(1,figsize=(10, 4))
    ax_sigma.set_title('Sigma')
    im = ax_sigma.imshow(sigma_matrix, cmap='viridis', 
            extent=(delay_increment_start/delay_chip, delay_increment_end/delay_chip, doppler_increment_start, doppler_increment_end),
            aspect="auto"
            )

    fig_ddm, ax_ddm = plt.subplots(1,figsize=(10, 4))
    ax_ddm.set_title('DDM')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    im = ax_ddm.imshow(ddm, cmap='jet', 
            extent=(delay_increment_start/delay_chip, delay_increment_end/delay_chip, doppler_increment_start, doppler_increment_end),
            aspect="auto"
            )

    # Image downscaling to desired resolution:
    # TODO: This is just an average of the pixels around the area
    # This is not valid, summation i srequired:
    # Di Simone > From a physical viewpoint, 
    # such an approach should call for summation instead of averaging
    # https://stackoverflow.com/questions/48121916/numpy-resize-rescale-image
    fig_ddm_rescaled, ax_ddm_rescaled = plt.subplots(1,figsize=(10, 4))
    ax_ddm_rescaled.set_title('DDM rescaled')
    rescaled_doppler_resolution = 500
    rescaled_delay_resolution_chips = 0.3
    ddm_rescaled = cv2.resize(ddm, 
            dsize=(
                int((delay_increment_end/delay_chip - delay_increment_start/delay_chip)/rescaled_delay_resolution_chips), 
                int((doppler_increment_end - doppler_increment_start)/rescaled_doppler_resolution)
                ), 
            interpolation=cv2.INTER_AREA
            )
    im = ax_ddm_rescaled.imshow(ddm_rescaled, cmap='jet', 
            extent=(delay_increment_start/delay_chip, delay_increment_end/delay_chip, doppler_increment_start, doppler_increment_end),
            aspect="auto"
            )

    fig_waf_delay, ax_waf_delay = plt.subplots(1,figsize=(10, 4))
    waf_delay_result = waf_delay(np.array(waf_delay_increment_values))**2
    ax_waf_delay.plot([i/delay_chip for i in waf_delay_increment_values], waf_delay_result)
    ax_waf_delay.set_title('waf_delay')

    fig_waf_frequency, ax_waf_frequency = plt.subplots(1,figsize=(10, 4))
    waf_frequency_result = waf_frequency(np.array(waf_doppler_increment_values))**2
    ax_waf_frequency.plot(waf_doppler_increment_values, waf_frequency_result)
    ax_waf_frequency.set_title('waf_freq')

    plt.show()

if __name__ == '__main__':
    main()

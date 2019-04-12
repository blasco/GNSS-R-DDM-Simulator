#!/usr/bin/env python

''' 
Jacobians computed as defined in equations 11 and 12. 
x(delay,doppler) and y(delay,doppler) form the solutions to equations 17 and 18 
generated with Maple 

    J. F. Marchan-Hernandez, A. Camps, N. Rodriguez-Alvarez, E.  Valencia, X.  
    Bosch-Lluis, and I. Ramos-Perez, “An Efficient Algorithm to the Simulation 
    of Delay–Doppler Maps of Reflected Global Navigation Satellite System 
    Signals,” IEEE Transactions on Geoscience and Remote Sensing, vol. 47, no. 
    8, pp.  2733–2740, Aug. 2009.  
'''

import numpy as np
from gnssr.utils import *

def delay_doppler_jacobian_1(delay, f_doppler, sim_config):
    delay_resolution = sim_config.delay_resolution
    doppler_resolution = sim_config.doppler_resolution

    j_11 = (1)*(
             x_delay_doppler_1(delay+delay_resolution/2, f_doppler, sim_config) - \
             x_delay_doppler_1(delay-delay_resolution/2, f_doppler, sim_config)
           )
    j_12 = (1)*( 
             x_delay_doppler_1(delay, f_doppler+doppler_resolution/2, sim_config) - \
             x_delay_doppler_1(delay, f_doppler-doppler_resolution/2, sim_config)
            )
    j_21 = (1)*(
             y_delay_doppler_1(delay+delay_resolution/2, f_doppler, sim_config) - \
             y_delay_doppler_1(delay-delay_resolution/2, f_doppler, sim_config)
            )
    j_22 = (1)*(
             y_delay_doppler_1(delay, f_doppler+doppler_resolution/2, sim_config) - \
             y_delay_doppler_1(delay, f_doppler-doppler_resolution/2, sim_config)
            )
    #return np.absolute(j_11*j_22 - j_12*j_21)
    J = j_11*j_22 - j_12*j_21
    np.place(J, J.imag != 0, 0)
    return np.absolute(J).astype(float)

def delay_doppler_jacobian_2(delay, f_doppler, sim_config):
    delay_resolution = sim_config.delay_resolution
    doppler_resolution = sim_config.doppler_resolution

    j_11 = (1/delay_resolution)* \
            (x_delay_doppler_2(delay+delay_resolution/2, f_doppler, sim_config) - \
             x_delay_doppler_2(delay-delay_resolution/2, f_doppler, sim_config))
    j_12 = (1/doppler_resolution)* \
            (x_delay_doppler_2(delay, f_doppler+doppler_resolution/2, sim_config) - \
             x_delay_doppler_2(delay, f_doppler-doppler_resolution/2, sim_config))
    j_21 = (1/delay_resolution)* \
            (y_delay_doppler_2(delay+delay_resolution/2, f_doppler, sim_config) - \
             y_delay_doppler_2(delay-delay_resolution/2, f_doppler, sim_config))
    j_22 = (1/doppler_resolution)* \
            (y_delay_doppler_2(delay, f_doppler+doppler_resolution/2, sim_config) - \
             y_delay_doppler_2(delay, f_doppler-doppler_resolution/2, sim_config))
    #return np.absolute(j_11*j_22 - j_12*j_21)
    J = j_11*j_22 - j_12*j_21
    np.place(J, J.imag != 0, 0)
    return np.absolute(J).astype(float)

def x_delay_doppler_1(delay, f_doppler, sim_config): 
    v_t = sim_config.v_t
    v_r = sim_config.v_r
    elevation = sim_config.elevation
    f_carrier = sim_config.f_carrier
    h_r = sim_config.h_r

    delay = delay.astype(complex)
    f_doppler = f_doppler.astype(complex)

    sol =  ((np.cos(elevation) ** 2 * v_t[1] * f_carrier + (np.sin(elevation) * v_t[2] * f_carrier + f_doppler * light_speed) * np.cos(elevation) - v_r[1] * f_carrier) * np.sqrt(-np.sin(elevation) ** 2 * v_r[0] ** 2 * f_carrier ** 2 * (-4 * np.cos(elevation) ** 3 * delay * h_r * v_t[1] * v_t[2] * f_carrier ** 2 * light_speed + f_carrier * (2 * h_r * light_speed * delay * f_carrier * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2]) * np.sin(elevation) + ((-v_t[2] + v_r[1] - v_t[1] - v_r[2]) * (v_t[2] + v_r[1] - v_t[1] + v_r[2]) * h_r ** 2 + light_speed ** 2 * delay ** 2 * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2])) * f_carrier - 4 * h_r * f_doppler * light_speed ** 2 * v_t[2] * delay) * np.cos(elevation) ** 2 - 2 * ((((v_r[2] + v_t[2]) * (v_r[1] - v_t[1]) * h_r ** 2 - light_speed ** 2 * v_t[1] * v_t[2] * delay ** 2) * f_carrier - 2 * h_r * f_doppler * light_speed ** 2 * v_t[1] * delay) * np.sin(elevation) + light_speed * (delay * (-2 * v_t[1] * v_t[2] + v_r[2] * (v_r[1] - v_t[1])) * h_r * f_carrier + ((v_r[1] - v_t[1]) * h_r ** 2 - light_speed ** 2 * v_t[1] * delay ** 2) * f_doppler)) * f_carrier * np.cos(elevation) + 2 * light_speed * (-h_r * delay * (v_r[0] ** 2 + v_r[1] ** 2 - v_r[2] * v_t[2] - v_t[2] ** 2) * f_carrier ** 2 + f_doppler * ((v_r[2] + v_t[2]) * h_r ** 2 + light_speed ** 2 * v_t[2] * delay ** 2) * f_carrier + h_r * f_doppler ** 2 * light_speed ** 2 * delay) * np.sin(elevation) + ((v_r[2] + v_t[2]) ** 2 * h_r ** 2 - light_speed ** 2 * delay ** 2 * (v_r[0] ** 2 + v_r[1] ** 2 - v_t[2] ** 2)) * f_carrier ** 2 + 2 * h_r * f_doppler * light_speed ** 2 * delay * (v_r[2] + 2 * v_t[2]) * f_carrier + f_doppler ** 2 * light_speed ** 2 * (delay ** 2 * light_speed ** 2 + h_r ** 2))) - v_r[0] ** 2 * f_carrier ** 2 * (-h_r * f_carrier * (v_r[1] - v_t[1]) * np.cos(elevation) ** 3 + (h_r * f_carrier * (v_r[2] + v_t[2]) * np.sin(elevation) + light_speed * (delay * v_t[2] * f_carrier + h_r * f_doppler)) * np.cos(elevation) ** 2 + f_carrier * (v_r[1] - v_t[1]) * (delay * np.sin(elevation) * light_speed + h_r) * np.cos(elevation) + (-h_r * f_carrier * (v_r[2] + v_t[2]) - f_doppler * light_speed ** 2 * delay) * np.sin(elevation) - light_speed * (delay * v_t[2] * f_carrier + h_r * f_doppler))) / np.sin(elevation) / v_r[0] / (f_carrier ** 2 * (v_t[1] ** 2 - v_t[2] ** 2) * np.cos(elevation) ** 4 + 2 * v_t[1] * f_carrier * (np.sin(elevation) * v_t[2] * f_carrier + f_doppler * light_speed) * np.cos(elevation) ** 3 + (2 * v_t[2] * np.sin(elevation) * f_carrier * f_doppler * light_speed + (-v_r[0] ** 2 - 2 * v_r[1] * v_t[1] + v_t[2] ** 2) * f_carrier ** 2 + f_doppler ** 2 * light_speed ** 2) * np.cos(elevation) ** 2 - 2 * v_r[1] * f_carrier * (np.sin(elevation) * v_t[2] * f_carrier + f_doppler * light_speed) * np.cos(elevation) + f_carrier ** 2 * (v_r[0] ** 2 + v_r[1] ** 2)) / f_carrier
    #np.place(sol, sol.imag > 0, np.nan)
    return sol

def x_delay_doppler_2(delay, f_doppler, sim_config): 
    v_t = sim_config.v_t
    v_r = sim_config.v_r
    elevation = sim_config.elevation
    f_carrier = sim_config.f_carrier
    h_r = sim_config.h_r

    delay = delay.astype(complex)
    f_doppler = f_doppler.astype(complex)

    sol =  ((-np.cos(elevation) ** 2 * v_t[1] * f_carrier + (-np.sin(elevation) * v_t[2] * f_carrier - f_doppler * light_speed) * np.cos(elevation) + v_r[1] * f_carrier) * np.sqrt(-np.sin(elevation) ** 2 * v_r[0] ** 2 * f_carrier ** 2 * (-4 * np.cos(elevation) ** 3 * delay * h_r * v_t[1] * v_t[2] * f_carrier ** 2 * light_speed + f_carrier * (2 * h_r * light_speed * delay * f_carrier * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2]) * np.sin(elevation) + ((-v_t[2] + v_r[1] - v_t[1] - v_r[2]) * (v_t[2] + v_r[1] - v_t[1] + v_r[2]) * h_r ** 2 + light_speed ** 2 * delay ** 2 * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2])) * f_carrier - 4 * h_r * f_doppler * light_speed ** 2 * v_t[2] * delay) * np.cos(elevation) ** 2 - 2 * ((((v_r[2] + v_t[2]) * (v_r[1] - v_t[1]) * h_r ** 2 - light_speed ** 2 * v_t[1] * v_t[2] * delay ** 2) * f_carrier - 2 * h_r * f_doppler * light_speed ** 2 * v_t[1] * delay) * np.sin(elevation) + light_speed * (delay * (-2 * v_t[1] * v_t[2] + v_r[2] * (v_r[1] - v_t[1])) * h_r * f_carrier + ((v_r[1] - v_t[1]) * h_r ** 2 - light_speed ** 2 * v_t[1] * delay ** 2) * f_doppler)) * f_carrier * np.cos(elevation) + 2 * light_speed * (-h_r * delay * (v_r[0] ** 2 + v_r[1] ** 2 - v_r[2] * v_t[2] - v_t[2] ** 2) * f_carrier ** 2 + f_doppler * ((v_r[2] + v_t[2]) * h_r ** 2 + light_speed ** 2 * v_t[2] * delay ** 2) * f_carrier + h_r * f_doppler ** 2 * light_speed ** 2 * delay) * np.sin(elevation) + ((v_r[2] + v_t[2]) ** 2 * h_r ** 2 - light_speed ** 2 * delay ** 2 * (v_r[0] ** 2 + v_r[1] ** 2 - v_t[2] ** 2)) * f_carrier ** 2 + 2 * h_r * f_doppler * light_speed ** 2 * delay * (v_r[2] + 2 * v_t[2]) * f_carrier + f_doppler ** 2 * light_speed ** 2 * (delay ** 2 * light_speed ** 2 + h_r ** 2))) - v_r[0] ** 2 * f_carrier ** 2 * (-h_r * f_carrier * (v_r[1] - v_t[1]) * np.cos(elevation) ** 3 + (h_r * f_carrier * (v_r[2] + v_t[2]) * np.sin(elevation) + light_speed * (delay * v_t[2] * f_carrier + h_r * f_doppler)) * np.cos(elevation) ** 2 + f_carrier * (v_r[1] - v_t[1]) * (delay * np.sin(elevation) * light_speed + h_r) * np.cos(elevation) + (-h_r * f_carrier * (v_r[2] + v_t[2]) - f_doppler * light_speed ** 2 * delay) * np.sin(elevation) - light_speed * (delay * v_t[2] * f_carrier + h_r * f_doppler))) / np.sin(elevation) / v_r[0] / (f_carrier ** 2 * (v_t[1] ** 2 - v_t[2] ** 2) * np.cos(elevation) ** 4 + 2 * v_t[1] * f_carrier * (np.sin(elevation) * v_t[2] * f_carrier + f_doppler * light_speed) * np.cos(elevation) ** 3 + (2 * v_t[2] * np.sin(elevation) * f_carrier * f_doppler * light_speed + (-v_r[0] ** 2 - 2 * v_r[1] * v_t[1] + v_t[2] ** 2) * f_carrier ** 2 + f_doppler ** 2 * light_speed ** 2) * np.cos(elevation) ** 2 - 2 * v_r[1] * f_carrier * (np.sin(elevation) * v_t[2] * f_carrier + f_doppler * light_speed) * np.cos(elevation) + f_carrier ** 2 * (v_r[0] ** 2 + v_r[1] ** 2)) / f_carrier
    #np.place(sol, sol.imag > 0, np.nan)
    return sol

def y_delay_doppler_1(delay, f_doppler, sim_config):
    v_t = sim_config.v_t
    v_r = sim_config.v_r
    elevation = sim_config.elevation
    f_carrier = sim_config.f_carrier
    h_r = sim_config.h_r

    delay = delay.astype(complex)
    f_doppler = f_doppler.astype(complex)

    sol =  (2 * np.cos(elevation) ** 4 * delay * v_t[1] * v_t[2] * f_carrier ** 2 * light_speed + f_carrier * (-light_speed * delay * f_carrier * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2]) * np.sin(elevation) + (v_t[2] ** 2 + v_r[2] * v_t[2] + v_t[1] * (v_r[1] - v_t[1])) * h_r * f_carrier + 2 * f_doppler * light_speed ** 2 * v_t[2] * delay) * np.cos(elevation) ** 3 + f_carrier * ((((v_r[1] - 2 * v_t[1]) * v_t[2] - v_r[2] * v_t[1]) * h_r * f_carrier - 2 * f_doppler * light_speed ** 2 * v_t[1] * delay) * np.sin(elevation) + light_speed * (-delay * f_carrier * (v_r[1] + 2 * v_t[1]) * v_t[2] + h_r * f_doppler * (v_r[1] - 2 * v_t[1]))) * np.cos(elevation) ** 2 + (-light_speed * (-delay * (v_r[0] ** 2 + v_r[1] * v_t[1] - v_t[2] ** 2) * f_carrier ** 2 + h_r * f_doppler * (v_r[2] + 2 * v_t[2]) * f_carrier + f_doppler ** 2 * light_speed ** 2 * delay) * np.sin(elevation) - (v_t[2] ** 2 + v_r[2] * v_t[2] + v_r[1] * (v_r[1] - v_t[1])) * h_r * f_carrier ** 2 - 2 * delay * v_t[2] * f_carrier * f_doppler * light_speed ** 2 - h_r * f_doppler ** 2 * light_speed ** 2) * np.cos(elevation) + (h_r * f_carrier * (v_r[2] + v_t[2]) + f_doppler * light_speed ** 2 * delay) * f_carrier * v_r[1] * np.sin(elevation) + delay * v_r[1] * v_t[2] * f_carrier ** 2 * light_speed + h_r * v_r[1] * f_carrier * f_doppler * light_speed + np.sqrt(-np.sin(elevation) ** 2 * v_r[0] ** 2 * f_carrier ** 2 * (-4 * np.cos(elevation) ** 3 * delay * h_r * v_t[1] * v_t[2] * f_carrier ** 2 * light_speed + f_carrier * (2 * h_r * light_speed * delay * f_carrier * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2]) * np.sin(elevation) + ((-v_t[2] + v_r[1] - v_t[1] - v_r[2]) * (v_t[2] + v_r[1] - v_t[1] + v_r[2]) * h_r ** 2 + light_speed ** 2 * delay ** 2 * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2])) * f_carrier - 4 * h_r * f_doppler * light_speed ** 2 * v_t[2] * delay) * np.cos(elevation) ** 2 - 2 * ((((v_r[2] + v_t[2]) * (v_r[1] - v_t[1]) * h_r ** 2 - light_speed ** 2 * v_t[1] * v_t[2] * delay ** 2) * f_carrier - 2 * h_r * f_doppler * light_speed ** 2 * v_t[1] * delay) * np.sin(elevation) + light_speed * (delay * (-2 * v_t[1] * v_t[2] + v_r[2] * (v_r[1] - v_t[1])) * h_r * f_carrier + ((v_r[1] - v_t[1]) * h_r ** 2 - light_speed ** 2 * v_t[1] * delay ** 2) * f_doppler)) * f_carrier * np.cos(elevation) + 2 * light_speed * (-h_r * delay * (v_r[0] ** 2 + v_r[1] ** 2 - v_r[2] * v_t[2] - v_t[2] ** 2) * f_carrier ** 2 + f_doppler * ((v_r[2] + v_t[2]) * h_r ** 2 + light_speed ** 2 * v_t[2] * delay ** 2) * f_carrier + h_r * f_doppler ** 2 * light_speed ** 2 * delay) * np.sin(elevation) + ((v_r[2] + v_t[2]) ** 2 * h_r ** 2 - light_speed ** 2 * delay ** 2 * (v_r[0] ** 2 + v_r[1] ** 2 - v_t[2] ** 2)) * f_carrier ** 2 + 2 * h_r * f_doppler * light_speed ** 2 * delay * (v_r[2] + 2 * v_t[2]) * f_carrier + f_doppler ** 2 * light_speed ** 2 * (delay ** 2 * light_speed ** 2 + h_r ** 2)))) / np.sin(elevation) / (f_carrier ** 2 * (v_t[1] ** 2 - v_t[2] ** 2) * np.cos(elevation) ** 4 + 2 * v_t[1] * f_carrier * (np.sin(elevation) * v_t[2] * f_carrier + f_doppler * light_speed) * np.cos(elevation) ** 3 + (2 * v_t[2] * np.sin(elevation) * f_carrier * f_doppler * light_speed + (-v_r[0] ** 2 - 2 * v_r[1] * v_t[1] + v_t[2] ** 2) * f_carrier ** 2 + f_doppler ** 2 * light_speed ** 2) * np.cos(elevation) ** 2 - 2 * v_r[1] * f_carrier * (np.sin(elevation) * v_t[2] * f_carrier + f_doppler * light_speed) * np.cos(elevation) + f_carrier ** 2 * (v_r[0] ** 2 + v_r[1] ** 2))
    #np.place(sol, sol.imag > 0, np.nan)
    return sol

def y_delay_doppler_2(delay, f_doppler, sim_config):
    v_t = sim_config.v_t
    v_r = sim_config.v_r
    elevation = sim_config.elevation
    f_carrier = sim_config.f_carrier
    h_r = sim_config.h_r

    delay = delay.astype(complex)
    f_doppler = f_doppler.astype(complex)

    sol =  (2 * np.cos(elevation) ** 4 * delay * v_t[1] * v_t[2] * f_carrier ** 2 * light_speed + f_carrier * (-light_speed * delay * f_carrier * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2]) * np.sin(elevation) + (v_t[2] ** 2 + v_r[2] * v_t[2] + v_t[1] * (v_r[1] - v_t[1])) * h_r * f_carrier + 2 * f_doppler * light_speed ** 2 * v_t[2] * delay) * np.cos(elevation) ** 3 + f_carrier * ((((v_r[1] - 2 * v_t[1]) * v_t[2] - v_r[2] * v_t[1]) * h_r * f_carrier - 2 * f_doppler * light_speed ** 2 * v_t[1] * delay) * np.sin(elevation) + light_speed * (-delay * f_carrier * (v_r[1] + 2 * v_t[1]) * v_t[2] + h_r * f_doppler * (v_r[1] - 2 * v_t[1]))) * np.cos(elevation) ** 2 + (-light_speed * (-delay * (v_r[0] ** 2 + v_r[1] * v_t[1] - v_t[2] ** 2) * f_carrier ** 2 + h_r * f_doppler * (v_r[2] + 2 * v_t[2]) * f_carrier + f_doppler ** 2 * light_speed ** 2 * delay) * np.sin(elevation) - (v_t[2] ** 2 + v_r[2] * v_t[2] + v_r[1] * (v_r[1] - v_t[1])) * h_r * f_carrier ** 2 - 2 * delay * v_t[2] * f_carrier * f_doppler * light_speed ** 2 - h_r * f_doppler ** 2 * light_speed ** 2) * np.cos(elevation) + (h_r * f_carrier * (v_r[2] + v_t[2]) + f_doppler * light_speed ** 2 * delay) * f_carrier * v_r[1] * np.sin(elevation) + delay * v_r[1] * v_t[2] * f_carrier ** 2 * light_speed + h_r * v_r[1] * f_carrier * f_doppler * light_speed - np.sqrt(-np.sin(elevation) ** 2 * v_r[0] ** 2 * f_carrier ** 2 * (-4 * np.cos(elevation) ** 3 * delay * h_r * v_t[1] * v_t[2] * f_carrier ** 2 * light_speed + f_carrier * (2 * h_r * light_speed * delay * f_carrier * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2]) * np.sin(elevation) + ((-v_t[2] + v_r[1] - v_t[1] - v_r[2]) * (v_t[2] + v_r[1] - v_t[1] + v_r[2]) * h_r ** 2 + light_speed ** 2 * delay ** 2 * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2])) * f_carrier - 4 * h_r * f_doppler * light_speed ** 2 * v_t[2] * delay) * np.cos(elevation) ** 2 - 2 * ((((v_r[2] + v_t[2]) * (v_r[1] - v_t[1]) * h_r ** 2 - light_speed ** 2 * v_t[1] * v_t[2] * delay ** 2) * f_carrier - 2 * h_r * f_doppler * light_speed ** 2 * v_t[1] * delay) * np.sin(elevation) + light_speed * (delay * (-2 * v_t[1] * v_t[2] + v_r[2] * (v_r[1] - v_t[1])) * h_r * f_carrier + ((v_r[1] - v_t[1]) * h_r ** 2 - light_speed ** 2 * v_t[1] * delay ** 2) * f_doppler)) * f_carrier * np.cos(elevation) + 2 * light_speed * (-h_r * delay * (v_r[0] ** 2 + v_r[1] ** 2 - v_r[2] * v_t[2] - v_t[2] ** 2) * f_carrier ** 2 + f_doppler * ((v_r[2] + v_t[2]) * h_r ** 2 + light_speed ** 2 * v_t[2] * delay ** 2) * f_carrier + h_r * f_doppler ** 2 * light_speed ** 2 * delay) * np.sin(elevation) + ((v_r[2] + v_t[2]) ** 2 * h_r ** 2 - light_speed ** 2 * delay ** 2 * (v_r[0] ** 2 + v_r[1] ** 2 - v_t[2] ** 2)) * f_carrier ** 2 + 2 * h_r * f_doppler * light_speed ** 2 * delay * (v_r[2] + 2 * v_t[2]) * f_carrier + f_doppler ** 2 * light_speed ** 2 * (delay ** 2 * light_speed ** 2 + h_r ** 2)))) / np.sin(elevation) / (f_carrier ** 2 * (v_t[1] ** 2 - v_t[2] ** 2) * np.cos(elevation) ** 4 + 2 * v_t[1] * f_carrier * (np.sin(elevation) * v_t[2] * f_carrier + f_doppler * light_speed) * np.cos(elevation) ** 3 + (2 * v_t[2] * np.sin(elevation) * f_carrier * f_doppler * light_speed + (-v_r[0] ** 2 - 2 * v_r[1] * v_t[1] + v_t[2] ** 2) * f_carrier ** 2 + f_doppler ** 2 * light_speed ** 2) * np.cos(elevation) ** 2 - 2 * v_r[1] * f_carrier * (np.sin(elevation) * v_t[2] * f_carrier + f_doppler * light_speed) * np.cos(elevation) + f_carrier ** 2 * (v_r[0] ** 2 + v_r[1] ** 2))
    #np.place(sol, sol.imag > 0, np.nan)
    return sol

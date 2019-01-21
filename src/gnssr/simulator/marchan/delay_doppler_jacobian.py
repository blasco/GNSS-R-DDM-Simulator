#!/usr/bin/env python

''' 
Facobians computed as defined in equations 11 and 12. 
x(delay,doppler) and y(delay,doppler) form the solutions to equations 17 and 18 
generated with Maple 

    J. F. Marchan-Hernandez, A. Camps, N. Rodriguez-Alvarez, E.  Valencia, X.  
    Bosch-Lluis, and I. Ramos-Perez, “An Efficient Algorithm to the Simulation 
    of Delay–Doppler Maps of Reflected Global Navigation Satellite System 
    Signals,” IEEE Transactions on Geoscience and Remote Sensing, vol. 47, no. 
    8, pp.  2733–2740, Aug. 2009.  
'''

import numpy as np
from problem_definition import *

def delay_doppler_jacobian_1(delay, f_doppler):
    try:
        j_11 = (1/delay_resolution)* \
                (x_delay_doppler_1(delay+delay_resolution/2, f_doppler) - \
                 x_delay_doppler_1(delay-delay_resolution/2, f_doppler))
        j_12 = (1/doppler_resolution)* \
                (x_delay_doppler_1(delay, f_doppler+doppler_resolution/2) - \
                 x_delay_doppler_1(delay, f_doppler-doppler_resolution/2))
        j_21 = (1/delay_resolution)* \
                (y_delay_doppler_1(delay+delay_resolution/2, f_doppler) - \
                 y_delay_doppler_1(delay-delay_resolution/2, f_doppler))
        j_22 = (1/doppler_resolution)* \
                (y_delay_doppler_1(delay, f_doppler+doppler_resolution/2) - \
                 y_delay_doppler_1(delay, f_doppler-doppler_resolution/2))
        return np.linalg.det(np.array([
                [j_11, j_12],
                [j_21, j_22]
                ], dtype=np.float64))

    except RuntimeWarning:
        return 0

def delay_doppler_jacobian_2(delay, f_doppler):
    try:
        j_11 = (1/delay_resolution)* \
                (x_delay_doppler_2(delay+delay_resolution/2, f_doppler) - \
                 x_delay_doppler_2(delay-delay_resolution/2, f_doppler))
        j_12 = (1/doppler_resolution)* \
                (x_delay_doppler_2(delay, f_doppler+doppler_resolution/2) - \
                 x_delay_doppler_2(delay, f_doppler-doppler_resolution/2))
        j_21 = (1/delay_resolution)* \
                (y_delay_doppler_2(delay+delay_resolution/2, f_doppler) - \
                 y_delay_doppler_2(delay-delay_resolution/2, f_doppler))
        j_22 = (1/doppler_resolution)* \
                (y_delay_doppler_2(delay, f_doppler+doppler_resolution/2) - \
                 y_delay_doppler_2(delay, f_doppler-doppler_resolution/2))
        return np.linalg.det(np.array([
                [j_11, j_12],
                [j_21, j_22]
                ], dtype=np.float64))

    except RuntimeWarning:
        return 0

def x_delay_doppler_1(delay, f_doppler): 
    sol =  ((np.cos(elevation) ** 2 * v_t[1] + (v_t[2] * np.sin(elevation) + f_doppler) * np.cos(elevation) - v_r[1]) * np.sqrt(-(-4 * np.cos(elevation) ** 3 * delay * h_r * v_t[1] * v_t[2] * light_speed + (2 * h_r * light_speed * delay * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2]) * np.sin(elevation) + (v_t[2] + v_r[1] - v_t[1] + v_r[2]) * (-v_t[2] + v_r[1] - v_t[1] - v_r[2]) * h_r ** 2 - 4 * delay * h_r * v_t[2] * f_doppler * light_speed + light_speed ** 2 * delay ** 2 * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2])) * np.cos(elevation) ** 2 + ((-2 * (v_r[2] + v_t[2]) * (v_r[1] - v_t[1]) * h_r ** 2 + 4 * delay * h_r * v_t[1] * f_doppler * light_speed + 2 * delay ** 2 * v_t[1] * v_t[2] * light_speed ** 2) * np.sin(elevation) - 2 * f_doppler * h_r ** 2 * (v_r[1] - v_t[1]) - 2 * light_speed * (-2 * v_t[2] * v_t[1] + v_r[2] * (v_r[1] - v_t[1])) * delay * h_r + 2 * f_doppler * light_speed ** 2 * v_t[1] * delay ** 2) * np.cos(elevation) + (2 * f_doppler * (v_r[2] + v_t[2]) * h_r ** 2 - 2 * light_speed * delay * (v_r[0] ** 2 + v_r[1] ** 2 - v_r[2] * v_t[2] - v_t[2] ** 2 - f_doppler ** 2) * h_r + 2 * delay ** 2 * v_t[2] * f_doppler * light_speed ** 2) * np.sin(elevation) + (v_r[2] ** 2 + 2 * v_r[2] * v_t[2] + v_t[2] ** 2 + f_doppler ** 2) * h_r ** 2 + 2 * light_speed * delay * f_doppler * (v_r[2] + 2 * v_t[2]) * h_r - light_speed ** 2 * delay ** 2 * (v_r[0] ** 2 + v_r[1] ** 2 - v_t[2] ** 2 - f_doppler ** 2)) * np.sin(elevation) ** 2 * v_r[0] ** 2) + (h_r * (v_r[1] - v_t[1]) * np.cos(elevation) ** 3 + (-h_r * (v_r[2] + v_t[2]) * np.sin(elevation) - delay * v_t[2] * light_speed - h_r * f_doppler) * np.cos(elevation) ** 2 - (v_r[1] - v_t[1]) * (delay * light_speed * np.sin(elevation) + h_r) * np.cos(elevation) + (light_speed * f_doppler * delay + v_r[2] * h_r + h_r * v_t[2]) * np.sin(elevation) + delay * v_t[2] * light_speed + h_r * f_doppler) * v_r[0] ** 2) / ((v_t[1] ** 2 - v_t[2] ** 2) * np.cos(elevation) ** 4 + 2 * v_t[1] * (v_t[2] * np.sin(elevation) + f_doppler) * np.cos(elevation) ** 3 + (2 * v_t[2] * np.sin(elevation) * f_doppler + f_doppler ** 2 - v_r[0] ** 2 - 2 * v_r[1] * v_t[1] + v_t[2] ** 2) * np.cos(elevation) ** 2 - 2 * v_r[1] * (v_t[2] * np.sin(elevation) + f_doppler) * np.cos(elevation) + v_r[0] ** 2 + v_r[1] ** 2) / np.sin(elevation) / v_r[0]
    if np.isnan(sol):
        return 0
    else:
        return sol

def x_delay_doppler_2(delay, f_doppler): 
    sol = ((-np.cos(elevation) ** 2 * v_t[1] + (-v_t[2] * np.sin(elevation) - f_doppler) * np.cos(elevation) + v_r[1]) * np.sqrt(-(-4 * np.cos(elevation) ** 3 * delay * h_r * v_t[1] * v_t[2] * light_speed + (2 * h_r * light_speed * delay * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2]) * np.sin(elevation) + (v_t[2] + v_r[1] - v_t[1] + v_r[2]) * (-v_t[2] + v_r[1] - v_t[1] - v_r[2]) * h_r ** 2 - 4 * delay * h_r * v_t[2] * f_doppler * light_speed + light_speed ** 2 * delay ** 2 * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2])) * np.cos(elevation) ** 2 + ((-2 * (v_r[2] + v_t[2]) * (v_r[1] - v_t[1]) * h_r ** 2 + 4 * delay * h_r * v_t[1] * f_doppler * light_speed + 2 * delay ** 2 * v_t[1] * v_t[2] * light_speed ** 2) * np.sin(elevation) - 2 * f_doppler * h_r ** 2 * (v_r[1] - v_t[1]) - 2 * light_speed * (-2 * v_t[2] * v_t[1] + v_r[2] * (v_r[1] - v_t[1])) * delay * h_r + 2 * f_doppler * light_speed ** 2 * v_t[1] * delay ** 2) * np.cos(elevation) + (2 * f_doppler * (v_r[2] + v_t[2]) * h_r ** 2 - 2 * light_speed * delay * (v_r[0] ** 2 + v_r[1] ** 2 - v_r[2] * v_t[2] - v_t[2] ** 2 - f_doppler ** 2) * h_r + 2 * delay ** 2 * v_t[2] * f_doppler * light_speed ** 2) * np.sin(elevation) + (v_r[2] ** 2 + 2 * v_r[2] * v_t[2] + v_t[2] ** 2 + f_doppler ** 2) * h_r ** 2 + 2 * light_speed * delay * f_doppler * (v_r[2] + 2 * v_t[2]) * h_r - light_speed ** 2 * delay ** 2 * (v_r[0] ** 2 + v_r[1] ** 2 - v_t[2] ** 2 - f_doppler ** 2)) * np.sin(elevation) ** 2 * v_r[0] ** 2) + (h_r * (v_r[1] - v_t[1]) * np.cos(elevation) ** 3 + (-h_r * (v_r[2] + v_t[2]) * np.sin(elevation) - delay * v_t[2] * light_speed - h_r * f_doppler) * np.cos(elevation) ** 2 - (v_r[1] - v_t[1]) * (delay * light_speed * np.sin(elevation) + h_r) * np.cos(elevation) + (light_speed * f_doppler * delay + v_r[2] * h_r + h_r * v_t[2]) * np.sin(elevation) + delay * v_t[2] * light_speed + h_r * f_doppler) * v_r[0] ** 2) / ((v_t[1] ** 2 - v_t[2] ** 2) * np.cos(elevation) ** 4 + 2 * v_t[1] * (v_t[2] * np.sin(elevation) + f_doppler) * np.cos(elevation) ** 3 + (2 * v_t[2] * np.sin(elevation) * f_doppler + f_doppler ** 2 - v_r[0] ** 2 - 2 * v_r[1] * v_t[1] + v_t[2] ** 2) * np.cos(elevation) ** 2 - 2 * v_r[1] * (v_t[2] * np.sin(elevation) + f_doppler) * np.cos(elevation) + v_r[0] ** 2 + v_r[1] ** 2) / np.sin(elevation) / v_r[0]
    if np.isnan(sol):
        return 0
    else:
        return sol

def y_delay_doppler_1(delay, f_doppler):
    sol = (2 * np.cos(elevation) ** 4 * delay * v_t[1] * v_t[2] * light_speed + (-light_speed * delay * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2]) * np.sin(elevation) + h_r * v_t[2] ** 2 + (2 * light_speed * f_doppler * delay + v_r[2] * h_r) * v_t[2] + h_r * v_t[1] * (v_r[1] - v_t[1])) * np.cos(elevation) ** 3 + (((v_r[1] - 2 * v_t[1]) * h_r * v_t[2] - v_t[1] * (2 * light_speed * f_doppler * delay + v_r[2] * h_r)) * np.sin(elevation) - light_speed * delay * (v_r[1] + 2 * v_t[1]) * v_t[2] + h_r * f_doppler * (v_r[1] - 2 * v_t[1])) * np.cos(elevation) ** 2 + ((-delay * v_t[2] ** 2 * light_speed - 2 * h_r * v_t[2] * f_doppler - h_r * v_r[2] * f_doppler + light_speed * delay * (v_r[0] ** 2 + v_r[1] * v_t[1] - f_doppler ** 2)) * np.sin(elevation) - h_r * v_t[2] ** 2 + (-2 * light_speed * f_doppler * delay - v_r[2] * h_r) * v_t[2] - h_r * (v_r[1] ** 2 - v_r[1] * v_t[1] + f_doppler ** 2)) * np.cos(elevation) + v_r[1] * (light_speed * f_doppler * delay + v_r[2] * h_r + h_r * v_t[2]) * np.sin(elevation) + delay * v_r[1] * v_t[2] * light_speed + h_r * v_r[1] * f_doppler + np.sqrt(-(-4 * np.cos(elevation) ** 3 * delay * h_r * v_t[1] * v_t[2] * light_speed + (2 * h_r * light_speed * delay * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2]) * np.sin(elevation) + (v_t[2] + v_r[1] - v_t[1] + v_r[2]) * (-v_t[2] + v_r[1] - v_t[1] - v_r[2]) * h_r ** 2 - 4 * delay * h_r * v_t[2] * f_doppler * light_speed + light_speed ** 2 * delay ** 2 * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2])) * np.cos(elevation) ** 2 + ((-2 * (v_r[2] + v_t[2]) * (v_r[1] - v_t[1]) * h_r ** 2 + 4 * delay * h_r * v_t[1] * f_doppler * light_speed + 2 * delay ** 2 * v_t[1] * v_t[2] * light_speed ** 2) * np.sin(elevation) - 2 * f_doppler * h_r ** 2 * (v_r[1] - v_t[1]) - 2 * light_speed * (-2 * v_t[2] * v_t[1] + v_r[2] * (v_r[1] - v_t[1])) * delay * h_r + 2 * f_doppler * light_speed ** 2 * v_t[1] * delay ** 2) * np.cos(elevation) + (2 * f_doppler * (v_r[2] + v_t[2]) * h_r ** 2 - 2 * light_speed * delay * (v_r[0] ** 2 + v_r[1] ** 2 - v_r[2] * v_t[2] - v_t[2] ** 2 - f_doppler ** 2) * h_r + 2 * delay ** 2 * v_t[2] * f_doppler * light_speed ** 2) * np.sin(elevation) + (v_r[2] ** 2 + 2 * v_r[2] * v_t[2] + v_t[2] ** 2 + f_doppler ** 2) * h_r ** 2 + 2 * light_speed * delay * f_doppler * (v_r[2] + 2 * v_t[2]) * h_r - light_speed ** 2 * delay ** 2 * (v_r[0] ** 2 + v_r[1] ** 2 - v_t[2] ** 2 - f_doppler ** 2)) * np.sin(elevation) ** 2 * v_r[0] ** 2)) / ((v_t[1] ** 2 - v_t[2] ** 2) * np.cos(elevation) ** 4 + 2 * v_t[1] * (v_t[2] * np.sin(elevation) + f_doppler) * np.cos(elevation) ** 3 + (2 * v_t[2] * np.sin(elevation) * f_doppler + f_doppler ** 2 - v_r[0] ** 2 - 2 * v_r[1] * v_t[1] + v_t[2] ** 2) * np.cos(elevation) ** 2 - 2 * v_r[1] * (v_t[2] * np.sin(elevation) + f_doppler) * np.cos(elevation) + v_r[0] ** 2 + v_r[1] ** 2) / np.sin(elevation)
    if np.isnan(sol):
        return 0
    else:
        return sol

def y_delay_doppler_2(delay, f_doppler):
    sol = (2 * np.cos(elevation) ** 4 * delay * v_t[1] * v_t[2] * light_speed + (-light_speed * delay * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2]) * np.sin(elevation) + h_r * v_t[2] ** 2 + (2 * light_speed * f_doppler * delay + v_r[2] * h_r) * v_t[2] + h_r * v_t[1] * (v_r[1] - v_t[1])) * np.cos(elevation) ** 3 + (((v_r[1] - 2 * v_t[1]) * h_r * v_t[2] - v_t[1] * (2 * light_speed * f_doppler * delay + v_r[2] * h_r)) * np.sin(elevation) - light_speed * delay * (v_r[1] + 2 * v_t[1]) * v_t[2] + h_r * f_doppler * (v_r[1] - 2 * v_t[1])) * np.cos(elevation) ** 2 + ((-delay * v_t[2] ** 2 * light_speed - 2 * h_r * v_t[2] * f_doppler - h_r * v_r[2] * f_doppler + light_speed * delay * (v_r[0] ** 2 + v_r[1] * v_t[1] - f_doppler ** 2)) * np.sin(elevation) - h_r * v_t[2] ** 2 + (-2 * light_speed * f_doppler * delay - v_r[2] * h_r) * v_t[2] - h_r * (v_r[1] ** 2 - v_r[1] * v_t[1] + f_doppler ** 2)) * np.cos(elevation) + v_r[1] * (light_speed * f_doppler * delay + v_r[2] * h_r + h_r * v_t[2]) * np.sin(elevation) + delay * v_r[1] * v_t[2] * light_speed + h_r * v_r[1] * f_doppler - np.sqrt(-(-4 * np.cos(elevation) ** 3 * delay * h_r * v_t[1] * v_t[2] * light_speed + (2 * h_r * light_speed * delay * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2]) * np.sin(elevation) + (v_t[2] + v_r[1] - v_t[1] + v_r[2]) * (-v_t[2] + v_r[1] - v_t[1] - v_r[2]) * h_r ** 2 - 4 * delay * h_r * v_t[2] * f_doppler * light_speed + light_speed ** 2 * delay ** 2 * (v_t[1] - v_t[2]) * (v_t[1] + v_t[2])) * np.cos(elevation) ** 2 + ((-2 * (v_r[2] + v_t[2]) * (v_r[1] - v_t[1]) * h_r ** 2 + 4 * delay * h_r * v_t[1] * f_doppler * light_speed + 2 * delay ** 2 * v_t[1] * v_t[2] * light_speed ** 2) * np.sin(elevation) - 2 * f_doppler * h_r ** 2 * (v_r[1] - v_t[1]) - 2 * light_speed * (-2 * v_t[2] * v_t[1] + v_r[2] * (v_r[1] - v_t[1])) * delay * h_r + 2 * f_doppler * light_speed ** 2 * v_t[1] * delay ** 2) * np.cos(elevation) + (2 * f_doppler * (v_r[2] + v_t[2]) * h_r ** 2 - 2 * light_speed * delay * (v_r[0] ** 2 + v_r[1] ** 2 - v_r[2] * v_t[2] - v_t[2] ** 2 - f_doppler ** 2) * h_r + 2 * delay ** 2 * v_t[2] * f_doppler * light_speed ** 2) * np.sin(elevation) + (v_r[2] ** 2 + 2 * v_r[2] * v_t[2] + v_t[2] ** 2 + f_doppler ** 2) * h_r ** 2 + 2 * light_speed * delay * f_doppler * (v_r[2] + 2 * v_t[2]) * h_r - light_speed ** 2 * delay ** 2 * (v_r[0] ** 2 + v_r[1] ** 2 - v_t[2] ** 2 - f_doppler ** 2)) * np.sin(elevation) ** 2 * v_r[0] ** 2)) / ((v_t[1] ** 2 - v_t[2] ** 2) * np.cos(elevation) ** 4 + 2 * v_t[1] * (v_t[2] * np.sin(elevation) + f_doppler) * np.cos(elevation) ** 3 + (2 * v_t[2] * np.sin(elevation) * f_doppler + f_doppler ** 2 - v_r[0] ** 2 - 2 * v_r[1] * v_t[1] + v_t[2] ** 2) * np.cos(elevation) ** 2 - 2 * v_r[1] * (v_t[2] * np.sin(elevation) + f_doppler) * np.cos(elevation) + v_r[0] ** 2 + v_r[1] ** 2) / np.sin(elevation)
    if np.isnan(sol):
        return 0
    else:
        return sol

print("d=1,f=500: {0}".format(format(x_delay_doppler_1(2/1.023e6, -2500),".8E")))
print("d=1,f=500: {0}".format(format(y_delay_doppler_1(2/1.023e6, -2500),".8E")))

print("d=1,f=500: {0}".format(format(x_delay_doppler_2(2/1.023e6, -2500),".8E")))
print("d=1,f=500: {0}".format(format(y_delay_doppler_2(2/1.023e6, -2500),".8E")))

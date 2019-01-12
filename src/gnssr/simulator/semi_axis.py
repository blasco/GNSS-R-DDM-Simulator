#!/usr/bin/env python

import numpy as np

def calculate_ellipse_semi_axes(h_t, h_r, elev_deg, T_c):
    elev = elev_deg*np.pi/180
    c = 299792458 # m/s
    rs = h_r/np.sin(elev)
    ts = h_t/np.sin(elev)
    a = (1/np.sin(elev))*((rs*ts*T_c*c)/(rs+ts))**(1/2)
    b = ((rs*ts*T_c*c)/(rs+ts))**(1/2)
    return a, b

# Altitude
h_t = 20000e3 # meters
h_r = 20e3 # meters
# Elevation
elev_deg = 42.2
# Delay
T_c = 2.44e-7*15 # s

print(calculate_ellipse_semi_axes(h_t, h_r, elev_deg, T_c))

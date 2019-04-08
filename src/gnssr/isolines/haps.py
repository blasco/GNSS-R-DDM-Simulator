#!/usr/bin/env python

from iso_lines import *
import matplotlib.pyplot as plt


haps_iso_lines = iso_lines()

# --- Parameters ---

# Altitude
haps_iso_lines.h_t = 20000e3 # meters
haps_iso_lines.h_r = 20e3 # meters

# Elevation
haps_iso_lines.elev_deg = 60 # deg

# Velocity
haps_iso_lines.v_tx = 2121 # m/s
haps_iso_lines.v_ty = 2121 # m/s
haps_iso_lines.v_tz = 5 # m/s

haps_iso_lines.v_rx = 25 # m/s
haps_iso_lines.v_ry = 25 # m/s
haps_iso_lines.v_rz = 25 # m/s

# Plotting Area
haps_iso_lines.extent_x0 =  -10e3 # meters
haps_iso_lines.extent_x1 =  10e3 # meters
haps_iso_lines.extent_y0 =  -10e3 # meters
haps_iso_lines.extent_y1 =  10e3 # meters

haps_iso_lines.linsapce_delta = 500

# -----------------

# Iso Delay values 
delay_start = 0 # C/A chips
delay_increment = 0.5# C/A chips
delay_end = 3 # C/A chips
delay_range = [delay_start, delay_increment, delay_end]

# Iso Doppler values
doppler_start = -70 # Hz
doppler_increment = 10 # Hz
doppler_end = 70 # Hz
doppler_range = [doppler_start, doppler_increment, doppler_end]

haps_iso_lines.plot_range(delay_range, doppler_range)
print(haps_iso_lines.calculate_ellipse_semi_axis())
plt.show()

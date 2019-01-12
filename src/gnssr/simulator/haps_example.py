#!/usr/bin/env python

from iso_lines import *
import matplotlib.pyplot as plt

haps_iso_lines = iso_lines()

# --- Parameters ---

# Altitude
haps_iso_lines.h_t = 20000e3 # meters
haps_iso_lines.h_r = 20e3 # meters

# Elevation
haps_iso_lines.elev_deg = 52.2 # deg

# Velocity
haps_iso_lines.v_tx = 2121 # m/s
haps_iso_lines.v_ty = 2121 # m/s
haps_iso_lines.v_tz = 5 # m/s

haps_iso_lines.v_rx = 5 # m/s
haps_iso_lines.v_ry = 5 # m/s
haps_iso_lines.v_rz = 5 # m/s

# Plotting Area
haps_iso_lines.extent_x0 =  -20e3 # meters
haps_iso_lines.extent_x1 =  20e3 # meters
haps_iso_lines.extent_y0 =  -20e3 # meters
haps_iso_lines.extent_y1 =  20e3 # meters

haps_iso_lines.linsapce_delta = 500

# -----------------

# Iso Delay values 
delay_start = 0 # micro secs
delay_increment = 2.44e-7*1e6 # micro secs
delay_end = 2.44e-7*1e6*15 # micro secs
delay_range = [delay_start, delay_increment, delay_end]

# Iso Doppler values
doppler_start = -3000 # Hz
doppler_increment = 500 # Hz
doppler_end = 3000 # Hz
doppler_range = [doppler_start, doppler_increment, doppler_end]


haps_iso_lines.plot_range(delay_range, doppler_range)
print(haps_iso_lines.calculate_ellipse_semi_axis())
plt.show()

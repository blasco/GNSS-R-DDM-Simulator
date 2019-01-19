#!/usr/bin/env python

from iso_lines import *
import matplotlib.pyplot as plt

tds_iso_lines = iso_lines()

# --- Parameters ---

# Altitude
tds_iso_lines.h_t = 20000e3 # meters
tds_iso_lines.h_r = 500e3 # meters

# Elevation
tds_iso_lines.elev_deg = 72.2

# Velocity
tds_iso_lines.v_tx = 2121
tds_iso_lines.v_ty = 2121
tds_iso_lines.v_tz = 5

tds_iso_lines.v_rx = 2210
tds_iso_lines.v_ry = 7299
tds_iso_lines.v_rz = 199

# Plotting Area
tds_iso_lines.extent_x0 =  -120e3 # meters
tds_iso_lines.extent_x1 =  120e3 # meters
tds_iso_lines.extent_y0 =  -120e3 # meters
tds_iso_lines.extent_y1 =  120e3 # meters

tds_iso_lines.linsapce_delta = 500

# -----------------

# Iso Delay values 
delay_start = 0 # C/A chips
delay_increment = 0.5 # C/A chips
delay_end = 15 # C/A chips
delay_range = [delay_start, delay_increment, delay_end]

# Iso Doppler values
doppler_start = -3000 # micro secs
doppler_increment = 500 # micro secs
doppler_end = 3000 # micro secs
doppler_range = [doppler_start, doppler_increment, doppler_end]

tds_iso_lines.plot_range(delay_range, doppler_range)
plt.show()

#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from gnssr.utils import *

class iso_lines:

    def __init__(self):
        # Altitude
        self.h_t = 20000e3 # meters
        self.h_r = 500e3 # meters

        # Elevation
        self.elev_deg = 72.2

        # Velocity
        self.v_tx = 2121
        self.v_ty = 2121
        self.v_tz = 5

        self.v_rx = 2210
        self.v_ry = 7299
        self.v_rz = 199

        # Plotting Area
        self.extent_x0 =  -20e3 # meters
        self.extent_x1 =  20e3 # meters
        self.extent_y0 =  -20e3 # meters
        self.extent_y1 =  20e3 # meters
        self.linsapce_delta = 500

    def z_eq(self, x, y):
        #R = 6373310.2078
        #return (R**2 -x**2 -y**2 )**(1/2) - R
        return 0

    def time_eq(self, x, y):
        self.elev = self.elev_deg*np.pi/180
        h_t = self.h_t 
        h_r = self.h_r
        elev = self.elev
        z_eq = self.z_eq
        c = 299792458 # m/s
        return (1/c)*( \
               (x**2 + (y-h_t/np.tan(elev))**2 + (z_eq(x,y)-h_t)**2)**(1/2) + \
               (x**2 + (-h_r/np.tan(elev) -y)**2 + (h_r - z_eq(x,y))**2)**(1/2) \
               )

    def calculate_ellipse_semi_axis(self):
        c = 299792458 # m/s
        T_c = 2.44e-7*15 # s
        rs = self.h_r/np.sin(self.elev)
        ts = self.h_t/np.sin(self.elev)
        a = (1/np.sin(self.elev))*((rs*ts*T_c*c)/(rs+ts))**(1/2)
        b = ((rs*ts*T_c*c)/(rs+ts))**(1/2)
        return a, b

    def time_inc_eq(self, x, y):
        return self.time_eq(x,y) - self.time_eq(0,0)

    def time_inc_eq_usec(self, x, y):
        return (self.time_eq(x,y) - self.time_eq(0,0))*1e6

    def time_inc_eq_chips(self, x, y):
        return (self.time_eq(x,y) - self.time_eq(0,0))*chips_per_second

    def doppler_eq(self, x, y):
        self.elev = self.elev_deg*np.pi/180
        h_t = self.h_t 
        h_r = self.h_r
        elev = self.elev
        z_eq = self.z_eq
        v_tx = self.v_tx
        v_ty = self.v_ty
        v_tz = self.v_tz
        v_rx = self.v_rx
        v_ry = self.v_ry
        v_rz = self.v_rz
        # GPS L1 center frequency
        c = 299792458 # m/s
        #f_c = 1575.42e6 # Hz 
        f_0 = 10.23e6 # Hz
        f_c = 154*f_0;
        return (f_c/c)*( \
                (v_tx*(x)  + v_ty*(y-h_t/np.tan(elev)) + v_tz*(z_eq(x,y)-h_t))  / (x**2 + (y-h_t/np.tan(elev))**2 + (z_eq(x,y)-h_t)**2)**(1/2) \
               -(v_rx*(-x) + v_ry*(-h_r/np.tan(elev)-y)  + v_rz*(h_r -z_eq(x,y))   ) / (x**2 + (-h_r/np.tan(elev) -y)**2 + (h_r - z_eq(x,y))**2)**(1/2) \
                )

    def doppler_inc_eq(self, x, y):
        return self.doppler_eq(x,y) - self.doppler_eq(0,0)

    def prepare_plot(self):
        # Grid
        self.X, self.Y = np.meshgrid(
                np.linspace(self.extent_x0, self.extent_x1, self.linsapce_delta), 
                np.linspace(self.extent_y0, self.extent_y1, self.linsapce_delta)
                )
        self.Z_time = self.time_inc_eq_chips(self.X,self.Y)
        self.Z_doppler = self.doppler_inc_eq(self.X,self.Y)

        self.fig_lines, self.ax_lines = plt.subplots(1,figsize=(10, 4))
        self.ax_lines.set_title('Iso-Delay and Iso-Doppler Lines')

        ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/1000))
        ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000))
        self.ax_lines.xaxis.set_major_formatter(ticks_x)
        self.ax_lines.yaxis.set_major_formatter(ticks_y)

        plt.xlabel('[km]')
        plt.ylabel('[km]')

    def plot(self):
        self.prepare_plot()
        # Iso-Delay
        c_time = self.ax_lines.contour(self.X, self.Y, self.Z_time, cmap='winter')  # Use this if no countours where found 
        self.fig_lines.colorbar(c_time, label='C/A chips', )
        # Iso-Doppler
        c_doppler = self.ax_lines.contour(self.X, self.Y, self.Z_doppler, cmap='autumn') # Use this if no countours where found 
        self.fig_lines.colorbar(c_doppler, label='Hz')

        plt.show(block=False)

    def plot_range(self, delay_range, doppler_range):
        delay_start, delay_increment, delay_end = delay_range
        doppler_start, doppler_increment, doppler_end = doppler_range
        self.prepare_plot()
        # Iso-Delay
        iso_delay_values = list(np.arange(delay_start, delay_end, delay_increment))
        c_time = self.ax_lines.contour(self.X, self.Y, self.Z_time, iso_delay_values, cmap='winter')
        self.fig_lines.colorbar(c_time, label='C/A chips', )
        # Iso-Doppler
        iso_doppler_values = list(np.arange(doppler_start, doppler_end, doppler_increment))
        c_doppler = self.ax_lines.contour(self.X, self.Y, self.Z_doppler, iso_doppler_values, cmap='autumn')
        self.fig_lines.colorbar(c_doppler, label='Hz')

        plt.show(block=False)

# Usage example
def main():
    # Iso Delay values 
    delay_start = 0 # micro secs
    delay_increment = 2.44e-7*1e6 # micro secs
    delay_end = 2.44e-7*1e6*15 # micro secs
    delay_range = [delay_start, delay_increment, delay_end]

    # Iso Doppler values
    doppler_start = -3000 # micro secs
    doppler_increment = 500 # micro secs
    doppler_end = 3000 # micro secs
    doppler_range = [doppler_start, doppler_increment, doppler_end]

    haps_iso = iso_lines()
    haps_iso.plot_range(delay_range, doppler_range)

if __name__ == '__main__':
    main()

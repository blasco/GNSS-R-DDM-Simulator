#!/usr/bin/env python

import cv2
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from gnssr.simulator.waf import *
from gnssr.simulator.simulation_configuration import *
from gnssr.utils import gps_ca_chips_per_second

def main():

    sim_config = simulation_configuration()
    sim_config.doppler_resolution = 20 # Hz
    sim_config.delay_resolution = 0.1/gps_ca_chips_per_second # seconds

    sim_config.delay_increment_start = -3/gps_ca_chips_per_second
    sim_config.delay_increment_end = 10/gps_ca_chips_per_second

    sim_config.doppler_increment_start = -3000
    sim_config.doppler_increment_end = 3000

    delay_increment_start = sim_config.delay_increment_start 
    delay_increment_end = sim_config.delay_increment_end 
    delay_resolution = sim_config.delay_resolution

    doppler_increment_start = sim_config.doppler_increment_start
    doppler_increment_end = sim_config.doppler_increment_end
    doppler_resolution = sim_config.doppler_resolution

    delay_chip = sim_config.delay_chip

    # WAF
    waf_delay_increment_values = list(np.arange(
        -delay_increment_end, 
        delay_increment_end, 
        delay_resolution
        ))
    waf_doppler_increment_values = list(np.arange(
        doppler_increment_start, 
        doppler_increment_end, 
        doppler_resolution
        ))
    waf_delay_grid, waf_doppler_grid = np.meshgrid(waf_delay_increment_values, waf_doppler_increment_values)

    waf_matrix = woodward_ambiguity_function(waf_delay_grid, waf_doppler_grid, sim_config)**2

    fig_waf, ax_waf = plt.subplots(1,figsize=(10, 4))
    ax_waf.set_title('WAF')
    im = ax_waf.imshow(waf_matrix, cmap='jet', 
            extent=(-delay_increment_end/delay_chip, delay_increment_end/delay_chip, doppler_increment_start, doppler_increment_end),
            aspect="auto"
            )
    cbar = fig_waf.colorbar(im)

    fig_waf_delay, ax_waf_delay = plt.subplots(1,figsize=(10, 4))
    waf_delay_result = waf_delay(np.array(waf_delay_increment_values), sim_config)**2
    ax_waf_delay.plot([i/delay_chip for i in waf_delay_increment_values], waf_delay_result)
    ax_waf_delay.fill_between([i/delay_chip for i in waf_delay_increment_values], 0, waf_delay_result)
    ax_waf_delay.set_title('waf_delay')

    fig_waf_frequency, ax_waf_frequency = plt.subplots(1,figsize=(10, 4))
    waf_frequency_result = waf_frequency(np.array(waf_doppler_increment_values), sim_config)**2
    ax_waf_frequency.plot(waf_doppler_increment_values, waf_frequency_result)
    ax_waf_frequency.fill_between(waf_doppler_increment_values, 0, waf_frequency_result, cmap='jet')
    ax_waf_frequency.set_title('waf_freq')

    fig_waf_frequency_2, ax_waf_frequency_2 = plt.subplots(1,figsize=(10, 4))
    xx=np.array(waf_delay_increment_values)/delay_chip
    yy=waf_delay_result

    path = Path(np.array([xx,yy]).transpose())
    patch = PathPatch(path, facecolor='none')
    plt.gca().add_patch(patch)

    im = plt.imshow(xx.reshape(yy.size,1),  cmap='jet',
            origin='lower',extent=[-delay_increment_end/delay_chip, delay_increment_end/delay_chip,0,1.05],aspect="auto", clip_path=patch, clip_on=True)
    plt.grid(linestyle='dotted')

    fig_waf_frequency_2, ax_waf_frequency_2 = plt.subplots(1,figsize=(10, 4))
    xx_doppler=np.array(waf_doppler_increment_values)
    yy_doppler=waf_frequency_result

    path = Path(np.array([xx_doppler,yy_doppler]).transpose())
    patch = PathPatch(path, facecolor='none')
    plt.gca().add_patch(patch)

    im = plt.imshow(xx_doppler.reshape(yy_doppler.size,1),  cmap='jet',
            origin='lower',extent=[-3000,3000,0,1.05],aspect="auto", clip_path=patch, clip_on=True)
    plt.grid(linestyle='dotted')

    plt.show()

if __name__ == '__main__':
    main()

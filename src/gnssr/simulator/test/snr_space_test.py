#!/usr/bin/env python

import numpy as np

import matplotlib.pyplot as plt

from gnssr.simulator.ddm import *
from gnssr.simulator.simulation_configuration import *

from gnssr.targets import *
from gnssr.tds.tds_data import *
from gnssr.tds.detection.find_targets import *
from gnssr.utils import *

import cv2

def main():

    sim_config = simulation_configuration()
    sim_config.doppler_resolution = 50

    delay_chip = sim_config.delay_chip
    # Really good
    #sim_config.u_10 = 4
    #sim_config.phi_0 = 35*np.pi/180
    #sim_config.u_10 = 4
    #sim_config.phi_0 = 38*np.pi/180
    sim_config.u_10 = 2
    sim_config.phi_0 = -83*np.pi/180

    file_root_name = 'raw/L1B/2015-04-01-H00'
    target = targets['hibernia']
    group = '000095'
    index = 415

    tds = tds_data(file_root_name, group, index)

    # Load TDS Geometry in simulation configuration

    number_of_delay_pixels = tds.metagrp.groups[tds.group].NumberOfDelayPixels
    number_of_doppler_pixels = tds.metagrp.groups[tds.group].NumberOfDopplerPixels

    delay_start = tds.calculate_delay_increment_chips(0)
    delay_end = tds.calculate_delay_increment_chips(number_of_delay_pixels-1)
    doppler_start = tds.calculate_doppler_increment(-np.floor(number_of_doppler_pixels/2))
    doppler_end = tds.calculate_doppler_increment(np.floor(number_of_doppler_pixels/2 - 0.5))

    delay_increment_start = sim_config.delay_increment_start 
    delay_increment_end = sim_config.delay_increment_end 
    delay_resolution = sim_config.delay_resolution

    sim_config.delay_increment_start = delay_start*delay_chip + 2*delay_resolution
    sim_config.delay_increment_end = delay_end*delay_chip
    sim_config.doppler_increment_start = -5000
    sim_config.doppler_increment_end = 5000

    doppler_increment_start = sim_config.doppler_increment_start
    doppler_increment_end = sim_config.doppler_increment_end
    doppler_resolution = sim_config.doppler_resolution

    r_sp, lat_sp, lon_sp = tds.find_sp();
    n_z = unit_vector(ellip_norm(r_sp))
    n_x = unit_vector(np.cross(n_z, r_sp-tds.r_t))
    n_y = unit_vector(np.cross(n_z, n_x))

    v_tx = np.dot(tds.v_t, n_x)
    v_ty = np.dot(tds.v_t, n_y)
    v_tz = np.dot(tds.v_t, n_z)

    v_rx = np.dot(tds.v_r, n_x)
    v_ry = np.dot(tds.v_r, n_y)
    v_rz = np.dot(tds.v_r, n_z)

    #sim_config.set_scenario_local_ref(
    #    h_r = np.dot((tds.r_r - r_sp), n_z),
    #    h_t = np.dot((tds.r_t - r_sp), n_z),
    #    elevation = angle_between(n_y, tds.r_t-r_sp),
    #    v_t = np.array([v_tx,v_ty,v_tz]),
    #    v_r = np.array([v_rx,v_ry,v_rz])
    #    )

    # DDM
    ddm_sim = simulate_ddm(sim_config)

    fig_ddm, ax_ddm = plt.subplots(1,figsize=(10, 4))
    plt.title('DDM simulation')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    contour = ax_ddm.imshow(ddm_sim, cmap='jet', 
            extent=(delay_start, delay_end, doppler_end, doppler_start), 
            aspect="auto"
            )
    cbar = fig_ddm.colorbar(contour, label='Power [Watt]')

    T_noise_receiver = 225
    k_b = 1.38e-23 # J/K
    y_noise = sim_config.coherent_integration_time*k_b*T_noise_receiver
    ddm_sim_snr = np.copy(10*np.log10(np.abs(ddm_sim)/y_noise))

    fig_ddm_snr, ax_ddm_snr = plt.subplots(1,figsize=(10, 4))
    plt.title('DDM SNR simulation')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    contour_snr = ax_ddm_snr.imshow(ddm_sim_snr, cmap='jet', 
            extent=(delay_start, delay_end, doppler_end, doppler_start), 
            aspect="auto"
            )
    cbar = fig_ddm_snr.colorbar(contour_snr, label='SNR [dB]')

    plt.show()

if __name__ == '__main__':
    main()

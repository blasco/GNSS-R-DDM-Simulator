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
    delay_chip = sim_config.delay_chip
    # Really good
    #sim_config.u_10 = 4
    #sim_config.phi_0 = 35*np.pi/180
    #sim_config.u_10 = 4
    #sim_config.phi_0 = 38*np.pi/180
    sim_config.u_10 = 3.9
    sim_config.phi_0 = 37*np.pi/180

    file_root_name = 'raw/L1B/2015-04-01-H00'
    target = targets['hibernia']
    group = '000095'
    index = 515

    tds = tds_data(file_root_name, group, index)

    # Sea clutter estimation. 
    p = target_processor();
    for i in range(index - 30, index + 5):
        ddm = tds.rootgrp.groups[group].variables['DDM'][i].data
        p.process_ddm(ddm)

    # Plot TDS DDM sea clutter

    datenum = tds.rootgrp.groups[tds.group].variables['IntegrationMidPointTime'][tds.index]

    string = str(datenum_to_pytime(float(datenum))) \
        + ' Lat: ' + "{0:.2f}".format(tds.lat_sp_tds) \
        + ' Lon: ' + "{0:.2f}".format(tds.lon_sp_tds)

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

    fig, ax = plt.subplots(1,figsize=(10, 4))
    ax.set_ylabel('Hz')
    ax.set_xlabel('C/A chips')
    im = ax.imshow(p.sea_clutter, cmap='jet', 
            extent=(delay_start, delay_end, doppler_end, doppler_start), 
            aspect=(number_of_doppler_pixels/number_of_delay_pixels)/np.abs(doppler_start/delay_start)
            )
    t = plt.text(0.01, 0.80, string, {'color': 'w', 'fontsize': 12}, transform=ax.transAxes)

    # Load TDS Geometry in simulation configuration
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

    sim_config.set_scenario_local_ref(
        h_r = np.dot((tds.r_r - r_sp), n_z),
        h_t = np.dot((tds.r_t - r_sp), n_z),
        elevation = angle_between(n_y, tds.r_t-r_sp),
        v_t = np.array([v_tx,v_ty,v_tz]),
        v_r = np.array([v_rx,v_ry,v_rz])
        )

    # DDM
    ddm_sim = normalize(simulate_ddm(sim_config)) 

    fig_ddm, ax_ddm = plt.subplots(1,figsize=(10, 4))
    ax_ddm.set_title('DDM simulator')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    im = ax_ddm.imshow(ddm_sim, cmap='jet', 
            extent=(delay_start, delay_end, doppler_end, doppler_start), 
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
    rescaled_doppler_resolution = tds.doppler_resolution
    rescaled_delay_resolution_chips = tds.time_delay_resolution/delay_chip
    ddm_rescaled = cv2.resize(ddm_sim, 
            dsize=(
                128, 
                20
                ), 
            interpolation=cv2.INTER_AREA
            ) 
    ddm_rescaled = ddm_rescaled + 0.02*np.max(ddm_rescaled)*np.random.rand(ddm_rescaled.shape[0],ddm_rescaled.shape[1])

    im = ax_ddm_rescaled.imshow(ddm_rescaled, cmap='jet', 
            extent=(delay_start, delay_end, doppler_end, doppler_start), 
            aspect=(number_of_doppler_pixels/number_of_delay_pixels)/np.abs(doppler_start/delay_start)
            )


    fig_diff, ax_diff = plt.subplots(1,figsize=(10, 4))
    ax_diff.set_title('DDM diff')
    plt.xlabel('C/A chips')
    plt.ylabel('Hz')
    ddm_diff = np.abs(ddm_rescaled-p.sea_clutter)
    im = ax_diff.imshow(ddm_diff, cmap='jet', 
            extent=(delay_start, delay_end, doppler_end, doppler_start), 
            aspect=(number_of_doppler_pixels/number_of_delay_pixels)/np.abs(doppler_start/delay_start)
            )

    plt.show()

if __name__ == '__main__':
    main()

#!/usr/bin/env python

from netCDF4 import Dataset
import os
import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
import matplotlib.patches as mpatches
from datetime import *

################################################# 
# Show plot with time and frequency delays
################################################# 

file_root_name = 'raw/L1B/2018-07-29-H18'
rootgrp = Dataset(os.path.join(os.environ['TDS_ROOT'], file_root_name+'-DDMs.nc') , "r", format="NETCDF4")
metagrp = Dataset(os.path.join(os.environ['TDS_ROOT'], file_root_name+'-metadata.nc') , "r", format="NETCDF4")
group = '000050'
index = 388

a = 6378137 #m
b = 6356752.314245 #m

def find_index_meta(group, index):
    datenum = rootgrp.groups[group].variables['IntegrationMidPointTime'][index]
    datenums_meta = np.array(metagrp.variables['IntegrationMidPointTime'])
    for i, datenum_meta in enumerate(datenums_meta):
        if datenum_meta == datenum:
            #print("Foud: {0} | {1}".format(i, datenum_meta))
            return i
            break

index_meta = find_index_meta(group,index)

def calculate_delay_increment_seconds(delay_pixel, group, index):
    delay_in_samples = delay_pixel * metagrp.groups[group].CodeDelaySpacingSamplesBetweenPixels
    delay_in_seconds = delay_in_samples / metagrp.groups[group].SamplingFrequency
    delay_of_tracking_point_in_seconds = (metagrp.groups[group].TrackingOffsetDelayNs - metagrp.groups[group].variables['SpecularPathRangeOffset'][index])*1e-9
    delay_increment_seconds = delay_in_seconds - delay_of_tracking_point_in_seconds
    return delay_increment_seconds

def calculate_delay_increment_chips(delay_pixel, group, index):
    chips_per_second = 1.023e6
    return calculate_delay_increment_seconds(delay_pixel, group, index)*chips_per_second

def calculate_doppler_increment(doppler_pixel, group, index):
    return doppler_pixel*metagrp.groups[group].DopplerResolution - metagrp.groups[group].TrackingOffsetDopplerHz

#x_col = list(map(lambda x: calculate_delay_increment_chips(x,group,index), np.arange(0,128)))
x_col = list(map(lambda x: calculate_delay_increment_seconds(x,group,index), np.arange(0,128)))
print(x_col[70])
y_col = list(map(lambda x: calculate_doppler_increment(x,group,index), np.arange(-10,10)))
x, y = np.meshgrid(x_col, y_col)

ddm = rootgrp.groups[group].variables['DDM'][index].data

datenum = rootgrp.groups[group].variables['IntegrationMidPointTime'][index]
datenum_meta = metagrp.variables['IntegrationMidPointTime'][index_meta]
print("datenum group: {0}\ndatenum meta: {1}".format(datenum, datenum_meta))
lat = metagrp.groups[group].variables['SpecularPointLat'][index]
lon = metagrp.groups[group].variables['SpecularPointLon'][index]

fig, ax = plt.subplots(1,figsize=(10, 4))
cs = ax.contourf(x, y, ddm, 30, cmap='viridis')
ax.set_title('DDM')
ax.grid(c='k', ls='-', alpha=0.3)

plt.show(block=False)

# https://en.wikipedia.org/wiki/IERS_Reference_Meridian: 
# It is also the reference meridian of the Global Positioning System (GPS) operated by the United States Department of Defense, and of WGS84 

w_earth = np.array([0, 0, 7.2921158553e-5]) # rad/sec

r_tx = metagrp.groups[group].variables['TransmitterPositionX'][index]
r_ty = metagrp.groups[group].variables['TransmitterPositionY'][index]
r_tz = metagrp.groups[group].variables['TransmitterPositionZ'][index]
r_t = np.array([r_tx,r_ty,r_tz])

v_tx = metagrp.groups[group].variables['TransmitterVelocityX'][index]
v_ty = metagrp.groups[group].variables['TransmitterVelocityY'][index]
v_tz = metagrp.groups[group].variables['TransmitterVelocityZ'][index]
v_t = np.array([v_tx,v_ty,v_tz])
# inertial velocity
v_t_i = v_t + np.cross(w_earth, r_t)

from datetime import *
from dateutil import tz

def datenum_to_pytime(matlab_datenum):
    python_datetime = datetime.fromordinal(int(matlab_datenum)) + timedelta(days=matlab_datenum%1) - timedelta(days = 366)
    return python_datetime.strftime('%Y-%m-%d %H:%M:%S')

# Make sure that they correspond to the same time
#datenum = rootgrp.groups[group].variables['IntegrationMidPointTime'][index]
#print(datenum_to_pytime(float(datenum)))
#datenum_meta = metagrp.variables['IntegrationMidPointTime'][index_meta]
#print(datenum_to_pytime(float(datenum_meta)))

r_rx = metagrp.variables['ReceiverPositionX'][index_meta]
r_ry = metagrp.variables['ReceiverPositionY'][index_meta]
r_rz = metagrp.variables['ReceiverPositionZ'][index_meta]
r_r = np.array([r_rx,r_ry,r_rz])

v_rx = metagrp.variables['ReceiverVelocityX'][index_meta]
v_ry = metagrp.variables['ReceiverVelocityY'][index_meta]
v_rz = metagrp.variables['ReceiverVelocityZ'][index_meta]
v_r = np.array([v_rx, v_ry, v_rz])
# inertial velocity
v_r_i = v_r + np.cross(w_earth, r_r)

r_spx = metagrp.groups[group].variables['SpecularPointPositionX'][index]
r_spy = metagrp.groups[group].variables['SpecularPointPositionY'][index]
r_spz = metagrp.groups[group].variables['SpecularPointPositionZ'][index]
r_sp_tds = np.array([r_spx,r_spy,r_spz])

i_sp = metagrp.groups[group].variables['SPIncidenceAngle'][index]

################################################# 
# Find specular point
################################################# 

def find_sp(r_sp_estimate):

    def ellip_norm(r):
        x,y,z = r
        return np.array([2*x/a**2,2*y/a**2,2*z/b**2])
        
    r_center = np.array(r_sp_estimate)
    for it in range(3):
        r_inc = 100e3/10**(it)
        r_step = 1e3/10**(it)

        x = np.arange(r_center[0]-r_inc, r_center[0]+r_inc, r_step)
        y = np.arange(r_center[1]-r_inc, r_center[1]+r_inc, r_step)
        xx, yy = np.meshgrid(x, y)
        z = b*(1-(xx/a)**2-(yy/a)**2)**(1/2)

        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #ax.set_aspect('equal')
        #ax.plot_surface(xx,yy,z)
        #set_axes_equal(ax)

        min = np.inf
        for x_i, x_val in enumerate(x):
            for y_i, y_val in enumerate(y):
                z_val = z[y_i,x_i]
                r = np.array([x_val, y_val, z_val])
                d = np.linalg.norm(r-r_t) + np.linalg.norm(r_r-r)
                if d < min:
                    min = d
                    r_min = r

        r_sol = r_min

        lon_computed = np.arctan2(r_sol[1],r_sol[0])*180/np.pi
        lat_computed = np.arcsin(abs(r_sol[2]/np.linalg.norm(r_sol)))*180/np.pi

        lat_sp = metagrp.groups[group].variables['SpecularPointLat'][index].data
        lon_sp = metagrp.groups[group].variables['SpecularPointLon'][index].data

        #print("TDS lat sp: {0} \nTDS lon sp: {1}\n".format(lat_sp, lon_sp) )
        #print("Computed lat sp: {0} \nComputed lon sp computd: {1}\n".format(lat_computed, lon_computed) )
        #print("TDS incidence angle {0:.2f} \nTDS reflection angle: {1:.2f}".format(\
        #    angle_between(ellip_norm(r_sp), (r_r-r_sp))*180/np.pi,\
        #    angle_between(ellip_norm(r_sp), (r_t-r_sp))*180/np.pi))
        #print("Computed incidence angle {0:.2f} \nComputed reflection angle: {1:.2f}\n".format(\
        #    angle_between(ellip_norm(r_sol), (r_r-r_sol))*180/np.pi,\
        #    angle_between(ellip_norm(r_sol), (r_t-r_sol))*180/np.pi))
        
        r_center = r_sol
    return r_center

r_sp = find_sp(r_sp_tds)

################################################# 
# Locate
################################################# 

from scipy.optimize import fsolve

def unit_vector(vector):
    """ Returns the unit vector of the vector. """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2' """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.dot(v1_u, v2_u))

def random_3d_vector():
    """
    Generates a random 3D unit vector (direction) with a uniform spherical distribution
    Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    :return:
    """
    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos( costheta )
    x = np.sin( theta) * np.cos( phi )
    y = np.sin( theta) * np.sin( phi )
    z = np.cos( theta )
    return np.array([x,y,z])


R = np.linalg.norm(r_sp)

# Speed of light
c = 299792458 #m/s
# GPS L1 center frequency
f_c = 1575.42e6 # Hz 

# TODO: plannar equation solver

def equations(p):
    x,y,z = p
    r = np.array([x,y,z])
    
    f = -2300 # Hz
    delay = 7e-6

    tau = c*delay
    tau_r = np.linalg.norm(r-r_t) + np.linalg.norm(r_r-r)
    tau_sp = np.linalg.norm(r_sp-r_t) + np.linalg.norm(r_r-r_sp)
    f1 = -tau + (tau_r - tau_sp)
    
    f_d_r = f_c/c*(np.dot(v_t,unit_vector(r-r_t)) - np.dot(v_r,unit_vector(r_r-r)))
    f_d_sp = f_c/c*(np.dot(v_t,unit_vector(r_sp-r_t)) - np.dot(v_r,unit_vector(r_r-r_sp)))
    f2 = -f + (f_d_r - f_d_sp)
    
    #f3 = R - np.linalg.norm(r)
    f3 = -1 + (x/a)**2 + (y/a)**2 + (z/b)**2 
    
    return (f1, f2, f3)

sol_1 = np.array(fsolve(equations, r_sp, xtol=1e-8, maxfev=1000, epsfcn=1e-24))
sol_2 = sol_1
k = 1
while np.linalg.norm(sol_1 - sol_2) < 0.5:
    r_0 =  r_sp + random_3d_vector()*1e3*k
    sol_2 = np.array(fsolve(equations,r_0))
    k = k+1

sol_3 = sol_2
k = 1
while np.linalg.norm(sol_3 - sol_2) < 0.5 or np.linalg.norm(sol_3 - sol_1) < 0.5:
    r_0 =  r_sp + random_3d_vector()*1e3*k
    sol_3 = np.array(fsolve(equations,r_0))
    k = k+1

def lon_lat(r):
    lon = np.arctan2(r[1],r[0])*180/np.pi
    lat = np.arcsin(abs(r[2]/np.linalg.norm(r)))*180/np.pi
    return lon,lat

lon_sp_tds, lat_sp_tds = lon_lat(r_sp_tds);

print("lat_sp_tds: {}".format(lat_sp_tds))
print("lon_sp_tds: {}".format(lon_sp_tds))

lon_sp, lat_sp = lon_lat(r_sp);

print("lat_sp: {}".format(lat_sp))
print("lon_sp: {}".format(lon_sp))

print('------------- sol 1')

r_sol = np.array(sol_1)

lon = np.arctan2(r_sol[1],r_sol[0])*180/np.pi
lat = np.arcsin(abs(r_sol[2]/np.linalg.norm(r_sol)))*180/np.pi

print("sol lat: {}".format(lat))
print("sol lon: {}".format(lon))

print('------------- sol 2')

r_sol = np.array(sol_2)

lon = np.arctan2(r_sol[1],r_sol[0])*180/np.pi
lat = np.arcsin(abs(r_sol[2]/np.linalg.norm(r_sol)))*180/np.pi

print("sol lat: {}".format(lat))
print("sol lon: {}".format(lon))

print('------------- sol 3')

r_sol = np.array(sol_3)

lon = np.arctan2(r_sol[1],r_sol[0])*180/np.pi
lat = np.arcsin(abs(r_sol[2]/np.linalg.norm(r_sol)))*180/np.pi

print("sol lat: {}".format(lat))
print("sol lon: {}".format(lon))

plt.show()

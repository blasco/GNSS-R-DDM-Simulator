#!/usr/bin/env python

import numpy as np
from shapely import geometry
import matplotlib.pyplot as plt

def find_intersection(contour1,contour2):
  p1 = contour1.collections[0].get_paths()[0]
  v1 = p1.vertices

  p2 = contour2.collections[0].get_paths()[0]
  v2 = p2.vertices

  poly1 = geometry.LineString(v1)
  poly2 = geometry.LineString(v2)

  intersection = poly1.intersection(poly2)

  return intersection

h = 500
elev = 60*np.pi/180

#c = 299792458 #m/s
#time_delay = 

def time_eq(x, y):
    return (x**2 +(y+h/np.tan(elev))**2 + h**2)**(1/2) - h/np.sin(elev) - y*np.cos(elev)

# GPS L1 center frequency
#f_c = 1575.42e6 # Hz 

v_ty = 1
v_tz = 1
v_rx = 15
v_ry = 15
v_rz = 15

def doppler_eq(x, y):
    return -v_ty*np.cos(elev) - v_tz*np.sin(elev) + (v_rx*x + v_ry*(y+h/np.tan(elev)) - v_rz*h)/(x**2 + (y+h/np.tan(elev))**2 + h**2)

extent_x0 =  -10e3
extent_x1 =  10e3
extent_y0 =  -10e3
extent_y1 =  10e3
linsapce_delta = 500
X, Y = np.meshgrid(
        np.linspace(extent_x0, extent_x1, linsapce_delta), 
        np.linspace(extent_y0, extent_y1, linsapce_delta)
        )
Z_time = time_eq(X,Y)
Z_doppler = doppler_eq(X,Y)
#c_doppler = ax_doppler.contour(x, y, z_doppler, [-1.3751])
#c_doppler = ax.contour(x, y, z_doppler, [-1.3751])

# Iso-Delay Contour
time_delay = 2500
fig_time, ax_time = plt.subplots(1,figsize=(10, 4))
ax_time.set_title('Eq Time')
cs_time = ax_time.contourf(X, Y, Z_time, 25, cmap='viridis')
fig_time.colorbar(cs_time)
c_time = ax_time.contour(X, Y, Z_time, [time_delay], colors='r', linestyles='dashdot')
ax_time.grid(c='k', ls='-', alpha=0.3)
plt.show(block=False)

# Iso-Delay Line
fig, ax = plt.subplots(1,figsize=(10, 4))
contour_time = ax.contour(X, Y, Z_time, [time_delay],
        colors='k', 
        linestyles='solid',
        extent=(extent_x0, extent_x1, extent_y0, extent_y1)
        )
plt.show(block=False)

# Iso-Doppler Contour
doppler_delay = -1.3750
fig_doppler, ax_doppler = plt.subplots(1,figsize=(10, 4))
ax_doppler.set_title('Eq doppler')
cs_doppler = ax_doppler.contourf(X, Y, Z_doppler, 25, cmap='viridis')
fig_doppler.colorbar(cs_doppler)
c_doppler = ax_doppler.contour(X, Y, Z_doppler, [doppler_delay], colors='r', linestyles='dashdot')
ax_doppler.grid(c='k', ls='-', alpha=0.3)
plt.show(block=False)

# Iso-Doppler Line
contour_doppler = ax.contour(X, Y, Z_doppler, [doppler_delay],
        colors='b', 
        linestyles='dashdot',
        extent=(extent_x0, extent_x1, extent_y0, extent_y1)
        )
plt.show(block=False)

# Intersections
intersections = find_intersection(contour_time, contour_doppler)

try:
    for i in intersections:
        ax.scatter(i.x, i.y)
except TypeError as te:
    print ('No intersections')

plt.show(block=False)

plt.show()


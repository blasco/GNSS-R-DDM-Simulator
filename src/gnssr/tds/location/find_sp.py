import numpy as np
from netCDF4 import Dataset

def find_sp(r_sp_estimate):
    a = 6378137 #m
    b = 6356752.314245 #m

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

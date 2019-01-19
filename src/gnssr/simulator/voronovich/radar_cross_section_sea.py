#!/usr/bin/env python
def rcs_sea(r):
    '''
    Radar Cross Section of the sea surface.
    Implements Equation 34
        [1]V. U. Zavorotny and A. G. Voronovich, “Scattering of GPS signals from 
        the ocean with wind remote sensing application,” IEEE Transactions on 
        Geoscience and Remote Sensing, vol. 38, no. 2, pp. 951–964, Mar. 2000.  
    '''
    q = np.linalg.norm(scattering_vector(r))
    q_z = np.linalg.norm(scattering_vector(r)[2])
    return np.pi*(fresnel_coefficient**2)*((q/q_z)**4)



def slope_probability_density_function(x):
    '''
    Implements Equation 4
        [1]Q. Yan and W. Huang, “GNSS-R Delay-Doppler Map Simulation Based on the 
        2004 Sumatra-Andaman Tsunami Event,” Journal of Sensors, vol. 2016, pp. 
        1–14, 2016.  
    '''
    # phi_0 is the angle between the up-down wind direction and the x-axis
    phi_0 = 0
    wind_rotation = np.array([
        [np.cos(phi_0), -np.sin(phi_0)],
        [np.sin(phi_0),  np.cos(phi_0)]
        ])
    covariance = np.array([
        [variance_upwind(wind_speed_10m_above_sea), 0],
        [0, variance_crosswind(wind_speed_10m_above_sea)]
        ])
    W = (wind_rotation.dot(covariance)).dot(np.transpose(wind_rotation)

    return 1/(2*np.pi*(numpy.linalg.det(W)**(1/2))) \
            *np.exp( \
                -1/2*(np.transpose(x).dot(np.linalg.inv(W))).dot(x) \
            )

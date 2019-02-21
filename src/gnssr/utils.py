import numpy as np
from datetime import *
from dateutil import tz

earth_semimajor_axis = 6378137 # meters
earth_semiminor_axis = 6356752.314245 # meters
earth_inertial_angular_speed = 7.2921158553e-5 # rad/sec
        
gps_ca_chips_per_second = 1.023e6

light_speed = 299792458.0 # m/s

def normalize(mat):
    return (mat - np.min(mat))/(np.max(mat)-np.min(mat))

def unit_vector(r):
    return r/np.linalg.norm(r)

def angle_between(v1, v2):
    """ 
    Returns the angle in radians between vectors 'v1' and 'v2' 
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.dot(v1_u, v2_u))

def ellip_norm(r):
    return np.array([2*r[0]/earth_semimajor_axis**2,2*r[1]/earth_semimajor_axis**2,2*r[2]/earth_semiminor_axis**2])

def datenum_to_pytime(matlab_datenum):
    """
    Formated using ISO 8601
    """
    python_datetime = datetime.fromordinal(int(matlab_datenum)) + timedelta(days=matlab_datenum%1) - timedelta(days = 366)
    return str(python_datetime.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]) + 'Z'
    #return python_datetime

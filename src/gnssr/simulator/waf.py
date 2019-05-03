import numpy as np

def woodward_ambiguity_function(delay, doppler, sim_config): 
    """ 
    The Woodward Ambiguity Function (WAF) can be understood as the impulse 
    response to the scattered signal from a single delay-Doppler cell. 
    
    It can be approximated by the product of a function dependent on the delay 
    increment and a function dependent on the doppler increment with respect to 
    the specular point.

    Implements equation 22:
        V. U. Zavorotny and A. G. Voronovich, “Scattering of GPS signals from the 
        ocean with wind remote sensing application,” IEEE Transactions on Geoscience and 
        Remote Sensing, vol. 38, no. 2, pp. 951–964, Mar. 2000.  

    Args:
        delay (numpy.ndarray with size (1,)): Delay increment.
        doppler (numpy.ndarray with size (1,)): Doppler increment.
        sim_config: Instance of simulation_configuration class.

    Returns:
        numpy.ndarray with size(1,).
    """ 
    return waf_delay(delay, sim_config) * np.absolute(waf_frequency(doppler, sim_config))

def waf_delay(delay, sim_config):
    """
    Delay increment component of the Woodward Ambiguity Function (WAF) decomposition.

    Implements equation 20:
        V. U. Zavorotny and A. G. Voronovich, “Scattering of GPS signals from the 
        ocean with wind remote sensing application,” IEEE Transactions on Geoscience and 
        Remote Sensing, vol. 38, no. 2, pp. 951–964, Mar. 2000.  

    -delay/coherent_integration_time is approximated by zero as explained by 
    the end of the paragraph below the referenced equation.

    Args:
        delay (numpy.ndarray with size(1,)): Delay increment.
        sim_config: Instance of simulation_configuration class.

    Returns:
        numpy.ndarray with size(1,).
    """
    delay_chip = sim_config.delay_chip
    return np.where(np.abs(delay) <= delay_chip, 
            1 - np.abs(delay/delay_chip), 
            0)

def waf_frequency(doppler, sim_config):
    """
    Doppler increment component of the Woodward Ambiguity Function (WAF) decomposition.

    Implements equation 21:
        V. U. Zavorotny and A. G. Voronovich, “Scattering of GPS signals from the 
        ocean with wind remote sensing application,” IEEE Transactions on Geoscience and 
        Remote Sensing, vol. 38, no. 2, pp. 951–964, Mar. 2000.  

    A coherent_integration_time factor is missing in the function defined in 
    the paper.

    Args:
        doppler (numpy.ndarray with size(1,)): Doppler increment.
        sim_config: Instance of simulation_configuration class.

    Returns:
        numpy.ndarray with size(1,).
    """
    coherent_integration_time = sim_config.coherent_integration_time
    x = np.pi*doppler*coherent_integration_time
    return np.absolute(
                np.where(np.abs(x) <= 0.2, 
                    1-(x**2)/6+(x**4)/120, # Taylor expansion around 0
                    np.sin(x)/(x)
                ) * \
                np.exp(-1j*x)
            )

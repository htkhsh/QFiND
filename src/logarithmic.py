import numpy as np
from scipy.integrate import quad
from init import setpara, const, opt
from specdens import sbeta  # Assuming sbeta is defined in specdens module

icm2ifs = const['icm2ifs']

def log_disc(Msp, Omegac, lamb=1.1):
    """
    Divide the frequency domain into logarithmically spaced intervals,
    compute coefficients 'g' and frequencies 'w' using numerical integration.

    Parameters:
    - Msp (int): Number of spectral points (must be even).
    - Omegac (float): Maximum frequency (cut-off frequency).

    Returns:
    - w (ndarray): Array of frequencies.
    - g (ndarray): Array of coefficients.
    """
    
    if Msp % 2 != 0:
        raise ValueError("The number of spectral points must be even.")
    
    # Initialize parameters
    Msp2 = Msp // 2  # Half the number of spectral points

    # Initialize arrays
    w = np.zeros(Msp)
    g = np.zeros(Msp)
    omegaj = np.zeros((2, Msp2))

    # Compute frequency intervals for negative and positive frequencies
    for j in range(Msp2):
        # Adjust indices for Python (0-based indexing)
        omegaj[0, j] = Omegac / lamb**(j + 1)  
        omegaj[1, j] = Omegac / lamb**(j)      

    # Define the integrand functions
    def integrand1(omega):
        return sbeta(omega)

    def integrand2(omega):
        return sbeta(omega) * omega

    # Compute 'g' and 'w' for negative frequencies
    for j in range(Msp2):
        # Integrate 'integrand1' from -omegaj[1, j] to -omegaj[0, j]
        res, _ = quad(integrand1, -omegaj[1, j], -omegaj[0, j], epsabs=1e-12)
        g[j] = res

        # Integrate 'integrand2' over the same interval
        res, _ = quad(integrand2, -omegaj[1, j], -omegaj[0, j], epsabs=1e-12)
        w[j] = res / g[j] if g[j] != 0 else 0  # Avoid division by zero

    # Compute 'g' and 'w' for positive frequencies
    for j in range(Msp2):
        # Integrate 'integrand1' from omegaj[0, j] to omegaj[1, j]
        res, _ = quad(integrand1, omegaj[0, j], omegaj[1, j], epsabs=1e-12)
        g[j + Msp2] = res

        # Integrate 'integrand2' over the same interval
        res, _ = quad(integrand2, omegaj[0, j], omegaj[1, j], epsabs=1e-12)
        idx = j + Msp2
        g_val = g[idx]
        w[idx] = res / g_val if g_val != 0 else 0  # Avoid division by zero

    return w, g

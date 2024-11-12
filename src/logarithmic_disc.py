import numpy as np
from scipy.integrate import quad
from init import setpara, const, opt
from specdens import sbeta  # Assuming sbeta is defined in specdens module

icm2ifs = const['icm2ifs']
M = opt['M']
N = opt['N']
N = int(N/2)
Omegac = opt['Omegac']
Omega_max = opt['Omega_max']
s = opt['s']
gamc = opt['gamc']
alpha = opt['alpha']
omega0 = gamc / N * (1.0 - np.exp(-Omegac / gamc))

def log_disc(Nsp):
    """
    Divide the frequency domain into logarithmically spaced intervals,
    compute coefficients 'g' and frequencies 'w' using numerical integration.

    Parameters:
    - Nsp (int): Number of spectral points (must be even).
    - Omegac (float): Maximum frequency (cut-off frequency).
    - stype (str): Spectral density type, used in sbeta function.
    - temperature (float): Temperature value, used in sbeta function.
    - bho (float): Beta * ħ (Planck's constant), used in sbeta function.

    Returns:
    - w (ndarray): Array of frequencies.
    - g (ndarray): Array of coefficients.
    """
    # Initialize parameters
    lamb = 1.1  # Logarithmic spacing factor
    Nsp2 = Nsp // 2  # Half the number of spectral points

    # Initialize arrays
    w = np.zeros(Nsp)
    g = np.zeros(Nsp)
    omegaj = np.zeros((2, Nsp2))

    # Compute frequency intervals for negative and positive frequencies
    for j in range(Nsp2):
        # Adjust indices for Python (0-based indexing)
        omegaj[0, j] = Omegac / lamb**(j + 1)  # Corresponds to omegaj(1, j) in Fortran
        omegaj[1, j] = Omegac / lamb**(j)      # Corresponds to omegaj(2, j) in Fortran

    # Define the integrand functions
    def integrand1(omega):
        return sbeta(omega)

    def integrand2(omega):
        return sbeta(omega) * omega

    # Compute 'g' and 'w' for negative frequencies
    for j in range(Nsp2):
        # Integrate 'integrand1' from -omegaj[1, j] to -omegaj[0, j]
        res, _ = quad(integrand1, -omegaj[1, j], -omegaj[0, j], epsabs=1e-12)
        g[j] = res

        # Integrate 'integrand2' over the same interval
        res, _ = quad(integrand2, -omegaj[1, j], -omegaj[0, j], epsabs=1e-12)
        w[j] = res / g[j] if g[j] != 0 else 0  # Avoid division by zero

    # Compute 'g' and 'w' for positive frequencies
    for j in range(Nsp2):
        # Integrate 'integrand1' from omegaj[0, j] to omegaj[1, j]
        res, _ = quad(integrand1, omegaj[0, j], omegaj[1, j], epsabs=1e-12)
        g[j + Nsp2] = res

        # Integrate 'integrand2' over the same interval
        res, _ = quad(integrand2, omegaj[0, j], omegaj[1, j], epsabs=1e-12)
        idx = j + Nsp2
        g_val = g[idx]
        w[idx] = res / g_val if g_val != 0 else 0  # Avoid division by zero

    return w, g

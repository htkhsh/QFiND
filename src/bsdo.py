import numpy as np
import sys
from init import setpara, const, opt
from scipy.integrate import quad
from specdens import sbeta
from lanczos import lanczos

f = open(str(sys.argv[1]), mode='r')
setpara(f)
icm2ifs = const['icm2ifs']
N_t = opt['N_t']
Tc = opt['Tc']
#Temp = opt['temperature']

def bsdo(N_w, Omega_min, Omega_max, Msp):
    """
    Compute the discretization of the spectral density using the BSDO method.

    Parameters:
    - N_w (int) : The number of discretization points.
    - Omega_min (float): The minimum frequency.
    - Omega_max (float): The maximum frequency.
        
    Returns:
    - k (ndarray) : The discretized frequencies.
    - zk (ndarray) : The squared first components of the eigenvectors (weights).
    """
    w = np.linspace(Omega_min, Omega_max, N_w) 
    j = sbeta(w)
    j[j < 0] = 0.0 

    wj = np.column_stack((w, j))

    # Compute discretized frequencies (wd) and weights (zd)
    wk, zk = orthpoly_discretization(Msp, wj)
    if Omega_min < 0:
        norm = quad(lambda w: sbeta(w), Omega_min, 0)[0] + quad(lambda w: sbeta(w), 0, Omega_max)[0]
    else:
        norm = quad(lambda w: sbeta(w), Omega_min, Omega_max)[0]
    zk = zk * norm

    return wk, zk


def orthpoly_discretization(N, wj):
    """
    Compute the discretization of orthogonal polynomials.

    Parameters:
    - N (int) : The number of discretization points.
    - wj (ndarray):
        An n x 2 array where:
        - wj[:,0] is the set of support points of the (quantum noise) spectral density (the frequencies).
        - wj[:,1] is the (quantum noise) spectral density.

    Returns:
    -wd (ndarray) : The eigenvalues (frequencies).
    -zd (ndarray) : The squared first components of the eigenvectors (weights).
    """
    # Ensure wj is a NumPy array
    wj = np.asarray(wj)
    
    # Compute the recurrence coefficients using the Lanczos algorithm
    ab = lanczos(N, wj)
    
    # Construct the tridiagonal matrix M
    diag_main = np.diag(ab[:, 0])
    off_diag_elements = np.sqrt(ab[1:N, 1])
    diag_upper = np.diag(off_diag_elements, k=1)
    diag_lower = np.diag(off_diag_elements, k=-1)
    M = diag_main + diag_upper + diag_lower
    
    # Compute the eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eigh(M)
    
    # Compute zd and wd
    zd = eigenvectors[0, :] ** 2
    wd = eigenvalues

    return wd, zd
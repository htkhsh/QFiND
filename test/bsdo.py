import numpy as np
import sys
from init import setpara, const, opt
from scipy.integrate import quad
from specdens import s_be, sbeta, sdens
from lanczos import lanczos

f = open(str(sys.argv[1]), mode='r')
setpara(f)
icm2ifs = const['icm2ifs']
N_w = opt['N_w']
Temp = opt['temperature']

def bsdo(N_w, Omega_min, Omega_max, Msp=N_w):
    w = np.linspace(Omega_min, Omega_max, N_w) 
    j = sbeta(w)
    j[j < 0] = 0.0  # Set negative values to zero

    wj = np.column_stack((w, j))

    # Compute discretized frequencies (wd) and weights (zd)
    wk, zk = orthpoly_discretization(Msp, wj)
    if Temp < 1e-12:
        norm = quad(lambda w: sbeta(w), Omega_min, Omega_max)[0]
    else:
        norm = quad(lambda w: sbeta(w), Omega_min, 0)[0] + quad(lambda w: sbeta(w), 0, Omega_max)[0]
    zk = zk * norm 

    return wk, zk

"""
def bsdo_be():
    if Temp < 1e-10:
        w = np.linspace(1e-15, Omegac, N_w) 
    else:
        w = np.linspace(-Omegac, Omegac, N_w) 

    j = s_be(w)
    j[j < 0] = 0.0  # Set negative values to zero

    wj = np.column_stack((w, j))

    # Compute discretized frequencies (wd) and weights (zd)
    wk, zk = orthpoly_discretization(Nbsdo, wj)
    if Temp < 1e-10:
        norm = quad(lambda w: s_be(w), 0, Omegac)[0]
    else:
        norm = quad(lambda w: s_be(w), -Omegac, 0)[0] + quad(lambda w: s_be(w), 0, Omegac)[0]
    zk = zk * norm

    return wk, zk


def bsdo_sd():
    if Temp < 1e-10:
        w = np.linspace(1e-15, Omegac, N_w) 
    else:
        w = np.linspace(-Omegac, Omegac, N_w) 

    j = sdens(w)
    j[j < 0] = 0.0  # Set negative values to zero

    wj = np.column_stack((w, j))

    # Compute discretized frequencies (wd) and weights (zd)
    wk, zk = orthpoly_discretization(Nbsdo, wj)
    if Temp < 1e-10:
        norm = quad(lambda w: sdens(w), 0, Omegac)[0]
    else:
        norm = quad(lambda w: sdens(w), -Omegac, 0)[0] + quad(lambda w: sdens(w), 0, Omegac)[0]
    zk = zk * norm

    return wk, zk
"""

def orthpoly_discretization(N, wj):
    """
    Compute the discretization of orthogonal polynomials.

    Parameters:
    N  : int
        The number of discretization points.
    wj : ndarray of shape (n, 2)
        An n x 2 array where:
        - wj[:,0] is the set of support points of the spectral density (the frequencies).
        - wj[:,1] is the spectral density.

    Returns:
    wd : ndarray
        The eigenvalues (frequencies).
    zd : ndarray
        The squared first components of the eigenvectors (weights).
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
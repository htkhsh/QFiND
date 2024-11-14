import numpy as np
import sys
from init import setpara, const, opt
from scipy.integrate import quad
from specdens import sbeta
from lanczos import lanczos

f = open(str(sys.argv[1]), mode='r')
setpara(f)
icm2ifs = const['icm2ifs']
M = opt['M']
N = opt['N']
Tc = opt['Tc']
Temp = opt['temperature']
Omegac = opt['Omegac']

def bsdo():
    if Temp < 1e-10:
        w = np.linspace(1e-15, Omegac, M) 
    else:
        w = np.linspace(-Omegac, Omegac, M) 
    j = sbeta(w)
    j[j < 0] = 0.0  # Set negative values to zero

    wj = np.column_stack((w, j))

    # Compute discretized frequencies (wd) and weights (zd)
    wk, zk = orthpoly_discretization(N, wj)
    if Temp < 1e-10:
        norm = quad(lambda w: sbeta(w), 0, Omegac)[0]
    else:
        norm = quad(lambda w: sbeta(w), -Omegac, 0)[0] + quad(lambda w: sbeta(w), 0, Omegac)[0]
    zk = zk*norm

    return wk, zk


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
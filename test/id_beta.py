import numpy as np
from init import const
from id_sub import id_freq_eps, id_freq_rank
from specdens import sbeta, s_be
from corrfunc import C_beta
from scipy.optimize import nnls
from bsdo import bsdo

icm2ifs = const['icm2ifs']
def edr_id_beta(N_tau, N_w, beta, omega_min, omega_max, eps, frank, rand=False):
    """
    Perform frequency estimation using interpolative decomposition (ID) and NNLS.

    Parameters:
    - N_tau (int): Number of imaginary time points.
    - N_w (int): Number of frequency points.
    - beta (float): Maximum imaginary time value (inverse temperature).
    - omega_min (float): Minimum frequency value.
    - omega_max (float): Maximum frequency value.
    - eps (float): Error tolerance for ID.
    - krank (int): Rank for the ID. If smaller than 1, ID uses error tolerance.

    Returns:
    - Nsp (int): Number of estimated frequencies.
    - w (ndarray): Estimated frequencies.
    - g (ndarray): Estimated coefficients (amplitudes).
    - krank (int): Updated rank after ID.
    """
    # Generate time mesh tf and frequency mesh wf
    tau, w = equispaced_mesh(N_tau,N_w,beta,omega_min,omega_max)

    # Create core matrix Kf
    f = create_integrand(tau,w)

    # Perform Interpolative Decomposition (ID)
    if frank < 1:
        frank, idx, B, err1 = id_freq_eps(f, eps, rand)
    else:
        idx, B, err1 = id_freq_rank(f, frank, rand)
        # krank remains the same
    
    print("Rank of f: ", frank)
    # Compute estimated frequencies wk
    wk = w[idx[:frank]]

    # Compute coefficients g using NNLS or another method
    zk, err2 = edr_coef(tau, B)

    ind = np.argsort(wk)
    wk = wk[ind]
    zk = zk[ind]

    Nsp = frank
    # Remove zero coefficients if any
    if np.min(zk) == 0.0:
        ind = np.where(zk > 0)[0]
        wk = wk[ind]
        zk = zk[ind]
        Nsp = len(wk)

    print("Number of sample points: ", Nsp)
    print("Error in ID: ", err1)
    print("Error in NNLS: ", err2)

    return Nsp, wk, zk, frank


def edr_coef(tau, B):
    """
    Use non-negative least squares (NNLS) to estimate coefficients g.

    Parameters:
    - tau (ndarray): Imaginary time array (size: N_tau)
    - frank (int): Approximation rank
    - B (ndarray): Input matrix (size: (2*N_t, frank))

    Returns:
    - g (ndarray): Estimated coefficients (size: frank)
    - err (float): Estimation error
    """
    N_tau = len(tau)
    c = np.zeros(N_tau)

    # Construct vector c with C_beta(tau_j)
    for j in range(N_tau):
        tau_j = tau[j]
        c[j] = C_beta(tau_j)

    # Solve the NNLS problem: minimize ||B * g - c|| subject to g >= 0
    g, err = nnls(B, c)

    return g, err


def equispaced_mesh(N_tau, N_w, beta, omega_min, omega_max):

    # Imaginary time grid (tau)
    tau = np.linspace(0, beta, N_tau)

    # Frequency grid (w)
    w = np.linspace(omega_min, omega_max, N_w)
    w = w * icm2ifs

    return tau, w


def create_integrand(tau, w):
    
    N_tau = len(tau)
    N_w = len(w)
    f = np.zeros((N_tau,N_w),dtype=float)
    for j in range(N_tau):
        for k in range(N_w):
            #f[j, k] = s_be(w[k],icm2ifs) * np.exp(w[k] * tau[j])
            f[j, k] = sbeta(w[k],icm2ifs) * np.exp(-w[k] * tau[j])
    
    return f
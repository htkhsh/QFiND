import numpy as np
from init import const
from id_sub import id_freq_eps, id_freq_rank
from specdens import sbeta
from corrfunc import S_exact, A_exact
from scipy.optimize import nnls

icm2ifs = const['icm2ifs']
def edr_id(N_t, N_w, tc, omega_min, omega_max, eps, frank, rand=False):
    """
    Perform frequency estimation using interpolative decomposition (ID) and NNLS.

    Parameters:
    - N_t (int): Number of time points.
    - N_w (int): Number of frequency points.
    - tc (float): Maximum time value.
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
    t, w = equispaced_mesh(N_t,N_w,tc,omega_min,omega_max)

    # Create core matrix Kf
    f = create_integrand(t,w)

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
    zk, err2 = edr_coef(t, B)

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


def edr_coef(t, B):
    """
    Use non-negative least squares (NNLS) to estimate coefficients g.

    Parameters:
    - t (ndarray): Time array (size: N_t)
    - frank (int): Approximation rank
    - B (ndarray): Input matrix (size: (2*N_t, frank))

    Returns:
    - g (ndarray): Estimated coefficients (size: frank)
    - err (float): Estimation error
    """
    N_t = len(t)
    c = np.zeros(2*N_t)

    # Construct vector c with S_exact(t) and A_exact(t)
    for i in range(N_t):
        ti = t[i]
        c[i] = S_exact(ti)
    for i in range(N_t):
        ti = t[i]
        c[N_t + i] = A_exact(ti)

    # Solve the NNLS problem: minimize ||B * g - c|| subject to g >= 0
    g, err = nnls(B, c)

    return g, err


def equispaced_mesh(N_t, N_w, tc, omega_min, omega_max):
    """
    Generate equispaced time and frequency grids.

    Parameters:
    - N_t (int): Number of time points.
    - N_w (int): Number of frequency points.
    - tc (float): Maximum time value.
    - omega_min (float): Minimum frequency value.
    - omega_max (float): Maximum frequency value.

    Returns:
    - t (ndarray): Time grid.
    - w (ndarray): Frequency grid.
    """
    # Time grid (t)
    t = np.linspace(0,tc,N_t)

    # Frequency grid (w)
    w = np.linspace(omega_min,omega_max,N_w)
    w = w * icm2ifs

    return t, w


def create_integrand(t, w):
    """
    Create the matrix f for the interpolative decomposition.

    Parameters:
    - t (ndarray): Time grid.
    - w (ndarray): Frequency grid.

    Returns:
    - f (ndarray): Matrix f with shape (2*N_t, N_w).
    """
    N_t = len(t)
    N_w = len(w)
    f = np.zeros((2*N_t,N_w),dtype=float)

    # Fill the first N_t rows of f with the real part of an integrand
    for i in range(N_t):
        for j in range(N_w):
            f[i, j] = sbeta(w[j],icm2ifs) * np.cos(w[j] * t[i])
    
    # Fill the next N_t rows of f with the imaginary part of an integrand
    for i in range(N_t, 2*N_t):
        for j in range(N_w):
            f[i, j] = -sbeta(w[j],icm2ifs) * np.sin(w[j] * t[i-N_t])
    
    return f






    
    


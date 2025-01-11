import numpy as np
from init import const
from id_sub import id_freq_eps, id_freq_rank
from specdens import sbeta
from corrfunc import S_mix, A_mix
from scipy.optimize import nnls

icm2ifs = const['icm2ifs']
def edr_id_mix(N_t, N_tau, N_w, tc, beta, omega_min, omega_max, eps, frank, id_in_time=False, rand=False):
    """
    Perform frequency estimation using interpolative decomposition (ID) and NNLS.

    Parameters:
    - N_t (int): Number of time points.
    - N_tau (int): Number of imaginary time points.
    - N_w (int): Number of frequency points.
    - tc (float): Maximum time value.
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
    t, tau, w = equispaced_mesh(N_t,N_tau,N_w,tc,beta,omega_min,omega_max)

    # Create core matrix Kf
    f = create_integrand(t,tau,w)
    
    if id_in_time:
        print("performing ID in time")
        frank_s = 5000
        fT = f.transpose()
        idx_s, B1, err_s = id_freq_rank(fT, frank_s, True)
        f = B1.transpose()

    # Perform Interpolative Decomposition (ID)
    print("performing ID in frequency")
    if frank < 1:
        frank, idx, B, err1 = id_freq_eps(f, eps, rand)
    else:
        idx, B, err1 = id_freq_rank(f, frank, rand)
        # krank remains the same
    
    print("Rank of f: ", frank)
    # Compute estimated frequencies wk
    wk = w[idx[:frank]]

    # Compute coefficients g using NNLS or another method
    if id_in_time:
        s = create_composite_time(t, tau)
        zk, err2 = edr_coef_s(s, idx_s, B)
    else:
        zk, err2 = edr_coef(t, tau, B)

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
    if id_in_time:
        print("Error in ID (time): ", err_s)
    print("Error in ID: ", err1)
    print("Error in NNLS: ", err2)

    return Nsp, wk, zk, frank


def edr_coef(t, tau, B):
    """
    Use non-negative least squares (NNLS) to estimate coefficients g.

    Parameters:
    - t (ndarray): Time array (size: N_t)
    - tau (ndarray): Imaginary time array (size: N_tau)
    - frank (int): Approximation rank
    - B (ndarray): Input matrix (size: (2*N_t, frank))

    Returns:
    - g (ndarray): Estimated coefficients (size: frank)
    - err (float): Estimation error
    """
    N_t = len(t)
    N_tau = len(tau)
    c = np.zeros(2*N_t*N_tau)

    # Construct vector c with S_mix(t) and A_mix(t)
    for i in range(N_t):
        for j in range(N_tau):
            t_i = t[i]
            tau_j = tau[j]
            c[i*N_tau+j] = S_mix(t_i, tau_j)
    for i in range(N_t):
        for j in range(N_tau):
            t_i = t[i]
            tau_j = tau[j]
            c[(i+N_t)*N_tau+j] = A_mix(t_i, tau_j)

    # Solve the NNLS problem: minimize ||B * g - c|| subject to g >= 0
    g, err = nnls(B, c)

    return g, err


def edr_coef_s(s, frank_s, idx_s, B):
    """
    Use non-negative least squares (NNLS) to estimate coefficients g.

    Parameters:
    - s (dict): Composite time array (size: N_s)
    - frank_s (int): Approximation rank
    - idx_s (ndarray): Indices of selected composite times (size: N_s)
    - B (ndarray): Input matrix (size: (N_s, frank))

    Returns:
    - g (ndarray): Estimated coefficients (size: frank)
    - err (float): Estimation error
    """
    N_s = len(s)
    c = np.zeros(frank_s)

    # Construct vector c with S_mix(t) and A_mix(t)
    for i in range(frank_s):
        t = s[idx_s[i]]['t']
        tau = s[idx_s[i]]['tau']
        if i < N_s:
            c[i] = S_mix(t, tau)
        else:
            c[i] = A_mix(t, tau)

    # Solve the NNLS problem: minimize ||B * g - c|| subject to g >= 0
    g, err = nnls(B, c)

    return g, err


def equispaced_mesh(N_t, N_tau, N_w, tc, beta, omega_min, omega_max):

    # Real time grid (t)
    t = np.linspace(0, tc, N_t)

    # Imaginary time grid (tau)
    tau = np.linspace(0, beta, N_tau)

    # Frequency grid (w)
    w = np.linspace(omega_min, omega_max, N_w)
    w = w * icm2ifs

    return t, tau, w

def create_composite_time(t, tau):
    """
    Create a composite time array.

    Parameters:
    - t (ndarray): Real time array.
    - tau (ndarray): Imaginary time array.

    Returns:
    - s (dict): Composite time.
    """

    N_t = len(t)
    N_tau = len(tau)
    # Dictionary of composite time (s)
    s = {}
    for i in range(N_t):
        for j in range(N_tau):
            s[i*N_tau+j] = {'t': t[i], 'tau': tau[j]}

    return s
    

def create_integrand(t, tau, w):
    
    N_t = len(t)
    N_tau = len(tau)
    N_w = len(w)
    f = np.zeros((2*N_t*N_tau,N_w),dtype=float)

    # Fill the first M rows of K with the real part of an integrand
    for i in range(N_t):
        for j in range(N_tau):
            for k in range(N_w):
                f[i*N_tau+j, k] = sbeta(w[k],icm2ifs) * np.cos(w[k] * t[i]) * np.exp(-w[k] * tau[j])
    
    # Fill the next M rows of K with the imaginary part of an integrand
    for i in range(N_t, 2*N_t):
        for j in range(N_tau):
            for k in range(N_w):
                f[i*N_tau+j, k] = sbeta(w[k],icm2ifs) * np.sin(w[k] * t[i-N_t]) * np.exp(-w[k] * tau[j])
    
    return f

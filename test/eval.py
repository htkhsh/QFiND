import numpy as np
from init import opt, const
from specdens import sbeta, s_be
from corrfunc import S_exact, A_exact, C_beta, S_mix, A_mix
from plot import plot_bcf, plot_bcf_float

hbar = const['hbar']
kb = const['kb']
icm2ifs = const['icm2ifs']
if opt['BCF_type'] == "Real" or opt['BCF_type'] == "Mix":
    N_t = opt['N_t']
    Tc = opt['Tc']
if opt['BCF_type'] == "Imag" or opt['BCF_type'] == "Mix":
    N_tau = opt['N_tau']
Temp = opt['temperature']
if Temp > 0.0:
    beta = hbar*1e15/kb/Temp
    
def calc_error(wk, zk):
    """
    Calculate the error between the approximated and exact correlation functions.

    Parameters:
    - Nsp (int): Number of spectral points.
    - wk (ndarray): Array of frequencies.
    - gk (ndarray): Array of coefficients.

    Returns:
    None
    """
    if opt['BCF_type'] == "Real":
        calc_error_t(wk, zk)
    elif opt['BCF_type'] == "Imag":
        calc_error_beta(wk, zk)
    elif opt['BCF_type'] == "Mix":
        calc_error_t(wk, zk)
        calc_error_beta(wk, zk)


def calc_error_t(wk, zk):
    t = np.linspace(0, Tc, N_t)
    approx = np.zeros(N_t, dtype=complex)
    exact = np.zeros(N_t, dtype=complex)
    error = np.zeros(N_t, dtype=complex)
    for i in range(N_t):
        ti = t[i]
        # Compute approximations\
        approx[i] = C_t(ti, zk, wk)
        # Compute exact values
        exact[i] = S_exact(ti) + 1j * A_exact(ti)
        # Compute error
        error[i] = (approx[i] - exact[i])
    # Print normalized errors
    norm = S_exact(0.0)
    normalized_max_error = np.max(np.abs(error)) / norm
    normalized_avg_error = np.sum(np.abs(error)) / (norm * N_t)
    print("Normalized maximum error:", normalized_max_error)
    print("Normalized average error:", normalized_avg_error)
    # Plot a BCF
    plot_bcf(t, exact, approx, error/norm)


def calc_error_beta(wk, zk):
    norm = C_beta(0)
    tau = np.linspace(0, beta, N_tau)
    approx = np.zeros(N_tau, dtype=float)
    exact = np.zeros(N_tau, dtype=float)
    error = np.zeros(N_tau, dtype=float)
    for i in range(N_tau):
        taui = tau[i]
        # Compute approximations\
        approx[i] = C_tau(taui, zk, wk)
        # Compute exact values
        exact[i] = C_beta(taui)
        # Compute error
        error[i] = (approx[i] - exact[i])
    # Print normalized errors
    norm = np.max(np.abs(exact))
    normalized_max_error = np.max(np.abs(error)) / norm
    normalized_avg_error = np.sum(np.abs(error)) / (norm * N_tau)
    print("Normalized maximum error:", normalized_max_error)
    print("Normalized average error:", normalized_avg_error)
    # Plot a BCF
    plot_bcf_float(tau, exact, approx, error/norm)


def C_t(t, zk, wk):
    """
    Compute the approximation of the correlation function.

    Parameters:
    - t (float): Time value.
    - Min (int): Maximum index (inclusive).
    - zk (ndarray): Array of coefficients (complex).
    - wk (ndarray): Array of frequencies (complex).

    Returns:
    - res (complex): Result of the approximation at time t.
    """

    N = len(wk)
    res = 0.0 + 0.0j
    for m in range(N): 
        res += zk[m] * sbeta(wk[m],icm2ifs) * np.exp(-1j * wk[m] * t)
    return res


def C_tau(tau, zk, wk):
    """
    Compute the approximation of the imaginary-time correlation function.

    Parameters:
    - tau (float): Imaginary time value.
    - rho (ndarray): Array of coefficients (complex).
    - alpha (ndarray): Array of frequencies (complex).

    Returns:
    - res (float): Result of the approximation at time t.
    """

    N = len(wk)
    res = 0.0
    for m in range(N): 
        res += zk[m] * s_be(wk[m],icm2ifs) * np.exp(wk[m] * tau)
    return res

def C_mix(t, tau, zk, wk):
    """
    Compute the approximation of the imaginary-time correlation function.

    Parameters:
    - t (float): Real time value.
    - tau (float): Imaginary time value.
    - rho (ndarray): Array of coefficients (real).
    - alpha (ndarray): Array of frequencies (complex).

    Returns:
    - res (float): Result of the approximation at time t.
    """

    N = len(wk)
    res = 0.0
    for m in range(N): 
        res += zk[m] * s_be(wk[m],icm2ifs) * np.exp(-wk[m] * (1j * t - tau))
    return res
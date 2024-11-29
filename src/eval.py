import numpy as np
from init import opt, const
from specdens import sbeta
from corrfunc import S_exact, A_exact
from plot import plot_bcf

N_t = opt['N_t']
Tc = opt['Tc']
icm2ifs = const['icm2ifs']
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
    
    t = np.linspace(0, Tc, N_t)
    approx = np.zeros(N_t, dtype=complex)
    exact = np.zeros(N_t, dtype=complex)
    error = np.zeros(N_t, dtype=complex)
    for i in range(N_t):
        ti = t[i]
        # Compute approximations
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


def C_t(t, zk, wk):
    """
    Compute the approximation of the correlation function.

    Parameters:
    - t (float): Time value.
    - gk (ndarray): Array of coefficients (complex).
    - wk (ndarray): Array of frequencies (complex).

    Returns:
    - res (complex): Result of the approximation at time t.
    """

    N = len(wk)
    res = 0.0 + 0.0j
    for m in range(N): 
        res += zk[m] * sbeta(wk[m], icm2ifs) * np.exp(-1j * wk[m] * t)
    return res
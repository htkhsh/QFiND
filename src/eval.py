import numpy as np
from init import opt, const
from specdens import sbeta
from corrfunc import S_exact, A_exact
from plot import plot_bcf

Nt = opt['Nt']
Tc = opt['Tc']
icm2ifs = const['icm2ifs']
def calc_error(wk, gk):
    """
    Calculate the error between the approximated and exact correlation functions.

    Parameters:
    - Nsp (int): Number of spectral points.
    - wk (ndarray): Array of frequencies.
    - gk (ndarray): Array of coefficients.

    Returns:
    None
    """
    
    c0 = S_exact(0.0)
    t = np.linspace(0, Tc, Nt)
    approx = np.zeros(Nt, dtype=complex)
    exact = np.zeros(Nt, dtype=complex)
    error = np.zeros(Nt, dtype=complex)
    for i in range(Nt):
        ti = t[i]
        # Compute approximations
        approx[i] = approximation(ti, gk, wk)
        # Compute exact values
        exact[i] = S_exact(ti) + 1j * A_exact(ti)
        # Compute error
        error[i] = (approx[i] - exact[i])

    # Print normalized errors
    normalized_max_error = np.max(np.abs(error)) / c0
    normalized_avg_error = np.sum(np.abs(error)) / (c0 * Nt)
    print("Normalized maximum error:", normalized_max_error)
    print("Normalized average error:", normalized_avg_error)

    # Plot a BCF
    plot_bcf(t, exact/c0, approx/c0, error/c0)


def approximation(t, rho, alpha):
    """
    Compute the approximation of the correlation function.

    Parameters:
    - t (float): Time value.
    - Min (int): Maximum index (inclusive).
    - rho (ndarray): Array of coefficients (complex).
    - alpha (ndarray): Array of frequencies (complex).

    Returns:
    - res (complex): Result of the approximation at time t.
    """

    N = len(alpha)
    res = 0.0 + 0.0j
    for m in range(N): 
        res += rho[m] * np.exp(-1j * alpha[m] * t)
    return res
import numpy as np
from specdens import sbeta
from corrfunc import S_exact, A_exact
from scipy.linalg.interpolative import interp_decomp, reconstruct_skel_matrix, reconstruct_matrix_from_id
from scipy.optimize import nnls

def edr_id(M, N, tc, omegac, eps, frank, rand=False):
    """
    Perform frequency estimation using interpolative decomposition (ID) and NNLS.

    Parameters:
    - M (int): Number of time points.
    - N (int): Number of frequency points.
    - tc (float): Maximum time value.
    - omegac (float): Maximum frequency value.
    - temp (float): Temperature parameter.
    - eps (float): Error tolerance for ID.
    - krank (int): Rank for the ID. If negative, ID uses error tolerance.

    Returns:
    - Nsp (int): Number of estimated frequencies.
    - w (ndarray): Estimated frequencies.
    - g (ndarray): Estimated coefficients (amplitudes).
    - krank (int): Updated rank after ID.
    """
    # Generate time mesh tf and frequency mesh wf
    t, w = equispaced_mesh(M,N,tc,omegac)

    # Create core matrix Kf
    f = create_integrand(M,N,t,w)

    #f = np.loadtxt('fmat.txt')
    # Perform Interpolative Decomposition (ID)
    if frank < 0:
        frank, idx, B, err1 = id_freq_eps(f, eps, rand)
    else:
        idx, B, err1 = id_freq_rank(f, frank, rand)
        # krank remains the same
    
    print(np.max(f))
    print(np.min(f))
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


def id_freq_eps(f, eps, rnd):
    """
    Perform interpolative decomposition on matrix K with error tolerance eps.

    Parameters:
    - K: ndarray
        Input matrix of shape (M2, N).
    - eps: float
        Error tolerance for the decomposition.

    Returns:
    - krank: int
        Rank of the approximation.
    - idx: ndarray
        Indices of the selected columns.
    - B: ndarray
        Matrix containing selected columns of K.
    - err: float
        Maximum absolute error between the original and approximated matrix.
    """

    f0 = f.copy()

    # Perform the interpolative decomposition
    frank, idx, proj = interp_decomp(f0, eps, rand=rnd)

    # Extract the selected columns to form matrix B
    B = reconstruct_skel_matrix(f0, frank, idx)

    # Reconstruct the approximated matrix K1
    f1 = reconstruct_matrix_from_id(B, idx, proj)

    # Compute the maximum absolute error
    err = np.max(np.abs(f1 - f))

    return frank, idx, B, err


def id_freq_rank(f, frank, rnd):
    """
    Perform interpolative decomposition on matrix K with specified rank krank.

    Parameters:
    - K: ndarray
        Input matrix of shape (M2, N).
    - krank: int
        Desired rank for the approximation.

    Returns:
    - idx: ndarray
        Indices of the selected columns.
    - B: ndarray
        Matrix containing selected columns of K.
    - err: float
        Maximum absolute error between the original and approximated matrix.
    """

    f0 = f.copy()

    # Perform the interpolative decomposition with specified rank
    idx, proj = interp_decomp(f0, frank, rand=rnd)

    # Extract the selected columns to form matrix B
    B = reconstruct_skel_matrix(f0, frank, idx)

    # Reconstruct the approximated matrix K1
    f1 = reconstruct_matrix_from_id(B, idx, proj)

    # Compute the maximum absolute error
    err = np.max(np.abs(f1 - f))

    return idx, B, err


def edr_coef(t, B):
    """
    Use non-negative least squares (NNLS) to estimate coefficients g.

    Parameters:
    - tf (ndarray): Time array (size: M)
    - krank (int): Approximation rank
    - B (ndarray): Input matrix (size: (2*M, krank))

    Returns:
    - g (ndarray): Estimated coefficients (size: krank)
    - err (float): Estimation error
    """
    M = len(t)
    c = np.zeros(2*M)

    # Construct vector c with S_exact(t) and A_exact(t)
    for i in range(M):
        ti = t[i]
        c[i] = S_exact(ti)
    for i in range(M):
        ti = t[i]
        c[M + i] = A_exact(ti)

    # Solve the NNLS problem: minimize ||B * g - c|| subject to g >= 0
    g, err = nnls(B, c)

    return g, err


def equispaced_mesh(M, N, tc, omegac):

    # Time grid (t)
    t = np.linspace(0,tc,M)

    # Frequency grid (w)
    w = np.linspace(-omegac,omegac,N)

    return t, w


def create_integrand(M, N, t, w):
    
    f = np.zeros((2*M,N),dtype=float)

    # Fill the first M rows of K with the real part of an integrand
    for j in range(N):
        for i in range(M):
            f[i, j] = sbeta(w[j]) * np.cos(w[j] * t[i])
    
    # Fill the next M rows of K with the imaginary part of an integrand
    for i in range(M, 2*M):
        for j in range(N):
            f[i, j] = -sbeta(w[j]) * np.sin(w[j] * t[i-M])
    
    return f






    
    


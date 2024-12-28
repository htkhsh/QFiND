import numpy as np
from scipy.linalg.interpolative import interp_decomp, reconstruct_skel_matrix, reconstruct_matrix_from_id#, reconstruct_interp_matrix

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
    #P = reconstruct_interp_matrix(idx, proj)

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
    #P = reconstruct_interp_matrix(idx, proj)

    # Reconstruct the approximated matrix K1
    f1 = reconstruct_matrix_from_id(B, idx, proj)

    # Compute the maximum absolute error
    err = np.max(np.abs(f1 - f))

    return idx, B, err

import numpy as np

def lanczos(N, xw):
    """
    Lanczos algorithm.

    Given the discrete inner product whose nodes are contained 
    in the first column, and whose weights are contained in the 
    second column, of the nx2 array xw, the call ab = lanczos(N, xw)
    generates the first N recurrence coefficients ab of the 
    corresponding discrete orthogonal polynomials. The N alpha-
    coefficients are stored in the first column, the N beta-
    coefficients in the second column, of the Nx2 array ab.

    The script is adapted from the routine RKPW in
    W.B. Gragg and W.J. Harrod, "The numerically stable 
    reconstruction of Jacobi matrices from spectral data", 
    Numer. Math. 44 (1984), 317-335.
    """
    Ncap = xw.shape[0]
    if N <= 0 or N > Ncap:
        raise ValueError('N out of range')
    p0 = np.zeros(Ncap + 1)
    p1 = np.zeros(Ncap + 1)
    p0[1:] = xw[:, 0]
    p1[1] = xw[0, 1]
    for n in range(1, Ncap):
        pn = xw[n, 1]
        gam = 1.0
        sig = 0.0
        t = 0.0
        xlam = xw[n, 0]
        for k in range(1, n + 2):
            rho = p1[k] + pn
            tmp = gam * rho
            tsig = sig
            if rho <= 0:
                gam = 1.0
                sig = 0.0
            else:
                gam = p1[k] / rho
                sig = pn / rho
            tk = sig * (p0[k] - xlam) - gam * t
            p0[k] = p0[k] - (tk - t)
            t = tk
            if sig <= 0:
                pn = tsig * p1[k]
            else:
                pn = (t ** 2) / sig
            tsig = sig
            p1[k] = tmp
    ab = np.column_stack((p0[1:N + 1], p1[1:N + 1]))
    return ab

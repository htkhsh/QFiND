import numpy as np
from scipy.optimize import least_squares
import sys
from init import setpara, const, opt
from scipy.integrate import quad
from specdens import spectral_density, sbeta, powerlaw_exp, sdens
from scipy.special import gammainc

f = open(str(sys.argv[1]), mode='r')
setpara(f)

icm2ifs = const['icm2ifs']
wmax_quad = opt['wmax_quad']
if opt['stype'] == 'PWR':
    s = opt['s']
    alpha = opt['alpha']
    gamc = opt['gamc']

def dos(omega):
    return sdens(omega)/omega

def int_dens(omega):
    if opt['stype'] == 'PWR':
        return gamc * gammainc(s,omega/gamc)
    else:
        return quad(lambda w: dos(w), 0.0, omega, epsabs=1e-12)[0]
 
def func(w):
    N = len(w)
    epsilon = (1.0 / float(N)) * int_dens(wmax_quad)
    F = np.zeros(N)
    for j in range(N):
        F[j] = int_dens(w[j]) - (j - 0.5) * epsilon
    return F


def mdm_ohmic(Msp, gamc, omegac):

    if Msp % 2 != 0:
        raise ValueError("Msp must be even.")

    Msp2 = Msp // 2  # Half the number of spectral points

    # Allocate arrays
    wj = np.zeros(Msp2)
    dos = np.zeros(Msp2)
    wk = np.zeros(Msp)
    zk = np.zeros(Msp)

    # Compute omega0
    omega0 = gamc / float(Msp2) * (1.0 - np.exp(-omegac / gamc))

    # Compute wj array
    for j in range(Msp2):
        # Adjust indices for Python (0-based indexing)
        numerator = float(Msp2)
        denominator = float(Msp2 - (j + 1) + 0.5)
        wj[j] = gamc * np.log(numerator / denominator)

    # Compute dos array
    wk = np.concatenate((-wj, wj))
    dos = np.exp(-np.abs(wk) / gamc) / omega0
    zk = 1.0 / dos

    return wk ,zk


def mdm_g(Msp, Omega_max):
    # Optimization options
    x0 = np.linspace(0, Omega_max, Msp)
    lb = np.zeros(Msp)
    ub = np.full(Msp, Omega_max)
    options = {'ftol': 1e-10, 'xtol': 1e-10, 'verbose': 2}

    # Solve the equation using least_squares
    res = least_squares(func, x0, bounds=(lb, ub), method='trf', **options)
    w = res.x

    # Calculate the error
    err = func(w)
    print(f'Maximum error: {np.max(np.abs(err))}')

    # Check convergence
    if not res.success:
        print('Solution did not converge:')
        print(res.message)
    else:
        print('Solution converged.')

    # Calculate zk
    norm = (1.0 / float(Msp)) * int_dens(Omega_max)
    print('norm:', norm)
    wk = np.concatenate((-w, w))
    zk = 1.0 / dos(np.abs(wk)) * norm 

    return wk, zk


import numpy as np
from scipy.integrate import quad
import sys
from init import setpara, const, opt
from specdens import spectral_density

f = open(str(sys.argv[1]), mode='r')
setpara(f)
hbar = const['hbar']
kb = const['kb']
Temp = opt['temperature']
stype = opt['stype']
Omega_max = opt['Omega_max']
bho = hbar*1e15/kb/Temp

# Real part of a BCF
def S_exact(t):
    if Temp == 0.0:
        def integrand_0K(w):
            res = spectral_density(stype, w) * np.cos(w * t) / np.pi
            return res
        res, err = quad(integrand_0K, 0.0, Omega_max, epsabs=1e-12)
    else:
        def integrand(w):
            res = spectral_density(stype, w) / np.tanh(0.5 * bho * w) * np.cos(w * t) / np.pi
            return res
        res, err = quad(integrand, 0.0, Omega_max, epsabs=1e-12)
    return res

# Imaginary part of a BCF
def A_exact(t):
    def integrand_A(w):
        res = -spectral_density(stype, w) * np.sin(w * t) / np.pi
        return res
    res, err = quad(integrand_A, 0.0, Omega_max, epsabs=1e-12)
    return res

import numpy as np
from scipy.integrate import quad
from scipy import constants as pc
import global_value as g
from my_function import setpara 
from specdens import spectral_density

setpara()
hbar = pc.h/2.0/np.pi
kb = pc.k
bho = hbar*1e15/kb/g.temperature

# Real part of a BCF
def S_exact(t):
    if g.temperature == 0.0:
        def integrand_0K(w):
            res = spectral_density(g.stype, w) * np.cos(w * t) / np.pi
            return res
        res, err = quad(integrand_0K, 0.0, g.wmax, epsabs=1e-12)
    else:
        def integrand(w):
            res = spectral_density(g.stype, w) / np.tanh(0.5 * bho * w) * np.cos(w * t) / np.pi
            return res
        res, err = quad(integrand, 0.0, g.wmax, epsabs=1e-12)
    return res

# Imaginary part of a BCF
def A_exact(t):
    def integrand_A(w):
        res = -spectral_density(g.stype, w) * np.sin(w * t) / np.pi
        return res
    res, err = quad(integrand_A, 0.0, g.wmax, epsabs=1e-12)
    return res

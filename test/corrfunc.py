import numpy as np
from scipy.integrate import quad
import sys
from init import setpara, const, opt
from specdens import spectral_density

f = open(str(sys.argv[1]), mode='r')
setpara(f)
hbar = const['hbar']
kb = const['kb']
icm2ifs = const['icm2ifs']
Temp = opt['temperature']
stype = opt['stype']
wmax_quad = opt['wmax_quad']
if Temp > 0.0:
    bho = hbar*1e15/kb/Temp

# Real part of a BCF
def S_exact(t):
    if Temp == 0.0:
        def integrand_0K(w):
            res = spectral_density(stype, w, icm2ifs) * np.cos(w * t) / np.pi
            return res
        res, err = quad(integrand_0K, 0.0, wmax_quad*icm2ifs, epsabs=1e-12)
    else:
        def integrand(w):
            res = spectral_density(stype, w, icm2ifs) / np.tanh(0.5 * bho * w) * np.cos(w * t) / np.pi
            return res
        res, err = quad(integrand, 0.0, wmax_quad*icm2ifs, epsabs=1e-12)
    return res

# Imaginary part of a BCF
def A_exact(t):
    def integrand(w):
        res = -spectral_density(stype, w, icm2ifs) * np.sin(w * t) / np.pi
        return res
    res, err = quad(integrand, 0.0, wmax_quad*icm2ifs, epsabs=1e-12)
    return res

def C_beta(tau):
    def integrand(w):
        res = spectral_density(stype, w, icm2ifs) * np.cosh(0.5 * bho * w - w * tau) / np.sinh(0.5 * bho * w) / np.pi
        return res
    res, err = quad(integrand, 0.0, wmax_quad*icm2ifs, epsabs=1e-12)
    return res

def S_mix(t, tau):
    def integrand(w):
        res = spectral_density(stype, w, icm2ifs) * np.cosh(0.5 * bho * w - w * tau) / np.sinh(0.5 * bho * w) * np.cos(w * t) / np.pi
        return res
    res, err = quad(integrand, 0.0, wmax_quad*icm2ifs, epsabs=1e-12)
    return res

def A_mix(t, tau):
    def integrand(w):
        res = spectral_density(stype, w, icm2ifs) * np.sinh(0.5 * bho * w - w * tau) / np.sinh(0.5 * bho * w) * np.sin(w * t) / np.pi
        return res
    res, err = quad(integrand, 0.0, wmax_quad*icm2ifs, epsabs=1e-12)
    return res

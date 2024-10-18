import numpy as np
from scipy import constants as pc
import global_value as g

def spectral_density(stype, omega):
    if stype == "PWR":
        res = powerlaw_exp(omega, g.s, g.alpha, g.gamc)
    elif stype == "TMn":
        res = tannor_meyer_n(omega, g.Omg, g.Gam, g.Lam)
    elif stype == "BOn":
        res = brownian(omega, g.Omg, g.Gam, g.Lam)
    else:
        res = 0.0 
    return res


def sbeta(omega):
    hbar = pc.h/2.0/np.pi
    kb = pc.k 
    bho = hbar*1e15/kb/g.temperature
    if g.temperature == 0.0:
        res = np.sign(omega) * spectral_density(g.stype, np.abs(omega)) / np.pi
    else:
        bho = 1.0
        res = (np.sign(omega) * spectral_density(g.stype, np.abs(omega)) *
               (1.0 / np.tanh(0.5 * bho * omega) + 1.0) / (2.0 * np.pi))
    
    return res


def powerlaw_exp(omega, s, alpha, gamc):
    res = np.pi * alpha * gamc**(1.0 - s) * omega**s * np.exp(-omega / gamc)
    return res


def tannor_meyer_n(omega, Omg, Gam, Lam):
    n = len(Omg)
    res = 0.0
    for i in range(n):
        p = 4.0 * Gam[i] * Lam[i] * (Omg[i]**2 + Gam[i]**2)
        deno = ((omega + Omg[i])**2 + Gam[i]**2) * ((omega - Omg[i])**2 + Gam[i]**2)
        res += p * omega / deno
    return res


def brownian(omega, Omg, Gam, Lam):
    n = len(Omg)
    res = 0.0
    for i in range(n):
        p = 2.0 * Gam[i] * Lam[i] * Omg[i]**2
        deno = (omega**2 - Omg[i]**2)**2 + Gam[i]**2 * omega**2
        res += p * omega / deno
    return res





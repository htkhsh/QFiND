import numpy as np
import sys
from init import setpara, const, opt

f = open(str(sys.argv[1]), mode='r')
setpara(f)
icm2ifs = const['icm2ifs']
hbar = const['hbar']
kb = const['kb'] 
Temp = opt['temperature']
stype = opt['stype']
bho = hbar*1e15/kb/Temp

def spectral_density(stype, omega, nrm=1.0):
    if stype == "PWR":
        res = powerlaw_exp(omega, opt['s'], opt['alpha'], opt['gamc']*nrm)
    elif stype == "TMn":
        res = tannor_meyer_n(omega, opt['Omg']*nrm, opt['Gam']*nrm, opt['Lam']*nrm)
    elif stype == "BOn":
        res = brownian(omega, opt['Omg']*nrm, opt['Gam']*nrm, opt['Lam']*nrm)
    else:
        res = 0.0 
    return res


def sbeta(omega, nrm=1.0):
    if Temp == 0.0:
        res = np.sign(omega) * spectral_density(stype, np.abs(omega), nrm) / np.pi
    else:
        res = (np.sign(omega) * spectral_density(stype, np.abs(omega), nrm) *
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





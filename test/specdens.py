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
if Temp > 0.0:
    bho = hbar*1e15/kb/Temp

def spectral_density(stype, omega, nrm=1.0):
    if stype == "PWR":
        res = powerlaw_exp(omega, opt['s'], opt['alpha'], opt['gamc']*nrm)
    elif stype == "TMn":
        res = tannor_meyer_n(omega, opt['Omg']*nrm, opt['Gam']*nrm, opt['Lam']*nrm)
    elif stype == "BOn":
        res = brownian(omega, opt['Omg']*nrm, opt['Gam']*nrm, opt['Lam']*nrm)
    elif stype == "AtCry":
        res = atcry(omega, nrm)
    else:
        res = 0.0 
    return res

def sdens(omega, nrm=1.0):
    res = spectral_density(stype, omega, nrm)
    return res

def sbeta(omega, nrm=1.0):
    if Temp == 0.0:
        res = np.sign(omega) * spectral_density(stype, np.abs(omega), nrm) / np.pi
    else:
        res = (np.sign(omega) * spectral_density(stype, np.abs(omega), nrm) *
               (1.0 / np.tanh(0.5 * bho*icm2ifs/nrm * omega) + 1.0) / (2.0 * np.pi))
    return res

def n_be(omega, nrm=1.0):
    res = 1.0 / (np.exp(bho*icm2ifs/nrm * omega) - 1.0)
    return res

def s_be(omega, nrm=1.0):
    if Temp == 0.0:
        res = np.sign(omega) * spectral_density(stype, np.abs(omega), nrm) / np.pi
    else:
        res = (np.sign(omega) * spectral_density(stype, np.abs(omega), nrm) *
               n_be(omega, nrm) / np.pi)
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

def atcry(w, nrm):
    jw = np.loadtxt('sd_atcry.txt')
    res = jw[:, 0] + jw[:, 1] * 1.0j
    pole = jw[:, 2] + jw[:, 3] * 1.0j
    res = res / icm2ifs * nrm
    pole = pole / icm2ifs * nrm

    w = np.array(w)
    p = np.zeros_like(w, dtype=float)

    mask = w >= 1.0*nrm
    w_filtered = w[mask]

    if w_filtered.size > 0:
        p_partial = np.zeros_like(w_filtered, dtype=float)
        for n in range(len(pole)):
            p_partial += (res[n] / (w_filtered - pole[n])).real
        p_partial *= 1.31591407196976
        p[mask] = p_partial

    return p



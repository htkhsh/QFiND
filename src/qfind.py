import numpy as np
import sys
from init import setpara, const, opt
from eval import calc_error
from specdens import sbeta
from id import edr_id
from bsdo import bsdo
from logarithmic import log_disc
from mdm import mdm_ohmic, mdm_g

# Parameter setting
f = open(str(sys.argv[1]), mode='r')
setpara(f)
icm2ifs = const['icm2ifs']
N_t = opt['N_t']
Tc = opt['Tc']

# Main routine
if opt['method'] == 'ID':
    Omega_min = opt['Omega_min']
    Omega_max = opt['Omega_max']
    N_w = opt['N_w']
    eps = opt['eps']
    frank = opt['frank']
    #rand = opt['rand']
    Msp, wk, zk, krank = edr_id(N_t, N_w, Tc, Omega_min, Omega_max, eps, frank)
elif opt['method'] == 'BSDO':
    Omega_min = opt['Omega_min']
    Omega_max = opt['Omega_max']
    N_w = opt['N_w']
    Msp = opt['Msp']
    wk, zk = bsdo(N_w, Omega_min, Omega_max, Msp)
    wk = wk * icm2ifs
    zk = zk * icm2ifs**2.0
    zk = zk / sbeta(wk, icm2ifs)
elif opt['method'] == 'LOG':
    Omega_max = opt['Omega_max']
    Msp = opt['Msp']
    wk, zk = log_disc(Msp, Omega_max)
    wk = wk * icm2ifs
    zk = zk * icm2ifs**2.0
    zk = zk / sbeta(wk, icm2ifs)
elif opt['method'] == 'MDM':
    Omega_max = opt['Omega_max']
    Msp = opt['Msp']
    if opt['stype'] == "PWR" and opt['s'] == 1.0:
        gamc = opt['gamc']
        wk, zk = mdm_ohmic(Msp, gamc, Omega_max)
    else:
        wk, zk = mdm_g(Msp, Omega_max)
    wk = wk * icm2ifs
    zk = zk * icm2ifs

# Error estimation and Visualization
calc_error(wk, zk)

# Write frequencies omega_k and coefficients g_k(beta) to a file
filename = "omega_g.txt"
with open(filename, 'w') as of:
    of.write("================================================================\n")
    of.write("         omega_k [cm^-1]          g_k [cm^-1]               \n")
    of.write("================================================================\n")
    for i in range(Msp):
        of.write('{:25.15e} {:25.15e}\n'.format(
            wk[i] / icm2ifs, np.sqrt(zk[i]*sbeta(wk[i], icm2ifs)) / icm2ifs**2.0
        ))
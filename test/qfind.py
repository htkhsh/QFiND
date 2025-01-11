import numpy as np
import sys
from init import setpara, const, opt
from eval import calc_error
from specdens import sbeta, s_be
from id import edr_id
from id_beta import edr_id_beta
from id_mix import edr_id_mix
from bsdo import bsdo, sdens

# Parameter setting
f = open(str(sys.argv[1]), mode='r')
setpara(f)
hbar = const['hbar']
kb = const['kb'] 
icm2ifs = const['icm2ifs']
N_w = opt['N_w']
Omega_min = opt['Omega_min']
Omega_max = opt['Omega_max']
if opt['BCF_type'] == "Real" or opt['BCF_type'] == "Mix":
    N_t = opt['N_t']
    Tc = opt['Tc']
    print("cutoff time: ", Tc)
if opt['BCF_type'] == "Imag" or opt['BCF_type'] == "Mix":
    N_tau = opt['N_tau']
Temp = opt['temperature']
if Temp > 0.0:
    beta = hbar*1e15/kb/Temp
else:
    beta = np.inf

print("beta: ", beta)

# Main routine
if opt['method'] == 'ID':
    eps = opt['eps']
    frank = opt['frank']
    rand = opt['rand']
    if opt['BCF_type'] == 'Real':
        Msp, wk, zk, krank = edr_id(N_t, N_w, Tc, Omega_min, Omega_max, eps, frank, rand)
    elif opt['BCF_type'] == 'Imag':
        Msp, wk, zk, krank = edr_id_beta(N_tau, N_w, beta, Omega_min, Omega_max, eps, frank, rand)
    elif opt['BCF_type'] == 'Mix':
        #Msp, wk, zk, krank = edr_id(N_t, N_w, Tc, Omega_min, Omega_max, eps, frank, rand)
        Msp, wk, zk, krank = edr_id_mix(N_t, N_tau, N_w, Tc, beta, Omega_min, Omega_max, eps, frank)
elif opt['method'] == 'BSDO':
    Msp = opt['Msp']
    wk, zk = bsdo(N_w, Omega_min, Omega_max, Msp)
    wk = wk * icm2ifs
    zk = zk * icm2ifs**2.0 
    zk = zk / sbeta(wk,icm2ifs)
"""
elif opt['method'] == 'BSDO-BE':
    wk, zk = bsdo_be()
    Msp = len(wk)
    wk = wk * icm2ifs
    zk = zk * icm2ifs**2.0 
    zk = zk / s_be(wk,icm2ifs)
elif opt['method'] == 'BSDO-SD':
    wk, zk = bsdo_sd()
    Msp = len(wk)
    wk = wk * icm2ifs
    zk = zk * icm2ifs**2.0 
    zk = zk / sdens(wk, icm2ifs)
"""

# Error estimation and Visualization
calc_error(wk, zk)

# Write frequencies omega_k and coefficients g_k(beta) to a file
filename = "omega_g.txt"
with open(filename, 'w') as of:
    of.write("================================================================\n")
    of.write("         omega_k [cm^-1]          z_k [cm^-1]                   \n")
    of.write("================================================================\n")
    for i in range(Msp):
        of.write('{:25.15e} {:25.15e}\n'.format(
            wk[i] / icm2ifs, zk[i] / icm2ifs
        ))
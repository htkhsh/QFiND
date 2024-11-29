import numpy as np
import sys
from init import setpara, const, opt
from eval import calc_error
from specdens import sbeta
from edr import edr_id
from bsdo import bsdo

# Parameter setting
f = open(str(sys.argv[1]), mode='r')
setpara(f)
icm2ifs = const['icm2ifs']
N_t = opt['N_t']
N_w = opt['N_w']
Tc = opt['Tc']
Omegac = opt['Omegac']

# Main routine
if opt['method'] == 'ID':
    eps = opt['eps']
    frank = opt['frank']
    rand = opt['rand']
    Nsp, wk, zk, krank = edr_id(N_t, N_w, Tc, Omegac*icm2ifs, eps, frank, rand)
if opt['method'] == 'BSDO':
    wk, zk = bsdo()
    Nsp = len(wk)
    wk = wk * icm2ifs
    zk = zk * icm2ifs**2.0
    zk = zk / sbeta(wk, icm2ifs)

# Error estimation and Visualization
calc_error(wk, zk)

# Write frequencies omega_k and coefficients g_k(beta) to a file
filename = "omega_g.txt"
with open(filename, 'w') as of:
    of.write("================================================================\n")
    of.write("         omega_k [cm^-1]          z_k [cm^-1]                   \n")
    of.write("================================================================\n")
    for i in range(Nsp):
        of.write('{:25.15e} {:25.15e}\n'.format(
            wk[i] / icm2ifs, zk[i] / icm2ifs
        ))
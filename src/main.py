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
    gk = zk * sbeta(wk,icm2ifs)
if opt['method'] == 'BSDO':
    wk, gk = bsdo()
    Nsp = len(wk)
    wk = wk * icm2ifs
    gk = gk * icm2ifs**2.0

# Error estimation and Visualization
calc_error(wk, gk)

# Write frequencies omega_k and coefficients g_k(beta) to a file
filename = "omega_g.txt"
with open(filename, 'w') as of:
    of.write("================================================================\n")
    of.write("         omega_k [cm^-1]          g_k(beta) [cm^-1]             \n")
    of.write("================================================================\n")
    for i in range(Nsp):
        of.write('{:25.15e} {:25.15e}\n'.format(
            wk[i] / icm2ifs, np.sqrt(gk[i]) / icm2ifs
        ))
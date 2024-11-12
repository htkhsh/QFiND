import numpy as np
import sys
from init import setpara, const, opt
from edr import edr_id
from eval import calc_error
from specdens import sbeta

# Parameter setting
f = open(str(sys.argv[1]), mode='r')
setpara(f)
icm2ifs = const['icm2ifs']
M = opt['M']
N = opt['N']
Tc = opt['Tc']
Omegac = opt['Omegac']
eps = opt['eps']
frank = opt['frank']
rand = opt['rand']
# Main routine
if opt['method'] == 'ID':
    Nsp, wk, zk, krank = edr_id(M, N, Tc, Omegac*icm2ifs, eps, frank, rand)

# Error estimation and Visualization
calc_error(wk, zk)

# Write frequencies omega_k and coefficients g_k(\beta) to a file
filename = "omega_g.txt"
with open(filename, 'w') as of:
    of.write("================================================================\n")
    of.write("       omega_k[cm^-1]      z_k*sbeta(omega_j)[(cm^-1)^2]        \n")
    of.write("================================================================\n")
    for i in range(Nsp):
        of.write('{:25.15e} {:25.15e}\n'.format(
            wk[i] / icm2ifs, zk[i] * sbeta(wk[i]) / (icm2ifs ** 2.0)
        ))
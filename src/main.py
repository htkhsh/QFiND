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
Nt = opt['Nt']
Nw = opt['Nw']
Tc = opt['Tc']
Omegac = opt['Omegac']
eps = opt['eps']
frank = opt['frank']
rand = opt['rand']
# Main routine
if opt['method'] == 'ID':
    Nsp, wk, zk, krank = edr_id(Nt, Nw, Tc, Omegac*icm2ifs, eps, frank, rand)
    gk = zk * sbeta(wk,icm2ifs)

# Error estimation and Visualization
calc_error(wk, gk)

# Write frequencies omega_k and coefficients g_k(\beta) to a file
filename = "omega_g.txt"
with open(filename, 'w') as of:
    of.write("================================================================\n")
    of.write("         omega_k [cm^-1]          g_k(\beta) [cm^-1]            \n")
    of.write("================================================================\n")
    for i in range(Nsp):
        of.write('{:25.15e} {:25.15e}\n'.format(
            wk[i] / icm2ifs, np.sqrt(gk[i]) / icm2ifs
        ))
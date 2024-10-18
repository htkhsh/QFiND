import numpy as np
import matplotlib.pyplot as plt
from my_function import setpara 
import global_value as g
from edr import edr_id
from eval import calc_error
from specdens import sbeta

# Parameter setting
setpara()

# Main routine
Nsp, wk, zk, krank = edr_id(g.M, g.N, g.tc, g.omegac, g.temperature, g.eps, g.krank)

# Error estimation and Plotting
calc_error(Nsp, wk, zk)

# Write frequencies omega_k and coefficients g_k(\beta) to a file
filename = "omega_g.txt"
with open(filename, 'w') as of:
    of.write("================================================================\n")
    of.write("       omega_k[cm^-1]      z_k*sbeta(omega_j)[(cm^-1)^2]        \n")
    of.write("================================================================\n")
    for i in range(Nsp):
        of.write('{:25.15e} {:25.15e}\n'.format(
            wk[i] / g.icm2ifs, zk[i] * sbeta(wk[i]) / (g.icm2ifs ** 2.0)
        ))
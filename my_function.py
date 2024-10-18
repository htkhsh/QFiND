import numpy as np
from scipy import constants as pc
import global_value as g

def setpara():
    g.icm2ifs = pc.c*1e2*2.0*np.pi*1e-15
    g.temperature = 300.0
    g.tc = 1000.0
    g.omegac = 300.0*g.icm2ifs
    g.wmax = 1000.0*g.icm2ifs
    g.M = 500
    g.N = 2000
    g.eps = 1e-2
    g.krank = -1
    # spectral density
    g.stype = "PWR"
    # PWR
    g.s = 1.0
    g.gamc = 53.0*g.icm2ifs
    g.alpha = 5.0/53.0
    # TMn and BOn
    g.Omg = np.array([])
    

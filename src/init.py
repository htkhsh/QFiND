import numpy as np
from scipy import constants as pc
global const
global opt

const = {}
const['icm2ifs'] = pc.c*1e2*2.0*np.pi*1e-15
const['hbar'] = pc.h/2.0/np.pi
const['kb'] = pc.k

opt = {}

def setpara(f):
    data = f.readlines()
    for line in data:
        # parse input, assign values to variables
        if "=" in line:
            key, value = line.split(" = ")
            opt[key.strip()] = value.strip()
        f.close()
    # Parse standard options
    opt['temperature'] = float(opt['temperature'])
    opt['Tc'] = float(opt['Tc'])
    opt['Omegac'] = float(opt['Omegac'])
    opt['Omega_max'] = float(opt['Omega_max'])
    opt['M'] = int(opt['M'])
    opt['N'] = int(opt['N'])
    opt['eps'] = float(opt['eps'])
    opt['frank'] = int(opt['frank'])
    opt['rand'] = bool(int(opt['rand']))
    # spectral density
    if opt['stype'] == "PWR":
        opt['s'] = float(opt['s'])
        opt['alpha'] = float(opt['alpha'])/float(opt['gamc'])
        opt['gamc'] = float(opt['gamc'])
    # TMn
    if opt['stype'] == "TMn" or opt['stype'] == "BOn":
        opt['Omg'] = np.array([float(num.strip()) for num in opt['Omg'].split(',')])*const['icm2ifs']
        opt['Gam'] = np.array([float(num.strip()) for num in opt['Gam'].split(',')])*const['icm2ifs']
        opt['Lam'] = np.array([float(num.strip()) for num in opt['Lam'].split(',')])*const['icm2ifs']